library(tidyverse)
library(Biostrings)


# Manual preprocessing ----------------------------------------------------

# Make sure that data/library/ contains the following files:
# * hcg_alpha.fasta - FASTA sequence of hCG alpha chain
# * hcg_beta.fasta  - FASTA sequence of hCG beta chain
# * glycan_nomenclature.csv - table with three columns
#                            `name` (glycan name),
#                            `type` (N or O), and
#                            `composition` (Byonic composition string)
#
# Convert the XLSX files in data/library/original/ to CSV format:
# * Delete empty column `MS alias`
# * Rename columns to `protein`, `sequence`, `glycan`, [sample names unchanged]
# 
# Generated files:
# * 150121 Byonic output N-glycan noNeuGc_noOx_noA4S1G3_v4 1miss.xlsx
#   --> byonic_results_N.csv
# * Byonic_output O_glycans_12_01_21.xlsx
#   --> byonic_results_O_56714.csv
#   |-> byonic_results_O_59433.csv



# Common data -------------------------------------------------------------

hcg_sequence <- 
  c("data/library/hcg_alpha.fasta", "data/library/hcg_beta.fasta") %>% 
  readAAStringSet() %>%
  as.character() %>%
  set_names(c("hcg_alpha", "hcg_beta"))

hcg_glycosylation_sites <- list(
  N = c(
    N52 = 52L,  # alpha
    N78 = 78L,  # alpha
    N13 = 13L,  # beta
    N30 = 30L   # beta
  ),
  O = c(
    S121 = 121L,  # all beta
    S127 = 127L,
    S130 = 130L,
    S132 = 132L,
    S138 = 138L,
    T140 = 140L
  )
)

name_regex <- regex(
  "
  (?: HexNac\\((\\d+)\\) )?
  (?: Hex   \\((\\d+)\\) )?
  (?: Neu5Ac\\((\\d+)\\) )?
  (?: Neu5Gc\\((\\d+)\\) )?
  (?: Fuc   \\((\\d+)\\) )?",
  comments = TRUE
)

hcg_glycan_nomenclature <-
  read_csv("data/library/glycan_nomenclature.csv") %>% 
  mutate(
    composition = composition %>%
      str_match(name_regex) %>%
      magrittr::set_colnames(
        c("match", "HexNAc", "Hex", "Neu5Ac", "Neu5Gc", "Fuc")
      ) %>%
      as_tibble() %>%
      select(-match) %>%
      mutate(across(everything(), ~replace_na(., "0") %>% as.integer()))
  )
hcg_glycan_nomenclature



# Functions ---------------------------------------------------------------

# # uncomment for testing purposes
# byonic_output <- "data/library/byonic_results_N.csv"
# protein_sequence <- hcg_sequence
# glycosylation_sites <- hcg_glycosylation_sites[["N"]]
# glycan_nomenclature <-  hcg_glycan_nomenclature

# prepares a glycan library
make_glycan_library <- function(byonic_output,
                                protein_sequence = hcg_sequence,
                                glycosylation_sites = hcg_glycosylation_sites,
                                glycan_nomenclature = hcg_glycan_nomenclature,
                                type = c("N", "O")) {
  type <- match.arg(type)
  glycosylation_sites <- glycosylation_sites[[type]]
  
  byonic_regex <- regex(
    "
    (?: HexNac\\((\\d+)\\) )?
    (?: Hex   \\((\\d+)\\) )?
    (?: NeuAc \\((\\d+)\\) )?
    (?: NeuGc \\((\\d+)\\) )?
    (?: Fuc   \\((\\d+)\\) )?",
    comments = TRUE
  )
  
  tidy_peptides <-
    read_csv(byonic_output) %>%
    fill(protein, sequence) %>%
    pivot_longer(
      !protein:glycan,
      names_to = "sample",
      values_to = "abundance"
    ) %>%
    extract(sample, into = c("batch", "sample"), regex = "BA(\\d+)_.+_(\\d+)$") %>% 
    replace_na(list(glycan = "none")) %>% 
    filter(!is.na(abundance))
  tidy_peptides
  
  peptides_with_sites <- 
    tidy_peptides %>% 
    mutate(
      protein = recode(
        protein,
        "sp|P01215|CGA 25-116" = "hcg_alpha",
        "tr|A0A0F7RQP8|CGB3 21-165" = "hcg_beta"
      ),
      location = str_locate(protein_sequence[protein], sequence) %>%
        as_tibble(),
    ) %>%
    unpack(location) %>%
    mutate(
      site = map2_chr(
        start,
        end,
        ~glycosylation_sites[between(glycosylation_sites, .x, .y)] %>%
          names() %>%
          str_c(collapse = "_")
      )
    )
  peptides_with_sites
  
  peptides_with_composition <-
    peptides_with_sites %>% 
    mutate(
      composition = glycan %>% 
        str_match_all(byonic_regex) %>%
        map_dfr(
          ~magrittr::set_colnames(
            .,
            c("match", "HexNAc", "Hex", "Neu5Ac", "Neu5Gc", "Fuc")
          ) %>%
            as_tibble() %>% 
            select(-match) %>% 
            mutate(across(everything(), as.integer)) %>% 
            summarise(across(everything(), sum, na.rm = TRUE))
        )
    )
  peptides_with_composition
  
  sites_with_abundance <- 
    peptides_with_composition %>% 
    group_by(batch, sample, protein, site, composition) %>%
    summarise(abundance = sum(abundance)) %>% 
    mutate(abundance = abundance / sum(abundance) * 100)
  sites_with_abundance

  sites_with_names <- 
    sites_with_abundance %>% 
    mutate(
      mofi_composition = composition %>%
        mutate(
          across(
            everything(),
            ~if_else(.x != 0L, str_glue("{.x} {cur_column()}"), NA_character_)
          )
        ) %>%
        unite(everything(), col = "mofi_string", sep = ", ", na.rm = TRUE) %>%
        deframe()
    ) %>% 
    left_join(
      glycan_nomenclature %>% filter(type == {{type}}),
      by = "composition"
    )
  sites_with_names
  
  sites_with_names %>% 
    select(
      Name = name, Composition = mofi_composition,
      Sites = site, Abundance = abundance
    )
}

# adds mean and sd for batches
make_final_dataset <- function(df) {
  df %>% 
    pivot_wider(names_from = sample, values_from = Abundance,
                names_prefix = "sample_") %>% 
    arrange(batch, protein, Sites) %>% 
    rowwise() %>% 
    mutate(
      mean = mean(c_across(starts_with("sample")), na.rm = TRUE),
      sd = sd(c_across(starts_with("sample")), na.rm = TRUE)
    )

}



# Analysis ----------------------------------------------------------------

samples <- tribble(
  ~infile, ~type, ~outfile,
  "byonic_results_N.csv", "N", "library_N.csv",
  "byonic_results_O_56714.csv", "O", "library_O_56714.csv",
  "byonic_results_O_59433.csv", "O", "library_O_59433.csv"
)

pwalk(
  samples,
  function(infile, type, outfile) {
    message("Analyzing ", infile)
    make_glycan_library(str_glue("data/library/{infile}"), type = type) %>% 
      make_final_dataset() %>% 
      write_csv(str_glue("output/library/{outfile}"))
  }
)

