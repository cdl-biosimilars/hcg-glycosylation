library(tidyverse)
library(fs)



# Manual preprocessing ----------------------------------------------------

# Make sure that data/fucsia/ contains glycan_properties.csv with columns
# `name` (glycan name)
# `sia_sites` (number of putative sialylation sites)
# `sia_occupied` (number of actually sialylated sites) and
# `fuc` (TRUE if core fucose is present, otherwise FALSE)
#
# Moreover, data/fucsia may contain any number of MoFi results for the hCG dimer
# (preprocessing: remove comment lines at the top of the files).
#
# The script also processes glycan libraries in output/library and MoFi results
# for hCG subunits in output/subunits.



# Common data -------------------------------------------------------------

glycan_types <-
  read_csv("data/fucsia/glycan_properties.csv", col_types = "ciil") %>% 
  mutate(
    fuc_sites = as.integer(!is.na(fuc)),
    fuc_occupied = as.integer(fuc)
  ) %>% 
  select(!fuc)
glycan_types



# Calculate degree of fucosylation ----------------------------------------

## 1. Intact level --------------------------------------------------------

get_dimer_fucsia <- function(f) {
  df <-
    read_csv(f) %>% 
    select(ID:`Hit Score`, `hCG alpha`, `hCG beta`) %>% 
    magrittr::set_colnames(c("id", "mass", "abundance", "hit", "hit_score",
                             "hcg_alpha", "hcg_beta")) %>% 
    mutate(across(c(mass, abundance, hit_score), as.numeric)) %>% 
    filter(!is.na(hit))
  
  mod_info <- 
    df %>% 
    unite(hcg_alpha, hcg_beta, col = "glycans", sep = "/") %>% 
    mutate(glycans = str_split(glycans, "/")) %>%
    unnest(glycans) %>%
    left_join(glycan_types, by = c(glycans = "name")) %>%
    group_by(id) %>%
    summarise(across(sia_sites:fuc_occupied, sum, na.rm = TRUE))
  mod_info
  
  df %>% 
    left_join(mod_info, by = "id") %>% 
    mutate(
      fuc_weight = abundance * hit_score * (fuc_sites > 0),
      fuc_weight_occupied = fuc_weight * fuc_occupied / fuc_sites,
      sia_weight = abundance * hit_score * (sia_sites > 0),
      sia_weight_occupied = sia_weight * sia_occupied / sia_sites
    ) %>%
    summarise(
      fuc = sum(fuc_weight_occupied, na.rm = TRUE) / sum(fuc_weight),
      sia = sum(sia_weight_occupied, na.rm = TRUE) / sum(sia_weight)
    ) %>% 
    ungroup()
}

dimer_results <- 
  dir_ls("data/fucsia", regexp = "Dimer") %>%
  map_dfr(get_dimer_fucsia, .id = "file") %>% 
  mutate(
    level = "dimer",
    batch = str_match(file, "BA(\\d+)_")[,2],
      protein = "hcg",
      sites = NA_character_
  ) %>%
  relocate(c(level, batch, protein, sites), .after = file)
dimer_results



## 2. Subunit level -------------------------------------------------------

get_subunit_fucsia <- function(f) {
  df <- read_csv(f)
  
  mod_info <- 
    df %>% 
    mutate(glycans = str_split(Name, "/")) %>%
    unnest(glycans) %>%
    left_join(glycan_types, by = c(glycans = "name")) %>%
    group_by(Name) %>%
    summarise(across(sia_sites:fuc_occupied, sum, na.rm = TRUE))
  mod_info
  
  df %>% 
    left_join(mod_info, by = "Name") %>%
    mutate(
      fuc_weight = `Fractional Abundance` * (fuc_sites > 0),
      fuc_weight_occupied = fuc_weight * fuc_occupied / fuc_sites,
      sia_weight = `Fractional Abundance` * (sia_sites > 0),
      sia_weight_occupied = sia_weight * sia_occupied / sia_sites
    ) %>%
    summarise(
      fuc = sum(fuc_weight_occupied, na.rm = TRUE) / sum(fuc_weight),
      sia = sum(sia_weight_occupied, na.rm = TRUE) / sum(sia_weight)
    ) %>% 
    ungroup()
}

subunit_results <- 
  dir_ls("output/subunits", regexp = "(Alpha|Beta).*\\.csv") %>% 
  map_dfr(get_subunit_fucsia, .id = "file") %>% 
  mutate(
    level = "subunit",
    batch = str_match(file, "BA(\\d+)_")[,2],
    protein = case_when(str_detect(file, "Alpha") ~ "hcg_alpha",
                        TRUE ~ "hcg_beta"),
    sites = NA_character_
  ) %>% 
  relocate(c(level, batch, protein, sites), .after = file)
subunit_results
  


## 3. Glycan library level ------------------------------------------------

get_library_fucsia <- function(f) {
  name_replacements <- c(
    "none" = "unmodified",
    "Core 1+Neu5Ac or GalNac-6S-3G or GalNac-3SG" = "1 x Core + 1 x S",
    "Core 1+2 Neu5Ac or GalNac-6S-3SG" = "1 x Core + 2 x S",
    "5 HexNAc, 6 Hex" = "A3G3",
    "6 HexNAc, 7 Hex, 1 Neu5Ac" = "A4S1G3",
    "6 HexNAc, 7 Hex, 4 Neu5Ac" = "A4S4",
    "4 HexNAc, 5 Hex, 2 Neu5Ac, 2 Fuc" = "A2S2F2",
    "Core 1" = "1 x Core",
    "2 HexNAc, 2 Hex" = "2 x Core",
    "2 HexNAc, 2 Hex, 1 Neu5Ac" = "2 x Core + 1 x S",
    "2 HexNAc, 2 Hex, 2 Neu5Ac" = "2 x Core + 2 x S",
    "2 HexNAc, 2 Hex, 3 Neu5Ac" = "2 x Core + 3 x S"
  )
  
  df <-
    read_csv(f) %>% 
    filter(!Sites %in% c("N13_N30", "N52_N30")) %>% 
    mutate(
      Name = Name %>%
        coalesce(Composition) %>%
        recode(!!!name_replacements)
    ) %>% 
    left_join(glycan_types, by = c(Name = "name"))
  
  df %>%
    group_by(batch, protein, sites = Sites) %>% 
    mutate(
      fuc_weight = mean * (fuc_sites > 0),
      fuc_weight_occupied = fuc_weight * fuc_occupied / fuc_sites,
      sia_weight = mean * (sia_sites > 0),
      sia_weight_occupied = sia_weight * sia_occupied / sia_sites
    ) %>%
    summarise(
      fuc = sum(fuc_weight_occupied, na.rm = TRUE) / sum(fuc_weight),
      sia = sum(sia_weight_occupied, na.rm = TRUE) / sum(sia_weight)
    ) %>% 
    ungroup() %>% 
    mutate(batch = str_glue("0{batch}"))
}

library_results <- 
  dir_ls("output/library", regexp = "library") %>% 
  map_dfr(get_library_fucsia, .id = "file") %>% 
  mutate(level = "library", batch = as.character(batch)) %>% 
  relocate(level, .after = file)
library_results



# Combine results and export ----------------------------------------------

bind_rows(list(dimer_results, subunit_results, library_results)) %>% 
  arrange(level, batch, protein, sites) %>% 
  filter(!(is.na(sites) & is.na(fuc) & is.na(sia))) %>%
  write_csv("output/fucsia/degree_fuc_sia.csv")
