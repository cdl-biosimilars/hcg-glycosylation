library(tidyverse)
library(fs)


# Manual preprocessing ----------------------------------------------------

# Make sure that data/desialylate/ contains the following files:
# * desialylated_glycan_names.csv - table with four columns
#                                   `name` (glycan name),
#                                   `HexNAc` (number of HexNAc residues),
#                                   `Hex` (number of hexose residues), and
#                                   `Fuc` (number of fucose residues)
#
# * any number of CSV files whose name starts with "Glycan" and which contain
#   viable MoFi glycan libraries (i.e., the following columns exist: `Name`,
#   `Composition`, `Sites`, and `Relative Abundance`)



# Common data -------------------------------------------------------------

nomenclature <- read_csv(
  "data/desialylate/desialylated_glycan_names.csv",
  col_types = "ciii"
)



# Functions ---------------------------------------------------------------

# # uncomment for testing purposes
# df <- read_csv("data/desialylate/Glycan Library_Ovitrelle_BA056714_Alpha_180121.csv")

desialylate <- function(df) {
  df %>% 
    mutate(
      sugars = str_match_all(Composition, "(\\d+) (\\w+)") %>% 
        map_dfr(
          ~magrittr::set_colnames(., c("match", "count", "residue")) %>% 
            as_tibble() %>%
            select(!match) %>% 
            pivot_wider(names_from = "residue", values_from = "count")
        ) %>%
        select(!any_of("NA")) %>%
        mutate(across(everything(), ~as.integer(.) %>% replace_na(0L)))
    ) %>% 
    unpack(sugars) %>%
    group_by(Sites, HexNAc, Hex, Fuc) %>%
    summarise(
      Name = str_c(Name, collapse = "_"),
      `Relative Abundance` = sum(`Relative Abundance`),
      .groups = "drop"
    ) %>%
    group_by(Sites) %>%
    mutate(
      `Relative Abundance` = `Relative Abundance` / sum(`Relative Abundance`) * 100,
      Composition = str_glue("{HexNAc} HexNAc, {Hex} Hex, {Fuc} Fuc")
    ) %>%
    left_join(nomenclature, by = c("HexNAc", "Hex", "Fuc")) %>%
    select(Name = name, Composition, Sites, `Relative Abundance`) %>%
    {.}
}



# Analysis ----------------------------------------------------------------

files <- dir_ls("data/desialylate/", type = "file", regex = "Glycan")

files %>% walk(
  function(file) {
    output_file <-
      file %>% 
      path_file() %>% 
      path_ext_remove() %>% 
      {str_glue("output/desialylate/{.}_desialylated.csv")}
    read_csv(file) %>%
      desialylate() %>%
      write_csv(output_file)
  }
)
