library(tidyverse)
library(fs)


# Manual preprocessing ----------------------------------------------------

# Make sure that data/subunits/ contains any number of CSV files:
# * format: exported MoFi results
# * file name must contain the string `sialidase` for desialylated samples
# * manually remove all comment lines at the beginning of each file,
#   i.e., those lines starting with `#`. (One column name also contains a
#   hash character, but is unquoted; hence, readr::read_csv(..., comment = "#")
#   will not work.)



# Functions ---------------------------------------------------------------

# uncomment for testing purposes
# file <- "data/subunits/Hit ox_Ovitrelle_BA056714_Alpha_180121.csv"
# desialylated <- FALSE

rank_glycoforms <- function(file, desialylated = FALSE) {
  if (desialylated)
    composition_glue <- "{Hex} Hex, {HexNAc} HexNAc, {Fuc} Fuc"
  else
    composition_glue <- "{Hex} Hex, {HexNAc} HexNAc, {Neu5Ac} Neu5Ac, {Fuc} Fuc"
  
  df <-
    read_csv(file) %>%
    filter(!is.na(Hit))
  
  df %>% 
    unite(matches("(N|O)\\d"), col = "name", sep = "/") %>% 
    mutate(
      score = as.numeric(`%`) * as.numeric(`Hit Score`),
      Composition = str_glue(composition_glue)
    ) %>%
    group_by(name) %>% 
    summarise(
      score = sum(score),
      Composition = first(Composition),
      .groups = "drop"
    ) %>%
    mutate(fractional_abundance = score / sum(score) * 100) %>%
    arrange(desc(fractional_abundance)) %>%
    mutate("Cumulative Abundance" = cumsum(fractional_abundance)) %>% 
    select(
      Name = name,
      Composition,
      "Fractional Abundance" = fractional_abundance,
      `Cumulative Abundance`
    )
}



# Analysis ----------------------------------------------------------------

dir_ls("data/subunits/", type = "file") %>% 
  walk(
    function(file) {
      output_file <-
        file %>% 
        path_file() %>% 
        path_ext_remove() %>% 
        {str_glue("output/subunits/{.}_glycoforms.csv")}
      
      desialylated <- str_detect(file, "sialidase")
      
      rank_glycoforms(file, desialylated) %>%
        write_csv(output_file)
    }
  )
