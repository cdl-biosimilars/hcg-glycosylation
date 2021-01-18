# Analysis of hCG data

This repository contains all R scripts that were used for hCG analysis.
Generally, required input files are stored in `data/[subfolder]`, while output files will be saved to `output/[subfolder]`.

Run the scripts in the following order

1. `make_glycan_library.R` creates MoFi-compatible glycan libraries from peptide mapping results exported from Byonic (subfolder `library`). 
2. `desialylate.R` creates a desialylated glycan library (subfolder `desialylate`).