# Analysis of hCG data

This repository contains all R scripts that were used for hCG analysis.
Generally, required input files are stored in `data/[subfolder]`, while output files will be saved to `output/[subfolder]`.

Run the scripts in the following order:

1. `make_glycan_library.R` (subfolder `library`) creates MoFi-compatible glycan libraries from peptide mapping results that were exported from Byonic. 
2. `desialylate.R` (subfolder `desialylate`) creates a desialylated glycan library.
3. `assemble_subunit_glycoforms.R` (subfolder `subunits`) extracts all glycoforms from MoFi results that were obtained for a subunit. Moreover, glycoforms are ranked by abundance, which is equal to the product of (a) relative abundance (peak height) in the deconvoluted mass spectrum (i.e., the value in column `%`) and (b) the hit score calculated by MoFi (column `Hit Score`).