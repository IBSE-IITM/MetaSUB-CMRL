# MetaSUB-CMRL
This repository contains the analysis code for the Chennai MetaSUB project. Both R and Python scripts are used for data analysis and plot creation. Please follow the sequence of scripts provided in this README to ensure that the required files are generated for smooth execution.

### Python environment setup
- Python = 3.11
- Use `requirments.txt` to install the required pakages. 
### Data folder structure
- Chennai_data : data of Chennai samples 
- MetaSUB_data.zip : data of MetaSUB (Extract it to folder `MetaSUB_data`)

### Taxonomic classification and comparative analysis of distinct microbial signatures
- `code/microbial_signatures_analysis.ipynb`

### Relative abundance plot
- `code/rel_abd_taxa_part_1_taxa_ordering.ipynb`
- `code/rel_abd_taxa_part_2_barplot.Rmd`

### species prevalance plot
- `code/species_prevalance_plot.R`
- `results/microbial_signatures/001/species_prevalence.csv`

### Metagenome-assembled genomes (MAGs) analysis
- `code/MAG_enrichment_analysis.ipynb`
- `data/Chennai_data/MAG_enrichment_analysis`
- `results/MAG_enrichment_analysis`
Here the results of enrichment analysis from Anvi'o tool is further analysed to extract the COG functionalities that are enriched in Chennai MAGs when compared against reference strain from NCBI. 


### AMR analysis
- `code/AMR_analysis.ipynb`
- `data/Chennai_data/AMR`

### Microbial diversity and composition across various surface types
- `code/surface_type_variation_part_1_Maaslin2.Rmd`
- `code/surface_type_variation_part_2.ipynb`
- `code/surface_type_variation_part_3_valcano_plot.Rmd`


