# Project: Nutitional and activity intervention on patients with cystic-fibrosis 

## Repository structure:

* data
Contain pre-procesed metadata and phyloseq object required for the analysis

* R
Contain the code to: 
  - Pre-process metadata: 0_Metadata_adjustment.R
  - Process raw sequencing data (ASV inference): 1_dada2_pipeline.R
  - Merge metadata and sequencing data: 2_phyloseq_preparation.R
  - Diversity analysis: 3_data_analysis.R

* figures
Contain final figures from the analysis

* tables
Contain tables with statistics 