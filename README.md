# Project: Diet and activity intervention on patients with cystic-fibrosis 

[[_TOC_]]
# Project status: 
- [X] planned
- [X] running
- [ ] manuscript writing
- [ ] submitted
- [ ] revision
- [ ] closed



## Deadlines


# Project description
(Unpdate later)

## Key questions
1. Are the excercise routine and a diet intervention impacting the microbiome of patients with cystic fibrosis?

    a. Can be observed in the faecal microbiota?

    b. Can be observed in the sputum micobiota?

# Team 
| Name | Affiliation | Responsibility | Position | Area |
| ------ | ------ | ------ | ------ | ------ |
| [Knoll, Rebecca](mailto:rebknoll@uni-mainz.de.ID?subject=SUBJECT%20CFProject) | University of Mainz | Data collection | PhD | Clinic |
| [Klopp, Jonas](mailto:jonklopp@uni-mainz.de.ID?subject=SUBJECT%20CFProject) | University of Mainz | Sequencing | PhD | Bioinformatics |
| [Markó, Lajos](mailto:lajosmarko@yahoo.com.ID?subject=SUBJECT%20CFProject) | ECRC, Charité, MDC | Data analysis | Postdoc | Clinic |
| [Jarquín-Díaz, Víctor Hugo](mailto:VictorHugo.JarquinDiaz@mdc-berlin.de.ID?subject=SUBJECT%20CFProject) | MDC | Data analysis | Postdoc | Bioinformatics |
| [Poplawska, Krystyna](mailto:krystyna.poplawska@unimedizin-mainz.de.ID?subject=SUBJECT%20CFProject) | University of Mainz | Project design, Supervision | PI | Clinic |
| [Forslund, Sofia](mailto:Sofia.Forslund@mdc-berlin.de.ID?subject=SUBJECT%20CFProject) | MDC | Data analysis, Supervision | PI | Bioinformatics |

# Data
* data types

| File or Directory | System | Path | short desciption |
| -------- | -------- | -------- | ---------| 
| Commed_Subset/ | VM    |/fast/AG_Forslund/Victor/CF_project/Seqdata/| all raw sequencing data  |


## Status

| Task | Status | Date | Responsible person|
| -------- | -------- | -------- | -------- |
| Data analysis | In progeress | 18.02.2021 | Víctor Hugo Jarquín-Díaz |

## Software

| Tool | Version | short description | 
| -------- | -------- | -------- |
| dada2 | 1.18 | Sequencing preprocessing | 


## Repository structure

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
