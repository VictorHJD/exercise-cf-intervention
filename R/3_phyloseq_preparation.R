##Cystic-Fibrosis Microbiome Project (Mainz)
##Phyloseq object preparation
##Víctor Hugo Jarquín-Díaz 15.02.2021


##Load matrices
  asvmat<- read.csv("~/CF_project/output/ASV_matrix.csv")
  rownames(asvmat)<-asvmat$X
  asvmat$X<- NULL
  taxamat<- read.csv("~/CF_project/output/Taxa_matrix.csv")
  rownames(taxamat)<-taxamat$X
  taxamat$X<- NULL
  dna<- readDNAStringSet("~/CF_project/output/ASV.fasta")
  ##Load sample data
  sample <- read.csv("~/CF_project/output//", dec=",", stringsAsFactors=FALSE)
  ##Add sample names used by Ulrike
  sample[ , "MDC_Names"] <- samplesMDC
  
  ##To phyloseq
  library(phyloseq)
  library(ggplot2)
  library(dplyr)
  
  #keep<- rownames(seqtab.nochim)
  ### Sample data includes those that didn't worked, so let's eliminate them 
  #samdata <- samdata[samdata$BeGenDiv_Name %in% keep, ]
  #rownames(samdata) <- samdata$sample_names
  
  ##To make Phyloseq object
  ##1) Use the ASV matrix and transform it to "OTU table" format
  asv<- otu_table(asvmat, taxa_are_rows = T)
  sample_names(asv)
  ##2) Use sample dataframe and transform it to "sample data" format
  sample<- sample_data(sample)
  sample_names(sample) <- sample_names(asv)
  ##3) Use taxa matrix and transform it to "tax table" format
  tax<-tax_table(as.matrix(taxamat))
  sample_names(tax)
  
  PS <- merge_phyloseq(asv, tax)
  
  ###Add phylogenetic tree
  require(ape)
  phylotree<- rtree(ntaxa(PS), rooted=TRUE, tip.label=taxa_names(PS))
  
  PS <- merge_phyloseq(asv, sample, tax, phylotree)
  
  table(sample$System, sample$Compartment) ## ---> README sample overview (previous filtering)
  
  saveRDS(PS, file="~/CF_project/output/PhyloSeqComp.Rds")
  #saveRDS(sample, file="~/CF_project/output/sample.Rds")

