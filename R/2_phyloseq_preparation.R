##Cystic-Fibrosis Microbiome Project (Mainz)
##Phyloseq object preparation
##Víctor Hugo Jarquín-Díaz 15.02.2021

##To phyloseq
library(phyloseq)
library(ggplot2)
library(dplyr)
library(Biostrings)

reRun<- FALSE

##Load data 
if(!exists("metadata.ps")){
  if(isTRUE(reRun)){
    source("R/0_Metadata_adjustment.R") ## Run the script at base directory of repository!   
  } else {
    metadata<- readRDS(file = "CF_project/exercise-cf-intervention/data/metadata_PS.rds") 
  }
}

##Load matrices
asvmat<- read.csv("~/CF_project/output/ASV_matrix.csv")
rownames(asvmat)<-asvmat$X
asvmat$X<- NULL

taxamat<- read.csv("~/CF_project/output/Taxa_matrix.csv")
rownames(taxamat)<-taxamat$X
taxamat$X<- NULL

dna<- readDNAStringSet("~/CF_project/output/ASV.fasta")

tree<- readRDS("~/CF_project/output/phylo_tree.rds")
  
##To make Phyloseq object

##1) Use the ASV matrix and transform it to "OTU table" format
##Correct sample names:
##Sample 10P3V3 is 10P2V3B after correction from Mainz data
##Sample 10P17V2 is 10P17V3A
asvmat%>%
  rename(X10P3V3 = "X10P2V3B")%>%
  rename(X10P17V2 = "X10P17V3A")->asvmat

asv<- otu_table(asvmat, taxa_are_rows = T)
sample_names(asv) <- gsub("X", "\\1", sample_names(asv)) ##eliminate the X at the beginning

##2) Use metadata dataframe and transform it to "sample data" format
sample<- sample_data(metadata.ps)
sample_names(sample) <- sample$SampleID 

##3) Use taxa matrix and transform it to "tax table" format
tax<-tax_table(as.matrix(taxamat))
sample_names(tax)

##4) Extract the phylogenetic tree
tree<- tree$tree
  
##5)Merge each element 
PS <- merge_phyloseq(asv, tax, sample, tree)
  
  saveRDS(PS, file="~/CF_project/exercise-cf-intervention/data/PhyloSeqComp.Rds")
  rm(asvmat, asv, dna, taxamat, tax, tree, sample) 
  