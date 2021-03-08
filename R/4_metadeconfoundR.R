##Cystic-Fibrosis Microbiome Project (Mainz)
##metadeconfoundR
##Víctor Hugo Jarquín-Díaz 18.02.2021

library(metadeconfoundR)
library(tidyverse)

##Load data

if(!exists("stool.microbiome")){
  stool.microbiome<- readRDS("CF_project/exercise-cf-intervention/data/Stool_rare_ASV.rds")
}

if(!exists("sputum.microbiome")){
  sputum.microbiome<- readRDS("CF_project/exercise-cf-intervention/data/Stool_rare_ASV.rds")
}

stool.metadata<- readRDS("CF_project/exercise-cf-intervention/data/Stool_rare_Metadata.rds")
sputum.metadata<- readRDS("CF_project/exercise-cf-intervention/data/Sput_rare_Metadata.rds")

##Subset ASV tables by Visit combination 

#1) Merge ASVs with visit 
##Take SampleID in the right order
SampleID<- rownames(stool.metadata)

##Transform metadata to dataframe 
x<- as.data.frame(stool.metadata)
y<- as.data.frame(stool.microbiome)

##Add IDs as column
x$SampleID<- SampleID
y$SampleID<- SampleID

x%>%
  select(SampleID, Visit)%>%
  left_join(y, by= "SampleID")%>%
  filter(Visit%in%c(1,2))-> y.V1V2

rownames(y.V1V2)<- y.V1V2$SampleID
y.V1V2$SampleID<- NULL
y.V1V2$Visit<- NULL

x[x$Visit %in% c(1,2), ]-> x.V1V2



