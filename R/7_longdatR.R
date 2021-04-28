##Cystic-Fibrosis Microbiome Project (Mainz)
##longdatR (Author: Chia-Yu Chen)
##Víctor Hugo Jarquín-Díaz 15.03.2021

library(MASS)
library(bestNormalize)
library(emmeans)
library(glmmTMB)
library(lme4)
library(lmtest)
library(orddom)
library(reshape2)
library(tidyverse)

##Source longdarR package 
source("longdat_vc1.0.1_for_distribution/R/longdat_vc1.0.1.R")

##Load data for analysis
if(!exists("stool.microbiome")){
  stool.microbiome<- readRDS("CF_project/exercise-cf-intervention/data/Stool_rare_ASV.rds")
}

if(!exists("sputum.microbiome")){
  sputum.microbiome<- readRDS("CF_project/exercise-cf-intervention/data/Sput_rare_ASV.rds")
}

if(!exists("stool.metadata")){
  stool.metadata<- readRDS("CF_project/exercise-cf-intervention/data/Stool_rare_Metadata.rds")
}

if(!exists("sputum.metadata")){
  sputum.metadata<- readRDS("CF_project/exercise-cf-intervention/data/Sput_rare_Metadata.rds")
}

## Input should be a data frame with patient in the first column 

##Change based on the sample type to work with
Stool<- TRUE
Sputum<- FALSE

if(Stool){
x<- as.data.frame(stool.metadata)
y<- as.data.frame(stool.microbiome)
}

if(Sputum){
  x<- as.data.frame(sputum.metadata)
  y<- as.data.frame(sputum.microbiome)
}

##Adjust metadata
x$Patient<- rownames(x)
x$Patient<- gsub("V\\d+", "\\1", x$Patient)
x$Patient<- gsub("^10P", "\\1", x$Patient)
x$Patient<- gsub("\\D+", "\\1", x$Patient)
x$Patient<- as.numeric(x$Patient)

x%>% 
  relocate(c(Patient, Visit))->x
x$SampleID<- rownames(x)

##All taxa should have round counts
y%>%mutate_if(is.numeric, round)->y
y$SampleID<- rownames(y)

##Remove [] from colum names
colnames(y)<- gsub("\\[|\\]", "", colnames(y))

x<- left_join(x, y, by= "SampleID")
rownames(x)<- x$SampleID
x$SampleID<- NULL

##Adjust type of variables

num.vars<- colnames(x[,c(3:6,8:25)])
fac.vars<- colnames(x[,c(1,2,7,27:69)])


## as.numeric alone will likely fail if stringsAsfactors is TRUE! 
x[, num.vars] <- apply(x[, num.vars], 2,
                              function (i) as.numeric(as.character(i)))
x[,fac.vars] <- lapply(x[fac.vars] , factor)

if(Sputum){
  ##Eliminate duplicates
  x$Duplicate<- rownames(x)
  x$Duplicate<- gsub("V\\d+", "\\1", x$Duplicate)
  x$Duplicate<- gsub("^10P", "\\1", x$Duplicate)
  x$Duplicate<- gsub("^[0-9]", "\\1", x$Duplicate)
  x$Duplicate<- gsub("^[0-9]", "\\1", x$Duplicate)
  
  x%>%
  mutate(Duplicate = case_when(Duplicate == ""  ~ FALSE,
                               Duplicate == "A" ~ TRUE))%>%
    dplyr::filter(Duplicate==FALSE)%>%
    dplyr::select(-Duplicate)->x
}

##Check for complete cases
x%>% 
  group_by(Patient)%>%
  arrange(Visit, .by_group = TRUE)%>%
  summarise(n(), .groups = "keep")%>%
  dplyr::rename(n = "n()")%>%
  filter(n == 3)-> Keep

Keep<- Keep$Patient
##Select just patients in Keep
x[x$Patient%in%Keep, ]-> x
 
##Remove columns with NAs
x<- x[ , colSums(is.na(x)) == 0]


if(Stool){
  write.table(x, file = "~/CF_project/output/longdat_Stool_input.txt", sep = "\t")
  setwd("CF_project/exercise-cf-intervention/")
  ##Run the pipeline
  longdat_disc(input = "/tables/longdat_Stool_input.txt", data_type = "count", 
               test_var = "Visit", variable_col = 70, fac_var = c(1,2,7,27:67), output_tag = "longdat_disc_Stool")
 
  ### Read the results
  result.stool <- read.table("longdat_disc_Stool_result_table.txt", header = T, row.names = 1, sep = "\t")
  confound.stool<- read.table("longdat_disc_Stool_confounders.txt", header = T, row.names = 1, sep = "\t")
}

if(Sputum){
  write.table(x, file = "~/CF_project/output/longdat_Sputum_input.txt", sep = "\t")
  setwd("CF_project/exercise-cf-intervention/")
  ##Run the pipeline
  longdat_disc(input = "/table/longdat_Sputum_input.txt", data_type = "count", 
               test_var = "Visit", variable_col = 57, fac_var = c(1,2,7,13:56), output_tag = "longdat_disc_Sputum")
  
  ### Read the results
  result.sputum <- read.table("longdat_disc_Sputum_result_table.txt", header = T, row.names = 1, sep = "\t")
  confound.sputum<- read.table("longdat_disc_Sputum_confounders.txt", header = T, row.names = 1, sep = "\t")
}