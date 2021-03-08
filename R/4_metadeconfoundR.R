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

#######################1) Stool###########################
##Take SampleID in the right order
SampleID<- rownames(stool.metadata)

##Transform metadata to dataframe 
x<- as.data.frame(stool.metadata)
y<- as.data.frame(stool.microbiome)

##Add IDs as column
x$SampleID<- SampleID
y$SampleID<- SampleID

##Subset V1 vs V2 
x%>%
  select(SampleID, Visit)%>%
  left_join(y, by= "SampleID")%>%
  filter(Visit%in%c(1,2))-> y.V1V2 ## Microbial ASV counts

rownames(y.V1V2)<- y.V1V2$SampleID
y.V1V2$SampleID<- NULL
y.V1V2$Visit<- NULL

x[x$Visit %in% c(1,2), ]-> x.V1V2 ## Metadata

##Transform to binary and re-structure to have in the wright order
##Visit will change to Status and 0= V1 and 1= V2
x.V1V2%>%
  mutate(Visit = case_when(Visit == 1  ~ 0,
                           Visit == 2 ~ 1))%>%
  mutate(sex = case_when(sex == 1  ~ 0,
                         sex == 2 ~ 1))%>%
  relocate(Visit)%>%
  dplyr::rename(Status= Visit)-> x.V1V2

x.V1V2$SampleID<- NULL

##Subset V2 vs V3 
x%>%
  select(SampleID, Visit)%>%
  left_join(y, by= "SampleID")%>%
  filter(Visit%in%c(2,3))-> y.V2V3 ## Microbial ASV counts

rownames(y.V2V3)<- y.V2V3$SampleID
y.V2V3$SampleID<- NULL
y.V2V3$Visit<- NULL

x[x$Visit %in% c(2, 3), ]-> x.V2V3 ## Metadata

##Transform to binary and re-structure to have in the wright order
##Visit will change to Status and 0= V2 and 1= V3
x.V2V3%>%
  mutate(Visit = case_when(Visit == 2  ~ 0,
                           Visit == 3 ~ 1))%>%
  mutate(sex = case_when(sex == 1  ~ 0,
                         sex == 2 ~ 1))%>%
  relocate(Visit)%>%
  dplyr::rename(Status= Visit)-> x.V2V3

x.V2V3$SampleID<- NULL

##Subset V1 vs V3 
x%>%
  select(SampleID, Visit)%>%
  left_join(y, by= "SampleID")%>%
  filter(Visit%in%c(1,3))-> y.V1V3 ## Microbial ASV counts

rownames(y.V1V3)<- y.V1V3$SampleID
y.V1V3$SampleID<- NULL
y.V1V3$Visit<- NULL

x[x$Visit %in% c(1,3), ]-> x.V1V3 ## Metadata

##Transform to binary and re-structure to have in the wright order
##Visit will change to Status and 0= V1 and 1= V3
x.V1V3%>%
  mutate(Visit = case_when(Visit == 1  ~ 0,
                           Visit == 3 ~ 1))%>%
  mutate(sex = case_when(sex == 1  ~ 0,
                         sex == 2 ~ 1))%>%
  relocate(Visit)%>%
  dplyr::rename(Status= Visit)-> x.V1V3

x.V1V3$SampleID<- NULL

##Use severity as case-control variable
x%>%
  mutate(sex = case_when(sex == 1  ~ 0,
                         sex == 2 ~ 1))%>%
  relocate(Severity)%>%
  dplyr::rename(Status= Severity)-> x.all

x.all$SampleID<- gsub("V\\d+", "\\1", x.all$SampleID)
x.all$SampleID<- gsub("^10P", "\\1", x.all$SampleID)
x.all$SampleID<- gsub("\\D+", "\\1", x.all$SampleID)
x.all$SampleID
x.all$SampleID<- as.numeric(x.all$SampleID)
x.all%>%
  rename(SampleID= "Patient")-> x.all

y.all<- y
y.all$SampleID<- NULL


###Let's run the pipeline!!!

MD.V1V2<- MetaDeconfound(featureMat = y.V1V2, metaMat = x.V1V2, nnodes = 25)
MD.V2V3<- MetaDeconfound(featureMat = y.V2V3, metaMat = x.V2V3, nnodes = 25)
MD.V1V3<- MetaDeconfound(featureMat = y.V1V3, metaMat = x.V1V3, nnodes = 25)
MD.all<- MetaDeconfound(featureMat = y.all, metaMat = x.all, nnodes = 25, randomVar = list("+ (1 | Patient)" , c("Patient")))

##Plot the results 

Cun.V1V2<- BuildHeatmap(MD.V1V2, cuneiform = T, coloring = 1)
Cun.V2V3<- BuildHeatmap(MD.V2V3, cuneiform = T, coloring = 1)
Cun.V1V3<- BuildHeatmap(MD.V1V3, cuneiform = T, coloring = 1)
Cun.all<- BuildHeatmap(MD.all, cuneiform = T, coloring = 1)

BuildHeatmap(MD.V1V2)
BuildHeatmap(MD.V2V3)
BuildHeatmap(MD.V1V3)
BuildHeatmap(MD.all)

##Edit cun plots

Cun.all+
  xlab("Variables")+
  ylab("ASVs Genus-level")+
  #geom_point(aes(color = cyl, size = , shape = gear))+
  labs(tag = "A)", caption = "Using all data with patient as random factor", 
       fill= "Effect Size \n (Cliff's Delta)", shape= "Confounding status")+
  theme(legend.key = element_rect(color = "black"))+
  theme(axis.text.x = element_text(size = 10, angle = 90))+
  scale_fill_continuous(guide = "colourbar") +
  scale_size(guide = "legend")-> A

Cun.V1V2+
  xlab("Variables")+
  ylab("ASVs Genus-level")+
  labs(title="MetaDeconfoundR summarizing coneiform plot (Stool V1 vs V2)",
       tag = "B)", caption = "Using all data with patient as random factor", 
       fill= "Effect Size \n (Cliff's Delta)", shape= "Confounding status")+
  theme(legend.key = element_rect(color = "black"))+
  theme(axis.text.x = element_text(size = 10, angle = 90, face="bold", color="black"),
        axis.text.y = element_text(size = 9, face="italic", color="black"))+
  scale_x_discrete(labels=c( "Dist"= "Distance","pFVC_Response" = "pFVC Response", "Polyethylenglykol_Movicol" = "Polyeth-Movicol", "Peak_power" = "Peak power",
                            "FFM_Luk" = "FFM (Lukaski %)", "FFM_Charatsi" = "FFM (Charatsi kg)",
                   "electrolyte_supp" = "Electrolyte Supp", "FFM_Response" = "FFM Response", "sex" = "Sex", 
                   "Sport_Response" = "Sport Response"))+
  scale_y_discrete(labels=c( "CAG-56"= "Firmicutes CAG-56","[Ruminococcus]_gnavus_group" = "Ruminococcus gnavus group", 
                             "DTU089" = "Ruminococcus DTU089", "[Clostridium]_innocuum_group" = "Clostridium innocuum group",
                             "[Eubacterium]_ventriosum_group" = "Eubacterium ventriosum group"))-> B

png("CF_project/exercise-cf-intervention/figures/Q4_Deconfound_Stool_V1V2.png", units = 'in', res = 300, width=10, height=8)
B
dev.off()
#######################2) Sputum####################
##Take SampleID in the right order
SampleID<- rownames(sputum.metadata)

##Transform metadata to dataframe 
x<- as.data.frame(sputum.metadata)
y<- as.data.frame(sputum.microbiome)

##Add IDs as column
x$SampleID<- SampleID
y$SampleID<- SampleID

##Subset V1 vs V2 
x%>%
  select(SampleID, Visit)%>%
  left_join(y, by= "SampleID")%>%
  filter(Visit%in%c(1,2))-> y.V1V2 ## Microbial ASV counts

rownames(y.V1V2)<- y.V1V2$SampleID
y.V1V2$SampleID<- NULL
y.V1V2$Visit<- NULL

x[x$Visit %in% c(1,2), ]-> x.V1V2 ## Metadata

##Transform to binary and re-structure to have in the wright order
##Visit will change to Status and 0= V1 and 1= V2
x.V1V2%>%
  mutate(Visit = case_when(Visit == 1  ~ 0,
                           Visit == 2 ~ 1))%>%
  mutate(sex = case_when(sex == 1  ~ 0,
                         sex == 2 ~ 1))%>%
  relocate(Visit)%>%
  dplyr::rename(Status= Visit)-> x.V1V2

x.V1V2$SampleID<- NULL

##Subset V2 vs V3 
x%>%
  select(SampleID, Visit)%>%
  left_join(y, by= "SampleID")%>%
  filter(Visit%in%c(2,3))-> y.V2V3 ## Microbial ASV counts

rownames(y.V2V3)<- y.V2V3$SampleID
y.V2V3$SampleID<- NULL
y.V2V3$Visit<- NULL

x[x$Visit %in% c(2, 3), ]-> x.V2V3 ## Metadata

##Transform to binary and re-structure to have in the wright order
##Visit will change to Status and 0= V2 and 1= V3
x.V2V3%>%
  mutate(Visit = case_when(Visit == 2  ~ 0,
                           Visit == 3 ~ 1))%>%
  mutate(sex = case_when(sex == 1  ~ 0,
                         sex == 2 ~ 1))%>%
  relocate(Visit)%>%
  dplyr::rename(Status= Visit)-> x.V2V3

x.V2V3$SampleID<- NULL

##Subset V1 vs V3 
x%>%
  select(SampleID, Visit)%>%
  left_join(y, by= "SampleID")%>%
  filter(Visit%in%c(1,3))-> y.V1V3 ## Microbial ASV counts

rownames(y.V1V3)<- y.V1V3$SampleID
y.V1V3$SampleID<- NULL
y.V1V3$Visit<- NULL

x[x$Visit %in% c(1,3), ]-> x.V1V3 ## Metadata

##Transform to binary and re-structure to have in the wright order
##Visit will change to Status and 0= V1 and 1= V3
x.V1V3%>%
  mutate(Visit = case_when(Visit == 1  ~ 0,
                           Visit == 3 ~ 1))%>%
  mutate(sex = case_when(sex == 1  ~ 0,
                         sex == 2 ~ 1))%>%
  relocate(Visit)%>%
  dplyr::rename(Status= Visit)-> x.V1V3

x.V1V3$SampleID<- NULL

##Use severity as case-control variable
x%>%
  mutate(sex = case_when(sex == 1  ~ 0,
                         sex == 2 ~ 1))%>%
  relocate(Severity)%>%
  dplyr::rename(Status= Severity)-> x.all

x.all$SampleID<- gsub("V\\d+", "\\1", x.all$SampleID)
x.all$SampleID<- gsub("^10P", "\\1", x.all$SampleID)
x.all$SampleID<- gsub("\\D+", "\\1", x.all$SampleID)
x.all$SampleID
x.all$SampleID<- as.numeric(x.all$SampleID)
x.all%>%
  rename(SampleID= "Patient")-> x.all

y.all<- y
y.all$SampleID<- NULL


###Let's run the pipeline!!!

MD.V1V2<- MetaDeconfound(featureMat = y.V1V2, metaMat = x.V1V2, nnodes = 25)
MD.V2V3<- MetaDeconfound(featureMat = y.V2V3, metaMat = x.V2V3, nnodes = 25)
MD.V1V3<- MetaDeconfound(featureMat = y.V1V3, metaMat = x.V1V3, nnodes = 25)
MD.all<- MetaDeconfound(featureMat = y.all, metaMat = x.all, nnodes = 25, randomVar = list("+ (1 | Patient)" , c("Patient")))

##Plot the results 

Cun.V1V2<- BuildHeatmap(MD.V1V2, cuneiform = T, coloring = 1)
Cun.V2V3<- BuildHeatmap(MD.V2V3, cuneiform = T, coloring = 1)
Cun.V1V3<- BuildHeatmap(MD.V1V3, cuneiform = T, coloring = 1)
Cun.all<- BuildHeatmap(MD.all, cuneiform = T, coloring = 1)

BuildHeatmap(MD.V1V2)
BuildHeatmap(MD.V2V3)
BuildHeatmap(MD.V1V3)
BuildHeatmap(MD.all)

##Edit cun plots

Cun.V1V2+
  xlab("Variables")+
  ylab("ASVs Genus-level")+
  #geom_point(aes(color = cyl, size = , shape = gear))+
  labs(tag = "A)", caption = "Using all data with patient as random factor", 
       fill= "Effect Size \n (Cliff's Delta)", shape= "Confounding status")+
  theme(legend.key = element_rect(color = "black"))+
  theme(axis.text.x = element_text(size = 10, angle = 90))+
  scale_fill_continuous(guide = "colourbar") +
  scale_size(guide = "legend")

  