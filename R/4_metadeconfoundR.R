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
  theme(axis.text.x = element_text(size = 10, angle = 90, face="bold", color="black"),
        axis.text.y = element_text(size = 9, face="italic", color="black"))+
  scale_x_discrete(labels=c( "Dist"= "Distance","pFVC_Response" = "pFVC Response", "Polyethylenglykol_Movicol" = "Polyeth-Movicol", 
                             "Peak_power" = "Peak power" , "FFM_Luk" = "FFM (Lukaski %)", "FFM_Charatsi" = "FFM (Charatsi kg)",
                             "electrolyte_supp" = "Electrolyte Supp", "FFM_Response" = "FFM Response", "sex" = "Sex", 
                             "Sport_Response" = "Sport Response", "extract_quant_ng_ul"= "DNA concentration (ng/µL)", "total_ng_DNA"= "Total DNA (ng)",
                             "Nutrition_supplementation"="Nutrition supp.", "heart_med"= "Heart medication", "DNAse_inh"= "DNAse inh"))+
  scale_y_discrete(labels=c( "CAG-56"= "Firmicutes CAG:56","[Ruminococcus]_gnavus_group" = "Ruminococcus gnavus group", 
                             "DTU089" = "Ruminococcus DTU:089", "[Clostridium]_innocuum_group" = "Clostridium innocuum group",
                             "[Eubacterium]_ventriosum_group" = "Eubacterium ventriosum group", 
                             "Erysipelotrichaceae_UCG-003"= "Erysipelotrichaceae UCG:003", 
                             "Lachnospiraceae_UCG-004"= "Lachnospiraceae UCG:004", "CAG-352"=  "Clostridium CAG:352",
                             "Family_XIII_UCG-001"= "Family XIII UCG:001", "Lachnospiraceae_FCS020_group"= "Lachnospiraceae FCS020 group",
                             "UCG-002"= "Ruminococcaceae UCG:002", "Lachnospiraceae_NK4A136_group"= "Lachnospiraceae NK4A136 group",
                             "NK4A214_group"= "Lachnospiraceae NK4A214 group",  "[Eubacterium]_hallii_group" = "Eubacterium hallii group", 
                             "Lachnospiraceae_UCG-009"= "Lachnospiraceae UCG:009", "Lachnospiraceae_UCG-010"= "Lachnospiraceae UCG:010",
                             "UCG-003"="Lachnospiraceae UCG:003"))-> A

Cun.V1V2+
  xlab("Variables")+
  ylab("ASVs Genus-level")+
  labs(title="MetaDeconfoundR summarizing coneiform plot (Stool Control: V1 vs Case: V2)",
       tag = "B)", caption = "Using V1V2 data", 
       fill= "Effect Size \n (Cliff's Delta)", shape= "Confounding status")+
  theme(legend.key = element_rect(color = "black"))+
  theme(axis.text.x = element_text(size = 10, angle = 90, face="bold", color="black"),
        axis.text.y = element_text(size = 9, face="italic", color="black"))+
  scale_x_discrete(labels=c( "Dist"= "Distance","pFVC_Response" = "pFVC Response", "Polyethylenglykol_Movicol" = "Polyeth-Movicol", 
                             "Peak_power" = "Peak power" , "FFM_Luk" = "FFM (Lukaski %)", "FFM_Charatsi" = "FFM (Charatsi kg)",
                   "electrolyte_supp" = "Electrolyte Supp", "FFM_Response" = "FFM Response", "sex" = "Sex", 
                   "Sport_Response" = "Sport Response"))+
  scale_y_discrete(labels=c( "CAG-56"= "Firmicutes CAG:56","[Ruminococcus]_gnavus_group" = "Ruminococcus gnavus group", 
                             "DTU089" = "Ruminococcus DTU:089", "[Clostridium]_innocuum_group" = "Clostridium innocuum group",
                             "[Eubacterium]_ventriosum_group" = "Eubacterium ventriosum group", 
                             "Erysipelotrichaceae_UCG-003"= "Erysipelotrichaceae UCG:003", 
                             "Lachnospiraceae_UCG-004"= "Lachnospiraceae UCG:004", "CAG-352"=  "Clostridium CAG:352",
                             "Family_XIII_UCG-001"= "Family XIII UCG:001"))-> B
Cun.V2V3+
  xlab("Variables")+
  ylab("ASVs Genus-level")+
  labs(title="MetaDeconfoundR summarizing coneiform plot (Stool Control: V2 vs Case: V3)",
       tag = "C)", caption = "Using V2V3 data", 
       fill= "Effect Size \n (Cliff's Delta)", shape= "Confounding status")+
  theme(legend.key = element_rect(color = "black"))+
  theme(axis.text.x = element_text(size = 10, angle = 90, face="bold", color="black"),
        axis.text.y = element_text(size = 9, face="italic", color="black"))+
  scale_x_discrete(labels=c( "Dist"= "Distance","pFVC_Response" = "pFVC Response", "Polyethylenglykol_Movicol" = "Polyeth-Movicol", 
                             "Peak_power" = "Peak power" , "FFM_Luk" = "FFM (Lukaski %)", "FFM_Charatsi" = "FFM (Charatsi kg)",
                             "electrolyte_supp" = "Electrolyte Supp", "FFM_Response" = "FFM Response", "sex" = "Sex", 
                             "Sport_Response" = "Sport Response", "extract_quant_ng_ul"= "DNA concentration (ng/µL)", 
                             "dna_quant_ng_ul"= "DNA sequenced (ng/µL)", "total_ng_DNA"= "Total DNA (ng)","anticholinergic_inh"= "Anticholinergic (inh)", "V02_B"= "Volume Oxygen (B)",
                             "Nutrition_supplementation"="Nutrition supp.", "heart_med"= "Heart medication", "DNAse_inh"= "DNAse (inh)"))+
  scale_y_discrete(labels=c( "CAG-56"= "Firmicutes CAG:56","[Ruminococcus]_gnavus_group" = "Ruminococcus gnavus group", 
                             "DTU089" = "Ruminococcus DTU:089", "[Clostridium]_innocuum_group" = "Clostridium innocuum group",
                             "[Eubacterium]_ventriosum_group" = "Eubacterium ventriosum group", 
                             "Erysipelotrichaceae_UCG-003"= "Erysipelotrichaceae UCG:003", 
                             "Lachnospiraceae_UCG-004"= "Lachnospiraceae UCG:004", "CAG-352"=  "Clostridium CAG:352",
                             "Family_XIII_UCG-001"= "Family XIII UCG:001", "Lachnospiraceae_FCS020_group"= "Lachnospiraceae FCS020 group",
                             "UCG-002"= "Ruminococcaceae UCG:002", "Lachnospiraceae_NK4A136_group"= "Lachnospiraceae NK4A136 group",
                             "NK4A214_group"= "Lachnospiraceae NK4A214 group",  "[Eubacterium]_hallii_group" = "Eubacterium hallii group", 
                             "Lachnospiraceae_UCG-009"= "Lachnospiraceae UCG:009", "Lachnospiraceae_UCG-010"= "Lachnospiraceae UCG:010",
                             "UCG-003"="Lachnospiraceae UCG:003", "Clostridium_sensu_stricto_1"= "Clostridium sensu stricto 1", 
                             "Incertae_Sedis"  ="Ruminococcaceae Incertae Sedis"))-> C
Cun.V1V3+
  xlab("Variables")+
  ylab("ASVs Genus-level")+
  labs(title="MetaDeconfoundR summarizing coneiform plot (Stool Control: V1 vs Case: V3)",
       tag = "D)", caption = "Using V1V3 data", 
       fill= "Effect Size \n (Cliff's Delta)", shape= "Confounding status")+
  theme(legend.key = element_rect(color = "black"))+
  theme(axis.text.x = element_text(size = 10, angle = 90, face="bold", color="black"),
        axis.text.y = element_text(size = 9, face="italic", color="black"))+
  scale_x_discrete(labels=c( "Dist"= "Distance","pFVC_Response" = "pFVC Response", "Polyethylenglykol_Movicol" = "Polyeth-Movicol", 
                             "Peak_power" = "Peak power" , "FFM_Luk" = "FFM (Lukaski %)", "FFM_Charatsi" = "FFM (Charatsi kg)",
                             "electrolyte_supp" = "Electrolyte Supp", "FFM_Response" = "FFM Response", "sex" = "Sex", 
                             "Sport_Response" = "Sport Response", "extract_quant_ng_ul"= "DNA concentration (ng/µL)", 
                             "dna_quant_ng_ul"= "DNA sequenced (ng/µL)", "total_ng_DNA"= "Total DNA (ng)","anticholinergic_inh"= "Anticholinergic (inh)", "V02_B"= "Volume Oxygen (B)",
                             "Nutrition_supplementation"="Nutrition supp.", "heart_med"= "Heart medication", "DNAse_inh"= "DNAse (inh)"))+
  scale_y_discrete(labels=c( "CAG-56"= "Firmicutes CAG:56","[Ruminococcus]_gnavus_group" = "Ruminococcus gnavus group", 
                             "DTU089" = "Ruminococcus DTU:089", "[Clostridium]_innocuum_group" = "Clostridium innocuum group",
                             "[Eubacterium]_ventriosum_group" = "Eubacterium ventriosum group", 
                             "Erysipelotrichaceae_UCG-003"= "Erysipelotrichaceae UCG:003", 
                             "Lachnospiraceae_UCG-004"= "Lachnospiraceae UCG:004", "CAG-352"=  "Clostridium CAG:352",
                             "Family_XIII_UCG-001"= "Family XIII UCG:001", "Lachnospiraceae_FCS020_group"= "Lachnospiraceae FCS020 group",
                             "UCG-002"= "Ruminococcaceae UCG:002", "Lachnospiraceae_NK4A136_group"= "Lachnospiraceae NK4A136 group",
                             "NK4A214_group"= "Lachnospiraceae NK4A214 group",  "[Eubacterium]_hallii_group" = "Eubacterium hallii group", 
                             "Lachnospiraceae_UCG-009"= "Lachnospiraceae UCG:009", "Lachnospiraceae_UCG-010"= "Lachnospiraceae UCG:010",
                             "UCG-003"="Lachnospiraceae UCG:003", "Clostridium_sensu_stricto_1"= "Clostridium sensu stricto 1", 
                             "Incertae_Sedis"  ="Ruminococcaceae Incertae Sedis"))-> D

png("CF_project/exercise-cf-intervention/figures/Q4_Deconfound_Stool_all_days.png", units = 'in', res = 300, width=10, height=8)
A
dev.off()
png("CF_project/exercise-cf-intervention/figures/Q4_Deconfound_Stool_V1V2.png", units = 'in', res = 300, width=10, height=8)
B
dev.off()
png("CF_project/exercise-cf-intervention/figures/Q4_Deconfound_Stool_V2V3.png", units = 'in', res = 300, width=10, height=8)
C
dev.off()
png("CF_project/exercise-cf-intervention/figures/Q4_Deconfound_Stool_V1V3.png", units = 'in', res = 300, width=10, height=8)
D
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

Cun.all+
  xlab("Variables")+
  ylab("ASVs Genus-level")+
  #geom_point(aes(color = cyl, size = , shape = gear))+
  labs(tag = "A)", caption = "Using all data with patient as random factor", 
       fill= "Effect Size \n (Cliff's Delta)", shape= "Confounding status")+
  theme(legend.key = element_rect(color = "black"))+
  theme(axis.text.x = element_text(size = 10, angle = 90, face="bold", color="black"),
        axis.text.y = element_text(size = 9, face="italic", color="black"))+
  scale_x_discrete(labels=c( "Dist"= "Distance","pFVC_Response" = "pFVC Response", "Polyethylenglykol_Movicol" = "Polyeth-Movicol", 
                             "Peak_power" = "Peak power" , "FFM_Luk" = "FFM (Lukaski %)", "FFM_Charatsi" = "FFM (Charatsi kg)",
                             "electrolyte_supp" = "Electrolyte Supp", "FFM_Response" = "FFM Response", "sex" = "Sex", 
                             "Sport_Response" = "Sport Response", "extract_quant_ng_ul"= "DNA concentration (ng/µL)", "total_ng_DNA"= "Total DNA (ng)",
                             "Nutrition_supplementation"="Nutrition supp.", "heart_med"= "Heart medication", "DNAse_inh"= "DNAse inh",
                             "antibiotics_inh"= "Antibiotics (inh)", "anticholinergic_inh"= "Anticholinergic (inh)", "Nutrition_Response"="Nutrition Response",
                             "Macrolides_po" = "Macrolides (oral)", "Steroids_nasal" = "Steroids (nasal)", "Steroids_po"= "Steroids (oral)", "V02_A"= "Volume Oxygen (A)"))+
  scale_y_discrete(labels=c( "CAG-56"= "Firmicutes CAG:56","[Ruminococcus]_gnavus_group" = "Ruminococcus gnavus group", 
                             "DTU089" = "Ruminococcus DTU:089", "[Clostridium]_innocuum_group" = "Clostridium innocuum group",
                             "[Eubacterium]_nodatum_group" = "Eubacterium nodatum group", "Candidatus_Saccharimonas" = "Candidatus Saccharimonas",
                             "Erysipelotrichaceae_UCG-003"= "Erysipelotrichaceae UCG:003", 
                             "Lachnospiraceae_UCG-004"= "Lachnospiraceae UCG:004", "CAG-352"=  "Clostridium CAG:352",
                             "Family_XIII_UCG-001"= "Family XIII UCG:001", "Lachnospiraceae_FCS020_group"= "Lachnospiraceae FCS020 group",
                             "UCG-002"= "Ruminococcaceae UCG:002", "Lachnospiraceae_NK4A136_group"= "Lachnospiraceae NK4A136 group",
                             "NK4A214_group"= "Lachnospiraceae NK4A214 group",  "[Eubacterium]_hallii_group" = "Eubacterium hallii group", 
                             "Lachnospiraceae_UCG-009"= "Lachnospiraceae UCG:009", "Lachnospiraceae_UCG-010"= "Lachnospiraceae UCG:010",
                             "UCG-003"="Lachnospiraceae UCG:003", "TM7x" ="Nanosynbacter lyticus TM7x"))-> A
 
Cun.V1V2+
  xlab("Variables")+
  ylab("ASVs Genus-level")+
  #geom_point(aes(color = cyl, size = , shape = gear))+
  labs(title="MetaDeconfoundR summarizing coneiform plot (Sputum Control: V1 vs Case: V2)",
       tag = "B)", caption = "Using all data with patient as random factor", 
       fill= "Effect Size \n (Cliff's Delta)", shape= "Confounding status")+
  theme(legend.key = element_rect(color = "black"))+
  theme(axis.text.x = element_text(size = 10, angle = 90, face="bold", color="black"),
        axis.text.y = element_text(size = 9, face="italic", color="black"))+
  scale_x_discrete(labels=c( "Dist"= "Distance","pFVC_Response" = "pFVC Response", "Polyethylenglykol_Movicol" = "Polyeth-Movicol", 
                             "Peak_power" = "Peak power" , "FFM_Luk" = "FFM (Lukaski %)", "FFM_Charatsi" = "FFM (Charatsi kg)",
                             "electrolyte_supp" = "Electrolyte Supp", "FFM_Response" = "FFM Response", "sex" = "Sex", 
                             "Sport_Response" = "Sport Response", "extract_quant_ng_ul"= "DNA concentration (ng/µL)", "total_ng_DNA"= "Total DNA (ng)",
                             "Nutrition_supplementation"="Nutrition supp.", "heart_med"= "Heart medication", "DNAse_inh"= "DNAse inh",
                             "antibiotics_inh"= "Antibiotics (inh)", "anticholinergic_inh"= "Anticholinergic (inh)", "Nutrition_Response"="Nutrition Response",
                             "Macrolides_po" = "Macrolides (oral)", "Steroids_nasal" = "Steroids (nasal)", "Steroids_po"= "Steroids (oral)", "V02_A"= "Volume Oxygen (A)"))+
  scale_y_discrete(labels=c( "CAG-56"= "Firmicutes CAG:56","[Ruminococcus]_gnavus_group" = "Ruminococcus gnavus group", 
                             "DTU089" = "Ruminococcus DTU:089", "[Clostridium]_innocuum_group" = "Clostridium innocuum group",
                             "[Eubacterium]_nodatum_group" = "Eubacterium nodatum group", "Candidatus_Saccharimonas" = "Candidatus Saccharimonas",
                             "Erysipelotrichaceae_UCG-003"= "Erysipelotrichaceae UCG:003", 
                             "Lachnospiraceae_UCG-004"= "Lachnospiraceae UCG:004", "CAG-352"=  "Clostridium CAG:352",
                             "Family_XIII_UCG-001"= "Family XIII UCG:001", "Lachnospiraceae_FCS020_group"= "Lachnospiraceae FCS020 group",
                             "UCG-002"= "Ruminococcaceae UCG:002", "Lachnospiraceae_NK4A136_group"= "Lachnospiraceae NK4A136 group",
                             "NK4A214_group"= "Lachnospiraceae NK4A214 group",  "[Eubacterium]_hallii_group" = "Eubacterium hallii group", 
                             "Lachnospiraceae_UCG-009"= "Lachnospiraceae UCG:009", "Lachnospiraceae_UCG-010"= "Lachnospiraceae UCG:010",
                             "UCG-003"="Lachnospiraceae UCG:003", "TM7x" ="Nanosynbacter lyticus TM7x"))-> B
Cun.V2V3+
  xlab("Variables")+
  ylab("ASVs Genus-level")+
  #geom_point(aes(color = cyl, size = , shape = gear))+
  labs(title="MetaDeconfoundR summarizing coneiform plot (Sputum Control: V2 vs Case: V3)",
       tag = "C)", caption = "Using all data with patient as random factor", 
       fill= "Effect Size \n (Cliff's Delta)", shape= "Confounding status")+
  theme(legend.key = element_rect(color = "black"))+
  theme(axis.text.x = element_text(size = 10, angle = 90, face="bold", color="black"),
        axis.text.y = element_text(size = 9, face="italic", color="black"))+
  scale_x_discrete(labels=c( "Dist"= "Distance","pFVC_Response" = "pFVC Response", "Polyethylenglykol_Movicol" = "Polyeth-Movicol", 
                             "Peak_power" = "Peak power" , "FFM_Luk" = "FFM (Lukaski %)", "FFM_Charatsi" = "FFM (Charatsi kg)",
                             "electrolyte_supp" = "Electrolyte Supp", "FFM_Response" = "FFM Response", "sex" = "Sex", 
                             "Sport_Response" = "Sport Response", "extract_quant_ng_ul"= "DNA concentration (ng/µL)", "total_ng_DNA"= "Total DNA (ng)",
                             "Nutrition_supplementation"="Nutrition supp.", "heart_med"= "Heart medication", "DNAse_inh"= "DNAse inh",
                             "antibiotics_inh"= "Antibiotics (inh)", "anticholinergic_inh"= "Anticholinergic (inh)", "Nutrition_Response"="Nutrition Response",
                             "Macrolides_po" = "Macrolides (oral)", "Steroids_nasal" = "Steroids (nasal)", "Steroids_po"= "Steroids (oral)", "V02_A"= "Volume Oxygen (A)"))+
  scale_y_discrete(labels=c( "CAG-56"= "Firmicutes CAG:56","[Ruminococcus]_gnavus_group" = "Ruminococcus gnavus group", 
                             "DTU089" = "Ruminococcus DTU:089", "[Clostridium]_innocuum_group" = "Clostridium innocuum group",
                             "[Eubacterium]_nodatum_group" = "Eubacterium nodatum group", "Candidatus_Saccharimonas" = "Candidatus Saccharimonas",
                             "Erysipelotrichaceae_UCG-003"= "Erysipelotrichaceae UCG:003", 
                             "Lachnospiraceae_UCG-004"= "Lachnospiraceae UCG:004", "CAG-352"=  "Clostridium CAG:352",
                             "Family_XIII_UCG-001"= "Family XIII UCG:001", "Lachnospiraceae_FCS020_group"= "Lachnospiraceae FCS020 group",
                             "UCG-002"= "Ruminococcaceae UCG:002", "Lachnospiraceae_NK4A136_group"= "Lachnospiraceae NK4A136 group",
                             "NK4A214_group"= "Lachnospiraceae NK4A214 group",  "[Eubacterium]_hallii_group" = "Eubacterium hallii group", 
                             "Lachnospiraceae_UCG-009"= "Lachnospiraceae UCG:009", "Lachnospiraceae_UCG-010"= "Lachnospiraceae UCG:010",
                             "UCG-003"="Lachnospiraceae UCG:003", "TM7x" ="Nanosynbacter lyticus TM7x"))-> C  

Cun.V1V3+
  xlab("Variables")+
  ylab("ASVs Genus-level")+
  #geom_point(aes(color = cyl, size = , shape = gear))+
  labs(title="MetaDeconfoundR summarizing coneiform plot (Sputum Control: V1 vs Case: V3)",
       tag = "D)", caption = "Using all data with patient as random factor", 
       fill= "Effect Size \n (Cliff's Delta)", shape= "Confounding status")+
  theme(legend.key = element_rect(color = "black"))+
  theme(axis.text.x = element_text(size = 10, angle = 90, face="bold", color="black"),
        axis.text.y = element_text(size = 9, face="italic", color="black"))+
  scale_x_discrete(labels=c( "Dist"= "Distance","pFVC_Response" = "pFVC Response", "Polyethylenglykol_Movicol" = "Polyeth-Movicol", 
                             "Peak_power" = "Peak power" , "FFM_Luk" = "FFM (Lukaski %)", "FFM_Charatsi" = "FFM (Charatsi kg)",
                             "electrolyte_supp" = "Electrolyte Supp", "FFM_Response" = "FFM Response", "sex" = "Sex", 
                             "Sport_Response" = "Sport Response", "extract_quant_ng_ul"= "DNA concentration (ng/µL)", "total_ng_DNA"= "Total DNA (ng)",
                             "Nutrition_supplementation"="Nutrition supp.", "heart_med"= "Heart medication", "DNAse_inh"= "DNAse inh",
                             "antibiotics_inh"= "Antibiotics (inh)", "anticholinergic_inh"= "Anticholinergic (inh)", "Nutrition_Response"="Nutrition Response",
                             "Macrolides_po" = "Macrolides (oral)", "Steroids_nasal" = "Steroids (nasal)",
                             "Steroids_po"= "Steroids (oral)", "V02_A"= "Volume Oxygen (A)", "kcal_kg_day"= "Diet (kcal/kg/day)"))+
  scale_y_discrete(labels=c( "CAG-56"= "Firmicutes CAG:56","[Ruminococcus]_gnavus_group" = "Ruminococcus gnavus group", 
                             "DTU089" = "Ruminococcus DTU:089", "[Clostridium]_innocuum_group" = "Clostridium innocuum group",
                             "[Eubacterium]_nodatum_group" = "Eubacterium nodatum group", "Candidatus_Saccharimonas" = "Candidatus Saccharimonas",
                             "Erysipelotrichaceae_UCG-003"= "Erysipelotrichaceae UCG:003", 
                             "Lachnospiraceae_UCG-004"= "Lachnospiraceae UCG:004", "CAG-352"=  "Clostridium CAG:352",
                             "Family_XIII_UCG-001"= "Family XIII UCG:001", "Lachnospiraceae_FCS020_group"= "Lachnospiraceae FCS020 group",
                             "UCG-002"= "Ruminococcaceae UCG:002", "Lachnospiraceae_NK4A136_group"= "Lachnospiraceae NK4A136 group",
                             "NK4A214_group"= "Lachnospiraceae NK4A214 group",  "[Eubacterium]_hallii_group" = "Eubacterium hallii group", 
                             "Lachnospiraceae_UCG-009"= "Lachnospiraceae UCG:009", "Lachnospiraceae_UCG-010"= "Lachnospiraceae UCG:010",
                             "UCG-003"="Lachnospiraceae UCG:003", "TM7x" ="Nanosynbacter lyticus TM7x"))-> D

png("CF_project/exercise-cf-intervention/figures/Q4_Deconfound_Sputum_all_days.png", units = 'in', res = 300, width=10, height=8)
A
dev.off()
png("CF_project/exercise-cf-intervention/figures/Q4_Deconfound_Sputum_V1V2.png", units = 'in', res = 300, width=10, height=8)
B
dev.off()
png("CF_project/exercise-cf-intervention/figures/Q4_Deconfound_Sputum_V2V3.png", units = 'in', res = 300, width=10, height=8)
C
dev.off()
png("CF_project/exercise-cf-intervention/figures/Q4_Deconfound_Sputum_V1V3.png", units = 'in', res = 300, width=10, height=8)
D
dev.off()