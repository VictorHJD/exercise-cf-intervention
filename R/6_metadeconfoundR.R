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
  sputum.microbiome<- readRDS("CF_project/exercise-cf-intervention/data/Sput_rare_ASV.rds")
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

##Adjust some variables
x%>%
  dplyr::mutate(sex = case_when(sex == 1  ~ 0,
                                sex == 2 ~ 1))-> x

##Subset V1 vs V2 
x%>%
  dplyr::select(SampleID, Visit)%>%
  left_join(y, by= "SampleID")%>%
  dplyr::filter(Visit%in%c(1,2))-> y.V1V2 ## Microbial ASV counts

rownames(y.V1V2)<- y.V1V2$SampleID
y.V1V2$SampleID<- NULL
y.V1V2$Visit<- NULL

x[x$Visit %in% c(1,2), ]-> x.V1V2 ## Metadata

##Transform to binary and re-structure to have in the wright order
##Visit will change to Status and 0= V1 and 1= V2
x.V1V2%>%
  dplyr::mutate(Visit = case_when(Visit == 1  ~ 0,
                           Visit == 2 ~ 1))%>%
  dplyr::relocate(Visit)%>%
  dplyr::rename(Status= Visit)-> x.V1V2  ###Without technical confunders

x.V1V2$SampleID<- NULL

##Subset V2 vs V3 
x%>%
  dplyr::select(SampleID, Visit)%>%
  left_join(y, by= "SampleID")%>%
  dplyr::filter(Visit%in%c(2,3))-> y.V2V3 ## Microbial ASV counts

rownames(y.V2V3)<- y.V2V3$SampleID
y.V2V3$SampleID<- NULL
y.V2V3$Visit<- NULL

x[x$Visit %in% c(2, 3), ]-> x.V2V3 ## Metadata

##Transform to binary and re-structure to have in the wright order
##Visit will change to Status and 0= V2 and 1= V3
x.V2V3%>%
  dplyr::mutate(Visit = case_when(Visit == 2  ~ 0,
                           Visit == 3 ~ 1))%>%
  dplyr::relocate(Visit)%>%
  dplyr::rename(Status= Visit)-> x.V2V3  ###Without technical confunders

x.V2V3$SampleID<- NULL

##Subset V1 vs V3 
x%>%
  dplyr::select(SampleID, Visit)%>%
  left_join(y, by= "SampleID")%>%
  dplyr::filter(Visit%in%c(1,3))-> y.V1V3 ## Microbial ASV counts

rownames(y.V1V3)<- y.V1V3$SampleID
y.V1V3$SampleID<- NULL
y.V1V3$Visit<- NULL

x[x$Visit %in% c(1,3), ]-> x.V1V3 ## Metadata

##Transform to binary and re-structure to have in the wright order
##Visit will change to Status and 0= V1 and 1= V3
x.V1V3%>%
  dplyr::mutate(Visit = case_when(Visit == 1  ~ 0,
                           Visit == 3 ~ 1))%>%
  dplyr::relocate(Visit)%>%
  dplyr::rename(Status= Visit)-> x.V1V3  ###Without technical confunders

x.V1V3$SampleID<- NULL

##Use severity as case-control variable
x$Mutation_severity<- NULL
x%>%
  dplyr::relocate(sex)%>%
  dplyr::rename(Status= Phenotype_severity)-> x.all

x.all$SampleID<- gsub("V\\d+", "\\1", x.all$SampleID)
x.all$SampleID<- gsub("^10P", "\\1", x.all$SampleID)
x.all$SampleID<- gsub("\\D+", "\\1", x.all$SampleID)
x.all$SampleID
x.all$SampleID<- as.numeric(x.all$SampleID)
x.all%>%
  dplyr::rename(Patient = SampleID)-> x.all ###Without technical confunders

###Add Antibiotic burden information 

antibioticB%>%
  dplyr::filter(material== "Stool")%>%
  dplyr::select(c(SampleID, AntibioticBurden_total, AntibioticBurden_iv))%>%
  dplyr::distinct()%>%
  column_to_rownames("SampleID")-> tmp2

tmp2$SampleID<-NULL

x.all<- cbind(x.all, tmp2)

y.all<- y
y.all$SampleID<- NULL

###Let's run the pipeline!!!

#MD.V1V2<- MetaDeconfound(featureMat = y.V1V2, metaMat = x.V1V2, nnodes = 25, randomVar = list("+ (1 | Patient_number)" , c("Patient_number")))
#MD.V2V3<- MetaDeconfound(featureMat = y.V2V3, metaMat = x.V2V3, nnodes = 25, randomVar = list("+ (1 | Patient_number)" , c("Patient_number")))
#MD.V1V3<- MetaDeconfound(featureMat = y.V1V3, metaMat = x.V1V3, nnodes = 25, randomVar = list("+ (1 | Patient_number)" , c("Patient_number")))
MD.all<- MetaDeconfound(featureMat = y.all, metaMat = x.all, nnodes = 25, randomVar = list("+ (1 | Patient)" , c("Patient")))

##Plot the results 

#Cun.V1V2<- BuildHeatmap(MD.V1V2, cuneiform = T, coloring = 1)
#Cun.V2V3<- BuildHeatmap(MD.V2V3, cuneiform = T, coloring = 1)
#Cun.V1V3<- BuildHeatmap(MD.V1V3, cuneiform = T, coloring = 1)
Cun.all<- BuildHeatmap(MD.all, cuneiform = T, coloring = 1)

Cun.all$data$metaVariable <- factor(Cun.all$data$metaVariable, 
                                    levels= c("sex","Length","Weight", "FFM_Charatsi", 
                                              "ppFEV1",  "Peak_power",  "Dist" ,  "DFr", 
                                              "DNAse_inh","heart_med", "Montelukast", "Number_iv_courses_priorstudy",
                                              "Number_iv_courses_duringstudy",
                                              "Nutrition_supplementation", "Polyethylenglykol_Movicol"))
                                               
##Edit cun plots
Cun.all+
  xlab("Variables")+
  ylab("ASVs Genus-level")+
  labs(title = "MetaDeconfoundR summarizing (Stool Microbiome)",tag = "A)",
       fill= "Effect Size \n (Cliff's Delta)", shape= "Confounding status")+
  theme(legend.key = element_rect(color = "black"))+
  theme(axis.text.x = element_text(size = 12, angle = 90, face="bold", color="black"),
        axis.title.x=element_blank(),
        axis.text.y = element_text(size = 12, face="italic", color="black"))+
  scale_x_discrete(labels=c( "Mutation_severity"= "Mutation severity", "Dist"= "Distance", "Polyethylenglykol_Movicol" = "Polyeth-Movicol", 
                             "Peak_power" = "Peak power" , "FFM_Charatsi" = "FFM (Charatsi kg)", "sex" = "Sex", "DFr" = "Fiber",
                             "Number_iv_courses_priorstudy" = "iv Antibiotic courses (prior)",
                             "Number_iv_courses_duringstudy" = "iv Antibiotic courses (during)",
                             "Nutrition_supplementation"="Nutrition supp.", "heart_med"= "Heart medication", "DNAse_inh"= "DNAse inh"))-> A
A <- A + theme(legend.position="none")

#Cun.V1V2+
#  xlab("Variables")+
#  ylab("ASVs Genus-level")+
#  labs(title="MetaDeconfoundR summarizing coneiform plot (Stool microbiome V1 vs V2)",
#       tag = "B)", caption = "Using V1V2 data", 
#       fill= "Effect Size \n (Cliff's Delta)", shape= "Confounding status")+
#  theme(legend.key = element_rect(color = "black"))+
#  theme(axis.text.x = element_text(size = 10, angle = 90, face="bold", color="black"),
#        axis.text.y = element_text(size = 9, face="italic", color="black"))+
#  scale_x_discrete(labels=c( "Dist"= "Distance","pFVC_Response" = "pFVC Response", "Polyethylenglykol_Movicol" = "Polyeth-Movicol", 
#                             "Peak_power" = "Peak power" , "FFM_Luk" = "FFM (Lukaski %)", "FFM_Charatsi" = "FFM (Charatsi kg)",
#                   "electrolyte_supp" = "Electrolyte Supp", "FFM_Response" = "FFM Response", "sex" = "Sex", 
#                   "Sport_Response" = "Sport Response"))+
#  scale_y_discrete(labels=c( "CAG-56"= "Firmicutes CAG:56","[Ruminococcus]_gnavus_group" = "Ruminococcus gnavus group", 
#                             "DTU089" = "Ruminococcus DTU:089", "[Clostridium]_innocuum_group" = "Clostridium innocuum group",
#                             "[Eubacterium]_ventriosum_group" = "Eubacterium ventriosum group", 
#                             "Erysipelotrichaceae_UCG-003"= "Erysipelotrichaceae UCG:003", 
#                             "Lachnospiraceae_UCG-004"= "Lachnospiraceae UCG:004", "CAG-352"=  "Clostridium CAG:352",
#                             "Family_XIII_UCG-001"= "Family XIII UCG:001"))-> B
#Cun.V2V3+
#  xlab("Variables")+
#  ylab("ASVs Genus-level")+
#  labs(title="MetaDeconfoundR summarizing coneiform plot (Stool microbiome V2 vs V3)",
#       tag = "C)", caption = "Using V2V3 data", 
#       fill= "Effect Size \n (Cliff's Delta)", shape= "Confounding status")+
#  theme(legend.key = element_rect(color = "black"))+
#  theme(axis.text.x = element_text(size = 10, angle = 90, face="bold", color="black"),
#        axis.text.y = element_text(size = 9, face="italic", color="black"))+
#  scale_x_discrete(labels=c( "Dist"= "Distance","pFVC_Response" = "pFVC Response", "Polyethylenglykol_Movicol" = "Polyeth-Movicol", 
#                             "Peak_power" = "Peak power" , "FFM_Luk" = "FFM (Lukaski %)", "FFM_Charatsi" = "FFM (Charatsi kg)",
#                             "electrolyte_supp" = "Electrolyte Supp", "FFM_Response" = "FFM Response", "sex" = "Sex", "DFr" = "Fiber",
#                             "Sport_Response" = "Sport Response", "extract_quant_ng_ul"= "DNA concentration (ng/µL)", 
#                             "dna_quant_ng_ul"= "DNA sequenced (ng/µL)", "total_ng_DNA"= "Total DNA (ng)","anticholinergic_inh"= "Anticholinergic (inh)", "V02_B"= "Volume Oxygen (B)",
#                             "Nutrition_supplementation"="Nutrition supp.", "heart_med"= "Heart medication", "DNAse_inh"= "DNAse (inh)"))+
#  scale_y_discrete(labels=c( "CAG-56"= "Firmicutes CAG:56","[Ruminococcus]_gnavus_group" = "Ruminococcus gnavus group", 
#                             "DTU089" = "Ruminococcus DTU:089", "[Clostridium]_innocuum_group" = "Clostridium innocuum group",
#                             "[Eubacterium]_ventriosum_group" = "Eubacterium ventriosum group", 
#                             "Erysipelotrichaceae_UCG-003"= "Erysipelotrichaceae UCG:003", 
#                             "Lachnospiraceae_UCG-004"= "Lachnospiraceae UCG:004", "CAG-352"=  "Clostridium CAG:352",
#                             "Family_XIII_UCG-001"= "Family XIII UCG:001", "Lachnospiraceae_FCS020_group"= "Lachnospiraceae FCS020 group",
#                             "UCG-002"= "Ruminococcaceae UCG:002", "Lachnospiraceae_NK4A136_group"= "Lachnospiraceae NK4A136 group",
#                             "NK4A214_group"= "Lachnospiraceae NK4A214 group",  "[Eubacterium]_hallii_group" = "Eubacterium hallii group", 
#                             "Lachnospiraceae_UCG-009"= "Lachnospiraceae UCG:009", "Lachnospiraceae_UCG-010"= "Lachnospiraceae UCG:010",
#                             "UCG-003"="Lachnospiraceae UCG:003", "Clostridium_sensu_stricto_1"= "Clostridium sensu stricto 1", 
#                             "Incertae_Sedis"  ="Ruminococcaceae Incertae Sedis"))-> C
#Cun.V1V3+
#  xlab("Variables")+
#  ylab("ASVs Genus-level")+
#  labs(title="MetaDeconfoundR summarizing coneiform plot (Stool microbiome V1 vs V3)",
#       tag = "D)", caption = "Using V1V3 data", 
#       fill= "Effect Size \n (Cliff's Delta)", shape= "Confounding status")+
#  theme(legend.key = element_rect(color = "black"))+
#  theme(axis.text.x = element_text(size = 10, angle = 90, face="bold", color="black"),
#        axis.text.y = element_text(size = 9, face="italic", color="black"))+
#  scale_x_discrete(labels=c( "Dist"= "Distance","pFVC_Response" = "pFVC Response", "Polyethylenglykol_Movicol" = "Polyeth-Movicol", 
#                             "Peak_power" = "Peak power" , "FFM_Luk" = "FFM (Lukaski %)", "FFM_Charatsi" = "FFM (Charatsi kg)",
#                             "electrolyte_supp" = "Electrolyte Supp", "FFM_Response" = "FFM Response", "sex" = "Sex", 
#                             "Sport_Response" = "Sport Response", "extract_quant_ng_ul"= "DNA concentration (ng/µL)", 
#                             "dna_quant_ng_ul"= "DNA sequenced (ng/µL)", "total_ng_DNA"= "Total DNA (ng)","anticholinergic_inh"= "Anticholinergic (inh)", "V02_B"= "Volume Oxygen (B)",
#                             "Nutrition_supplementation"="Nutrition supp.", "heart_med"= "Heart medication", "DNAse_inh"= "DNAse (inh)"))+
#  scale_y_discrete(labels=c( "CAG-56"= "Firmicutes CAG:56","[Ruminococcus]_gnavus_group" = "Ruminococcus gnavus group", 
#                             "DTU089" = "Ruminococcus DTU:089", "[Clostridium]_innocuum_group" = "Clostridium innocuum group",
#                             "[Eubacterium]_ventriosum_group" = "Eubacterium ventriosum group", 
#                             "Erysipelotrichaceae_UCG-003"= "Erysipelotrichaceae UCG:003", 
#                             "Lachnospiraceae_UCG-004"= "Lachnospiraceae UCG:004", "CAG-352"=  "Clostridium CAG:352",
#                             "Family_XIII_UCG-001"= "Family XIII UCG:001", "Lachnospiraceae_FCS020_group"= "Lachnospiraceae FCS020 group",
#                             "UCG-002"= "Ruminococcaceae UCG:002", "Lachnospiraceae_NK4A136_group"= "Lachnospiraceae NK4A136 group",
#                             "NK4A214_group"= "Lachnospiraceae NK4A214 group",  "[Eubacterium]_hallii_group" = "Eubacterium hallii group", 
#                             "Lachnospiraceae_UCG-009"= "Lachnospiraceae UCG:009", "Lachnospiraceae_UCG-010"= "Lachnospiraceae UCG:010",
#                             "UCG-003"="Lachnospiraceae UCG:003", "Clostridium_sensu_stricto_1"= "Clostridium sensu stricto 1", 
#                             "Incertae_Sedis"  ="Ruminococcaceae Incertae Sedis"))-> D

#png("CF_project/exercise-cf-intervention/figures/Q4_Deconfound_Stool_all_days.png", units = 'in', res = 300, width=10, height=8)
#A
#dev.off()
#png("CF_project/exercise-cf-intervention/figures/Q4_Deconfound_Stool_V1V2.png", units = 'in', res = 300, width=10, height=8)
#B
#dev.off()
#png("CF_project/exercise-cf-intervention/figures/Q4_Deconfound_Stool_V2V3.png", units = 'in', res = 300, width=10, height=8)
#C
#dev.off()
#png("CF_project/exercise-cf-intervention/figures/Q4_Deconfound_Stool_V1V3.png", units = 'in', res = 300, width=10, height=8)
#D
#dev.off()
#######################2) Sputum####################
##Take SampleID in the right order
SampleID<- rownames(sputum.metadata)

##Transform metadata to dataframe 
x<- as.data.frame(sputum.metadata)
y<- as.data.frame(sputum.microbiome)

##Adjust some variables
x%>%
  dplyr::mutate(sex = case_when(sex == 1  ~ 0,
                         sex == 2 ~ 1))-> x
##Add IDs as column
x$SampleID<- SampleID
y$SampleID<- SampleID

##Subset V1 vs V2 
x%>%
  dplyr::select(SampleID, Visit)%>%
  left_join(y, by= "SampleID")%>%
  dplyr::filter(Visit%in%c(1,2))-> y.V1V2 ## Microbial ASV counts

rownames(y.V1V2)<- y.V1V2$SampleID
y.V1V2$SampleID<- NULL
y.V1V2$Visit<- NULL

x[x$Visit %in% c(1,2), ]-> x.V1V2 ## Metadata

##Transform to binary and re-structure to have in the wright order
##Visit will change to Status and 0= V1 and 1= V2
x.V1V2%>%
  dplyr::mutate(Visit = case_when(Visit == 1  ~ 0,
                           Visit == 2 ~ 1))%>%
  dplyr::relocate(Visit)%>%
  dplyr::rename(Status= Visit)-> x.V1V2  ###Without technical confunders

x.V1V2$SampleID<- NULL

##Subset V2 vs V3 
x%>%
  dplyr::select(SampleID, Visit)%>%
  left_join(y, by= "SampleID")%>%
  dplyr::filter(Visit%in%c(2,3))-> y.V2V3 ## Microbial ASV counts

rownames(y.V2V3)<- y.V2V3$SampleID
y.V2V3$SampleID<- NULL
y.V2V3$Visit<- NULL

x[x$Visit %in% c(2, 3), ]-> x.V2V3 ## Metadata

##Transform to binary and re-structure to have in the wright order
##Visit will change to Status and 0= V2 and 1= V3
x.V2V3%>%
  dplyr::mutate(Visit = case_when(Visit == 2  ~ 0,
                           Visit == 3 ~ 1))%>%
  dplyr::relocate(Visit)%>%
  dplyr::rename(Status= Visit)-> x.V2V3  ###Without technical confunders

x.V2V3$SampleID<- NULL

##Subset V1 vs V3 
x%>%
  dplyr::select(SampleID, Visit)%>%
  left_join(y, by= "SampleID")%>%
  dplyr::filter(Visit%in%c(1,3))-> y.V1V3 ## Microbial ASV counts

rownames(y.V1V3)<- y.V1V3$SampleID
y.V1V3$SampleID<- NULL
y.V1V3$Visit<- NULL

x[x$Visit %in% c(1,3), ]-> x.V1V3 ## Metadata

##Transform to binary and re-structure to have in the wright order
##Visit will change to Status and 0= V1 and 1= V3
x.V1V3%>%
  dplyr::mutate(Visit = case_when(Visit == 1  ~ 0,
                           Visit == 3 ~ 1))%>%
  dplyr::relocate(Visit)%>%
  dplyr::rename(Status= Visit)-> x.V1V3  ###Without technical confunders

x.V1V3$SampleID<- NULL

##Use severity as case-control variable
x$Mutation_severity<- NULL

x%>%
  dplyr::relocate(Phenotype_severity)%>%
  dplyr::rename(Status= Phenotype_severity)-> x.all

x.all$SampleID<- gsub("V\\d+", "\\1", x.all$SampleID)
x.all$SampleID<- gsub("^10P", "\\1", x.all$SampleID)
x.all$SampleID<- gsub("\\D+", "\\1", x.all$SampleID)
x.all$SampleID
x.all$SampleID<- as.numeric(x.all$SampleID)
x.all%>%
  dplyr::rename(Patient = SampleID)-> x.all ###Without technical confunders

###Add Antibiotic burden information 

antibioticB%>%
  dplyr::filter(material== "Sputum")%>%
  dplyr::select(c(SampleID, AntibioticBurden_total, AntibioticBurden_iv))%>%
  dplyr::distinct()%>%
  column_to_rownames("SampleID")-> tmp2

tmp2$SampleID<-NULL

x.all<- cbind(x.all, tmp2)

y.all<- y
y.all$SampleID<- NULL


###Let's run the pipeline!!!

#MD.V1V2<- MetaDeconfound(featureMat = y.V1V2, metaMat = x.V1V2, nnodes = 25, randomVar = list("+ (1 | Patient_number)" , c("Patient_number")))
#MD.V2V3<- MetaDeconfound(featureMat = y.V2V3, metaMat = x.V2V3, nnodes = 25, randomVar = list("+ (1 | Patient_number)" , c("Patient_number")))
#MD.V1V3<- MetaDeconfound(featureMat = y.V1V3, metaMat = x.V1V3, nnodes = 25, randomVar = list("+ (1 | Patient_number)" , c("Patient_number")))
#MD.V1V3<- MetaDeconfound(featureMat = y.V1V3, metaMat = x.V1V3, nnodes = 25)
MD.all<- MetaDeconfound(featureMat = y.all, metaMat = x.all, nnodes = 25, randomVar = list("+ (1 | Patient)" , c("Patient")))

##Plot the results 

#Cun.V1V2<- BuildHeatmap(MD.V1V2, cuneiform = T, coloring = 1)
#Cun.V2V3<- BuildHeatmap(MD.V2V3, cuneiform = T, coloring = 1)
#Cun.V1V3<- BuildHeatmap(MD.V1V3, cuneiform = T, coloring = 1)
Cun.all<- BuildHeatmap(MD.all, cuneiform = T, coloring = 1)

Cun.all$data$metaVariable <- factor(Cun.all$data$metaVariable, 
                                    levels= c("sex", "Length", "Weight", "FFM_Charatsi",
                                              "ppFVC","Peak_power", "CHO", "Lipids",
                                              "Macrolides_po", "Montelukast", "Nutrition_supplementation", 
                                              "PPI", "Steroids_inh", "Pseudomonas_status", 
                                              "AntibioticBurden_total", "AntibioticBurden_iv"))

  
##Edit cun plots
Cun.all+
  xlab("Variables")+
  ylab("ASVs Genus-level")+
  labs(title = "MetaDeconfoundR summarizing (Sputum Microbiome)",
       fill= "Effect Size \n (Cliff's Delta)", shape= "Confounding status")+
  theme(legend.key = element_rect(color = "black"))+
  theme(axis.text.x = element_text(size = 12, angle = 90, face="bold", color="black"),
        axis.title.x=element_blank(),
        axis.text.y = element_text(size = 12, face="italic", color="black"))+
  scale_x_discrete(labels=c("Peak_power" = "Peak power" , "FFM_Charatsi" = "FFM (Charatsi kg)",
                              "sex" = "Sex", "Pseudomonas_status"= "Pseudomonas culture",
                            "AntibioticBurden_total" = "Total Antibiotic burden",
                            "AntibioticBurden_iv" = "iv Antibiotic burden",
                             "Nutrition_supplementation"="Nutrition supp.", "Steroids_inh"= "Steroids (inh)",
                             "Macrolides_po" = "Macrolides (oral)"))-> B

plot<-ggarrange(A,  B,  ncol=1, nrow=2, common.legend = T, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q6_Deconfound_Stool_Sputum.png", plot = plot, width = 15, height = 25, dpi = 600)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q6_Deconfound_Stool_Sputum.pdf", plot = plot, width = 15, height = 25, dpi = 600)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q6_Deconfound_Stool_Sputum.svg", plot = plot, width = 15, height = 25, dpi = 600)

##Store them separately for crafting
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q6_Deconfound_Stool.png", plot = A, width = 10, height = 10)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q6_Deconfound_Stool.pdf", plot = A, width = 10, height = 10)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q6_Deconfound_Stool.svg", plot = A, width = 10, height = 10)


ggsave(file = "CF_project/exercise-cf-intervention/figures/Q6_Deconfound_Sputum.png", plot = B, width = 10, height = 10)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q6_Deconfound_Sputum.pdf", plot = B, width = 10, height = 10)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q6_Deconfound_Sputum.svg", plot = B, width = 10, height = 10)

#Cun.V1V2+
#  xlab("Variables")+
#  ylab("ASVs Genus-level")+
  #geom_point(aes(color = cyl, size = , shape = gear))+
#  labs(title="MetaDeconfoundR summarizing coneiform plot (Sputum V1 vs V2)",
#       tag = "B)", caption = "Using V1-V2 data", 
#       fill= "Effect Size \n (Cliff's Delta)", shape= "Confounding status")+
#  theme(legend.key = element_rect(color = "black"))+
#  theme(axis.text.x = element_text(size = 10, angle = 90, face="bold", color="black"),
#        axis.text.y = element_text(size = 9, face="italic", color="black"))+
#  scale_x_discrete(labels=c( "Dist"= "Distance","pFVC_Response" = "pFVC Response", "Polyethylenglykol_Movicol" = "Polyeth-Movicol", 
#                             "Peak_power" = "Peak power" , "FFM_Luk" = "FFM (Lukaski %)", "FFM_Charatsi" = "FFM (Charatsi kg)",
#                             "electrolyte_supp" = "Electrolyte Supp", "FFM_Response" = "FFM Response", "sex" = "Sex", 
#                             "Sport_Response" = "Sport Response", "extract_quant_ng_ul"= "DNA concentration (ng/µL)", "total_ng_DNA"= "Total DNA (ng)",
#                             "Nutrition_supplementation"="Nutrition supp.", "heart_med"= "Heart medication", "DNAse_inh"= "DNAse inh",
#                             "antibiotics_inh"= "Antibiotics (inh)", "anticholinergic_inh"= "Anticholinergic (inh)", "Nutrition_Response"="Nutrition Response",
#                             "Macrolides_po" = "Macrolides (oral)", "Steroids_nasal" = "Steroids (nasal)", "Steroids_po"= "Steroids (oral)", "V02_A"= "Volume Oxygen (A)"))+
#  scale_y_discrete(labels=c( "CAG-56"= "Firmicutes CAG:56","[Ruminococcus]_gnavus_group" = "Ruminococcus gnavus group", 
#                             "DTU089" = "Ruminococcus DTU:089", "[Clostridium]_innocuum_group" = "Clostridium innocuum group",
#                             "[Eubacterium]_nodatum_group" = "Eubacterium nodatum group", "Candidatus_Saccharimonas" = "Candidatus Saccharimonas",
#                             "Erysipelotrichaceae_UCG-003"= "Erysipelotrichaceae UCG:003", 
#                             "Lachnospiraceae_UCG-004"= "Lachnospiraceae UCG:004", "CAG-352"=  "Clostridium CAG:352",
#                             "Family_XIII_UCG-001"= "Family XIII UCG:001", "Lachnospiraceae_FCS020_group"= "Lachnospiraceae FCS020 group",
#                             "UCG-002"= "Ruminococcaceae UCG:002", "Lachnospiraceae_NK4A136_group"= "Lachnospiraceae NK4A136 group",
#                             "NK4A214_group"= "Lachnospiraceae NK4A214 group",  "[Eubacterium]_hallii_group" = "Eubacterium hallii group", 
#                             "Lachnospiraceae_UCG-009"= "Lachnospiraceae UCG:009", "Lachnospiraceae_UCG-010"= "Lachnospiraceae UCG:010",
#                             "UCG-003"="Lachnospiraceae UCG:003", "TM7x" ="Nanosynbacter lyticus TM7x"))-> B
#Cun.V2V3+
#  xlab("Variables")+
#  ylab("ASVs Genus-level")+
#  #geom_point(aes(color = cyl, size = , shape = gear))+
#  labs(title="MetaDeconfoundR summarizing coneiform plot (Sputum V2 vs V3)",
#       tag = "C)", caption = "Using V2-V3 data", 
#       fill= "Effect Size \n (Cliff's Delta)", shape= "Confounding status")+
#  theme(legend.key = element_rect(color = "black"))+
#  theme(axis.text.x = element_text(size = 10, angle = 90, face="bold", color="black"),
#        axis.text.y = element_text(size = 9, face="italic", color="black"))+
#  scale_x_discrete(labels=c( "Dist"= "Distance","pFVC_Response" = "pFVC Response", "Polyethylenglykol_Movicol" = "Polyeth-Movicol", 
#                             "Peak_power" = "Peak power" , "FFM_Luk" = "FFM (Lukaski %)", "FFM_Charatsi" = "FFM (Charatsi kg)",
#                             "electrolyte_supp" = "Electrolyte Supp", "FFM_Response" = "FFM Response", "sex" = "Sex", "Pseudomonas_status"= "Psudomonas culture",
#                             "Sport_Response" = "Sport Response", "extract_quant_ng_ul"= "DNA concentration (ng/µL)", "total_ng_DNA"= "Total DNA (ng)",
#                             "Nutrition_supplementation"="Nutrition supp.", "heart_med"= "Heart medication", "DNAse_inh"= "DNAse inh",
#                             "antibiotics_inh"= "Antibiotics (inh)", "anticholinergic_inh"= "Anticholinergic (inh)", "Nutrition_Response"="Nutrition Response",
#                             "Macrolides_po" = "Macrolides (oral)", "Steroids_nasal" = "Steroids (nasal)", "Steroids_po"= "Steroids (oral)", "V02_A"= "Volume Oxygen (A)"))+
#  scale_y_discrete(labels=c( "CAG-56"= "Firmicutes CAG:56","[Ruminococcus]_gnavus_group" = "Ruminococcus gnavus group", 
#                             "DTU089" = "Ruminococcus DTU:089", "[Clostridium]_innocuum_group" = "Clostridium innocuum group",
#                             "[Eubacterium]_nodatum_group" = "Eubacterium nodatum group", "Candidatus_Saccharimonas" = "Candidatus Saccharimonas",
#                             "Erysipelotrichaceae_UCG-003"= "Erysipelotrichaceae UCG:003", 
#                             "Lachnospiraceae_UCG-004"= "Lachnospiraceae UCG:004", "CAG-352"=  "Clostridium CAG:352",
#                             "Family_XIII_UCG-001"= "Family XIII UCG:001", "Lachnospiraceae_FCS020_group"= "Lachnospiraceae FCS020 group",
#                             "UCG-002"= "Ruminococcaceae UCG:002", "Lachnospiraceae_NK4A136_group"= "Lachnospiraceae NK4A136 group",
#                             "NK4A214_group"= "Lachnospiraceae NK4A214 group",  "[Eubacterium]_hallii_group" = "Eubacterium hallii group", 
#                             "Lachnospiraceae_UCG-009"= "Lachnospiraceae UCG:009", "Lachnospiraceae_UCG-010"= "Lachnospiraceae UCG:010",
#                             "UCG-003"="Lachnospiraceae UCG:003", "TM7x" ="Nanosynbacter lyticus TM7x"))-> C  #

#Cun.V1V3+
#  xlab("Variables")+
#  ylab("ASVs Genus-level")+
#  #geom_point(aes(color = cyl, size = , shape = gear))+
#  labs(title="MetaDeconfoundR summarizing coneiform plot (Sputum V1 vs V3)",
#       tag = "D)", caption = "Using V1-V3 data", 
#       fill= "Effect Size \n (Cliff's Delta)", shape= "Confounding status")+
#  theme(legend.key = element_rect(color = "black"))+
#  theme(axis.text.x = element_text(size = 10, angle = 90, face="bold", color="black"),
#        axis.text.y = element_text(size = 9, face="italic", color="black"))+
#  scale_x_discrete(labels=c( "Dist"= "Distance","pFVC_Response" = "pFVC Response", "Polyethylenglykol_Movicol" = "Polyeth-Movicol", 
#                             "Peak_power" = "Peak power" , "FFM_Luk" = "FFM (Lukaski %)", "FFM_Charatsi" = "FFM (Charatsi kg)",
#                             "electrolyte_supp" = "Electrolyte Supp", "FFM_Response" = "FFM Response", "sex" = "Sex", 
#                             "Sport_Response" = "Sport Response", "extract_quant_ng_ul"= "DNA concentration (ng/µL)", "total_ng_DNA"= "Total DNA (ng)",
#                             "Nutrition_supplementation"="Nutrition supp.", "heart_med"= "Heart medication", "DNAse_inh"= "DNAse inh",
#                             "antibiotics_inh"= "Antibiotics (inh)", "anticholinergic_inh"= "Anticholinergic (inh)", "Nutrition_Response"="Nutrition Response",
#                             "Macrolides_po" = "Macrolides (oral)", "Steroids_nasal" = "Steroids (nasal)", "Pseudomonas_status"= "Pseudomonas culture",
#                             "Steroids_po"= "Steroids (oral)", "V02_A"= "Volume Oxygen (A)", "kcal_kg_day"= "Diet (kcal/kg/day)"))+
#  scale_y_discrete(labels=c( "CAG-56"= "Firmicutes CAG:56","[Ruminococcus]_gnavus_group" = "Ruminococcus gnavus group", 
#                             "DTU089" = "Ruminococcus DTU:089", "[Clostridium]_innocuum_group" = "Clostridium innocuum group",
#                             "[Eubacterium]_nodatum_group" = "Eubacterium nodatum group", "Candidatus_Saccharimonas" = "Candidatus Saccharimonas",
#                            "Erysipelotrichaceae_UCG-003"= "Erysipelotrichaceae UCG:003", 
#                             "Lachnospiraceae_UCG-004"= "Lachnospiraceae UCG:004", "CAG-352"=  "Clostridium CAG:352",
#                             "Family_XIII_UCG-001"= "Family XIII UCG:001", "Lachnospiraceae_FCS020_group"= "Lachnospiraceae FCS020 group",
#                             "UCG-002"= "Ruminococcaceae UCG:002", "Lachnospiraceae_NK4A136_group"= "Lachnospiraceae NK4A136 group",
#                             "NK4A214_group"= "Lachnospiraceae NK4A214 group",  "[Eubacterium]_hallii_group" = "Eubacterium hallii group", 
#                             "Lachnospiraceae_UCG-009"= "Lachnospiraceae UCG:009", "Lachnospiraceae_UCG-010"= "Lachnospiraceae UCG:010",
#                             "UCG-003"="Lachnospiraceae UCG:003", "TM7x" ="Nanosynbacter lyticus TM7x"))-> D#

#png("CF_project/exercise-cf-intervention/figures/Q4_Deconfound_Sputum_all_days.png", units = 'in', res = 300, width=10, height=8)
#A
#dev.off()
#png("CF_project/exercise-cf-intervention/figures/Q4_Deconfound_Sputum_V1V2.png", units = 'in', res = 300, width=10, height=8)
#B
#dev.off()
#png("CF_project/exercise-cf-intervention/figures/Q4_Deconfound_Sputum_V2V3.png", units = 'in', res = 300, width=10, height=8)
#C
#dev.off()
#png("CF_project/exercise-cf-intervention/figures/Q4_Deconfound_Sputum_V1V3.png", units = 'in', res = 300, width=10, height=8)
#D
#dev.off()


