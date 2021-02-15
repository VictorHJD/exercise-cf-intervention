##Cystic-Fibrosis Microbiome Project (Mainz)
##Metadata preparation (Try to organize the data provided by the team in Mainz)
##Víctor Hugo Jarquín-Díaz 15.02.2021

library(vroom)
library(tidyverse)
library(readxl)

##Load the data
tech<- read_tsv("~/CF_project/Metadata/technical_Metadata.tsv") ##Contains technical data from DNA (Mainz)
lung<- read_excel("~/CF_project/Metadata/LungFunction_PhysicalFitness.xlsx") ##Contains lung function data
nutri<- read_excel("~/CF_project/Metadata/Nutrition_BodyComposition.xlsx") ##Contains nutritional data
clinic<- read_excel("~/CF_project/Metadata/KlinDaten080221.xlsx") ##Contains clinical data (Antibiotics, Bacterial isolates, etc)
resp<- read_excel("~/CF_project/Metadata/Responder_nonResponder.xlsx") ##Response to intervantions data (ask how were assigned the categories)


##Select useful data and uniform columns
##1) Tech data
tech%>%
  select(Comed_token,X.SampleID,Position, 
         material, extract_quant_ng_ul, total_ng_DNA, dna_quant_ng_ul, Benzoase, Comed_id, Comed_visit)%>%
  rename(SampleID= 2, Patient_number= 9, Visit= 10)%>%
  mutate(Patient_number = paste0("P", Patient_number))%>%
  mutate(Visit = paste0("V", Visit))%>%
  group_by(Patient_number)-> tech
  
##2)Lung function data
##Re-shape dataframe
lung%>%
  rename(Patient_number= 1, sex=2)%>%
  select(Patient_number, sex, age, ppFEV1_V1, ppFEV1_V2, ppFEV1_V3)%>%
  gather(ppFEV1_Visit, ppFEV1, ppFEV1_V1:ppFEV1_V3)%>%
  separate(ppFEV1_Visit, c("Test", "Visit"))%>%
  mutate(Comed_token= paste0(Patient_number, Visit))%>%
  select(Patient_number, sex, age, Visit, Comed_token, ppFEV1)-> tmp1

lung%>%
  rename(Patient_number= 1)%>%
  select(Patient_number, DistanzV1:DistanzV3)%>%
  gather(Dist_Visit, Dist, DistanzV1:DistanzV3)%>%
  mutate(Visit = case_when(Dist_Visit == "DistanzV1"  ~ "V1",
                           Dist_Visit == "DistanzV2" ~ "V2",
                           Dist_Visit == "DistanzV3" ~ "V3"))%>%
  mutate(Comed_token= paste0(Patient_number, Visit))%>%
  select(Comed_token, Dist)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")

lung%>%
  rename(Patient_number= 1)%>%
  select(Patient_number, `Peakpower (W/kg)V1`,`Peakpower (W/kg)V2`, `Peakpower (W/kg)V3`)%>%
  gather(Peak_power_Visit, Peak_power, `Peakpower (W/kg)V1`:`Peakpower (W/kg)V3`)%>%
  mutate(Visit = case_when(Peak_power_Visit == "Peakpower (W/kg)V1"  ~ "V1",
                           Peak_power_Visit == "Peakpower (W/kg)V2" ~ "V2",
                           Peak_power_Visit == "Peakpower (W/kg)V3" ~ "V3"))%>%
  mutate(Comed_token= paste0(Patient_number, Visit))%>%
  select(Comed_token, Peak_power)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")

lung%>%
  rename(Patient_number= 1)%>%
  select(Patient_number, `VO2peak (ml/min/kg)V1`:`VO2peak (ml/min/kg)V3`)%>%
  gather(V02_Visit, V02_A, `VO2peak (ml/min/kg)V1`:`VO2peak (ml/min/kg)V3`)%>%
  mutate(Visit = case_when(V02_Visit == "VO2peak (ml/min/kg)V1"  ~ "V1",
                           V02_Visit == "VO2peak (ml/min/kg)V2" ~ "V2",
                           V02_Visit == "VO2peak (ml/min/kg)V3" ~ "V3"))%>%
  mutate(Comed_token= paste0(Patient_number, Visit))%>%
  select(Comed_token, V02_A)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")

lung%>%
  rename(Patient_number= 1)%>%
  select(Patient_number, `VO2peak (l/min)V1`:`VO2peak (l/min)V3`)%>%
  gather(V02_Visit, V02_B, `VO2peak (l/min)V1`:`VO2peak (l/min)V3`)%>%
  mutate(Visit = case_when(V02_Visit == "VO2peak (l/min)V1"  ~ "V1",
                           V02_Visit == "VO2peak (l/min)V2" ~ "V2",
                           V02_Visit == "VO2peak (l/min)V3" ~ "V3"))%>%
  mutate(Comed_token= paste0(Patient_number, Visit))%>%
  select(Comed_token, V02_B)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")

lung%>%
  rename(Patient_number= 1)%>%
  select(Patient_number, BMIV1:BMIV3)%>%
  gather(BMI_Visit, BMI, BMIV1:BMIV3)%>%
  mutate(Visit = case_when(BMI_Visit == "BMIV1"  ~ "V1",
                           BMI_Visit == "BMIV2" ~ "V2",
                           BMI_Visit == "BMIV3" ~ "V3"))%>%
  mutate(Comed_token= paste0(Patient_number, Visit))%>%
  select(Comed_token, BMI)-> tmp2

lung<- left_join(tmp1, tmp2, by= "Comed_token") ##Index table usable in R
rm(tmp1, tmp2)

##3) Nutritional data
##Re-shape dataframe
nutri%>%
  rename(Patient_number= 1)%>%
  mutate(Patient_number = paste0("P", Patient_number))%>%
  select(-c(V1_BMI:V3_BMI))-> nutri

nutri%>%
  select(Patient_number, `FFM_V1_Charatsi(kg)`:`FFM_V3_Charatsi(kg)`)%>%
  gather(FFM_Visit, FFM_Charatsi, `FFM_V1_Charatsi(kg)`:`FFM_V3_Charatsi(kg)`)%>%
  mutate(Visit = case_when(FFM_Visit == "FFM_V1_Charatsi(kg)"  ~ "V1",
                           FFM_Visit == "FFM_V2_Charatsi(kg)" ~ "V2",
                           FFM_Visit == "FFM_V3_Charatsi(kg)" ~ "V3"))%>%
  mutate(Comed_token= paste0(Patient_number, Visit))%>%
  select(Patient_number, Visit, Comed_token, FFM_Charatsi)-> tmp1

