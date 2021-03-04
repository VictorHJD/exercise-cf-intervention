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
clinic<- read.csv("~/CF_project/Metadata/sample_data_indexed_medication.csv", sep = ";") ##Contains clinical data (Antibiotics, Bacterial isolates, etc)
genetics<- read_excel("~/CF_project/Metadata/KlinDaten080221.xlsx")
resp<- read_excel("~/CF_project/Metadata/Responder_nonResponder.xlsx") ##Response to interventions data (ask how were assigned the categories)

##Pre-organized data with adjustments to samples:
##Sample P3V3 stool was incorrectly assigned, it is actually P2V3 stool
##Sample P17V2 stool, it is P17V3 stool

data.mainz<- read_tsv("~/CF_project/Metadata/Sample_Metadata_combine_rebecca.csv")

##Select useful data and uniform columns
data.mainz%>%
  select(c(1,2,5,7,11,13,16,17,27, 55:75))%>%
  rename(SampleID= 2, Visit= 10, sex= 11)%>%
  mutate(Visit = paste0("V", Visit))%>%
  mutate(Patient_number= Comed_token)%>%
  mutate(Patient_number = gsub("V\\d+", "\\1", basename(Patient_number)))%>%
  group_by(Patient_number)->data.mainz


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
##Eliminate BMI that was already included in the previous table (this value actually can be computed based on wight and lenght)
nutri%>%
  rename(Patient_number= 1)%>%
  mutate(Patient_number = paste0("P", Patient_number))%>%
  select(-c(V1_BMI:V3_BMI))-> nutri
##Fat free mass (Charatsi)
nutri%>%
  select(Patient_number, `FFM_V1_Charatsi(kg)`:`FFM_V3_Charatsi(kg)`)%>%
  gather(FFM_Visit, FFM_Charatsi, `FFM_V1_Charatsi(kg)`:`FFM_V3_Charatsi(kg)`)%>%
  mutate(Visit = case_when(FFM_Visit == "FFM_V1_Charatsi(kg)"  ~ "V1",
                           FFM_Visit == "FFM_V2_Charatsi(kg)" ~ "V2",
                           FFM_Visit == "FFM_V3_Charatsi(kg)" ~ "V3"))%>%
  mutate(Comed_token= paste0(Patient_number, Visit))%>%
  select(Patient_number, Visit, Comed_token, FFM_Charatsi)-> tmp1
##Lenght (cm)
nutri%>%
  select(Patient_number, V1_length_cm:V3_length_cm)%>%
  gather(Length_Visit, Length, V1_length_cm:V3_length_cm)%>%
  mutate(Visit = case_when(Length_Visit == "V1_length_cm"  ~ "V1",
                           Length_Visit == "V2_length_cm" ~ "V2",
                           Length_Visit == "V3_length_cm" ~ "V3"))%>%
  mutate(Comed_token= paste0(Patient_number, Visit))%>%
  select(Comed_token, Length)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")
##Weight (kg)
nutri%>%
  select(Patient_number, V1_weight_kg:V3_weight_kg)%>%
  gather(Weight_Visit, Weight, V1_weight_kg:V3_weight_kg)%>%
  mutate(Visit = case_when(Weight_Visit == "V1_weight_kg"  ~ "V1",
                           Weight_Visit == "V2_weight_kg" ~ "V2",
                           Weight_Visit == "V3_weight_kg" ~ "V3"))%>%
  mutate(Comed_token= paste0(Patient_number, Visit))%>%
  select(Comed_token, Weight)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")
##Fat free mass (Lukaski) bioelectrical impedance measurements
nutri%>%
  select(Patient_number, `V1_FFM_%(Lukaski)`:`V3_FFM_%(Lukaski)`)%>%
  gather(FFM_Luk_Visit, FFM_Luk,`V1_FFM_%(Lukaski)`:`V3_FFM_%(Lukaski)`)%>%
  mutate(Visit = case_when(FFM_Luk_Visit == "V1_FFM_%(Lukaski)"  ~ "V1",
                           FFM_Luk_Visit == "V2_FFM_%(Lukaski)" ~ "V2",
                           FFM_Luk_Visit == "V3_FFM_%(Lukaski)" ~ "V3"))%>%
  mutate(Comed_token= paste0(Patient_number, Visit))%>%
  select(Comed_token, FFM_Luk)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")
##Average of kcal
nutri%>%
  select(Patient_number, `V2 kcal Durchschnitt`:`V3 kcal Durchschnitt`)%>%
  gather(kcal_Visit, kcal_Avg,`V2 kcal Durchschnitt`:`V3 kcal Durchschnitt`)%>%
  mutate(Visit = case_when(kcal_Visit == "V1 kcal Durchschnitt"  ~ "V1",
                           kcal_Visit == "V2 kcal Durchschnitt" ~ "V2",
                           kcal_Visit == "V3 kcal Durchschnitt" ~ "V3"))%>%
  mutate(Comed_token= paste0(Patient_number, Visit))%>%
  select(Comed_token, kcal_Avg)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")
##Kcal body weight ratio per day
nutri%>%
  select(Patient_number, `V2_kcal_kg(BW)/d`:`V3_kcal_kg(BW)/d`)%>%
  gather(kcal_Visit, kcal_kg_day,`V2_kcal_kg(BW)/d`:`V3_kcal_kg(BW)/d`)%>%
  mutate(Visit = case_when(kcal_Visit == "V1_kcal_kg(BW)/d"  ~ "V1",
                           kcal_Visit == "V2_kcal_kg(BW)/d" ~ "V2",
                           kcal_Visit == "V3_kcal_kg(BW)/d" ~ "V3"))%>%
  mutate(Comed_token= paste0(Patient_number, Visit))%>%
  select(Comed_token, kcal_kg_day)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")
##Protein (%)
nutri%>%
  select(Patient_number, `V2_Eiweiß(%)`:`V3_Eiweiß(%)`)%>%
  gather(Prot_Visit, Protein,`V2_Eiweiß(%)`:`V3_Eiweiß(%)`)%>%
  mutate(Visit = case_when(Prot_Visit == "V1_Eiweiß(%)"  ~ "V1",
                           Prot_Visit == "V2_Eiweiß(%)" ~ "V2",
                           Prot_Visit == "V3_Eiweiß(%)" ~ "V3"))%>%
  mutate(Comed_token= paste0(Patient_number, Visit))%>%
  select(Comed_token, Protein)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")
##Lipids (%)
nutri%>%
  select(Patient_number,  `V2_Fett(%)`:`V3_Fett(%)`)%>%
  gather(Lip_Visit, Lipids,`V2_Fett(%)`:`V3_Fett(%)`)%>%
  mutate(Visit = case_when(Lip_Visit == "V1_Fett(%)"  ~ "V1",
                           Lip_Visit == "V2_Fett(%)" ~ "V2",
                           Lip_Visit == "V3_Fett(%)" ~ "V3"))%>%
  mutate(Comed_token= paste0(Patient_number, Visit))%>%
  select(Comed_token, Lipids)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")
##Carbohydrates (%CHO)
nutri%>%
  select(Patient_number, `V2_Kohlenhydrate(%)`:`V3_Kohlenhydrate(%)`)%>%
  gather(CHO_Visit, CHO,`V2_Kohlenhydrate(%)`:`V3_Kohlenhydrate(%)`)%>%
  mutate(Visit = case_when(CHO_Visit == "V1_Kohlenhydrate(%)"  ~ "V1",
                           CHO_Visit == "V2_Kohlenhydrate(%)" ~ "V2",
                           CHO_Visit == "V3_Kohlenhydrate(%)" ~ "V3"))%>%
  mutate(Comed_token= paste0(Patient_number, Visit))%>%
  select(Comed_token, CHO)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")
##Dietary fiber (%DFr)
nutri%>%
  select(Patient_number, c(`V2_Ballaststoffe(%)`,`V3_Ballaststoffe(%)`))%>%
  gather(DFr_Visit, DFr,c(`V2_Ballaststoffe(%)`,`V3_Ballaststoffe(%)`))%>%
  mutate(Visit = case_when(DFr_Visit == "V1_Ballaststoffe(%)"  ~ "V1",
                           DFr_Visit == "V2_Ballaststoffe(%)" ~ "V2",
                           DFr_Visit == "V3_Ballaststoffe(%)" ~ "V3"))%>%
  mutate(Comed_token= paste0(Patient_number, Visit))%>%
  select(Comed_token, DFr)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")
##Alcohol (%EtOH)
nutri%>%
  select(Patient_number, c(`V2_Alkohol(%)`,`V3_Alkohol(%)`))%>%
  gather(EtOH_Visit, EtOH,c(`V2_Alkohol(%)`,`V3_Alkohol(%)`))%>%
  mutate(Visit = case_when(EtOH_Visit == "V1_Alkohol(%)"  ~ "V1",
                           EtOH_Visit == "V2_Alkohol(%)" ~ "V2",
                           EtOH_Visit == "V3_Alkohol(%)" ~ "V3"))%>%
  mutate(Comed_token= paste0(Patient_number, Visit))%>%
  select(Comed_token, EtOH)-> tmp2

nutri<- left_join(tmp1, tmp2, by= "Comed_token")
rm(tmp1, tmp2)

##4)Genetic data
##Re-shape dataframe
##Select not redundant tables
genetics%>%
  dplyr::rename(Patient_number= 1)%>%
  select(c(Patient_number, Mutation))-> genotype

genotype%>%
  mutate(Severity = case_when(Mutation == "F508del/F508del"  ~ "Severe",
                              Mutation != "F508del/F508del" ~ "Mild"))-> genotype

rm(genetics)

##5) Medication data 
clinic%>%
  select(c(SampleID,antibiotics_inh:Cystagon))-> clinic

##Merge all the info
##Merge Nutritional data with Lung data

lung%>%
  select(-c(Patient_number,Visit))%>%
  left_join(nutri, lung, by= "Comed_token")-> tmp1

tmp1%>%
  select(-c(Patient_number,Visit))-> tmp1

setdiff(tmp1$Comed_token, tech$Comed_token)

left_join(tech, tmp1, by="Comed_token")-> metadata ##This is partial (clinical info not yet included)

rm(tmp1)

##Add genotype 
left_join(metadata, genotype, by="Patient_number")-> metadata

left_join(metadata, clinic, by="SampleID")-> metadata
saveRDS(metadata, "CF_project/exercise-cf-intervention/data/metadata_indexed.rds")
write.csv(metadata, "~/CF_project/exercise-cf-intervention/data/metadata_indexed.csv")
