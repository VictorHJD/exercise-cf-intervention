##Cystic-Fibrosis Microbiome Project (Mainz)
##Metadata preparation (Try to organize the data provided by the team in Mainz)
##Víctor Hugo Jarquín-Díaz 15.02.2021

library(vroom)
library(tidyverse)
library(readxl)

##Load the data
##Pre-organized data with adjustments to samples:
##Sample P3V3 stool was incorrectly assigned, it is actually P2V3 stool
##Sample P17V2 stool, it is P17V3 stool

data.mainz<- read_tsv("~/CF_project/Metadata/Sample_Metadata_combine_rebecca.csv") ##Contains technical data from DNA (Mainz)

#tech<- read_tsv("~/CF_project/Metadata/technical_Metadata.tsv") 
lung<- read_excel("~/CF_project/Metadata/LungFunction_PhysicalFitness.xlsx") ##Contains lung function data
nutri<- read_excel("~/CF_project/Metadata/Nutrition_BodyComposition.xlsx") ##Contains nutritional data
clinic<- read.csv("~/CF_project/Metadata/sample_data_indexed_medication.csv", sep = ";") ##Contains clinical data (Antibiotics, Bacterial isolates, etc)
genetics<- read_excel("~/CF_project/Metadata/KlinDaten080221.xlsx")
resp<- read_excel("~/CF_project/Metadata/Responder_nonResponder.xlsx") ##Response to interventions data (ask how were assigned the categories)
classifier<- read.csv("~/CF_project/Metadata/sample_data_indexed_classifier.csv") ##Classification phenotic, genotipic, severity, Pseudomonas, Sport
self.train<-  read.csv("CF_project/Metadata/IPAQ_self_reported_training.csv")
bacteria.cult<- read_excel("CF_project/Metadata/ClinicalColonization.xlsx")

##Select useful data and uniform columns
##1) Tech data

data.mainz%>%
  dplyr::select(c(1,2,5,7,11,13,16,17,27, 55:57))%>%
  dplyr::rename(SampleID= X.SampleID, Visit= Comed_visit, sex= "sex(m1_w2)")%>%
  dplyr::mutate(Visit = paste0("V", Visit))%>%
  dplyr::mutate(Patient_number= Comed_token)%>%
  dplyr::mutate(Patient_number = gsub("V\\d+", "\\1", basename(Patient_number)))%>%
  dplyr::  group_by(Patient_number)%>%
  dplyr::select(c(1:10,13))->data.mainz


##2)Lung function data
##Re-shape dataframe
lung%>%
  dplyr::rename(Patient_number= ID, sex= `sex (m1_w2)`)%>%
  dplyr::select(Patient_number, sex, age, ppFEV1_V1, ppFEV1_V2, ppFEV1_V3)%>%
  gather(ppFEV1_Visit, ppFEV1, ppFEV1_V1:ppFEV1_V3)%>%
  separate(ppFEV1_Visit, c("Test", "Visit"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, Visit))%>%
  dplyr::select(Patient_number, sex, age, Visit, Comed_token, ppFEV1)-> tmp1

lung%>%
  dplyr::rename(Patient_number= ID)%>%
  dplyr::select(Patient_number, ppFVC_V1:ppFVC_V3)%>%
  gather(ppFVC_Visit, ppFVC, ppFVC_V1:ppFVC_V3)%>%
  dplyr::mutate(Visit = case_when(ppFVC_Visit == "ppFVC_V1"  ~ "V1",
                                  ppFVC_Visit == "ppFVC_V2" ~ "V2",
                                  ppFVC_Visit == "ppFVC_V3" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, Visit))%>%
  dplyr::select(Comed_token, ppFVC)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")

lung%>%
  dplyr::rename(Patient_number= ID)%>%
  dplyr::select(Patient_number, DistanzV1:DistanzV3)%>%
  gather(Dist_Visit, Dist, DistanzV1:DistanzV3)%>%
  dplyr::mutate(Visit = case_when(Dist_Visit == "DistanzV1"  ~ "V1",
                                  Dist_Visit == "DistanzV2" ~ "V2",
                                  Dist_Visit == "DistanzV3" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, Visit))%>%
  dplyr::select(Comed_token, Dist)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")

lung%>%
  dplyr::rename(Patient_number= ID)%>%
  dplyr::select(Patient_number, `Peakpower (W/kg)V1`,`Peakpower (W/kg)V2`, `Peakpower (W/kg)V3`)%>%
  gather(Peak_power_Visit, Peak_power, `Peakpower (W/kg)V1`:`Peakpower (W/kg)V3`)%>%
  dplyr::mutate(Visit = case_when(Peak_power_Visit == "Peakpower (W/kg)V1"  ~ "V1",
                           Peak_power_Visit == "Peakpower (W/kg)V2" ~ "V2",
                           Peak_power_Visit == "Peakpower (W/kg)V3" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, Visit))%>%
  dplyr::select(Comed_token, Peak_power)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")

lung%>%
  dplyr::rename(Patient_number= ID)%>%
  dplyr::select(Patient_number, `VO2peak (ml/min/kg)V1`:`VO2peak (ml/min/kg)V3`)%>%
  gather(V02_Visit, V02_A, `VO2peak (ml/min/kg)V1`:`VO2peak (ml/min/kg)V3`)%>%
  dplyr::mutate(Visit = case_when(V02_Visit == "VO2peak (ml/min/kg)V1"  ~ "V1",
                           V02_Visit == "VO2peak (ml/min/kg)V2" ~ "V2",
                           V02_Visit == "VO2peak (ml/min/kg)V3" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, Visit))%>%
  dplyr::select(Comed_token, V02_A)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")

lung%>%
  dplyr::rename(Patient_number= ID)%>%
  dplyr::select(Patient_number, `VO2peak (l/min)V1`:`VO2peak (l/min)V3`)%>%
  gather(V02_Visit, V02_B, `VO2peak (l/min)V1`:`VO2peak (l/min)V3`)%>%
  dplyr::mutate(Visit = case_when(V02_Visit == "VO2peak (l/min)V1"  ~ "V1",
                           V02_Visit == "VO2peak (l/min)V2" ~ "V2",
                           V02_Visit == "VO2peak (l/min)V3" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, Visit))%>%
  dplyr::select(Comed_token, V02_B)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")

lung%>%
  dplyr::rename(Patient_number= ID)%>%
  dplyr::select(Patient_number, BMIV1:BMIV3)%>%
  gather(BMI_Visit, BMI, BMIV1:BMIV3)%>%
  dplyr::mutate(Visit = case_when(BMI_Visit == "BMIV1"  ~ "V1",
                           BMI_Visit == "BMIV2" ~ "V2",
                           BMI_Visit == "BMIV3" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, Visit))%>%
  dplyr::select(Comed_token, BMI)-> tmp2

lung<- left_join(tmp1, tmp2, by= "Comed_token") ##Index table usable in R
rm(tmp1, tmp2)

##3) Nutritional data
##Re-shape dataframe
##Eliminate BMI that was already included in the previous table (this value actually can be computed based on wight and lenght)
nutri%>%
  dplyr::rename(Patient_number= Patient)%>%
  dplyr::mutate(Patient_number = paste0("P", Patient_number))%>%
  dplyr::select(-c(V1_BMI:V3_BMI))-> nutri
##Fat free mass (Charatsi)
nutri%>%
  dplyr::select(Patient_number, `FFM_V1_Charatsi(kg)`:`FFM_V3_Charatsi(kg)`)%>%
  gather(FFM_Visit, FFM_Charatsi, `FFM_V1_Charatsi(kg)`:`FFM_V3_Charatsi(kg)`)%>%
  dplyr::mutate(Visit = case_when(FFM_Visit == "FFM_V1_Charatsi(kg)"  ~ "V1",
                           FFM_Visit == "FFM_V2_Charatsi(kg)" ~ "V2",
                           FFM_Visit == "FFM_V3_Charatsi(kg)" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, Visit))%>%
  dplyr::select(Patient_number, Visit, Comed_token, FFM_Charatsi)-> tmp1
##Lenght (cm)
nutri%>%
  dplyr::select(Patient_number, V1_length_cm:V3_length_cm)%>%
  gather(Length_Visit, Length, V1_length_cm:V3_length_cm)%>%
  dplyr::mutate(Visit = case_when(Length_Visit == "V1_length_cm"  ~ "V1",
                           Length_Visit == "V2_length_cm" ~ "V2",
                           Length_Visit == "V3_length_cm" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, Visit))%>%
  dplyr::select(Comed_token, Length)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")
##Weight (kg)
nutri%>%
  dplyr::select(Patient_number, V1_weight_kg:V3_weight_kg)%>%
  gather(Weight_Visit, Weight, V1_weight_kg:V3_weight_kg)%>%
  dplyr::mutate(Visit = case_when(Weight_Visit == "V1_weight_kg"  ~ "V1",
                           Weight_Visit == "V2_weight_kg" ~ "V2",
                           Weight_Visit == "V3_weight_kg" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, Visit))%>%
  dplyr::select(Comed_token, Weight)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")
##Fat free mass (Lukaski) bioelectrical impedance measurements
nutri%>%
  dplyr::select(Patient_number, `V1_FFM_%(Lukaski)`:`V3_FFM_%(Lukaski)`)%>%
  gather(FFM_Luk_Visit, FFM_Luk,`V1_FFM_%(Lukaski)`:`V3_FFM_%(Lukaski)`)%>%
  dplyr::mutate(Visit = case_when(FFM_Luk_Visit == "V1_FFM_%(Lukaski)"  ~ "V1",
                           FFM_Luk_Visit == "V2_FFM_%(Lukaski)" ~ "V2",
                           FFM_Luk_Visit == "V3_FFM_%(Lukaski)" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, Visit))%>%
  dplyr::select(Comed_token, FFM_Luk)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")
##Average of kcal
nutri%>%
  dplyr::select(Patient_number, `V2 kcal Durchschnitt`:`V3 kcal Durchschnitt`)%>%
  gather(kcal_Visit, kcal_Avg,`V2 kcal Durchschnitt`:`V3 kcal Durchschnitt`)%>%
  dplyr::mutate(Visit = case_when(kcal_Visit == "V1 kcal Durchschnitt"  ~ "V1",
                           kcal_Visit == "V2 kcal Durchschnitt" ~ "V2",
                           kcal_Visit == "V3 kcal Durchschnitt" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, Visit))%>%
  dplyr::select(Comed_token, kcal_Avg)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")
##Kcal body weight ratio per day
nutri%>%
  dplyr::select(Patient_number, `V2_kcal_kg(BW)/d`:`V3_kcal_kg(BW)/d`)%>%
  gather(kcal_Visit, kcal_kg_day,`V2_kcal_kg(BW)/d`:`V3_kcal_kg(BW)/d`)%>%
  dplyr::mutate(Visit = case_when(kcal_Visit == "V1_kcal_kg(BW)/d"  ~ "V1",
                           kcal_Visit == "V2_kcal_kg(BW)/d" ~ "V2",
                           kcal_Visit == "V3_kcal_kg(BW)/d" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, Visit))%>%
  dplyr::select(Comed_token, kcal_kg_day)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")
##Protein (%)
nutri%>%
  dplyr::select(Patient_number, `V2_Eiweiß(%)`:`V3_Eiweiß(%)`)%>%
  gather(Prot_Visit, Protein,`V2_Eiweiß(%)`:`V3_Eiweiß(%)`)%>%
  dplyr::mutate(Visit = case_when(Prot_Visit == "V1_Eiweiß(%)"  ~ "V1",
                           Prot_Visit == "V2_Eiweiß(%)" ~ "V2",
                           Prot_Visit == "V3_Eiweiß(%)" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, Visit))%>%
  dplyr::select(Comed_token, Protein)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")
##Lipids (%)
nutri%>%
  dplyr::select(Patient_number,  `V2_Fett(%)`:`V3_Fett(%)`)%>%
  gather(Lip_Visit, Lipids,`V2_Fett(%)`:`V3_Fett(%)`)%>%
  dplyr::mutate(Visit = case_when(Lip_Visit == "V1_Fett(%)"  ~ "V1",
                           Lip_Visit == "V2_Fett(%)" ~ "V2",
                           Lip_Visit == "V3_Fett(%)" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, Visit))%>%
  dplyr::select(Comed_token, Lipids)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")
##Carbohydrates (%CHO)
nutri%>%
  dplyr::select(Patient_number, `V2_Kohlenhydrate(%)`:`V3_Kohlenhydrate(%)`)%>%
  gather(CHO_Visit, CHO,`V2_Kohlenhydrate(%)`:`V3_Kohlenhydrate(%)`)%>%
  dplyr::mutate(Visit = case_when(CHO_Visit == "V1_Kohlenhydrate(%)"  ~ "V1",
                           CHO_Visit == "V2_Kohlenhydrate(%)" ~ "V2",
                           CHO_Visit == "V3_Kohlenhydrate(%)" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, Visit))%>%
  dplyr::select(Comed_token, CHO)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")
##Dietary fiber (%DFr)
nutri%>%
  dplyr::select(Patient_number, c(`V2_Ballaststoffe(%)`,`V3_Ballaststoffe(%)`))%>%
  gather(DFr_Visit, DFr,c(`V2_Ballaststoffe(%)`,`V3_Ballaststoffe(%)`))%>%
  dplyr::mutate(Visit = case_when(DFr_Visit == "V1_Ballaststoffe(%)"  ~ "V1",
                           DFr_Visit == "V2_Ballaststoffe(%)" ~ "V2",
                           DFr_Visit == "V3_Ballaststoffe(%)" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, Visit))%>%
  dplyr::select(Comed_token, DFr)-> tmp2

tmp1<- left_join(tmp1, tmp2, by= "Comed_token")
##Alcohol (%EtOH)
nutri%>%
  dplyr::select(Patient_number, c(`V2_Alkohol(%)`,`V3_Alkohol(%)`))%>%
  gather(EtOH_Visit, EtOH,c(`V2_Alkohol(%)`,`V3_Alkohol(%)`))%>%
  dplyr::mutate(Visit = case_when(EtOH_Visit == "V1_Alkohol(%)"  ~ "V1",
                           EtOH_Visit == "V2_Alkohol(%)" ~ "V2",
                           EtOH_Visit == "V3_Alkohol(%)" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, Visit))%>%
  dplyr::select(Comed_token, EtOH)-> tmp2

nutri<- left_join(tmp1, tmp2, by= "Comed_token")
rm(tmp1, tmp2)

##4)Severity classification
##Re-shape dataframe
##Select not redundant tables
classifier%>%
  dplyr::select(c(SampleID, Benzoase, Phenotype_severity, Pseudomonas_status, Sport_responderVO2max.5ml.min.kg, Mutation_severity))%>%
  dplyr::rename(Sport_Response= 5)-> severity

rm(classifier)
#genetics%>%
#  dplyr::rename(Patient_number= 1)%>%
#  select(c(Patient_number, Mutation))-> genotype

#genotype%>%
#  mutate(Severity = case_when(Mutation == "F508del/F508del"  ~ "Severe",
#                              Mutation != "F508del/F508del" ~ "Mild"))-> genotype

#rm(genetics)

##5) Medication data 
clinic%>%
  dplyr::select(c(SampleID,antibiotics_inh:Cystagon))-> clinic

##6) Responder/non Responder data
resp%>%
  dplyr::select(1:5)%>%
  dplyr::rename(Patient_number= 1, Nutrition_Response= 2, FFM_Response= 3, Sport_Response= 4, pFVC_Response= 5)%>%
  dplyr::select(1,2,3,5)-> resp

##Merge all the info
##Merge Nutritional data with Lung data
lung%>%
  dplyr::select(-c(Patient_number,Visit))%>%
  left_join(nutri, lung, by= "Comed_token")-> tmp1

write.csv(tmp1, "~/CF_project/exercise-cf-intervention/data/metadata_body_lung.csv", row.names = F) #--> #--> For body composition
saveRDS(tmp1, "CF_project/exercise-cf-intervention/data/metadata_body_lung.rds") #--> For body composition

tmp1%>%
  dplyr::select(-c(Patient_number,Visit))-> tmp1

setdiff(tmp1$Comed_token, data.mainz$Comed_token)

left_join(tmp1, data.mainz, by="Comed_token")-> metadata ##This is partial (clinical info not yet included)
left_join(data.mainz, tmp1,  by="Comed_token")-> metadata.ps ##This is partial (for phyloseq object)
rm(tmp1)

##Add response data 
left_join(metadata, resp, by="Patient_number")-> metadata
left_join(metadata.ps, resp, by="Patient_number")-> metadata.ps
##Add medication data
setdiff(metadata$SampleID, clinic$SampleID)
setdiff(metadata.ps$SampleID, clinic$SampleID)

##In clinic:
##SampleID 10P3V3 change to 10P2V3B in metadata
##SampleID 10P17V2 change to 10P17V3A in metadata 
##Make adjustment before mearging

clinic%>%
  dplyr::mutate(SampleID=replace(SampleID, SampleID=="10P3V3", "10P2V3B"))%>%
  dplyr::mutate(SampleID=replace(SampleID, SampleID=="10P17V2", "10P17V3A"))-> clinic

##Check that this is fixed
setdiff(metadata$SampleID, clinic$SampleID)
setdiff(metadata.ps$SampleID, clinic$SampleID)

##Now join medication data
left_join(metadata, clinic, by="SampleID")-> metadata
left_join(metadata.ps, clinic, by="SampleID")-> metadata.ps
##Add genotype 
left_join(metadata, severity, by="SampleID")-> metadata
left_join(metadata.ps, severity, by="SampleID")-> metadata.ps

metadata$FFM_Luk<- NULL
metadata$Nutrition_Response<- NULL
metadata$FFM_Response<- NULL
metadata$pFVC_Response<- NULL

metadata.ps$FFM_Luk<- NULL
metadata.ps$Nutrition_Response<- NULL
metadata.ps$FFM_Response<- NULL
metadata.ps$pFVC_Response<- NULL

write.csv(metadata, "~/CF_project/exercise-cf-intervention/data/metadata_indexed.csv", row.names = F)
saveRDS(metadata, "CF_project/exercise-cf-intervention/data/metadata_indexed.rds")

write.csv(metadata.ps, "~/CF_project/exercise-cf-intervention/data/metadata_PS.csv", row.names = F)
saveRDS(metadata.ps, "CF_project/exercise-cf-intervention/data/metadata_PS.rds")

rm(clinic, data.mainz, lung, nutri, severity, resp, genetics)

##Adjust self train data and save it. 

self.train%>%
  dplyr::select(Patient_number, V1_Train_Week, V2_Train_Week, V3_Train_Week)%>%
  gather(Train_Week_Visit, Train_Week, V1_Train_Week:V3_Train_Week)%>%
  dplyr::mutate(Visit = case_when(Train_Week_Visit == "V1_Train_Week"  ~ "V1",
                                  Train_Week_Visit == "V2_Train_Week" ~ "V2",
                                  Train_Week_Visit == "V3_Train_Week" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, "_", Visit))%>%
  dplyr::select(Comed_token, Train_Week)-> tmp2

self.train%>%
  dplyr::select(Patient_number, V1_Sessions_Week, V2_Sessions_Week, V3_Sessions_Week)%>%
  gather(Sessions_Week_Visit, Sessions_Week, V1_Sessions_Week:V3_Sessions_Week)%>%
  dplyr::mutate(Visit = case_when(Sessions_Week_Visit == "V1_Sessions_Week"  ~ "V1",
                                  Sessions_Week_Visit == "V2_Sessions_Week" ~ "V2",
                                  Sessions_Week_Visit == "V3_Sessions_Week" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, "_", Visit))%>%
  dplyr::select(Comed_token, Sessions_Week)-> tmp3

left_join(tmp2, tmp3, by="Comed_token")%>%
  separate(Comed_token, c("Patient_number", "Visit"))-> self.train

saveRDS(self.train, "CF_project/exercise-cf-intervention/data/self.train.rds")

##Adjust Bacterial culture data 
##Remove spaces in colnames 
colnames(bacteria.cult)<- gsub(" ", "_", colnames(bacteria.cult))

##Pseudomonas  
bacteria.cult%>%
  dplyr::rename(Patient_number= ID, V1_Pseudomonas_aeruginosa_MDR= `V1_Pseudomonas_aeruginosa_3/4MRGN`,
                V2_Pseudomonas_aeruginosa_MDR = `V2_Pseudomonas_aeruginosa_3/4MRGN`,
                V3_Pseudomonas_aeruginosa_MDR = `V3_Pseudomonas_aeruginosa_3/4MRGN`)%>%
  dplyr::select(Patient_number, V1_Pseudomonas_aeruginosa:V1_Pseudomonas_aeruginosa_MDR, 
                V2_Pseudomonas_aeruginosa:V2_Pseudomonas_aeruginosa_MDR,
                V3_Pseudomonas_aeruginosa:V3_Pseudomonas_aeruginosa_MDR)-> tmp

tmp[,2:10]<- lapply(tmp[,2:10], function(x) as.numeric(as.character(x)))

tmp%>%
  dplyr::mutate(V1_Pseudomonas_aeruginosa_HV = case_when(V1_Pseudomonas_aeruginosa_mucoid == 0 & V1_Pseudomonas_aeruginosa_MDR == 0  ~ 0,
                                                         V1_Pseudomonas_aeruginosa_mucoid == 1 & V1_Pseudomonas_aeruginosa_MDR == 0 ~ 1,
                                                         V1_Pseudomonas_aeruginosa_mucoid == 0 & V1_Pseudomonas_aeruginosa_MDR == 1 ~ 1,
                                                         V1_Pseudomonas_aeruginosa_mucoid == 1 & V1_Pseudomonas_aeruginosa_MDR == 1 ~ 1))%>%
  dplyr::mutate(V2_Pseudomonas_aeruginosa_HV = case_when(V2_Pseudomonas_aeruginosa_mucoid == 0 & V2_Pseudomonas_aeruginosa_MDR == 0  ~ 0,
                                                         V2_Pseudomonas_aeruginosa_mucoid == 1 & V2_Pseudomonas_aeruginosa_MDR == 0 ~ 1,
                                                         V2_Pseudomonas_aeruginosa_mucoid == 0 & V2_Pseudomonas_aeruginosa_MDR == 1 ~ 1,
                                                         V2_Pseudomonas_aeruginosa_mucoid == 1 & V2_Pseudomonas_aeruginosa_MDR == 1 ~ 1))%>%
  dplyr::mutate(V3_Pseudomonas_aeruginosa_HV = case_when(V3_Pseudomonas_aeruginosa_mucoid == 0 & V3_Pseudomonas_aeruginosa_MDR == 0  ~ 0,
                                                         V3_Pseudomonas_aeruginosa_mucoid == 1 & V3_Pseudomonas_aeruginosa_MDR == 0 ~ 1,
                                                         V3_Pseudomonas_aeruginosa_mucoid == 0 & V3_Pseudomonas_aeruginosa_MDR == 1 ~ 1,
                                                         V3_Pseudomonas_aeruginosa_mucoid == 1 & V3_Pseudomonas_aeruginosa_MDR == 1 ~ 1))%>%
  dplyr::select(Patient_number, V1_Pseudomonas_aeruginosa, V1_Pseudomonas_aeruginosa_HV, 
                V2_Pseudomonas_aeruginosa,V2_Pseudomonas_aeruginosa_HV,
                V3_Pseudomonas_aeruginosa, V3_Pseudomonas_aeruginosa_HV)%>%
  dplyr::mutate(V1_Pseudomonas_aeruginosa = case_when(V1_Pseudomonas_aeruginosa == 0 & V1_Pseudomonas_aeruginosa_HV == 0  ~ 0,
                                                      V1_Pseudomonas_aeruginosa == 1 & V1_Pseudomonas_aeruginosa_HV == 0 ~ 1,
                                                      V1_Pseudomonas_aeruginosa == 0 & V1_Pseudomonas_aeruginosa_HV == 1 ~ 1,
                                                      V1_Pseudomonas_aeruginosa == 1 & V1_Pseudomonas_aeruginosa_HV == 1 ~ 1))%>%
  dplyr::mutate(V2_Pseudomonas_aeruginosa = case_when(V2_Pseudomonas_aeruginosa == 0 & V2_Pseudomonas_aeruginosa_HV == 0  ~ 0,
                                                      V2_Pseudomonas_aeruginosa == 1 & V2_Pseudomonas_aeruginosa_HV == 0 ~ 1,
                                                      V2_Pseudomonas_aeruginosa == 0 & V2_Pseudomonas_aeruginosa_HV == 1 ~ 1,
                                                      V2_Pseudomonas_aeruginosa == 1 & V2_Pseudomonas_aeruginosa_HV == 1 ~ 1))%>%
  dplyr::mutate(V3_Pseudomonas_aeruginosa = case_when(V3_Pseudomonas_aeruginosa == 0 & V3_Pseudomonas_aeruginosa_HV == 0  ~ 0,
                                                      V3_Pseudomonas_aeruginosa == 1 & V3_Pseudomonas_aeruginosa_HV == 0 ~ 1,
                                                      V3_Pseudomonas_aeruginosa == 0 & V3_Pseudomonas_aeruginosa_HV == 1 ~ 1,
                                                      V3_Pseudomonas_aeruginosa == 1 & V3_Pseudomonas_aeruginosa_HV == 1 ~ 1))-> tmp
tmp%>%
  dplyr::select(Patient_number, V1_Pseudomonas_aeruginosa, 
                V2_Pseudomonas_aeruginosa,
                V3_Pseudomonas_aeruginosa)%>%
  gather(Pseudomonas_Visit, Pseudomonas_culture, c(V1_Pseudomonas_aeruginosa, V2_Pseudomonas_aeruginosa, V3_Pseudomonas_aeruginosa))%>%
  dplyr::mutate(Visit = case_when(Pseudomonas_Visit == "V1_Pseudomonas_aeruginosa"  ~ "V1",
                                  Pseudomonas_Visit == "V2_Pseudomonas_aeruginosa" ~ "V2",
                                  Pseudomonas_Visit == "V3_Pseudomonas_aeruginosa" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, "_", Visit))%>%
  dplyr::select(Comed_token, Pseudomonas_culture)%>%
  distinct()-> tmp1

tmp%>%
  dplyr::select(Patient_number, V1_Pseudomonas_aeruginosa_HV, 
                V2_Pseudomonas_aeruginosa_HV,
                V3_Pseudomonas_aeruginosa_HV)%>%
  gather(Pseudomonas_Visit, Pseudomonas_HV, c(V1_Pseudomonas_aeruginosa_HV, V2_Pseudomonas_aeruginosa_HV, V3_Pseudomonas_aeruginosa_HV))%>%
  dplyr::mutate(Visit = case_when(Pseudomonas_Visit == "V1_Pseudomonas_aeruginosa_HV"  ~ "V1",
                                  Pseudomonas_Visit == "V2_Pseudomonas_aeruginosa_HV" ~ "V2",
                                  Pseudomonas_Visit == "V3_Pseudomonas_aeruginosa_HV" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, "_", Visit))%>%
  dplyr::select(Comed_token, Pseudomonas_HV)%>%
  distinct()%>%
  left_join(tmp1, by= "Comed_token")-> tmp

##Staphylococcus
bacteria.cult%>%
  dplyr::rename(Patient_number= ID)%>%
  dplyr::select(Patient_number, V1_Staphylococcus_aureus, V2_Staphylococcus_aureus,
               V3_Staphylococcus_aureus)-> tmp1

tmp1[,2:4]<- lapply(tmp1[,2:4], function(x) as.numeric(as.character(x)))

tmp1%>%
  gather(Staphylococcus_Visit, Staphylococcus_aureus, c(V1_Staphylococcus_aureus, V2_Staphylococcus_aureus, V3_Staphylococcus_aureus))%>%
  dplyr::mutate(Visit = case_when(Staphylococcus_Visit == "V1_Staphylococcus_aureus"  ~ "V1",
                                  Staphylococcus_Visit == "V2_Staphylococcus_aureus" ~ "V2",
                                  Staphylococcus_Visit == "V3_Staphylococcus_aureus" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, "_", Visit))%>%
  dplyr::select(Comed_token, Staphylococcus_aureus)%>%
  left_join(tmp, by= "Comed_token")-> tmp

##Stenotrophomonas_maltophilia
bacteria.cult%>%
  dplyr::rename(Patient_number= ID)%>%
  dplyr::select(Patient_number, V1_Stenotrophomonas_maltophilia, V2_Stenotrophomonas_maltophilia,
                V3_Stenotrophomonas_maltophilia)-> tmp1

tmp1[,2:4]<- lapply(tmp1[,2:4], function(x) as.numeric(as.character(x)))

tmp1%>%
  gather(Stenotrophomonas_maltophilia_Visit, Stenotrophomonas_maltophilia, c(V1_Stenotrophomonas_maltophilia, 
                                                                             V2_Stenotrophomonas_maltophilia, 
                                                                             V3_Stenotrophomonas_maltophilia))%>%
  dplyr::mutate(Visit = case_when(Stenotrophomonas_maltophilia_Visit == "V1_Stenotrophomonas_maltophilia"  ~ "V1",
                                  Stenotrophomonas_maltophilia_Visit == "V2_Stenotrophomonas_maltophilia" ~ "V2",
                                  Stenotrophomonas_maltophilia_Visit == "V3_Stenotrophomonas_maltophilia" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, "_", Visit))%>%
  dplyr::select(Comed_token, Stenotrophomonas_maltophilia)%>%
  left_join(tmp, by= "Comed_token")-> tmp

##Candida
bacteria.cult%>%
  dplyr::rename(Patient_number= ID)%>%
  dplyr::select(Patient_number, V1_Candida_spp, V2_Candida_spp,
                V3_Candida_spp)-> tmp1

tmp1[,2:4]<- lapply(tmp1[,2:4], function(x) as.numeric(as.character(x)))

tmp1%>%
  gather(Candida_Visit, Candida_spp, c(V1_Candida_spp, 
                                       V2_Candida_spp, 
                                       V3_Candida_spp))%>%
  dplyr::mutate(Visit = case_when(Candida_Visit == "V1_Candida_spp"  ~ "V1",
                                  Candida_Visit == "V2_Candida_spp" ~ "V2",
                                  Candida_Visit == "V3_Candida_spp" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, "_", Visit))%>%
  dplyr::select(Comed_token, Candida_spp)%>%
  left_join(tmp, by= "Comed_token")-> tmp

##Aspergillus ---> Ask how they discarded potential environmental contaminants?
bacteria.cult%>%
  dplyr::rename(Patient_number= ID)%>%
  dplyr::select(Patient_number, V1_Aspergillus, V2_Aspergillus,
                V3_Aspergillus)-> tmp1

tmp1[,2:4]<- lapply(tmp1[,2:4], function(x) as.numeric(as.character(x)))

tmp1%>%
  gather(Aspergillus_Visit, Aspergillus_spp, c(V1_Aspergillus, 
                                       V2_Aspergillus, 
                                       V3_Aspergillus))%>%
  dplyr::mutate(Visit = case_when(Aspergillus_Visit == "V1_Aspergillus"  ~ "V1",
                                  Aspergillus_Visit == "V2_Aspergillus" ~ "V2",
                                  Aspergillus_Visit == "V3_Aspergillus" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, "_", Visit))%>%
  dplyr::select(Comed_token, Aspergillus_spp)%>%
  left_join(tmp, by= "Comed_token")-> tmp


##Burkholderia_spp
bacteria.cult%>%
  dplyr::rename(Patient_number= ID)%>%
  dplyr::select(Patient_number, V1_Burkholderia_spp, V2_Burkholderia_spp,
                V3_Burkholderia_spp)-> tmp1

tmp1[,2:4]<- lapply(tmp1[,2:4], function(x) as.numeric(as.character(x)))

tmp1%>%
  gather(Burkholderia_spp_Visit, Burkholderia_spp, c(V1_Burkholderia_spp, 
                                               V2_Burkholderia_spp, 
                                               V3_Burkholderia_spp))%>%
  dplyr::mutate(Visit = case_when(Burkholderia_spp_Visit == "V1_Burkholderia_spp"  ~ "V1",
                                  Burkholderia_spp_Visit == "V2_Burkholderia_spp" ~ "V2",
                                  Burkholderia_spp_Visit == "V3_Burkholderia_spp" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, "_", Visit))%>%
  dplyr::select(Comed_token, Burkholderia_spp)%>%
  left_join(tmp, by= "Comed_token")-> tmp

##Haemophilus influenzae
bacteria.cult%>%
  dplyr::rename(Patient_number= ID, V1_Haemophilus_influenzae = `V1_Haemophilus influenzae`,
              V2_Haemophilus_influenzae = `V2_Haemophilus influenzae`,
              V3_Haemophilus_influenzae = `V3_Haemophilus influenzae`)%>%
  dplyr::select(Patient_number, V1_Haemophilus_influenzae, 
                V2_Haemophilus_influenzae,
                V3_Haemophilus_influenzae)-> tmp1
tmp1[,2:4]<- lapply(tmp1[,2:4], function(x) as.numeric(as.character(x)))

tmp1%>%
  gather(Haemophilus_influenzae_Visit, Haemophilus_influenzae, c(V1_Haemophilus_influenzae, 
                                               V2_Haemophilus_influenzae, 
                                               V3_Haemophilus_influenzae))%>%
  dplyr::mutate(Visit = case_when(Haemophilus_influenzae_Visit == "V1_Haemophilus_influenzae"  ~ "V1",
                                  Haemophilus_influenzae_Visit == "V2_Haemophilus_influenzae" ~ "V2",
                                  Haemophilus_influenzae_Visit == "V3_Haemophilus_influenzae" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, "_", Visit))%>%
  dplyr::select(Comed_token, Haemophilus_influenzae)%>%
  left_join(tmp, by= "Comed_token")-> tmp

##Achromobacter_xylosoxidans
bacteria.cult%>%
  dplyr::rename(Patient_number= ID)%>%
  dplyr::select(Patient_number, V1_Achromobacter_xylosoxidans, 
                V2_Achromobacter_xylosoxidans,
                V3_Achromobacter_xylosoxidans)-> tmp1

tmp1[,2:4]<- lapply(tmp1[,2:4], function(x) as.numeric(as.character(x)))

tmp1%>%
  gather(Achromobacter_xylosoxidans_Visit, Achromobacter_xylosoxidans, c(V1_Achromobacter_xylosoxidans, 
                                               V2_Achromobacter_xylosoxidans, 
                                               V3_Achromobacter_xylosoxidans))%>%
  dplyr::mutate(Visit = case_when(Achromobacter_xylosoxidans_Visit == "V1_Achromobacter_xylosoxidans"  ~ "V1",
                                  Achromobacter_xylosoxidans_Visit == "V2_Achromobacter_xylosoxidans" ~ "V2",
                                  Achromobacter_xylosoxidans_Visit == "V3_Achromobacter_xylosoxidans" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, "_", Visit))%>%
  dplyr::select(Comed_token, Achromobacter_xylosoxidans)%>%
  left_join(tmp, by= "Comed_token")-> tmp

##Trichosporon inkin
bacteria.cult%>%
  dplyr::rename(Patient_number= ID)%>%
  dplyr::select(Patient_number, V1_Trichosporon_inkin, 
                V2_Trichosporon_inkin,
                V3_Trichosporon_inkin)-> tmp1

tmp1[,2:4]<- lapply(tmp1[,2:4], function(x) as.numeric(as.character(x)))

tmp1%>%
  gather(Trichosporon_inkin_Visit, Trichosporon_inkin, c(V1_Trichosporon_inkin, 
                                               V2_Trichosporon_inkin, 
                                               V3_Trichosporon_inkin))%>%
  dplyr::mutate(Visit = case_when(Trichosporon_inkin_Visit == "V1_Trichosporon_inkin"  ~ "V1",
                                  Trichosporon_inkin_Visit == "V2_Trichosporon_inkin" ~ "V2",
                                  Trichosporon_inkin_Visit == "V3_Trichosporon_inkin" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, "_", Visit))%>%
  dplyr::select(Comed_token, Trichosporon_inkin)%>%
  left_join(tmp, by= "Comed_token")-> tmp


##Mycobacterium
bacteria.cult%>%
  dplyr::rename(Patient_number= ID, V1_Mycobacterium_abcessus = `V1_Mycobacterium abcessus abcessus`,
                V2_Mycobacterium_abcessus = `V2_Mycobacterium abcessus abcessus`,
                V3_Mycobacterium_abcessus = `V3_Mycobacterium abcessus abcessus`)%>%
  dplyr::select(Patient_number, V1_Mycobacterium_abcessus, 
                V2_Mycobacterium_abcessus,
                V3_Mycobacterium_abcessus)-> tmp1

tmp1[,2:4]<- lapply(tmp1[,2:4], function(x) as.numeric(as.character(x)))

tmp1%>%
  gather(Mycobacterium_abcessus_Visit, Mycobacterium_abcessus, c(V1_Mycobacterium_abcessus, 
                                               V2_Mycobacterium_abcessus, 
                                               V3_Mycobacterium_abcessus))%>%
  dplyr::mutate(Visit = case_when(Mycobacterium_abcessus_Visit == "V1_Mycobacterium_abcessus"  ~ "V1",
                                  Mycobacterium_abcessus_Visit == "V2_Mycobacterium_abcessus" ~ "V2",
                                  Mycobacterium_abcessus_Visit == "V3_Mycobacterium_abcessus" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, "_", Visit))%>%
  dplyr::select(Comed_token, Mycobacterium_abcessus)%>%
  left_join(tmp, by= "Comed_token")-> tmp

##Enterobacteriaceae not included because screened just on V1 and V2 
##Fungi
bacteria.cult%>%
  dplyr::rename(Patient_number= ID)%>%
  dplyr::select(Patient_number, V1_fungi, 
                V2_fungi,
                V3_fungi)-> tmp1

tmp1[,2:4]<- lapply(tmp1[,2:4], function(x) as.numeric(as.character(x)))

tmp1%>%
  gather(Fungi_Visit, Other_Fungi, c(V1_fungi, V2_fungi, V3_fungi))%>%
  dplyr::mutate(Visit = case_when(Fungi_Visit == "V1_fungi"  ~ "V1",
                                  Fungi_Visit == "V2_fungi" ~ "V2",
                                  Fungi_Visit == "V3_fungi" ~ "V3"))%>%
  dplyr::mutate(Comed_token= paste0(Patient_number, "_", Visit))%>%
  dplyr::select(Comed_token, Other_Fungi)%>%
  left_join(tmp, by= "Comed_token")-> bacteria.cult

saveRDS(bacteria.cult, "CF_project/exercise-cf-intervention/data/bacteria.culture.data.rds")
rm(tmp, tmp1)