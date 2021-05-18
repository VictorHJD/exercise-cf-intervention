##Cystic-Fibrosis Microbiome Project (Mainz)
##Lung funtion and body composition by patient
##Víctor Hugo Jarquín-Díaz 18.05.2021

##Load libraries
library(ggpubr)
library(tidyverse)

tmp1<- readRDS("CF_project/exercise-cf-intervention/data/metadata_body_lung.rds") #--> For body composition

###Antrophometric data and lung function 
##Define colors for patients

pal.CF<- c("P1"="#1B9E77","P2"= "#D95F02","P3"= "#7570B3","P4"= "#E7298A","P5"= "#66A61E",
           "P6"="#E6AB02","P7"= "#A6761D","P8"= "#666666","P9"= "#A6CEE3","P10"= "#1F78B4",
           "P11"= "#B2DF8A","P12"= "#33A02C","P13"= "#FB9A99","P14"="#E31A1C","P15"= "#FDBF6F",
           "P16"= "#FF7F00","P17"= "#CAB2D6","P18"= "#6A3D9A","P19"= "#FFFF99")

##Extract required data
tmp1%>%
  mutate(Patient_number = 
           fct_relevel(Patient_number, 
                       "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                       "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  group_by(Patient_number)%>%
  distinct(Visit, .keep_all = TRUE)%>%
  dplyr::select(c(Patient_number, Visit, 
                  ppFEV1, ppFVC, BMI, FFM_Charatsi))%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  filter(!is.na(ppFEV1)) -> df

##Plot
ggplot(df, aes(x= Visit, y= ppFEV1))+
  geom_boxplot(alpha= 0.5, color= "black")+
  geom_point(shape=21, size=3, aes(fill= Patient_number), color= "black")+
  scale_fill_manual(values = pal.CF)+
  xlab("Visit")+
  ylab("Lung function (ppFEV1)")+
  labs(fill= "Patient number")+
  theme_bw()+ 
  scale_x_discrete(labels=c("V1" = paste0("V1 (N=", sum(df$Visit=="V1"), ")"), 
                            "V2" = paste0("V2 (N=", sum(df$Visit=="V2"), ")"),
                            "V3" = paste0("V3 (N=", sum(df$Visit=="V3"), ")")))+
  theme(text = element_text(size=16))-> A

 
ggplot(df, aes(x= Visit, y= ppFVC))+
  geom_boxplot(alpha= 0.5, color= "black")+
  geom_point(shape=21, size=3, aes(fill= Patient_number), color= "black")+
  scale_fill_manual(values = pal.CF)+
  xlab("Visit")+
  labs(fill= "Patient number")+
  ylab("Lung function (ppFVC)")+
  theme_bw()+ 
  scale_x_discrete(labels=c("V1" = paste0("V1 (N=", sum(df$Visit=="V1"), ")"), 
                            "V2" = paste0("V2 (N=", sum(df$Visit=="V2"), ")"),
                            "V3" = paste0("V3 (N=", sum(df$Visit=="V3"), ")")))+
  theme(text = element_text(size=16))-> B

ggplot(df, aes(x= Visit, y= BMI))+
  geom_boxplot(alpha= 0.5, color= "black")+
  geom_point(shape=21, size=3, aes(fill= Patient_number), color= "black")+
  scale_fill_manual(values = pal.CF)+
  xlab("Visit")+
  labs(fill= "Patient number")+
  ylab("Body Mass Index (BMI)")+
  theme_bw()+ 
  scale_x_discrete(labels=c("V1" = paste0("V1 (N=", sum(df$Visit=="V1"), ")"), 
                            "V2" = paste0("V2 (N=", sum(df$Visit=="V2"), ")"),
                            "V3" = paste0("V3 (N=", sum(df$Visit=="V3"), ")")))+
  theme(text = element_text(size=16))-> C

ggplot(df, aes(x= Visit, y= FFM_Charatsi))+
  geom_boxplot(alpha= 0.5, color= "black")+
  geom_point(shape=21, size=3, aes(fill= Patient_number), color= "black")+
  scale_fill_manual(values = pal.CF)+
  xlab("Visit")+
  labs(fill= "Patient number")+
  ylab("Fat-free mass (FFM Charatsi kg)")+
  theme_bw()+ 
  scale_x_discrete(labels=c("V1" = paste0("V1 (N=", sum(df$Visit=="V1" & !is.na(df$FFM_Charatsi)), ")"), 
                            "V2" = paste0("V2 (N=", sum(df$Visit=="V2"), ")"),
                            "V3" = paste0("V3 (N=", sum(df$Visit=="V3"), ")")))+
  theme(text = element_text(size=16))-> D

E<-ggarrange(A,B,C,D, ncol=2, nrow=2, common.legend = TRUE, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q0_Lung_Body.pdf", plot = E, width = 10, height = 8)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q0_Lung_Body.png", plot = E, width = 10, height = 8)

rm(A,B,C,D,E)