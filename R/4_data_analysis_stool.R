##Cystic-Fibrosis Microbiome Project (Mainz)
##Data analysis: Stool samples (Alpha, Beta diversity and differential abundance)
##Víctor Hugo Jarquín-Díaz 28.04.2021

##To phyloseq
library(phyloseq)
library(microbiome)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(data.table)
library(knitr)
library(grid)
library(gridExtra)
library(RColorBrewer)
##For analysis with linear models
library("lme4")
library("lmtest")
library("rcompanion")
library("car")
library("merTools")
library(sjPlot)
library(sjlabelled)
library(sjmisc)

reRun<- FALSE 

##Load data 
##Run from the root of the repo at  ~/CF_project/exercise-cf-intervention/

if(!exists("PS4.stool")){
  if(isTRUE(reRun)){
    source("~/CF_project/exercise-cf-intervention/R/3_data_analysis.R") ## Run the script at base directory of repository!   
  } 
}

##Color palette for patients ##
pal.CF<- c((brewer.pal(n = 8, name = "Dark2")), (brewer.pal(n = 11, name = "Paired")) )

pal.CF<- c("P1"="#1B9E77","P2"= "#D95F02","P3"= "#7570B3","P4"= "#E7298A","P5"= "#66A61E",
           "P6"="#E6AB02","P7"= "#A6761D","P8"= "#666666","P9"= "#A6CEE3","P10"= "#1F78B4",
           "P11"= "#B2DF8A","P12"= "#33A02C","P13"= "#FB9A99","P14"="#E31A1C","P15"= "#FDBF6F",
           "P16"= "#FF7F00","P17"= "#CAB2D6","P18"= "#6A3D9A","P19"= "#FFFF99")
##Stool#####################
sdt%>%
  dplyr::filter(material=="Stool")%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))-> sdt.stool
##Richness
##Plot 
sdt%>%
  dplyr::filter(material=="Stool")%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Visit, y= chao1))+
  geom_violin(color= "black", outlier.colour = "white", trim = F)+
  geom_point(shape=21, size=3, aes(fill= Patient_number), color= "black", alpha= 0.5)+
  xlab("Visit")+
  ylab("Richness (Chao1 Index)")+
  geom_line(aes(group = Patient_number), color= "gray")+
  scale_fill_manual(values = pal.CF)+
  labs(tag= "A)", fill= "Patient number")+
  theme_classic()+
  theme(text = element_text(size=16))-> A

sdt%>%
  dplyr::filter(material=="Sputum")%>%
  dplyr::filter(Benzoase==1)%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Visit, y= chao1))+
  geom_violin(color= "black", outlier.colour = "white", trim = F)+
  geom_point(shape=21, size=3, aes(fill= Patient_number), color= "black", alpha= 0.5)+
  xlab("Visit")+
  ylab("Richness (Chao1 Index)")+
  geom_line(aes(group = Patient_number), color= "gray")+
  scale_fill_manual(values = pal.CF)+
  labs(tag= "C)", fill= "Patient number")+
  theme_classic()+
  theme(text = element_text(size=16))-> C

##Shannon diversity 
##Plot 
sdt%>%
  dplyr::filter(material=="Stool")%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Visit, y= diversity_shannon))+
  geom_violin(color= "black", outlier.colour = "white", trim = F)+
  geom_point(shape=21, size=3, aes(fill= Patient_number), color= "black", alpha= 0.5)+
  xlab("Visit")+
  ylab("Diversity (Shannon Index)")+
  geom_line(aes(group = Patient_number), color= "gray")+
  scale_fill_manual(values = pal.CF)+
  labs(tag= "B)", fill= "Patient number")+
  theme_classic()+
  theme(text = element_text(size=16))-> B

sdt%>%
  dplyr::filter(material=="Sputum")%>%
  dplyr::filter(Benzoase==1)%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Visit, y= diversity_shannon))+
  geom_violin(color= "black", outlier.colour = "white", trim = F)+
  geom_point(shape=21, size=3, aes(fill= Patient_number), color= "black", alpha= 0.5)+
  xlab("Visit")+
  ylab("Diversity (Shannon Index)")+
  geom_line(aes(group = Patient_number), color= "gray")+
  scale_fill_manual(values = pal.CF)+
  labs(tag= "D)", fill= "Patient number")+
  theme_classic()+
  theme(text = element_text(size=16))-> D

E<- ggarrange(A, B, C, D, ncol=2, nrow=2, common.legend = TRUE, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Alpha_Material.pdf", plot = E, width = 10, height = 8)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Alpha_Material.png", plot = E, width = 10, height = 8)
rm(A,B,C,D,E)

###Diversity and lung function 
sdt%>%
  dplyr::filter(material=="Stool")%>%
  mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= diversity_shannon, y= ppFEV1, shape= Visit))+
  geom_point(size=3, aes(fill= Visit), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  geom_smooth(method=lm, se = T, aes(color= Visit))+
  theme_bw()+
  labs(tag= "A)")+
  xlab("Alpha diveristy (Shannon Index)")+
  ylab("Lung function (ppFEV1)")+
  stat_cor(method = "spearman", label.x = 2, label.y = 30)+ # Add sperman`s correlation coefficient
  theme(text = element_text(size=16), legend.position = "none")+
  facet_grid(rows = vars(Visit))-> A

sdt%>%
  dplyr::filter(material=="Stool")%>%
  mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= diversity_shannon, y= ppFEV1))+
  geom_point(size=3, aes(fill= Patient_number, shape= Visit), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  scale_fill_manual(values = pal.CF)+
  geom_smooth(method=lm, se = T, color= "black")+
  theme_bw()+
  labs(tag= "B)")+
  labs(fill = "Patient")+
  labs(shape = "Visit")+
  guides(fill = guide_legend(override.aes=list(shape=c(21)), ncol = 6), shape= guide_legend(nrow = 3))+
  xlab("Alpha diveristy (Shannon Index)")+
  ylab("Lung function (ppFEV1)")+
  stat_cor(method = "spearman", label.x = 2, label.y = 30)+ # Add sperman`s correlation coefficient
  theme(text = element_text(size=16), legend.position="bottom", legend.box = "horizontal")-> B

C<- grid.arrange(A,B, heights = c(3, 2))

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Alpha_Lung_Stool.pdf", plot = C, width = 10, height = 10)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Alpha_Lung_Stool.png", plot = C, width = 10, height = 10)

rm(A,B,C)

lf.model.stool1 <- lm(ppFEV1 ~ diversity_shannon, data= sdt.stool) ##Null
lf.model.stool2 <- lm(ppFEV1 ~ diversity_shannon * Visit, data= sdt.stool) ##Full

car::Anova(lf.model.stool1, type=3) 
car::Anova(lf.model.stool2, type=3) 

lrtest(lf.model.stool1, lf.model.stool2)

lf.model.stool.lsm <-
  lsmeans::lsmeans(lf.model.stool2,
                   pairwise~diversity_shannon:Visit,
                   adjust="fdr")

lf.model.stool.lsm$contrasts

##Bray-Curtis
BC_dist<- phyloseq::distance(PS4.stool,
                             method="bray", weighted=F)
ordination<- ordinate(PS4.stool,
                      method="PCoA", distance= BC_dist)

plot_ordination(PS4.stool, ordination)+ 
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(fill= Patient_number, shape= Visit), color= "black")+
  scale_shape_manual(values = c(21, 24, 22))+
  scale_fill_manual(values = pal.CF)+
  #geom_point(color= "black", size= 1.5)+
  labs(title = "Bray-Curtis dissimilariy",tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(fill = "Patient")+
  labs(shape = "Visit")+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> D

###Phenotype severity
plot_ordination(PS4.stool, ordination)+ 
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(fill= Patient_number, shape= as.factor(Phenotype_severity)), color= "black")+
  scale_shape_manual(values = c(25, 24), labels = c("Low severity", "High severity"))+
  scale_fill_manual(values = pal.CF)+
  labs(title = "Bray-Curtis dissimilariy stool",tag= "A)")+
  stat_ellipse(aes(color = as.factor(Phenotype_severity)))+
  scale_color_manual(values=c("#8A9045FF", "#800000FF"))+
  theme_bw()+
  theme(text = element_text(size=16))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= F)+
  labs(fill = "Patient")+
  labs(shape = "Phenotype severity")+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> PCo.Sev.Stool

saveRDS(PCo.Sev.Stool, "CF_project/exercise-cf-intervention/data/PCo.Sev.Stool.rds")

##Stratified for Patient number 
tmp<- row.names(PS4.stool@sam_data)

tmp<- sdt[rownames(sdt)%in%tmp, ]

BC.test.stool<- vegan::adonis(BC_dist~ Phenotype_severity + Mutation_severity + sex + age +  Visit + BMI,
                              permutations = 999, data = tmp, na.action = F, strata = tmp$Patient_number)

kable(BC.test.stool$aov.tab)

## Differences are not linked to severity phenotype or genotype. 
##Extract pairwise distances per patient
BC_dist.stool<- as.matrix(BC_dist)
tmp1<- cbind(sdt.stool, BC_dist.stool)

##V1 vs V2
tmp1%>%
  filter(Visit== "V1")%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  dplyr::select(-c(rownames(tmp1[tmp1$Visit=="V1",])))%>%
  dplyr::select(-c(rownames(tmp1[tmp1$Visit=="V3",])))-> V1vsV2.stool


x<- c(V1vsV2.stool["10P1V1","10P1V2"], V1vsV2.stool["10P3V1A","10P3V2A"], 
      V1vsV2.stool["10P4V1A","10P4V2A"], V1vsV2.stool["10P6V1A","10P6V2A"],
      V1vsV2.stool["10P7V1A","10P7V2"], V1vsV2.stool["10P8V1A","10P8V2A"], V1vsV2.stool["10P9V1B","10P9V2A"],
      V1vsV2.stool["10P13V1","10P13V2"], V1vsV2.stool["10P14V1A","10P14V2A"], V1vsV2.stool["10P15V1A","10P15V2A"],
      V1vsV2.stool["10P16V1","10P16V2"])

y<- c("P1", "P3", "P4", "P6", "P7", "P8", "P9", "P13", "P14", "P15", "P16")

tmp2<- data.frame(x,y)
tmp2[,3]<- "V1_V2"
colnames(tmp2)<- c("BC_dist", "Patient_number", "Group")

##V1 vs V3
tmp1%>%
  filter(Visit== "V1")%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  dplyr::select(-c(rownames(tmp1[tmp1$Visit=="V1",])))%>%
  dplyr::select(-c(rownames(tmp1[tmp1$Visit=="V2",])))-> V1vsV3.stool

x<- c(V1vsV3.stool["10P3V1A","10P3V3A"], 
      V1vsV3.stool["10P4V1A","10P4V3A"], V1vsV3.stool["10P6V1A","10P6V3A"],
      V1vsV3.stool["10P7V1A","10P7V3A"],  V1vsV3.stool["10P9V1B","10P9V3A"],
      V1vsV3.stool["10P14V1A","10P14V3A"], V1vsV3.stool["10P15V1A","10P15V3A"],
      V1vsV3.stool["10P17V1","10P17V3A"], V1vsV3.stool["10P18V1A","10P18V3A"])

y<- c("P3", "P4", "P6", "P7", "P9", "P14", "P15", "P17","P18")

tmp3<- data.frame(x,y)
tmp3[,3]<- "V1_V3"
colnames(tmp3)<- c("BC_dist", "Patient_number", "Group")

##V2 vs V3
tmp1%>%
  filter(Visit== "V2")%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  dplyr::select(-c(rownames(tmp1[tmp1$Visit=="V1",])))%>%
  dplyr::select(-c(rownames(tmp1[tmp1$Visit=="V2",])))-> V2vsV3.stool

x<- c(V2vsV3.stool["10P2V2B","10P2V3B"], V2vsV3.stool["10P3V2A","10P3V3A"], 
      V2vsV3.stool["10P4V2A","10P4V3A"], V2vsV3.stool["10P6V2A","10P6V3A"],
      V2vsV3.stool["10P7V2","10P7V3A"],  V2vsV3.stool["10P9V2A","10P9V3A"],
      V2vsV3.stool["10P14V2A","10P14V3A"], V2vsV3.stool["10P15V2A","10P15V3A"])

y<- c("P2", "P3", "P4", "P6", "P7", "P9", "P14", "P15")

tmp4<- data.frame(x,y)
tmp4[,3]<- "V2_V3"
colnames(tmp4)<- c("BC_dist", "Patient_number", "Group")

##rowbind the 3 dataframes

BC_dist.stool<- bind_rows(tmp2, tmp3, tmp4)
rm(tmp2, tmp3, tmp4)

##From metadata extract classifiers and severity status 
metadata%>%
  dplyr::filter(material=="Stool")%>%
  group_by(Patient_number)%>%
  distinct(Patient_number, .keep_all = TRUE)%>%
  dplyr::select(c(Patient_number, Phenotype_severity, 
                  Pseudomonas_status, Sport_Response, Mutation_severity))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))-> tmp2

BC_dist.stool%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  left_join(tmp2, by="Patient_number")-> BC_dist.stool

rm(tmp2)

##Add training information
training%>%
  dplyr::filter(material=="Stool")%>%
  group_by(Patient_number)%>%
  distinct(Patient_number, .keep_all = TRUE)%>%
  dplyr::select(c(Patient_number, Mean_MET_V1V2:Percentage_Trainingsweeks_n52))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))-> training.stool

training.stool%>%
  dplyr::select(-c(Mean_Trainingtime_womissingvalues_V1V2, Mean_Trainingfrequency_womissingvalues_V1V2))%>%
  gather(key = "Measurment", value = "Value",
         Mean_MET_V1V2:Percentage_Trainingsweeks_n52)%>%
  separate(Measurment, c("A", "Measurment", "Group"))-> tmp3

tmp3$Group<- gsub("n52", "V1V3", basename(tmp3$Group))

tmp3%>%
  dplyr::filter(A!= "Percentage")%>%
  dplyr::select(c(Patient_number, Measurment, Group, Value))%>%
  dplyr::mutate(Measurment= as.factor(Measurment))%>%
  spread(key = "Measurment", value = "Value")%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  dplyr::mutate(ID= paste0(Patient_number, Group))%>%
  dplyr::select(-c(Patient_number, Group))-> tmp3

tmp3$Patient_number<- NULL
tmp3$ID<- gsub('(V\\d+)(V\\d+)$', '\\1_\\2', basename(tmp3$ID))

BC_dist.stool%>%
  dplyr::mutate(ID= paste0(Patient_number, Group))%>%
  left_join(tmp3, by="ID")-> BC_dist.stool

##Extract ppFEV1 and estimate deltas per Visit combination 
metadata%>%
  dplyr::select(c(Comed_token, ppFEV1))%>%
  dplyr::mutate(Comed_token= gsub("^(.*)V", "\\1_V", Comed_token))%>%
  separate(Comed_token, c("Patient_number", "Visit"))%>%
  group_by(Patient_number)%>%
  dplyr::select(c(Patient_number, Visit, ppFEV1))%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  dplyr::mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  dplyr::distinct()%>%
  spread(Visit, ppFEV1)%>%
  dplyr::mutate(V1_V2= V1 - V2)%>%
  dplyr::mutate(V1_V3= V1 - V3)%>%
  dplyr::mutate(V2_V3= V2 - V3)%>%
  dplyr::select(c(Patient_number, V1_V2, V2_V3, V1_V3))%>%
  pivot_longer(!Patient_number, names_to = "Group", values_to = "ppFEV1")%>%
  dplyr::mutate(ID= paste0(Patient_number, Group))%>%
  dplyr::ungroup()%>%
  dplyr::select(c(ID, ppFEV1))-> tmp2

BC_dist.stool%>%
  left_join(tmp2, by="ID")-> BC_dist.stool

##Extract ppFVC and estimate deltas per Visit
metadata%>%
  dplyr::select(c(Comed_token, ppFVC))%>%
  dplyr::mutate(Comed_token= gsub("^(.*)V", "\\1_V", Comed_token))%>%
  separate(Comed_token, c("Patient_number", "Visit"))%>%
  group_by(Patient_number)%>%
  dplyr::select(c(Patient_number, Visit, ppFVC))%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  dplyr::mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  dplyr::distinct()%>%
  spread(Visit, ppFVC)%>%
  dplyr::mutate(V1_V2= V1 - V2)%>%
  dplyr::mutate(V1_V3= V1 - V3)%>%
  dplyr::mutate(V2_V3= V2 - V3)%>%
  dplyr::select(c(Patient_number, V1_V2, V2_V3, V1_V3))%>%
  pivot_longer(!Patient_number, names_to = "Group", values_to = "ppFVC")%>%
  dplyr::mutate(ID= paste0(Patient_number, Group))%>%
  dplyr::ungroup()%>%
  dplyr::select(c(ID, ppFVC))-> tmp2

BC_dist.stool%>%
  left_join(tmp2, by="ID")-> BC_dist.stool

##Add a time between visits (Overall for know but ask values per patient per period)
BC_dist.stool%>%
  dplyr::mutate(Months= case_when(Group == "V1_V3" ~ 12,
                                  Group == "V1_V2" ~ 3,
                                  Group == "V2_V3" ~ 19))-> BC_dist.stool

##Is visit impacting differences in composition by patient? 
BC_dist.stool%>% 
  wilcox_test(BC_dist ~ Group)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Group")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "~/CF_project/exercise-cf-intervention/tables/Q3_Sample_Visit_BC.csv")

BC_dist.stool%>%
  wilcox_effsize(BC_dist ~ Group)

##Plot 
BC_dist.stool%>%
  ggplot(aes(x= Group, y= BC_dist))+
  geom_boxplot(color="black", alpha= 0.5)+
  geom_point(shape=21, position=position_jitter(0.2), size=3, aes(fill= Patient_number), color= "black")+
  xlab("Visit")+
  ylab("Bray-Curtis dissimilarity")+
  labs(tag= "B)", caption = get_pwc_label(stats.test))+
  scale_fill_manual(values = pal.CF)+
  theme_classic()+
  theme(text = element_text(size=16), legend.position = "none")+
  stat_pvalue_manual(stats.test, hide.ns = F,label = "{p.adj}{p.adj.signif}")->E

f<-ggarrange(D, E, ncol=1, nrow=2, common.legend = TRUE, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Stool.pdf", plot = f, width = 10, height = 8)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Stool.png", plot = f, width = 10, height = 8)

rm(D,E,f)

##Correlation with training 
##Frequency
BC_dist.stool%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Trainingfrequency, y= BC_dist))+
  geom_point(size=2.5, aes(shape= Group, fill= Patient_number), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  xlab("Mean training frequency per period")+
  ylab("Bray-Curtis dissimilarity")+
  labs(tag= "A)")+
  theme_classic()+
  theme(text = element_text(size=16))+
  scale_fill_manual(values = pal.CF)+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(fill = "Patient")+
  labs(shape = "Visit period")+
  geom_smooth(method = lm, se=FALSE)-> A

##Time
BC_dist.stool%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Trainingtime, y= BC_dist))+
  geom_point(size=2.5, aes(shape= Group, fill= Patient_number), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  xlab("Mean training time per period")+
  ylab("Bray-Curtis dissimilarity")+
  labs(tag= "B)")+
  theme_classic()+
  theme(text = element_text(size=16))+
  scale_fill_manual(values = pal.CF)+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(fill = "Patient")+
  labs(shape = "Visit period")+
  geom_smooth(method = lm, se=FALSE)-> B

##ppFEV1
BC_dist.stool%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= ppFEV1, y= BC_dist))+
  geom_point(position=position_jitter(0.2), size=2.5, aes(shape= Group, fill= Patient_number), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  xlab("Difference in ppFEV1 between visits")+
  ylab("Bray-Curtis dissimilarity")+
  labs(tag= "C)")+
  theme_classic()+
  scale_fill_manual(values = pal.CF)+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(fill = "Patient")+
  labs(shape = "Visit period")+
  theme(text = element_text(size=16))+
  geom_smooth(method = lm, se=FALSE)-> C

##ppFVC
BC_dist.stool%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= ppFVC, y= BC_dist))+
  geom_point(size=2.5, aes(shape= Group, fill= Patient_number), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  xlab("Difference in ppFVC between visits")+
  ylab("Bray-Curtis dissimilarity")+
  labs(tag= "D)")+
  theme_classic()+
  scale_fill_manual(values = pal.CF)+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(fill = "Patient")+
  labs(shape = "Visit period")+
  theme(text = element_text(size=16))+
  geom_smooth(method = lm, se=FALSE)-> D

plot<-ggarrange(A, B, C, D, ncol=2, nrow=2, common.legend = TRUE, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Stool_Training.pdf", plot = plot, width = 10, height = 12)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Stool_Training.png", plot = plot, width = 10, height = 12)

rm(A,B,C,D, plot)

###Mixed effect models 
BC_dist.stool%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))-> BC_dist.stool
##Check for complete cases
BC_dist.stool%>%
  dplyr::select(ppFVC, Trainingfrequency, Trainingtime, ppFEV1, BC_dist, Patient_number, Group)%>%
  dplyr::filter(complete.cases(.))-> tmp

##qqPlots (Check whther our variables are normaly distributed)
qqPlot(BC_dist.stool$BC_dist)
qqPlot(BC_dist.stool$Trainingfrequency)
qqPlot(BC_dist.stool$Trainingtime)
qqPlot(BC_dist.stool$ppFEV1)
qqPlot(BC_dist.stool$ppFVC)

myLRTsignificanceFactors <- function(modnull, mod2, mod3, mod4, mod5, mod6, mod7, mod8, modFull){
  return(list(signif2 = lrtest(modnull, mod2),
              signif3 = lrtest(modnull, mod3),
              signif4 = lrtest(modnull, mod4),
              signif5 = lrtest(modnull, mod5),
              signif6 = lrtest(modnull, mod6),
              signif7 = lrtest(modnull, mod7),
              signif8 = lrtest(modnull, mod8),
              signifFull = lrtest(modnull, modFull)))
}

##Based on the experimental design we should have a nested model with Intervisit as main grouping factor and patient 
##However, not all the grouping factors are complete so we should treat it as a crossed model 
##with patient and intervisit group as individual random effects

tr0<- lmer(BC_dist ~ (1 | Patient_number) + (1| Group), data = tmp) ##Null model
tr1<- lmer(BC_dist ~ Trainingfrequency + (1 | Patient_number)+ (1| Group), data = tmp)
tr2<- lmer(BC_dist ~ Trainingtime + (1 | Patient_number)+ (1| Group), data = tmp)
tr3<- lmer(BC_dist ~ ppFEV1 + (1 | Patient_number)+ (1| Group), data = tmp)
tr4<- lmer(BC_dist ~ ppFVC + (1 | Patient_number)+ (1| Group), data = tmp)
tr5<- lmer(BC_dist ~ Trainingfrequency + Trainingtime + (1 | Patient_number)+ (1| Group), data = tmp)
tr6<- lmer(BC_dist ~ Trainingfrequency + Trainingtime + ppFEV1 + (1 | Patient_number)+ (1| Group), data = tmp)
tr7<- lmer(BC_dist ~ Trainingfrequency + Trainingtime + ppFVC + (1 | Patient_number)+ (1| Group), data = tmp)
tr8<- lmer(BC_dist ~ Trainingfrequency + Trainingtime + ppFVC + ppFEV1 + (1 | Patient_number)+ (1| Group), data = tmp)

myLRTsignificanceFactors(modnull =tr0, tr1, tr2, tr3, tr4, tr5, tr6, tr7, tr8)

##Each factor ad predictive value to the model 
##What happen with interactions 
tr9<-glmer(BC_dist ~ Trainingfrequency*Trainingtime*ppFVC + (1 | Patient_number) + (1| Group), data = BC_dist.stool)
summary(tr9)

tr10<-lmer(BC_dist ~ Trainingfrequency*Trainingtime*ppFEV1*ppFVC + (1 | Patient_number) + (1| Group), data = BC_dist.stool)
summary(tr10) ##Full model

lrtest(tr9, tr10)

A<- plotREsim(REsim(tr10))  ## plot the interval estimates
A$data%>%
  dplyr::mutate(groupID = fct_relevel(groupID, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))-> A$data
A+
  geom_point(shape= 21, size=2.5, aes(fill= groupID), color= "black")+
  scale_fill_manual(values = pal.CF)+
  xlab(label = NULL)+
  ylab(label = "Estimates")+
  labs(title = NULL, tag= "A)", fill= "Patient number")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        text = element_text(size=16))-> A

##PLot model  
B<- plot_model(tr10, p.adjust = "BH", vline.color = "gray",
               axis.labels = c( "Trainingfrequency"= "Frequency",
                                "Trainingtime" = "Time", 
                                "Trainingfrequency:Trainingtime" = "Frequency*Time", 
                                "Trainingfrequency:ppFEV1" = "Frequency*ppFEV1",
                                "Trainingtime:ppFEV1" = "Time*ppFEV1", 
                                "Trainingfrequency:ppFVC"= "Frequency*ppFVC", 
                                "Trainingtime:ppFVC"= "Time*ppFVC", 
                                "Trainingfrequency:Trainingtime:ppFEV1"= "Frequency*Time)*ppFEV1",
                                "Trainingfrequency:Trainingtime:ppFVC"= "(Frequency*Time)*ppFVC",
                                "ppFEV1:ppFVC"= "ppFEV1*ppFVC",
                                "Trainingfrequency:ppFEV1:ppFVC"= "(Frequency*ppFEV1)*ppFVC", 
                                "Trainingtime:ppFEV1:ppFVC"= "(Time*ppFEV1)*ppFVC",
                                "Trainingfrequency:Trainingtime:ppFEV1:ppFVC"= "(Frequency*Time*ppFEV1)*ppFVC"))+
  scale_y_continuous(limits = c(-0.4, 0.4))+
  geom_point(shape= 21, size=2.5, aes(fill= group), color= "black")+
  labs(title = NULL, tag= "B)")+
  theme_classic()+
  theme(text = element_text(size=16))
  
C<-ggarrange(A, B, ncol=1, nrow=2)

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Stool_Training_Effect_Ranges.png", plot = C, width = 8, height = 8)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Stool_Training_Effect_Ranges.pdf", plot = C, width = 8, height = 8)

rm(A, B, C)

##Test predictive value of BC differences within patient to lung function measurements Delta ppFEV1 between visits

##Mixed effect models
##with patient as random effect
tr0<- lmer(ppFVC ~ (1 | Patient_number)+ (1| Group), data = tmp) ##Null model
tr1<- lmer(ppFVC ~ Trainingfrequency + (1 | Patient_number) + (1| Group), data = tmp)
tr2<- lmer(ppFVC ~ Trainingtime + (1 | Patient_number) + (1| Group), data = tmp)
tr3<- lmer(ppFVC ~ ppFEV1 + (1 | Patient_number) + (1| Group), data = tmp)
tr4<- lmer(ppFVC ~ BC_dist + (1 | Patient_number) + (1| Group), data = tmp)
tr5<- lmer(ppFVC ~ Trainingfrequency + Trainingtime + (1 | Patient_number)+ (1| Group), data = tmp)
tr6<- lmer(ppFVC ~ Trainingfrequency + Trainingtime + ppFEV1 + (1 | Patient_number)+ (1| Group), data = tmp)
tr7<- lmer(ppFVC ~ Trainingfrequency + Trainingtime + BC_dist + (1 | Patient_number)+ (1| Group), data = tmp)
tr8<- lmer(ppFVC ~ Trainingfrequency + Trainingtime + BC_dist + ppFEV1 + (1 | Patient_number)+ (1| Group), data = tmp)

myLRTsignificanceFactors(modnull =tr0, tr1, tr2, tr3, tr4, tr5, tr6, tr7, tr8)

##The most significant predictors are ppFEV1 (model 3) and BC dissimilarity (model 4) (Check interaction)
##BC alone is a better predictor
lrtest(tr3, tr4)
##What about their interaction?
tr9<- lmer(ppFVC ~ ppFEV1*BC_dist + (1 | Patient_number) + (1| Group), data = tmp)
summary(tr9)

Model.Stool<- plot_model(tr9, p.adjust = "BH", vline.color = "gray", show.values = TRUE, value.offset = .3,
               axis.labels = c( "Trainingfrequency"= "Frequency",
                                "Trainingtime" = "Time", 
                                "BC_dist"="Bray-Curtis",
                                "Trainingfrequency:Trainingtime" = "Frequency*Time", 
                                "Trainingfrequency:ppFEV1" = "Frequency*ppFEV1",
                                "Trainingtime:ppFEV1" = "Time*ppFEV1", 
                                "Trainingfrequency:BC_dist"= "Frequency*Bray-Curtis", 
                                "Trainingtime:BC_dist"= "Time*Bray-Curtis", 
                                "Trainingfrequency:Trainingtime:ppFEV1"= "Frequency*Time)*ppFEV1",
                                "Trainingfrequency:Trainingtime:BC_dist"= "(Frequency*Time)*Bray-Curtis",
                                "ppFEV1:BC_dist"= "ppFEV1*Bray-Curtis",
                                "Trainingfrequency:ppFEV1:BC_dist"= "(Frequency*ppFEV1)*Bray-Curtis", 
                                "Trainingtime:ppFEV1:BC_dist"= "(Time*ppFEV1)*Bray-Curtis",
                                "Trainingfrequency:Trainingtime:ppFEV1:BC_dist"= "(Frequency*Time*ppFEV1)*Bray-Curtis"))+
  scale_y_continuous(limits = c(-30, 30))+
  geom_point(shape= 21, size=2.5, aes(fill= group), color= "black")+
  labs(title = NULL, tag= "A)")+
  theme_classic()+
  theme(text = element_text(size=16))


###Does differences in bacterial composition within patient predict severity status 
BC_dist.stool%>%
  dplyr::mutate(Phenotype_severity = case_when(Phenotype_severity == 2  ~ 1,
                                               Phenotype_severity == 1 ~ 0))%>%
  dplyr::mutate(Mutation_severity = case_when(Mutation_severity == 2  ~ 1,
                                              Mutation_severity == 1 ~ 0))-> BC_dist.stool

##Logistic regression 
log.model1 <-glmer(Phenotype_severity ~Trainingfrequency + Trainingtime +
                     ppFVC + ppFEV1 + BC_dist + (1 | Patient_number), data = BC_dist.stool, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)##Full model

summary(log.model1)$coef

log.model2 <- glmer(Phenotype_severity ~Trainingfrequency + Trainingtime +
                      ppFVC + ppFEV1 + (1 | Patient_number), data = BC_dist.stool, family = binomial, 
                    control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)##Without BC dissimilarity model

summary(log.model2)$coef

lrtest(log.model1, log.model2)  ##Bray-Curtis dissimilarity do not add predictive power for phenotype prediction

BC_dist.stool%>%
  ggplot(aes(ppFEV1, Phenotype_severity)) +
  geom_point(size=2.5, aes(shape= Group, fill= Patient_number), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  scale_fill_manual(values = pal.CF)+
  labs(x = "Difference in ppFEV1 between visits",
       y = "Probability of severe CF phenotype", tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(fill = "Patient")+
  labs(shape = "Visit interval")+
  geom_smooth(method = "glm", method.args = list(family = "binomial"))

log.model3 <-glmer(Sport_Response ~Trainingfrequency + Trainingtime +
                     ppFVC + ppFEV1 + BC_dist + (1 | Patient_number), data = BC_dist.stool, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)##Full model

summary(log.model3)$coef

log.model4 <- glmer(Sport_Response ~Trainingfrequency + Trainingtime +
                      ppFVC + ppFEV1 + (1 | Patient_number), data = BC_dist.stool, family = binomial, 
                    control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)##Without BC dissimilarity model

summary(log.model4)$coef

lrtest(log.model3, log.model4) 

##Microbial composition do not add predictive power response to sports intervention
###Naive correlation with nutritional and respiratory activity
##Glom by genus
PS.stool.Gen<- tax_glom(PS4.stool, "Genus", NArm = T)

##Adjust ASV table for merging with taxa information
otu<- PS.stool.Gen@otu_table
otu<-as.data.frame(otu)
otu<- rownames_to_column(otu, var = "ASV")

##Select just the genus 
tax<- PS.stool.Gen@tax_table
tax<-as.data.frame(tax)
tax%>%
  dplyr::select(Genus, Phylum)-> tax
tax$Genus<-gsub(" ", "_", basename(tax$Genus))
tax<- rownames_to_column(tax, var = "ASV")

##Use genus as rownames
stool.microbiome<- plyr::join(otu, tax, by= "ASV")
stool.microbiome$ASV<- NULL
rownames(stool.microbiome)<- paste0(stool.microbiome$Phylum, "-", stool.microbiome$Genus)
stool.microbiome$Genus<- NULL
stool.microbiome$Phylum<- NULL

##Transpose dataframe so samples are rows 
stool.microbiome<- t(stool.microbiome)
#saveRDS(stool.microbiome, "~/CF_project/exercise-cf-intervention/data/Stool_rare_ASV.rds")#--> For MetadeconfoundR
##Select useful metrics
y<-sdt.stool

y<- y[,c(31:95)] ##Eliminate non-numeric variables

##Make an adjustment to make visit and patient numeric
y$Visit<- as.numeric(gsub("V", "\\1", y$Visit))
y$Patient_number<- as.numeric(gsub("P", "\\1", y$Patient_number))

##Transform severity to binary 
y%>%
  dplyr::mutate(Phenotype_severity = case_when(Phenotype_severity == 2  ~ 1,
                                        Phenotype_severity == 1 ~ 0))%>%
  dplyr::mutate(Mutation_severity = case_when(Mutation_severity == 2  ~ 1,
                                       Mutation_severity == 1 ~ 0))-> y


y <- as.matrix(y) # Metadata (39 samples x 65 respiratory and nutritional variables)
#saveRDS(y, "~/CF_project/exercise-cf-intervention/data/Stool_rare_Metadata.rds")#--> For MetadeconfoundR

# Cross correlate data sets, output is also available in a handy table format
###Let's run this with metadeconfoundR (other script for that)

###Deseq2 analysis
library("DESeq2")

##Get raw data to run this
PS3.stool<- subset_samples(PS3, material%in%c("Stool"))
PS3.stool<- tax_glom(PS3.stool, "Genus")
##Adjustment make phenotype and genotype as factor 
PS3.stool@sam_data$Phenotype_severity <- as.factor(PS3.stool@sam_data$Phenotype_severity)
PS3.stool@sam_data$Mutation_severity <- as.factor(PS3.stool@sam_data$Mutation_severity)
##Severity classification based on:
#genotype
deseq.severity<- phyloseq_to_deseq2(PS3.stool, ~ Mutation_severity)

# calculate geometric means prior to estimate size factors
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans<- apply(counts(deseq.severity), 1, gm_mean)
deseq.severity<- estimateSizeFactors(deseq.severity, geoMeans = geoMeans)
deseq.severity<- DESeq(deseq.severity, fitType="local")

ac.res <- results(deseq.severity)

##Select cut-off value
sigtab <- cbind(as(ac.res,"data.frame"), as(tax_table(PS3.stool)[rownames(ac.res), ], "matrix"))
sigtab$Species<- NULL
head(sigtab,20)

##Volcano plot to detect differential taxa in stool microbiome between severity 
#The significantly deferentially abundant genes are the ones found upper-left and upper-right corners
##Add a column to the data specifying if they are highly (positive) or lowly abundant (negative)
## Considering the comparison D21 vs D0

# add a column of Non-significant
sigtab$AbundLev <- "NS"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "High" 
sigtab$AbundLev[sigtab$log2FoldChange > 0.6 & sigtab$padj < 0.001] <- "High"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
sigtab$AbundLev[sigtab$log2FoldChange < -0.6 & sigtab$padj < 0.001] <- "Low"

#Organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
#plot adding up all layers we have seen so far
sigtab%>%
  dplyr::mutate(Genus=  case_when(Genus == "CAG-56"  ~ "Firmicutes CAG:56",
                                  Genus == "[Eubacterium] ruminantium group" ~ "Eubacterium ruminantium group",
                                  Genus == "UCG-004" ~ "Lachnospiraceae UCG:004",
                                  TRUE ~ Genus))%>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(shape=21, size=3, alpha= 0.5, aes(fill= AbundLev), color= "black")+
  ggrepel::geom_text_repel(aes(col=AbundLev, label=Genus)) +
  scale_fill_manual(values=c("#8A9045FF", "#800000FF", "#767676FF")) +
  scale_color_manual(values=c("#8A9045FF", "#800000FF", "#767676FF")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.001), col="black", linetype= "dashed") +
  labs(tag= "A)", x= "log2 Fold change", y= "-Log10 (p Adjusted)", fill= "Genus\nabundance")+
  theme_bw()+
  theme(text = element_text(size=16))+
  guides(color= F)-> A

##Save table 
#write.csv(sigtab, "~/CF_project/exercise-cf-intervention/tables/Q5_DeSeq2_Abund_Stool_Severity.csv") #Genotype
#write.csv(sigtab, "~/CF_project/exercise-cf-intervention/tables/Q5_DeSeq2_Abund_Stool_Gen_Severity.csv") #Phenotype

##Phenotype
deseq.severity<- phyloseq_to_deseq2(PS3.stool, ~ Phenotype_severity)
# calculate geometric means prior to estimate size factors
geoMeans<- apply(counts(deseq.severity), 1, gm_mean)
deseq.severity<- estimateSizeFactors(deseq.severity, geoMeans = geoMeans)
deseq.severity<- DESeq(deseq.severity, fitType="local")

ac.res <- results(deseq.severity)

##Select cut-off value
sigtab <- cbind(as(ac.res,"data.frame"), as(tax_table(PS3.stool)[rownames(ac.res), ], "matrix"))
sigtab$Species<- NULL
head(sigtab,20)

# add a column of Non-significant
sigtab$AbundLev <- "NS"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "High" 
sigtab$AbundLev[sigtab$log2FoldChange > 0.6 & sigtab$padj < 0.001] <- "High"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
sigtab$AbundLev[sigtab$log2FoldChange < -0.6 & sigtab$padj < 0.001] <- "Low"

#Organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
#plot adding up all layers we have seen so far
sigtab%>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(shape=21, size=3, alpha= 0.5, aes(fill= AbundLev), color= "black")+
  ggrepel::geom_text_repel(aes(col=AbundLev, label=Genus)) +
  scale_fill_manual(values=c("#8A9045FF", "#800000FF", "#767676FF")) +
  scale_color_manual(values=c("#8A9045FF", "#800000FF", "#767676FF")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.001), col="black", linetype= "dashed") +
  labs(tag= "B)", x= "log2 Fold change", y= "-Log10 (p Adjusted)", fill= "Genus\nabundance")+
  theme_bw()+
  theme(text = element_text(size=16))+
  guides(color= F)-> B

C<-ggarrange(A, B, ncol=1, nrow=2, common.legend = T, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q5_DefAbund_Stool_Severity.png", plot = C, width = 8, height = 8)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q5_DefAbund_Stool_Severity.pdf", plot = C, width = 8, height = 8)

rm(A, B, C)

###Analysis by time points
###Make Phloseq subsets to run the analysis
##Adjustment make phenotype and genotype as factor 
PS3.stool@sam_data$Visit <- as.factor(PS3.stool@sam_data$Visit)
##V1V2
PS3.stool12<- subset_samples(PS3.stool, Visit%in%c("V1", "V2"))
##V2V3
PS3.stool23<- subset_samples(PS3.stool, Visit%in%c("V2", "V3"))
##V1V3
PS3.stool13<- subset_samples(PS3.stool, Visit%in%c("V1", "V3"))

##Visit 2 vs 1
deseq.visit<- phyloseq_to_deseq2(PS3.stool12, ~ Visit)

# calculate geometric means prior to estimate size factors
geoMeans<- apply(counts(deseq.visit), 1, gm_mean)
deseq.visit<- estimateSizeFactors(deseq.visit, geoMeans = geoMeans)
deseq.visit<- DESeq(deseq.visit, fitType="local")

ac.res <- results(deseq.visit)

##Select cut-off value
sigtab <- cbind(as(ac.res,"data.frame"), as(tax_table(PS3.stool12)[rownames(ac.res), ], "matrix"))
sigtab$Species<- NULL
head(sigtab,20)

# add a column of Non-significant
sigtab$AbundLev <- "NS"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "High" 
sigtab$AbundLev[sigtab$log2FoldChange > 0.6 & sigtab$padj < 0.001] <- "High"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
sigtab$AbundLev[sigtab$log2FoldChange < -0.6 & sigtab$padj < 0.001] <- "Low"

#Organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
#plot adding up all layers we have seen so far
sigtab%>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(shape=21, size=3, alpha= 0.5, aes(fill= AbundLev), color= "black")+
  ggrepel::geom_text_repel(aes(col=AbundLev, label=Genus)) +
  scale_fill_manual(values=c("#8A9045FF", "#800000FF", "#767676FF")) +
  scale_color_manual(values=c("#8A9045FF", "#800000FF", "#767676FF")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.001), col="black", linetype= "dashed") +
  labs(tag= "A)", x= "log2 Fold change", y= "-Log10 (p Adjusted)", fill= "Genus abundance")+
  theme_bw()+
  theme(text = element_text(size=16), legend.position = "top")+
  guides(color= F)-> A

##Extract the legend from A to use it later as a common legend 
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(A)

##Remove legend from A 
A <- A + theme(legend.position="none")

##Visit 3 vs 2
deseq.visit<- phyloseq_to_deseq2(PS3.stool23, ~ Visit)

# calculate geometric means prior to estimate size factors
geoMeans<- apply(counts(deseq.visit), 1, gm_mean)
deseq.visit<- estimateSizeFactors(deseq.visit, geoMeans = geoMeans)
deseq.visit<- DESeq(deseq.visit, fitType="local")

ac.res <- results(deseq.visit)

##Add taxa information
sigtab <- cbind(as(ac.res,"data.frame"), as(tax_table(PS3.stool23)[rownames(ac.res), ], "matrix"))
sigtab$Species<- NULL
head(sigtab,20)

# add a column of Non-significant
sigtab$AbundLev <- "NS"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "High" 
sigtab$AbundLev[sigtab$log2FoldChange > 0.6 & sigtab$padj < 0.001] <- "High"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
sigtab$AbundLev[sigtab$log2FoldChange < -0.6 & sigtab$padj < 0.001] <- "Low"

#Organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
#plot adding up all layers we have seen so far
sigtab%>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(shape=21, size=3, alpha= 0.5, aes(fill= AbundLev), color= "black")+
  ggrepel::geom_text_repel(aes(col=AbundLev, label=Genus)) +
  scale_fill_manual(values=c("#800000FF", "#767676FF")) +
  scale_color_manual(values=c("#800000FF", "#767676FF")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.001), col="black", linetype= "dashed") +
  labs(tag= "B)", x= "log2 Fold change", y= "-Log10 (p Adjusted)", fill= "Genus\nabundance")+
  theme_bw()+
  theme(text = element_text(size=16), legend.position = "none")+
  guides(color= F)-> B

##Visit 3 vs 1
deseq.visit<- phyloseq_to_deseq2(PS3.stool13, ~ Visit)

# calculate geometric means prior to estimate size factors
geoMeans<- apply(counts(deseq.visit), 1, gm_mean)
deseq.visit<- estimateSizeFactors(deseq.visit, geoMeans = geoMeans)
deseq.visit<- DESeq(deseq.visit, fitType="local")

ac.res <- results(deseq.visit)

##Add taxa information
sigtab <- cbind(as(ac.res,"data.frame"), as(tax_table(PS3.stool13)[rownames(ac.res), ], "matrix"))
sigtab$Species<- NULL
head(sigtab,20)

# add a column of Non-significant
sigtab$AbundLev <- "NS"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "High" 
sigtab$AbundLev[sigtab$log2FoldChange > 0.6 & sigtab$padj < 0.001] <- "High"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
sigtab$AbundLev[sigtab$log2FoldChange < -0.6 & sigtab$padj < 0.001] <- "Low"

#Organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
#plot adding up all layers we have seen so far
sigtab%>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(shape=21, size=3, alpha= 0.5, aes(fill= AbundLev), color= "black")+
  ggrepel::geom_text_repel(aes(col=AbundLev, label=Genus)) +
  scale_fill_manual(values=c("#800000FF", "#767676FF")) +
  scale_color_manual(values=c("#800000FF", "#767676FF")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.001), col="black", linetype= "dashed") +
  labs(tag= "C)", x= "log2 Fold change", y= "-Log10 (p Adjusted)", fill= "Genus\nabundance")+
  theme_bw()+
  theme(text = element_text(size=16), legend.position = "none")+
  guides(color= F)-> C

D<- grid.arrange(legend, A,B,C, nrow=4, heights=c(0.5, 2.5, 2.5, 2.5))

##Save just when the three objects are in the environment
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q5_DefAbund_Stool_Visit.png", plot = D, width = 8, height = 8)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q5_DefAbund_Stool_Visit.pdf", plot = D, width = 8, height = 8)

rm(A,B,C,D)

##Response to sports classification:
#sports
PS3.stool@sam_data$Sport_Response <- as.factor(PS3.stool@sam_data$Sport_Response)
deseq.sports<- phyloseq_to_deseq2(PS3.stool, ~ Sport_Response)

geoMeans<- apply(counts(deseq.sports), 1, gm_mean)
deseq.sports<- estimateSizeFactors(deseq.sports, geoMeans = geoMeans)
deseq.sports<- DESeq(deseq.sports, fitType="local")

ac.res <- results(deseq.sports)

##Select cut-off value
sigtab <- cbind(as(ac.res,"data.frame"), as(tax_table(PS3.stool)[rownames(ac.res), ], "matrix"))
sigtab$Species<- NULL
head(sigtab,20)

##Volcano plot to detect differential taxa in stool microbiome between severity 
#The significantly deferentially abundant genes are the ones found upper-left and upper-right corners
##Add a column to the data specifying if they are highly (positive) or lowly abundant (negative)
## Considering the comparison D21 vs D0

# add a column of Non-significant
sigtab$AbundLev <- "NS"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "High" 
sigtab$AbundLev[sigtab$log2FoldChange > 0.6 & sigtab$padj < 0.001] <- "High"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
sigtab$AbundLev[sigtab$log2FoldChange < -0.6 & sigtab$padj < 0.001] <- "Low"

#Organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
#plot adding up all layers we have seen so far
sigtab%>%
  dplyr::mutate(Genus=  case_when(Genus == "CAG-56"  ~ "Firmicutes CAG:56",
                                  Genus == "[Eubacterium] eligens group" ~ "Eubacterium eligens group",
                                  Genus == "UCG-004" ~ "Lachnospiraceae UCG:004",
                                  TRUE ~ Genus))%>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(shape=21, size=3, alpha= 0.5, aes(fill= AbundLev), color= "black")+
  ggrepel::geom_text_repel(aes(col=AbundLev, label=Genus)) +
  scale_fill_manual(values=c("#800000FF", "#767676FF")) +
  scale_color_manual(values=c("#800000FF", "#767676FF")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.001), col="black", linetype= "dashed") +
  labs(tag= "A)", x= "log2 Fold change", y= "-Log10 (p Adjusted)", fill= "Genus\nabundance")+
  theme_bw()+
  theme(text = element_text(size=16))+
  guides(color= F)-> A

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q5_DefAbund_Stool_Sports.png", plot = A, width = 8, height = 8)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q5_DefAbund_Stool_Sports.pdf", plot = A, width = 8, height = 8)

rm(A)

##Create biom format object for PICRUSt2
require("biomformat")
asvmat.rare<- as.matrix(PS3.stool@otu_table)
biom.tmp<- make_biom(asvmat.rare, matrix_element_type = "int")
write_biom(biom.tmp,"CF_project/exercise-cf-intervention/data/biom_stool.biom") ##Good biom for test

##Select sequences from the ASV in PS3.stool
library(Biostrings)
dna<- readDNAStringSet( "~/CF_project/output/ASV.fasta", format = "fasta")
keep <- data.frame(name = rownames(asvmat.rare))
names(dna)
dna.stool<- dna[keep$name]
writeXStringSet(dna.stool, "CF_project/exercise-cf-intervention/data/Stool_ASV.fasta") #-> For Picrust2
