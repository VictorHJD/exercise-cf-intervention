##Cystic-Fibrosis Microbiome Project (Mainz)
##Data analysis: Sputum samples (Beta diversity and differential abundance)
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

##Load data 
##Run from the root of the repo at  ~/CF_project/exercise-cf-intervention/

if(!exists("PS4.sputum")){
  if(isTRUE(reRun)){
    source("~/CF_project/exercise-cf-intervention/R/3_data_analysis.R") ## Run the script at base directory of repository!   
  } 
}

training<- read_tsv("~/CF_project/Metadata/sample_data_indexed_training.csv")

##Color palette for patients ##
pal.CF<- c((brewer.pal(n = 8, name = "Dark2")), (brewer.pal(n = 11, name = "Paired")) )

pal.CF<- c("P1"="#1B9E77","P2"= "#D95F02","P3"= "#7570B3","P4"= "#E7298A","P5"= "#66A61E",
           "P6"="#E6AB02","P7"= "#A6761D","P8"= "#666666","P9"= "#A6CEE3","P10"= "#1F78B4",
           "P11"= "#B2DF8A","P12"= "#33A02C","P13"= "#FB9A99","P14"="#E31A1C","P15"= "#FDBF6F",
           "P16"= "#FF7F00","P17"= "#CAB2D6","P18"= "#6A3D9A","P19"= "#FFFF99")

######Sputum###################
sdt%>%
  dplyr::filter(Benzoase==1)%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))-> sdt.sputum
##Bray-Curtis
BC_dist<- phyloseq::distance(PS4.sput,
                             method="bray", weighted=F)
ordination<- ordinate(PS4.sput,
                      method="PCoA", distance= BC_dist)

plot_ordination(PS4.sput, ordination, shape= "Visit")+ 
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
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> A

##Stratified for Patient number 
BC.test.sputum<- vegan::adonis(BC_dist~ Phenotype_severity+ Mutation_severity + sex + age +  Visit + BMI,
                               permutations = 999, data = sdt.sputum, na.action = F, strata = sdt.sputum$Patient_number)

##BMI significant predictor explaining 4% of the variation

##Extract pairwise distances per patient
BC_dist.sputum<- as.matrix(BC_dist)
tmp1<- cbind(sdt.sputum, BC_dist.sputum)

##V1 vs V2
tmp1%>%
  filter(Visit== "V1")%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  dplyr::select(-c(rownames(tmp1[tmp1$Visit=="V1",])))%>%
  dplyr::select(-c(rownames(tmp1[tmp1$Visit=="V3",])))-> V1vsV2.sputum

x<- c(V1vsV2.sputum["10P2V1","10P2V2A"], 
      V1vsV2.sputum["10P3V1","10P3V2"], V1vsV2.sputum["10P4V1","10P4V2"],
      V1vsV2.sputum["10P6V1","10P6V2"], V1vsV2.sputum["10P8V1","10P8V2"], 
      V1vsV2.sputum["10P9V1","10P9V2"],
      V1vsV2.sputum["10P10V1","10P10V2"], V1vsV2.sputum["10P11V1","10P11V2"], V1vsV2.sputum["10P14V1","10P14V2"],
      V1vsV2.sputum["10P15V1","10P15V2"], V1vsV2.sputum["10P18V1","10P18V2"])

y<- c("P2", "P3", "P4", "P6", "P8", "P9","P10", "P11", "P14", "P15", "P18")

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
  dplyr::select(-c(rownames(tmp1[tmp1$Visit=="V2",])))-> V1vsV3.sputum

x<- c(V1vsV3.sputum["10P2V1","10P2V3"], 
      V1vsV3.sputum["10P4V1","10P4V3"], V1vsV3.sputum["10P6V1","10P6V3"],
      V1vsV3.sputum["10P7V1","10P7V3"],  V1vsV3.sputum["10P9V1","10P9V3"],
      V1vsV3.sputum["10P14V1","10P14V3"], V1vsV3.sputum["10P15V1","10P15V3"],
      V1vsV3.sputum["10P18V1","10P18V3"])

y<- c("P2", "P4", "P6", "P7", "P9","P14", "P15", "P18")

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
  dplyr::select(-c(rownames(tmp1[tmp1$Visit=="V2",])))-> V2vsV3.sputum

x<- c(V2vsV3.sputum["10P2V2","10P2V3"], 
      V2vsV3.sputum["10P4V2","10P4V3"], V2vsV3.sputum["10P6V2","10P6V3"],
      V2vsV3.sputum["10P9V2","10P9V3"],  V2vsV3.sputum["10P14V2","10P14V3"],
      V2vsV3.sputum["10P15V2","10P15V3"], V2vsV3.sputum["10P18V2","10P18V3"])

y<- c("P2", "P4", "P6", "P9", "P14", "P15", "P18")

tmp4<- data.frame(x,y)
tmp4[,3]<- "V2_V3"
colnames(tmp4)<- c("BC_dist", "Patient_number", "Group")

##rowbind the 3 dataframes

BC_dist.sputum<- bind_rows(tmp2, tmp3, tmp4)

metadata%>%
  dplyr::filter(Benzoase==1)%>%
  group_by(Patient_number)%>%
  distinct(Patient_number, .keep_all = TRUE)%>%
  dplyr::select(c(Patient_number, Phenotype_severity, 
                  Pseudomonas_status, Sport_Response, Mutation_severity))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))-> tmp2

BC_dist.sputum%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  left_join(tmp2, by="Patient_number")-> BC_dist.sputum

##Add training information
training%>%
  dplyr::filter(Benzoase==1)%>%
  group_by(Patient_number)%>%
  distinct(Patient_number, .keep_all = TRUE)%>%
  dplyr::select(c(Patient_number, Mean_MET_V1V2:Percentage_Trainingsweeks_n52))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))-> training.sputum

training.sputum%>%
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

BC_dist.sputum%>%
  dplyr::mutate(ID= paste0(Patient_number, Group))%>%
  left_join(tmp3, by="ID")-> BC_dist.sputum

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

BC_dist.sputum%>%
  left_join(tmp2, by="ID")-> BC_dist.sputum

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

BC_dist.sputum%>%
  left_join(tmp2, by="ID")-> BC_dist.sputum

##Add a time between visits (Overall for know but ask values per patient per period)
BC_dist.sputum%>%
  dplyr::mutate(Months= case_when(Group == "V1_V3" ~ 12,
                                  Group == "V1_V2" ~ 3,
                                  Group == "V2_V3" ~ 19))-> BC_dist.sputum
#### By visit
BC_dist.sputum%>% 
  wilcox_test(BC_dist ~ Group)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Group")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "~/CF_project/exercise-cf-intervention/tables/Q3_Sample_Visit_BC_Sputum.csv")

BC_dist.sputum%>%
  wilcox_effsize(BC_dist ~ Group)

##Plot 
BC_dist.sputum%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Group, y= BC_dist))+
  geom_boxplot(color="black", alpha= 0.5)+
  geom_point(shape=21, position=position_jitter(0.2), size=3, aes(fill= Patient_number), color= "black")+
  xlab("Visit")+
  ylab("Bray-Curtis dissimilarity")+
  labs(tag= "B)", caption = get_pwc_label(stats.test))+
  scale_fill_manual(values = pal.CF)+
  theme_classic()+
  theme(text = element_text(size=16), legend.position = "none")+
  stat_pvalue_manual(stats.test, hide.ns = F,label = "{p.adj}{p.adj.signif}")->B

C<-ggarrange(A, B, ncol=1, nrow=2, common.legend = TRUE, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Sputum.pdf", plot = C, width = 10, height = 8)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Sputum.png", plot = C, width = 10, height = 8)

rm(A,B,C)

##Correlation with training 
##Frequency
BC_dist.sputum%>%
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
BC_dist.sputum%>%
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
BC_dist.sputum%>%
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

##ppFEV1
BC_dist.sputum%>%
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

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Sputum_Training.pdf", plot = plot, width = 10, height = 12)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Sputum_Training.png", plot = plot, width = 10, height = 12)

rm(A,B,C,D)

###Mixed effect models 
##Check for complete cases
#BC_dist.sputum%>% 
#  group_by(Patient_number)%>%
#  arrange(Group, .by_group = TRUE)%>%
#  summarise(n(), .groups = "keep")%>%
#  dplyr::rename(n = "n()")%>%
#  filter(n == 3)-> Keep

#Keep<- Keep$Patient_number

##Select just patients in Keep
#BC_dist.sputum[BC_dist.sputum$Patient_number%in%Keep, ]-> x

BC_dist.sputum%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))-> BC_dist.sputum
##qqPlots
qqPlot(BC_dist.sputum$BC_dist)
qqPlot(BC_dist.sputum$Trainingfrequency)
qqPlot(BC_dist.sputum$Trainingtime)

##GLM (Training weeks was omited due to a lot of NAs)
tr0<- glm(BC_dist ~ 1, data = BC_dist.sputum, na.action = na.exclude) ##Null model
tr1<- glm(BC_dist ~ Trainingfrequency, data =BC_dist.sputum, na.action = na.exclude)
tr2<- glm(BC_dist ~ Trainingtime, data = BC_dist.sputum, na.action = na.exclude)
tr3<- glm(BC_dist ~ Trainingfrequency + Trainingtime, data = BC_dist.sputum, na.action = na.exclude)
tr4<- glm(BC_dist ~ Trainingfrequency*Trainingtime , data =BC_dist.sputum, na.action = na.exclude) ##Full model

##Comparisons between models
lrtest(tr1, tr2) ##Significant difference Training Frequency seems to be a better predictor for microbial differences among visits
lrtest(tr1, tr3) ##Not significant difference 
lrtest(tr2, tr3) ##Not significant difference 
lrtest(tr1, tr4) ##No difference between model with training frequency alone and its interaction with time
lrtest(tr2, tr4) ##No difference between model with training time alone and its interaction with frequency
lrtest(tr3, tr4) ##No difference between model with training time and frequency and their interaction 

##Mixed effect models
##with patient as random effect
tr5<-lmer(BC_dist ~ Trainingfrequency * Trainingtime + (1 | Patient_number), data = BC_dist.sputum)
summary(tr5)
confint(tr5) ##Confidence interval for each fixed effect
ranef(tr5)$Patient_number ##Estimates for random effect
predictInterval(tr5)  ## for various model predictions, possibly with new data
REsim(tr5) ## mean, median and sd of the random effect estimates

lrtest(tr4, tr5)

A<- plotREsim(REsim(tr5))  ## plot the interval estimates
A$data%>%
  dplyr::mutate(groupID = fct_relevel(groupID, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))-> A$data
A+
  geom_point(aes(color= groupID))+
  xlab(label = NULL)+
  labs(title = NULL, color = "Patiente \n number")+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(), text = element_text(size=16))-> A

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Sputum_Training_Effect_Ranges.png", plot = A, width = 8, height = 8)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Sputum_Training_Effect_Ranges.pdf", plot = A, width = 8, height = 8)

rm(A)

##Nutritional, pFVC, FFM response with patient random effect
##Adjust to binary 1 responders = Improvement, 0 No responders= Stable, decreased or not detected
BC_dist.sputum%>%
  dplyr::mutate(Nutrition_Response = case_when(Nutrition_Response == 1  ~ 1,
                                               Nutrition_Response %in% c(0, 2,3) ~ 0))%>%
  dplyr::mutate(FFM_Response = case_when(FFM_Response == 1  ~ 1,
                                         FFM_Response %in% c(0, 2,3) ~ 0))%>%
  dplyr::mutate(pFVC_Response = case_when(pFVC_Response == 1  ~ 1,
                                          pFVC_Response %in% c(0, 2,3) ~ 0))-> BC_dist.sputum

##Does the nutrition response is linked to changes in microbial composition?
p1<- glm(Nutrition_Response ~ BC_dist, data = BC_dist.sputum, na.action = na.exclude, 
         family = binomial()) ## NS

##Does improvement on fat free mass response is linked to changes in microbial composition?
p2<- glm(FFM_Response ~ BC_dist, data = BC_dist.sputum, na.action = na.exclude,
         family = binomial()) ## NS

##Does improvement in predicted forced vital capacity linked to changes in microbiota or training?
p3<- glm(pFVC_Response ~ BC_dist, data = BC_dist.sputum, na.action = na.exclude,
         family = binomial())
p4<- glm(pFVC_Response ~ Trainingfrequency, data = BC_dist.sputum, na.action = na.exclude,
         family = binomial())
p5<- glm(pFVC_Response ~ Trainingtime, data = BC_dist.sputum, na.action = na.exclude,
         family = binomial())
p6<- glm(pFVC_Response ~ BC_dist + Trainingfrequency, data = BC_dist.sputum, na.action = na.exclude,
         family = binomial())
p7<- glm(pFVC_Response ~ BC_dist + Trainingtime, data = BC_dist.sputum, na.action = na.exclude,
         family = binomial())
p8<- glm(pFVC_Response ~ Trainingfrequency + Trainingtime, data = BC_dist.sputum, na.action = na.exclude,
         family = binomial())
p9<- glm(pFVC_Response ~ BC_dist + Trainingfrequency + Trainingtime, data = BC_dist.sputum, na.action = na.exclude,
         family = binomial())
p10<- glm(pFVC_Response ~ BC_dist*Trainingfrequency*Trainingtime, data = BC_dist.sputum, na.action = na.exclude,
          family = binomial())

##Glom by genus
PS.sputum.Gen<- tax_glom(PS4.sput, "Genus", NArm = T)

##Adjust ASV table for merging with taxa information
otu<- PS.sputum.Gen@otu_table
otu<-as.data.frame(otu)
otu<- rownames_to_column(otu, var = "ASV")

##Select just the genus 
tax<- PS.sputum.Gen@tax_table
tax<-as.data.frame(tax)
tax%>%
  dplyr::select(Genus)-> tax
tax$Genus<-gsub(" ", "_", basename(tax$Genus))
tax<- rownames_to_column(tax, var = "ASV")

##Use genus as rownames
sputum.microbiome<- plyr::join(otu, tax, by= "ASV")
sputum.microbiome$ASV<- NULL
rownames(sputum.microbiome)<- sputum.microbiome$Genus
sputum.microbiome$Genus<- NULL

##Transpose dataframe so samples are rows 
sputum.microbiome<- t(sputum.microbiome)
#saveRDS(sputum.microbiome, "~/CF_project/exercise-cf-intervention/data/Sput_rare_ASV.rds") #--> For MetadeconfoundR

x <- log10(sputum.microbiome+1) # ASV Log10 (39 samples x 205 genera)

##Select useful metrics
y<-sdt.sputum

y<- y[,c(26:27,29:98)] ##Eliminate non-numeric variables

##Make an adjustment to visit and patient to make it numeric
y$Visit<- as.numeric(gsub("V", "\\1", y$Visit))
y$Patient_number<- as.numeric(gsub("P", "\\1", y$Patient_number))

##Transform severity to binary 
y%>%
  mutate(Phenotype_severity = case_when(Phenotype_severity == 2  ~ 1,
                                        Phenotype_severity == 1 ~ 0))%>%
  mutate(Mutation_severity = case_when(Mutation_severity == 2  ~ 1,
                                       Mutation_severity == 1 ~ 0))-> y


y <- as.matrix(y) # Metadata (39 samples x 72 technical, respiratory and nutritional variables)
#saveRDS(y, "~/CF_project/exercise-cf-intervention/data/Sput_rare_Metadata.rds") #--> For MetadeconfoundR

# Cross correlate data sets, output is also available in a handy table format
correlation.table <- associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
kable(head(correlation.table))

heat(correlation.table, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05) 

ggplot(correlation.table, aes(X1, X2, group=X2)) + 
  geom_tile(aes(fill = Correlation)) +
  geom_text(aes(fill = correlation.table$Correlation, label = round(correlation.table$Correlation, 1)), size = 5) +
  scale_fill_gradientn("Correlation", 
                       breaks = seq(from = -1, to = 1,  by = 0.25), 
                       colours = c("blue", "white", "red"), 
                       limits = c(-1, 1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(x = "", y = "", tag = "B)")

###Deseq2 analysis
##Get raw data to run this
PS2.sputum<- subset_samples(PS2, Benzoase%in%c(1))
PS2.sputum<- tax_glom(PS2.sputum, "Genus")

##Severity classification based on:
#genotype
deseq.severity<- phyloseq_to_deseq2(PS2.sputum, ~ Mutation_severity)

##Phenotype
#deseq.severity<- phyloseq_to_deseq2(PS2.sputum, ~ Phenotype_severity)

# calculate geometric means prior to estimate size factors
geoMeans<- apply(counts(deseq.severity), 1, gm_mean)
deseq.severity<- estimateSizeFactors(deseq.severity, geoMeans = geoMeans)
deseq.severity<- DESeq(deseq.severity, fitType="local")

ac.res <- results(deseq.severity)
##Remove NA
ac.res <- ac.res[order(ac.res$padj, na.last=NA), ]

##Select cut-off value
alpha <- 0.05
Bac.sigtab <- ac.res[(ac.res$padj < alpha), ]

Bac.sigtab <- cbind(as(Bac.sigtab, "data.frame"), as(tax_table(PS2.sputum)[rownames(Bac.sigtab), ], "matrix"))

##Adjust value
Bac.posigtab <- Bac.sigtab[Bac.sigtab[, "log2FoldChange"] > -7, ]
Bac.posigtab <- Bac.posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

##Save table 
#write.csv(Bac.posigtab, "~/CF_project/exercise-cf-intervention/tables/Q5_DeSeq2_Abund_Sputum_Gen_Severity.csv") #Genotype
#write.csv(Bac.posigtab, "~/CF_project/exercise-cf-intervention/tables/Q5_DeSeq2_Abund_Sputum_Phen_Severity.csv") #Phenotype

###Negative Binomial in Microbiome Differential Abundance Testing (plot)
ggplot(na.omit(Bac.posigtab), aes(x=reorder(Genus, -log2FoldChange), y=log2FoldChange, color=Phylum)) +
  geom_point(shape=21, position=position_jitter(0.2), size=3, aes(fill= Phylum), color= "black") + 
  scale_fill_brewer(palette = "Set1")+
  coord_flip()+
  geom_hline(aes(yintercept = 0), color = "gray70", size = 0.6)+
  xlab("ASVs Genus-level")+
  ylab("Mild <-- Log-2-Fold-Change --> Severe")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5), text = element_text(size=16), 
        axis.text.y = element_text(face="italic", color="black"))+
  labs(tag = "B)")-> B

#png("CF_project/exercise-cf-intervention/figures/Q5_DefAbund_Sputum_Severity.png", units = 'in', res = 300, width=20, height=10)
ggarrange(A, B, ncol=2, nrow=1, common.legend = F, legend="right")
#dev.off()

###Analysis by time points
###Make Phloseq subsets to run the analysis
##V1V2
PS2.sputum12<- subset_samples(PS2.sputum, Visit%in%c("V1", "V2"))
##V2V3
PS2.sputum23<- subset_samples(PS2.sputum, Visit%in%c("V2", "V3"))
##V1V3
PS2.sputum13<- subset_samples(PS2.sputum, Visit%in%c("V1", "V3"))

deseq.visit<- phyloseq_to_deseq2(PS2.sputum13, ~ Visit)

geoMeans<- apply(counts(deseq.visit), 1, gm_mean)
deseq.visit<- estimateSizeFactors(deseq.visit, geoMeans = geoMeans)
deseq.visit<- DESeq(deseq.visit, fitType="local")

ac.res <- results(deseq.visit)
##Remove NA
ac.res <- ac.res[order(ac.res$padj, na.last=NA), ]

##Select cut-off value
alpha <- 0.05
Bac.sigtab <- ac.res[(ac.res$padj < alpha), ]
##Adjust based on the visit combination
Bac.sigtab <- cbind(as(ac.res, "data.frame"), as(tax_table(PS2.sputum13)[rownames(ac.res), ], "matrix")) ##Change based on the visit

##Adjust value
Bac.posigtab <- Bac.sigtab[Bac.sigtab[, "log2FoldChange"] > -30, ]
Bac.posigtab <- Bac.posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

#write.csv(Bac.posigtab, "~/CF_project/exercise-cf-intervention/tables/Q5_DeSeq2_Abund_Sputum_V1V2.csv") 
#write.csv(Bac.posigtab, "~/CF_project/exercise-cf-intervention/tables/Q5_DeSeq2_Abund_Sputum_V2V3.csv") 
#write.csv(Bac.posigtab, "~/CF_project/exercise-cf-intervention/tables/Q5_DeSeq2_Abund_Sputum_V1V3.csv") 

### IN sputum nothing was significant between days by DeSeq2