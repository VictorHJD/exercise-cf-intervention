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

dom.palette<- c( "Streptococcus"= "#925E9FFF", "Staphylococcus"= "#008B45FF", 
                 "Stenotrophomonas"= "#B09C85FF", "Pseudomonas" = "#ED0000FF")
######Sputum###################
sdt%>%
  dplyr::filter(Benzoase==1)%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))-> sdt.sputum

###Diversity and lung function 
sdt%>%
  dplyr::filter(material=="Sputum")%>%
  dplyr::filter(Benzoase==1)%>%
  mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  ggplot(aes(x= diversity_shannon, y= ppFEV1, shape= Visit))+
  geom_point(size=3, aes(fill= Patient_number), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  geom_smooth(method=lm, se = F, color= "black")+
  scale_fill_manual(values = pal.CF)+
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=70),alpha=0.05,fill="grey")+
  geom_hline(yintercept=70, linetype="dashed", color = "red")+
  theme_bw()+
  labs(tag= "A)")+
  xlab("Alpha diveristy (Shannon Index)")+
  ylab("Lung function (ppFEV1)")+
  stat_cor(method = "spearman", label.x = 0, label.y = 85)+
  theme(text = element_text(size=16), legend.position = "none")+
  facet_grid(rows = vars(Visit))-> A

sdt%>%
  dplyr::filter(material=="Sputum")%>%
  dplyr::filter(Benzoase==1)%>%
  mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  ggplot(aes(x= diversity_shannon, y= ppFEV1))+
  geom_point(size=3, aes(fill= Patient_number, shape= Visit), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  geom_smooth(method=lm, se = F, color= "black")+
  scale_fill_manual(values = pal.CF)+
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=70),alpha=0.01,fill="grey")+
  geom_hline(yintercept=70, linetype="dashed", color = "red")+
  theme_bw()+
  labs(tag= "B)")+
  labs(fill = "Patient")+
  labs(shape = "Visit")+
  guides(fill = guide_legend(override.aes=list(shape=c(21)), ncol = 6), shape= guide_legend(nrow = 3))+
  xlab("Alpha diveristy (Shannon Index)")+
  ylab("Lung function (ppFEV1)")+
  stat_cor(method = "spearman", label.x = 3, label.y = 30)+ # Add sperman`s correlation coefficient
  theme(text = element_text(size=16), legend.position="bottom", legend.box = "horizontal")-> B

A+
  xlab(NULL)-> A

C<- grid.arrange(A,B, heights = c(3, 2))

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Alpha_Lung_Sputum.pdf", plot = C, width = 8, height = 10, dpi = 600)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Alpha_Lung_Sputum.png", plot = C, width = 8, height = 10, dpi = 600)

rm(A,B,C)

lf.model.sputum1 <- lm(ppFEV1 ~ diversity_shannon, data= sdt.sputum)
lf.model.sputum2 <- lm(ppFEV1 ~ diversity_shannon * Visit, data= sdt.sputum)

car::Anova(lf.model.sputum1, type=3) 
car::Anova(lf.model.sputum2, type=3) 

lf.model.sputum.lsm <-
  lsmeans::lsmeans(lf.model.sputum2,
                   pairwise~diversity_shannon:Visit,
                   adjust="fdr")

lf.model.sputum.lsm$contrasts

###Dominance and lung function 
sdt%>%
  dplyr::filter(material=="Sputum")%>%
  dplyr::filter(Benzoase==1)%>%
  mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= dominance_dbp, y= ppFEV1, shape= Visit))+
  geom_point(size=3, aes(fill= Patient_number), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  geom_smooth(method=lm, se = F, color= "black")+
  scale_fill_manual(values = pal.CF)+
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=70),alpha=0.05,fill="grey")+
  geom_hline(yintercept=70, linetype="dashed", color = "red")+
  theme_bw()+
  labs(tag= "A)")+
  xlab("Dominance (Berger-Parker Index)")+
  ylab("Lung function (ppFEV1)")+
  stat_cor(method = "spearman", label.x = 0, label.y = 85)+ # Add sperman`s correlation coefficient
  theme(text = element_text(size=16), legend.position = "none")+
  facet_grid(rows = vars(Visit))-> A

sdt%>%
  dplyr::filter(material=="Sputum")%>%
  dplyr::filter(Benzoase==1)%>%
  mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= dominance_dbp, y= ppFEV1))+
  geom_point(size=3, aes(fill= Patient_number, shape= Visit), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  geom_smooth(method=lm, se = F, color= "black")+
  scale_fill_manual(values = pal.CF)+
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=70),alpha=0.01,fill="grey")+
  geom_hline(yintercept=70, linetype="dashed", color = "red")+
  theme_bw()+
  labs(tag= "B)")+
  labs(fill = "Patient")+
  labs(shape = "Visit")+
  guides(fill = guide_legend(override.aes=list(shape=c(21)), ncol = 6), shape= guide_legend(nrow = 3))+
  xlab("Dominance (Berger-Parker Index)")+
  ylab("Lung function (ppFEV1)")+
  stat_cor(method = "spearman", label.x = 0, label.y = 30)+ # Add sperman`s correlation coefficient
  theme(text = element_text(size=16), legend.position="bottom", legend.box = "horizontal")-> B

C<- grid.arrange(A,B, heights = c(3, 2))

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Dominance_Lung_Sputuml.pdf", plot = C, width = 10, height = 10)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Dominance_Lung_Sputum.png", plot = C, width = 10, height = 10)

rm(A,B,C)

##Barplot by sample 
gen.sputum<- count.high.genus(x = PS4.sput.Gen, num = 10)

sdt.sputum%>%
  rownames_to_column()%>%
  dplyr::select(c(1, 32:33))%>%
  dplyr::rename(SampleID = rowname)%>%
  left_join(gen.sputum, by="SampleID")-> gen.sputum

#set color palette to accommodate the number of genera
length(unique(gen.sputum$Genus))

#plot
gen.sputum%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x=Patient_number, y=Abundance, fill=Genus))+ 
  geom_bar(aes(), stat="identity", position="stack") + 
  facet_wrap(~Visit, scales= "free_x", nrow=1)+
  scale_fill_manual(values=tax.palette) + 
  guides(fill=guide_legend(nrow=4))+
  theme_bw()+
  labs(tag= "B)")+
  ylab("Relative abundance (%)")+
  xlab("Patient number")+
  theme(legend.position="bottom")# -> B

#C<- grid.arrange(A,B, heights = c(3, 3))

###Colonizers and satellites
Type.sputum<- count.genus(x = PS4.sput.Gen, num = 4.5)

sdt.sputum%>%
  rownames_to_column()%>%
  dplyr::select(c(1, 32:33))%>%
  dplyr::rename(SampleID = rowname)%>%
  left_join(Type.sputum, by="SampleID")-> Type.sputum

saveRDS(Type.sputum, "CF_project/exercise-cf-intervention/data/Type.sputum.rds")

Type.sputum%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x=Patient_number, y=Abundance, fill=Type))+ 
  geom_bar(aes(), stat="identity", position="stack") + 
  facet_wrap(~Visit, scales= "free_x", nrow=1)+
  guides(fill=guide_legend(nrow=5))+
  theme_bw()+
  labs(tag= "A)")+
  ylab("Relative abundance (%)")+
  xlab("Patient number")

#ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Composition_Genus.pdf", plot = C, width = 10, height = 12)
#ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Composition_Genus.png", plot = C, width = 10, height = 12)

#rm(A,B,C)
##Comparison V1 vs V3 (Individual look)
sdt%>%
  dplyr::filter(material=="Sputum")%>%
  dplyr::filter(Benzoase==1)%>%
  mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  dplyr::filter(Visit!="V2")%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= dominance_dbp, y= ppFEV1))+
  geom_point(size=3, aes(fill= Patient_number, shape= Visit), color= "black")+
  scale_shape_manual(values = c(21, 24))+ 
  scale_fill_manual(values = pal.CF)+
  geom_smooth(method=lm, se = F, aes(color= Patient_number))+
  scale_color_manual(values = pal.CF)+
  theme_bw()+
  labs(fill = "Patient")+
  labs(shape = "Visit")+
  guides(fill = guide_legend(override.aes=list(shape=c(21)), ncol = 6), shape= guide_legend(nrow = 3), color= "none")+
  xlab("Dominance (Berger-Parker Index)")+
  ylab("Lung function (ppFEV1)")+
  theme(text = element_text(size=16), legend.position="bottom", legend.box = "horizontal")

##Get top taxa per patient
top.sputum<- find.top.asv(PS4.sput, "Genus", 1)
top.sputum$Species<- NULL
##Add sample data
sdt.sputum%>%
  rownames_to_column()%>%
  dplyr::select(c(1, 32:95))%>%
  dplyr::rename(SampleID = rowname)%>%
  left_join(top.sputum, by="SampleID")-> top.sputum

top.sputum %>% 
  count(Genus)-> tmp

top.sputum%>%
  left_join(tmp, by="Genus")-> top.sputum

##Plot Lung function and dominant bugs
top.sputum%>%
  mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(y= reorder(Genus, n), x= ppFEV1))+
  geom_point(shape= 21, aes(fill= Patient_number, size = Abundance), color= "black", alpha= 0.75)+
  scale_fill_manual(values = pal.CF)+
  scale_size(range = c(.1, 10))+
  theme_bw()+
  labs(fill = "Patient", size = "Rel. abund (%)", tag = "A)")+
  labs()+
  guides(fill = guide_legend(override.aes=list(shape=c(21)), ncol = 6, size= 10), size= guide_legend(nrow = 2), color= "none")+
  facet_wrap(~Visit, scales= "free_x", nrow=1)+
  scale_x_continuous(breaks = c(20, 40, 60, 80, 100),limits = c(20, 100))+
  ylab("Dominant bacterial genus")+
  xlab("Lung function (ppFEV1)")+
  geom_vline(xintercept=70, linetype="dashed", color = "red")+
  theme(text = element_text(size=16), legend.position="bottom", legend.box = "horizontal",
        axis.text.y = element_text(size = 9, face="italic", color="black")) -> A

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Dominant_Lungfunct_Sputum.pdf", plot = A, width = 10, height = 8, dpi = 600)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Dominant_Lungfunct_Sputum.png", plot = A, width = 10, height = 8, dpi = 600)

##Make a heatmap
tmp<- as.data.frame(otu_table(tax_glom(PS4.sput, "Genus")))

tmp1<- count.genus(PS4.sput.Gen, 1)

tmp1%>%
  ungroup()%>%
  dplyr::select(c(Genus, ASV))%>%
  unique()-> tmp1

#tmp<- tmp[rownames(tmp) %in% tmp1$ASV, ] ##Subset just dominant genus for all samples

tmp<- Rel.abund_fun(tmp) ##Transform into relative abundance 

tmp%>%
  rownames_to_column()%>%
  dplyr::rename(ASV= rowname)%>%
  left_join(tmp1, by="ASV")%>%
  unite(ASV_Genus, c("ASV", "Genus"))%>%
  column_to_rownames(var = "ASV_Genus")->tmp
##Eliminate rows with zero
tmp<- tmp[apply(tmp[,-1], 1, function(x) !all(x<0.9)),]

library("pheatmap")
library(dendextend)
###In order to add the annotations in good order, 
#it is necessary to have the same order in the intersection matrix and in the annotation table
P.clust <- hclust(dist(t(tmp)), method = "ward.D") ##Dendogram

as.dendrogram(P.clust) %>%
  plot(horiz = T)

P.col <- cutree(tree = P.clust, k = 2)
P.col  <- data.frame(cluster = ifelse(test = P.col  == 1, yes = "cluster 1", no = "cluster 2"))
P.col$SampleID <- rownames(P.col)
P.col$SampleID<- NULL
P.col<- cbind(P.col, top.sputum)

col_groups <- P.col %>%
  mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  dplyr::select(c(SampleID, Patient_number,Genus)) ##Here It is possible to add the other characteristics

row.names(col_groups)<- col_groups$SampleID

col_groups$SampleID<- NULL

colour_groups <- list(Patient_number= pal.CF, Genus= dom.palette)

sput.heatmap <- pheatmap(tmp, cluster_rows = F, cluster_cols = T,
                          color = colorRampPalette(c("white","#832424FF"))(50), #"#3A3A98FF",
                          border_color = NA,
                          annotation_col = col_groups, 
                          annotation_colors = colour_groups,
                          show_rownames = T,
                          show_colnames = F)

ggsave(file = "CF_project/exercise-cf-intervention/figures/Heatmap_Sputum.pdf", plot = sput.heatmap, width = 10, height = 8, dpi = 600)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Heatmap_Sputum.png", plot = sput.heatmap, width = 10, height = 8, dpi = 600)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Heatmap_Sputum.svg", plot = sput.heatmap, width = 10, height = 8, dpi = 600)

##Bray-Curtis
BC_dist<- phyloseq::distance(PS4.sput,
                             method="bray", weighted=F)
ordination<- ordinate(PS4.sput,
                      method="PCoA", distance= BC_dist)

plot_ordination(PS4.sput, ordination, shape= "Visit")+ 
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(fill= Patient_number, shape= Visit), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+
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

###Phenotype severity
plot_ordination(PS4.sput, ordination)+ 
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(fill= Patient_number, shape= as.factor(Phenotype_severity)), color= "black")+
  scale_shape_manual(values = c(25, 24), labels = c("Low severity", "High severity"))+
  scale_fill_manual(values = pal.CF)+
  labs(title = "Bray-Curtis dissimilariy sputum",tag= "B)")+
  stat_ellipse(aes(color = as.factor(Phenotype_severity)))+
  scale_color_manual(values=c("#8A9045FF", "#800000FF"))+
  theme_bw()+
  theme(text = element_text(size=16))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= F)+
  labs(fill = "Patient")+
  labs(shape = "Phenotype severity")+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> PCo.Sev.Sputum

PCo.Sev.Stool<- readRDS("CF_project/exercise-cf-intervention/data/PCo.Sev.Stool.rds")

plot<-ggarrange(PCo.Sev.Stool, PCo.Sev.Sputum, ncol=1, nrow=2, common.legend = TRUE, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Severity.pdf", plot = plot, width = 10, height = 8)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Severity.png", plot = plot, width = 10, height = 8)

rm(PCo.Sev.Sputum, PCo.Sev.Stool, plot)

##Stratified for Patient number 
BC.test.sputum<- vegan::adonis(BC_dist~ Phenotype_severity+ Mutation_severity + sex + age +  Visit + BMI,
                               permutations = 999, data = sdt.sputum, na.action = F, strata = sdt.sputum$Patient_number)

kable(BC.test.sputum$aov.tab)
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

###Add Antibiotic intake information 

antibiotic<- read.csv("~/CF_project/Metadata/sample_data_indexed_antibiotics.csv")

antibiotic%>%
  dplyr::select(c(Comed_token, Number_antibioticCourses_priorstudystart, Number_antibioticCourses_duringstudy,
                  Number_iv_courses_priorstudy, Number_iv_courses_duringstudy))%>%
  dplyr::mutate(Comed_token= gsub("^(.*)V", "\\1_V", Comed_token))%>%
  separate(Comed_token, c("Patient_number", "Visit"))%>%
  dplyr::select(c(Patient_number, Number_antibioticCourses_priorstudystart, Number_antibioticCourses_duringstudy,
                  Number_iv_courses_priorstudy, Number_iv_courses_duringstudy))%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  unique()-> tmp2

BC_dist.sputum%>%
  left_join(tmp2, by="Patient_number")-> BC_dist.sputum

##Add total antibiotic burden
antibioticB<- read.csv("~/CF_project/Metadata/antibioticBurden.csv")

antibioticB%>%
  dplyr::select(c(Comed_token, Antibiotic_until15daysbefore, AntibioticBurden_total, AntibioticBurden_iv))%>%
  dplyr::distinct()%>%
  #There is a duplicated P7V1, ask Rebecca which one is the correct one, for now eliminate one line 38 seems to be the problem
  dplyr::slice(-38)%>%
  dplyr::mutate(Comed_token= gsub("^(.*)V", "\\1_V", Comed_token))%>%
  separate(Comed_token, c("Patient_number", "Visit"))%>%
  group_by(Patient_number)%>%
  dplyr::select(c(Patient_number, Visit, AntibioticBurden_total))%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  dplyr::mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  dplyr::distinct()%>%
  spread(Visit, AntibioticBurden_total)%>%
  dplyr::mutate(V1_V2= V1 - V2)%>%
  dplyr::mutate(V1_V3= V1 - V3)%>%
  dplyr::mutate(V2_V3= V2 - V3)%>%
  dplyr::select(c(Patient_number, V1_V2, V2_V3, V1_V3))%>%
  pivot_longer(!Patient_number, names_to = "Group", values_to = "AntibioticBurden_total")%>%
  dplyr::mutate(ID= paste0(Patient_number, Group))%>%
  dplyr::ungroup()%>%
  dplyr::select(c(ID, AntibioticBurden_total))-> tmp2

BC_dist.sputum%>%
  left_join(tmp2, by="ID")-> BC_dist.sputum

##Add total IV antibiotic burden
antibioticB%>%
  dplyr::select(c(Comed_token, Antibiotic_until15daysbefore, AntibioticBurden_total, AntibioticBurden_iv))%>%
  dplyr::distinct()%>%
  #There is a duplicated P7V1, ask Rebecca which one is the correct one, for now eliminate one line 38 seems to be the problem
  dplyr::slice(-38)%>%
  dplyr::mutate(Comed_token= gsub("^(.*)V", "\\1_V", Comed_token))%>%
  separate(Comed_token, c("Patient_number", "Visit"))%>%
  group_by(Patient_number)%>%
  dplyr::select(c(Patient_number, Visit, AntibioticBurden_iv))%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  dplyr::mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  dplyr::distinct()%>%
  spread(Visit, AntibioticBurden_iv)%>%
  dplyr::mutate(V1_V2= V1 - V2)%>%
  dplyr::mutate(V1_V3= V1 - V3)%>%
  dplyr::mutate(V2_V3= V2 - V3)%>%
  dplyr::select(c(Patient_number, V1_V2, V2_V3, V1_V3))%>%
  pivot_longer(!Patient_number, names_to = "Group", values_to = "AntibioticBurden_iv")%>%
  dplyr::mutate(ID= paste0(Patient_number, Group))%>%
  dplyr::ungroup()%>%
  dplyr::select(c(ID, AntibioticBurden_iv))-> tmp2

BC_dist.sputum%>%
  left_join(tmp2, by="ID")-> BC_dist.sputum

##Add a time between visits (Overall for know but ask values per patient per period)
BC_dist.sputum%>%
  dplyr::mutate(Months= case_when(Group == "V1_V3" ~ 12,
                                  Group == "V1_V2" ~ 3,
                                  Group == "V2_V3" ~ 19))-> BC_dist.sputum

saveRDS(BC_dist.sputum, "~/CF_project/exercise-cf-intervention/data/BC_dist.sputum.rds")

##Is antibiotic burden differences impacting in composition by patient? 
### Linear model test
require("lmtest")
require("lme4")
test<- BC_dist.sputum[!is.na(BC_dist.sputum$AntibioticBurden_total),]

print(summary (lmer (data = test, rank (BC_dist) ~  AntibioticBurden_total + 
                       (1 | Patient_number) + (1 | Group), REML = F)))

##Nested model for AntibioticBurden_total
pABXburden<- lrtest (lmer (data = test, rank (BC_dist) ~ AntibioticBurden_total + (1 | Patient_number) + (1 | Group), REML = F),
                     lmer (data = test, rank (BC_dist) ~ (1 | Patient_number) + (1 | Group), REML = F))$'Pr(>Chisq)' [2]

##How much variance is explained by each?
mm.ABX <- lmer (data = test, rank (BC_dist) ~ AntibioticBurden_total  + (1 | Patient_number) + (1 | Group), REML = F)
varianceTable <- as.data.frame(anova (mm.ABX))
varianceTable$VarExplained <- varianceTable$`Sum Sq` / sum (resid (mm.ABX)^2)
varianceTable$Variable <- rownames(varianceTable)
varianceTable[2, ] <- c(rep(1), (1 - sum(varianceTable$VarExplained)), "Residuals")
varianceTable$VarExplained <- as.numeric(varianceTable$VarExplained)
varianceTable$VarLabels <- scales::percent(varianceTable$VarExplained)
print(varianceTable)

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
  scale_x_discrete(labels=c("V1_V2" =  "V1 to V2", 
                            "V2_V3" = "V2 to V3",
                            "V1_V3" = "V1 to V3"))+
  stat_pvalue_manual(stats.test, hide.ns = F, step.increase = 0.05,
                     tip.length = 0, label = "{p.adj} {p.adj.signif}")->B

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

##Correlation with Antibiotic intake 
BC_dist.sputum%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Number_antibioticCourses_priorstudystart, y= BC_dist, color= Group))+
  geom_point(size=2.5, aes(shape= Group, fill= Patient_number), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  xlab("Number of antibiotic courses prior study start")+
  ylab("Bray-Curtis dissimilarity (Sputum microbiome)")+
  labs(tag= "A)")+
  theme_classic()+
  geom_smooth(method=lm, se = F)+
  stat_regline_equation(label.x = 2, label.y = 1)+ 
  stat_cor(label.x = 2, label.y = 0.9, 
           aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+
  scale_fill_manual(values = pal.CF)+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= "none")+
  labs(fill = "Patient")+
  labs(shape = "Visit period")+
  facet_wrap(~Group)+
  theme(text = element_text(size=16))-> A

BC_dist.sputum%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Number_antibioticCourses_duringstudy, y= BC_dist, color= Group))+
  geom_point(size=2.5, aes(shape= Group, fill= Patient_number), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  xlab("Number of antibiotic courses during study start")+
  ylab("Bray-Curtis dissimilarity (Sputum microbiome)")+
  labs(tag= "B)")+
  theme_classic()+
  geom_smooth(method=lm, se = F)+
  stat_regline_equation(label.x = 2, label.y = 1)+ 
  stat_cor(label.x = 2, label.y = 0.9, 
           aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+
  scale_fill_manual(values = pal.CF)+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= "none")+
  labs(fill = "Patient")+
  labs(shape = "Visit period")+
  facet_wrap(~Group)+
  theme(text = element_text(size=16))-> B

BC_dist.sputum%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Number_iv_courses_priorstudy, y= BC_dist, color= Group))+
  geom_point(size=2.5, aes(shape= Group, fill= Patient_number), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  xlab("Number of iv antibiotic courses prior study start")+
  ylab("Bray-Curtis dissimilarity (Sputum microbiome)")+
  labs(tag= "C)")+
  theme_classic()+
  geom_smooth(method=lm, se = F)+
  stat_regline_equation(label.x = 0.5, label.y = 1)+ 
  stat_cor(label.x = 0.5, label.y = 0.9, 
           aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+
  scale_fill_manual(values = pal.CF)+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= "none")+
  labs(fill = "Patient")+
  labs(shape = "Visit period")+
  facet_wrap(~Group)+
  theme(text = element_text(size=16))-> C

BC_dist.sputum%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Number_iv_courses_duringstudy, y= BC_dist, color= Group))+
  geom_point(size=2.5, aes(shape= Group, fill= Patient_number), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  xlab("Number of iv antibiotic courses during study")+
  ylab("Bray-Curtis dissimilarity (Sputum microbiome)")+
  labs(tag= "D)")+
  theme_classic()+
  geom_smooth(method=lm, se = F)+
  stat_regline_equation(label.x = 0.5, label.y = 0.8)+ 
  stat_cor(label.x = 0.5, label.y = 0.75, 
           aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+
  scale_fill_manual(values = pal.CF)+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= "none")+
  labs(fill = "Patient")+
  labs(shape = "Visit period")+
  facet_wrap(~Group)+
  theme(text = element_text(size=16))-> D

plot1<-ggarrange(A, B, ncol=1, nrow=2, common.legend = TRUE, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Sputum_Oral_Antibiotics.pdf", plot = plot1, width = 10, height = 12)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Sputum_Oral_Antibiotics.png", plot = plot1, width = 10, height = 12)

plot2<-ggarrange(C, D, ncol=1, nrow=2, common.legend = TRUE, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Sputum_iv_Antibiotics.pdf", plot = plot2, width = 10, height = 12)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Sputum_iv_Antibiotics.png", plot = plot2, width = 10, height = 12)

rm(A,B, C, D, plot1, plot2)

##
BC_dist.sputum%>%
  dplyr::filter(Group=="V1_V3")%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Number_antibioticCourses_duringstudy, y= BC_dist))+
  geom_point(size=2.5, shape= 21, aes(fill= Patient_number), color= "black")+
  xlab("Number of antibiotic courses during study")+
  ylab("Bray-Curtis dissimilarity V1 - V3 \n (Sputum microbiome)")+
  labs(tag= "A)")+
  theme_classic()+
  geom_smooth(method=lm, se = T, color= "black")+
  stat_regline_equation(label.x = 3, label.y = 0.9)+ 
  stat_cor(label.x = 3, label.y = 0.8, 
           aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+
  scale_fill_manual(values = pal.CF)+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(fill = "Patient")+
  labs(shape = "Visit period")+
  theme(text = element_text(size=16))-> plot1

BC_dist.sputum%>%
  dplyr::filter(Group=="V1_V3")%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Number_iv_courses_duringstudy, y= BC_dist))+
  geom_point(size=2.5, shape= 21, aes(fill= Patient_number), color= "black")+
  xlab("Number of iv antibiotic courses during study")+
  ylab("Bray-Curtis dissimilarity V1 - V3 \n (Sputum microbiome)")+
  labs(tag= "B)")+
  theme_classic()+
  geom_smooth(method=lm, se = T, color= "black")+
  stat_regline_equation(label.x = 0.5, label.y = 1.1)+ 
  stat_cor(label.x = 0.5, label.y = 1, 
           aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+
  scale_fill_manual(values = pal.CF)+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(fill = "Patient")+
  labs(shape = "Visit period")+
  theme(text = element_text(size=16))-> plot2

plot1<-ggarrange(plot1, plot2, ncol=1, nrow=2, common.legend = TRUE, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Sputum_V1V3_Antibiotics.pdf", plot = plot1, width = 10, height = 12)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Sputum_V1V3_Antibiotics.png", plot = plot1, width = 10, height = 12)

###Matrix for Mantel test

BC_sputum<- as.matrix(BC_dist)

##Change row and colnames for a uniform system between sample type Patient:Visit
rownames(BC_sputum)<-  gsub("10P", "P", rownames(BC_sputum))
rownames(BC_sputum)<-  gsub("A", "\\1", rownames(BC_sputum))

colnames(BC_sputum)<-  gsub("10P", "P", colnames(BC_sputum))
colnames(BC_sputum)<-  gsub("A", "\\1", colnames(BC_sputum))

saveRDS(BC_sputum, "~/CF_project/exercise-cf-intervention/data/BC_sputum_matrix.rds")

###Mixed effect models 
BC_dist.sputum%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))-> BC_dist.sputum
BC_dist.sputum%>%
  dplyr::select(ppFVC, Trainingfrequency, Trainingtime, ppFEV1, BC_dist, Patient_number, Group, 
                Number_antibioticCourses_priorstudystart, Number_antibioticCourses_duringstudy, 
                Number_iv_courses_priorstudy, Number_iv_courses_duringstudy)%>%
  dplyr::filter(complete.cases(.))-> tmp

colnames(tmp)<- c("ppFVC", "Trainingfrequency", "Trainingtime", "ppFEV1", "BC_dist", "Patient_number", "Group", 
                  "ab_prior", "ab_during", "iv_prior", "iv_during")
##qqPlots
qqPlot(BC_dist.sputum$BC_dist)
qqPlot(BC_dist.sputum$Trainingfrequency)
qqPlot(BC_dist.sputum$Trainingtime)
qqPlot(BC_dist.sputum$ppFEV1)
qqPlot(BC_dist.sputum$ppFVC)

##Model selection do it with glm
full.model<- glm(BC_dist ~ Trainingfrequency*Trainingtime*ppFEV1*ppFVC*ab_prior*ab_during*iv_prior*iv_during, data = tmp) ##Full model
# Stepwise regression model
step.model <- MASS::stepAIC(full.model, direction = "both", 
                            trace = FALSE)
summary(step.model)
step.model.df<- as.data.frame(coef(summary(step.model)))

step.model.df%>%
  mutate(p.adj = p.adjust(`Pr(>|t|)`, method='BH')) %>%
  add_significance()%>%
  rownames_to_column()-> step.model.df

write.csv(step.model.df, "~/CF_project/exercise-cf-intervention/tables/Q2_BC_Sports_Lung_Sputum.csv", row.names = F)

##Antibiotics
full.model<- glm(BC_dist ~ ab_prior*ab_during*iv_prior*iv_during, data = tmp) ##Full model
# Stepwise regression model
step.model <- MASS::stepAIC(full.model, direction = "both", 
                            trace = FALSE)
summary(step.model)
step.model.df<- as.data.frame(coef(summary(step.model)))

step.model.df%>%
  mutate(p.adj = p.adjust(`Pr(>|t|)`, method='BH')) %>%
  add_significance()%>%
  rownames_to_column()-> step.model.df

write.csv(step.model.df, "~/CF_project/exercise-cf-intervention/tables/Q2_BC_Sputum_antibiotics.csv", row.names = F)

##Mixed effect models
##with patient and intervisit group as individual random effects

tr0<- lmer(BC_dist ~ (1 | Patient_number), data = tmp) ##Null model
tr1<-lmer(BC_dist ~ Trainingfrequency*Trainingtime + (1 | Patient_number), data = tmp)
tr2<-lmer(BC_dist ~ ppFEV1*ppFVC + (1 | Patient_number), data = tmp)

lrtest(tr1, tr2)
##Each factor add predictive value to the model 
##What happen with interactions 
tr3<-lmer(BC_dist ~ Trainingfrequency*Trainingtime*ppFEV1*ppFVC + (1 | Patient_number), data = BC_dist.sputum)
summary(tr3) ##Full model

lrtest(tr3, tr1)
lrtest(tr3, tr2)

A<- plotREsim(REsim(tr3))  ## plot the interval estimates
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
B<- plot_model(tr3, p.adjust = "BH", vline.color = "gray",
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

C<-ggarrange(B, A, ncol=1, nrow=2)

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Sputum_Training_Effect_Ranges.png", plot = C, width = 8, height = 8)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Sputum_Training_Effect_Ranges.pdf", plot = C, width = 8, height = 8)

rm(A, B, C)

##Test predictive value of BC differences within patient to lung function measurements Delta ppFEV1/ppFVC between visits
##Model selection do it with glm
full.model<- glm(ppFEV1 ~ Trainingfrequency*Trainingtime*BC_dist, data = tmp) ##Full model
# Stepwise regression model
step.model <- MASS::stepAIC(full.model, direction = "both", 
                            trace = FALSE)
summary(step.model)

##Model selection do it with glm
full.model<- glm(ppFVC ~ Trainingfrequency*Trainingtime*BC_dist, data = tmp) ##Full model
# Stepwise regression model
step.model <- MASS::stepAIC(full.model, direction = "both", 
                            trace = FALSE)
summary(step.model)

step.model.df<- as.data.frame(coef(summary(step.model))) ##The full model is the best!!!

step.model.df%>%
  mutate(p.adj = p.adjust(`Pr(>|t|)`, method='BH')) %>%
  add_significance()%>%
  rownames_to_column()-> step.model.df

write.csv(step.model.df, "~/CF_project/exercise-cf-intervention/tables/Q2_Lung_BC_Sports_Sputum.csv", row.names = F)

### Likelihood ratio test
tr1<- lmer(ppFVC ~ Trainingfrequency*Trainingtime + (1 | Patient_number), data = tmp) ##Null model
tr2<- lmer(ppFEV1 ~ Trainingfrequency + (1 | Patient_number), data = tmp) ##Null model
tr3<- lmer(ppFEV1 ~ Trainingfrequency + BC_dist + (1 | Patient_number), data = tmp) ##Best model
tr4<-lmer(ppFVC ~ Trainingfrequency*Trainingtime*BC_dist + (1 | Patient_number), data = tmp)
summary(tr4) ## Best model

lrtest(tr1, tr4)
lrtest(tr2, tr3)

###Does differences in bacterial composition within patient predict severity status 
BC_dist.sputum%>%
  dplyr::mutate(Phenotype_severity = case_when(Phenotype_severity == 2  ~ 1,
                                               Phenotype_severity == 1 ~ 0))%>%
  dplyr::mutate(Mutation_severity = case_when(Mutation_severity == 2  ~ 1,
                                              Mutation_severity == 1 ~ 0))-> BC_dist.sputum

##Logistic regression 
log.model1 <-glmer(Phenotype_severity ~Trainingfrequency + Trainingtime +
                     ppFVC + ppFEV1 + BC_dist + (1 | Patient_number), data = BC_dist.sputum, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)##Full model

summary(log.model1)$coef

log.model2 <- glmer(Phenotype_severity ~Trainingfrequency + Trainingtime +
                      ppFVC + ppFEV1 + (1 | Patient_number), data = BC_dist.sputum, family = binomial, 
                    control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)##Without BC dissimilarity model

summary(log.model2)$coef

lrtest(log.model1, log.model2)  ##Bray-Curtis dissimilarity do not add predictive power for phenotype prediction

BC_dist.sputum%>%
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
                     ppFVC + ppFEV1 + BC_dist + (1 | Patient_number), data = BC_dist.sputum, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)##Full model

summary(log.model3)$coef

log.model4 <- glmer(Sport_Response ~Trainingfrequency + Trainingtime +
                      ppFVC + ppFEV1 + (1 | Patient_number), data = BC_dist.sputum, family = binomial, 
                    control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)##Without BC dissimilarity model

summary(log.model4)$coef

lrtest(log.model3, log.model4) 

##Microbial composition do not add predictive power response to sports intervention 
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
  dplyr::select(Genus, Phylum)-> tax
tax$Genus<-gsub(" ", "_", basename(tax$Genus))
tax<- rownames_to_column(tax, var = "ASV")

##Use genus as rownames
sputum.microbiome<- plyr::join(otu, tax, by= "ASV")
sputum.microbiome$ASV<- NULL
rownames(sputum.microbiome)<- paste0(sputum.microbiome$Phylum, "-", sputum.microbiome$Genus)
sputum.microbiome$Genus<- NULL
sputum.microbiome$Phylum<- NULL

##Transpose dataframe so samples are rows 
sputum.microbiome<- t(sputum.microbiome)
#saveRDS(sputum.microbiome, "~/CF_project/exercise-cf-intervention/data/Sput_rare_ASV.rds") #--> For MetadeconfoundR

##Select useful metrics
y<-sdt.sputum

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
#saveRDS(y, "~/CF_project/exercise-cf-intervention/data/Sput_rare_Metadata.rds")#--> For MetadeconfoundR

# Cross correlate data sets, output is also available in a handy table format
###Let's run this with metadeconfoundR (other script for that)

###Deseq2 analysis
##Get raw data to run this
PS3.sputum<- subset_samples(PS3, Benzoase%in%c(1))
PS3.sputum<- tax_glom(PS3.sputum, "Genus")

##Adjustment make phenotype and genotype as factor 
PS3.sputum@sam_data$Phenotype_severity <- as.factor(PS3.sputum@sam_data$Phenotype_severity)
PS3.sputum@sam_data$Mutation_severity <- as.factor(PS3.sputum@sam_data$Mutation_severity)
##Severity classification based on:
#genotype
deseq.severity<- phyloseq_to_deseq2(PS3.sputum, ~ Mutation_severity)

# calculate geometric means prior to estimate size factors
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans<- apply(counts(deseq.severity), 1, gm_mean)
deseq.severity<- estimateSizeFactors(deseq.severity, geoMeans = geoMeans)
deseq.severity<- DESeq(deseq.severity, fitType="local")

ac.res <- results(deseq.severity)

##Select cut-off value
sigtab <- cbind(as(ac.res,"data.frame"), as(tax_table(PS3.sputum)[rownames(ac.res), ], "matrix"))
sigtab$Species<- NULL
head(sigtab,20)

##Volcano plot to detect differential taxa in sputum microbiome between severity 
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
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(shape=21, size=3, alpha= 0.5, aes(fill= AbundLev), color= "black")+
  ggrepel::geom_text_repel(aes(col=AbundLev, label=Genus)) +
  scale_fill_manual(values=c("#8A9045FF", "#767676FF")) +
  scale_color_manual(values=c("#8A9045FF", "#767676FF")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.001), col="black", linetype= "dashed") +
  labs(tag= "A)", x= "log2 Fold change", y= "-Log10 (p Adjusted)", fill= "Genus\nabundance")+
  theme_bw()+
  theme(text = element_text(size=16))+
  guides(color= F)-> A

##Save table 
#write.csv(sigtab, "~/CF_project/exercise-cf-intervention/tables/Q5_DeSeq2_Abund_sputum_Severity.csv") #Genotype
#write.csv(sigtab, "~/CF_project/exercise-cf-intervention/tables/Q5_DeSeq2_Abund_sputum_Gen_Severity.csv") #Phenotype

##Phenotype
deseq.severity<- phyloseq_to_deseq2(PS3.sputum, ~ Phenotype_severity)
# calculate geometric means prior to estimate size factors
geoMeans<- apply(counts(deseq.severity), 1, gm_mean)
deseq.severity<- estimateSizeFactors(deseq.severity, geoMeans = geoMeans)
deseq.severity<- DESeq(deseq.severity, fitType="local")

ac.res <- results(deseq.severity)

##Select cut-off value
sigtab <- cbind(as(ac.res,"data.frame"), as(tax_table(PS3.sputum)[rownames(ac.res), ], "matrix"))
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
  dplyr::mutate(Genus=  case_when(Genus == "TM7x" ~ "Nanosynbacter lyticus TM7x",
                                  TRUE ~ Genus))%>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(shape=21, size=3, alpha= 0.5, aes(fill= AbundLev), color= "black")+
  ggrepel::geom_text_repel(aes(col=AbundLev, label=Genus)) +
  scale_fill_manual(values=c("#767676FF")) +
  scale_color_manual(values=c("#767676FF")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.001), col="black", linetype= "dashed") +
  labs(tag= "B)", x= "log2 Fold change", y= "-Log10 (p Adjusted)", fill= "Genus\nabundance")+
  theme_bw()+
  theme(text = element_text(size=16))+
  guides(color= F)-> B

C<-ggarrange(A, B, ncol=1, nrow=2, common.legend = T, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q5_DefAbund_sputum_Severity.png", plot = C, width = 8, height = 8)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q5_DefAbund_sputum_Severity.pdf", plot = C, width = 8, height = 8)

rm(A, B, C)

###Analysis by time points
###Make Phloseq subsets to run the analysis
##Adjustment make phenotype and genotype as factor 
PS3.sputum@sam_data$Visit <- as.factor(PS3.sputum@sam_data$Visit)
##V1V2
PS3.sputum12<- subset_samples(PS3.sputum, Visit%in%c("V1", "V2"))
##V2V3
PS3.sputum23<- subset_samples(PS3.sputum, Visit%in%c("V2", "V3"))
##V1V3
PS3.sputum13<- subset_samples(PS3.sputum, Visit%in%c("V1", "V3"))

##Visit 2 vs 1
deseq.visit<- phyloseq_to_deseq2(PS3.sputum12, ~ Visit)

# calculate geometric means prior to estimate size factors
geoMeans<- apply(counts(deseq.visit), 1, gm_mean)
deseq.visit<- estimateSizeFactors(deseq.visit, geoMeans = geoMeans)
deseq.visit<- DESeq(deseq.visit, fitType="local")

ac.res <- results(deseq.visit)

##Select cut-off value
sigtab <- cbind(as(ac.res,"data.frame"), as(tax_table(PS3.sputum12)[rownames(ac.res), ], "matrix"))
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
  scale_fill_manual(values=c("#767676FF")) +
  scale_color_manual(values=c("#767676FF")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.001), col="black", linetype= "dashed") +
  labs(tag= "A)", x= "log2 Fold change", y= "-Log10 (p Adjusted)", fill= "Genus abundance")+
  theme_bw()+
  theme(text = element_text(size=16), legend.position = "top")+
  guides(color= F)-> A

A <- A + theme(legend.position="none")

##Visit 3 vs 2
deseq.visit<- phyloseq_to_deseq2(PS3.sputum23, ~ Visit)

# calculate geometric means prior to estimate size factors
geoMeans<- apply(counts(deseq.visit), 1, gm_mean)
deseq.visit<- estimateSizeFactors(deseq.visit, geoMeans = geoMeans)
deseq.visit<- DESeq(deseq.visit, fitType="local")

ac.res <- results(deseq.visit)

##Add taxa information
sigtab <- cbind(as(ac.res,"data.frame"), as(tax_table(PS3.sputum23)[rownames(ac.res), ], "matrix"))
sigtab$Species<- NULL
head(sigtab,20)

# add a column of Non-significant
sigtab$AbundLev <- "NS"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "High" 
sigtab$AbundLev[sigtab$log2FoldChange > 0.6 & sigtab$padj < 0.001] <- "High"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
sigtab$AbundLev[sigtab$log2FoldChange < -0.6 & sigtab$padj < 0.001] <- "Low"

#Organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
#plot adding up all layers we have seen so far Actinomyces sp. oral taxon 848 str. F0332
sigtab%>%
  dplyr::mutate(Genus=  case_when(Genus == "F0332" ~ "Actinomyces F0332",
                                  TRUE ~ Genus))%>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(shape=21, size=3, alpha= 0.5, aes(fill= AbundLev), color= "black")+
  ggrepel::geom_text_repel(aes(col=AbundLev, label=Genus)) +
  scale_fill_manual(values=c("#8A9045FF", "#767676FF")) +
  scale_color_manual(values=c("#8A9045FF","#767676FF")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.001), col="black", linetype= "dashed") +
  labs(tag= "B)", x= "log2 Fold change", y= "-Log10 (p Adjusted)", fill= "Genus abundance")+
  theme_bw()+
  theme(text = element_text(size=16), legend.position = "top")+
  guides(color= F)-> B

##Extract the legend from A to use it later as a common legend 
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(B)

##Remove legend from B
B <- B + theme(legend.position="none")

##Visit 3 vs 1
deseq.visit<- phyloseq_to_deseq2(PS3.sputum13, ~ Visit)

# calculate geometric means prior to estimate size factors
geoMeans<- apply(counts(deseq.visit), 1, gm_mean)
deseq.visit<- estimateSizeFactors(deseq.visit, geoMeans = geoMeans)
deseq.visit<- DESeq(deseq.visit, fitType="local")

ac.res <- results(deseq.visit)

##Add taxa information
sigtab <- cbind(as(ac.res,"data.frame"), as(tax_table(PS3.sputum13)[rownames(ac.res), ], "matrix"))
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
  dplyr::mutate(Genus=  case_when(Genus == "F0332" ~ "Actinomyces F0332",
                                  TRUE ~ Genus))%>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(shape=21, size=3, alpha= 0.5, aes(fill= AbundLev), color= "black")+
  ggrepel::geom_text_repel(aes(col=AbundLev, label=Genus)) +
  scale_fill_manual(values=c("#767676FF")) +
  scale_color_manual(values=c("#767676FF")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.001), col="black", linetype= "dashed") +
  labs(tag= "C)", x= "log2 Fold change", y= "-Log10 (p Adjusted)", fill= "Genus\nabundance")+
  theme_bw()+
  theme(text = element_text(size=16), legend.position = "none")+
  guides(color= F)-> C

D<- grid.arrange(legend, A,B,C, nrow=4, heights=c(0.5, 2.5, 2.5, 2.5))

##Save just when the three objects are in the environment
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q5_DefAbund_sputum_Visit.png", plot = D, width = 8, height = 8)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q5_DefAbund_sputum_Visit.pdf", plot = D, width = 8, height = 8)

rm(A,B,C, D)
##Response to sports classification:
#sports
PS3.sputum@sam_data$Sport_Response <- as.factor(PS3.sputum@sam_data$Sport_Response)
deseq.sports<- phyloseq_to_deseq2(PS3.sputum, ~ Sport_Response)

geoMeans<- apply(counts(deseq.sports), 1, gm_mean)
deseq.sports<- estimateSizeFactors(deseq.sports, geoMeans = geoMeans)
deseq.sports<- DESeq(deseq.sports, fitType="local")

ac.res <- results(deseq.sports)

##Select cut-off value
sigtab <- cbind(as(ac.res,"data.frame"), as(tax_table(PS3.sputum)[rownames(ac.res), ], "matrix"))
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
  scale_fill_manual(values=c("#767676FF")) +
  scale_color_manual(values=c("#767676FF")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.001), col="black", linetype= "dashed") +
  labs(tag= "B)", x= "log2 Fold change", y= "-Log10 (p Adjusted)", fill= "Genus\nabundance")+
  theme_bw()+
  theme(text = element_text(size=16))+
  guides(color= F)-> B

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q5_DefAbund_Sputum_Sports.png", plot = B, width = 8, height = 8)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q5_DefAbund_Sputum_Sports.pdf", plot = B, width = 8, height = 8)

rm(B)
##Correlations with responders for sputum and stool
##Lactobacillus 
tmp1<- as.data.frame(sputum.microbiome[,"Firmicutes-Lactobacillus"])
colnames(tmp1)<- c("Firmicutes-Lactobacillus")
sdt.sputum<- cbind(sdt.sputum, tmp1)

tmp1<- as.data.frame(stool.microbiome[,"Firmicutes-Lactobacillus"])
colnames(tmp1)<- c("Firmicutes-Lactobacillus")
sdt.stool<- cbind(sdt.stool, tmp1)

##Pseudomonas
tmp2<- as.data.frame(sputum.microbiome[,"Proteobacteria-Pseudomonas"])
colnames(tmp2)<- c("Proteobacteria-Pseudomonas")
sdt.sputum<- cbind(sdt.sputum, tmp2)

tmp2<- as.data.frame(stool.microbiome[,"Proteobacteria-Pseudomonas"])
colnames(tmp2)<- c("Proteobacteria-Pseudomonas")
sdt.stool<- cbind(sdt.stool, tmp2)

##Logistic regression 
###Stool
sdt.stool%>%
  dplyr::mutate(Phenotype_severity = case_when(Phenotype_severity == 2  ~ 1,
                                               Phenotype_severity == 1 ~ 0))%>%
  dplyr::mutate(Mutation_severity = case_when(Mutation_severity == 2  ~ 1,
                                              Mutation_severity == 1 ~ 0))-> sdt.stool

log.model1 <- glm(Phenotype_severity ~ `Proteobacteria-Pseudomonas`
                  , data = sdt.stool, family = binomial)
summary(log.model1)$coef

sdt.stool%>%
  ggplot(aes(`Proteobacteria-Pseudomonas`/1e6, Phenotype_severity)) +
  geom_point(size=3, aes(fill= Patient_number, shape= Visit), color= "black")+
  scale_shape_manual(values = c(21, 24, 22))+
  scale_fill_manual(values = pal.CF)+
  labs(x = "Psudomonas relative abundance",
       y = "Probability of severe CF phenotype", tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(fill = "Patient")+
  labs(shape = "Visit")+
  geom_smooth(method = "glm", method.args = list(family = "binomial"))-> A

sdt.stool%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  dplyr::group_by(Visit)%>%
  ggplot(aes(Visit, `Proteobacteria-Pseudomonas`/1e6)) +
  geom_line(aes(group = Patient_number), color= "gray")+
  geom_boxplot(color="black", alpha= 0.5)+
  scale_y_log10()+
  geom_point(size=3, aes(fill= Patient_number, shape= Visit), color= "black")+
  scale_shape_manual(values = c(21, 24, 22))+
  xlab("Visit")+
  ylab("Pseudomonas relative abundance")+
  labs(tag= "B)")+
  scale_fill_manual(values = pal.CF)+
  theme_classic()+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(fill = "Patient")+
  theme(text = element_text(size=16))-> B

sdt.stool%>%
  mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  dplyr::group_by(Visit)%>%
  ggplot(aes(Visit, `Firmicutes-Lactobacillus`/1e6)) +
  geom_line(aes(group = Patient_number), color= "gray")+
  geom_boxplot(color="black", alpha= 0.5)+
  scale_y_log10()+
  geom_point(size=3, aes(fill= Patient_number, shape= Visit), color= "black")+
  scale_shape_manual(values = c(21, 24, 22))+
  xlab("Visit")+
  ylab("Lactobacillus relative abundance")+
  labs(tag= "C)")+
  scale_fill_manual(values = pal.CF)+
  theme_classic()+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(fill = "Patient")+
  theme(text = element_text(size=16))-> C

###Sputum
sdt.sputum%>%
  dplyr::mutate(Phenotype_severity = case_when(Phenotype_severity == 2  ~ 1,
                                               Phenotype_severity == 1 ~ 0))%>%
  dplyr::mutate(Mutation_severity = case_when(Mutation_severity == 2  ~ 1,
                                              Mutation_severity == 1 ~ 0))-> sdt.sputum

log.model2 <- glm(Phenotype_severity ~ `Proteobacteria-Pseudomonas`
                    , data = sdt.sputum, family = binomial)
summary(log.model2)$coef

sdt.sputum%>%
  ggplot(aes(`Proteobacteria-Pseudomonas`/1e6, Phenotype_severity)) +
  geom_point(size=3, aes(fill= Patient_number, shape= Visit), color= "black")+
  scale_shape_manual(values = c(21, 24, 22))+
  scale_fill_manual(values = pal.CF)+
  labs(x = "Pseudomonas relative abundance",
       y = "Probability of severe CF phenotype", tag= "D)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(fill = "Patient")+
  labs(shape = "Visit")+
  geom_smooth(method = "glm", method.args = list(family = "binomial"))-> D

sdt.sputum%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  dplyr::group_by(Visit)%>%
  ggplot(aes(Visit, `Proteobacteria-Pseudomonas`/1e6)) +
  geom_line(aes(group = Patient_number), color= "gray")+
  geom_boxplot(color="black", alpha= 0.5)+
  scale_y_log10()+
  geom_point(size=3, aes(fill= Patient_number, shape= Visit), color= "black")+
  scale_shape_manual(values = c(21, 24, 22))+
  xlab("Visit")+
  ylab("Pseudomonas relative abundance")+
  labs(tag= "E)")+
  scale_fill_manual(values = pal.CF)+
  theme_classic()+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(fill = "Patient")+
  theme(text = element_text(size=16))-> E

sdt.sputum%>%
  mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  dplyr::group_by(Visit)%>%
  ggplot(aes(Visit, `Firmicutes-Lactobacillus`/1e6)) +
  geom_line(aes(group = Patient_number), color= "gray")+
  geom_boxplot(color="black", alpha= 0.5)+
  scale_y_log10()+
  geom_point(size=3, aes(fill= Patient_number, shape= Visit), color= "black")+
  scale_shape_manual(values = c(21, 24, 22))+
  xlab("Visit")+
  ylab("Lactobacillus relative abundance")+
  labs(tag= "F)")+
  scale_fill_manual(values = pal.CF)+
  theme_classic()+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(fill = "Patient")+
  theme(text = element_text(size=16))-> f
 
plot<-ggarrange(A, D, B,  E, C, f, ncol=2, nrow=3, common.legend = T, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q6_LogReg_Stool_Sputum.png", plot = plot, width = 10, height = 13)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q6_LogReg_Stool_Sputum.pdf", plot = plot, width = 10, height = 13)

rm(A, B, C, D, E, f)

###Pseudomonas and lung function 
sdt.sputum%>%
  mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  dplyr::group_by(Visit)%>%
  dplyr::mutate(Pseudomonas_RA= (`Proteobacteria-Pseudomonas`/1E6)*100)%>%
  ggplot(aes(x=Pseudomonas_RA, y= ppFEV1, shape= Visit))+
  geom_point(size=3, aes(fill= Patient_number, shape= Visit), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  geom_smooth(method=lm, se = T,aes(color= Visit))+
  theme_bw()+
  scale_fill_manual(values = pal.CF)+
  labs(tag= "A)")+
  xlab("Pseudomonas relative abundance")+
  ylab("Lung function (ppFEV1)")+
  geom_hline(yintercept=70, linetype="dashed", color = "red")+
  stat_cor(method = "spearman", label.x = 60, label.y = 80)+
  theme(text = element_text(size=16), legend.position = "none")+
  facet_grid(rows = vars(Visit))-> A

sdt.sputum%>%
  mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  dplyr::group_by(Visit)%>%
  dplyr::mutate(Pseudomonas_RA= (`Proteobacteria-Pseudomonas`/1E6)*100)%>%
  ggplot(aes(x= Pseudomonas_RA, y= ppFEV1))+
  geom_point(size=3, aes(fill= Patient_number, shape= Visit), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  scale_fill_manual(values = pal.CF)+
  geom_smooth(method=lm, se = T, color= "black")+
  theme_bw()+
  labs(tag= "B)")+
  labs(fill = "Patient")+
  labs(shape = "Visit")+
  geom_hline(yintercept=70, linetype="dashed", color = "red")+
  guides(fill = guide_legend(override.aes=list(shape=c(21)), ncol = 6), shape= guide_legend(nrow = 3))+
  xlab("Pseudomonas relative abundance")+
  ylab("Lung function (ppFEV1)")+
  stat_cor(method = "spearman", label.x = 60, label.y = 80)+ # Add sperman`s correlation coefficient
  theme(text = element_text(size=16), legend.position="bottom", legend.box = "horizontal")-> B

C<- grid.arrange(A,B, heights = c(3, 2))

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Lung_Pseudomonas_Sputum.pdf", plot = C, width = 10, height = 10)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Lung_Pseudomonas_Sputum.png", plot = C, width = 10, height = 10)

rm(A,B,C)

sdt.sputum%>%
  mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  dplyr::group_by(Visit)%>%
  dplyr::mutate(Pseudomonas_RA= (`Proteobacteria-Pseudomonas`/1E6)*100) -> tmp

model.Pseu.sputum1 <- lm(ppFEV1 ~ Pseudomonas_RA, data= tmp)
model.Pseu.sputum2 <- lm(ppFEV1 ~ Pseudomonas_RA * Visit, data= tmp)

car::Anova(model.Pseu.sputum1, type=3) 
car::Anova(model.Pseu.sputum2, type=3) 

model.Pseu.sputum.lsm <-
  lsmeans::lsmeans(model.Pseu.sputum2,
                   pairwise~Pseudomonas_RA:Visit,
                   adjust="fdr")

lf.model.sputum.lsm$contrasts


##Create biom format object for PICRUSt2
require("biomformat")
asvmat.rare<- as.matrix(PS3.sputum@otu_table)
biom.tmp<- make_biom(asvmat.rare, matrix_element_type = "int")
write_biom(biom.tmp,"CF_project/exercise-cf-intervention/data/biom_sputum.biom") ##Good biom for test

##Select sequences from the ASV in PS3.stool
library(Biostrings)
dna<- readDNAStringSet( "~/CF_project/output/ASV.fasta", format = "fasta")
keep <- data.frame(name = rownames(asvmat.rare))
names(dna)
dna.sputum<- dna[keep$name]
writeXStringSet(dna.sputum, "CF_project/exercise-cf-intervention/data/Sputum_ASV.fasta") #-> For Picrust2
