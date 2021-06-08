##Cystic-Fibrosis Microbiome Project (Mainz)
##Data analysis: Picrust2 results Stool and Sputum samples
##Víctor Hugo Jarquín-Díaz 28.04.2021

###Analysis of predicted metagenomes generated in PICRUSt2
library(lme4)
library(lmtest)
library(ggplot2)
library(reshape2)
library(vegan)
library(gtools)
library(ggpubr)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(ALDEx2)
library(tidyverse)
library(viridis)
library("car")
library("merTools")
library(sjPlot)
library(sjlabelled)
library(sjmisc)


##Scaling abundances function
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

##1) Sample data
##For this analysis we need
#training
#sdt.sputum
#sdt.stool
#top.sputum
#top.stool

##2)Predicted metagenomes

##Enzyme classification (EC) abundances per sample (Count table)
PredMet.Sput<- read.table("CF_project/picrust2_sputum_out/EC_metagenome_out/pred_metagenome_unstrat.tsv", header = T, sep = "\t")
PredMet.Stool<- read.table("CF_project/picrust2_stool_out/EC_metagenome_out/pred_metagenome_unstrat.tsv", header = T, sep = "\t")

##3) Predicted KEGG onthology (KO)
##How the ASVs contribute to KEGG onthology (KO) abundances in each sample
PredKO.Sput<- read.table("CF_project/picrust2_sputum_out/KO_metagenome_out/pred_metagenome_unstrat.tsv", header = T, sep = "\t")
PredKO.Stool<- read.table("CF_project/picrust2_stool_out/KO_metagenome_out/pred_metagenome_unstrat.tsv", header = T, sep = "\t")

##4)Predicted pathway functions
##Metabolic pathway abundances per sample
PredPath.Sput<- read.table("CF_project/picrust2_sputum_out/pathways_out/path_abun_unstrat.tsv", header = T, sep = "\t")
PredPath.Stool<- read.table("CF_project/picrust2_stool_out/pathways_out/path_abun_unstrat.tsv", header = T, sep = "\t")

##Adjust tables
#Make functions/pathways the row names 
##Sputum
###Pathways
rownames(PredPath.Sput)<-PredPath.Sput$pathway
PredPath.Sput$pathway<- NULL
PredPath.Sput<- as.matrix(PredPath.Sput)
PredPath.Sput<- otu_table(PredPath.Sput, taxa_are_rows = T)
sample_names(PredPath.Sput) <- gsub("X", "\\1", sample_names(PredPath.Sput)) ##Adjust col names to match SampleID

###KO functions
rownames(PredKO.Sput)<-PredKO.Sput$function.
PredKO.Sput$function.<- NULL
PredKO.Sput<- as.matrix(PredKO.Sput)
PredKO.Sput<- otu_table(PredKO.Sput, taxa_are_rows = T)
sample_names(PredKO.Sput) <- gsub("X", "\\1", sample_names(PredKO.Sput))

##Stool
###Pathway
rownames(PredPath.Stool)<-PredPath.Stool$pathway
PredPath.Stool$pathway<- NULL
PredPath.Stool<- as.matrix(PredPath.Stool)
PredPath.Stool<- otu_table(PredPath.Stool, taxa_are_rows = T)
sample_names(PredPath.Stool) <- gsub("X", "\\1", sample_names(PredPath.Stool)) ##Adjust col names to match SampleID

###KO functions
rownames(PredKO.Stool)<-PredKO.Stool$function.
PredKO.Stool$function.<- NULL
PredKO.Stool<- as.matrix(PredKO.Stool)
PredKO.Stool<- otu_table(PredKO.Stool, taxa_are_rows = T)
sample_names(PredKO.Stool) <- gsub("X", "\\1", sample_names(PredKO.Stool))

#Prepare sample data 
##Sputum
###Training information
training%>%
  dplyr::filter(Benzoase==1)%>%
  dplyr::select(c(SampleID, Mean_MET_V1V2:Percentage_Trainingsweeks_n52))%>%
  column_to_rownames(var = "SampleID")-> x

###Sample information
sdt.sputum%>%
  dplyr::mutate(Phenotype_severity = case_when(Phenotype_severity == 2  ~ 1,
                                               Phenotype_severity == 1 ~ 0))%>%
  dplyr::mutate(Mutation_severity = case_when(Mutation_severity == 2  ~ 1,
                                              Mutation_severity == 1 ~ 0))%>%
  dplyr::mutate(sex = case_when(sex == 1  ~ 0,
                                sex == 2 ~ 1))%>%
  cbind(x)-> y

###Dominant taxa 
top.sputum%>%
  dplyr::select(c(SampleID, ASV:n))%>%
  column_to_rownames(var = "SampleID")%>%
  cbind(y)-> tmp.sputum

y<- sample_data(tmp.sputum)
##Stool
###Training information
training%>%
  dplyr::filter(material=="Stool")%>%
  dplyr::select(c(SampleID, Mean_MET_V1V2:Percentage_Trainingsweeks_n52))%>%
  column_to_rownames(var = "SampleID")-> x

###Sample information
sdt.stool%>%
  dplyr::mutate(Phenotype_severity = case_when(Phenotype_severity == 2  ~ 1,
                                               Phenotype_severity == 1 ~ 0))%>%
  dplyr::mutate(Mutation_severity = case_when(Mutation_severity == 2  ~ 1,
                                              Mutation_severity == 1 ~ 0))%>%
  dplyr::mutate(sex = case_when(sex == 1  ~ 0,
                                sex == 2 ~ 1))%>%
  cbind(x)-> y

###Dominant taxa 
top.stool%>%
  dplyr::select(c(SampleID, ASV:n))%>%
  column_to_rownames(var = "SampleID")%>%
  cbind(y)-> tmp.stool

z<- sample_data(tmp.stool)

#Make Phyloseq objects for pathways and KO
##Sputum
###Pathways
PS.Path.sputum <- merge_phyloseq(PredPath.Sput, y)
###KO
PS.KO.sputum <- merge_phyloseq(PredKO.Sput, y)

##Stool
###Pathways
PS.Path.stool <- merge_phyloseq(PredPath.Stool, z)
###KO
PS.KO.stool <- merge_phyloseq(PredKO.Stool, z)

rm(PredKO.Sput, PredKO.Stool, PredPath.Sput, PredPath.Stool, x, y, z)

##################Generate PCAs#####################
#Sputum
##KO
##Transform dataset 
PS.KO.sputum.clr <- microbiome::transform(PS.KO.sputum, "clr") #Centered log ratio transformation
Ord.KO.sputum.clr <- phyloseq::ordinate(PS.KO.sputum.clr, "RDA") #principal components analysis
#Examine eigenvalues and % prop. variance explained
head(Ord.KO.sputum.clr$CA$eig)
sapply(Ord.KO.sputum.clr$CA$eig[1:6], function(x) x / sum(Ord.KO.sputum.clr$CA$eig))

##KOs contributing into PC1 and PC2
ind.coord <- data.frame(Ord.KO.sputum.clr$CA$v)
sdev_ind <- apply(ind.coord, 1, sd)
ind_cont_PCA1 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC1^2 /sdev_ind))))
ind_cont_PCA1_top <- ind_cont_PCA1 %>% 
  rownames_to_column("otu") %>% 
  filter(PCA >= 0.0015)%>% 
  column_to_rownames("otu")
ind_cont_PCA2 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC2^2 /sdev_ind))))
ind_cont_PCA2_top <- ind_cont_PCA2 %>% 
  rownames_to_column("otu") %>% 
  filter(PCA >= 0.0015)%>% 
  column_to_rownames("otu")

ind_cont_PCA_top.sputum <- rbind(ind_cont_PCA1_top, ind_cont_PCA2_top)
sum(ind_cont_PCA1_top$PCA) / sum(ind_cont_PCA1$PCA)
nrow(ind_cont_PCA1_top)
##303 KO contribute for the 30% of the variation in PC1
sum(ind_cont_PCA2_top$PCA) / sum(ind_cont_PCA2$PCA)
nrow(ind_cont_PCA2_top)
##122 KO contribute for the 18.9% of the variation in PC2

ind_cont_PCA_top.sputum%>%
  rownames_to_column(var = "KO")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 25)-> ind_cont_PCA_top.sputum

###Plot PCA
plot_ordination(PS.KO.sputum.clr, ordination = Ord.KO.sputum.clr)+ 
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(fill= Genus, shape= Visit), color= "black")+
  scale_shape_manual(values = c(21, 24, 22))+
  scale_fill_manual(values = tax.palette)+
  labs(title = "PCA (centered-log ratio KO prediction)",tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(fill = "Dominant taxa")+
  labs(shape = "Visit")+
  xlab(paste0("PC 1 [", round(Ord.KO.sputum.clr$CA$eig[1] / sum(Ord.KO.sputum.clr$CA$eig)*100, digits = 2), "%]"))+
  ylab(paste0("PC 2 [", round(Ord.KO.sputum.clr$CA$eig[2] / sum(Ord.KO.sputum.clr$CA$eig)*100, digits = 2), "%]"))-> A

## Bray Curtis
BC_dist<- phyloseq::distance(PS.KO.sputum, method="bray", weighted=F)
Ord.KO.sputum <- ordinate(PS.KO.sputum, method="PCoA", distance="bray") 

plot_ordination(PS.KO.sputum, ordination = Ord.KO.sputum)+ 
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(fill= Genus, shape= Visit), color= "black")+
  scale_shape_manual(values = c(21, 24, 22))+
  scale_fill_manual(values = tax.palette)+
  labs(title = "Bray-Curtis dissimilariy (KO prediction)",tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(fill = "Dominant taxa")+
  labs(shape = "Visit")+
  xlab(paste0("PCo 1 [", round(Ord.KO.sputum$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(Ord.KO.sputum$values[2,2]*100, digits = 2), "%]")) -> B

BC.KO.sputum<- vegan::adonis(BC_dist~ Phenotype_severity + Mutation_severity + Genus + sex + age +  Visit + BMI,
                              permutations = 999, data = tmp.sputum, na.action = F, strata = tmp.sputum$Patient_number)

kable(BC.KO.sputum$aov.tab)#-> Report 

C<- ggarrange(A, B, ncol=1, nrow=2, common.legend = TRUE, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q7_KO_pred_ord_sputum.pdf", plot = C, width = 8, height = 10)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q7_KO_pred_ord_sputum.png", plot = C, width = 8, height = 10)

rm(A,B,C)

##Pathways
PS.Path.sputum.clr <- microbiome::transform(PS.Path.sputum, "clr")  
Ord.Path.sputum.clr <- phyloseq::ordinate(PS.Path.sputum.clr, "RDA")
#Examine eigenvalues and % prop. variance explained
head(Ord.Path.sputum.clr$CA$eig)
sapply(Ord.Path.sputum.clr$CA$eig[1:6], function(x) x / sum(Ord.Path.sputum.clr$CA$eig))

ind.coord.pw <- data.frame(Ord.Path.sputum.clr$CA$v)
sdev_ind.pw <- apply(ind.coord.pw, 1, sd)
ind_cont_PCA1.pw <- data.frame(PCA = (100*(1 / nrow(ind.coord.pw)*(ind.coord.pw$PC1^2 /sdev_ind.pw))))
ind_cont_PCA1_top.pw <- ind_cont_PCA1.pw %>% 
  rownames_to_column("otu") %>% 
  filter(PCA >= 0.008)%>% 
  column_to_rownames("otu")
ind_cont_PCA2.pw <- data.frame(PCA = (100*(1 / nrow(ind.coord.pw)*(ind.coord.pw$PC2^2 /sdev_ind.pw))))
ind_cont_PCA2_top.pw <- ind_cont_PCA2.pw %>% 
  rownames_to_column("otu") %>% 
  filter(PCA >= 0.02)%>% 
  column_to_rownames("otu")
ind_cont_PCA_top.pw.sputum <- rbind(ind_cont_PCA1_top.pw, ind_cont_PCA2_top.pw)

ind_cont_PCA_top.pw.sputum%>%
  rownames_to_column(var = "Pathway")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 26)-> ind_cont_PCA_top.pw.sputum

##Plot
plot_ordination(PS.Path.sputum.clr, ordination = Ord.Path.sputum.clr)+
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(fill= Genus, shape= Visit), color= "black")+
  scale_shape_manual(values = c(21, 24, 22))+
  scale_fill_manual(values = tax.palette)+
  labs(title = "PCA (centered-log ratio pathway prediction)",tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(fill = "Dominant taxa")+
  labs(shape = "Visit")+
  xlab(paste0("PC 1 [", round(Ord.Path.sputum.clr$CA$eig[1] / sum(Ord.Path.sputum.clr$CA$eig)*100, digits = 2), "%]"))+
  ylab(paste0("PC 2 [", round(Ord.Path.sputum.clr$CA$eig[2] / sum(Ord.Path.sputum.clr$CA$eig)*100, digits = 2), "%]"))-> A

## Bray Curtis
BC_dist<- phyloseq::distance(PS.Path.sputum, method="bray", weighted=F)
Ord.path.sputum <- ordinate(PS.Path.sputum, method="PCoA", distance="bray") 

plot_ordination(PS.Path.sputum, ordination = Ord.path.sputum)+ 
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(fill= Genus, shape= Visit), color= "black")+
  scale_shape_manual(values = c(21, 24, 22))+
  scale_fill_manual(values = tax.palette)+
  labs(title = "Bray-Curtis dissimilariy (Pathway prediction)",tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(fill = "Dominant taxa")+
  labs(shape = "Visit")+
  xlab(paste0("PCo 1 [", round(Ord.path.sputum$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(Ord.path.sputum$values[2,2]*100, digits = 2), "%]")) -> B

BC.path.sputum<- vegan::adonis(BC_dist~ Phenotype_severity + Mutation_severity + Genus + sex + age +  Visit + BMI,
                             permutations = 999, data = tmp.sputum, na.action = F, strata = tmp.sputum$Patient_number)

kable(BC.path.sputum$aov.tab)#-> Report 

C<- ggarrange(A, B, ncol=1, nrow=2, common.legend = TRUE, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q7_path_pred_ord_sputum.pdf", plot = C, width = 8, height = 10)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q7_path_pred_ord_sputum.png", plot = C, width = 8, height = 10)

####Heat map pathways#####
##Matrix Pathways sample
PathM<- PS.Path.sputum.clr@otu_table
top<- ind_cont_PCA_top.pw.sputum$Pathway

PathM<- PathM[rownames(PathM) %in% top]

P.clust <- hclust(dist(t(PathM)), method = "complete") ##Dendogram

as.dendrogram(P.clust) %>%
  plot(horiz = T)

P.col <- cutree(tree = P.clust, k = 2)
P.col  <- data.frame(cluster = ifelse(test = P.col  == 1, yes = "cluster 1", no = "cluster 2"))
P.col$SampleID <- rownames(P.col)

P.col<- cbind(P.col, tmp.sputum)

col_groups <- P.col %>%
  mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  dplyr::select(c(SampleID, Patient_number, Visit, ppFEV1, Genus)) ##Here It is possible to add the other characteristics

row.names(col_groups)<- col_groups$SampleID

col_groups$SampleID<- NULL

colour_groups <- list(Patient_number= pal.CF%in%tmp.sputum$Patient_number, Genus= c("Streptococcus"= "#925E9FFF","Staphylococcus"= "#008B45FF", 
                                                       "Stenotrophomonas"= "#B09C85FF", "Pseudomonas" = "#ED0000FF"))

sputum.heatmap <- pheatmap(PathM, cluster_rows = F, cluster_cols = T,
                          color = colorRampPalette(c("white","#832424FF"))(100), #"#3A3A98FF",
                          border_color = NA,
                          annotation_col = col_groups, 
                          annotation_colors = colour_groups,
                          show_rownames = T,
                          show_colnames = F)

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q7_path_pred_heatmap_sputum.pdf", plot = sputum.heatmap, width = 8, height = 10)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q7_path_pred_heatmap_sputum.png", plot = sputum.heatmap, width = 8, height = 10)
