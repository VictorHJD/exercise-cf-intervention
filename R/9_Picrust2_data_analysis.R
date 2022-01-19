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
##Descriptions
PredKO.Sput.des<- read.table("CF_project/picrust2_sputum_out/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv", header = T, sep = "\t")
PredKO.Stool.des<- read.table("CF_project/picrust2_stool_out/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv", header = T, sep = "\t")

##4)Predicted pathway functions
##Metabolic pathway abundances per sample
PredPath.Sput<- read.table("CF_project/picrust2_sputum_out/pathways_out/path_abun_unstrat.tsv", header = T, sep = "\t")
PredPath.Stool<- read.table("CF_project/picrust2_stool_out/pathways_out/path_abun_unstrat.tsv", header = T, sep = "\t")
##Descriptions
PredPath.Sput.des<- read.table("CF_project/picrust2_sputum_out/pathways_out/path_abun_unstrat_descrip.tsv", header = T, sep = "\t")
PredPath.Stool.des<- read.table("CF_project/picrust2_stool_out/pathways_out/path_abun_unstrat_descrip.tsv", header = T, sep = "\t")


##Adjust tables
#Make functions/pathways the row names 
##Sputum
###Pathways
rownames(PredPath.Sput)<-PredPath.Sput$pathway
PredPath.Sput$pathway<- NULL
PredPath.Sput<- as.matrix(PredPath.Sput)
PredPath.Sput<- otu_table(PredPath.Sput, taxa_are_rows = T)
sample_names(PredPath.Sput) <- gsub("X", "\\1", sample_names(PredPath.Sput)) ##Adjust col names to match SampleID
###Descriptions
PredPath.Sput.des%>% 
  dplyr::select(c(pathway,description))%>%
  dplyr::rename(Pathway= pathway)->PredPath.Sput.des

###KO functions
rownames(PredKO.Sput)<-PredKO.Sput$function.
PredKO.Sput$function.<- NULL
PredKO.Sput<- as.matrix(PredKO.Sput)
PredKO.Sput<- otu_table(PredKO.Sput, taxa_are_rows = T)
sample_names(PredKO.Sput) <- gsub("X", "\\1", sample_names(PredKO.Sput))
###Descriptions
PredKO.Sput.des%>% 
  dplyr::select(c(function.,description))%>%
  dplyr::rename(KO= function.)->PredKO.Sput.des

##Stool
###Pathway
rownames(PredPath.Stool)<-PredPath.Stool$pathway
PredPath.Stool$pathway<- NULL
PredPath.Stool<- as.matrix(PredPath.Stool)
PredPath.Stool<- otu_table(PredPath.Stool, taxa_are_rows = T)
sample_names(PredPath.Stool) <- gsub("X", "\\1", sample_names(PredPath.Stool)) ##Adjust col names to match SampleID
###Descriptions
PredPath.Stool.des%>% 
  dplyr::select(c(pathway,description))%>%
  dplyr::rename(Pathway= pathway)->PredPath.Stool.des

###KO functions
rownames(PredKO.Stool)<-PredKO.Stool$function.
PredKO.Stool$function.<- NULL
PredKO.Stool<- as.matrix(PredKO.Stool)
PredKO.Stool<- otu_table(PredKO.Stool, taxa_are_rows = T)
sample_names(PredKO.Stool) <- gsub("X", "\\1", sample_names(PredKO.Stool))
###Descriptions
PredKO.Stool.des%>% 
  dplyr::select(c(function.,description))%>%
  dplyr::rename(KO= function.)->PredKO.Stool.des

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

##################Generate Sputum PCAs#####################
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
  slice_head(n = 25)%>%
  left_join(PredKO.Sput.des, by= "KO")-> ind_cont_PCA_top.sputum

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
  slice_head(n = 26)%>%
  left_join(PredPath.Sput.des, by= "Pathway")-> ind_cont_PCA_top.pw.sputum

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

rm(A, B, C)
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

####Based on the dominance make new PCo for sputum 
##Bray-Curtis
PS4.sput.2<- PS4.sput
PS4.sput.2@sam_data<- sample_data(tmp.sputum)

BC_dist<- phyloseq::distance(PS4.sput.2,
                             method="bray", weighted=F)
ordination<- ordinate(PS4.sput.2,
                      method="PCoA", distance= BC_dist)

##Permanova
BC.test.sputum2<- vegan::adonis(BC_dist~ Phenotype_severity+ Mutation_severity + Genus + sex + age +  Visit + BMI,
                               permutations = 999, data = tmp.sputum, na.action = F, strata = tmp.sputum$Patient_number)

kable(BC.test.sputum2$aov.tab)

## Calculate multivariate dispersion (aka distance to the centroid)
mvd<- vegan::betadisper(BC_dist, tmp.sputum$Genus, type = "centroid")
mvd.perm<- vegan::permutest(mvd, permutations = 999)

plot(mvd)
anova(mvd)
boxplot(mvd)
plot(TukeyHSD(mvd))

##Extract centroids and vectors 
centroids<-data.frame(grps=rownames(mvd$centroids), data.frame(mvd$centroids))
vectors<-data.frame(group=mvd$group,data.frame(mvd$vectors))

seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("Genus","v.PCoA1","v.PCoA2","PCoA1","PCoA2")
seg.data$Genus<- NULL
##Add Line information for each sample 
tmp.sputum%>%
  dplyr::select(c(Genus, Visit, Patient_number, Phenotype_severity))%>%
  dplyr::mutate(Phenotype_severity = as.factor(Phenotype_severity))%>%
  rownames_to_column(var = "SampleID")%>%
  dplyr::group_by(Genus)%>%
  unique()%>%
  column_to_rownames(var = "SampleID")%>%
  cbind(seg.data)-> seg.data


##Plot BC
ggplot() + 
  geom_point(data=centroids[,1:4], aes(x=PCoA1,y=PCoA2, color= grps, group=grps),size=4, shape= 4) + 
  geom_line(data=centroids[,1:4], aes(x=PCoA1,y=PCoA2), linetype = 2, color=" gray")+
  geom_segment(data=centroids[,1:4], aes(x=centroids[1,2], y=centroids[1,3], xend=centroids[2,2], yend=centroids[2,3]), linetype = 2, color=" gray")+ 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2,shape=Phenotype_severity, fill= Genus),size=3) +
  #geom_point(size=3, aes(fill= Patient_number, shape= as.factor(Phenotype_severity)), color= "black")+
  scale_shape_manual(values = c(25, 24), labels = c("Low (ppFEV1 > 70% at V1)", "High (ppFEV1 < 70% at V1)"))+
  #scale_shape_manual(values = c(21, 24, 22))+
  scale_fill_manual(values = tax.palette)+
  scale_color_manual(values = tax.palette)+
  stat_ellipse(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2,color= Genus))+
  labs(title = "Bray-Curtis dissimilarity",tag= "A)")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= F)+
  theme_bw()+
  theme(text = element_text(size=16))+
  labs(fill = "Dominant taxa")+
  #labs(shape = "Visit")+
  labs(shape = "Phenotype severity")+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> A

# Horizontal marginal boxplot - to appear at the top of the chart
seg.data%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  wilcox_test(v.PCoA1 ~ Genus)%>%
  adjust_pvalue(method = "BH") %>%
  add_significance()%>%
  add_xy_position(x = "Genus")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "~/CF_project/exercise-cf-intervention/tables/Q2_PCoA1_Dom_tax_Sput.csv")

seg.data%>%
  ggplot(aes(x=Genus, y= v.PCoA1))+
  geom_boxplot(aes(fill= Genus), color= "black", alpha= 0.5)+
  #geom_jitter(position = position_jitter(width = 0.5),size=3, aes(shape=Visit, fill= Genus)) +
  #scale_shape_manual(values = c(21, 24, 22))+
  scale_fill_manual(values = tax.palette)+
  scale_color_manual(values = tax.palette)+
  theme_bw()+
  stat_pvalue_manual(stats.test, hide.ns = T,
                     tip.length = 0, label = "{p.adj.signif}", coord.flip = T)+
  coord_flip()+ 
  theme(axis.text = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none",
        plot.margin = unit(c(1, 0.2, -0.5, 0.5), "lines"))-> xplot

# Vertical marginal boxplot - to appear at the right of the chart
seg.data%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  wilcox_test(v.PCoA2 ~ Genus)%>%
  adjust_pvalue(method = "BH") %>%
  add_significance()%>%
  add_xy_position(x = "Genus")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "~/CF_project/exercise-cf-intervention/tables/Q2_PCoA2_Dom_tax_Sput.csv")

seg.data%>%
  ggplot(aes(x=Genus, y= v.PCoA2))+
  geom_boxplot(aes(fill= Genus), color= "black", alpha= 0.5)+
  #geom_jitter(position = position_jitter(width = 0.5),size=3, aes(shape=Visit, fill= Genus)) +
  #scale_shape_manual(values = c(21, 24, 22))+
  scale_fill_manual(values = tax.palette)+
  scale_color_manual(values = tax.palette)+
  theme_bw()+
  stat_pvalue_manual(stats.test, hide.ns = T,
                     tip.length = 0, label = "{p.adj.signif}")+
  theme(axis.text = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none",
        plot.margin = unit(c(0.2, 1, 0.5, -0.5), "lines"))->  yplot 

require(cowplot)
p1<- insert_xaxis_grob(A, xplot, grid::unit(.25, "null"), position = "top")
p2<- insert_yaxis_grob(p1, yplot, grid::unit(.25, "null"), position = "right")
A2<- ggdraw(p2)

#### By visit
BC_dist.sputum%>% 
  wilcox_test(BC_dist ~ Group)%>%
  adjust_pvalue(method = "BH") %>%
  add_significance()%>%
  add_xy_position(x = "Group")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "~/CF_project/exercise-cf-intervention/tables/Q3_Sample_Visit_BC_Sputum.csv")

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
  labs(tag= "B)", caption = get_pwc_label(stats.test), fill = "Patient")+
  scale_fill_manual(values = pal.CF)+
  theme_classic()+
  guides(fill = guide_legend(override.aes=list(shape=c(21)), ncol = 6, size= 10), color= "none")+
  theme(text = element_text(size=16), legend.position="bottom", legend.box = "horizontal")+
  scale_x_discrete(labels=c("V1_V2" =  "V1 to V2", 
                            "V2_V3" = "V2 to V3",
                            "V1_V3" = "V1 to V3"))+
  stat_pvalue_manual(stats.test, hide.ns = F, step.increase = 0.05,
                     tip.length = 0, label = "{p.adj} {p.adj.signif}")->B

C<- ggarrange(A, B, ncol=1, nrow=2, common.legend = F)

ggsave(file = "CF_project/exercise-cf-intervention/figures/Fig_3_Sput.pdf", plot = C, width = 10, height = 12)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Fig_3_Sput.png", plot = C, width = 10, height = 12)

C<- ggarrange(A2, B, ncol=1, nrow=2, common.legend = F)

ggsave(file = "CF_project/exercise-cf-intervention/figures/Fig_3.2_Sput.pdf", plot = C, width = 10, height = 12)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Fig_3.2_Sput.png", plot = C, width = 10, height = 12)

##Extract distance to the centroid for each samples 
centroid.dist<- data.frame(mvd$distances)

centroid.dist%>%
  dplyr::rename(cent.dist= mvd.distances)%>%
  cbind(seg.data)-> seg.data

##Plot distances comparison 
seg.data%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  wilcox_test(cent.dist ~ Genus)%>%
  adjust_pvalue(method = "BH") %>%
  add_significance()%>%
  add_xy_position(x = "Genus")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "~/CF_project/exercise-cf-intervention/tables/Q2_Dist_Centroid_Dom_tax_Sput.csv")

seg.data%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Genus, y= cent.dist))+
  geom_boxplot(color="black", alpha= 0.5)+
  geom_point(position=position_jitter(0.2), size=3, aes(fill= Patient_number, shape= Visit), color= "black")+
  scale_shape_manual(values = c(21, 24, 22))+
  xlab("Dominant genus")+
  ylab("Bray-Curtis dissimilarity \n (Distance to centroid)")+
  labs(tag= "B)", caption = get_pwc_label(stats.test), fill = "Patient")+
  guides(fill = guide_legend(override.aes=list(shape=c(21)), size= 10), color= "none")+
  scale_fill_manual(values = pal.CF)+
  theme_classic()+
  theme(text = element_text(size=16), legend.position="bottom", legend.box = "horizontal")+
  stat_pvalue_manual(stats.test, hide.ns = T, step.increase = 0.05,
                     tip.length = 0, label = "{p.adj} {p.adj.signif}")-> B

#C<- ggarrange(A, B, ncol=1, nrow=2, common.legend = F)

#ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Dominant_Sput.pdf", plot = C, width = 10, height = 11)
#ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Dominant_Sput.png", plot = C, width = 10, height = 11)

#C<- ggarrange(A2, B, ncol=1, nrow=2, common.legend = F)

#ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Dominant_Sput_2.pdf", plot = C, width = 10, height = 12)
#ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Dominant_Sput_2.png", plot = C, width = 10, height = 12)

rm(A, A2, B, C)

##Dominant taxa for sputum estimate change in relative abundance 
##Staphylococcus
tmp1<- as.data.frame(sputum.microbiome[,"Firmicutes-Staphylococcus"])
colnames(tmp1)<- c("Firmicutes-Staphylococcus")
tmp1%>%
  dplyr::mutate(`Firmicutes-Staphylococcus`= round((`Firmicutes-Staphylococcus`/1e6)*100, digits = 2))-> tmp1

##Pseudomonas
tmp2<- as.data.frame(sputum.microbiome[,"Proteobacteria-Pseudomonas"])
colnames(tmp2)<- c("Proteobacteria-Pseudomonas")
tmp2%>%
  dplyr::mutate(`Proteobacteria-Pseudomonas`= round((`Proteobacteria-Pseudomonas`/1e6)*100, digits = 2))-> tmp2

##Streptococcus
tmp3<- as.data.frame(sputum.microbiome[,"Firmicutes-Streptococcus"])
colnames(tmp3)<- c("Firmicutes-Streptococcus")
tmp3%>%
  dplyr::mutate(`Firmicutes-Streptococcus`= round((`Firmicutes-Streptococcus`/1e6)*100, digits = 2))-> tmp3

##Stenotrophomonas
tmp4<- as.data.frame(sputum.microbiome[,"Proteobacteria-Stenotrophomonas"])
colnames(tmp4)<- c("Proteobacteria-Stenotrophomonas")
tmp4%>%
  dplyr::mutate(`Proteobacteria-Stenotrophomonas`= round((`Proteobacteria-Stenotrophomonas`/1e6)*100, digits = 2))-> tmp4

tmp.sputum%>%
  cbind(tmp1)%>%
  cbind(tmp2)%>%
  cbind(tmp3)%>%
  cbind(tmp4)-> tmp.sputum

##Estimate deltas of abundance per Visit combination 
tmp.sputum%>%
  dplyr::select(c(Patient_number, Visit, `Firmicutes-Staphylococcus`))%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  dplyr::mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  dplyr::distinct()%>%
  spread(Visit, `Firmicutes-Staphylococcus`)%>%
  dplyr::mutate(V1_V2= V1 - V2)%>%
  dplyr::mutate(V1_V3= V1 - V3)%>%
  dplyr::mutate(V2_V3= V2 - V3)%>%
  dplyr::select(c(Patient_number, V1_V2, V2_V3, V1_V3))%>%
  pivot_longer(!Patient_number, names_to = "Group", values_to = "Firmicutes-Staphylococcus")%>%
  dplyr::mutate(ID= paste0(Patient_number, Group))%>%
  dplyr::ungroup()%>%
  dplyr::select(c(ID, `Firmicutes-Staphylococcus`))-> tmp1

BC_dist.sputum%>%
  left_join(tmp1, by="ID")-> BC_dist.sputum

tmp.sputum%>%
  dplyr::select(c(Patient_number, Visit, `Proteobacteria-Pseudomonas`))%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  dplyr::mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  dplyr::distinct()%>%
  spread(Visit, `Proteobacteria-Pseudomonas`)%>%
  dplyr::mutate(V1_V2= V1 - V2)%>%
  dplyr::mutate(V1_V3= V1 - V3)%>%
  dplyr::mutate(V2_V3= V2 - V3)%>%
  dplyr::select(c(Patient_number, V1_V2, V2_V3, V1_V3))%>%
  pivot_longer(!Patient_number, names_to = "Group", values_to = "Proteobacteria-Pseudomonas")%>%
  dplyr::mutate(ID= paste0(Patient_number, Group))%>%
  dplyr::ungroup()%>%
  dplyr::select(c(ID, `Proteobacteria-Pseudomonas`))-> tmp2

BC_dist.sputum%>%
  left_join(tmp2, by="ID")-> BC_dist.sputum

tmp.sputum%>%
  dplyr::select(c(Patient_number, Visit, `Firmicutes-Streptococcus`))%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  dplyr::mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  dplyr::distinct()%>%
  spread(Visit, `Firmicutes-Streptococcus`)%>%
  dplyr::mutate(V1_V2= V1 - V2)%>%
  dplyr::mutate(V1_V3= V1 - V3)%>%
  dplyr::mutate(V2_V3= V2 - V3)%>%
  dplyr::select(c(Patient_number, V1_V2, V2_V3, V1_V3))%>%
  pivot_longer(!Patient_number, names_to = "Group", values_to = "Firmicutes-Streptococcus")%>%
  dplyr::mutate(ID= paste0(Patient_number, Group))%>%
  dplyr::ungroup()%>%
  dplyr::select(c(ID, `Firmicutes-Streptococcus`))-> tmp3

BC_dist.sputum%>%
  left_join(tmp3, by="ID")-> BC_dist.sputum

tmp.sputum%>%
  dplyr::select(c(Patient_number, Visit, `Proteobacteria-Stenotrophomonas`))%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  dplyr::mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  dplyr::distinct()%>%
  spread(Visit, `Proteobacteria-Stenotrophomonas`)%>%
  dplyr::mutate(V1_V2= V1 - V2)%>%
  dplyr::mutate(V1_V3= V1 - V3)%>%
  dplyr::mutate(V2_V3= V2 - V3)%>%
  dplyr::select(c(Patient_number, V1_V2, V2_V3, V1_V3))%>%
  pivot_longer(!Patient_number, names_to = "Group", values_to = "Proteobacteria-Stenotrophomonas")%>%
  dplyr::mutate(ID= paste0(Patient_number, Group))%>%
  dplyr::ungroup()%>%
  dplyr::select(c(ID, `Proteobacteria-Stenotrophomonas`))-> tmp4

BC_dist.sputum%>%
  left_join(tmp4, by="ID")-> BC_dist.sputum

##Model time 
BC_dist.sputum%>%
  dplyr::select(ppFVC, Trainingfrequency, Trainingtime, ppFEV1, BC_dist, Patient_number, Group, 
                `Firmicutes-Staphylococcus`,`Proteobacteria-Pseudomonas`, 
                `Firmicutes-Streptococcus`, `Proteobacteria-Stenotrophomonas`)%>%
  dplyr::filter(complete.cases(.))-> tmp

BC_dist.sputum%>%
  dplyr::filter(`Firmicutes-Staphylococcus`!= 0)%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= (`Firmicutes-Staphylococcus`)*-1, y= (ppFEV1)*-1))+
  geom_point(size=3, aes(fill= Patient_number, shape= Group), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  geom_smooth(method=lm, se = F, color= "black")+
  scale_fill_manual(values = pal.CF)+
  #geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=70),alpha=0.01,fill="grey")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  geom_vline(xintercept=0, linetype="dashed", color = "red")+
  theme_bw()+
  labs(tag= "B)")+
  labs(fill = "Patient")+
  labs(shape = "Visit")+
  guides(fill = guide_legend(override.aes=list(shape=c(21)), ncol = 6), shape= guide_legend(nrow = 3))+
  xlab("Change in Staphylococcus relative abundance")+
  ylab("Change in lung function (Δ ppFEV1)")+
  theme(text = element_text(size=16), legend.position="bottom", legend.box = "horizontal")
  
##Model selection do it with glm
tmp%>%
  dplyr::select(Trainingfrequency, Trainingtime, ppFEV1, 
                `Firmicutes-Staphylococcus`,`Proteobacteria-Pseudomonas`, 
                `Firmicutes-Streptococcus`, `Proteobacteria-Stenotrophomonas`)-> tmp.ppFEV1

full.model<- glm(ppFEV1 ~ ., data = tmp.ppFEV1) ##Full model
# Stepwise regression model
step.model <- MASS::stepAIC(full.model, direction = "both", 
                            trace = FALSE)
summary(step.model)
step.model.df<- as.data.frame(coef(summary(step.model)))

step.model.df%>%
  mutate(p.adj = p.adjust(`Pr(>|t|)`, method='BH')) %>%
  add_significance()%>%
  rownames_to_column()-> step.model.df

write.csv(step.model.df, "~/CF_project/exercise-cf-intervention/tables/Q2_Model_ppFEV1_Dominant_taxa_Sputum.csv", row.names = F)

##Mixed effect models
##with patient and intervisit group as individual random effects

tr0<- lmer(ppFEV1 ~ Trainingfrequency + (1 | Patient_number), data = tmp) ##Null model
tr1<- lmer(ppFEV1 ~ Trainingfrequency + BC_dist + (1 | Patient_number), data = tmp) ##Null model
tr2<-lmer(ppFEV1 ~ Trainingfrequency + `Firmicutes-Staphylococcus` + (1 | Patient_number), data = tmp)

lrtest(tr0, tr1)
lrtest(tr1, tr2)
lrtest(tr0, tr2)

##Model selection do it with glm
tmp%>%
  dplyr::select(Trainingfrequency, Trainingtime, ppFVC, 
                `Firmicutes-Staphylococcus`,`Proteobacteria-Pseudomonas`, 
                `Firmicutes-Streptococcus`, `Proteobacteria-Stenotrophomonas`)-> tmp.ppFVC

full.model<- glm(ppFVC ~ ., data = tmp.ppFVC) ##Full model
# Stepwise regression model
step.model <- MASS::stepAIC(full.model, direction = "both", 
                            trace = FALSE)
summary(step.model)
step.model.df<- as.data.frame(coef(summary(step.model)))

step.model.df%>%
  mutate(p.adj = p.adjust(`Pr(>|t|)`, method='BH')) %>%
  add_significance()%>%
  rownames_to_column()-> step.model.df

write.csv(step.model.df, "~/CF_project/exercise-cf-intervention/tables/Q2_Model_ppFVC_Dominant_taxa_Sputum.csv", row.names = F)

##Mixed effect models
##with patient and intervisit group as individual random effects

tr0<- lmer(ppFVC ~ Trainingfrequency + (1 | Patient_number), data = tmp) ##Null model
tr1<- lmer(ppFVC ~ Trainingfrequency + BC_dist + (1 | Patient_number), data = tmp) ##Null model
tr2<-lmer(ppFVC ~ Trainingfrequency + `Firmicutes-Staphylococcus` + (1 | Patient_number), data = tmp)

lrtest(tr0, tr2)
lrtest(tr1, tr2)

##################Generate Stool PCAs#####################
##KO
##Transform dataset 
PS.KO.stool.clr <- microbiome::transform(PS.KO.stool, "clr") #Centered log ratio transformation
Ord.KO.stool.clr <- phyloseq::ordinate(PS.KO.stool.clr, "RDA") #principal components analysis
#Examine eigenvalues and % prop. variance explained
head(Ord.KO.stool.clr$CA$eig)
sapply(Ord.KO.stool.clr$CA$eig[1:6], function(x) x / sum(Ord.KO.stool.clr$CA$eig))

##KOs contributing into PC1 and PC2
ind.coord <- data.frame(Ord.KO.stool.clr$CA$v)
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

ind_cont_PCA_top.stool <- rbind(ind_cont_PCA1_top, ind_cont_PCA2_top)
sum(ind_cont_PCA1_top$PCA) / sum(ind_cont_PCA1$PCA)
nrow(ind_cont_PCA1_top)
##488 KO contribute for the 47.3% of the variation in PC1
sum(ind_cont_PCA2_top$PCA) / sum(ind_cont_PCA2$PCA)
nrow(ind_cont_PCA2_top)
##90 KO contribute for the 20.4% of the variation in PC2

ind_cont_PCA_top.stool%>%
  rownames_to_column(var = "KO")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 25)%>%
  left_join(PredKO.Stool.des, by= "KO")-> ind_cont_PCA_top.stool

###Plot PCA
plot_ordination(PS.KO.stool.clr, ordination = Ord.KO.stool.clr)+ 
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
  xlab(paste0("PC 1 [", round(Ord.KO.stool.clr$CA$eig[1] / sum(Ord.KO.stool.clr$CA$eig)*100, digits = 2), "%]"))+
  ylab(paste0("PC 2 [", round(Ord.KO.stool.clr$CA$eig[2] / sum(Ord.KO.stool.clr$CA$eig)*100, digits = 2), "%]"))-> A

## Bray Curtis
BC_dist<- phyloseq::distance(PS.KO.stool, method="bray", weighted=F)
Ord.KO.stool <- ordinate(PS.KO.stool, method="PCoA", distance="bray") 

plot_ordination(PS.KO.stool, ordination = Ord.KO.stool)+ 
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
  xlab(paste0("PCo 1 [", round(Ord.KO.stool$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(Ord.KO.stool$values[2,2]*100, digits = 2), "%]")) -> B

BC.KO.stool<- vegan::adonis(BC_dist~ Phenotype_severity + Mutation_severity + Genus + sex + age +  Visit + BMI,
                             permutations = 999, data = tmp.stool, na.action = F, strata = tmp.stool$Patient_number)

kable(BC.KO.stool$aov.tab)#-> Report 

C<- ggarrange(A, B, ncol=1, nrow=2, common.legend = TRUE, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q7_KO_pred_ord_stool.pdf", plot = C, width = 8, height = 10)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q7_KO_pred_ord_stool.png", plot = C, width = 8, height = 10)

rm(A,B,C)

##Pathways
PS.Path.stool.clr <- microbiome::transform(PS.Path.stool, "clr")  
Ord.Path.stool.clr <- phyloseq::ordinate(PS.Path.stool.clr, "RDA")
#Examine eigenvalues and % prop. variance explained
head(Ord.Path.stool.clr$CA$eig)
sapply(Ord.Path.stool.clr$CA$eig[1:6], function(x) x / sum(Ord.Path.stool.clr$CA$eig))

ind.coord.pw <- data.frame(Ord.Path.stool.clr$CA$v)
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
ind_cont_PCA_top.pw.stool <- rbind(ind_cont_PCA1_top.pw, ind_cont_PCA2_top.pw)

ind_cont_PCA_top.pw.stool%>%
  rownames_to_column(var = "Pathway")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 26)%>%
  left_join(PredPath.Sput.des, by= "Pathway")-> ind_cont_PCA_top.pw.stool

##Plot
plot_ordination(PS.Path.stool.clr, ordination = Ord.Path.stool.clr)+
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
  xlab(paste0("PC 1 [", round(Ord.Path.stool.clr$CA$eig[1] / sum(Ord.Path.stool.clr$CA$eig)*100, digits = 2), "%]"))+
  ylab(paste0("PC 2 [", round(Ord.Path.stool.clr$CA$eig[2] / sum(Ord.Path.stool.clr$CA$eig)*100, digits = 2), "%]"))-> A

## Bray Curtis
BC_dist<- phyloseq::distance(PS.Path.stool, method="bray", weighted=F)
Ord.path.stool <- ordinate(PS.Path.stool, method="PCoA", distance="bray") 

plot_ordination(PS.Path.stool, ordination = Ord.path.stool)+ 
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
  xlab(paste0("PCo 1 [", round(Ord.path.stool$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(Ord.path.stool$values[2,2]*100, digits = 2), "%]")) -> B

BC.path.stool<- vegan::adonis(BC_dist~ Phenotype_severity + Mutation_severity + Genus + sex + age +  Visit + BMI,
                               permutations = 999, data = tmp.stool, na.action = F, strata = tmp.stool$Patient_number)

kable(BC.path.stool$aov.tab)#-> Report 

C<- ggarrange(A, B, ncol=1, nrow=2, common.legend = TRUE, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q7_path_pred_ord_stool.pdf", plot = C, width = 8, height = 10)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q7_path_pred_ord_stool.png", plot = C, width = 8, height = 10)

rm(A, B, C)
####Heat map pathways#####
##Matrix Pathways sample
PathM<- PS.Path.stool.clr@otu_table
top<- ind_cont_PCA_top.pw.stool$Pathway

PathM<- PathM[rownames(PathM) %in% top]

P.clust <- hclust(dist(t(PathM)), method = "complete") ##Dendogram

as.dendrogram(P.clust) %>%
  plot(horiz = T)

P.col <- cutree(tree = P.clust, k = 2)
P.col  <- data.frame(cluster = ifelse(test = P.col  == 1, yes = "cluster 1", no = "cluster 2"))
P.col$SampleID <- rownames(P.col)

P.col<- cbind(P.col, tmp.stool)

col_groups <- P.col %>%
  mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  dplyr::select(c(SampleID, Patient_number, Visit, ppFEV1, Genus)) ##Here It is possible to add the other characteristics

row.names(col_groups)<- col_groups$SampleID

col_groups$SampleID<- NULL

colour_groups <- list(Patient_number= pal.CF, Genus= c("Prevotella"= "#3C5488FF", "Blautia" = "#AD002AFF","Bacteroides" = "#00A087FF",  
                                                                                   "Bifidobacterium" = "#E64B35FF", "Faecalibacterium"= "#8491B4FF", "Alistipes" = "#FAFD7CFF",
                                                                                   "Lactobacillus"=  "#631879FF", "Enterococcus" = "#E18727FF", "Butyricicoccus" = "#FFDC91FF", 
                                                                                  "Clostridium sensu stricto 1"= "#00468BFF"))
stool.heatmap <- pheatmap(PathM, cluster_rows = F, cluster_cols = T,
                           color = colorRampPalette(c("white","#832424FF"))(100), #"#3A3A98FF",
                           border_color = NA,
                           annotation_col = col_groups, 
                           annotation_colors = colour_groups,
                           show_rownames = T,
                           show_colnames = F)

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q7_path_pred_heatmap_stool.pdf", plot = stool.heatmap, width = 8, height = 10)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q7_path_pred_heatmap_stool.png", plot = stool.heatmap, width = 8, height = 10)

####Based on the dominance make new PCo for stool 
##Bray-Curtis
PS4.stool.2<- PS4.stool
top.stool%>%
  column_to_rownames("SampleID")-> tmp

##Add antibiotic burden
antibioticB%>%
  dplyr::filter(material== "Stool")%>%
  dplyr::select(c(SampleID, AntibioticBurden_total, AntibioticBurden_iv))%>%
  dplyr::distinct()%>%
  column_to_rownames("SampleID")-> ABX.stool
  
tmp<- cbind(tmp, ABX.stool)

PS4.stool.2@sam_data<- sample_data(tmp)

BC_dist<- phyloseq::distance(PS4.stool.2,
                             method="bray", weighted=F)
ordination<- ordinate(PS4.stool.2,
                      method="PCoA", distance= BC_dist)

##Permanova
BC.test.stool2<- vegan::adonis2(BC_dist~ Phenotype_severity+ Genus + AntibioticBurden_total + sex + age +  Visit + BMI,
                                permutations = 999, data = tmp, na.action = F)

kable(BC.test.stool2$aov.tab)

## Calculate multivariate dispersion (aka distance to the centroid)
mvd<- vegan::betadisper(BC_dist, tmp.stool$Genus, type = "centroid")
mvd.perm<- vegan::permutest(mvd, permutations = 999)

plot(mvd)
anova(mvd)
boxplot(mvd)
plot(TukeyHSD(mvd))

##Extract centroids and vectors 
centroids<-data.frame(grps=rownames(mvd$centroids), data.frame(mvd$centroids))
vectors<-data.frame(group=mvd$group,data.frame(mvd$vectors))

seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("Genus","v.PCoA1","v.PCoA2","PCoA1","PCoA2")
seg.data$Genus<- NULL
##Add Line information for each sample 
tmp.stool%>%
  dplyr::select(c(Genus, Visit, Patient_number, Phenotype_severity))%>%
  dplyr::mutate(Phenotype_severity = as.factor(Phenotype_severity))%>%
  rownames_to_column(var = "SampleID")%>%
  dplyr::group_by(Genus)%>%
  unique()%>%
  column_to_rownames(var = "SampleID")%>%
  cbind(seg.data)-> seg.data

##Plot BC
ggplot() + 
  #geom_point(data=centroids[,1:4], aes(x=PCoA1,y=PCoA2, color= grps, group=grps),size=4, shape= 4) + 
  #geom_line(data=centroids[,1:4], aes(x=PCoA1,y=PCoA2), linetype = 2, color=" gray")+
  #geom_segment(data=centroids[,1:4], aes(x=centroids[1,2], y=centroids[1,3], xend=centroids[2,2], yend=centroids[2,3]), linetype = 2, color=" gray")+ 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2,shape=Phenotype_severity, fill= Genus),size=3) +
  #geom_point(size=3, aes(fill= Patient_number, shape= as.factor(Phenotype_severity)), color= "black")+
  scale_shape_manual(values = c(25, 24), labels = c("Low (ppFEV1 > 70% at V1)", "High (ppFEV1 < 70% at V1)"))+
  #scale_shape_manual(values = c(21, 24, 22))+
  scale_fill_manual(values = tax.palette)+
  scale_color_manual(values = tax.palette)+
  #stat_ellipse(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2,color= Genus))+
  labs(title = "Bray-Curtis dissimilarity",tag= "A)")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= F)+
  theme_bw()+
  theme(text = element_text(size=16))+
  labs(fill = "Dominant taxa")+
  labs(shape = "Phenotype severity")+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> A

#ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Dominant_Stool.pdf", plot = A, width = 10, height = 11)
#ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Dominant_Stool.png", plot = A, width = 10, height = 11)

##Is visit impacting differences in composition by patient? 
BC_dist.stool%>% 
  wilcox_test(BC_dist ~ Group)%>%
  adjust_pvalue(method = "BH") %>%
  add_significance()%>%
  add_xy_position(x = "Group")-> stats.test

##Save statistical analysis
#x <- stats.test
#x$groups<- NULL
#write.csv(x, "~/CF_project/exercise-cf-intervention/tables/Q3_Sample_Visit_BC.csv")

##Plot 
BC_dist.stool%>%
  ggplot(aes(x= Group, y= BC_dist))+
  geom_boxplot(color="black", alpha= 0.5)+
  geom_point(shape=21, position=position_jitter(0.2), size=3, aes(fill= Patient_number), color= "black")+
  xlab("Visit")+
  ylab("Bray-Curtis dissimilarity")+
  labs(tag= "B)", caption = get_pwc_label(stats.test), fill = "Patient")+
  scale_fill_manual(values = pal.CF)+
  theme_classic()+
  guides(fill = guide_legend(override.aes=list(shape=c(21)), ncol = 8, size= 10), color= "none")+
  theme(text = element_text(size=16), legend.position="bottom", legend.box = "horizontal")+
  scale_x_discrete(labels=c("V1_V2" =  "V1 to V2", 
                            "V2_V3" = "V2 to V3",
                            "V1_V3" = "V1 to V3"))+
  stat_pvalue_manual(stats.test, hide.ns = F, step.increase = 0.05,
                     tip.length = 0, label = "{p.adj} {p.adj.signif}")->B

C<- ggarrange(A, B, ncol=1, nrow=2)

ggsave(file = "CF_project/exercise-cf-intervention/figures/Fig4_Stool.pdf", plot = C, width = 10, height = 12)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Fig4_Stool.png", plot = C, width = 10, height = 12)

rm(A, B, C)

##Dominant taxa for stool estimate change in relative abundance 
##Bacteroides
tmp1<- as.data.frame(stool.microbiome[,"Bacteroidota-Bacteroides"])
colnames(tmp1)<- c("Bacteroidota-Bacteroides")
tmp1%>%
  dplyr::mutate(`Bacteroidota-Bacteroides`= round((`Bacteroidota-Bacteroides`/1e6)*100, digits = 2))-> tmp1

##Prevotella
tmp2<- as.data.frame(stool.microbiome[,"Bacteroidota-Prevotella"])
colnames(tmp2)<- c("Bacteroidota-Prevotella")
tmp2%>%
  dplyr::mutate(`Bacteroidota-Prevotella`= round((`Bacteroidota-Prevotella`/1e6)*100, digits = 2))-> tmp2

##Lactobacillus
tmp3<- as.data.frame(stool.microbiome[,"Firmicutes-Lactobacillus"])
colnames(tmp3)<- c("Firmicutes-Lactobacillus")
tmp3%>%
  dplyr::mutate(`Firmicutes-Lactobacillus`= round((`Firmicutes-Lactobacillus`/1e6)*100, digits = 2))-> tmp3

##Firmicutes
tmp4<- as.data.frame(stool.microbiome[,"Firmicutes-Blautia"])
colnames(tmp4)<- c("Firmicutes-Blautia")
tmp4%>%
  dplyr::mutate(`Firmicutes-Blautia`= round((`Firmicutes-Blautia`/1e6)*100, digits = 2))-> tmp4

##Bifidobacterium
tmp5<- as.data.frame(stool.microbiome[,"Firmicutes-Bifidobacterium"])
colnames(tmp5)<- c("Firmicutes-Bifidobacterium")
tmp5%>%
  dplyr::mutate(`Firmicutes-Bifidobacterium`= round((`Firmicutes-Bifidobacterium`/1e6)*100, digits = 2))-> tmp5

tmp.stool%>%
  cbind(tmp1)%>%
  cbind(tmp2)%>%
  cbind(tmp3)%>%
  cbind(tmp4)-> tmp.stool

##Estimate deltas of abundance per Visit combination 
tmp.stool%>%
  dplyr::select(c(Patient_number, Visit, `Firmicutes-Staphylococcus`))%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  dplyr::mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  dplyr::distinct()%>%
  spread(Visit, `Firmicutes-Staphylococcus`)%>%
  dplyr::mutate(V1_V2= V1 - V2)%>%
  dplyr::mutate(V1_V3= V1 - V3)%>%
  dplyr::mutate(V2_V3= V2 - V3)%>%
  dplyr::select(c(Patient_number, V1_V2, V2_V3, V1_V3))%>%
  pivot_longer(!Patient_number, names_to = "Group", values_to = "Firmicutes-Staphylococcus")%>%
  dplyr::mutate(ID= paste0(Patient_number, Group))%>%
  dplyr::ungroup()%>%
  dplyr::select(c(ID, `Firmicutes-Staphylococcus`))-> tmp1

BC_dist.stool%>%
  left_join(tmp1, by="ID")-> BC_dist.stool

tmp.stool%>%
  dplyr::select(c(Patient_number, Visit, `Proteobacteria-Pseudomonas`))%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  dplyr::mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  dplyr::distinct()%>%
  spread(Visit, `Proteobacteria-Pseudomonas`)%>%
  dplyr::mutate(V1_V2= V1 - V2)%>%
  dplyr::mutate(V1_V3= V1 - V3)%>%
  dplyr::mutate(V2_V3= V2 - V3)%>%
  dplyr::select(c(Patient_number, V1_V2, V2_V3, V1_V3))%>%
  pivot_longer(!Patient_number, names_to = "Group", values_to = "Proteobacteria-Pseudomonas")%>%
  dplyr::mutate(ID= paste0(Patient_number, Group))%>%
  dplyr::ungroup()%>%
  dplyr::select(c(ID, `Proteobacteria-Pseudomonas`))-> tmp2

BC_dist.stool%>%
  left_join(tmp2, by="ID")-> BC_dist.stool

tmp.stool%>%
  dplyr::select(c(Patient_number, Visit, `Firmicutes-Streptococcus`))%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  dplyr::mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  dplyr::distinct()%>%
  spread(Visit, `Firmicutes-Streptococcus`)%>%
  dplyr::mutate(V1_V2= V1 - V2)%>%
  dplyr::mutate(V1_V3= V1 - V3)%>%
  dplyr::mutate(V2_V3= V2 - V3)%>%
  dplyr::select(c(Patient_number, V1_V2, V2_V3, V1_V3))%>%
  pivot_longer(!Patient_number, names_to = "Group", values_to = "Firmicutes-Streptococcus")%>%
  dplyr::mutate(ID= paste0(Patient_number, Group))%>%
  dplyr::ungroup()%>%
  dplyr::select(c(ID, `Firmicutes-Streptococcus`))-> tmp3

BC_dist.stool%>%
  left_join(tmp3, by="ID")-> BC_dist.stool

tmp.stool%>%
  dplyr::select(c(Patient_number, Visit, `Proteobacteria-Stenotrophomonas`))%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  dplyr::mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  dplyr::distinct()%>%
  spread(Visit, `Proteobacteria-Stenotrophomonas`)%>%
  dplyr::mutate(V1_V2= V1 - V2)%>%
  dplyr::mutate(V1_V3= V1 - V3)%>%
  dplyr::mutate(V2_V3= V2 - V3)%>%
  dplyr::select(c(Patient_number, V1_V2, V2_V3, V1_V3))%>%
  pivot_longer(!Patient_number, names_to = "Group", values_to = "Proteobacteria-Stenotrophomonas")%>%
  dplyr::mutate(ID= paste0(Patient_number, Group))%>%
  dplyr::ungroup()%>%
  dplyr::select(c(ID, `Proteobacteria-Stenotrophomonas`))-> tmp4

BC_dist.stool%>%
  left_join(tmp4, by="ID")-> BC_dist.stool

##Model time 
BC_dist.stool%>%
  dplyr::select(ppFVC, Trainingfrequency, Trainingtime, ppFEV1, BC_dist, Patient_number, Group, 
                `Firmicutes-Staphylococcus`,`Proteobacteria-Pseudomonas`, 
                `Firmicutes-Streptococcus`, `Proteobacteria-Stenotrophomonas`)%>%
  dplyr::filter(complete.cases(.))-> tmp

##Model selection do it with glm
tmp%>%
  dplyr::select(Trainingfrequency, Trainingtime, ppFEV1, 
                `Firmicutes-Staphylococcus`,`Proteobacteria-Pseudomonas`, 
                `Firmicutes-Streptococcus`, `Proteobacteria-Stenotrophomonas`)-> tmp.ppFEV1

full.model<- glm(ppFEV1 ~ ., data = tmp.ppFEV1) ##Full model
# Stepwise regression model
step.model <- MASS::stepAIC(full.model, direction = "both", 
                            trace = FALSE)
summary(step.model)
step.model.df<- as.data.frame(coef(summary(step.model)))

step.model.df%>%
  mutate(p.adj = p.adjust(`Pr(>|t|)`, method='BH')) %>%
  add_significance()%>%
  rownames_to_column()-> step.model.df

write.csv(step.model.df, "~/CF_project/exercise-cf-intervention/tables/Q2_Model_ppFEV1_Dominant_taxa_stool.csv", row.names = F)

##Mixed effect models
##with patient and intervisit group as individual random effects

tr0<- lmer(ppFEV1 ~ Trainingfrequency + (1 | Patient_number), data = tmp) ##Null model
tr1<- lmer(ppFEV1 ~ Trainingfrequency + BC_dist + (1 | Patient_number), data = tmp) ##Null model
tr2<-lmer(ppFEV1 ~ Trainingfrequency + `Firmicutes-Staphylococcus` + (1 | Patient_number), data = tmp)

lrtest(tr0, tr1)
lrtest(tr1, tr2)
lrtest(tr0, tr2)

##Model selection do it with glm
tmp%>%
  dplyr::select(Trainingfrequency, Trainingtime, ppFVC, 
                `Firmicutes-Staphylococcus`,`Proteobacteria-Pseudomonas`, 
                `Firmicutes-Streptococcus`, `Proteobacteria-Stenotrophomonas`)-> tmp.ppFVC

full.model<- glm(ppFVC ~ ., data = tmp.ppFVC) ##Full model
# Stepwise regression model
step.model <- MASS::stepAIC(full.model, direction = "both", 
                            trace = FALSE)
summary(step.model)
step.model.df<- as.data.frame(coef(summary(step.model)))

step.model.df%>%
  mutate(p.adj = p.adjust(`Pr(>|t|)`, method='BH')) %>%
  add_significance()%>%
  rownames_to_column()-> step.model.df

write.csv(step.model.df, "~/CF_project/exercise-cf-intervention/tables/Q2_Model_ppFVC_Dominant_taxa_stool.csv", row.names = F)

##Mixed effect models
##with patient and intervisit group as individual random effects

tr0<- lmer(ppFVC ~ Trainingfrequency + (1 | Patient_number), data = tmp) ##Null model
tr1<- lmer(ppFVC ~ Trainingfrequency + BC_dist + (1 | Patient_number), data = tmp) ##Null model
tr2<-lmer(ppFVC ~ Trainingfrequency + `Firmicutes-Staphylococcus` + (1 | Patient_number), data = tmp)

lrtest(tr0, tr2)
lrtest(tr1, tr2)