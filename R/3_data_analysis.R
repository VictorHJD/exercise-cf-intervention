##Cystic-Fibrosis Microbiome Project (Mainz)
##Data analysis
##Víctor Hugo Jarquín-Díaz 18.02.2021

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

##Load data 
if(!exists("PS")){
  if(isTRUE(reRun)){
    source("R/2_phyloseq_preparation.R") ## Run the script at base directory of repository!   
  } else {
    metadata<- readRDS(file = "~/CF_project/exercise-cf-intervention/data/PhyloSeqComp.Rds") ##New annotation SILVA
  }
}

##Have a look into the data
summarize_phyloseq(PS) ##Non-Normalized

## HOW many ASVs for off-target eukaryotes and archaea
table(tax_table(PS)[, "Kingdom"], exclude = NULL) ## ---> README results summary

## HOW many reads for off-target eukaryotes and archaea
by((otu_table(PS)), tax_table(PS)[, "Kingdom"], sum) ## --->  README results summary

##Negative controls are way bellow the samples, let's remove them

PS1 <- prune_samples(sample_sums(PS)>190, PS)

##Filtering 
hist(sample_sums(PS1))

##1) Taxa filtering: Remove "NA" ASVs
PS2<- subset_taxa(PS1, !is.na(Phylum) & !Phylum %in% c(""))

##2) Taxa filtering: Remove low prevalent taxa
##Create a prevalence dataframe 
Prevdf<- apply(X = otu_table(PS2),
               MARGIN = 1,
               FUN = function(x){sum(x > 0)})

##Add taxonomy and total read counts to this data.frame
Prevdf<- data.frame(Prevalence = Prevdf,
                    TotalAbundance = taxa_sums(PS2),
                    tax_table(PS2))

plyr::ddply(Prevdf, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})

ggplot(Prevdf, aes(TotalAbundance, Prevalence / nsamples(PS2),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Log 10 Total Reads") + ylab("Prevalence [Prop. of Samples]") +
  theme_bw()+
  facet_wrap(~Phylum) + theme(legend.position="none")

##3) Transform to even sampling depth
## Rarefy without replacement
PS3<- rarefy_even_depth(PS2, rngseed=2020, sample.size=min(sample_sums(PS2)), replace=F)
##readcount(PS3)

## Merge ASVs that have the same taxonomy at a certain taxonomic rank (in this case Phylum and Family)
PS.Fam<-  tax_glom(PS3, "Family", NArm = F)
summarize_phyloseq(PS.Fam)

PS.Gen<-  tax_glom(PS3, "Genus", NArm = T)
summarize_phyloseq(PS.Gen)

PS.Phy<-  tax_glom(PS3, "Phylum", NArm = F)
summarize_phyloseq(PS.Phy)

plot_bar(PS.Phy, fill="Phylum") + facet_wrap(~material, scales= "free_x", nrow=1)

##Alpha diversity (not-rarefied)
##Estimate global indicators
alphaDiv <-microbiome::alpha(PS3, index = "all")
kable(head(alphaDiv))

alphaDiv%>%
  rownames_to_column()%>%
  rename(rowname = "SampleID")->tmp1

as_tibble(sample)->tmp2

tmp1<-inner_join(tmp1, tmp2, by="SampleID")
rownames(tmp1)<- tmp1$SampleID
tmp1$SampleID<- NULL
sdt<- tmp1
rm(tmp1,tmp2)

###Alpha diversity 
##Q1: Differences in alpha diversity among days between sample types
##Richness
sdt%>% 
  mutate(Visit = fct_relevel(Visit, 
                                   "V1", "V2", "V3"))%>%
  dplyr::group_by(Visit)%>%
  wilcox_test(chao1 ~ material)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Visit")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "~/CF_project/exercise-cf-intervention/tables/Q1_Sample_Visit_Chao1.csv")

sdt%>% 
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  dplyr::group_by(Visit)%>%
  wilcox_effsize(chao1 ~ material)

##Plot 
sdt%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  dplyr::group_by(Visit)%>%
  ggplot(aes(x= Visit, y= chao1))+
  geom_boxplot(aes(color= material), alpha= 0.5)+
  geom_point(shape=21, position=position_jitter(0.2), size=3, aes(fill= material), color= "black")+
  xlab("Visit")+
  ylab("Richness (Chao1 Index)")+
  labs(tag= "A)", caption = get_pwc_label(stats.test))+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -2, hide.ns = F,label = "{p.adj}{p.adj.signif}")-> A

##Shannon diversity 
sdt%>% 
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  dplyr::group_by(Visit)%>%
  wilcox_test(diversity_shannon ~ material)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Visit")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "~/CF_project/exercise-cf-intervention/tables/Q1_Sample_Visit_Shannon.csv")

sdt%>% 
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  dplyr::group_by(Visit)%>%
  wilcox_effsize(diversity_shannon ~ material)

##Plot 
sdt%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  dplyr::group_by(Visit)%>%
  ggplot(aes(x= Visit, y= diversity_shannon))+
  geom_boxplot(aes(color= material), alpha= 0.5)+
  geom_point(shape=21, position=position_jitter(0.2), size=3, aes(fill= material), color= "black")+
  xlab("Visit")+
  ylab("Diversity (Shannon Index)")+
  labs(tag= "B)", caption = get_pwc_label(stats.test))+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_pvalue_manual(stats.test, hide.ns = F,label = "{p.adj}{p.adj.signif}")-> B

##Beta diversity
##Remove ASVs that do not show appear more than 5 times in more than half the samples
wh0 <- genefilter_sample(PS3, filterfun_sample(function(x) x > 5), A=0.05*nsamples(PS3))
PS4<- prune_taxa(wh0, PS3)

##Transform to an even sample size
##PS4<- transform_sample_counts(PS4, function(x) 1E6 * x/sum(x))
readcount(PS4)

patient<- get_variable(PS4, "Patient_number")
patient<- fct_relevel(patient, "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18")
sample_data(PS4)$Patient_number <- patient

##Weighted Unifrac
wunifrac_dist<- phyloseq::distance(PS4,
                                   method="unifrac", weighted=T)
ordination<- ordinate(PS4,
                      method="PCoA", distance=wunifrac_dist)
plot_ordination(PS4, ordination, shape= "material")+ 
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(color= Patient_number))+
  geom_point(color= "black", size= 1.5)+
  labs(title = "Weighted UniFrac",tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  labs(colour = "Patient")+
  labs(shape = "Sample type")

vegan::adonis(wunifrac_dist~ Severity + Patient_number + material + Visit,
              permutations = 999, data = sdt)

##Bray-Curtis
BC_dist<- phyloseq::distance(PS4,
                                   method="bray", weighted=T)
ordination<- ordinate(PS4,
                      method="PCoA", distance= BC_dist)
plot_ordination(PS4, ordination, shape= "material")+ 
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(color= Patient_number))+
  #geom_point(color= "black", size= 1.5)+
  labs(title = "Bray-Curtis dissimilariy",tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  labs(colour = "Patient")+
  labs(shape = "Sample type")+
  xlab("PCo1 (26.9%)")+
  ylab("PCo2 (14.7%)")-> C

vegan::adonis(BC_dist~ Severity + Patient_number + material + Visit,
              permutations = 999, data = sdt)

png("CF_project/exercise-cf-intervention/figures/Q1_Alpha_div_General.png", units = 'in', res = 300, width=14, height=14)
grid.arrange(A, B)
dev.off()

png("CF_project/exercise-cf-intervention/figures/Q1_Beta_div_General.png", units = 'in', res = 300, width=14, height=14)
C
dev.off()

rm(A,B,C)
##Q2: Analysis by sample type
PS4.stool<- subset_samples(PS4, material%in%c("Stool"))
PS4.sput<-  subset_samples(PS4, material%in%c("Sputum"))

plot_bar(PS4.stool, fill="Phylum")+ 
  facet_wrap(~Visit, scales= "free_x", nrow=1)+
  labs(tag= "A)")-> A

plot_bar(PS4.sput, fill="Phylum")+ 
  facet_wrap(~Visit, scales= "free_x", nrow=1)+
  labs(tag= "B)")-> B

png("CF_project/exercise-cf-intervention/figures/Q1_Composition_General.png", units = 'in', res = 300, width=10, height=8)
grid.arrange(A, B)
dev.off()
rm(A,B)

##Stool#####################

##Plot 
sdt%>%
  dplyr::filter(material=="Stool")%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Visit, y= chao1))+
  geom_boxplot(color= "black", outlier.colour = "white")+
  geom_point(shape=21, size=3, aes(fill= Patient_number), color= "black")+
  xlab("Visit")+
  ylab("Richness (Chao1 Index)")+
  geom_line(aes(group = Patient_number), color= "gray")+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))

##Shannon diversity 
sdt%>% 
  dplyr::filter(material=="Stool")%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  wilcox_test(diversity_shannon ~ Patient_number)%>%
  add_significance()%>%
  add_xy_position(x = "Visit")

##Plot 
sdt%>%
  dplyr::filter(material=="Stool")%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Visit, y= diversity_shannon))+
  geom_boxplot(color= "black", outlier.colour = "white")+
  geom_point(shape=21, size=3, aes(fill= Patient_number), color= "black")+
  xlab("Visit")+
  ylab("Diversity (Shannon Index)")+
  geom_line(aes(group = Patient_number), color= "gray")+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))

##Bray-Curtis
BC_dist<- phyloseq::distance(PS4.stool,
                             method="bray", weighted=F)
ordination<- ordinate(PS4.stool,
                      method="PCoA", distance= BC_dist)
plot_ordination(PS4.stool, ordination, shape= "Visit")+ 
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(color= Patient_number))+
  #geom_point(color= "black", size= 1.5)+
  labs(title = "Bray-Curtis dissimilariy",tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  labs(colour = "Patient")+
  labs(shape = "Visit")+
  xlab("PCo1 (15.6%)")+
  ylab("PCo2 (11.0%)")-> D

ordination$values[1,2]
ordination$values[2,2]

sdt%>%
  dplyr::filter(material=="Stool")-> tmp1

BC.test<- vegan::adonis(BC_dist~ Severity + Patient_number + sex + age +  Visit + BMI,
              permutations = 999, data = tmp1, na.action = F)

##Patient is the only significant factor and explains 64% of the variation.

png("CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Stool.png", units = 'in', res = 300, width=10, height=8)
D
dev.off()
rm(A,B,C,D)

## Differences are more linked to patient rather than Visit... but difficult to assess. 
##Extract pairwise distances per patient
BC_dist.stool<- as.matrix(BC_dist)
nm<- rownames(BC_dist.stool)
tmp1<- subset(BC_dist.stool, grepl("V1", nm))
tmp2<- tmp1[,grepl("V2", colnames(tmp1))]

###Correlation with nutritional and respiratory activity
##Glom by genus
PS.stool.Gen<-  tax_glom(PS4.stool, "Genus", NArm = T)

##Adjust ASV table for merging with taxa information
otu<- PS.stool.Gen@otu_table
otu<-as.data.frame(otu)
otu<- rownames_to_column(otu, var = "ASV")

##Select just the genus 
tax<- PS.stool.Gen@tax_table
tax<-as.data.frame(tax)
tax%>%
  select(Genus)-> tax
tax$Genus<-gsub(" ", "_", basename(tax$Genus))
tax<- rownames_to_column(tax, var = "ASV")

##Use genus as rownames
stool.microbiome<- plyr::join(otu, tax, by= "ASV")
stool.microbiome$ASV<- NULL
rownames(stool.microbiome)<- stool.microbiome$Genus
stool.microbiome$Genus<- NULL

##Transpose dataframe so samples are rows 
stool.microbiome<- t(stool.microbiome)

x <- log10(stool.microbiome+1) # ASV Log10 (39 samples x 122 genera)

##Select useful metrics
y<-as.data.frame(sdt.stool)

y<- y[,c(5:7,11:29)]
 
y <- as.matrix(y) # Metadata (39 samples x 22 technical, respiratory and nutritional variables)

# Cross correlate data sets
correlations <- associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)

# Or, alternatively, the same output is also available in a handy table format
correlation.table <- associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
kable(head(correlation.table))

heat(correlation.table, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05) 

ggplot(correlation.table, aes(X1, X2, group=X1)) + 
  geom_tile(aes(fill = Correlation)) +
  geom_text(aes(fill = correlation.table$Correlation, label = round(correlation.table$Correlation, 1)), size = 5) +
  scale_fill_gradientn("Spearman's \n Correlation", 
                              breaks = seq(from = -1, to = 1,  by = 0.25), 
                              colours = c("blue", "white", "red"), 
                              limits = c(-1, 1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(x = "", y = "", tag = "A)")-> A

png("CF_project/exercise-cf-intervention/figures/Q3_Correlation_Stool.png", units = 'in', res = 300, width=10, height=8)
A
dev.off()

## Bray-Curtis dissimilarity estimation among samples 
foo.matrix<- as.matrix(stool.microbiome)
foo.braycurt<- vegan::vegdist(foo.matrix, method = "bray")
as.matrix(foo.braycurt)

###Using pheatmap to include annotations 
foo.clust <- hclust(dist(foo.braycurt), method = "complete") ##Dendogram
require(dendextend)
as.dendrogram(foo.clust) %>%
  plot(horiz = TRUE)

foo.col <- cutree(tree = foo.clust, k = 2)
foo.col  <- data.frame(cluster = ifelse(test = foo.col  == 1, yes = "cluster 1", no = "cluster 2"))
foo.col$SampleID <- rownames(foo.col)

y<-as.data.frame(sdt.stool)

foo.col <- left_join(foo.col, y, by="SampleID", sort= F)

col_groups <- foo.col %>%
  select("SampleID","Patient_number","Visit","sex") ##Here It is possible to add the expected size 

row.names(col_groups)<- col_groups$SampleID

col_groups$SampleID<- NULL

colour_groups <- list( Visit= c("V1"= "#E3DAC9", "V2"= "pink","V3"= "#440154FF"),
                       Patient_number= c("P1"= "#8DD3C7","P2"= "#009999" ,"P3"= "#FFFFB3", "P4"= "#BEBADA", "P5"= "#FB8072", "P6"= "#80B1D3",
                                 "P7"="#FDB462", "P8"= "#B3DE69", "P9"="#FC4E07","P10"= "#FCCDE5", "P11"= "#D9D9D9",
                                 "P12"="#D95F02", "P13"= "#7570B3", "P14"= "#E7298A", "P15"= "#A6761D", "P16"= "#66A61E", 
                                 "P17"= "#E6AB02", "P18"= "#D0FF14"))
require(pheatmap)
require(viridis)
BCheatmap.stool <- pheatmap(foo.braycurt, 
                      color = viridis(100),
                      border_color = NA,
                      annotation_col = col_groups, 
                      annotation_colors = colour_groups,
                      show_rownames = F,
                      show_colnames = F,
                      main= "Bray-Curtis dissimilarity among stool samples")
png("CF_project/exercise-cf-intervention/figures/Q2_BCHeatmap_Stool.png", units = 'in', res = 300, width=10, height=8)
BCheatmap.stool
dev.off()

library("DESeq2"); packageVersion("DESeq2")

tmp2<- sdt.stool[,c(2,9)]
tmp2<- left_join(genotype, tmp2, "Patient_number")
tmp2[,2:4]-> tmp2
tmp2<- as.data.frame(tmp2)
rownames(tmp2)<- tmp2$SampleID
tmp2$SampleID<- NULL

tmp2<- cbind(sdt.stool, tmp2)
rownames(tmp2) <- sample(sample_names(PS3.stool))
sample_data(PS3.stool)$Severity <- tmp2$Severity

deseq.severity<- phyloseq_to_deseq2(PS3.stool, ~ Severity)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans<- apply(counts(deseq.severity), 1, gm_mean)
deseq.severity<- estimateSizeFactors(deseq.severity, geoMeans = geoMeans)
deseq.severity<- DESeq(deseq.severity, test="LRT", fitType="parametric", reduced= ~ 1)

ac.res <- results(deseq.severity, cooksCutoff = FALSE)
alpha <- 0.05
Bac.sigtab <- deseq.severity[which(deseq.severity$padj < alpha), ]
Bac.sigtab <- cbind(as(Bac.sigtab, "data.frame"), as(tax_table(PS.PA.genus)[rownames(Bac.sigtab), ], "matrix"))
rownames(Bac.sigtab) <- NULL
Bac.sigtab

# Phylum order
x <- tapply(Bac.sigtab$log2FoldChange, Bac.sigtab$Phylum, function(x) max(x))
x <- sort(x, TRUE)
Bac.sigtab$Phylum <- factor(as.character(Bac.sigtab$Phylum), levels=names(x))
# Genus order
x <- tapply(Bac.sigtab$log2FoldChange, Bac.sigtab$Genus, function(x) max(x))
x <- sort(x, TRUE)
Bac.sigtab$Genus <- factor(as.character(Bac.sigtab$Genus), levels=names(x))

ggplot(Bac.sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) +
  geom_point(shape=21, position=position_jitter(0.2), size=3, aes(fill= Phylum), color= "black") + 
  scale_fill_brewer(palette = "Set1")+
  coord_flip()+
  geom_hline(aes(yintercept = 0), color = "gray70", size = 0.6)+
  xlab("Bacteria Genus")+
  ylab("Ascaris <-- Log-2-Fold-Change --> Jejunum")+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5), text = element_text(size=16))+
  theme_bw()

######Sputum###################
##Bray-Curtis
BC_dist<- phyloseq::distance(PS3.sput,
                             method="bray", weighted=F)
ordination<- ordinate(PS3.sput,
                      method="PCoA", distance= BC_dist)
plot_ordination(PS3.sput, ordination, shape= "Visit")+ 
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(color= Patient_number))+
  #geom_point(color= "black", size= 1.5)+
  labs(title = "Bray-Curtis dissimilariy",tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  labs(colour = "Patient")+
  labs(shape = "Visit")+
  xlab("PCo1 (41.1%)")+
  ylab("PCo2 (20%)")-> A

ordination$values[1,2]
ordination$values[2,2]

sdt%>%
  dplyr::filter(material=="Sputum")-> tmp1

tmp1<- left_join(tmp1, genotype, "Patient_number")

vegan::adonis(BC_dist~ Patient_number + Severity + Visit + BMI + sex + age,
              permutations = 999, data = tmp1, na.action = F)

##Patient is the only significant factor and explains 79.8% of the variation.


png("CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Sputum.png", units = 'in', res = 300, width=10, height=8)
A
dev.off()
rm(A)

##Glom by genus
PS.sputum.Gen<- tax_glom(PS3.sput, "Genus", NArm = T)

##Adjust ASV table for merging with taxa information
otu<- PS.sputum.Gen@otu_table
otu<-as.data.frame(otu)
otu<- rownames_to_column(otu, var = "ASV")

##Select just the genus 
tax<- PS.sputum.Gen@tax_table
tax<-as.data.frame(tax)
tax%>%
  select(Genus)-> tax
tax$Genus<-gsub(" ", "_", basename(tax$Genus))
tax<- rownames_to_column(tax, var = "ASV")

##Use genus as rownames
sputum.microbiome<- plyr::join(otu, tax, by= "ASV")
sputum.microbiome$ASV<- NULL
rownames(sputum.microbiome)<- sputum.microbiome$Genus
sputum.microbiome$Genus<- NULL

##Transpose dataframe so samples are rows 
sputum.microbiome<- t(sputum.microbiome)

x <- log10(sputum.microbiome+1) # ASV Log10 (39 samples x 122 genera)

##Select useful metrics
sdt.sputum<- sample_data(PS3.sput)

y<-as.data.frame(sdt.sputum)

y<- y[,c(5:7,11:29)]

y <- as.matrix(y) # Metadata (39 samples x 22 technical, respiratory and nutritional variables)

# Cross correlate data sets
correlations <- associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)

# Or, alternatively, the same output is also available in a handy table format
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
  labs(x = "", y = "", tag = "B)")-> B

png("CF_project/exercise-cf-intervention/figures/Q3_Correlation_Sputum.png", units = 'in', res = 300, width=10, height=8)
B
dev.off()

## Bray-Curtis dissimilarity estimation among samples 
foo.matrix<- as.matrix(sputum.microbiome)
foo.braycurt<- vegan::vegdist(foo.matrix, method = "bray")
as.matrix(foo.braycurt)
###Using pheatmap to include annotations 
foo.clust <- hclust(dist(foo.braycurt), method = "complete") ##Dendogram
require(dendextend)
as.dendrogram(foo.clust) %>%
  plot(horiz = TRUE)

foo.col <- cutree(tree = foo.clust, k = 2)
foo.col  <- data.frame(cluster = ifelse(test = foo.col  == 1, yes = "cluster 1", no = "cluster 2"))
foo.col$SampleID <- rownames(foo.col)

y<-as.data.frame(sdt.sputum)

foo.col <- left_join(foo.col, y, by="SampleID", sort= F)

col_groups <- foo.col %>%
  select("SampleID","Patient_number","Visit","sex") ##Here It is possible to add the expected size 

row.names(col_groups)<- col_groups$SampleID

col_groups$SampleID<- NULL

colour_groups <- list( Visit= c("V1"= "#E3DAC9", "V2"= "pink","V3"= "#440154FF"),
                       Patient_number= c("P1"= "#8DD3C7","P2"= "#009999" ,"P3"= "#FFFFB3", "P4"= "#BEBADA", "P5"= "#FB8072", "P6"= "#80B1D3",
                                         "P7"="#FDB462", "P8"= "#B3DE69", "P9"="#FC4E07","P10"= "#FCCDE5", "P11"= "#D9D9D9",
                                         "P12"="#D95F02", "P13"= "#7570B3", "P14"= "#E7298A", "P15"= "#A6761D", "P16"= "#66A61E", 
                                         "P17"= "#E6AB02", "P18"= "#D0FF14"))
require(pheatmap)
require(viridis)
BCheatmap.sputum <- pheatmap(foo.braycurt, 
                            color = viridis(100),
                            border_color = NA,
                            annotation_col = col_groups, 
                            annotation_colors = colour_groups,
                            show_rownames = F,
                            show_colnames = F,
                            main= "Bray-Curtis dissimilarity among sputum samples")
png("CF_project/exercise-cf-intervention/figures/Q2_BCHeatmap_Sputum.png", units = 'in', res = 300, width=10, height=8)
BCheatmap.sputum
dev.off()
