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
##PS3<- rarefy_even_depth(PS2, rngseed=2020, sample.size=min(sample_sums(PS2)), replace=F)
##readcount(PS3)

## Merge ASVs that have the same taxonomy at a certain taxonomic rank (in this case Phylum and Family)
PS.Fam<-  tax_glom(PS2, "Family", NArm = F)
summarize_phyloseq(PS.Fam)

PS.Gen<-  tax_glom(PS2, "Genus", NArm = T)
summarize_phyloseq(PS.Gen)

PS.Phy<-  tax_glom(PS2, "Phylum", NArm = F)
summarize_phyloseq(PS.Phy)

plot_bar(PS.Phy, fill="Phylum") + facet_wrap(~material, scales= "free_x", nrow=1)

##Alpha diversity (not-rarefied)
##Estimate global indicators
alphaDiv <-microbiome::alpha(PS2, index = "all")
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
wh0 <- genefilter_sample(PS2, filterfun_sample(function(x) x > 5), A=0.05*nsamples(PS2))
PS3<- prune_taxa(wh0, PS2)

##Transform to an even sample size
PS3<- transform_sample_counts(PS3, function(x) 1E6 * x/sum(x))

#phylum.sum <- tapply(taxa_sums(PS3), tax_table(PS3)[, "Phylum"], sum, na.rm=TRUE)
#top5phyla <- names(sort(phylum.sum, TRUE))[1:5]
#PS3<- prune_taxa((tax_table(PS3)[, "Phylum"] %in% top5phyla), PS3)

patient<- get_variable(PS3, "Patient_number")
patient<- fct_relevel(patient, "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18")
sample_data(PS3)$Patient_number <- patient

##Weighted Unifrac
wunifrac_dist<- phyloseq::distance(PS3,
                                   method="unifrac", weighted=T)
ordination<- ordinate(PS3,
                      method="PCoA", distance=wunifrac_dist)
plot_ordination(PS3, ordination, shape= "material")+ 
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(color= Patient_number))+
  geom_point(color= "black", size= 1.5)+
  labs(title = "Weighted UniFrac",tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  labs(colour = "Patient")+
  labs(shape = "Sample type")

vegan::adonis(wunifrac_dist~ Patient_number + material + Visit,
              permutations = 999, data = sdt)

##Bray-Curtis
BC_dist<- phyloseq::distance(PS3,
                                   method="bray", weighted=T)
ordination<- ordinate(PS3,
                      method="PCoA", distance= BC_dist)
plot_ordination(PS3, ordination, shape= "material")+ 
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(color= Patient_number))+
  #geom_point(color= "black", size= 1.5)+
  labs(title = "Bray-Curtis dissimilariy",tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  labs(colour = "Patient")+
  labs(shape = "Sample type")+
  xlab("PCo1 (26.6%)")+
  ylab("PCo2 (14.6%)")-> c

vegan::adonis(BC_dist~ Patient_number + material + Visit,
              permutations = 999, data = sdt)

png("CF_project/exercise-cf-intervention/figures/Q1_Alpha_div_General.png", units = 'in', res = 300, width=14, height=14)
grid.arrange(A, B)
dev.off()

png("CF_project/exercise-cf-intervention/figures/Q1_Beta_div_General.png", units = 'in', res = 300, width=14, height=14)
C
dev.off()

##Q2: Analysis by sample type
PS3.stool<- subset_samples(PS3, material%in%c("Stool"))
PS3.sput<-  subset_samples(PS3, material%in%c("Sputum"))

plot_bar(PS3.stool, fill="Phylum")+ 
  facet_wrap(~Visit, scales= "free_x", nrow=1)+
  labs(tag= "A)")-> A

plot_bar(PS3.sput, fill="Phylum")+ 
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
BC_dist<- phyloseq::distance(PS3.stool,
                             method="bray", weighted=F)
ordination<- ordinate(PS3.stool,
                      method="PCoA", distance= BC_dist)
plot_ordination(PS3.stool, ordination, shape= "Visit")+ 
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(color= Patient_number))+
  #geom_point(color= "black", size= 1.5)+
  labs(title = "Bray-Curtis dissimilariy",tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  labs(colour = "Patient")+
  labs(shape = "Visit")+
  xlab("PCo1 (15.3%)")+
  ylab("PCo2 (10.6%)")-> D

ordination$values[1,2]
ordination$values[2,2]

png("CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Stool.png", units = 'in', res = 300, width=10, height=8)
D
dev.off()
rm(A,B,C,D)

## Differences are more linked to patient rather than Visit... but difficult to assess. 
sdt.stool<- sample_data(PS3.stool)

mvd.Stool<- vegan::betadisper(BC_dist, sdt.stool$Patient_number, type = "centroid")
vegan::permutest(mvd.Stool, permutations = 999)
anova(mvd.Stool)
plot(mvd.Stool)
boxplot(mvd.Stool)
tuky<- TukeyHSD(mvd.Stool)

tuky<- as.data.frame(tuky$group)

tuky%>%
  rename("p adj" = "p.adj")%>%
  filter(p.adj <= 0.05)%>%
  arrange(desc(diff))%>%
  rownames_to_column(var = "Pat_comp")%>%
  ggplot(aes(diff, Pat_comp))+
  geom_boxplot()

plot(TukeyHSD(mvd.Stool))

###Correlation with nutritional and respiratory activity
##Glom by genus
PS.stool.Gen<-  tax_glom(PS3.stool, "Genus", NArm = T)

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
