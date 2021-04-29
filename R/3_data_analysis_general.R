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
##For analysis with linear models
library("lme4")
library("lmtest")
library("rcompanion")
library("car")
library("merTools")

reRun<- FALSE

##Load data 
if(!exists("PS")){
  if(isTRUE(reRun)){
    source("R/2_phyloseq_preparation.R") ## Run the script at base directory of repository!   
  } else {
    PS<- readRDS(file = "~/CF_project/exercise-cf-intervention/data/PhyloSeqComp.Rds") ##New annotation SILVA
  }
}

training<- read_tsv("~/CF_project/Metadata/sample_data_indexed_training.csv")

##Discard some variables after discussion with Rebecca 
metadata<- readRDS(file = "~/CF_project/exercise-cf-intervention/data/metadata_indexed.rds")

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

#3) Remove low prevalent ASVs (This data is used for alpha diversity and differential abundance analysis)
##Remove ASVs that do not show appear more than 5 times in more than 10% of the samples
wh0 <- genefilter_sample(PS2, filterfun_sample(function(x) x > 5), A=0.01*nsamples(PS2))
PS3<- prune_taxa(wh0, PS2)
hist(readcount(PS3))

##4) Transform to even sampling depth
## Rarefy without replacement to the min sequencing depth 
#PS4<- rarefy_even_depth(PS3, rngseed=2020, sample.size=min(sample_sums(PS3)), replace=F)
##Normalization transformation to an even sample size
PS4<- transform_sample_counts(PS3, function(x) 1E6 * x/sum(x))
readcount(PS4)
##Tiny adjustment
patient<- get_variable(PS4, "Patient_number")
patient<- fct_relevel(patient, "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18")
sample_data(PS4)$Patient_number <- patient

##Count ASVs
asv.sample<- as.data.frame(PS3@otu_table)
test<- data.frame()
for (i in 1:ncol(asv.sample)) {
  asv<- data.frame()
  asv[1,1]<- sum(asv.sample[,i]!=0)
  rownames(asv)<- paste0("Sample", i)
  test <- rbind(test, asv) ### Join all the "individual" data frames into the final data frame 
}

asv.sample<- test
rm(asv, test)

summary(asv.sample$V1) #--> For the results of ASVs per sample

## Merge ASVs that have the same taxonomy at a certain taxonomic rank (in this case Phylum and Family)
PS.Fam<-  tax_glom(PS4, "Family", NArm = F)
summarize_phyloseq(PS.Fam)

PS.Gen<-  tax_glom(PS4, "Genus", NArm = T)
summarize_phyloseq(PS.Gen)

PS.Phy<-  tax_glom(PS4, "Phylum", NArm = F)
summarize_phyloseq(PS.Phy)

plot_bar(PS.Phy, fill="Phylum") + facet_wrap(~material, scales= "free_x", nrow=1)

##Alpha diversity (not-rarefied)
##Estimate global indicators
alphaDiv <-microbiome::alpha(PS3, index = "all")
kable(head(alphaDiv))

alphaDiv%>%
  rownames_to_column()%>%
  dplyr::rename(SampleID = rowname)->tmp1

as_tibble(sample_data(PS4))->tmp2

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

C<-ggarrange(A, B, ncol=1, nrow=2, common.legend = TRUE, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Alpha_div_General.pdf", plot = C, width = 10, height = 8)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Alpha_div_General.png", plot = C, width = 10, height = 8)
rm(A,B,C)

##Beta diversity
##Bray-Curtis
BC_dist<- phyloseq::distance(PS4, method="bray", weighted=T)
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
  xlab("PCo1 (24.8%)")+
  ylab("PCo2 (14.0%)")

##Q2: Analysis by sample type
PS4.stool<- subset_samples(PS4, material%in%c("Stool"))
PS4.stool.Phy<-  tax_glom(PS4.stool, "Phylum", NArm = F)

##Subset those samples with Benzose treatment (Correct Sputum samples)
PS4.sput<-  subset_samples(PS4, Benzoase==1)
PS4.sput.Phy<-  tax_glom(PS4.sput, "Phylum", NArm = F)

plot_bar(PS4.stool.Phy, x = "Patient_number", fill="Phylum")+ 
  facet_wrap(~Visit, scales= "free_x", nrow=1)+
  xlab("Patient number")+
  labs(tag= "A)")-> A

plot_bar(PS4.sput.Phy, x = "Patient_number",fill="Phylum")+ 
  facet_wrap(~Visit, scales= "free_x", nrow=1)+
  xlab("Patient number")+
  labs(tag= "B)")-> B

C<-ggarrange(A, B, ncol=1, nrow=2, common.legend = TRUE, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Composition_General.pdf", plot = C, width = 10, height = 8)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Composition_General.png", plot = C, width = 10, height = 8)
rm(A,B,C)

##Remove what is not need it in the next script to clean the global environment
rm(PS, PS1, PS2, Prevdf, wh0, asv.sample)