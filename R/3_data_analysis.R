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
metadata$FFM_Luk<- NULL
metadata$Nutrition_Response<- NULL
metadata$FFM_Response<- NULL
metadata$pFVC_Response<- NULL

##Also in the PS object
sample_data(PS)$FFM_Luk<- NULL
sample_data(PS)$Nutrition_Response<- NULL
sample_data(PS)$FFM_Response<- NULL
sample_data(PS)$pFVC_Response<- NULL

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

#3) Remove low prevalent ASVs
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

summary(asv.sample$V1)

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

vegan::adonis(BC_dist~ Severity + sex + age +  Visit + BMI,
              permutations = 999, data = sdt, na.action = F, strata = sdt$Patient_number)

#png("CF_project/exercise-cf-intervention/figures/Q1_Alpha_div_General.png", units = 'in', res = 300, width=14, height=14)
grid.arrange(A, B)
#dev.off()

##Q2: Analysis by sample type
PS4.stool<- subset_samples(PS4, material%in%c("Stool"))
PS4.stool.Phy<-  tax_glom(PS4.stool, "Phylum", NArm = F)
#PS4.sput<-  subset_samples(PS4, material%in%c("Sputum"))

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

##Stool#####################
sdt%>%
  dplyr::filter(material=="Stool")%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))-> sdt.stool

sdt.stool%>%
  wilcox_test(diversity_shannon ~ Visit)%>%
  adjust_pvalue(method = "fdr") %>%
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
  ggplot(aes(x= Visit, y= chao1))+
  geom_violin(color= "black", outlier.colour = "white", trim = F)+
  geom_point(shape=21, size=3, aes(fill= Patient_number), color= "black", alpha= 0.5)+
  xlab("Visit")+
  ylab("Richness (Chao1 Index)")+
  geom_line(aes(group = Patient_number), color= "gray")+
  labs(tag= "A)")+
  theme_bw()+
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
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16))-> C

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
  geom_violin(color= "black", outlier.colour = "white", trim = F)+
  geom_point(shape=21, size=3, aes(fill= Patient_number), color= "black", alpha= 0.5)+
  xlab("Visit")+
  ylab("Diversity (Shannon Index)")+
  geom_line(aes(group = Patient_number), color= "gray")+
  labs(tag= "B)")+
  theme_bw()+
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
  labs(tag= "D)")+
  theme_bw()+
  theme(text = element_text(size=16))-> D

E<- ggarrange(A, B, C, D, ncol=2, nrow=2, common.legend = TRUE, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Alpha_Material.pdf", plot = E, width = 10, height = 8)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Alpha_Material.png", plot = E, width = 10, height = 8)
rm(A,B,C,D,E)

##Try linear models 
p0<- lm(chao1~ 1, data = sdt.stool, na.action = na.exclude)
p1<- lm(chao1~ Visit + Patient_number, data = sdt.stool, na.action = na.exclude)
p2<- lm(chao1~ Visit, data = sdt.stool, na.action = na.exclude)
p3<- lm(chao1~ Patient_number, data = sdt.stool, na.action = na.exclude)

lrtest(p0, p1)

p0<- lm(diversity_shannon~ 1, data = sdt.stool, na.action = na.exclude)
p1<- lm(diversity_shannon~ Visit + Patient_number, data = sdt.stool, na.action = na.exclude)
p2<- lm(diversity_shannon~ Visit, data = sdt.stool, na.action = na.exclude)
p3<- lm(diversity_shannon~ Patient_number, data = sdt.stool, na.action = na.exclude)

lrtest(p0, p1)

rm(p0, p1, p2, p3)

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
  xlab(paste0("PCo1 ", round(ordination$values[1,2]*100, digits = 2)))+
  ylab(paste0("PCo2 ", round(ordination$values[2,2]*100, digits = 2)))-> D

##Stratified for Patient number 
BC.test.stool<- vegan::adonis(BC_dist~ Phenotype_severity + Mutation_severity + sex + age +  Visit + BMI,
              permutations = 999, data = sdt.stool, na.action = F, strata = sdt.stool$Patient_number)

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

##From metadata extract responders and severity status 
metadata%>%
  dplyr::filter(material=="Stool")%>%
  group_by(Patient_number)%>%
  distinct(Patient_number, .keep_all = TRUE)%>%
  dplyr::select(c(Patient_number, Nutrition_Response, FFM_Response, pFVC_Response, Phenotype_severity, 
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
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))-> training

training%>%
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
  theme_bw()+
  theme(text = element_text(size=16))+
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
  geom_point(position=position_jitter(0.2), size=2.5, aes(shape= Group, fill= Patient_number), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  xlab("Training frequency")+
  ylab("Bray-Curtis dissimilarity")+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  geom_smooth(method = lm, se=FALSE)-> A

##Time
BC_dist.stool%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Trainingtime, y= BC_dist))+
  geom_point(position=position_jitter(0.2), size=2.5, aes(shape= Group, fill= Patient_number), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  xlab("Training time")+
  ylab("Bray-Curtis dissimilarity")+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  geom_smooth(method = lm, se=FALSE)-> B

##Weeks
BC_dist.stool%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Trainingsweeks, y= BC_dist))+
  geom_point(shape= 22, position=position_jitter(0.2), size=2.5, aes(fill= Patient_number), color= "black")+
  xlab("Training weeks")+
  ylab("Bray-Curtis dissimilarity")+
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  geom_smooth(method = lm, se=FALSE) -> C

D<-ggarrange(A, B, C, ncol=1, nrow=3, common.legend = TRUE, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Stool_Training.pdf", plot = D, width = 10, height = 12)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Stool_Training.png", plot = D, width = 10, height = 12)

rm(A,B,C,D)

###Mixed effect models 
##Check for complete cases
#BC_dist.stool%>% 
#  group_by(Patient_number)%>%
#  arrange(Group, .by_group = TRUE)%>%
#  summarise(n(), .groups = "keep")%>%
#  dplyr::rename(n = "n()")%>%
#  filter(n == 3)-> Keep

#Keep<- Keep$Patient_number

##Select just patients in Keep
#BC_dist.stool[BC_dist.stool$Patient_number%in%Keep, ]-> x

BC_dist.stool%>%
dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                           "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                           "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))-> BC_dist.stool
##qqPlots
qqPlot(BC_dist.stool$BC_dist)
qqPlot(BC_dist.stool$Trainingfrequency)
qqPlot(BC_dist.stool$Trainingtime)

##GLM (Training weeks was omited due to a lot of NAs)
tr0<- glm(BC_dist ~ 1, data = BC_dist.stool, na.action = na.exclude) ##Null model
tr1<- glm(BC_dist ~ Trainingfrequency, data =BC_dist.stool, na.action = na.exclude)
tr2<- glm(BC_dist ~ Trainingtime, data = BC_dist.stool, na.action = na.exclude)
tr3<- glm(BC_dist ~ Trainingfrequency + Trainingtime, data = BC_dist.stool, na.action = na.exclude)
tr4<- glm(BC_dist ~ Trainingfrequency*Trainingtime , data =BC_dist.stool, na.action = na.exclude) ##Full model

##Comparisons between models
lrtest(tr1, tr2) ##Significant difference Training Frequency seems to be a better predictor for microbial differences among visits
lrtest(tr1, tr3) ##Not significant difference 
lrtest(tr2, tr3) ##Not significant difference 
lrtest(tr1, tr4) ##No difference between model with training frequency alone and its interaction with time
lrtest(tr2, tr4) ##No difference between model with training time alone and its interaction with frequency
lrtest(tr3, tr4) ##No difference between model with training time and frequency and their interaction 

##Mixed effect models
##with patient as random effect
tr5<-lmer(BC_dist ~ Trainingfrequency * Trainingtime + (1 | Patient_number), data = BC_dist.stool)
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

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Stool_Training_Effect_Ranges.png", plot = A, width = 8, height = 8)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Stool_Training_Effect_Ranges.pdf", plot = A, width = 8, height = 8)

rm(A)

##Nutritional, pFVC, FFM response with patient random effect
##Adjust to binary 1 responders = Improvement, 0 No responders= Stable, decreased or not detected
BC_dist.stool%>%
  dplyr::mutate(Nutrition_Response = case_when(Nutrition_Response == 1  ~ 1,
                                        Nutrition_Response %in% c(0, 2,3) ~ 0))%>%
  dplyr::mutate(FFM_Response = case_when(FFM_Response == 1  ~ 1,
                                         FFM_Response %in% c(0, 2,3) ~ 0))%>%
  dplyr::mutate(pFVC_Response = case_when(pFVC_Response == 1  ~ 1,
                                          pFVC_Response %in% c(0, 2,3) ~ 0))-> BC_dist.stool

##Does the nutrition response is linked to changes in microbial composition?
p1<- glm(Nutrition_Response ~ BC_dist, data = BC_dist.stool, na.action = na.exclude, 
         family = binomial()) ## NS

##Does improvement on fat free mass response is linked to changes in microbial composition?
p2<- glm(FFM_Response ~ BC_dist, data = BC_dist.stool, na.action = na.exclude,
         family = binomial()) ## NS

##Does improvement in predicted forced vital capacity linked to changes in microbiota or training?
p3<- glm(pFVC_Response ~ BC_dist, data = BC_dist.stool, na.action = na.exclude,
         family = binomial())
p4<- glm(pFVC_Response ~ Trainingfrequency, data = BC_dist.stool, na.action = na.exclude,
         family = binomial())
p5<- glm(pFVC_Response ~ Trainingtime, data = BC_dist.stool, na.action = na.exclude,
         family = binomial())
p6<- glm(pFVC_Response ~ BC_dist + Trainingfrequency, data = BC_dist.stool, na.action = na.exclude,
         family = binomial())
p7<- glm(pFVC_Response ~ BC_dist + Trainingtime, data = BC_dist.stool, na.action = na.exclude,
         family = binomial())
p8<- glm(pFVC_Response ~ Trainingfrequency + Trainingtime, data = BC_dist.stool, na.action = na.exclude,
         family = binomial())
p9<- glm(pFVC_Response ~ BC_dist + Trainingfrequency + Trainingtime, data = BC_dist.stool, na.action = na.exclude,
         family = binomial())
p10<- glm(pFVC_Response ~ BC_dist*Trainingfrequency*Trainingtime, data = BC_dist.stool, na.action = na.exclude,
          family = binomial())

###Naive correlation with nutritional and respiratory activity
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
  dplyr::select(Genus)-> tax
tax$Genus<-gsub(" ", "_", basename(tax$Genus))
tax<- rownames_to_column(tax, var = "ASV")

##Use genus as rownames
stool.microbiome<- plyr::join(otu, tax, by= "ASV")
stool.microbiome$ASV<- NULL
rownames(stool.microbiome)<- stool.microbiome$Genus
stool.microbiome$Genus<- NULL

##Transpose dataframe so samples are rows 
stool.microbiome<- t(stool.microbiome)
#saveRDS(stool.microbiome, "~/CF_project/exercise-cf-intervention/data/Stool_rare_ASV.rds")#--> For MetadeconfoundR

x <- log10(stool.microbiome+1) # ASV Log10 (39 samples x 205 genera)

##Select useful metrics
y<-sdt.stool

y<- y[,c(26:27,29:98)] ##Eliminate non-numeric variables

##Make an adjustment to make visit and patient numeric
y$Visit<- as.numeric(gsub("V", "\\1", y$Visit))
y$Patient_number<- as.numeric(gsub("P", "\\1", y$Patient_number))

##Transform severity to binary 
y%>%
  mutate(Phenotype_severity = case_when(Phenotype_severity == 2  ~ 1,
                                        Phenotype_severity == 1 ~ 0))%>%
  mutate(Mutation_severity = case_when(Mutation_severity == 2  ~ 1,
                                        Mutation_severity == 1 ~ 0))-> y


y <- as.matrix(y) # Metadata (39 samples x 68 technical, respiratory and nutritional variables)
#saveRDS(y, "~/CF_project/exercise-cf-intervention/data/Stool_rare_Metadata.rds")#--> For MetadeconfoundR

# Cross correlate data sets, output is also available in a handy table format
correlation.table <- associate(y, x, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
kable(head(correlation.table))

heat((correlation.table), "X1", "X2", fill = "Correlation", star = "p.adj")

ggplot(correlation.table, aes(X1, X2, group=X1)) + 
  geom_tile(aes(fill = Correlation)) +
  #geom_text(aes(fill = correlation.table$Correlation, label = round(correlation.table$Correlation, 1)), size = 5) +
  scale_fill_gradientn("Spearman's \n Correlation", 
                              breaks = seq(from = -1, to = 1,  by = 0.25), 
                              colours = c("blue", "white", "red"), 
                              limits = c(-1, 1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(x = "", y = "", tag = "A)")

###Let's run this with metadeconfoundR (other script for that)

###Deseq2 analysis
library("DESeq2"); packageVersion("DESeq2")

##Get raw data to run this
PS2.stool<- subset_samples(PS2, material%in%c("Stool"))
PS2.stool<- tax_glom(PS2.stool, "Genus")

##Severity classification based on:
#genotype
deseq.severity<- phyloseq_to_deseq2(PS2.stool, ~ Mutation_severity)

##Phenotype
#deseq.severity<- phyloseq_to_deseq2(PS2.stool, ~ Phenotype_severity)

# calculate geometric means prior to estimate size factors
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans<- apply(counts(deseq.severity), 1, gm_mean)
deseq.severity<- estimateSizeFactors(deseq.severity, geoMeans = geoMeans)
deseq.severity<- DESeq(deseq.severity, fitType="local")

ac.res <- results(deseq.severity)
##Remove NA
ac.res <- ac.res[order(ac.res$padj, na.last=NA), ]

##Select cut-off value
alpha <- 0.05
Bac.sigtab <- ac.res[(ac.res$padj < alpha), ]

Bac.sigtab <- cbind(as(Bac.sigtab, "data.frame"), as(tax_table(PS2.stool)[rownames(Bac.sigtab), ], "matrix"))

##Adjust value
Bac.posigtab <- Bac.sigtab[Bac.sigtab[, "log2FoldChange"] > -10, ]
Bac.posigtab <- Bac.posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

##Save table 
#write.csv(Bac.posigtab, "~/CF_project/exercise-cf-intervention/tables/Q5_DeSeq2_Abund_Stool_Severity.csv") #Genotype
#write.csv(Bac.posigtab, "~/CF_project/exercise-cf-intervention/tables/Q5_DeSeq2_Abund_Stool_Phen_Severity.csv") #Phenotype

###Negative Binomial in Microbiome Differential Abundance Testing (plot)
ggplot(na.omit(Bac.posigtab), aes(x=reorder(Genus, -log2FoldChange), y=log2FoldChange, color=Phylum)) +
  geom_point(shape=21, position=position_jitter(0.2), size=3, aes(fill= Phylum), color= "black") + 
  scale_fill_brewer(palette = "Set1")+
  coord_flip()+
  geom_hline(aes(yintercept = 15), color = "gray70", size = 0.6)+
  xlab("ASVs Genus-level")+
  ylab("Mild <-- Log-2-Fold-Change --> Severe")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5), text = element_text(size=16), 
        axis.text.y = element_text(face="italic", color="black"))+
  labs(tag = "A)")-> A #Change to B when Phenotype

#png("CF_project/exercise-cf-intervention/figures/Q5_DefAbund_Stool_Severity.png", units = 'in', res = 300, width=20, height=10)
ggarrange(A, B, ncol=2, nrow=1, common.legend = F, legend="right")
#dev.off()

###Analysis by time points
###Make Phloseq subsets to run the analysis
##V1V2
PS2.stool12<- subset_samples(PS2.stool, Visit%in%c("V1", "V2"))
##V2V3
PS2.stool23<- subset_samples(PS2.stool, Visit%in%c("V2", "V3"))
##V1V3
PS2.stool13<- subset_samples(PS2.stool, Visit%in%c("V1", "V3"))

deseq.visit<- phyloseq_to_deseq2(PS2.stool13, ~ Visit)

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
Bac.sigtab <- cbind(as(Bac.sigtab, "data.frame"), as(tax_table(PS2.stool12)[rownames(Bac.sigtab), ], "matrix")) ##Change based on the visit

##Adjust value
Bac.posigtab <- Bac.sigtab[Bac.sigtab[, "log2FoldChange"] > -30, ]
Bac.posigtab <- Bac.posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

#write.csv(Bac.posigtab, "~/CF_project/exercise-cf-intervention/tables/Q5_DeSeq2_Abund_Stool_V1V2.csv") 
#write.csv(Bac.posigtab, "~/CF_project/exercise-cf-intervention/tables/Q5_DeSeq2_Abund_Stool_V2V3.csv") 
#write.csv(Bac.posigtab, "~/CF_project/exercise-cf-intervention/tables/Q5_DeSeq2_Abund_Stool_V1V3.csv") 

###Negative Binomial in Microbiome Differential Abundance Testing (plot)
ggplot(na.omit(Bac.posigtab), aes(x=reorder(Genus, -padj), y=log2FoldChange, color=Phylum)) +
  geom_point(shape=21, position=position_jitter(0.2), size=3, aes(fill= Phylum), color= "black") + 
  scale_fill_brewer(palette = "Set1")+
  coord_flip()+
  geom_hline(aes(yintercept = 0), color = "gray70", size = 0.6)+
  xlab("ASVs Genus-level")+
  ylab("V1 <-- Log-2-Fold-Change --> V3")+ ##Change this based on the visit
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5), text = element_text(size=16), 
        axis.text.y = element_text(face="italic", color="black"))+
  labs(tag = "C)")->C

##Save just when the three objects are in the environment
#png("CF_project/exercise-cf-intervention/figures/Q5_DefAbund_Stool_Visit.png", units = 'in', res = 300, width=10, height=13)
ggarrange(A, B, C, ncol=1, nrow=3, common.legend = F, legend="right")
#dev.off()

rm(A,B,C)

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
  geom_point(size=3, aes(color= Patient_number))+
  #geom_point(color= "black", size= 1.5)+
  labs(title = "Bray-Curtis dissimilariy",tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  labs(colour = "Patient")+
  labs(shape = "Visit")+
  xlab(paste0("PCo1 ", round(ordination$values[1,2]*100, digits = 2)))+
  ylab(paste0("PCo2 ", round(ordination$values[2,2]*100, digits = 2)))-> A

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
  dplyr::select(c(Patient_number, Nutrition_Response, FFM_Response, pFVC_Response, Phenotype_severity, 
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
  theme_bw()+
  theme(text = element_text(size=16))+
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
  geom_point(position=position_jitter(0.2), size=2.5, aes(shape= Group, fill= Patient_number), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  xlab("Training frequency")+
  ylab("Bray-Curtis dissimilarity")+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  geom_smooth(method = lm, se=FALSE)-> A

##Time
BC_dist.sputum%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Trainingtime, y= BC_dist))+
  geom_point(position=position_jitter(0.2), size=2.5, aes(shape= Group, fill= Patient_number), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  xlab("Training time")+
  ylab("Bray-Curtis dissimilarity")+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  geom_smooth(method = lm, se=FALSE)-> B

##Weeks
BC_dist.sputum%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Trainingsweeks, y= BC_dist))+
  geom_point(shape= 22, position=position_jitter(0.2), size=2.5, aes(fill= Patient_number), color= "black")+
  xlab("Training weeks")+
  ylab("Bray-Curtis dissimilarity")+
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  geom_smooth(method = lm, se=FALSE) -> C

D<-ggarrange(A, B, C, ncol=1, nrow=3, common.legend = TRUE, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Sputum_Training.pdf", plot = D, width = 10, height = 12)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Sputum_Training.png", plot = D, width = 10, height = 12)

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