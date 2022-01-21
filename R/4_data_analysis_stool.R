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

##Functions 
##Find dominant taxa per samples
find.top.asv <- function(x, taxa, num){
  require(phyloseq)
  require(magrittr)
  
  top.taxa <- tax_glom(x,taxa)
  otu <- as(otu_table(top.taxa), "matrix")
  # transpose if necessary
  if(taxa_are_rows(top.taxa)){otu <- t(otu)}
  otu <- otu_table(otu, taxa_are_rows = F)
  tax <- tax_table(top.taxa)
  # Coerce to data.frame
  n <- as.data.frame(tax)
  n%>%
    rownames_to_column()%>%
    dplyr::rename(ASV = rowname)-> n
  
  j1 <- apply(otu,1,sort,index.return=T, decreasing=T) # modifying which.max to return a list of sorted index
  j2 <- lapply(j1,'[[',"x") # select for Names
  
  m <- data.frame(unlist(j2))
  m%>%
    rownames_to_column()%>%
    dplyr::filter(unlist.j2.!=0)%>%
    separate(rowname, c("SampleID", "ASV"))%>%
    dplyr::group_by(SampleID)%>%
    slice_max(order_by = unlist.j2., n = num)%>%
    dplyr::rename(Abundance = unlist.j2.)%>%
    dplyr::mutate(Abundance = (Abundance/1E6)*100)%>%
    left_join(n, by="ASV")->m
  
  rm(top.taxa, otu, tax, j1, j2, n)
  return(m)
}

##Get data frame for bar plot at genus level 
count.high.genus <- function(x, num){
  require(phyloseq)
  require(magrittr)
  #x is a phyloseq object glomed to Genus
  #num is the threshold of Relative abundance desired 
  otu <- as(otu_table(x), "matrix")
  # transpose if necessary
  if(taxa_are_rows(x)){otu <- t(otu)}
  otu <- otu_table(otu, taxa_are_rows = F)
  tax <- tax_table(x)
  # Coerce to data.frame
  n <- as.data.frame(tax)
  n%>%
    rownames_to_column()%>%
    dplyr::rename(ASV = rowname)-> n
  
  j1 <- apply(otu,1,sort,index.return=T, decreasing=T) # modifying which.max to return a list of sorted index
  j2 <- lapply(j1,'[[',"x") # select for Names
  
  m <- data.frame(unlist(j2))
  
  m%>%
    rownames_to_column()%>%
    dplyr::filter(unlist.j2.!=0)%>%
    separate(rowname, c("SampleID", "ASV"))%>%
    dplyr::group_by(SampleID)%>%
    dplyr::rename(Abundance = unlist.j2.)%>%
    dplyr::mutate(Abundance = (Abundance/1E6)*100)%>%
    left_join(n, by="ASV")%>%
    mutate(Main_taxa= Abundance>= num)%>%
    dplyr::mutate(Genus= case_when(Main_taxa== FALSE ~ "Taxa less represented", TRUE ~ as.character(Genus)))%>%
    arrange(SampleID, desc(Genus))->m
  
  m$Genus[is.na(m$Genus)]<- "Unassigned" ##Change NA's into Unassigned 
  m$Species<- NULL
  
  rm(otu, tax, j1, j2, n)
  return(m)
}

count.genus <- function(x, num){
  require(phyloseq)
  require(magrittr)
  #x is a phyloseq object glomed to Genus
  #num is the threshold of Relative abundance desired 
  otu <- as(otu_table(x), "matrix")
  # transpose if necessary
  if(taxa_are_rows(x)){otu <- t(otu)}
  otu <- otu_table(otu, taxa_are_rows = F)
  tax <- tax_table(x)
  # Coerce to data.frame
  n <- as.data.frame(tax)
  n%>%
    rownames_to_column()%>%
    dplyr::rename(ASV = rowname)-> n
  
  j1 <- apply(otu,1,sort,index.return=T, decreasing=T) # modifying which.max to return a list of sorted index
  j2 <- lapply(j1,'[[',"x") # select for Names
  
  m <- data.frame(unlist(j2))
  
  m%>%
    rownames_to_column()%>%
    dplyr::filter(unlist.j2.!=0)%>%
    separate(rowname, c("SampleID", "ASV"))%>%
    dplyr::group_by(SampleID)%>%
    dplyr::rename(Abundance = unlist.j2.)%>%
    dplyr::mutate(Abundance = (Abundance/1E6)*100)%>%
    left_join(n, by="ASV")%>%
    mutate(Main_taxa= Abundance>= num)%>%
    dplyr::mutate(Type= case_when(Main_taxa== FALSE ~ "Satellites", TRUE ~ "Colonizers"))%>%
    arrange(SampleID, desc(Genus))->m
  
  m$Genus[is.na(m$Genus)]<- "Unassigned" ##Change NA's into Unassigned 
  m$Species<- NULL
  
  rm(otu, tax, j1, j2, n)
  return(m)
}

##Transform abundance into relative abundance
Rel.abund_fun <- function(df){
  df2 <- sapply(df, function(x) (x/1E6)*100)  
  colnames(df2) <- colnames(df)
  rownames(df2) <- rownames(df)
  df2<- as.data.frame(df2)
  return(df2)
}

##Color palette for patients ##
pal.CF<- c((brewer.pal(n = 8, name = "Dark2")), (brewer.pal(n = 11, name = "Paired")) )

pal.CF<- c("P1"="#1B9E77","P2"= "#D95F02","P3"= "#7570B3","P4"= "#E7298A","P5"= "#66A61E",
           "P6"="#E6AB02","P7"= "#A6761D","P8"= "#666666","P9"= "#A6CEE3","P10"= "#1F78B4",
           "P11"= "#B2DF8A","P12"= "#33A02C","P13"= "#FB9A99","P14"="#E31A1C","P15"= "#FDBF6F",
           "P16"= "#FF7F00","P17"= "#CAB2D6","P18"= "#6A3D9A","P19"= "#FFFF99")

##For taxa 
tax.palette<- c("Taxa less represented" = "#767676FF",  "Unassigned"="lightgray", "Prevotella"= "#3C5488FF", "Blautia" = "#AD002AFF",
                "Bacteroides" = "#00A087FF",  "Bifidobacterium" = "#E64B35FF", "Subdoligranulum"= "#F39B7FFF", "Faecalibacterium"= "#8491B4FF",   
                "Collinsella"= "#CD534CFF", "Alistipes" = "#FAFD7CFF", "Holdemanella"= "#7E6148FF","Lactobacillus"=  "#631879FF",
                "Ruminococcus"= "#BC3C29FF","Escherichia-Shigella" = "#0072B5FF", "Enterococcus" = "#E18727FF", "Roseburia"= "#E762D7FF",                  
                "Acidaminococcus"= "#7876B1FF", "Intestinibacter"="#6F99ADFF", "Butyricicoccus" = "#FFDC91FF", 
                "[Ruminococcus] gnavus group" = "#EE4C97FF","Clostridium sensu stricto 1"= "#00468BFF", "Agathobacter" = "#0099B4FF" , 
                "Lachnoclostridium"= "#42B540FF", "Pediococcus"= "#925E9FFF", "Streptococcus"= "#925E9FFF", "Staphylococcus"= "#008B45FF", 
                "Rothia" ="#FDAF91FF","Alloprevotella"= "#4DBBD5FF", "Veillonella"= "#3B4992FF", "Stenotrophomonas"= "#B09C85FF", 
                "Porphyromonas"= "#BB0021FF", "Granulicatella"= "#A20056FF","Gemella" = "#0073C2FF", 
                "Achromobacter"= "#EFC000FF" , "Pseudomonas" = "#ED0000FF", "Actinomyces" = "#91D1C2FF", 
                "Haemophilus" ="#7AA6DCFF"  ,  "Lachnoanaerobaculum"= "#003C67FF",   "Fusobacterium"= "#8F7700FF", 
                "Campylobacter"= "#A73030FF",  "Neisseria"= "#3B3B3BFF")

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
  labs(tag= "B)", fill= "Patient number")+
  theme_classic()+
  theme(text = element_text(size=16))-> B

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
  labs(tag= "C)", fill= "Patient number")+
  theme_classic()+
  theme(text = element_text(size=16))-> C

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

#E<- ggarrange(A, B, C, D, ncol=2, nrow=2, common.legend = TRUE, legend="right")

#ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Alpha_Material.pdf", plot = E, width = 10, height = 8)
#ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Alpha_Material.png", plot = E, width = 10, height = 8)
#rm(A,B,C,D,E)

##Dominance 
sdt%>%
  dplyr::filter(material=="Stool")%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Visit, y= dominance_dbp))+
  geom_violin(color= "black", outlier.colour = "white", trim = F)+
  geom_point(shape=21, size=3, aes(fill= Patient_number), color= "black", alpha= 0.5)+
  xlab("Visit")+
  ylab("Dominance (Berger-Parker Index)")+
  geom_line(aes(group = Patient_number), color= "gray")+
  scale_fill_manual(values = pal.CF)+
  labs(tag= "E)", fill= "Patient number")+
  theme_classic()+
  theme(text = element_text(size=16))-> E

sdt%>%
  dplyr::filter(material=="Sputum")%>%
  dplyr::filter(Benzoase==1)%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Visit, y= dominance_dbp))+
  geom_violin(color= "black", outlier.colour = "white", trim = F)+
  geom_point(shape=21, size=3, aes(fill= Patient_number), color= "black", alpha= 0.5)+
  xlab("Visit")+
  ylab("Dominance (Berger-Parker Index)")+
  geom_line(aes(group = Patient_number), color= "gray")+
  scale_fill_manual(values = pal.CF)+
  labs(tag= "F)", fill= "Patient number")+
  theme_classic()+
  theme(text = element_text(size=16))-> F1

##Evenness
sdt%>%
  dplyr::filter(material=="Stool")%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Visit, y= evenness_pielou))+
  geom_violin(color= "black", outlier.colour = "white", trim = F)+
  geom_point(shape=21, size=3, aes(fill= Patient_number), color= "black", alpha= 0.5)+
  xlab("Visit")+
  ylab("Evenness (Pielou Index)")+
  geom_line(aes(group = Patient_number), color= "gray")+
  scale_fill_manual(values = pal.CF)+
  labs(tag= "G)", fill= "Patient number")+
  theme_classic()+
  theme(text = element_text(size=16))-> G

sdt%>%
  dplyr::filter(material=="Sputum")%>%
  dplyr::filter(Benzoase==1)%>%
  mutate(Visit = fct_relevel(Visit, 
                             "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Visit, y= evenness_pielou))+
  geom_violin(color= "black", outlier.colour = "white", trim = F)+
  geom_point(shape=21, size=3, aes(fill= Patient_number), color= "black", alpha= 0.5)+
  xlab("Visit")+
  ylab("Evenness (Pielou Index)")+
  geom_line(aes(group = Patient_number), color= "gray")+
  scale_fill_manual(values = pal.CF)+
  labs(tag= "H)", fill= "Patient number")+
  theme_classic()+
  theme(text = element_text(size=16))-> H

I<- ggarrange(A, B, C, D, E, F1, G, H, ncol = 2, nrow = 4, common.legend = TRUE, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Alpha_Material_Final.pdf", plot = I, width = 10, height = 13, dpi = 600)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Alpha_Material_Final.png", plot = I, width = 10, height = 13, dpi = 600)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Alpha_Material_Final.svg", plot = I, width = 10, height = 13, dpi = 600)

rm(A, B, C, D, E, F1, G, H, I)

##Barplot by sample 
gen.stool<- count.high.genus(x = PS4.stool.Gen, num = 10)

sdt.stool%>%
  rownames_to_column()%>%
  dplyr::select(c(1, 32:33))%>%
  dplyr::rename(SampleID = rowname)%>%
  left_join(gen.stool, by="SampleID")-> gen.stool

#plot
gen.stool%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x=Patient_number, y=Abundance, fill=Genus))+ 
  geom_bar(aes(), stat="identity", position="stack") + 
  facet_wrap(~Visit, scales= "free_x", nrow=1)+
  scale_fill_manual(values=tax.palette) + 
  guides(fill=guide_legend(nrow=5))+
  theme_bw()+
  labs(tag= "A)")+
  ylab("Relative abundance (%)")+
  xlab("Patient number")+
  theme(legend.position="bottom")#-> A

###Colonizers and satellites
Type.stool<- count.genus(x = PS4.stool.Gen, num = 4.5)

sdt.stool%>%
  rownames_to_column()%>%
  dplyr::select(c(1, 32:33))%>%
  dplyr::rename(SampleID = rowname)%>%
  left_join(Type.stool, by="SampleID")-> Type.stool

saveRDS(Type.stool, "CF_project/exercise-cf-intervention/data/Type.stool.rds")

Type.stool%>%
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

###Diversity and lung function 
sdt%>%
  dplyr::filter(material=="Stool")%>%
  mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= diversity_shannon, y= ppFEV1, shape= Visit))+
  geom_point(size=3, aes(fill= Patient_number), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=70),alpha=0.05,fill="grey")+
  geom_smooth(method=lm, se = F, color= "black")+
  scale_fill_manual(values = pal.CF)+
  theme_bw()+
  labs(tag= "A)")+
  xlab("Alpha diveristy (Shannon Index)")+
  ylab("Lung function (ppFEV1)")+
  geom_hline(yintercept=70, linetype="dashed", color = "red")+
  stat_cor(method = "spearman", label.x = 2, label.y = 85)+ # Add sperman`s correlation coefficient
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
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=70),alpha=0.01,fill="grey")+
  geom_smooth(method=lm, se = F, color= "black")+
  theme_bw()+
  labs(tag= "B)")+
  labs(fill = "Patient")+
  labs(shape = "Visit")+
  guides(fill = guide_legend(override.aes=list(shape=c(21)), ncol = 6), shape= guide_legend(nrow = 3))+
  xlab("Alpha diveristy (Shannon Index)")+
  ylab("Lung function (ppFEV1)")+
  geom_hline(yintercept=70, linetype="dashed", color = "red")+
  stat_cor(method = "spearman", label.x = 2, label.y = 30)+ # Add sperman`s correlation coefficient
  theme(text = element_text(size=16), legend.position="bottom", legend.box = "horizontal")-> B

A+
  xlab(NULL)-> A

C<- grid.arrange(A,B, heights = c(3, 2))

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Alpha_Lung_Stool.pdf", plot = C, width = 8, height = 10)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Alpha_Lung_Stool.png", plot = C, width = 8, height = 10)

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

###Dominance and lung function 
sdt%>%
  dplyr::filter(material=="Stool")%>%
  mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= dominance_dbp, y= ppFEV1, shape= Visit))+
  geom_point(size=3, aes(fill= Visit), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  geom_smooth(method=lm, se = T, aes(color= Visit))+
  theme_bw()+
  labs(tag= "A)")+
  xlab("Dominance (Berger-Parker Index)")+
  ylab("Lung function (ppFEV1)")+
  geom_hline(yintercept=70, linetype="dashed", color = "red")+
  stat_cor(method = "spearman", label.x = 0, label.y = 30)+ # Add sperman`s correlation coefficient
  theme(text = element_text(size=16), legend.position = "none")+
  facet_grid(rows = vars(Visit))-> A

sdt%>%
  dplyr::filter(material=="Stool")%>%
  mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= dominance_dbp, y= ppFEV1))+
  geom_point(size=3, aes(fill= Patient_number, shape= Visit), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  scale_fill_manual(values = pal.CF)+
  geom_smooth(method=lm, se = T, color= "black")+
  theme_bw()+
  labs(tag= "B)")+
  labs(fill = "Patient")+
  labs(shape = "Visit")+
  guides(fill = guide_legend(override.aes=list(shape=c(21)), ncol = 6), shape= guide_legend(nrow = 3))+
  xlab("Dominance (Berger-Parker Index)")+
  ylab("Lung function (ppFEV1)")+
  geom_hline(yintercept=70, linetype="dashed", color = "red")+
  stat_cor(method = "spearman", label.x = 0, label.y = 30)+ # Add sperman`s correlation coefficient
  theme(text = element_text(size=16), legend.position="bottom", legend.box = "horizontal")-> B

C<- grid.arrange(A,B, heights = c(3, 2))

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Dominance_Lung_Stool.pdf", plot = C, width = 10, height = 10)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Dominance_Lung_Stool.png", plot = C, width = 10, height = 10)

rm(A,B,C)

###Extract dominant taxa per sample
top.stool<- find.top.asv(PS4.stool, "Genus", 1)
top.stool$Species<- NULL

##Add sample data
sdt.stool%>%
  rownames_to_column()%>%
  dplyr::select(c(1, 3, 6, 10, 14, 32:95))%>%
  dplyr::rename(SampleID = rowname)%>%
  left_join(top.stool, by="SampleID")-> top.stool

top.stool %>% 
  count(Genus)-> tmp

top.stool%>%
  left_join(tmp, by="Genus")-> top.stool

##Plot Lung function and dominant bugs
top.stool%>%
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
  scale_x_continuous(breaks = c(20, 40, 60, 80, 100),limits = c(20, 100))+
  facet_wrap(~Visit, scales= "free_x", nrow=1)+
  ylab("Dominant bacterial genus")+
  xlab("Lung function (ppFEV1)")+
  geom_vline(xintercept=70, linetype="dashed", color = "red")+
  theme(text = element_text(size=16), legend.position="bottom", legend.box = "horizontal",
        axis.text.y = element_text(size = 9, face="italic", color="black"))-> A

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Dominant_Lungfunct_Stool.pdf", plot = A, width = 10, height = 8)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q1_Dominant_Lungfunct_Stool.png", plot = A, width = 10, height = 8)

##Make a heatmap
tmp<- as.data.frame(otu_table(tax_glom(PS4.stool, "Genus")))

tmp1<- count.genus(PS4.stool.Gen, 1)

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

library("pheatmap")
library(dendextend)
###In order to add the annotations in good order, 
#it is necessary to have the same order in the intersection matrix and in the annotation table
P.clust <- hclust(dist(t(tmp)), method = "complete") ##Dendogram

as.dendrogram(P.clust) %>%
  plot(horiz = T)

P.col <- cutree(tree = P.clust, k = 2)
P.col  <- data.frame(cluster = ifelse(test = P.col  == 1, yes = "cluster 1", no = "cluster 2"))
P.col$SampleID <- rownames(P.col)
P.col$SampleID<- NULL
P.col<- cbind(P.col, top.stool)

col_groups <- P.col %>%
  mutate(Visit = fct_relevel(Visit, "V1", "V2", "V3"))%>%
  mutate(Patient_number = fct_relevel(Patient_number, 
                                      "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                      "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  dplyr::select(c(SampleID, Patient_number, Visit, ppFEV1, Genus)) ##Here It is possible to add the other characteristics

row.names(col_groups)<- col_groups$SampleID

col_groups$SampleID<- NULL

colour_groups <- list(Patient_number= pal.CF, Genus= tax.palette)

stool.heatmap <- pheatmap(tmp, cluster_rows = F, cluster_cols = T,
                        color = colorRampPalette(c("white","#832424FF"))(100), #"#3A3A98FF",
                        border_color = NA,
                        annotation_col = col_groups, 
                        annotation_colors = colour_groups,
                        show_rownames = T,
                        show_colnames = F)

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

BC.test.stool<- vegan::adonis2(BC_dist~ Phenotype_severity + sex + age +  Visit + BMI,
                              permutations = 999, data = tmp, na.action = F)

kable(BC.test.stool$aov.tab)

## Differences are not linked to severity phenotype or genotype. 
##Check for complete cases
sdt.stool%>% 
  group_by(Patient_number)%>%
  arrange(Visit, .by_group = TRUE)%>%
  summarise(n(), .groups = "keep")%>%
  dplyr::rename(n = "n()")%>%
  filter(n == 3)-> Keep

Keep<- Keep$Patient_number

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

###Add Antibiotic intake information 

antibiotic<- read.csv("~/CF_project/Metadata/sample_data_indexed_antibiotics_1.csv")

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

BC_dist.stool%>%
  left_join(tmp2, by="Patient_number")-> BC_dist.stool

##Add a time between visits (Overall for know but ask values per patient per period)
BC_dist.stool%>%
  dplyr::mutate(Months= case_when(Group == "V1_V3" ~ 12,
                                  Group == "V1_V2" ~ 3,
                                  Group == "V2_V3" ~ 19))-> BC_dist.stool

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

BC_dist.stool%>%
  left_join(tmp2, by="ID")-> BC_dist.stool

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

BC_dist.stool%>%
  left_join(tmp2, by="ID")-> BC_dist.stool

saveRDS(BC_dist.stool, "~/CF_project/exercise-cf-intervention/data/BC_dist.stool.rds")

##Is visit impacting differences in composition by patient? 
BC_dist.stool%>% 
  wilcox_test(BC_dist ~ Group)%>%
  adjust_pvalue(method = "BH") %>%
  add_significance()%>%
  add_xy_position(x = "Group")-> stats.test ## Not significant 

BC_dist.stool%>%
  wilcox_effsize(BC_dist ~ Group)

##Plot 
BC_dist.stool%>%
  ggplot(aes(x= Group, y= BC_dist))+
  geom_boxplot(color="black", alpha= 0.5)+
  geom_point(shape=21, position=position_jitter(0.2), size=3, aes(fill= Patient_number), color= "black")+
  xlab("Visit period")+
  ylab("Bray-Curtis dissimilarity")+
  labs(tag= "B)", caption = get_pwc_label(stats.test))+
  scale_fill_manual(values = pal.CF)+
  theme_classic()+
  theme(text = element_text(size=16), legend.position = "none")+
  scale_x_discrete(labels=c("V1_V2" =  "V1 to V2", 
                            "V2_V3" = "V2 to V3",
                            "V1_V3" = "V1 to V3"))->E

f<-ggarrange(D, E, ncol=1, nrow=2, common.legend = TRUE, legend="right") ##Final figure at script 9

#ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Stool.pdf", plot = f, width = 10, height = 8)
#ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Stool.png", plot = f, width = 10, height = 8)

rm(D,E,f)

##Is antibiotic burden differences impacting in composition by patient? 
### Linear model test
require("lmtest")
require("lme4")
test<- BC_dist.stool[!is.na(BC_dist.stool$AntibioticBurden_total),]

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

#Percentage of Variance explained
require("ggrepel")
pVarExpl <- ggplot (data = varianceTable) +
  geom_bar(aes (x = "", y = VarExplained, fill = Variable), width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c( "#DDCC77", "#CC6677","#117733","white"),
                    labels = c( "ABX iv burden", "ABX burden", "ABX prior Visit", "Residuals"))+
  theme_void() +
  geom_text_repel(aes(x=1.3, y = VarExplained/2, label=VarLabels))

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

##Correlation with Antibiotic intake 
BC_dist.stool%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Number_antibioticCourses_priorstudystart, y= BC_dist, color= Group))+
  geom_point(size=2.5, aes(shape= Group, fill= Patient_number), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  xlab("Number of antibiotic courses prior study start")+
  ylab("Bray-Curtis dissimilarity (Stool microbiome)")+
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

BC_dist.stool%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Number_antibioticCourses_duringstudy, y= BC_dist, color= Group))+
  geom_point(size=2.5, aes(shape= Group, fill= Patient_number), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  xlab("Number of antibiotic courses during study")+
  ylab("Bray-Curtis dissimilarity (Stool microbiome)")+
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

plot<-ggarrange(A, B, ncol=1, nrow=2, common.legend = TRUE, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Stool_Antibiotics.pdf", plot = plot, width = 10, height = 12)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Stool_Antibiotics.png", plot = plot, width = 10, height = 12)

rm(A,B, plot)

BC_dist.stool%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Number_iv_courses_priorstudy, y= BC_dist, color= Group))+
  geom_point(size=2.5, aes(shape= Group, fill= Patient_number), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  xlab("Number of iv antibiotic courses prior study start")+
  ylab("Bray-Curtis dissimilarity (Stool microbiome)")+
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

BC_dist.stool%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Number_iv_courses_duringstudy, y= BC_dist, color= Group))+
  geom_point(size=2.5, aes(shape= Group, fill= Patient_number), color= "black")+
  scale_shape_manual(values = c(21, 22, 24))+ 
  xlab("Number of iv antibiotic courses during study")+
  ylab("Bray-Curtis dissimilarity (Stool microbiome)")+
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

plot2<-ggarrange(C, D, ncol=1, nrow=2, common.legend = TRUE, legend="right")

ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Stool_iv_Antibiotics.pdf", plot = plot2, width = 10, height = 12)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Stool_iv_Antibiotics.png", plot = plot2, width = 10, height = 12)

###
BC_dist.stool%>%
  dplyr::filter(Group=="V1_V3")%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))%>%
  ggplot(aes(x= Number_antibioticCourses_duringstudy, y= BC_dist))+
  geom_point(size=2.5, shape= 21, aes(fill= Patient_number), color= "black")+
  xlab("Number of antibiotic courses during study")+
  ylab("Bray-Curtis dissimilarity V1 - V3 \n (Stool microbiome)")+
  labs(tag= "A)")+
  theme_classic()+
  geom_smooth(method=lm, se = T, color= "black")+
  stat_regline_equation(label.x = 2, label.y = 0.85)+ 
  stat_cor(label.x = 2, label.y = 0.8, 
           aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+
  scale_fill_manual(values = pal.CF)+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(fill = "Patient")+
  labs(shape = "Visit period")+
  theme(text = element_text(size=16))-> plot


ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Stool_V1V3_Antibiotics.pdf", plot = plot, width = 10, height = 12)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Q2_Beta_div_Stool_V1V3_Antibiotics.png", plot = plot, width = 10, height = 12)


###Matrix for Mantel test

BC_stool<- as.matrix(BC_dist)

##Change row and colnames for a uniform system between sample type Patient:Visit
rownames(BC_stool)<-  gsub("10P", "P", rownames(BC_stool))
rownames(BC_stool)<-  gsub("A", "\\1", rownames(BC_stool))

colnames(BC_stool)<-  gsub("10P", "P", colnames(BC_stool))
colnames(BC_stool)<-  gsub("A", "\\1", colnames(BC_stool))

saveRDS(BC_stool, "~/CF_project/exercise-cf-intervention/data/BC_stool_matrix.rds")

###Mixed effect models 
BC_dist.stool%>%
  dplyr::mutate(Patient_number = fct_relevel(Patient_number, 
                                             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                                             "P10", "P11", "P12", "P13", "P14","P15", "P16", "P17", "P18"))-> BC_dist.stool
##Check for complete cases
BC_dist.stool%>%
  dplyr::select(ppFVC, Trainingfrequency, Trainingtime, ppFEV1, BC_dist, Patient_number, Group, 
                Number_antibioticCourses_priorstudystart, Number_antibioticCourses_duringstudy, 
                Number_iv_courses_priorstudy, Number_iv_courses_duringstudy)%>%
  dplyr::filter(complete.cases(.))-> tmp

colnames(tmp)<- c("ppFVC", "Trainingfrequency", "Trainingtime", "ppFEV1", "BC_dist", "Patient_number", "Group", 
                  "ab_prior", "ab_during", "iv_prior", "iv_during")

##qqPlots (Check whther our variables are normaly distributed)
qqPlot(BC_dist.stool$BC_dist)
qqPlot(BC_dist.stool$Trainingfrequency)
qqPlot(BC_dist.stool$Trainingtime)
qqPlot(BC_dist.stool$ppFEV1)
qqPlot(BC_dist.stool$ppFVC)

##Based on the experimental design we should have a nested model with Intervisit as main grouping factor and patient 
##However, not all the grouping factors are complete so we should treat it as a crossed model 
##with patient and intervisit group as individual random effects

##Model selection do it with glm
full.model<- glm(BC_dist ~ Trainingfrequency*Trainingtime*ppFEV1*ppFVC, data = tmp) ##Full model
# Stepwise regression model
step.model <- MASS::stepAIC(full.model, direction = "both", 
                            trace = FALSE)

summary(step.model)
step.model.df<- as.data.frame(coef(summary(step.model))) ##The full model is the best!!!

step.model.df%>%
  mutate(p.adj = p.adjust(`Pr(>|t|)`, method='BH')) %>%
  add_significance()%>%
  rownames_to_column()-> step.model.df

write.csv(step.model.df, "~/CF_project/exercise-cf-intervention/tables/Q2_BC_Sports_Lung_Stool.csv", row.names = F)

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

write.csv(step.model.df, "~/CF_project/exercise-cf-intervention/tables/Q2_BC_Stool_antibiotics.csv", row.names = F)

## From model selection step ppFVC:ppFEV1 interaction are the best explanatory variables for BC dissimilarities
tr0<- lmer(BC_dist ~ (1 | Patient_number), data = tmp) ##Null model

tr2<-lmer(BC_dist ~ Trainingfrequency*Trainingtime+ (1 | Patient_number), data = tmp)
summary(tr2) ##Training

tr3<-lmer(BC_dist ~ ab_prior + ab_during + iv_prior + iv_during + 
            ab_prior:ab_during + (1 | Patient_number), data = tmp)
summary(tr3) ##Antibiotics

lrtest(tr0, tr3)

##Each factor ad predictive value to the model 
##What happen with interactions 
summary(glm(BC_dist ~ Trainingfrequency*Trainingtime*ppFEV1*ppFVC, data = tmp))
tr4<-lmer(BC_dist ~ Trainingfrequency*Trainingtime*ppFEV1*ppFVC + (1 | Patient_number), data = tmp)
summary(tr4) ##"Best" model

step.model.df<- as.data.frame(coef(summary(tr4))) ##T

step.model.df%>%
  mutate(p.adj = p.adjust(`Pr(>|t|)`, method='BH')) %>%
  add_significance()%>%
  rownames_to_column()-> step.model.df

write.csv(step.model.df, "~/CF_project/exercise-cf-intervention/tables/Q2_GLMM_BC_Sports_Lung_Stool.csv", row.names = F)

lrtest(tr4, tr2)
lrtest(tr4, tr3)

A<- plotREsim(REsim(tr4))  ## plot the interval estimates
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
B<- plot_model(tr4, p.adjust = "BH", vline.color = "gray",
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
##Model selection do it with glm
full.model<- glm(ppFEV1 ~ Trainingfrequency*Trainingtime*BC_dist, data = tmp) ##Full model
# Stepwise regression model
step.model <- MASS::stepAIC(full.model, direction = "both", 
                            trace = FALSE)
summary(step.model)

full.model<- glm(ppFVC ~ Trainingfrequency*Trainingtime*BC_dist, data = tmp) ##Full model
# Stepwise regression model
step.model <- MASS::stepAIC(full.model, direction = "both", 
                            trace = FALSE)
summary(step.model)

best.model<- glm(ppFVC ~ Trainingfrequency+ BC_dist + Trainingfrequency*BC_dist, data = tmp)
summary(best.model)

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
