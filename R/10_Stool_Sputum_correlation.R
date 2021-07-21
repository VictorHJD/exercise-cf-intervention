##Cystic-Fibrosis Microbiome Project (Mainz)
##Data analysis: Mantel test Stool and Sputum samples
##Víctor Hugo Jarquín-Díaz 28.04.2021

library(vegan)

##Load data

BC_sputum <- readRDS("../vjarqui/CF_project/exercise-cf-intervention/data/BC_sputum_matrix.rds")
BC_stool<- readRDS("../vjarqui/CF_project/exercise-cf-intervention/data/BC_stool_matrix.rds")

##Check that both have same samples in same order 
##1) Intersection
sput_stool_inter <- intersect(colnames(BC_sputum), colnames(BC_stool))

##2) Subset
MSm<- BC_sputum[sput_stool_inter, sput_stool_inter]
MSl<- BC_stool[sput_stool_inter, sput_stool_inter]

##3) Mantel test 
sput_stool <-  mantel(MSm, MSl, method = "spearman", permutations = 9999, na.rm = TRUE)
sput_stool

#Mantel statistic based on Spearman's rank correlation rho 

#Call:
#mantel(xdis = MSm, ydis = MSl, method = "spearman", permutations = 9999,      na.rm = TRUE) 

#Mantel statistic r: 0.2759 
#      Significance: 3e-04 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#0.110 0.139 0.165 0.192 
#Permutation: free
#Number of permutations: 9999

###Network
library(phyloseq)
library(igraph)
library(ggnetwork)

##Stool 
PS4.stool.2<- PS4.stool
PS4.stool.2@sam_data<- sample_data(tmp.stool)

##Use bray-curtis dissimilarity distance
net.stool <- make_network(PS4.stool.2, type="samples", distance="bray",
                          max.dist = 0.8, keep.isolates=T)
  
set.seed(2020)
plot_network(net.stool, PS4.stool.2, line_weight=0.3, color = NULL, shape = "Visit",
             label= "Patient_number", line_alpha = 0.1,  hjust = 2.5)+ 
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(fill= Genus, shape= Visit), color= "black")+
  scale_shape_manual(values = c(21, 22, 24), labels = c("Visit 1", "Visit 2", "Visit 3"))+
  scale_fill_manual(values = tax.palette)+
  labs(tag= "A)")+
  theme(text = element_text(size=16))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= F)+
  labs(fill = "Dominant taxa")+
  labs(shape = "Visit")-> A

##Reduce maximum distance to distinguish relevant clusters
net.stool.gen <- make_network(PS4.stool.Gen, type="taxa", distance="bray",
                             max.dist = 0.5, keep.isolates=F)

set.seed(2020)
plot_network(net.stool.gen, PS4.stool.Gen, line_weight=0.5, type = "taxa", color = NULL,
             line_alpha = 0.1,  hjust = 1.5)+ 
  theme(aspect.ratio=0.5)+
  geom_point(shape= 21, size=3, aes(fill= Phylum), color= "black")+
  labs(tag= "B)")+
  theme(text = element_text(size=16), legend.position="none")-> C

##Testing relation between graph and factors:
#Compute the null distribution and p-value for the test

library("phyloseqGraphTest")
gt.stool<- graph_perm_test(PS4.stool.2, "Genus", distance="bray",
                           type = c("mst", "knn", "threshold.value",
                                    "threshold.nedges"),  nperm=1000)
gt.stool$pval
plot_permutations(gt.stool)

##Sputum
PS4.sputum.2<- PS4.sput
PS4.sputum.2@sam_data<- sample_data(tmp.sputum)

##Use bray-curtis dissimilarity distance
net.sput <- make_network(PS4.sputum.2, type="samples", distance="bray",
                          max.dist = 0.8, keep.isolates=T)

set.seed(2020)
plot_network(net.sput, PS4.sputum.2, line_weight=0.3, color = NULL, shape = "Visit",
             label= "Patient_number", line_alpha = 0.1,  hjust = 2.5)+ 
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(fill= Genus, shape= Visit), color= "black")+
  scale_shape_manual(values = c(21, 22, 24), labels = c("Visit 1", "Visit 2", "Visit 3"))+
  scale_fill_manual(values = tax.palette)+
  labs(tag= "B)")+
  theme(text = element_text(size=16))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= F)+
  labs(fill = "Dominant taxa")+
  labs(shape = "Visit")-> B

##Taxa level
net.sput.gen <- make_network(PS4.sput.Gen, type="taxa", distance="bray",
                         max.dist = 0.7, keep.isolates=F)

set.seed(2020)
plot_network(net.sput.gen, PS4.sput.Gen, line_weight=0.5, type = "taxa", color = NULL,
              line_alpha = 0.1,  hjust = 1.5)+ 
  theme(aspect.ratio=0.5)+
  geom_point(shape= 21, size=3, aes(fill= Phylum), color= "black")+
  labs(tag= "A)")+
  theme(text = element_text(size=16), legend.position="bottom")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= F)+
  labs(fill = "Phylum")-> D

##Tests
gt.sput<- graph_perm_test(PS4.sputum.2, "Genus", distance="bray",
                           type = c("mst", "knn", "threshold.value",
                                    "threshold.nedges"),  nperm=1000)
gt.sput$pval
plot_permutations(gt.sput)


##Save
plot1<- ggarrange(A, B, ncol=2, nrow=1)

ggsave(file = "CF_project/exercise-cf-intervention/figures/Fig5_Patient_Network.pdf", plot = plot1, width = 14, height = 5, dpi = 400)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Fig5_Patient_Network.png", plot = plot1, width = 14, height = 5, dpi = 400)

plot2<- ggarrange(D, C, ncol=2, nrow=1, common.legend = T)

ggsave(file = "CF_project/exercise-cf-intervention/figures/Fig6_ASV_Network.pdf", plot = plot2, width = 14, height = 6, dpi = 400)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Fig6_ASV_Network.png", plot = plot2, width = 14, height = 6, dpi = 400)

rm(A,B, C, D, plot1, plot2)
###Line plot to compare change of total proportion of satellites vs colonizer during study
##Get the sum of each group by patient by Visit

Type.sputum%>%
  dplyr::group_by(Patient_number, Visit, Type)%>%
  dplyr::summarise(sum= sum(Abundance))%>%
  dplyr::group_by(Visit, Type)%>%
  dplyr::summarise(mean = mean(sum), sd= sd(sum), n = n())%>%
  dplyr::mutate(se = sd / sqrt(n),
                lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
                upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se)%>%
  dplyr::mutate(upper.ran = mean + (2*sd), 
                lower.ran = mean - (2*sd))%>%
  ggplot(aes(x=Visit, y=mean, group = Type, color=Type))+ 
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), width=.1) +
  geom_line(aes(linetype=Type)) + 
  geom_point(shape=21, aes(fill= Type))+
  labs(x="Visit", y = "Relative abundance (%)", tag = "A)", color= "Group of taxa")+
  guides(color = guide_legend(override.aes=list(shape=c(21), fill= c('#999999','#E69F00')), 
                              nrow = 1, size= 10),linetype= "none", fill= "none")+
  theme_classic() + 
  scale_color_manual(values=c('#999999','#E69F00'))+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  theme(text = element_text(size=16))-> A

Type.stool%>%
  dplyr::group_by(Patient_number, Visit, Type)%>%
  dplyr::summarise(sum= sum(Abundance))%>%
  dplyr::group_by(Visit, Type)%>%
  dplyr::summarise(mean = mean(sum), sd= sd(sum), n = n())%>%
  dplyr::mutate(se = sd / sqrt(n),
                lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
                upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se)%>%
  dplyr::mutate(upper.ran = mean + (2*sd), 
                lower.ran = mean - (2*sd))%>%
  ggplot(aes(x=Visit, y=mean, group = Type, color=Type))+ 
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), width=.1) +
  geom_line(aes(linetype=Type)) + 
  geom_point(shape=21, aes(fill= Type))+
  labs(x="Visit", y = "Relative abundance (%)", tag = "B)", color= "Group of taxa")+
  guides(color = guide_legend(override.aes=list(shape=c(21), fill= c('#999999','#E69F00')), 
                              nrow = 1, size= 10),linetype= "none", fill= "none")+
  theme_classic() + 
  scale_color_manual(values=c('#999999','#E69F00'))+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  theme(text = element_text(size=16))-> B

##Save
plot1<- ggarrange(A, B, ncol=1, nrow=2, legend = "top", common.legend = T)

ggsave(file = "CF_project/exercise-cf-intervention/figures/Fig6_Colonizer_Satellites.pdf", plot = plot1, width = 10, height = 8, dpi = 400)
ggsave(file = "CF_project/exercise-cf-intervention/figures/Fig6_Colonizer_Satellites.png", plot = plot1, width = 10, height = 8, dpi = 400)
