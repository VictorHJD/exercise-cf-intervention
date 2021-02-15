##Cystic-Fibrosis Microbiome Project (Mainz)
##dada2 pipeline v1.18
##Víctor Hugo Jarquín-Díaz 08.02.2021

library("dada2")

Refil<- FALSE ##Turn to TRUE just if filter parameters are changed 
Taxanot<- FALSE ##Turn to TRUE just in case new taxonomic assignment is required


if(Refil){
  ###Load files 
  path <- "~/CF_project/Seqdata/Commed_Subset/" # CHANGE ME to the directory containing the fastq files after unzipping.
  list.files(path)
  fastqFiles <- list.files(path, pattern=".fastq$", full.names=TRUE) #take all fastaq files from the folder 
  fastqF <- grep("_R1_001.fastq", fastqFiles, value = TRUE) #separate the forward reads
  fastqR <- grep("_R2_001.fastq", fastqFiles, value = TRUE) #separate the reverse reads 
  
  samples <- gsub("_\\d+_L001_R1_001.fastq", "\\1", basename(fastqF))
  
  ###Check quality of the reads 
  plotQualityProfile(fastqF[10:11])
  plotQualityProfile(fastqR[10:11])
  
  filt_path <- "~/CF_project/tmp/filtered" ##Filtering
  
  #Pipeline filtration 
  if(!file_test("-d", filt_path)) dir.create(filt_path)
  filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq"))
  names(filtFs) <- samples
  filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq"))
  names(filtRs) <- samples
  
  out <- filterAndTrim(fastqF, filtFs, fastqR, filtRs, truncLen=c(245, 240), ##This gives just 10bp overlap, Prev. conditions were 240,250
                       maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft = c(19, 21), ##Remove primers
                       #'16s_V4': ('GTGCCAGCMGCCGCGGTAA','GACTACHVGGGTWTCTAATCC'),#forward:515f, reverse:806r
                       compress=TRUE, multithread=TRUE)
  head(out)

  ###Learning errors
  errF <- learnErrors(filtFs, multithread=TRUE)
  #100351684 total bases in 444034 reads from 15 samples will be used for learning the error rates.
  errR <- learnErrors(filtRs, multithread=TRUE)
  #102335634 total bases in 467286 reads from 16 samples will be used for learning the error rates.
  
  plotErrors(errF, nominalQ=TRUE)
  plotErrors(errR, nominalQ=TRUE)
  
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  derepRs <- derepFastq(filtRs, verbose=TRUE)
  # Name the derep-class objects by the sample names
  names(derepFs) <- samples
  names(derepRs) <- samples
  
  ##Sample inference 
  dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
  dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
  
  ##Merge paired reads
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
  # Inspect the merger data.frame from the first sample
  head(mergers[[1]])
  
  ##Construction of sequence table 
  seqtab <- makeSequenceTable(mergers)
  dim(seqtab)
  
  # Inspect distribution of sequence lengths
  table(nchar(getSequences(seqtab))) ##Everything looks good but we have some short fragments 
  plot(table(nchar(getSequences(seqtab))))
  
  ##Since our expected amplicon size is 440bp with primers, Let's make an in silico cut
  seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 249:253] ##New filtering gives amplicons ~251bp expected as amplicon
  plot(table(nchar(getSequences(seqtab2))))
  
  ##Remove chimeras 
  ##Use just the data with the expected size
  seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
  dim(seqtab.nochim)
  sum(seqtab.nochim)/sum(seqtab)
  ##With new filtering conditions 53.5% of the reads passed 
  saveRDS(seqtab.nochim, "~/CF_project/tmp/seqtab_final.rds") ##Final sequence table without chimeras
}

if(Taxanot){
  seqtab.nochim<- readRDS("~/CF_project/tmp/seqtab_final.rds")
  
  ##Create an ASV matrix with ASV as rows and samples as columns
  asvmat <- t(seqtab.nochim) #Removing sequence rownames for display only
  asv.names <- as.data.frame(rownames(asvmat))
  rownames(asv.names) <- paste0("ASV", 1:nrow(asv.names))
  rownames(asvmat)<-NULL
  rownames(asvmat) <- paste0("ASV", 1:nrow(asvmat))
  #colnames(asvmat) <- paste0("Sample", 1:ncol(asvmat)) #keep original sample ID
  head(asvmat)
  write.csv(asvmat, "~/CF_project/output/ASV_matrix.csv")
  
  ##Get count of ASVs detected by sample
  asv.sample<- as.data.frame(asvmat)
  test<- data.frame()
  for (i in 1:ncol(asv.sample)) {
    asv<- data.frame()
    asv[1,1]<- sum(asv.sample[,i]!=0)
    rownames(asv)<- paste0("Sample", i)
    test <- rbind(test, asv) ### Join all the "individual" data frames into the final data frame 
  }
  asv.sample<- as.matrix(test)
  colnames(asv.sample)<- "ASVs_dada2"
  rm(test,asv, i)
  
  ###Track reads through the pipeline 
  getN <- function(x) sum(getUniques(x))
  track <- cbind(sapply(fastqF, getN), sapply(fastqR, getN), sapply(filtFs, getN), sapply(filtRs, getN), sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
  #If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
  track<-track[,c(1,3,5:8)]
  colnames(track) <- c("input_dada2", "filtered_dada2", "denoisedF_dada2", "denoisedR_dada2", "merged_dada2", "nonchim_dada2")
  rownames(track) <- samples
  rownames(asv.sample) <- samples
  track<-cbind(track, asv.sample)
  head(track)
  saveRDS(track, "~/CF_project/output/Track_DADA2.rds") ##Final track of dada2 pipeline
  
  ##Taxonomic annotation using naive Bayesian classifier from dada2 with SILVA db version 138
  taxa <- assignTaxonomy(seqtab.nochim, "~/SILVA_db/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
  taxa <- addSpecies(taxa, "~/SILVA_db/silva_species_assignment_v138.fa.gz")
  saveRDS(taxa, "~/CF_project/output/tax_final.rds") ##Final taxonomy
  
  ##Create an Taxa matrix with ASV as rows and taxonomic level as columns
  taxamat <- taxa # Removing sequence rownames for display only
  rownames(taxamat) <- NULL
  rownames(taxamat) <- paste0("ASV", 1:nrow(taxamat))
  head(taxamat)
  write.csv(taxamat, "~/CF_project/output/Taxa_matrix.csv")
  
  ##Transform seqtab.nochim to fasta file 
  library(DECIPHER); packageVersion("DECIPHER")
  ##Create a DNAString set from the ASVs
  dna <- DNAStringSet(getSequences(seqtab.nochim))
  names(dna)<- paste0("ASV", 1:length(dna)) ##Give short names to each sequence
  writeXStringSet(dna, "~/CF_project/output/ASV.fasta") ##Export fasta seq
  
  ##Create a Phylogenetic tree
  library(phangorn)
  Align16S<- AlignSeqs(dna, anchor= NA, verbose= FALSE) ##Alignment
  
  phangAlign16S <- phyDat(as(Align16S, "matrix"), type="DNA")
  dm16S <- dist.ml(phangAlign16S) ## Distance matrix
  treeNJ16S <- NJ(dm16S) # Note, tip order != sequence order
  plot(treeNJ16S, type= "unrooted", use.edge.length= FALSE, no.margin= TRUE, show.tip.label= TRUE) ##Neighbor-Joining tree
  fit1 <- pml(treeNJ16S, data=phangAlign16S)
  fitGTR16S <- update(fit1, k=4, inv=0.2)
  fitGTR16S <- optim.pml(fitGTR16S, model="GTR", optInv=TRUE, optGamma=TRUE,
                         rearrangement = "stochastic", control = pml.control(trace = 0))
  plot(fitGTR16S, type= "unrooted", use.edge.length= FALSE, no.margin= TRUE, show.tip.label= FALSE)
  
  tree<- fitGTR16S$tree
  
  saveRDS(fitGTR16S, "~/CF_project/output/phylo_tree.rds")
  
}
  
  