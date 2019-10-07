library(dada2)
library(ape)
packageVersion("dada2")
library(dplyr)
library(ggplot2)
library(ggpubr)
library(DECIPHER)
library(phyloseq)
library(phangorn)
library(openxlsx)
library(picante)
#source("https://bioconductor.org/biocLite.R")
#workflow reference https://compbiocore.github.io/metagenomics-workshop/assets/DADA2_tutorial.html
#biocLite("GenomeInfoDbData")
#biocLite("DECIPHER")
#biocLite("phyloseq")
#biocLite("phangorn")

###Save Working
save.image("/Users/arianaliu1214/Desktop/Microb_counting/DNA/Diversity coding /DNA sequence for 14 cruises/18s all/14c_HNF/18s_dada2_phyloseq.rdata")
# Close R, Re-open R
load("/Users/arianaliu1214/Desktop/Microb_counting/DNA/Diversity coding /DNA sequence for 14 cruises/18s all/14c_HNF/18s_dada2")
###Save individual
saveRDS(seqtab.nochim, "/Users/arianaliu1214/Desktop/Microb_counting/DNA/Diversity coding /DNA sequence for 14 cruises/18s all/14c_HNF/seqtab.nochim.18s.rds")
# Close and re-open R
seqtab <- readRDS("/Users/ariana/Desktop/Microb_counting/DNA/Diversity coding /DNA sequence for 14 cruises/18s all/14c_HNF/seqtab.nochim.18s.rds")




path= "/Users/arianaliu1214/Desktop/Microb_counting/DNA/Diversity coding /DNA sequence for 14 cruises/18s all/14c_HNF/"
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(c(list.files(path, pattern="_R1_001.fastq", full.names = TRUE),list.files(path, pattern="_R1_001.clean.fastq", full.names = TRUE)))
fnRs <- sort(c(list.files(path, pattern="_R2_001.fastq", full.names = TRUE),list.files(path, pattern="_R2_001.clean.fastq", full.names = TRUE)))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
# sample.names <- paste(sapply(strsplit(basename(fnFs), "_"), `[`, 1),"_",sapply(strsplit(basename(fnFs), "_"), `[`, 4),sep="")
# sample.names <-substr(basename(fnFs),7,19)
sample.names <- paste(sapply(strsplit(basename(fnFs), "_"), `[`, 1)) 
#plot different runs for strick quality control
plotQualityProfile(fnFs[c(31,39,8,82,96)])+geom_hline(yintercept=30, linetype="dashed", color = "blue")
plotQualityProfile(fnRs[c(3,39,8,82,96)])+geom_hline(yintercept=30, linetype="dashed", color = "blue")

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#Cut Q30 and +- 10 for high quailty
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(270,170), trimLeft=c(20,18),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]


mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
### give a range of sequences

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(371,390)]
table(nchar(getSequences(seqtab2)))
plot(table(nchar(getSequences(seqtab2))))
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))

sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Save OTU table
#write.xlsx(seqtab.nochim,file = "18s_seqtab_nochim(seq_no).xlsx")

########### END Trim DNA ##########################

########### Start Taxa ##########################
### change to silva 132
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa", multithread=TRUE)
taxa <- addSpecies(taxa, "silva_species_assignment_v132.fa")

#system.time({taxa<-addSpecies(taxa, "silva_species_assignment_v132.fa",verbose=T)
#colnames(taxa)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
#unname(taxa)})

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
### save taxa table
#write.table(taxa, "taxa_seq_18s.txt",quote=FALSE,sep="\t")
#write.table(taxa.print, "taxa_print_18s.txt",quote=FALSE,sep="\t")

#####################################
#####Construct Phylogenetic Tree#####
#####################################

#Extract sequences from DADA2 output
seqtab.nochim<-read.xlsx("18s_seqtab_nochim(seq_no).xlsx",sheet=1,colNames = T,rowNames = T)
sequences<-getSequences(t(seqtab.nochim))
names(sequences)<-sequences
#Run Sequence Alignment (MSA) using DECIPHER
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)
#Change sequence alignment output into a phyDat structure
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
#Create distance matrix
dm <- dist.ml(phang.align)
#Preform Neighbor joining
treeNJ <- NJ(dm) # Note, tip order != sequence order
#Internal maximum likelihood
fit = pml(treeNJ, data=phang.align)
#negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
#Save tree
#tree18s = phy_tree(treeNJ)
#ape::write.tree(tree18s, "/Users/arianaliu1214/Desktop/Microb_counting/DNA/Diversity coding /DNA sequence for 14 cruises/18s all/14c_HNF/treeNJ_18s.tree")
#=========== End ===============================#



#####################################
#####    Import to phyloseq     #####
#####################################
#Reference https://vaulot.github.io/tutorials/Phyloseq_tutorial.html
###For OTU and tax, use "as.matrix" to avoid error in phyloseq step
###(Pre)Processing Data
# otu<-read.xlsx("HNF_OTUtaxa(20190827).xlsx",sheet=1,colNames = T,rowNames = T)
# tax<-read.xlsx("HNF_OTUtaxa(20190827).xlsx",sheet=2,colNames = T,rowNames = T)
# treeNJ<-read.tree("treeNJ_18s.tree")

otu<- as.matrix(t(seqtab.nochim), taxa_are_rows=T,colNames = TRUE) #ok
tax<- as.matrix(taxa)
identical(rownames(tax), rownames(otu))
sample.names <- read.xlsx("14c_group_info.xlsx",colNames = TRUE)
rownames(sample.names)<-sample.names[,2]
identical(row.names(sample.names),colnames(otu))
#if FALSE: write.xlsx(colnames(otu),"column_names_otu.xlsx")
#rownames(otu)<-c(paste0("seq",1:6882)),Change sequence info to seq.numbers

###Mapping file 1st column name must
#samdf<-seqtab.nochim[,1]
#map <- import_qiime_sample_data("~/DADA2_Tutorial/Tutorial_Map.txt")
#identical(rownames(seqtab.nochim), rownames(sample.names))---if TRUE continue
map <- read.xlsx("14c_group_info.xlsx",colNames = T,rowNames = T)# environment factors or grouping factors of samples
# mapfile <- system.file("extdata", "master_map.txt", package = "phyloseq")
#map <- import_qiime_sample_data(read.csv2("14c_group_info.csv",header=F,row.names=F,col.names = T,sep=","))
map <- map[!duplicated(map$SampleID),] 
rownames(map)<-colnames(otu)
identical(row.names(map),colnames(otu)) #double check to be TRUE
ps <- phyloseq(otu_table(otu, taxa_are_rows=T),phy_tree(treeNJ),tax_table(tax),sample_names(sample.names),sample_data(map))
#Merge PhyloSeq object with map, merge_phyloseq(ps,map)
ps

#Currently, the phylogenetic tree is not rooted Though it is not necessary here, 
#you will need to root the tree if you want to calculate any phylogeny based diversity metrics (like Unifrac)
set.seed(711)
phy_tree(ps) <- root(phy_tree(ps), sample(taxa_names(ps), 1), resolve.root = TRUE)
is.rooted(phy_tree(ps))


###Trim taxa###
HNF_species_seq<-read.xlsx("HNF_taxaNamelist.xlsx",sheet = 3,colNames = TRUE)
rownames(HNF_species_seq)<-HNF_species_seq[,1]
HNF_taxa=read.xlsx("HNF_taxaNamelist.xlsx",sheet=2)
rank_names(ps)

ps_trim<- subset_taxa(ps,rownames(tax_table(ps))%in%rownames(HNF_species_seq))
'%ni%' <- Negate('%in%')
ps_trim<-subset_samples(ps_trim,rownames(sam_data(ps_trim))%ni%"S18F12Or2Cr2023St01SC40A") #Due to low sequence count
ps_trim
# ex1 = prune_taxa(myTaxa, GlobalPatterns)

#Rarefy by phyloseq
#Reference: from dada2 benjjneb
#http://joey711.github.io/phyloseq-demo/phyloseq-demo.html
#Plot original sample size
rarecurve(t(otu_table(ps_trim)), step=50, cex=0.5)
# rarefy without replacement
ps.rarefied = rarefy_even_depth(ps_trim, rngseed=1, sample.size=0.9*min(sample_sums(ps_trim)), replace=F)

rarecurve(t(otu_table(ps.rarefied)), step=50, cex=0.5)
sample_sums(ps.rarefied)#Check sum of each sample site is same

#Add pd from picante package to phyloseq estimate function
estimate_richness2<-function (physeq, split = TRUE, measures = NULL) 
{
  if (!any(otu_table(physeq) == 1)) {
    warning("The data you have provided does not have\n", 
            "any singletons. This is highly suspicious. Results of richness\n", 
            "estimates (for example) are probably unreliable, or wrong, if you have already\n", 
            "trimmed low-abundance taxa from the data.\n", "\n", 
            "We recommended that you find the un-trimmed data and retry.")
  }
  if (!split) {
    OTU <- taxa_sums(physeq)
  }
  else if (split) {
    OTU <- as(otu_table(physeq), "matrix")
    if (taxa_are_rows(physeq)) {
      OTU <- t(OTU)
    }
  }
  renamevec = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", 
                "InvSimpson", "Fisher")
  names(renamevec) <- c("S.obs", "S.chao1", "S.ACE", "shannon", 
                        "simpson", "invsimpson", "fisher")
  if (is.null(measures)) {
    measures = as.character(renamevec)
  }
  if (any(measures %in% names(renamevec))) {
    measures[measures %in% names(renamevec)] <- renamevec[names(renamevec) %in% 
                                                            measures]
  }
  if (!any(measures %in% renamevec)) {
    stop("None of the `measures` you provided are supported. Try default `NULL` instead.")
  }
  outlist = vector("list")
  estimRmeas = c("Chao1", "Observed", "ACE")
  if (any(estimRmeas %in% measures)) {
    outlist <- c(outlist, list(t(data.frame(estimateR(OTU)))))
  }
  if ("Shannon" %in% measures) {
    outlist <- c(outlist, list(shannon = diversity(OTU, index = "shannon")))
  }
  if ("Simpson" %in% measures) {
    outlist <- c(outlist, list(simpson = diversity(OTU, index = "simpson")))
  }
  if ( "FaithPD" %in% measures){
    outlist <- c(outlist, list(FaithPD = t(picante::pd(samp = OTU, tree = phy_tree(physeq), include.root = F))[1,] ))
  }
  if ("InvSimpson" %in% measures) {
    outlist <- c(outlist, list(invsimpson = diversity(OTU, 
                                                      index = "invsimpson")))
  }
  if ("Fisher" %in% measures) {
    fisher = tryCatch(fisher.alpha(OTU, se = TRUE), warning = function(w) {
      warning("phyloseq::estimate_richness: Warning in fisher.alpha(). See `?fisher.fit` or ?`fisher.alpha`. Treat fisher results with caution")
      suppressWarnings(fisher.alpha(OTU, se = TRUE)[, c("alpha", 
                                                        "se")])
    })
    if (!is.null(dim(fisher))) {
      colnames(fisher)[1:2] <- c("Fisher", "se.fisher")
      outlist <- c(outlist, list(fisher))
    }
    else {
      outlist <- c(outlist, Fisher = list(fisher))
    }
  }
  out = do.call("cbind", outlist)
  namechange = intersect(colnames(out), names(renamevec))
  colnames(out)[colnames(out) %in% namechange] <- renamevec[namechange]
  colkeep = sapply(paste0("(se\\.){0,}", measures), grep, colnames(out), 
                   ignore.case = TRUE)
  out = out[, sort(unique(unlist(colkeep))), drop = FALSE]
  out <- as.data.frame(out)
  return(out)
}


###Alpha Diversity###
div_a<-estimate_richness2(ps.rarefied,split = TRUE,c("Observed", "Chao1","Shannon", "Simpson","FaithPD"))
#Save diversity index
write.xlsx(div_a,file = "HNF_a_diversity_rar(20190830).xlsx",colNames=T,rowNames=T)


#Beta Diversity
ord.nmds.bray<- ordinate(ps_trim,method="NMDS",distance = "bray")

p1<-plot_ordination(ps_trim,ord.nmds.bray,color="Cruise",title="Cruise",label = "Cruise")
p2<-plot_ordination(ps1,ord.nmds.bray,color="Station",title="Station",label = "Cruise")
p3<-plot_ordination(ps1,ord.nmds.bray,color="Q4temp",title="Temperature Quartile",label = "Cruise")+
  scale_color_manual(values = c("#99CCFF","#6699FF","#0033CC","#003366"))
nmds_plot<-ggarrange(p1,p2,p3, labels = c("A", "B","C"),legend = "right",ncol=3,nrow=1)
annotate_figure(nmds_plot,top = text_grob("Bray NMDS-18s", face = "bold", size = 14))


# CAP ordinate
# http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
# Remove data points with missing metadata
erie_not_na <- ps1 %>%
  subset_samples(
    !is.na(q0_obs_Bac)&!is.na(q1_obs_Bac)&!is.na(q2_obs_Bac)&
      !is.na(prey_biomass)&!is.na(ABB)&!is.na(HBB)&
      !is.na(temperature)&!is.na(salinity)&
      !is.na(nitrite)&!is.na(nitrate)&!is.na(phosphate))

bray_not_na <- phyloseq::distance(physeq = erie_not_na, method = "bray")

cap_ord <- ordinate(
  physeq = erie_not_na, 
  method = "CAP",
  distance = bray_not_na,
  formula = ~ q1_obs_Bac + q2_obs_Bac + ABB + HBB + temperature +salinity+nitrite+nitrate+phosphate
)

# CAP plot
cap_plot <- plot_ordination(
  physeq = erie_not_na, 
  ordination = cap_ord, 
  color = "Cruise", 
  axes = c(1,2)) + 
  aes(shape = Station) + 
  geom_point(aes(colour = Cruise), size = 4) + 
  geom_point(colour = "grey90", size = 1.5) 


# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
anova(cap_ord)
cap_plot<-cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  )
annotate_figure(cap_plot,top = text_grob("CAP-HNF", face = "bold", size = 14),bottom  = text_grob(paste("anova model p=",anova(cap_ord)[1,4]), face = "bold", size = 8))
