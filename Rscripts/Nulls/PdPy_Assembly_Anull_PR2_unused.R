###############################################################################################
##### Loading packages ########################################################################
###############################################################################################
if (!require(tidyverse)) {
  install.packages("tidyverse", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(tidyverse)
}else{library(tidyverse)}

if (!require(vegan)) {
  install.packages("vegan", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(vegan)
}else{library(vegan)}

if (!require(ape)) {
  install.packages("ape", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(ape)
}else{library(ape)}

if (!require(parallel)) {
  install.packages("vegan", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(parallel)
}else{library(parallel)}

if (!require(SpadeR)) {
  install.packages("SpadeR", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(SpadeR)
}else{library(SpadeR)}

if (!require(iNEXT)) {
  install.packages("iNEXT", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(iNEXT)
}else{library(iNEXT)}

if (!require(picante)) {
  install.packages("picante", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(picante)
}else{library(picante)}

if (!require(geiger)) {
  install.packages("geiger", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(geiger)
}else{library(geiger)}

if (!require(GUniFrac)) {
  install.packages("GUniFrac", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(GUniFrac)
}else{library(GUniFrac)}

if (!require(phytools)) {
  install.packages("phytools", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(phytools)
}else{library(phytools)}

if (!require(ecodist)) {
  install.packages("ecodist", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(ecodist)
}else{library(ecodist)}

if (!require(ade4)) {
  install.packages("ade4", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(ade4)
}else{library(ade4)}

if (!require(abind)) {
  install.packages("abind", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(abind)
}else{library(abind)}
###############################################################################################
##### Loading packages ########################################################################
###############################################################################################

###############################################################################################
##### Loading data ############################################################################
###############################################################################################
Bac_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_Bac_seqXst_PR2_3.csv",
                                       sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
Bac_ra_comm <- Bac_comm / rowSums(Bac_comm)
Bac_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_Bac_treeNJ_PR2_3.tree")

HNF_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_HNF_seqXst_PR2_3.csv",
                                       sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
HNF_ra_comm <- HNF_comm / rowSums(HNF_comm)
HNF_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_HNF_treeNJ_PR2_3.tree")

NF_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_prot_seqXst_PR2_3.csv",
                                      sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
NF_ra_comm <- NF_comm / rowSums(NF_comm)
NF_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_prot_treeNJ_PR2_3.tree")
###############################################################################################
##### Loading data ############################################################################
###############################################################################################

###############################################################################################
##### Phylogenetic turnover / deterministic vs stochastic processes_Ampd #######################
###############################################################################################

##### Nanoflagellate phylogenetic turnover ##########################################
ini <- Sys.time()
numCores <- detectCores()
numCores

cl <- makeCluster(numCores - 24)

clusterEvalQ(cl, {
  library(vegan)
  library(tidyverse)
  library(picante)
  NF_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_prot_seqXst_PR2_3.csv", 
                                        sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
  NF_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_prot_treeNJ_PR2_3.tree")
})

NFPhylo_null_func <- function(x){
  # set up parameters and data
  community = NF_comm
  phylo = NF_phylo
  # performing comdistnt to calculate MNTD
  as.matrix(mpd(community, cophenetic(tipShuffle(phylo)), abundance.weighted = TRUE))
}
nsim.list <- sapply(1:999, list)
test <- parLapply(cl, nsim.list, NFPhylo_null_func)
test[[1000]] <- as.matrix(mpd(NF_comm, cophenetic(NF_phylo), abundance.weighted = TRUE))

Sys.time() - ini

NFPhylo_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(row.names(NF_comm)) %>%
  rename(obs = X1000)
write.table(NFPhylo_null, file = "/home/zac422/Desktop/OSCAR/Nulls/NF_Ampd_null_3.csv",
            sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
##### Nanoflagellate phylogenetic turnover ##########################################

###############################################################################################
##### Phylogenetic turnover / deterministic vs stochastic processes_Ampd ######################
###############################################################################################

###############################################################################################
##### Phylogenetic turnover / deterministic vs stochastic processes_Amntd #####################
###############################################################################################

##### Nanoflagellate phylogenetic turnover ##########################################
ini <- Sys.time()
numCores <- detectCores()
numCores

cl <- makeCluster(numCores - 24)

clusterEvalQ(cl, {
  library(vegan)
  library(tidyverse)
  library(picante)
  NF_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_prot_seqXst_PR2_3.csv", 
                                        sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
  NF_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_prot_treeNJ_PR2_3.tree")
})

NFPhylo_null_func <- function(x){
  # set up parameters and data
  community = NF_comm
  phylo = NF_phylo
  # performing comdistnt to calculate MNTD
  as.matrix(mntd(community, cophenetic(tipShuffle(phylo)), abundance.weighted = TRUE))
}

nsim.list <- sapply(1:999, list)
test <- parLapply(cl, nsim.list, NFPhylo_null_func)
test[[1000]] <- as.matrix(mntd(NF_comm, cophenetic(NF_phylo), abundance.weighted = TRUE))

Sys.time() - ini

NFPhylo_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(row.names(NF_comm)) %>%
  rename(obs = X1000)
write.table(NFPhylo_null, file = "/home/zac422/Desktop/OSCAR/Nulls/NF_Amntd_null_3.csv",
            sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
##### Nanoflagellate phylogenetic turnover ##########################################

###############################################################################################
##### Phylogenetic turnover / deterministic vs stochastic processes_Amntd #####################
###############################################################################################