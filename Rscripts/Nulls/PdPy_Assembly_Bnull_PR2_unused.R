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
##### Phylogenetic turnover / deterministic vs stochastic processes_Bmpd ######################
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
  as.matrix(comdist(community, cophenetic(tipShuffle(phylo)), abundance.weighted = TRUE))
}
nsim.list <- sapply(1:999, list)
test <- parLapply(cl, nsim.list, NFPhylo_null_func)
test[[1000]] <- as.matrix(comdist(NF_comm, cophenetic(NF_phylo), abundance.weighted = TRUE))

Sys.time() - ini

NFPhylo_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(expand.grid(row.names(NF_comm), row.names(NF_comm))) %>%
  rename(obs = X1000)
write.table(NFPhylo_null, file = "/home/zac422/Desktop/OSCAR/Nulls/NF_Bmpd_null_3.csv",
            sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
##### Nanoflagellate phylogenetic turnover ##########################################

###############################################################################################
##### Phylogenetic turnover / deterministic vs stochastic processes_Bmpd ######################
###############################################################################################

###############################################################################################
##### Phylogenetic turnover / deterministic vs stochastic processes_Bmntd #####################
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
  as.matrix(comdistnt(community, cophenetic(tipShuffle(phylo)), abundance.weighted = TRUE))
}
nsim.list <- sapply(1:999, list)
test <- parLapply(cl, nsim.list, NFPhylo_null_func)
test[[1000]] <- as.matrix(comdistnt(NF_comm, cophenetic(NF_phylo), abundance.weighted = TRUE))

Sys.time() - ini

NFPhylo_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(expand.grid(row.names(NF_comm), row.names(NF_comm))) %>%
  rename(obs = X1000)
write.table(NFPhylo_null, file = "/home/zac422/Desktop/OSCAR/Nulls/NF_Bmntd_null_3.csv",
            sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
##### Nanoflagellate phylogenetic turnover ##########################################

###############################################################################################
##### Phylogenetic turnover / deterministic vs stochastic processes_Bmntd #####################
###############################################################################################

###############################################################################################
##### OTU turnover / Dispersal limitation vs homogeneous dispersal ############################
###############################################################################################

##### Bacteria OTU turnover #########################################################
numCores <- detectCores()
numCores

cl <- makeCluster(numCores - 24)

clusterEvalQ(cl, {
  library(vegan)
  library(tidyverse)
  library(picante)
  Bac_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_Bac_seqXst_PR2_3.csv",
                                         sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
})

BacOTU_null_func <- function(x){
  # set up parameters and data
  community = Bac_comm
  n.mod = "independentswap"
  method = "chao"
  # performing comdistnt to calculate MNTD
  as.matrix(vegdist(randomizeMatrix(community, null.model = n.mod), method = method))
}

nsim.list <- sapply(1:999, list)
test <- parLapply(cl, nsim.list, BacOTU_null_func)
test[[1000]] <- as.matrix(vegdist(Bac_comm, null.model = "independentswap", method = "chao"))

BacOTU_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(expand.grid(row.names(Bac_comm), row.names(Bac_comm))) %>%
  rename(obs = X1000)
write.table(BacOTU_null, file = "/home/zac422/Desktop/OSCAR/Bac_chao_null_PR2.csv",
            sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
##### Bacteria OTU turnover #########################################################

##### Nanoflagellate OTU turnover ###################################################
numCores <- detectCores()
numCores

cl <- makeCluster(numCores - 24)

clusterEvalQ(cl, {
  library(vegan)
  library(tidyverse)
  library(picante)
  NF_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS/sECS_prot_seqXst_PR2_new.csv", 
                                        sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
})


NFOTU_null_func <- function(x){
  # set up parameters and data
  community = NF_comm
  n.mod = "independentswap"
  method = "chao"
  # performing comdistnt to calculate MNTD
  as.matrix(vegdist(randomizeMatrix(community, null.model = n.mod), method = method))
}

nsim.list <- sapply(1:999, list)
test <- parLapply(cl, nsim.list, NFOTU_null_func)
test[[1000]] <- as.matrix(vegdist(NF_comm, null.model = "independentswap", method = "chao"))

NFOTU_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(expand.grid(row.names(NF_comm), row.names(NF_comm))) %>%
  rename(obs = X1000)
write.table(NFOTU_null, file = "/home/zac422/Desktop/OSCAR/NF_chao_null_PR2.csv",
            sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
##### Nanoflagellate OTU turnover ###################################################

##### Hetero-trophic Nanoflagellate OTU turnover ####################################
numCores <- detectCores()
numCores

cl <- makeCluster(numCores - 24)

clusterEvalQ(cl, {
  library(vegan)
  library(tidyverse)
  library(picante)
  HNF_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS/sECS_HNF_seqXst_PR2_new.csv",
                                         sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
})

HNFOTU_null_func <- function(x){
  # set up parameters and data
  community = HNF_comm
  n.mod = "independentswap"
  method = "chao"
  # performing comdistnt to calculate MNTD
  as.matrix(vegdist(randomizeMatrix(community, null.model = n.mod), method = method))
}

nsim.list <- sapply(1:999, list)
test <- parLapply(cl, nsim.list, HNFOTU_null_func)
test[[1000]] <- as.matrix(vegdist(HNF_comm, null.model = "independentswap", method = "chao"))

HNFOTU_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(expand.grid(row.names(HNF_comm), row.names(HNF_comm))) %>%
  rename(obs = X1000)
write.table(HNFOTU_null, file = "/home/zac422/Desktop/OSCAR/HNF_chao_null_PR2.csv",
            sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
##### Hetero-trophic Nanoflagellate OTU turnover ####################################

###############################################################################################
##### OTU turnover / Dispersal limitation vs homogeneous dispersal ############################
###############################################################################################

