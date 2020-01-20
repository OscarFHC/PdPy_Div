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
Bac_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/16s_seqXst.csv",
                                       sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
Bac_ra_comm <- Bac_comm / rowSums(Bac_comm)
Bac_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/treeNJ_16s.tree")

NF_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/prot_seqXst_PR2.csv", sep = ",", 
                                      header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
NF_ra_comm <- NF_comm / rowSums(NF_comm)
NF_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/prot_treeNJ_PR2.tree")

HNF_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/prot_seqXst_PR2.csv", sep = ",", 
                                       header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
HNF_ra_comm <- HNF_comm / rowSums(HNF_comm)
HNF_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/HNF_treeNJ_PR2.tree")


Vars <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_Vars.csv", sep = ",", 
                   header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
###############################################################################################
##### Loading data ############################################################################
###############################################################################################

###############################################################################################
##### OTU turnover / Dispersal limitation vs homogeneous dispersal ############################
###############################################################################################

##### Bacteria OTU turnover #########################################################
ini <- Sys.time()
numCores <- detectCores()
numCores

cl <- makeCluster(numCores)

clusterEvalQ(cl, {
  library(vegan)
  library(tidyverse)
  library(picante)
  Bac_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/16s_seqXst.csv", sep = ",", 
                                         header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
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

Sys.time() - ini

BacOTU_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(expand.grid(row.names(Bac_comm), row.names(Bac_comm))) %>%
  rename(obs = X1000)
write.table(BacOTU_null, file = "D:/Research/PdPy_Div_Results/Bac_chao_null.csv",
            sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
##### Bacteria OTU turnover #########################################################

##### Nanoflagellate OTU turnover ###################################################
ini <- Sys.time()
numCores <- detectCores()
numCores

cl <- makeCluster(numCores)

clusterEvalQ(cl, {
  library(vegan)
  library(tidyverse)
  library(picante)
  NF_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/prot_seqXst_PR2.csv", sep = ",", 
                                        header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
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

Sys.time() - ini

NFOTU_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(expand.grid(row.names(NF_comm), row.names(NF_comm))) %>%
  rename(obs = X1000)
write.table(NFOTU_null, file = "D:/Research/PdPy_Div_Results/NF_chao_null_PR2.csv",
            sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
##### Nanoflagellate OTU turnover ###################################################

##### Hetero-trophic Nanoflagellate OTU turnover ####################################
ini <- Sys.time()
numCores <- detectCores()
numCores

cl <- makeCluster(numCores)

clusterEvalQ(cl, {
  library(vegan)
  library(tidyverse)
  library(picante)
  HNF_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/HNF_seqXst_PR2.csv", sep = ",", 
                                        header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
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

Sys.time() - ini

HNFOTU_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(expand.grid(row.names(HNF_comm), row.names(HNF_comm))) %>%
  rename(obs = X1000)
write.table(HNFOTU_null, file = "D:/Research/PdPy_Div_Results/HNF_chao_null_PR2.csv",
            sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
##### Hetero-trophic Nanoflagellate OTU turnover ####################################

###############################################################################################
##### OTU turnover / Dispersal limitation vs homogeneous dispersal ############################
###############################################################################################
