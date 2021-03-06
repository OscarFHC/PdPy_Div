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
###############################################################################################
##### Loading data ############################################################################
###############################################################################################

###############################################################################################
##### Phylogenetic turnover / deterministic vs stochastic processes_Bmpd ######################
###############################################################################################

##### Bacteria phylogenetic turnover ################################################
ini <- Sys.time()
numCores <- detectCores()
numCores

cl <- makeCluster(numCores - 24)

clusterEvalQ(cl, {
  library(vegan)
  #library(tidyverse)
  library(picante)
  Bac_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_Bac_seqXst_PR2_3.csv",
                                         sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
  Bac_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_Bac_treeNJ_PR2_3.tree")
  
})

BacPhylo_null_func <- function(x){
  #set up parameters and data
  community = Bac_comm
  phylo = Bac_phylo
  # performing comdistnt to calculate MNTD
  as.matrix(comdist(community, cophenetic(tipShuffle(phylo)), abundance.weighted = TRUE))
}
nsim.list <- sapply(1:999, list)
test <- parLapply(cl, nsim.list, BacPhylo_null_func)
test[[1000]] <- as.matrix(comdist(Bac_comm, cophenetic(Bac_phylo), abundance.weighted = TRUE))

Sys.time() - ini

BacPhylo_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(expand.grid(row.names(Bac_comm), row.names(Bac_comm))) %>%
  rename(obs = X1000)
write.table(BacPhylo_null, file = "/home/zac422/Desktop/OSCAR/Nulls/Bac_Bmpd_null_3.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
##### Bacteria phylogenetic turnover ################################################

###############################################################################################
##### Phylogenetic turnover / deterministic vs stochastic processes_Bmpd ######################
###############################################################################################

###############################################################################################
##### Phylogenetic turnover / deterministic vs stochastic processes_Bmntd #####################
###############################################################################################

##### Bacteria phylogenetic turnover ################################################
ini <- Sys.time()
numCores <- detectCores()
numCores

cl <- makeCluster(numCores - 24)

clusterEvalQ(cl, {
  library(vegan)
  #library(tidyverse)
  library(picante)
  Bac_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_Bac_seqXst_PR2_3.csv",
                                         sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
  Bac_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_Bac_treeNJ_PR2_3.tree")
})

BacPhylo_null_func <- function(x){
  #set up parameters and data
  community = Bac_comm
  phylo = Bac_phylo
  # performing comdistnt to calculate MNTD
  as.matrix(comdistnt(community, cophenetic(tipShuffle(phylo)), abundance.weighted = TRUE))
}
nsim.list <- sapply(1:999, list)
test <- parLapply(cl, nsim.list, BacPhylo_null_func)
test[[1000]] <- as.matrix(comdistnt(Bac_comm, cophenetic(Bac_phylo), abundance.weighted = TRUE))

Sys.time() - ini

BacPhylo_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(expand.grid(row.names(Bac_comm), row.names(Bac_comm))) %>%
  rename(obs = X1000)
write.table(BacPhylo_null, file = "/home/zac422/Desktop/OSCAR/Nulls/Bac_Bmntd_null_3.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
##### Bacteria phylogenetic turnover ################################################

###############################################################################################
##### Phylogenetic turnover / deterministic vs stochastic processes_Bmntd #####################
###############################################################################################

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
###############################################################################################
##### Loading data ############################################################################
###############################################################################################

###############################################################################################
##### Phylogenetic turnover / deterministic vs stochastic processes_Ampd #######################
###############################################################################################

##### Bacteria phylogenetic turnover ################################################
ini <- Sys.time()
numCores <- detectCores()
numCores
cl <- makeCluster(numCores - 24)

clusterEvalQ(cl, {
  library(vegan)
  #library(tidyverse)
  library(picante)
  Bac_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_Bac_seqXst_PR2_3.csv",
                                         sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
  Bac_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_Bac_treeNJ_PR2_3.tree")
})

BacPhylo_null_func <- function(x){
  #set up parameters and data
  community = Bac_comm
  phylo = Bac_phylo
  # performing comdistnt to calculate MNTD
  as.matrix(mpd(Bac_comm, cophenetic(tipShuffle(Bac_phylo)), abundance.weighted = TRUE))
}
nsim.list <- sapply(1:999, list)
test <- parLapply(cl, nsim.list, BacPhylo_null_func)
test[[1000]] <- as.matrix(mpd(Bac_comm, cophenetic(Bac_phylo), abundance.weighted = TRUE))

Sys.time() - ini

BacPhylo_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(row.names(Bac_comm)) %>%
  rename(obs = X1000)
write.table(BacPhylo_null, file = "/home/zac422/Desktop/OSCAR/Nulls/Bac_Ampd_null_3.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
##### Bacteria phylogenetic turnover ################################################

###############################################################################################
##### Phylogenetic turnover / deterministic vs stochastic processes_Ampd ######################
###############################################################################################

###############################################################################################
##### Phylogenetic turnover / deterministic vs stochastic processes_Amntd #####################
###############################################################################################

##### Bacteria phylogenetic turnover ################################################
ini <- Sys.time()
numCores <- detectCores()
numCores
cl <- makeCluster(numCores - 24)

clusterEvalQ(cl, {
  library(vegan)
  #library(tidyverse)
  library(picante)
  Bac_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_Bac_seqXst_PR2_3.csv",
                                         sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
  Bac_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_Bac_treeNJ_PR2_3.tree")
})

BacPhylo_null_func <- function(x){
  #set up parameters and data
  community = Bac_comm
  phylo = Bac_phylo
  # performing comdistnt to calculate MNTD
  as.matrix(mntd(Bac_comm, cophenetic(tipShuffle(Bac_phylo)), abundance.weighted = TRUE))
}

nsim.list <- sapply(1:999, list)
test <- parLapply(cl, nsim.list, BacPhylo_null_func)
test[[1000]] <- as.matrix(mntd(Bac_comm, cophenetic(Bac_phylo), abundance.weighted = TRUE))

Sys.time() - ini

BacPhylo_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(row.names(Bac_comm)) %>%
  rename(obs = X1000)
write.table(BacPhylo_null, file = "/home/zac422/Desktop/OSCAR/Nulls/Bac_Amntd_null_3.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
##### Bacteria phylogenetic turnover ################################################

###############################################################################################
##### Phylogenetic turnover / deterministic vs stochastic processes_Amntd #####################
###############################################################################################
