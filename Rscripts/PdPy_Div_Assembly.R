###############################################################################################
##### Loading packages ########################################################################
###############################################################################################
if (!require(rmarkdown)) {
  install.packages("rmarkdown", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(rmarkdown)
}else{library(rmarkdown)}

if (!require(knitr)) {
  install.packages("knitr", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(knitr)
}else{library(knitr)}

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

# install.packages("iNEXT")
# 
# ## install the latest version from github
# install.packages('devtools')
# library(devtools)
# install_github('JohnsonHsieh/iNEXT')
# library(iNEXT)

if (!require(lavaan)) {
  install.packages("lavaan", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(lavaan)
}else{library(lavaan)}

if (!require(nlme)) {
  install.packages("nlme", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(nlme)
}else{library(nlme)}

if (!require(piecewiseSEM)) {
  install.packages("piecewiseSEM", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(piecewiseSEM)
}else{library(piecewiseSEM)}

if (!require(brms)) {
  install.packages("brms", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(brms)
}else{library(brms)}

if (!require(blavaan)) {
  install.packages("blavaan", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(blavaan)
}else{library(blavaan)}

if (!require(rstanarm)) {
  install.packages("rstanarm", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(rstanarm)
}else{library(rstanarm)}

library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

if (!require(loo)) {
  install.packages("loo", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(loo)
}else{library(loo)}

if (!require(abind)) {
  install.packages("abind", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(abind)
}else{library(abind)}

if (!require(ggplot2)) {
  install.packages("ggplot2", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(ggplot2)
}else{library(ggplot2)}

if (!require(viridis)) {
  install.packages("viridis", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(viridis)
}else{library(viridis)}

if (!require(GGally)) {
  install.packages("GGally", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(GGally)
}else{library(GGally)}

if (!require(cowplot)) {
  install.packages("cowplot", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(cowplot)
}else{library(cowplot)}
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

NF_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/18s_seqXst.csv", sep = ",", 
                                      header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
NF_ra_comm <- NF_comm / rowSums(NF_comm)
HNF_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/HNF_seqXst.csv", sep = ",", 
                                       header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
HNF_ra_comm <- HNF_comm / rowSums(HNF_comm)
NF_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/treeNJ_18s.tree")

Vars <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_Vars.csv", sep = ",", 
                   header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
###############################################################################################
##### Loading data ############################################################################
###############################################################################################

###############################################################################################
##### Phylogenetic turnover / deterministic vs stochastic processes ###########################
###############################################################################################

##### Bacteria phylogenetic turnover ################################################
numCores <- detectCores()
numCores

cl <- makeCluster(numCores)

clusterEvalQ(cl, {
  library(vegan)
  library(tidyverse)
  library(picante)
  Bac_comm <- as.data.frame(t(read.table(file = "D:/Research/PdPy_Div/data/16s_seqXst.csv", sep = ",", 
                                         header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
  Bac_phylo<- read.tree(file = "D:/Research/PdPy_Div/data/treeNJ_16s.tree")
  
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

BacPhylo_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(expand.grid(row.names(Bac_comm), row.names(Bac_comm))) %>%
  rename(obs = X1000)
write.table(BacOTU_null, file = "D:/Research/PdPy_Div_Results/Bac_MNTD_null.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
##### Bacteria phylogenetic turnover ################################################

##### Nanoflagellate phylogenetic turnover ##########################################
numCores <- detectCores()
numCores

cl <- makeCluster(numCores)

clusterEvalQ(cl, {
  library(vegan)
  library(tidyverse)
  library(picante)
  NF_comm <- as.data.frame(t(read.table(file = "D:/Research/PdPy_Div/data/18s_seqXst.csv", sep = ",", 
                                        header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
  NF_phylo<- read.tree(file = "D:/Research/PdPy_Div/data/treeNJ_18s.tree")
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

NFPhylo_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(expand.grid(row.names(NF_comm), row.names(NF_comm))) %>%
  rename(obs = X1000)
write.table(NFPhylo_null, file = "D:/Research/PdPy_Div_Results/NF_MNTD_null.csv",
            sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
##### Nanoflagellate phylogenetic turnover ##########################################

##### Hetero-trophic Nanoflagellate phylogenetic turnover ###########################
numCores <- detectCores()
numCores

cl <- makeCluster(numCores)

clusterEvalQ(cl, {
  library(vegan)
  library(tidyverse)
  library(picante)
  
  HNF_comm <- as.data.frame(t(read.table(file = "D:/Research/PdPy_Div/data/HNF_seqXst.csv", sep = ",", 
                                        header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
  HNF_phylo<- read.tree(file = "D:/Research/PdPy_Div/data/treeNJ_18s.tree")
})

HNFPhylo_null_func <- function(x){
  # set up parameters and data
  community = HNF_comm
  phylo = HNF_phylo
  # performing comdistnt to calculate MNTD
  as.matrix(comdistnt(community, cophenetic(tipShuffle(phylo)), abundance.weighted = TRUE))
}
# ini <- Sys.time()
nsim.list <- sapply(1:999, list)
test <- parLapply(cl, nsim.list, HNFPhylo_null_func)
#Sys.time() - ini
test[[1000]] <- as.matrix(comdistnt(HNF_comm, cophenetic(HNF_phylo), abundance.weighted = TRUE))

HNFPhylo_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(expand.grid(row.names(NF_comm), row.names(NF_comm))) %>%
  rename(obs = X1000)
write.table(HNFPhylo_null, file = "D:/Research/PdPy_Div_Results/HNF_MNTD_null.csv", 
           sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
##### Hetero-trophic Nanoflagellate phylogenetic turnover ###########################

###############################################################################################
##### Phylogenetic turnover / deterministic vs stochastic processes ###########################
###############################################################################################




# OTU turnover / Dispersal limitation vs homogeneous dispersal

```{r, Bacteria OTU turnover}
numCores <- detectCores()
numCores

cl <- makeCluster(numCores)

clusterEvalQ(cl, {
  library(vegan)
  library(tidyverse)
  library(picante)
  Bac_comm <- as.data.frame(t(read.table(file = "D:/Research/PdPy_Div/data/16s_seqXst.csv", sep = ",", 
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

BacOTU_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(expand.grid(row.names(Bac_comm), row.names(Bac_comm))) %>%
  rename(obs = X1000)
# write.table(BacOTU_null, file = "D:/Research/PdPy_Div/Results/Bac_chao_null.csv", 
#             sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
```

```{r, Nanoflagellate OTU turnover}
numCores <- detectCores()
numCores

cl <- makeCluster(numCores)

clusterEvalQ(cl, {
  library(vegan)
  library(tidyverse)
  library(picante)
  NF_comm <- as.data.frame(t(read.table(file = "D:/Research/PdPy_Div/data/18s_seqXst.csv", sep = ",", 
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

NFOTU_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(expand.grid(row.names(NF_comm), row.names(NF_comm))) %>%
  rename(obs = X1000)
write.table(NFOTU_null, file = "D:/Research/PdPy_Div/Results/NF_chao_null.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
```

```{r, Hetero-trophic Nanoflagellate OTU turnover}
numCores <- detectCores()
numCores

cl <- makeCluster(numCores)

clusterEvalQ(cl, {
  library(vegan)
  library(tidyverse)
  library(picante)
  NF_comm <- as.data.frame(t(read.table(file = "D:/Research/PdPy_Div/data/18s_seqXst.csv", sep = ",", 
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

NFOTU_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(expand.grid(row.names(NF_comm), row.names(NF_comm))) %>%
  rename(obs = X1000)
write.table(NFOTU_null, file = "D:/Research/PdPy_Div/Results/NF_chao_null.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
```