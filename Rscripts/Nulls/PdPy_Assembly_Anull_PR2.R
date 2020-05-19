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

if (!require(mgcv)) {
  install.packages("mgcv", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(mgcv)
}else{library(mgcv)}
###############################################################################################
##### Loading packages ########################################################################
###############################################################################################

###############################################################################################
##### Loading data ############################################################################
###############################################################################################
Bac_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_4/sECS_Bac_seqXst_PR2_4.csv",
                                       sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
Bac_ra_comm <- Bac_comm / rowSums(Bac_comm)
Bac_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_4/sECS_Bac_treeNJ_PR2_4.tree")

HNF_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_4/sECS_HNF_seqXst_PR2_4.csv",
                                       sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
HNF_ra_comm <- HNF_comm / rowSums(HNF_comm)
HNF_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_4/sECS_HNF_treeNJ_PR2_4.tree")

# NF_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_prot_seqXst_PR2_3.csv",
#                                       sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
# NF_ra_comm <- NF_comm / rowSums(NF_comm)
# NF_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_prot_treeNJ_PR2_3.tree")
###############################################################################################
##### Loading data ############################################################################
###############################################################################################

###############################################################################################
##### Predator-prey diversity correlation under neutral processes #############################
###############################################################################################
# Bac_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_Bac_seqXst_PR2_3.csv",
#                                        sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
# HNF_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_HNF_seqXst_PR2_3.csv",
#                                        sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
# Bac_A <- iNEXT(t(Bac_comm), q = 0, datatype = "abundance", size = max(colSums(HNF_comm)) + 100000)$AsyEst %>% 
#   select(Site, Diversity, Estimator) %>% 
#   spread(Diversity, Estimator) %>%
#   rename(Bac_q0 = "Species richness", Bac_q1 = "Shannon diversity", Bac_q2 = "Simpson diversity") %>%
#   mutate(Site = rownames(Bac_comm))
# HNF_A <- iNEXT(t(HNF_comm), q = 0, datatype = "abundance", size = max(colSums(HNF_comm)) + 100000)$AsyEst %>% 
#   select(Site, Diversity, Estimator) %>% 
#   spread(Diversity, Estimator) %>%
#   rename(HNF_q0 = "Species richness", HNF_q1 = "Shannon diversity", HNF_q2 = "Simpson diversity") %>%
#   mutate(Site = rownames(HNF_comm))
# 
# HNF_Bac_A <- Bac_A %>%
#   inner_join(HNF_A, by = c("Site" = "Site")) %>%
#   mutate(ln.Bac_q0 = log(Bac_q0),
#          ln.HNF_q0 = log(HNF_q0),
#          ln.Bac_q1 = log(Bac_q1),
#          ln.HNF_q1 = log(HNF_q1),
#          ln.Bac_q2 = log(Bac_q2),
#          ln.HNF_q2 = log(HNF_q2))
# Bac_comm <- Bac_comm[which(rownames(Bac_comm) %in% HNF_Bac_A$Site), ]
# q0coef <- summary(gam(ln.HNF_q0 ~ ln.Bac_q0, data = HNF_Bac_A))$p.coef[2]
# q1coef <- summary(gam(ln.HNF_q1 ~ ln.Bac_q1, data = HNF_Bac_A))$p.coef[2]
# q2coef <- summary(gam(ln.HNF_q2 ~ ln.Bac_q2, data = HNF_Bac_A))$p.coef[2]
# 
# Bac_SpList <- colnames(Bac_comm)
# Bac_CommSize <- rowSums(Bac_comm)
# Bac_PoolProb <- colSums(Bac_comm) / sum(Bac_comm)
# HNF_SpList <- colnames(HNF_comm)
# HNF_CommSize <- rowSums(HNF_comm)
# HNF_PoolProb <- colSums(HNF_comm) / sum(HNF_comm)
# 
# 
# 
# for (i in 1:999){
#   # ini <- Sys.time()
#   Bac_RandComm <- as.data.frame(matrix(0, nrow(Bac_comm), ncol(Bac_comm)))
#   rownames(Bac_RandComm) <- rownames(Bac_comm)
#   colnames(Bac_RandComm) <- colnames(Bac_comm)
#   HNF_RandComm <- as.data.frame(matrix(0, nrow(HNF_comm), ncol(HNF_comm)))
#   rownames(HNF_RandComm) <- rownames(HNF_comm)
#   colnames(HNF_RandComm) <- colnames(HNF_comm)
#   
#   for (j in 1:nrow(Bac_RandComm)){
#     Bac <- sample(Bac_SpList, Bac_CommSize[j], prob = Bac_PoolProb, replace = TRUE) 
#     HNF <- sample(HNF_SpList, HNF_CommSize[j], prob = HNF_PoolProb, replace = TRUE) 
#     
#     for (m in 1:ncol(Bac_RandComm)){
#       Bac_RandComm[j, m] <- length(which(Bac == colnames(Bac_RandComm)[m]))
#     }
#     
#     for (n in 1:ncol(HNF_RandComm)){
#       HNF_RandComm[j, n] <- length(which(HNF == colnames(HNF_RandComm)[n]))
#     }
#   }
#   
#   Bac_RandA <- iNEXT(t(Bac_RandComm), q = 0, datatype = "abundance", size = max(colSums(Bac_RandComm)) + 100000)$AsyEst %>% 
#     select(Site, Diversity, Estimator) %>% 
#     spread(Diversity, Estimator) %>%
#     rename(Bac_q0 = "Species richness", Bac_q1 = "Shannon diversity", Bac_q2 = "Simpson diversity") %>%
#     mutate(Site = rownames(Bac_RandComm))
#   HNF_RandA <- iNEXT(t(HNF_RandComm), q = 0, datatype = "abundance", size = max(colSums(HNF_RandComm)) + 100000)$AsyEst %>% 
#     select(Site, Diversity, Estimator) %>% 
#     spread(Diversity, Estimator) %>%
#     rename(HNF_q0 = "Species richness", HNF_q1 = "Shannon diversity", HNF_q2 = "Simpson diversity") %>%
#     mutate(Site = rownames(HNF_RandComm))
#   
#   HNF_Bac_RandA <- Bac_RandA %>%
#     inner_join(HNF_RandA, by = c("Site" = "Site")) %>%
#     mutate(ln.Bac_q0 = log(Bac_q0),
#            ln.HNF_q0 = log(HNF_q0),
#            ln.Bac_q1 = log(Bac_q1),
#            ln.HNF_q1 = log(HNF_q1),
#            ln.Bac_q2 = log(Bac_q2),
#            ln.HNF_q2 = log(HNF_q2))
#   
#   q0coef <- c(q0coef, summary(gam(ln.HNF_q0 ~ ln.Bac_q0, data = HNF_Bac_RandA))$p.coef[2])
#   q1coef <- c(q1coef, summary(gam(ln.HNF_q1 ~ ln.Bac_q1, data = HNF_Bac_RandA))$p.coef[2])
#   q2coef <- c(q2coef, summary(gam(ln.HNF_q2 ~ ln.Bac_q2, data = HNF_Bac_RandA))$p.coef[2])
#   # Sys.time() - ini
# }
# 
# write.table(cbind(q0coef, q1coef, q2coef), file = "D:/Research/PdPy_Div_Results/nulls_PR2_4/Neutral_ACorr_null_4.csv", 
#             sep = ",", col.names = TRUE, row.names = FALSE)
###############################################################################################
##### Predator-prey diversity correlation under neutral processes #############################
###############################################################################################

###############################################################################################
##### Predator-prey diversity correlation under neutral processes #############################
###############################################################################################
Bac_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_4/sECS_Bac_seqXst_PR2_4.csv",
                                       sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
HNF_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_4/sECS_HNF_seqXst_PR2_4.csv",
                                       sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
Bac_A <- iNEXT(t(Bac_comm), q = 0, datatype = "abundance", size = max(colSums(HNF_comm)) + 100000)$AsyEst %>%
  select(Site, Diversity, Estimator) %>%
  spread(Diversity, Estimator) %>%
  rename(Bac_q0 = "Species richness", Bac_q1 = "Shannon diversity", Bac_q2 = "Simpson diversity") %>%
  mutate(Site = rownames(Bac_comm))
HNF_A <- iNEXT(t(HNF_comm), q = 0, datatype = "abundance", size = max(colSums(HNF_comm)) + 100000)$AsyEst %>%
  select(Site, Diversity, Estimator) %>%
  spread(Diversity, Estimator) %>%
  rename(HNF_q0 = "Species richness", HNF_q1 = "Shannon diversity", HNF_q2 = "Simpson diversity") %>%
  mutate(Site = rownames(HNF_comm))

HNF_Bac_A <- Bac_A %>%
  inner_join(HNF_A, by = c("Site" = "Site")) %>%
  mutate(ln.Bac_q0 = log(Bac_q0),
         ln.HNF_q0 = log(HNF_q0),
         ln.Bac_q1 = log(Bac_q1),
         ln.HNF_q1 = log(HNF_q1),
         ln.Bac_q2 = log(Bac_q2),
         ln.HNF_q2 = log(HNF_q2))
Bac_comm <- Bac_comm[which(rownames(Bac_comm) %in% HNF_Bac_A$Site), ]
q0coef <- cor.test(HNF_Bac_A$ln.HNF_q0, HNF_Bac_A$ln.Bac_q0)$estimate
q1coef <- cor.test(HNF_Bac_A$ln.HNF_q1, HNF_Bac_A$ln.Bac_q1)$estimate
q2coef <- cor.test(HNF_Bac_A$ln.HNF_q2, HNF_Bac_A$ln.Bac_q2)$estimate

# now randomly shuffle HNF community
for (i in 1:999){
  HNF_RandComm <- HNF_comm[sample(nrow(HNF_comm)), ]
  
  HNF_RandA <- iNEXT(t(HNF_RandComm), q = 0, datatype = "abundance", size = max(colSums(HNF_RandComm)) + 100000)$AsyEst %>%
    select(Site, Diversity, Estimator) %>%
    spread(Diversity, Estimator) %>%
    rename(HNF_q0 = "Species richness", HNF_q1 = "Shannon diversity", HNF_q2 = "Simpson diversity") %>%
    mutate(Site = rownames(HNF_comm))
  
  HNF_Bac_RandA <- Bac_A %>%
    inner_join(HNF_RandA, by = c("Site" = "Site")) %>%
    mutate(ln.Bac_q0 = log(Bac_q0),
           ln.HNF_q0 = log(HNF_q0),
           ln.Bac_q1 = log(Bac_q1),
           ln.HNF_q1 = log(HNF_q1),
           ln.Bac_q2 = log(Bac_q2),
           ln.HNF_q2 = log(HNF_q2))
  
  q0coef <- c(q0coef, cor.test(HNF_Bac_RandA$ln.HNF_q0, HNF_Bac_RandA$ln.Bac_q0)$estimate)
  q1coef <- c(q1coef, cor.test(HNF_Bac_RandA$ln.HNF_q1, HNF_Bac_RandA$ln.Bac_q1)$estimate)
  q2coef <- c(q2coef, cor.test(HNF_Bac_RandA$ln.HNF_q2, HNF_Bac_RandA$ln.Bac_q2)$estimate)
}

write.table(cbind(q0coef, q1coef, q2coef), file = "D:/Research/PdPy_Div/Anulls_PR2_4/Neutral_ACorr_null.csv",
            sep = ",", col.names = TRUE, row.names = FALSE)
###############################################################################################
##### Predator-prey diversity correlation under neutral processes #############################
###############################################################################################
  

###############################################################################################
##### Phylogenetic turnover / deterministic vs stochastic processes_Ampd ######################
###############################################################################################

##### Bacteria phylogenetic turnover ################################################
ini <- Sys.time()
numCores <- detectCores()
numCores
cl <- makeCluster(numCores - 2)

clusterEvalQ(cl, {
  library(vegan)
  #library(tidyverse)
  library(picante)
  Bac_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_4/sECS_Bac_seqXst_PR2_4.csv",
                                         sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
  Bac_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_4/sECS_Bac_treeNJ_PR2_4.tree")
})

BacPhylo_null_func <- function(x){
  #set up parameters and data
  community = Bac_comm
  phylo = Bac_phylo
  # performing comdistnt to calculate MPD
  as.matrix(mpd(Bac_comm, cophenetic(tipShuffle(Bac_phylo)), abundance.weighted = TRUE))
}
nsim.list <- sapply(1:999, list)
test <- parLapply(cl, nsim.list, BacPhylo_null_func)
test[[1000]] <- as.matrix(mpd(Bac_comm, cophenetic(Bac_phylo), abundance.weighted = TRUE))

Sys.time() - ini

BacPhylo_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(row.names(Bac_comm)) %>%
  rename(obs = X1000)
write.table(BacPhylo_null, file = "/home/zac422/Desktop/OSCAR/Nulls/Bac_Ampd_null_4.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
##### Bacteria phylogenetic turnover ################################################

##### Hetero-trophic Nanoflagellate phylogenetic turnover ###########################
ini <- Sys.time()
numCores <- detectCores()
numCores

cl <- makeCluster(numCores - 2)

clusterEvalQ(cl, {
  library(vegan)
  library(tidyverse)
  library(picante)
  HNF_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_4/sECS_HNF_seqXst_PR2_4.csv",
                                         sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
  HNF_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_4/sECS_HNF_treeNJ_PR2_4.tree")
})

HNFPhylo_null_func <- function(x){
  # set up parameters and data
  community = HNF_comm
  phylo = HNF_phylo
  # performing comdistnt to calculate MPD
  as.matrix(mpd(community, cophenetic(tipShuffle(phylo)), abundance.weighted = TRUE))
}
nsim.list <- sapply(1:999, list)
test <- parLapply(cl, nsim.list, HNFPhylo_null_func)
test[[1000]] <- as.matrix(mpd(HNF_comm, cophenetic(HNF_phylo), abundance.weighted = TRUE))

Sys.time() - ini

HNFPhylo_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(row.names(HNF_comm)) %>%
  rename(obs = X1000)
write.table(HNFPhylo_null, file = "D:/Research/PdPy_Div/Anulls_PR2_4/HNF_Ampd_null_4.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
##### Hetero-trophic Nanoflagellate phylogenetic turnover ###########################

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
  Bac_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_4/sECS_Bac_seqXst_PR2_4.csv",
                                         sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
  Bac_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_4/sECS_Bac_treeNJ_PR2_4.tree")
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
write.table(BacPhylo_null, file = "/home/zac422/Desktop/OSCAR/Nulls/Bac_Amntd_null_4.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
##### Bacteria phylogenetic turnover ################################################

##### Hetero-trophic Nanoflagellate phylogenetic turnover ###########################
ini <- Sys.time()
numCores <- detectCores()
numCores

cl <- makeCluster(numCores - 2)

clusterEvalQ(cl, {
  library(vegan)
  library(tidyverse)
  library(picante)
  HNF_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_4/sECS_HNF_seqXst_PR2_4.csv",
                                         sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
  HNF_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_4/sECS_HNF_treeNJ_PR2_4.tree")
})

HNFPhylo_null_func <- function(x){
  # set up parameters and data
  community = HNF_comm
  phylo = HNF_phylo
  # performing comdistnt to calculate MNTD
  as.matrix(mntd(community, cophenetic(tipShuffle(phylo)), abundance.weighted = TRUE))
}
nsim.list <- sapply(1:999, list)
test <- parLapply(cl, nsim.list, HNFPhylo_null_func)
test[[1000]] <- as.matrix(mntd(HNF_comm, cophenetic(HNF_phylo), abundance.weighted = TRUE))

Sys.time() - ini

HNFPhylo_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(row.names(HNF_comm)) %>%
  rename(obs = X1000)
write.table(HNFPhylo_null, file = "D:/Research/PdPy_Div/Anulls_PR2_4/HNF_Amntd_null_4.csv", 
           sep = ",", col.names = TRUE, row.names = FALSE)
stopCluster(cl)
##### Hetero-trophic Nanoflagellate phylogenetic turnover ###########################

###############################################################################################
##### Phylogenetic turnover / deterministic vs stochastic processes_Amntd #####################
###############################################################################################
