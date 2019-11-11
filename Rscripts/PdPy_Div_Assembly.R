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
##### Loading nulls ###########################################################################
###############################################################################################

###############################################################################################
##### Loading nulls ###########################################################################
###############################################################################################


