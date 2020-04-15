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

if (!require(mgcv)) {
  install.packages("mgcv", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(mgcv)
}else{library(mgcv)}

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

if (!require(semTools)) {
  install.packages("semTools", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(semTools)
}else{library(semTools)}

if (!require(psych)) {
  install.packages("psych", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(psych)
}else{library(psych)}

if (!require(plspm)) {
  install.packages("plspm", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(plspm)
}else{library(plspm)}

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

if (!require(lmodel2)) {
  install.packages("lmodel2", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(lmodel2)
}else{library(lmodel2)}

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
##### Loading functions #######################################################################
###############################################################################################
cor_fun <- function(data, mapping, method="pearson", ndp=2, sz=5, stars=TRUE, ...){
  
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  corr <- cor.test(x, y, method = method)
  est <- corr$estimate
  lb.size <- 6#sz* abs(est) 
  
  if(stars){
    stars <- c("***", "**", "*", "")[findInterval(corr$p.value, c(0, 0.001, 0.01, 0.05, 1))]
    lbl <- paste0(round(est, ndp), stars)
  }else{
    lbl <- round(est, ndp)
  }
  
  ggplot(data = data, mapping = mapping) + 
    annotate("text", x = mean(x, na.rm = TRUE), y = mean(y, na.rm = TRUE), label = lbl, size = lb.size, ...)+
    theme(panel.grid = element_blank())
}

fit_fun <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method = loess, fill = "red", color = "red", ...) +
    geom_smooth(method = lm, fill = "blue", color = "blue", ...)
  p
}
###############################################################################################
##### Loading functions #######################################################################
###############################################################################################

###############################################################################################
##### Loading data ############################################################################
###############################################################################################
Bac_comm <- as.data.frame(t(read.table(
  file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_Bac_seqXst_PR2_3.csv", 
  sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
Bac_ra_comm <- Bac_comm / rowSums(Bac_comm)
Bac_phylo<- read.tree(
  file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_Bac_treeNJ_PR2_3.tree")

HNF_comm <- as.data.frame(t(read.table(
  file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_HNF_seqXst_PR2_3.csv", 
  sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
HNF_ra_comm <- HNF_comm / rowSums(HNF_comm)
HNF_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_HNF_treeNJ_PR2_3.tree")

Vars <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_Vars.csv", sep = ",", 
                   header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
###############################################################################################
##### Loading data ############################################################################
###############################################################################################

###############################################################################################
##### Loading nulls ###########################################################################
###############################################################################################
Bac_Amntd_null <- read.table(file = "D:/Research/PdPy_Div_Results/nulls_PR2_3/Bac_Amntd_null_3.csv", sep = ",", 
                             header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
Bac_Ampd_null <- read.table(file = "D:/Research/PdPy_Div_Results/nulls_PR2_3/Bac_Ampd_null_3.csv", sep = ",", 
                            header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
Bac_Bmntd_null <- read.table(file = "D:/Research/PdPy_Div_Results/nulls_PR2_3/Bac_Bmntd_null_3.csv", sep = ",", 
                             header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
Bac_Bmpd_null <- read.table(file = "D:/Research/PdPy_Div_Results/nulls_PR2_3/Bac_Bmpd_null_3.csv", sep = ",", 
                            header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
Bac_Amntd <- Bac_Amntd_null %>% 
  select(c(obs, Site)) %>%
  mutate(Amntd_null_mean = apply(Bac_Amntd_null[, !names(Bac_Amntd_null) %in% c("obs", "Site")], 1, mean),
         Amntd_null_sd = apply(Bac_Amntd_null[, !names(Bac_Amntd_null) %in% c("obs", "Site")], 1, sd),
         Bac_Amntd_select = (obs - Amntd_null_mean) / Amntd_null_sd,
         Bac_Amntd_select_p = pnorm(-abs(Bac_Amntd_select), 0, 1))
Bac_Ampd <- Bac_Ampd_null %>% 
  select(c(obs, Site)) %>%
  mutate(Ampd_null_mean = apply(Bac_Ampd_null[, !names(Bac_Ampd_null) %in% c("obs", "Site")], 1, mean),
         Ampd_null_sd = apply(Bac_Ampd_null[, !names(Bac_Ampd_null) %in% c("obs", "Site")], 1, sd),
         Bac_Ampd_select = (obs - Ampd_null_mean) / Ampd_null_sd,
         Bac_Ampd_select_p = pnorm(-abs(Bac_Ampd_select), 0, 1))
Bac_Bmntd <- Bac_Bmntd_null %>% 
  select(c(obs, Var1, Var2)) %>%
  mutate(Bmntd_null_mean = apply(Bac_Bmntd_null[, !names(Bac_Bmntd_null) %in% c("obs", "Var1", "Var2")], 1, mean),
         Bmntd_null_sd = apply(Bac_Bmntd_null[, !names(Bac_Bmntd_null) %in% c("obs", "Var1", "Var2")], 1, sd),
         Bac_Bmntd_select_strength = (obs - Bmntd_null_mean) / Bmntd_null_sd,
         Bac_Bmntd_select_p = pnorm(-abs(Bac_Bmntd_select_strength), 0, 1))
Bac_Bmpd <- Bac_Bmpd_null %>% 
  select(c(obs, Var1, Var2)) %>%
  mutate(Bmpd_null_mean = apply(Bac_Bmpd_null[, !names(Bac_Bmpd_null) %in% c("obs", "Var1", "Var2")], 1, mean),
         Bmpd_null_sd = apply(Bac_Bmpd_null[, !names(Bac_Bmpd_null) %in% c("obs", "Var1", "Var2")], 1, sd),
         Bac_Bmpd_select_strength = (obs - Bmpd_null_mean) / Bmpd_null_sd,
         Bac_Bmpd_select_p = pnorm(-abs(Bac_Bmpd_select_strength), 0, 1))

HNF_Amntd_null <- read.table(file = "D:/Research/PdPy_Div_Results/nulls_PR2_3/HNF_Amntd_null_3.csv", sep = ",", 
                             header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
HNF_Ampd_null <- read.table(file = "D:/Research/PdPy_Div_Results/nulls_PR2_3/HNF_Ampd_null_3.csv", sep = ",", 
                            header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
HNF_Bmntd_null <- read.table(file = "D:/Research/PdPy_Div_Results/nulls_PR2_3/HNF_Bmntd_null_3.csv", sep = ",", 
                             header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
HNF_Bmpd_null <- read.table(file = "D:/Research/PdPy_Div_Results/nulls_PR2_3/HNF_Bmpd_null_3.csv", sep = ",", 
                            header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
HNF_Amntd <- HNF_Amntd_null %>% 
  select(c(obs, Site)) %>%
  mutate(Amntd_null_mean = apply(HNF_Amntd_null[, !names(HNF_Amntd_null) %in% c("obs", "Site")], 1, mean),
         Amntd_null_sd = apply(HNF_Amntd_null[, !names(HNF_Amntd_null) %in% c("obs", "Site")], 1, sd),
         HNF_Amntd_select = (obs - Amntd_null_mean) / Amntd_null_sd,
         HNF_Amntd_select_p = pnorm(-abs(HNF_Amntd_select), 0, 1))
HNF_Ampd <- HNF_Ampd_null %>% 
  select(c(obs, Site)) %>%
  mutate(Ampd_null_mean = apply(HNF_Ampd_null[, !names(HNF_Ampd_null) %in% c("obs", "Site")], 1, mean),
         Ampd_null_sd = apply(HNF_Ampd_null[, !names(HNF_Ampd_null) %in% c("obs", "Site")], 1, sd),
         HNF_Ampd_select = (obs - Ampd_null_mean) / Ampd_null_sd,
         HNF_Ampd_select_p = pnorm(-abs(HNF_Ampd_select), 0, 1))
HNF_Bmntd <- HNF_Bmntd_null %>% select(c(obs, Var1, Var2)) %>%
  mutate(Bmntd_null_mean = apply(HNF_Bmntd_null[, !names(HNF_Bmntd_null) %in% c("obs", "Var1", "Var2")], 1, mean),
         Bmntd_null_sd = apply(HNF_Bmntd_null[, !names(HNF_Bmntd_null) %in% c("obs", "Var1", "Var2")], 1, sd),
         HNF_Bmntd_select_strength = (obs - Bmntd_null_mean) / Bmntd_null_sd,
         HNF_Bmntd_select_p = pnorm(-abs(HNF_Bmntd_select_strength), 0, 1))
HNF_Bmpd <- HNF_Bmpd_null %>% select(c(obs, Var1, Var2)) %>%
  mutate(Bmpd_null_mean = apply(HNF_Bmpd_null[, !names(HNF_Bmpd_null) %in% c("obs", "Var1", "Var2")], 1, mean),
         Bmpd_null_sd = apply(HNF_Bmpd_null[, !names(HNF_Bmpd_null) %in% c("obs", "Var1", "Var2")], 1, sd),
         HNF_Bmpd_select_strength = (obs - Bmpd_null_mean) / Bmpd_null_sd,
         HNF_Bmpd_select_p = pnorm(-abs(HNF_Bmpd_select_strength), 0, 1))
###############################################################################################
##### Loading nulls ###########################################################################
###############################################################################################

###############################################################################################
##### Phylogenetic signal #####################################################################
###############################################################################################
##### Bacteria ##########
# Bac_Niche <- as.data.frame(matrix(0, ncol(Bac_comm), 8)) %>%
#   #mutate(Bac_ID = names(Bac_comm)) %>%
#   rename(Temp = V1, Sal = V2, PAR = V3, NO2 = V4, NO3 = V5, DIN = V6, PO3 = V7, Chla = V8) 
# 
# for(i in 1:ncol(Bac_comm)){ 
#   Bac_Niche[i, "Temp"] = sum(Bac_comm[,names(Bac_comm)[i]] * Vars$Temp) / sum(Bac_comm[,names(Bac_comm)[i]])
#   Bac_Niche[i, "Sal"] = sum(Bac_comm[,names(Bac_comm)[i]] * Vars$Sal) / sum(Bac_comm[,names(Bac_comm)[i]])
#   Bac_Niche[i, "PAR"] = sum(Bac_comm[,names(Bac_comm)[i]] * Vars$PAR) / sum(Bac_comm[,names(Bac_comm)[i]])
#   Bac_Niche[i, "NO2"] = sum(Bac_comm[,names(Bac_comm)[i]] * Vars$NO2) / sum(Bac_comm[,names(Bac_comm)[i]])
#   Bac_Niche[i, "NO3"] = sum(Bac_comm[,names(Bac_comm)[i]] * Vars$NO3) / sum(Bac_comm[,names(Bac_comm)[i]])
#   Bac_Niche[i, "DIN"] = sum(Bac_comm[,names(Bac_comm)[i]] * Vars$DIN) / sum(Bac_comm[,names(Bac_comm)[i]])
#   Bac_Niche[i, "PO3"] = sum(Bac_comm[,names(Bac_comm)[i]] * Vars$PO3) / sum(Bac_comm[,names(Bac_comm)[i]])
#   Bac_Niche[i, "Chla"] = sum(Bac_comm[,names(Bac_comm)[i]] * Vars$Chla) / sum(Bac_comm[,names(Bac_comm)[i]])
# }
# 
# Bac_Niche <- as.data.frame(scale(Bac_Niche, center = TRUE, scale = TRUE))
# Bac_Niche_Dist <- vegdist(Bac_Niche, method = "euclidean")
# Bac_Phylo_Dist <- cophenetic(Bac_phylo)
# Bac_PhySig <- mantel.correlog(Bac_Niche_Dist, Bac_Phylo_Dist, n.class = 50)
# Bac_PhySig_fig <- plot(Bac_PhySig)
# ggsave(plot(Bac_PhySig), file = "D:/Research/PdPy_Div_Results/p_Bac_PhySig.png", dpi = 600) 
##### Bacteria ##########
 
##### HNF ##########
# HNF_noZeroComm <- HNF_comm[, which(colSums(HNF_comm) != 0)]
# HNF_Niche <- as.data.frame(matrix(0, ncol(HNF_noZeroComm), 8)) %>%
#   # mutate(HNF_ID = names(HNF_comm)) %>%
#   rename(Temp = V1, Sal = V2, PAR = V3, NO2 = V4, NO3 = V5, DIN = V6, PO3 = V7, Chla = V8) 
# Vars_HNF <- Vars[Vars$SampleID %in% row.names(HNF_noZeroComm),]
# 
# for(i in 1:ncol(HNF_noZeroComm)){ 
#   HNF_Niche[i, "Temp"] = sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]] * Vars_HNF$Temp) / sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]])
#   HNF_Niche[i, "Sal"] = sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]] * Vars_HNF$Sal) / sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]])
#   HNF_Niche[i, "PAR"] = sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]] * Vars_HNF$PAR) / sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]])
#   HNF_Niche[i, "NO2"] = sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]] * Vars_HNF$NO2) / sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]])
#   HNF_Niche[i, "NO3"] = sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]] * Vars_HNF$NO3) / sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]])
#   HNF_Niche[i, "DIN"] = sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]] * Vars_HNF$DIN) / sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]])
#   HNF_Niche[i, "PO3"] = sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]] * Vars_HNF$PO3) / sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]])
#   HNF_Niche[i, "Chla"] = sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]] * Vars_HNF$Chla) / sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]])
# }
# 
# HNF_Niche <- as.data.frame(scale(HNF_Niche, center = TRUE, scale = TRUE))
# HNF_Niche_Dist <- vegdist(HNF_Niche, method = "euclidean")
# HNF_Phylo_Dist <- cophenetic(drop.tip(HNF_phylo, which(colSums(HNF_comm) == 0)))
# HNF_PhySig <- mantel.correlog(HNF_Niche_Dist, HNF_Phylo_Dist, n.class = 50)
# HNF_PhySig_fig <- plot(HNF_PhySig)
# ggsave(plot(HNF_PhySig), file = "D:/Research/PdPy_Div_Results/p_HNF_PhySig.png", dpi = 600) 
##### HNF ##########
###############################################################################################
##### Phylogenetic signal #####################################################################
###############################################################################################

###############################################################################################
##### sampling maps ###########################################################################
###############################################################################################
St_sECS = read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/raw/St_Location.csv", 
                     header = T, fill = T, sep = ",")
library(ggmap)
library(mapproj)
bbox <- ggmap::make_bbox(Lon, Lat, St_sECS, f = 0.1)
map <- get_map(
  location = bbox #c(lon=-85, lat=44) # google search string
  #, zoom = 7 # larger is closer
  #, maptype = "terrian" # map type
  , source = "google"
)

map = ggmap(map) + 
  geom_point(data = St_sECS[1:6,], aes(x = Lon, y = Lat), size = 4) + 
  labs(x="Longitude", y="Latitude") + 
  geom_text(data = St_sECS[1:6,], aes(x = Lon, y = Lat, label = Station), hjust = -0.4, vjust = -0.3, size = 6) + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16))
# ggsave(map, file = "D:/Research/PdPy_Div_Results/Figs/PR2_3/p_map.png", dpi = 600) 
###############################################################################################
##### sampling maps ###########################################################################
###############################################################################################

###############################################################################################
##### Preping data ############################################################################
###############################################################################################
Bac_A <- iNEXT(t(Bac_comm), q = 0, datatype = "abundance", size = max(colSums(Bac_comm)) + 100000)$AsyEst %>% 
  select(Site, Diversity, Estimator) %>% 
  spread(Diversity, Estimator) %>%
  rename(Bac_q0 = "Species richness", Bac_q1 = "Shannon diversity", Bac_q2 = "Simpson diversity") %>%
  mutate(Site = rownames(Bac_comm))

HNF_A <- iNEXT(t(HNF_comm), q = 0, datatype = "abundance", size = max(colSums(HNF_comm)) + 100000)$AsyEst %>% 
  select(Site, Diversity, Estimator) %>% 
  spread(Diversity, Estimator) %>%
  rename(HNF_q0 = "Species richness", HNF_q1 = "Shannon diversity", HNF_q2 = "Simpson diversity") %>%
  mutate(Site = rownames(HNF_comm))

Bac_Bmntd_selec <- Bac_Bmntd %>%
  mutate(Cr_V1 = substr(Var1, start = 1, stop = 9),
         Cr_V2 = substr(Var2, start = 1, stop = 9)) %>%
  filter(Cr_V1 == Cr_V2) %>%
  #filter(Bac_select_p < 0.05) %>%
  group_by(Var2) %>%
  summarize(Bac_Bmntd_select = mean(Bac_Bmntd_select_strength, na.rm = TRUE))

Bac_Bmpd_selec <- Bac_Bmpd %>%
  mutate(Cr_V1 = substr(Var1, start = 1, stop = 9),
         Cr_V2 = substr(Var2, start = 1, stop = 9)) %>%
  filter(Cr_V1 == Cr_V2) %>%
  #filter(Bac_select_p < 0.05) %>%
  group_by(Var2) %>%
  summarize(Bac_Bmpd_select = mean(Bac_Bmpd_select_strength, na.rm = TRUE))

HNF_Bmntd_selec <- HNF_Bmntd %>%
  mutate(Cr_V1 = substr(Var1, start = 1, stop = 9),
         Cr_V2 = substr(Var2, start = 1, stop = 9)) %>%
  filter(Cr_V1 == Cr_V2) %>%
  #filter(HNF_select_p < 0.05) %>%
  group_by(Var2) %>%
  summarize(HNF_Bmntd_select = mean(HNF_Bmntd_select_strength, na.rm = TRUE))

HNF_Bmpd_selec <- HNF_Bmpd %>%
  mutate(Cr_V1 = substr(Var1, start = 1, stop = 9),
         Cr_V2 = substr(Var2, start = 1, stop = 9)) %>%
  filter(Cr_V1 == Cr_V2) %>%
  #filter(HNF_select_p < 0.05) %>%
  group_by(Var2) %>%
  summarize(HNF_Bmpd_select = mean(HNF_Bmpd_select_strength, na.rm = TRUE))

HNF_Bac_A <- Bac_A %>%
  inner_join(Bac_Amntd, by = c("Site" = "Site")) %>%
  inner_join(Bac_Ampd, by = c("Site" = "Site")) %>%
  #inner_join(Bac_Bmntd_selec, by = c("Site" = "Var2")) %>%
  #inner_join(Bac_Bmpd_selec, by = c("Site" = "Var2")) %>%
  
  inner_join(HNF_A, by = c("Site" = "Site")) %>%
  inner_join(HNF_Amntd, by = c("Site" = "Site")) %>%
  inner_join(HNF_Ampd, by = c("Site" = "Site")) %>%
  #inner_join(HNF_Bmntd_selec, by = c("Site" = "Var2")) %>%
  #inner_join(HNF_Bmpd_selec, by = c("Site" = "Var2")) %>%
  
  inner_join(Vars, by = c("Site" = "SampleID")) %>%
  filter(!is.na(NF_Biom)) %>%
  mutate(ln.Bac_q0 = log(Bac_q0),
         ln.HNF_q0 = log(HNF_q0),
         ln.Bac_q1 = log(Bac_q1),
         ln.HNF_q1 = log(HNF_q1),
         ln.Bac_q2 = log(Bac_q2),
         ln.HNF_q2 = log(HNF_q2),
         ln.Bac_Biom = log(Bac_Biom),
         ln.HNF_Biom = log(HNF_Biom),
         ln.Temp = log(Temp),
         ln.Sal = log(Sal),
         ln.PAR = log(PAR),
         ln.NO2 = log(NO2 + 0.0001),
         ln.NO3 = log(NO3 + 0.0001),
         ln.DIN = log(DIN + 0.0001),
         ln.PO3 = log(PO3 + 0.0001), 
         ln.Chla = log(Chla + 0.00001))
HNF_Bac_A <- as.data.frame(HNF_Bac_A)
head(HNF_Bac_A)
###############################################################################################
##### Preping data ############################################################################
###############################################################################################

###############################################################################################
##### Envi exploratory factor analyses and  pair-wise correlations ############################
###############################################################################################

##### exploratory factor analyses on environmental data ##########
### plotting
p_Envi_pairs <- HNF_Bac_A %>%
  ggpairs(columns = c("Temp", "Sal", "PAR", "NO2", "NO3", "DIN", "PO3", "Chla"),
          columnLabels = c("Temperature", "Salinity", "PAR", "Nitrite", "Nitrate", "TN", "TP", "Chla"),
          upper = list(continuous = cor_fun),
          lower = list(continuous = fit_fun)) +
  theme(strip.text.x = element_text(color = "black", size = 14),
        strip.text.y = element_text(angle = 45, color = "black", size = 14))
p_Envi_pairs
# ggsave(p_Envi_pairs, file = "D:/Research/PdPy_Div_Results/p_Envi_pairs.jpeg", 
#        dpi = 600, width = 34, height = 28, units = "cm")

Envi <- HNF_Bac_A[, c("Temp", "Sal", "PAR", "DIN", "PO3", "Chla")] #, "NO2", "NO3"

fa <- fa.parallel(Envi, fm = "mle", fa = 'fa')

fa1 <- fa(Envi, nfactors = 1, rotate = "varimax", fm = "mle", n.iter = 1000)
print(fa1)
fa2 <- fa(Envi, nfactors = 2, rotate = "varimax", fm = "mle", n.iter = 1000)
print(fa2)
fa3 <- fa(Envi, nfactors = 3, rotate = "varimax", fm = "mle", n.iter = 1000)
print(fa3)
fa4 <- fa(Envi, nfactors = 4, rotate = "varimax", fm = "mle", n.iter = 1000)
print(fa4)
fa5 <- fa(Envi, nfactors = 5, rotate = "varimax", fm = "mle", n.iter = 1000)
print(fa5)
fa.diagram(fa3)
##### exploratory factor analyses on environmental data ##########

##### Pair-wise plot of bio-variables ##########
p_Adiv_AAssemb_pairs <- HNF_Bac_A %>%
  ggpairs(columns = c("ln.Bac_q0", "ln.HNF_q0", "ln.Bac_q1", "ln.HNF_q1", "ln.Bac_q2", "ln.HNF_q2",
                      "Bac_Amntd_select", "Bac_Ampd_select", "HNF_Amntd_select", "HNF_Ampd_select"),
          columnLabels = c("Bacteria\nspecies\nrichness", "HNF\nspecies\nrichness",
                           "Bacteria\nShannon\ndiversity", "HNF\nShannon\ndiversity", 
                           "Bacteria\nSimpson\ndiversity", "HNF\nSimpson\ndiversity",
                           "Bacteria\n\U03B1NTI", "Bacteria\n\U03B1MPTI", "HNF\n\U03B1NTI", "HNF\n\U03B1MPTI"),
          #mapping = ggplot2::aes(colour = Cruise),
          upper = list(continuous = cor_fun),
          lower = list(continuous = fit_fun)) +
  theme(strip.text.x = element_text(color = "black", size = 14),
        strip.text.y = element_text(angle = 45, color = "black", size = 14))
p_Adiv_AAssemb_pairs
ggsave(p_Adiv_AAssemb_pairs, file = "D:/Research/PdPy_Div_Results/Figs/PR2_3/p_Adiv_AAssemb_pairs.png",
       dpi = 600, width = 34, height = 28, units = "cm")

# zooming 
p_Adiv_AAssemb_zoompairs <- HNF_Bac_A %>%
  ggpairs(columns = c("ln.Bac_q1", "ln.HNF_q1", "Bac_Ampd_select", "HNF_Ampd_select"),
          columnLabels = c("Bacteria\nspecies\nrichness", "HNF\nSimpson\ndiversity",
                           "Bacteria\n\U03B1MPTI", "HNF\n\U03B1MPTI"),
          #mapping = ggplot2::aes(colour = Cruise),
          upper = list(continuous = cor_fun),
          lower = list(continuous = fit_fun)) +
  theme(strip.text.x = element_text(color = "black", size = 14),
        strip.text.y = element_text(angle = 45, color = "black", size = 14))
p_Adiv_AAssemb_zoompairs

ggsave(p_Adiv_AAssemb_zoompairs, file = "D:/Research/PdPy_Div_Results/Figs/PR2_3/p_Adiv_AAssemb_zoompairs.png",
       dpi = 600, width = 34, height = 28, units = "cm")

##### Pair-wise plot of bio-variables ##########
###############################################################################################
##### Envi exploratory factor analyses and  pair-wise correlations ############################
###############################################################################################

###############################################################################################
##### Assemnbly processes force visualization  ################################################
###############################################################################################
force.labs <- c("Bacteria \U03B1MPTI", "HNF \U03B1MPTI")
names(force.labs) <- c("Bac_Ampd_select", "HNF_Ampd_select")

p_force <- HNF_Bac_A %>% 
  select(Bac_Ampd_select, HNF_Ampd_select) %>%
  gather(key = "force", value = "strength") %>%
  ggplot(aes(strength)) + 
    geom_histogram(aes(y =..density..),
                   col = "red", 
                   fill = "green", 
                   alpha = 0.2) + 
    geom_density(col = 2) +
    facet_wrap(~ force, ncol = 2, dir = "v", labeller = labeller(force = force.labs)) +  
  theme(
    strip.text.x = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16)
  )
p_force
ggsave(p_force, file = "D:/Research/PdPy_Div_Results/Figs/PR2_new/p_force.png",
       dpi = 600, width = 34, height = 28, units = "cm")

###############################################################################################
##### Assemnbly processes force visualization  ################################################
###############################################################################################

###############################################################################################
##### Simple HNF and Bac alppha diversity relationship  #######################################
###############################################################################################
##### without selection #####
### linear model
# Statistical check on non-linear relationship between Bac and HNF alpha diversity
gam0 <- gam(ln.HNF_q1 ~ ln.Bac_q1, data = HNF_Bac_A)
gam1 <- gam(ln.HNF_q1 ~ s(ln.Bac_q1), data = HNF_Bac_A)
anova(gam0, gam1, test = "F")

gam_select <- gam(ln.HNF_q1 ~ ln.Bac_q1*Bac_Ampd_select + ln.Bac_q1*HNF_Ampd_select, data = HNF_Bac_A)
summary(gam_select)
lme_select <- lme(ln.HNF_q1 ~ ln.Bac_q1*Bac_Ampd_select + ln.Bac_q1*HNF_Amntd_select, random = ( ~1 | Cruise), data = HNF_Bac_A)
summary(lme_select)
### plotting (without selection)
p_HNFq1_Bacq1 <- HNF_Bac_A %>% 
  select(ln.Bac_q1, ln.HNF_q1, Bac_Ampd_select, HNF_Ampd_select, Cruise) %>%
  ggplot(aes(x = ln.Bac_q1, y = ln.HNF_q1)) + 
    geom_point(size = 3, aes(color = Cruise)) + 
    geom_smooth(formula = y ~ x, method = "lm", se = TRUE) + 
    # geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "dotted") + 
    labs(x = expression("Log[ Bacteria Shannon diversity (Hill number = 1) ]"),
         y = expression("Log[ HNF Shannon diversity (Hill number = 1) ]")) + 
    theme(
      strip.text.x = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 16)
    )
p_HNFq1_Bacq1
ggsave(p_HNFq1_Bacq1, file = "D:/Research/PdPy_Div_Results/Figs/PR2_3/p_ADiv_HNFq1_Bacq1.png",
       dpi = 600, width = 34, height = 28, units = "cm")
##### without selection #####

##### add selection and dispersal #####
### linear model
HNFq1_Bacq1_0 <- glm(ln.HNF_q1 ~ ln.Bac_q1, data = HNF_Bac_A)
HNFq1_Bacq1_1 <- lme(ln.HNF_q1 ~ ln.Bac_q1, random = ( ~1 | Cruise), data = HNF_Bac_A, method = "ML")

anova(HNFq1_Bacq1_0, HNFq1_Bacq1_1, test = "F")
AIC(HNFq1_Bacq1_0, HNFq1_Bacq1_1)

summary(HNFq1_Bacq1_1)

# Plotting HNF and Bac alppha diversity relationship with selection as color code
community.labs <- c("Colored by selection on bacteria community", "Colored by selection on HNF community")
names(community.labs) <- c("Bac_Ampd_select", "HNF_Ampd_select")

p_ADiv_Select <- HNF_Bac_A %>% 
  select(ln.Bac_q1, ln.HNF_q1, Bac_Ampd_select, HNF_Ampd_select) %>%
  gather(key = "community", value = "Selection", -c(ln.Bac_q1, ln.HNF_q1)) %>%
  ggplot(aes(x = ln.Bac_q1, y = ln.HNF_q1, color = Selection)) + 
    geom_point(size = 3) + 
    geom_smooth(formula = y ~ x, method = "lm", se = TRUE) + 
    #geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "dotted") + 
    facet_grid(~ community, labeller = labeller(community = community.labs)) +  
    scale_colour_viridis(alpha = 0.7) + 
    labs(x = expression("Log[ Bacteria Shannon diversity (Hill number = 1) ]"),
         y = expression("Log[ HNF Shannon diversity (Hill number = 1) ]"),
         colour = expression(paste("\u03B1", "MPTI"))) + 
    theme(
      strip.text.x = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 16)
    )
p_ADiv_Select
ggsave(p_ADiv_Select, file = "D:/Research/PdPy_Div_Results/Figs/PR2_3/p_ADiv_HNFq1_Bacq1_Ampd.png",
       dpi = 600, width = 34, height = 28, units = "cm")
##### add selection and dispersal #####
###############################################################################################
##### Simple HNF and Bac alppha diversity relationship  #######################################
###############################################################################################

###############################################################################################
##### Testing hypothesis: HNFq0 -> Bac selection -> Bacq0 & Bacq0 -> Bac selection ->  HNFq0 ##
###############################################################################################
##### HNFq1 -> Bac selection -> Bacq1 ##########
### linear model testing
BacS_HNFq1.0 <- lm(Bac_Ampd_select ~ ln.HNF_q1, data = HNF_Bac_A)
BacS_HNFq1.Cr <- lme(Bac_Ampd_select ~ ln.HNF_q1, random = ~ 1 | Cruise, data = HNF_Bac_A, method = "ML")
BacS_HNFq1.Season <- lme(Bac_Amntd_select ~ ln.HNF_q1, random = ~ 1 | Season, data = HNF_Bac_A, method = "ML")
AIC(BacS_HNFq1.1, BacS_HNFq1.Cr, BacS_HNFq1.Season)
summary(BacS_HNFq1.Cr)
summary(gam(Bac_Amntd_select ~ s(ln.HNF_q2, bs = "ts"), data = HNF_Bac_A))

Bacq1_BacS.0 <- lm(ln.Bac_q1 ~ Bac_Ampd_select, data = HNF_Bac_A)
Bacq1_BacS.Cr <- lme(ln.Bac_q1 ~ Bac_Ampd_select, random = ~ 1 | Cruise, data = HNF_Bac_A, method = "ML")
Bacq1_BacS.Season <- lme(ln.Bac_q1 ~ Bac_Ampd_select, random = ~ 1 | Season, data = HNF_Bac_A, method = "ML")
AIC(Bacq1_BacS.0, Bacq1_BacS.Cr, Bacq1_BacS.Season)
summary(Bacq1_BacS.Cr)
summary(gam(ln.Bac_q1 ~ s(Bac_Ampd_select, bs = "ts"), data = HNF_Bac_A))

### plotting
p_BacSelect_Bacq1 <- HNF_Bac_A %>% 
  select(ln.Bac_q1, ln.HNF_q1, Bac_Ampd_select, HNF_Ampd_select, Cruise) %>%
  ggplot(aes(x = Bac_Ampd_select, y = ln.Bac_q1)) + 
    geom_point(aes(color = Cruise), size = 3) + 
    geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "solid") + 
    #geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "dotted") + 
    scale_colour_viridis(alpha = 0.7, discrete=TRUE) + 
    labs(x = expression("Deterministic assembly processes (\U03B1MPTI) of Bacteria community "),
         y = expression(atop("Log[ Bacterial Shannon diversity", "(Hill number = 1) ]"))) + 
    theme(axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 0, b = 0, l = 20)),
          axis.title.x = element_text(size = 20),
          axis.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 16))

p_HNFq1_BacSelect <- HNF_Bac_A %>% 
  select(ln.Bac_q1, ln.HNF_q1, Bac_Ampd_select, HNF_Ampd_select, Cruise) %>%
  ggplot(aes(x = ln.HNF_q1, y = Bac_Ampd_select)) + 
    geom_point(aes(color = Cruise), size = 3) + 
    geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "solid") + 
    #geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "solid") + 
    scale_colour_viridis(alpha = 0.7, discrete=TRUE) + 
    labs(x = expression("Log[ HNF Shannon diversity (Hill number = 1) ]"),
         y = expression(atop("Deterministic assembly processes (\U03B1MPTI)", "of Bacteria community"))) + 
    theme(axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 0, b = 0, l = 20)),
          axis.title.x = element_text(size = 20),
          axis.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 16))

legend <- get_legend(
  p_HNFq1_BacSelect + theme(legend.box.margin = margin(0, 0, 0, 12))
)
p_HNFq1_BacSelect_Bacq1 <- plot_grid(
  plot_grid(p_BacSelect_Bacq1 + theme(legend.position = "none"),
            p_HNFq1_BacSelect + theme(legend.position = "none"),
            ncol = 1, labels = "AUTO", hjust = 0),
  legend, rel_widths = c(3, .4))
p_HNFq1_BacSelect_Bacq1
ggsave(p_HNFq1_BacSelect_Bacq1, file = "D:/Research/PdPy_Div_Results/Figs/PR2_3/p_HNFq1_BacAmpd_Bacq1.png",
       dpi = 600, width = 36, height = 42, units = "cm")
##### HNFq1 -> Bac selection -> Bacq1 ##########

##### Bacq1 -> HNF selection -> HNFq1 ##########
### linear model testing
HNFS_Bacq1.0 <- lm(HNF_Ampd_select ~ ln.Bac_q1, data = HNF_Bac_A)
HNFS_Bacq1.Cr <- lme(HNF_Ampd_select ~ ln.Bac_q1, random = ~ 1 | Cruise, data = HNF_Bac_A, method = "ML")
HNFS_Bacq1.Season <- lme(HNF_Ampd_select ~ ln.Bac_q1, random = ~ 1 | Season, data = HNF_Bac_A, method = "ML")
AIC(HNFS_Bacq1.0, HNFS_Bacq1.Cr, HNFS_Bacq1.Season)
summary(HNFS_Bacq1.Cr)
summary(gam(HNF_Ampd_select ~ s(ln.Bac_q1, bs = "ts"), data = HNF_Bac_A))

HNFq1_HNFS.0 <- lm(ln.HNF_q1 ~ HNF_Amntd_select, data = HNF_Bac_A)
HNFq1_HNFS.Cr <- lme(ln.HNF_q1 ~ HNF_Amntd_select, random = ~ 1 | Cruise, data = HNF_Bac_A, method = "ML")
HNFq1_HNFS.Season <- lme(ln.HNF_q1 ~ HNF_Amntd_select, random = ~ 1 | Season, data = HNF_Bac_A, method = "ML")
AIC(HNFq1_HNFS.0, HNFq1_HNFS.Cr, HNFq1_HNFS.Season)
summary(HNFq1_HNFS.Cr)
summary(gam(ln.HNF_q1 ~ s(HNF_Ampd_select, bs = "ts"), data = HNF_Bac_A))

### plotting
p_HNFSelect_HNFq1 <- HNF_Bac_A %>% 
  select(ln.Bac_q1, ln.HNF_q1, Bac_Ampd_select, HNF_Ampd_select, Cruise) %>%
  ggplot(aes(x = HNF_Ampd_select, y = ln.HNF_q1)) + 
    geom_point(aes(color = Cruise), size = 3) + 
    geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "solid") + 
    #geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "dotted") + 
    scale_colour_viridis(alpha = 0.7, discrete=TRUE) + 
    labs(x = expression("Deterministic assembly processes (\U03B1MPTI) of HNF community "),
         y = expression(atop("Log[ HNF Shannon diversity", "(Hill number = 1) ]"))) + 
    theme(axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 0, b = 0, l = 20)),
          axis.title.x = element_text(size = 20),
          axis.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 16))

p_Bacq1_HNFSelect <- HNF_Bac_A %>% 
  select(ln.Bac_q1, ln.HNF_q1, Bac_Ampd_select, HNF_Ampd_select, Cruise) %>%
  ggplot(aes(x = ln.Bac_q1, y = HNF_Ampd_select)) + 
  geom_point(aes(color = Cruise), size = 3) + 
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "dashed") + 
  #geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "dotted") + 
  scale_colour_viridis(alpha = 0.7, discrete=TRUE) + 
  labs(x = expression("Log[ Bacterial Shannon diversity (Hill number = 1) ]"),
       y = expression(atop("Deterministic assembly processes (\U03B1MPTI)", "of HNF community "))) + 
  theme(axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 0, b = 0, l = 20)),
        axis.title.x = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))

legend <- get_legend(
  p_Bacq1_HNFSelect + theme(legend.box.margin = margin(0, 0, 0, 12))
)
p_Bacq1_HNFSelect_HNFq1 <- plot_grid(
  plot_grid(p_HNFSelect_HNFq1 + theme(legend.position = "none"),
            p_Bacq1_HNFSelect + theme(legend.position = "none"),
            ncol = 1, labels = "AUTO", hjust = 0),
  legend, rel_widths = c(3, .4))
p_Bacq1_HNFSelect_HNFq1
ggsave(p_Bacq1_HNFSelect_HNFq1, file = "D:/Research/PdPy_Div_Results/Figs/PR2_3/p_Bacq1_HNFAmpd_HNFq1.png",
       dpi = 600, width = 36, height = 42, units = "cm")
##### Bacq0 -> HNF selection ->  HNFq0 ##########

##### Path model analysis : Bac_q0 vs HNF_q0 ##########
##### Step 1: no random effects #####
Bacq1_HNFq1_mod1.0 <- '
  # regressions
    ln.Bac_q1 ~ Bac_Ampd_select + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_q1 ~ HNF_Ampd_select + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    Bac_Ampd_select ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    HNF_Ampd_select ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
  
  # correlated error
    ln.Bac_q1 ~~ ln.HNF_q1
'

Bacq1_HNFq1_mod1.1 <- '
  # regressions
    ln.Bac_q1 ~ Bac_Ampd_select + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_q1 ~ HNF_Ampd_select + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    Bac_Ampd_select ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    HNF_Ampd_select ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.Bac_Biom ~ ln.Bac_q1 + ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_Biom ~ ln.HNF_q1 + ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
  
  # correlated error
    ln.Bac_q1 ~~ ln.HNF_q1
    ln.Bac_Biom ~~ ln.HNF_Biom
'

Bacq1_HNFq1_mod1.2 <- '
  # regressions
    ln.Bac_q1 ~ Bac_Ampd_select + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_q1 ~ HNF_Ampd_select + ln.Bac_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    Bac_Ampd_select ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    HNF_Ampd_select ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.Bac_Biom ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_Biom ~ ln.HNF_q1 + ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
  
  # correlated error  
    ln.Bac_q1 ~~ ln.HNF_q1
    ln.Bac_Biom ~~ ln.HNF_Biom
'

Bacq1_HNFq1_mod1.3 <- '
  # regressions
    ln.Bac_q1 ~ Bac_Ampd_select + ln.HNF_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_q1 ~ HNF_Ampd_select + ln.Bac_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    Bac_Ampd_select ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    HNF_Ampd_select ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.Bac_Biom ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_Biom ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
  
  # correlated error
    ln.Bac_q1 ~~ ln.HNF_q1
    ln.Bac_Biom ~~ ln.HNF_Biom
'

Bacq1_HNFq1_mod1.4 <- '
  # regressions
    ln.Bac_q1 ~ Bac_Ampd_select + ln.HNF_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_q1 ~ HNF_Ampd_select + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    Bac_Ampd_select ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    HNF_Ampd_select ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.Bac_Biom ~ ln.Bac_q1 + ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_Biom ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
  
  # correlated error
    ln.Bac_q1 ~~ ln.HNF_q1
    ln.Bac_Biom ~~ ln.HNF_Biom
'

Bacq1_HNFq1_mod1.5 <- '
  # regressions
    ln.Bac_q1 ~ Bac_Ampd_select + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_q1 ~ HNF_Ampd_select + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    Bac_Ampd_select ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    HNF_Ampd_select ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.Bac_Biom ~ ln.Bac_q1 + ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_Biom ~ ln.HNF_q1 + ln.Bac_q1 + ln.Bac_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    
  # correlated error
    ln.Bac_q1 ~~ ln.HNF_q1
'

Bacq1_HNFq1_mod1.6 <- '
  # regressions
    ln.Bac_q1 ~ Bac_Ampd_select + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_q1 ~ HNF_Ampd_select + ln.Bac_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    Bac_Ampd_select ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    HNF_Ampd_select ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.Bac_Biom ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_Biom ~ ln.HNF_q1 + ln.Bac_q1 + ln.Bac_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    
  # correlated error
    ln.Bac_q1 ~~ ln.HNF_q1
'

Bacq1_HNFq1_mod1.7 <- '
  # regressions
    ln.Bac_q1 ~ Bac_Ampd_select + ln.HNF_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_q1 ~ HNF_Ampd_select + ln.Bac_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    Bac_Ampd_select ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    HNF_Ampd_select ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.Bac_Biom ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_Biom ~ ln.HNF_q1 + ln.Bac_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    
  # correlated error
    ln.Bac_q1 ~~ ln.HNF_q1
'

Bacq1_HNFq1_mod1.8 <- '
  # regressions
    ln.Bac_q1 ~ Bac_Ampd_select + ln.HNF_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_q1 ~ HNF_Ampd_select + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    Bac_Ampd_select ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    HNF_Ampd_select ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.Bac_Biom ~ ln.Bac_q1 + ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_Biom ~ ln.HNF_q1 + ln.Bac_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    
  # correlated error
    ln.Bac_q1 ~~ ln.HNF_q1
'

Bacq1_HNFq1_mod1.9 <- '
  # regressions
    ln.Bac_q1 ~ Bac_Ampd_select + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_q1 ~ HNF_Ampd_select + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    Bac_Ampd_select ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    HNF_Ampd_select ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.Bac_Biom ~ ln.Bac_q1 + ln.HNF_q1 + ln.HNF_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_Biom ~ ln.HNF_q1 + ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    
  # correlated error
    ln.Bac_q1 ~~ ln.HNF_q1
'

Bacq1_HNFq1_mod1.10 <- '
  # regressions
    ln.Bac_q1 ~ Bac_Ampd_select + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_q1 ~ HNF_Ampd_select + ln.Bac_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    Bac_Ampd_select ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    HNF_Ampd_select ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.Bac_Biom ~ ln.Bac_q1 + ln.HNF_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_Biom ~ ln.HNF_q1 + ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    
  # correlated error
    ln.Bac_q1 ~~ ln.HNF_q1
'

Bacq1_HNFq1_mod1.11 <- '
  # regressions
    ln.Bac_q1 ~ Bac_Ampd_select + ln.HNF_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_q1 ~ HNF_Ampd_select + ln.Bac_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    Bac_Ampd_select ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    HNF_Ampd_select ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.Bac_Biom ~ ln.Bac_q1 + ln.HNF_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_Biom ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    
  # correlated error
    ln.Bac_q1 ~~ ln.HNF_q1
'

Bacq1_HNFq1_mod1.12 <- '
  # regressions
    ln.Bac_q1 ~ Bac_Ampd_select + ln.HNF_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_q1 ~ HNF_Ampd_select + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    Bac_Ampd_select ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    HNF_Ampd_select ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.Bac_Biom ~ ln.Bac_q1 + ln.HNF_q1 + ln.HNF_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_Biom ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    
  # correlated error
    ln.Bac_q1 ~~ ln.HNF_q1
'
Bacq1_HNFq1_lavaan1.0 <- sem(Bacq1_HNFq1_mod1.0, data = HNF_Bac_A)
Bacq1_HNFq1_lavaan1.1 <- sem(Bacq1_HNFq1_mod1.1, data = HNF_Bac_A)
Bacq1_HNFq1_lavaan1.2 <- sem(Bacq1_HNFq1_mod1.2, data = HNF_Bac_A)
Bacq1_HNFq1_lavaan1.3 <- sem(Bacq1_HNFq1_mod1.3, data = HNF_Bac_A)
Bacq1_HNFq1_lavaan1.4 <- sem(Bacq1_HNFq1_mod1.4, data = HNF_Bac_A)
Bacq1_HNFq1_lavaan1.5 <- sem(Bacq1_HNFq1_mod1.5, data = HNF_Bac_A)
Bacq1_HNFq1_lavaan1.6 <- sem(Bacq1_HNFq1_mod1.6, data = HNF_Bac_A)
Bacq1_HNFq1_lavaan1.7 <- sem(Bacq1_HNFq1_mod1.7, data = HNF_Bac_A)
Bacq1_HNFq1_lavaan1.8 <- sem(Bacq1_HNFq1_mod1.8, data = HNF_Bac_A)
Bacq1_HNFq1_lavaan1.9 <- sem(Bacq1_HNFq1_mod1.9, data = HNF_Bac_A)
Bacq1_HNFq1_lavaan1.10 <- sem(Bacq1_HNFq1_mod1.10, data = HNF_Bac_A)
Bacq1_HNFq1_lavaan1.11 <- sem(Bacq1_HNFq1_mod1.11, data = HNF_Bac_A)
Bacq1_HNFq1_lavaan1.12 <- sem(Bacq1_HNFq1_mod1.12, data = HNF_Bac_A)

AICstep1 <- AIC(Bacq1_HNFq1_lavaan1.0, Bacq1_HNFq1_lavaan1.1, Bacq1_HNFq1_lavaan1.2, Bacq1_HNFq1_lavaan1.3,
                Bacq1_HNFq1_lavaan1.4, Bacq1_HNFq1_lavaan1.5, Bacq1_HNFq1_lavaan1.6, Bacq1_HNFq1_lavaan1.7,
                Bacq1_HNFq1_lavaan1.8, Bacq1_HNFq1_lavaan1.9, Bacq1_HNFq1_lavaan1.10, Bacq1_HNFq1_lavaan1.11,
                Bacq1_HNFq1_lavaan1.12)
AICstep1 <- AICstep1 %>% cbind(row.names(AICstep1)) %>%
  arrange(AIC)
AICstep1

Bacq1_HNFq1_mod1.3 <- '
  # regressions
    ln.Bac_q1 ~ Bac_Ampd_select + ln.HNF_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_q1 ~ HNF_Ampd_select + ln.Bac_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    Bac_Ampd_select ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    HNF_Ampd_select ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.Bac_Biom ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
    ln.HNF_Biom ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3
  
  # correlated error
    ln.Bac_q1 ~~ ln.HNF_q1
    ln.Bac_Biom ~~ ln.HNF_Biom
'
Bacq1_HNFq1_lavaan1.3 <- sem(Bacq1_HNFq1_mod1.3, data = HNF_Bac_A)#, se = "bootstrap")
summary(Bacq1_HNFq1_lavaan1.3, fit.measures = TRUE)
##### Step 1: no random effects #####
##### Step 2 : include grouping variables (random effects) and environmental variables #####
Bacq1_HNFq1_psem2.0 <- psem(
  lm(ln.Bac_q1 ~ Bac_Ampd_select + ln.HNF_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3,
     data = HNF_Bac_A),
  lm(ln.HNF_q1 ~ HNF_Ampd_select + ln.Bac_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3,
     data = HNF_Bac_A),
  lm(Bac_Ampd_select ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3,
    data = HNF_Bac_A),
  # lm(HNF_Ampd_select ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3,
  #    data = HNF_Bac_A),
  # lm(ln.Bac_Biom ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3,
  #    data = HNF_Bac_A),
  lm(ln.HNF_Biom ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3,
    data = HNF_Bac_A),
  
  ln.Bac_q1 %~~% ln.HNF_q1,
  ln.Bac_Biom %~~% ln.HNF_Biom
)
summary(Bacq1_HNFq1_psem2.0)
Bacq1_HNFq1_psem2.Cr <- psem(
  lme(ln.Bac_q1 ~ Bac_Ampd_select + ln.HNF_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3,
      random = ~ 1 | Cruise, data = HNF_Bac_A, method = "ML"),
  lme(ln.HNF_q1 ~ HNF_Ampd_select + ln.Bac_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3,
      random = ~ 1 | Cruise, data = HNF_Bac_A, method = "ML"),
  lme(Bac_Ampd_select ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3,
      random = ~ 1 | Cruise, data = HNF_Bac_A, method = "ML"),
  # lme(HNF_Ampd_select ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3,
  #     random = ~ 1 | Cruise, data = HNF_Bac_A),
  # lme(ln.Bac_Biom ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3,
  #     random = ~ 1 | Cruise, data = HNF_Bac_A),
  # lme(ln.HNF_Biom ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3,
  #     random = ~ 1 | Cruise, data = HNF_Bac_A, method = "ML"),
  
  ln.Bac_q1 %~~% ln.HNF_q1
  # ln.Bac_Biom %~~% ln.HNF_Biom
)
summary(Bacq1_HNFq1_psem2.Cr, fit.measures = TRUE)
summary(lme(HNF_Ampd_select ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3,
    random = ~ 1 | Cruise, data = HNF_Bac_A))
##### Path model analysis : Bac_q0 vs HNF_q0 ##########
###############################################################################################
##### Testing hypothesis: HNFq0 -> Bac selection -> Bacq0 & Bacq0 -> Bac selection ->  HNFq0 ##
###############################################################################################
