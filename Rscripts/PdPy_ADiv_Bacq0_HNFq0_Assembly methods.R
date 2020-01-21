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
Bac_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/16s_seqXst.csv",
                                       sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
Bac_ra_comm <- Bac_comm / rowSums(Bac_comm)
Bac_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/treeNJ_16s.tree")

NF_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/prot_seqXst_PR2.csv", 
                                      sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
NF_ra_comm <- NF_comm / rowSums(NF_comm)
NF_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/prot_treeNJ_PR2.tree")

HNF_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/HNF_seqXst_PR2.csv",
                                       sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
HNF_ra_comm <- HNF_comm / rowSums(HNF_comm)
HNF_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/HNF_treeNJ_PR2.tree")

Vars <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_Vars.csv", sep = ",", 
                   header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
###############################################################################################
##### Loading data ############################################################################
###############################################################################################

###############################################################################################
##### Loading nulls ###########################################################################
###############################################################################################
Bac_Amntd_null <- read.table(file = "D:/Research/PdPy_Div_Results/Bac_Amntd_null.csv", sep = ",", 
                             header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
Bac_Ampd_null <- read.table(file = "D:/Research/PdPy_Div_Results/Bac_Ampd_null.csv", sep = ",", 
                            header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
Bac_Bmntd_null <- read.table(file = "D:/Research/PdPy_Div_Results/Bac_Bmntd_null.csv", sep = ",", 
                             header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
Bac_Bmpd_null <- read.table(file = "D:/Research/PdPy_Div_Results/Bac_Bmpd_null.csv", sep = ",", 
                            header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
Bac_Chao_null <- read.table(file = "D:/Research/PdPy_Div_Results/Bac_Chao_null.csv", sep = ",", 
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
Bac_BDiv_Chao <- Bac_Chao_null %>% 
  select(c(obs, Var1, Var2)) %>%
  mutate(Chao_null_mean = apply(Bac_Chao_null[, !names(Bac_Chao_null) %in% c("obs", "Var1", "Var2")], 1, mean),
         Chao_null_sd = apply(Bac_Chao_null[, !names(Bac_Chao_null) %in% c("obs", "Var1", "Var2")], 1, sd),
         Bac_disp_strength = (obs - Chao_null_mean) / Chao_null_sd,
         Bac_disp_p = pnorm(Bac_disp_strength, 0, 1))

HNF_Amntd_null <- read.table(file = "D:/Research/PdPy_Div_Results/HNF_Amntd_null_PR2.csv", sep = ",", 
                             header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
HNF_Ampd_null <- read.table(file = "D:/Research/PdPy_Div_Results/HNF_Ampd_null_PR2.csv", sep = ",", 
                            header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
HNF_Bmntd_null <- read.table(file = "D:/Research/PdPy_Div_Results/HNF_Bmntd_null_PR2.csv", sep = ",", 
                             header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
HNF_Bmpd_null <- read.table(file = "D:/Research/PdPy_Div_Results/HNF_Bmpd_null_PR2.csv", sep = ",", 
                            header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
HNF_Chao_null <- read.table(file = "D:/Research/PdPy_Div_Results/HNF_Chao_null_PR2.csv", sep = ",", 
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
HNF_BDiv_Chao <- HNF_Chao_null %>% 
  select(c(obs, Var1, Var2)) %>%
  mutate(Chao_null_mean = apply(HNF_Chao_null[, !names(HNF_Chao_null) %in% c("obs", "Var1", "Var2")], 1, mean),
         Chao_null_sd = apply(HNF_Chao_null[, !names(HNF_Chao_null) %in% c("obs", "Var1", "Var2")], 1, sd),
         HNF_disp_strength = (obs - Chao_null_mean) / Chao_null_sd,
         HNF_disp_p = pnorm(HNF_disp_strength, 0, 1))
###############################################################################################
##### Loading nulls ###########################################################################
###############################################################################################

###############################################################################################
##### Envi exploratory factor analyses and  pair-wise correlations ############################
###############################################################################################
##### Preping data ##########
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

Bac_disp <- Bac_BDiv_Chao %>%
  mutate(Cr_V1 = substr(Var1, start = 1, stop = 9),
         Cr_V2 = substr(Var2, start = 1, stop = 9)) %>%
  filter(Cr_V1 == Cr_V2) %>%
  #filter(Bac_select_p < 0.05) %>%
  group_by(Var2) %>%
  summarize(Bac_disp = mean(Bac_disp_strength, na.rm = TRUE))

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

HNF_disp <- HNF_BDiv_Chao %>%
  mutate(Cr_V1 = substr(Var1, start = 1, stop = 9),
         Cr_V2 = substr(Var2, start = 1, stop = 9)) %>%
  filter(Cr_V1 == Cr_V2) %>%
  #filter(Bac_select_p < 0.05) %>%
  group_by(Var2) %>%
  summarize(HNF_disp = mean(HNF_disp_strength, na.rm = TRUE))

HNF_Bac_A <- Bac_A %>%
  inner_join(Bac_Amntd, by = c("Site" = "Site")) %>%
  inner_join(Bac_Ampd, by = c("Site" = "Site")) %>%
  inner_join(Bac_Bmntd_selec, by = c("Site" = "Var2")) %>%
  inner_join(Bac_Bmpd_selec, by = c("Site" = "Var2")) %>%
  inner_join(Bac_disp, by = c("Site" = "Var2")) %>%
  
  inner_join(HNF_A, by = c("Site" = "Site")) %>%
  inner_join(HNF_Amntd, by = c("Site" = "Site")) %>%
  inner_join(HNF_Ampd, by = c("Site" = "Site")) %>%
  inner_join(HNF_Bmntd_selec, by = c("Site" = "Var2")) %>%
  inner_join(HNF_Bmpd_selec, by = c("Site" = "Var2")) %>%
  inner_join(HNF_disp, by = c("Site" = "Var2")) %>%
  
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
##### Preping data ##########