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
  file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS/sECS_Bac_seqXst_PR2_new.csv", 
  sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
Bac_ra_comm <- Bac_comm / rowSums(Bac_comm)
Bac_phylo<- read.tree(
  file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS/sECS_Bac_treeNJ_PR2_new.tree")

HNF_comm <- as.data.frame(t(read.table(
  file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS/sECS_HNF_seqXst_PR2_new.csv", 
  sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
HNF_ra_comm <- HNF_comm / rowSums(HNF_comm)
HNF_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS/sECS_HNF_treeNJ_PR2_new.tree")

NF_comm <- as.data.frame(t(read.table(
  file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS/sECS_prot_seqXst_PR2_new.csv", 
  sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
NF_ra_comm <- NF_comm / rowSums(NF_comm)
NF_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS/sECS_prot_treeNJ_PR2_new.tree")

Vars <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_Vars.csv", sep = ",", 
                   header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
###############################################################################################
##### Loading data ############################################################################
###############################################################################################

###############################################################################################
##### Loading nulls ###########################################################################
###############################################################################################
Bac_Amntd_null <- read.table(file = "D:/Research/PdPy_Div_Results/PR2_new/Bac_Amntd_null.csv", sep = ",", 
                             header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
Bac_Ampd_null <- read.table(file = "D:/Research/PdPy_Div_Results/PR2_new/Bac_Ampd_null.csv", sep = ",", 
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

HNF_Amntd_null <- read.table(file = "D:/Research/PdPy_Div_Results/PR2_new/HNF_Amntd_null_PR2.csv", sep = ",", 
                             header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
HNF_Ampd_null <- read.table(file = "D:/Research/PdPy_Div_Results/PR2_new/HNF_Ampd_null_PR2.csv", sep = ",", 
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
                           "Bacteria\naNTI", "Bacteria\naMPTI", "HNF\naNTI", "HNF\naMPTI"),
          #mapping = ggplot2::aes(colour = Cruise),
          upper = list(continuous = cor_fun),
          lower = list(continuous = fit_fun)) +
  theme(strip.text.x = element_text(color = "black", size = 14),
        strip.text.y = element_text(angle = 45, color = "black", size = 14))
p_Adiv_AAssemb_pairs
ggsave(p_Adiv_AAssemb_pairs, file = "D:/Research/PdPy_Div_Results/Figs/PR2_new/p_Adiv_AAssemb_pairs.png",
       dpi = 600, width = 34, height = 28, units = "cm")

# zooming 
p_Adiv_AAssemb_zoompairs <- HNF_Bac_A %>%
  ggpairs(columns = c("ln.Bac_q0", "ln.HNF_q2", "Bac_Amntd_select", "HNF_Amntd_select"),
          columnLabels = c("Bacteria\nspecies\nrichness", "HNF\nSimpson\ndiversity",
                           "Bacteria\naNTI", "HNF\naNTI"),
          #mapping = ggplot2::aes(colour = Cruise),
          upper = list(continuous = cor_fun),
          lower = list(continuous = fit_fun)) +
  theme(strip.text.x = element_text(color = "black", size = 14),
        strip.text.y = element_text(angle = 45, color = "black", size = 14))
p_Adiv_AAssemb_zoompairs

ggsave(p_Adiv_AAssemb_zoompairs, file = "D:/Research/PdPy_Div_Results/Figs/PR2_new/p_Adiv_AAssemb_zoompairs.png",
       dpi = 600, width = 34, height = 28, units = "cm")

##### Pair-wise plot of bio-variables ##########
###############################################################################################
##### Envi exploratory factor analyses and  pair-wise correlations ############################
###############################################################################################

###############################################################################################
##### Assemnbly processes and dispersal force visualization  ##################################
###############################################################################################

force.labs <- c("Bacteria \U03B2NTI", "HNF \U03B2NTI", "Dispersal force on Bacteria", "Dispersal force on HNF")
names(force.labs) <- c("Bac_Amntd_select", "HNF_Amntd_select", "Bac_disp", "HNF_disp")

p_force <- HNF_Bac_A %>% 
  select(Bac_Amntd_select, HNF_Amntd_select, Bac_disp, HNF_disp) %>%
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
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16)
  )
p_force
ggsave(p_force, file = "D:/Research/PdPy_Div_Results/Figs/PR2_new/p_force.png",
       dpi = 600, width = 34, height = 28, units = "cm")

###############################################################################################
##### Assemnbly processes and dispersal force visualization  ##################################
###############################################################################################

###############################################################################################
##### Simple HNF and Bac alppha diversity relationship  #######################################
###############################################################################################
##### without selection #####
### linear model
# Statistical check on non-linear relationship between Bac and HNF alpha diversity
gam0 <- gam(ln.HNF_q2 ~ ln.Bac_q0, data = HNF_Bac_A)
gam1 <- gam(ln.HNF_q2 ~ s(ln.Bac_q0), data = HNF_Bac_A)
anova(gam0, gam1, test = "F")

gam_select <- gam(ln.HNF_q2 ~ ln.Bac_q0*Bac_Amntd_select + ln.Bac_q0*HNF_Amntd_select + Bac_disp + HNF_disp, data = HNF_Bac_A)
summary(gam_select)
lme_select <- lme(ln.HNF_q2 ~ ln.Bac_q0*Bac_Amntd_select + ln.Bac_q0*HNF_Amntd_select + Bac_disp + HNF_disp, 
                  random = ( ~1 | Cruise), data = HNF_Bac_A)
summary(lme_select)
### plotting (without selection)
p_HNFq2_Bacq0 <- HNF_Bac_A %>% 
  select(ln.Bac_q0, ln.HNF_q2, Bac_Amntd_select, HNF_Amntd_select) %>%
  ggplot(aes(x = ln.Bac_q0, y = ln.HNF_q2)) + 
  geom_point(size = 3) + 
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE) + 
  geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "dotted") + 
  scale_colour_viridis(alpha = 0.7) + 
  labs(x = expression("Log[ Bacteria species richness (Hill number = 0) ]"),
       y = expression("Log[ HNF Simpson diversity (Hill number = 2) ]")) + 
  theme(
    strip.text.x = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16)
  )
p_HNFq2_Bacq0
ggsave(p_HNFq2_Bacq0, file = "D:/Research/PdPy_Div_Results/Figs/PR2_new/p_ADiv_HNFq2_Bacq0.png",
       dpi = 600, width = 34, height = 28, units = "cm")
##### without selection #####

##### add selection and dispersal #####
### linear model
HNFq2_Bacq0_0 <- glm(ln.HNF_q2 ~ ln.Bac_q0, data = HNF_Bac_A)
HNFq2_Bacq0_1 <- glm(ln.HNF_q2 ~ ln.Bac_q0 + HNF_Amntd_select, data = HNF_Bac_A)
HNFq2_Bacq0_2 <- glm(ln.HNF_q2 ~ ln.Bac_q0 * HNF_Amntd_select + Bac_disp + HNF_disp, data = HNF_Bac_A)
anova(HNFq2_Bacq0_0, HNFq2_Bacq0_1, HNFq2_Bacq0_2, test = "F")
summary(HNFq2_Bacq0_0)

# Plotting HNF and Bac alppha diversity relationship with selection as color code
community.labs <- c("Colored by selection on bacteria community", "Colored by selection on HNF community")
names(community.labs) <- c("Bac_Amntd_select", "HNF_Amntd_select")

p_ADiv_Select <- HNF_Bac_A %>% 
  select(ln.Bac_q0, ln.HNF_q2, Bac_Amntd_select, HNF_Amntd_select) %>%
  gather(key = "community", value = "Selection", -c(ln.Bac_q0, ln.HNF_q2)) %>%
  ggplot(aes(x = ln.Bac_q0, y = ln.HNF_q2, color = Selection)) + 
    geom_point(size = 3) + 
    geom_smooth(formula = y ~ x, method = "lm", se = TRUE) + 
    geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "dotted") + 
    facet_grid(~ community, labeller = labeller(community = community.labs)) +  
    scale_colour_viridis(alpha = 0.7) + 
    labs(x = expression("Log[ Bacteria species richness (Hill number = 0) ]"),
         y = expression("Log[ HNF Simpson diversity (Hill number = 2) ]"),
         colour = expression(paste("\U03B2", "NTI"))) + 
    theme(
      strip.text.x = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 16)
    )
p_ADiv_Select
ggsave(p_ADiv_Select, file = "D:/Research/PdPy_Div_Results/Figs/PR2_new/p_ADiv_HNFq2_Bacq0_Amntd.png",
       dpi = 600, width = 34, height = 28, units = "cm")
##### add selection and dispersal #####
###############################################################################################
##### Simple HNF and Bac alppha diversity relationship  #######################################
###############################################################################################

###############################################################################################
##### Testing hypothesis: HNFq0 -> Bac selection -> Bacq0 & Bacq0 -> Bac selection ->  HNFq0 ##
###############################################################################################
##### HNFq2 -> Bac selection -> Bacq0 ##########
### linear model testing
BacS_HNFq2.0 <- lm(Bac_Amntd_select ~ ln.HNF_q2, data = HNF_Bac_A)
BacS_HNFq2.Cr <- lme(Bac_Amntd_select ~ ln.HNF_q2, random = ~ 1 | Cruise, data = HNF_Bac_A)
BacS_HNFq2.Season <- lme(Bac_Amntd_select ~ ln.HNF_q2, random = ~ 1 | Season, data = HNF_Bac_A)
BacS_HNFq2.St <- lme(Bac_Amntd_select ~ ln.HNF_q2, random = ~ 1 | Station, data = HNF_Bac_A)
AIC(BacS_HNFq2.0, BacS_HNFq2.Cr, BacS_HNFq2.Season, BacS_HNFq2.St)
summary(BacS_HNFq2.0)
summary(gam(Bac_Amntd_select ~ s(ln.HNF_q2, bs = "ts"), data = HNF_Bac_A))

Bacq0_BacS.0 <- lm(ln.Bac_q0 ~ Bac_Amntd_select + Bac_disp, data = HNF_Bac_A)
Bacq0_BacS.Cr <- lme(ln.Bac_q0 ~ Bac_Amntd_select + Bac_disp, random = ~ 1 | Cruise, data = HNF_Bac_A)
Bacq0_BacS.Season <- lme(ln.Bac_q0 ~ Bac_Amntd_select + Bac_disp, random = ~ 1 | Season, data = HNF_Bac_A)
Bacq0_BacS.St <- lme(ln.Bac_q0 ~ Bac_Amntd_select + Bac_disp, random = ~ 1 | Station, data = HNF_Bac_A)
AIC(Bacq0_BacS.0, Bacq0_BacS.Cr, Bacq0_BacS.Season, Bacq0_BacS.St)
summary(Bacq0_BacS.0)
summary(gam(ln.Bac_q0 ~ s(Bac_Amntd_select, bs = "ts"), data = HNF_Bac_A))
anova(Bacq0_BacS.0, gam(ln.Bac_q0 ~ s(Bac_Amntd_select, bs = "ts"), data = HNF_Bac_A))
### plotting
p_HNFq2_BacSelect <- HNF_Bac_A %>% 
  select(ln.Bac_q0, ln.HNF_q2, Bac_Amntd_select, HNF_Amntd_select, Cruise) %>%
  ggplot(aes(x = ln.HNF_q2, y = Bac_Amntd_select)) + 
  geom_point(aes(color = Cruise), size = 3) + 
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "solid") + 
  geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "solid") + 
  scale_colour_viridis(alpha = 0.7, discrete=TRUE) + 
  labs(x = expression("Log[ HNF Simpson diversity (Hill number = 2) ]"),
       y = expression("Deterministic assembly processes ( \U03B2NTI) of Bacteria community ")) + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))

p_BacSelect_Bacq0 <- HNF_Bac_A %>% 
  select(ln.Bac_q0, ln.HNF_q2, Bac_Amntd_select, HNF_Amntd_select, Cruise) %>%
  ggplot(aes(x = Bac_Amntd_select, y = ln.Bac_q0)) + 
  geom_point(aes(color = Cruise), size = 3) + 
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "solid") + 
  geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "dotted") + 
  scale_colour_viridis(alpha = 0.7, discrete=TRUE) + 
  labs(x = expression("Deterministic assembly processes ( \U03B2NTI) of Bacteria community "),
       y = expression("Log[ Bacteria species richness (Hill number = 0) ]")) + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))

legend <- get_legend(
  p_HNFq2_BacSelect + theme(legend.box.margin = margin(0, 0, 0, 12))
)
p_HNFq2_BacSelect_Bacq0 <- plot_grid(
  plot_grid(p_HNFq2_BacSelect + theme(legend.position = "none"),
            p_BacSelect_Bacq0 + theme(legend.position = "none")),
  legend, rel_widths = c(3, .4))
p_HNFq2_BacSelect_Bacq0
ggsave(p_HNFq2_BacSelect_Bacq0, file = "D:/Research/PdPy_Div_Results/Figs/PR2_new/p_HNFq2_BacAmntd_Bacq0.png",
       dpi = 600, width = 56, height = 28, units = "cm")
##### HNFq2 -> Bac selection -> Bacq0 ##########

##### Bacq0 -> HNF selection ->  HNFq2 ##########
### linear model testing
HNFS_Bacq0.0 <- lm(HNF_Amntd_select ~ ln.Bac_q0, data = HNF_Bac_A)
HNFS_Bacq0.Cr <- lme(HNF_Amntd_select ~ ln.Bac_q0, random = ~ 1 | Cruise, data = HNF_Bac_A)
HNFS_Bacq0.Season <- lme(HNF_Amntd_select ~ ln.Bac_q0, random = ~ 1 | Season, data = HNF_Bac_A)
HNFS_Bacq0.St <- lme(HNF_Amntd_select ~ ln.Bac_q0, random = ~ 1 | Station, data = HNF_Bac_A)
AIC(HNFS_Bacq0.0, HNFS_Bacq0.Cr, HNFS_Bacq0.Season, HNFS_Bacq0.St)
summary(HNFS_Bacq0.0)
summary(gam(HNF_Amntd_select ~ s(ln.Bac_q0, bs = "ts"), data = HNF_Bac_A))

HNFq2_HNFS.0 <- lm(ln.HNF_q2 ~ HNF_Amntd_select + HNF_disp, data = HNF_Bac_A)
HNFq2_HNFS.Cr <- lme(ln.HNF_q2 ~ HNF_Amntd_select + HNF_disp, random = ~ 1 | Cruise, data = HNF_Bac_A)
HNFq2_HNFS.Season <- lme(ln.HNF_q2 ~ HNF_Amntd_select + HNF_disp, random = ~ 1 | Season, data = HNF_Bac_A)
HNFq2_HNFS.St <- lme(ln.HNF_q2 ~ HNF_Amntd_select + HNF_disp, random = ~ 1 | Station, data = HNF_Bac_A)
AIC(HNFq2_HNFS.0, HNFq2_HNFS.Cr, HNFq2_HNFS.Season, HNFq2_HNFS.St)
summary(HNFq2_HNFS.0)
summary(gam(ln.HNF_q2 ~ s(HNF_Amntd_select, bs = "ts"), data = HNF_Bac_A))
AIC(HNFq2_HNFS.0, gam(ln.HNF_q2 ~ s(HNF_Amntd_select, bs = "ts"), data = HNF_Bac_A))

### plotting
p_Bacq0_HNFSelect <- HNF_Bac_A %>% 
  select(ln.Bac_q0, ln.HNF_q0, Bac_Amntd_select, HNF_Amntd_select, Cruise) %>%
  ggplot(aes(x = ln.Bac_q0, y = HNF_Amntd_select)) + 
  geom_point(aes(color = Cruise), size = 3) + 
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "dotted") + 
  geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "dotted") + 
  scale_colour_viridis(alpha = 0.7, discrete=TRUE) + 
  labs(x = expression("Log[ Bacteria species richness (Hill number = 0) ]"),
       y = expression("Deterministic assembly processes ( \U03B2NTI) of HNF community ")) + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))

p_HNFSelect_HNFq2 <- HNF_Bac_A %>% 
  select(ln.Bac_q0, ln.HNF_q2, Bac_Amntd_select, HNF_Amntd_select, Cruise) %>%
  ggplot(aes(x = HNF_Amntd_select, y = ln.HNF_q2)) + 
  geom_point(aes(color = Cruise), size = 3) + 
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "dotted") + 
  geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "solid") + 
  scale_colour_viridis(alpha = 0.7, discrete=TRUE) + 
  labs(x = expression("Deterministic assembly processes ( \U03B2NTI) of HNF community "),
       y = expression("Log[ HNF Simpson diversity (Hill number = 2) ]")) + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))

legend <- get_legend(
  p_Bacq0_HNFSelect + theme(legend.box.margin = margin(0, 0, 0, 12))
)
p_Select_Bacq0_HNFq2 <- plot_grid(
  plot_grid(p_Bacq0_HNFSelect + theme(legend.position="none"),
            p_HNFSelect_HNFq2 + theme(legend.position="none")),
  legend, rel_widths = c(3, .4))
p_Select_Bacq0_HNFq2
ggsave(p_Select_Bacq0_HNFq2, file = "D:/Research/PdPy_Div_Results/Figs/PR2_new/p_Bacq0_HNFAmntd_HNFq2.png",
       dpi = 600, width = 56, height = 28, units = "cm")
##### Bacq0 -> HNF selection ->  HNFq0 ##########

##### Path model analysis : Bac_q0 vs HNF_q0 ##########
##### Step 1: no random effects #####
Bacq0_HNFq2_mod1.0 <- '
  # regressions
    ln.Bac_q0 ~ Bac_Amntd_select + Bac_disp 
    ln.HNF_q2 ~ HNF_Amntd_select + HNF_disp
    Bac_Amntd_select ~ ln.HNF_q2
    HNF_Amntd_select ~ ln.Bac_q0
  
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q2
'

Bacq0_HNFq2_mod1.1 <- '
  # regressions
    ln.Bac_q0 ~ Bac_Amntd_select + Bac_disp
    ln.HNF_q2 ~ HNF_Amntd_select + HNF_disp
    Bac_Amntd_select ~ ln.HNF_q2
    HNF_Amntd_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0 + ln.HNF_q2
    ln.HNF_Biom ~ ln.HNF_q2 + ln.Bac_q0
  
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q2
    ln.Bac_Biom ~~ ln.HNF_Biom
'

Bacq0_HNFq2_mod1.2 <- '
  # regressions
    ln.Bac_q0 ~ Bac_Amntd_select + Bac_disp
    ln.HNF_q2 ~ HNF_Amntd_select + HNF_disp + ln.Bac_Biom
    Bac_Amntd_select ~ ln.HNF_q2
    HNF_Amntd_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0
    ln.HNF_Biom ~ ln.HNF_q2 + ln.Bac_q0
  
  # correlated error  
    ln.Bac_q0 ~~ ln.HNF_q2
    ln.Bac_Biom ~~ ln.HNF_Biom
'

Bacq0_HNFq2_mod1.3 <- '
  # regressions
    ln.Bac_q0 ~ Bac_Amntd_select + Bac_disp + ln.HNF_Biom
    ln.HNF_q2 ~ HNF_Amntd_select + HNF_disp + ln.Bac_Biom
    Bac_Amntd_select ~ ln.HNF_q2
    HNF_Amntd_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0
    ln.HNF_Biom ~ ln.HNF_q2
  
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q2
    ln.Bac_Biom ~~ ln.HNF_Biom
'

Bacq0_HNFq2_mod1.4 <- '
  # regressions
    ln.Bac_q0 ~ Bac_Amntd_select + Bac_disp + ln.HNF_Biom
    ln.HNF_q2 ~ HNF_Amntd_select + HNF_disp
    Bac_Amntd_select ~ ln.HNF_q2
    HNF_Amntd_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0 + ln.HNF_q2
    ln.HNF_Biom ~ ln.HNF_q2
  
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q2
    ln.Bac_Biom ~~ ln.HNF_Biom
'

Bacq0_HNFq2_mod1.5 <- '
  # regressions
    ln.Bac_q0 ~ Bac_Amntd_select + Bac_disp
    ln.HNF_q2 ~ HNF_Amntd_select + HNF_disp
    Bac_Amntd_select ~ ln.HNF_q2
    HNF_Amntd_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0 + ln.HNF_q2
    ln.HNF_Biom ~ ln.HNF_q2 + ln.Bac_q0 + ln.Bac_Biom
    
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q2
'

Bacq0_HNFq2_mod1.6 <- '
  # regressions
    ln.Bac_q0 ~ Bac_Amntd_select + Bac_disp
    ln.HNF_q2 ~ HNF_Amntd_select + HNF_disp + ln.Bac_Biom
    Bac_Amntd_select ~ ln.HNF_q2
    HNF_Amntd_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0
    ln.HNF_Biom ~ ln.HNF_q2 + ln.Bac_q0 + ln.Bac_Biom
    
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q2
'

Bacq0_HNFq2_mod1.7 <- '
  # regressions
    ln.Bac_q0 ~ Bac_Amntd_select + Bac_disp + ln.HNF_Biom
    ln.HNF_q2 ~ HNF_Amntd_select + HNF_disp + ln.Bac_Biom
    Bac_Amntd_select ~ ln.HNF_q2
    HNF_Amntd_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0
    ln.HNF_Biom ~ ln.HNF_q2 + ln.Bac_Biom
    
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q2
'

Bacq0_HNFq2_mod1.8 <- '
  # regressions
    ln.Bac_q0 ~ Bac_Amntd_select + Bac_disp + ln.HNF_Biom
    ln.HNF_q2 ~ HNF_Amntd_select + HNF_disp
    Bac_Amntd_select ~ ln.HNF_q2
    HNF_Amntd_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0 + ln.HNF_q2
    ln.HNF_Biom ~ ln.HNF_q2 + ln.Bac_Biom
    
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q2
'

Bacq0_HNFq2_mod1.9 <- '
  # regressions
    ln.Bac_q0 ~ Bac_Amntd_select + Bac_disp
    ln.HNF_q2 ~ HNF_Amntd_select + HNF_disp
    Bac_Amntd_select ~ ln.HNF_q2
    HNF_Amntd_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0 + ln.HNF_q2 + ln.HNF_Biom
    ln.HNF_Biom ~ ln.HNF_q2 + ln.Bac_q0
    
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q2
'

Bacq0_HNFq2_mod1.10 <- '
  # regressions
    ln.Bac_q0 ~ Bac_Amntd_select + Bac_disp
    ln.HNF_q2 ~ HNF_Amntd_select + HNF_disp + ln.Bac_Biom
    Bac_Amntd_select ~ ln.HNF_q2
    HNF_Amntd_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0 + ln.HNF_Biom
    ln.HNF_Biom ~ ln.HNF_q2 + ln.Bac_q0
    
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q2
'

Bacq0_HNFq2_mod1.11 <- '
  # regressions
    ln.Bac_q0 ~ Bac_Amntd_select + Bac_disp + ln.HNF_Biom
    ln.HNF_q2 ~ HNF_Amntd_select + HNF_disp + ln.Bac_Biom
    Bac_Amntd_select ~ ln.HNF_q2
    HNF_Amntd_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0 + ln.HNF_Biom
    ln.HNF_Biom ~ ln.HNF_q2
    
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q2
'

Bacq0_HNFq2_mod1.12 <- '
  # regressions
    ln.Bac_q0 ~ Bac_Amntd_select + Bac_disp + ln.HNF_Biom
    ln.HNF_q2 ~ HNF_Amntd_select + HNF_disp
    Bac_Amntd_select ~ ln.HNF_q2
    HNF_Amntd_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0 + ln.HNF_q2 + ln.HNF_Biom
    ln.HNF_Biom ~ ln.HNF_q2
    
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q2
'
Bacq0_HNFq2_lavaan1.0 <- sem(Bacq0_HNFq2_mod1.0, data = HNF_Bac_A)#, se = "bootstrap")
Bacq0_HNFq2_lavaan1.1 <- sem(Bacq0_HNFq2_mod1.1, data = HNF_Bac_A)
Bacq0_HNFq2_lavaan1.2 <- sem(Bacq0_HNFq2_mod1.2, data = HNF_Bac_A)
Bacq0_HNFq2_lavaan1.3 <- sem(Bacq0_HNFq2_mod1.3, data = HNF_Bac_A)
Bacq0_HNFq2_lavaan1.4 <- sem(Bacq0_HNFq2_mod1.4, data = HNF_Bac_A)
Bacq0_HNFq2_lavaan1.5 <- sem(Bacq0_HNFq2_mod1.5, data = HNF_Bac_A)
Bacq0_HNFq2_lavaan1.6 <- sem(Bacq0_HNFq2_mod1.6, data = HNF_Bac_A)
Bacq0_HNFq2_lavaan1.7 <- sem(Bacq0_HNFq2_mod1.7, data = HNF_Bac_A)
Bacq0_HNFq2_lavaan1.8 <- sem(Bacq0_HNFq2_mod1.8, data = HNF_Bac_A)
Bacq0_HNFq2_lavaan1.9 <- sem(Bacq0_HNFq2_mod1.9, data = HNF_Bac_A)
Bacq0_HNFq2_lavaan1.10 <- sem(Bacq0_HNFq2_mod1.10, data = HNF_Bac_A)
Bacq0_HNFq2_lavaan1.11 <- sem(Bacq0_HNFq2_mod1.11, data = HNF_Bac_A)
Bacq0_HNFq2_lavaan1.12 <- sem(Bacq0_HNFq2_mod1.12, data = HNF_Bac_A)

AICstep1 <- AIC(Bacq0_HNFq2_lavaan1.0, Bacq0_HNFq2_lavaan1.1, Bacq0_HNFq2_lavaan1.2, Bacq0_HNFq2_lavaan1.3,
                Bacq0_HNFq2_lavaan1.4, Bacq0_HNFq2_lavaan1.5, Bacq0_HNFq2_lavaan1.6, Bacq0_HNFq2_lavaan1.7,
                Bacq0_HNFq2_lavaan1.8, Bacq0_HNFq2_lavaan1.9, Bacq0_HNFq2_lavaan1.10, Bacq0_HNFq2_lavaan1.11,
                Bacq0_HNFq2_lavaan1.12)
AICstep1 <- AICstep1 %>% cbind(row.names(AICstep1)) %>%
  arrange(AIC)
AICstep1

Bacq0_HNFq2_mod1.1 <- '
  # regressions
    ln.Bac_q0 ~ Bac_Amntd_select + Bac_disp
    ln.HNF_q2 ~ HNF_Amntd_select + HNF_disp
    Bac_Amntd_select ~ ln.HNF_q2
    HNF_Amntd_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0 + ln.HNF_q2
    ln.HNF_Biom ~ ln.HNF_q2 + ln.Bac_q0
  
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q2
    ln.Bac_Biom ~~ ln.HNF_Biom
'
Bacq0_HNFq2_lavaan1.1 <- sem(Bacq0_HNFq2_mod1.1, data = HNF_Bac_A)
summary(Bacq0_HNFq2_lavaan1.1, fit.measures = TRUE)
summary(Bacq0_HNFq2_lavaan1.0, fit.measures = TRUE)
##### Step 1: no random effects #####
##### Step 2 : include grouping variables (random effects) and environmental variables #####
Bacq0_HNFq0_psem2.0 <- psem(
  lm(ln.Bac_q0 ~ Bac_select + Bac_disp + ln.HNF_Biom, # + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      data = HNF_Bac_A),
  lm(ln.HNF_q0 ~ HNF_select + HNF_disp + ln.Bac_Biom, # + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      data = HNF_Bac_A),
  lm(Bac_select ~ ln.HNF_q0 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla,
      data = HNF_Bac_A),
  # lm(HNF_select ~ ln.Bac_q0 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla,
  #     data = HNF_Bac_A),

  # lm(ln.Bac_Biom ~ ln.Bac_q0 + ln.HNF_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla,
  #     data = HNF_Bac_A),
  lm(ln.HNF_Biom ~ ln.HNF_q0+ ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla,
      data = HNF_Bac_A),
  
  ln.Bac_q0 %~~% ln.HNF_q0
)
summary(Bacq0_HNFq0_psem2.0)
Bacq0_HNFq0_psem2.Cr <- psem(
  lme(ln.Bac_q0 ~ Bac_select + Bac_disp + ln.HNF_Biom, # + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Cruise, data = HNF_Bac_A),
  lme(ln.HNF_q0 ~ HNF_select + HNF_disp + ln.Bac_Biom, # + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Cruise, data = HNF_Bac_A),
  lme(Bac_select ~ ln.HNF_q0 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Cruise, data = HNF_Bac_A),
  # lme(HNF_select ~ ln.Bac_q0 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
  #     random = ~ 1 | Cruise, data = HNF_Bac_A),
  # 
  # lme(ln.Bac_Biom ~ ln.Bac_q0 + ln.HNF_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
  #    random = ~ 1 | Cruise, data = HNF_Bac_A),
  lme(ln.HNF_Biom ~ ln.HNF_q0 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Cruise, data = HNF_Bac_A),
  
  ln.Bac_q0 %~~% ln.HNF_q0
)
summary(Bacq0_HNFq0_psem2.Cr)
Bacq0_HNFq0_psem2.Season <- psem(
  lme(ln.Bac_q0 ~ Bac_select + Bac_disp + ln.HNF_Biom, #  + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Season, data = HNF_Bac_A),
  lme(ln.HNF_q0 ~ HNF_select + HNF_disp + ln.Bac_Biom, #  + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Season, data = HNF_Bac_A),
  lme(Bac_select ~ ln.HNF_q0 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Season, data = HNF_Bac_A),
  # lme(HNF_select ~ ln.Bac_q0 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
  #     random = ~ 1 | Season, data = HNF_Bac_A),
  # 
  # lme(ln.Bac_Biom ~ ln.Bac_q0 + ln.HNF_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla,
  #     random = ~ 1 | Season, data = HNF_Bac_A),
  lme(ln.HNF_Biom ~ ln.HNF_q0+ ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Season, data = HNF_Bac_A),
  
  ln.Bac_q0 %~~% ln.HNF_q0
)
summary(Bacq0_HNFq0_psem2.Season)
Bacq0_HNFq0_psem2.St <- psem(
  lme(ln.Bac_q0 ~ Bac_select + Bac_disp + ln.HNF_Biom, #  + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Station, data = HNF_Bac_A),
  lme(ln.HNF_q0 ~ HNF_select + HNF_disp + ln.Bac_Biom, #  + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Station, data = HNF_Bac_A),
  lme(Bac_select ~ ln.HNF_q0 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Station, data = HNF_Bac_A),
  # lme(HNF_select ~ ln.Bac_q0 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
  #     random = ~ 1 | Station, data = HNF_Bac_A),
  # 
  # lme(ln.Bac_Biom ~ ln.Bac_q0 + ln.HNF_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla,
  #     random = ~ 1 | Station, data = HNF_Bac_A),
  lme(ln.HNF_Biom ~ ln.HNF_q0 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Station, data = HNF_Bac_A),
  
  ln.Bac_q0 %~~% ln.HNF_q0
)
summary(Bacq0_HNFq0_psem2.St)

anova(Bacq0_HNFq0_psem2.0, Bacq0_HNFq0_psem2.Cr, Bacq0_HNFq0_psem2.Season, Bacq0_HNFq0_psem2.St, test = "F")
anova(Bacq0_HNFq0_psem2.0, Bacq0_HNFq0_psem2.Cr)
summary(Bacq0_HNFq0_psem2.Cr, fit.measures = TRUE)
##### Path model analysis : Bac_q0 vs HNF_q0 ##########
###############################################################################################
##### Testing hypothesis: HNFq0 -> Bac selection -> Bacq0 & Bacq0 -> Bac selection ->  HNFq0 ##
###############################################################################################
