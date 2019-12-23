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
Bac_MNTD_null <- read.table(file = "D:/Research/PdPy_Div_Results/Bac_MNTD_null.csv", sep = ",", 
                            header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
Bac_MNTD <- Bac_MNTD_null %>% 
  select(c(obs, Var1, Var2)) %>%
  mutate(MNTD_null_mean = apply(Bac_MNTD_null[, !names(Bac_MNTD_null) %in% c("obs", "Var1", "Var2")], 1, mean),
         MNTD_null_sd = apply(Bac_MNTD_null[, !names(Bac_MNTD_null) %in% c("obs", "Var1", "Var2")], 1, sd),
         Bac_select_strength = (obs - MNTD_null_mean) / MNTD_null_sd,
         Bac_select_p = pnorm(-abs(Bac_select_strength), 0, 1))

Bac_Chao_null <- read.table(file = "D:/Research/PdPy_Div_Results/Bac_Chao_null.csv", sep = ",", 
                            header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
Bac_BDiv_Chao <- Bac_Chao_null %>% 
  select(c(obs, Var1, Var2)) %>%
  mutate(Chao_null_mean = apply(Bac_Chao_null[, !names(Bac_Chao_null) %in% c("obs", "Var1", "Var2")], 1, mean),
         Chao_null_sd = apply(Bac_Chao_null[, !names(Bac_Chao_null) %in% c("obs", "Var1", "Var2")], 1, sd),
         Bac_disp_strength = (obs - Chao_null_mean) / Chao_null_sd,
         Bac_disp_p = pnorm(Bac_disp_strength, 0, 1))

HNF_MNTD_null <- read.table(file = "D:/Research/PdPy_Div_Results/HNF_MNTD_null_PR2.csv", sep = ",", 
                            header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
HNF_MNTD <- HNF_MNTD_null %>% select(c(obs, Var1, Var2)) %>%
  mutate(MNTD_null_mean = apply(HNF_MNTD_null[, !names(HNF_MNTD_null) %in% c("obs", "Var1", "Var2")], 1, mean),
         MNTD_null_sd = apply(HNF_MNTD_null[, !names(HNF_MNTD_null) %in% c("obs", "Var1", "Var2")], 1, sd),
         HNF_select_strength = (obs - MNTD_null_mean) / MNTD_null_sd,
         HNF_select_p = pnorm(-abs(HNF_select_strength), 0, 1))
HNF_Chao_null <- read.table(file = "D:/Research/PdPy_Div_Results/HNF_Chao_null_PR2.csv", sep = ",", 
                            header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
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
##### Alpha level analyses ####################################################################
###############################################################################################
##### Preping data ##########
Bac_A <- iNEXT(t(Bac_comm), q = 0, datatype = "abundance", size = max(colSums(Bac_comm)) + 100000)$AsyEst %>% 
  select(Site, Diversity, Estimator) %>% 
  spread(Diversity, Estimator) %>%
  rename(Bac_SR = "Species richness", Bac_Shannon = "Shannon diversity", Bac_Simpson = "Simpson diversity") %>%
  mutate(Site = rownames(Bac_comm))

HNF_A <- iNEXT(t(HNF_comm), q = 0, datatype = "abundance", size = max(colSums(HNF_comm)) + 100000)$AsyEst %>% 
  select(Site, Diversity, Estimator) %>% 
  spread(Diversity, Estimator) %>%
  rename(HNF_SR = "Species richness", HNF_Shannon = "Shannon diversity", HNF_Simpson = "Simpson diversity") %>%
  mutate(Site = rownames(HNF_comm))

Bac_selec <- Bac_MNTD %>%
  filter(Bac_select_p < 0.05) %>%
  group_by(Var2) %>%
  summarize(Bac_select = mean(Bac_select_strength, na.rm = TRUE))

HNF_selec <- HNF_MNTD %>%
  filter(HNF_select_p < 0.05) %>%
  group_by(Var2) %>%
  summarize(HNF_select = mean(HNF_select_strength, na.rm = TRUE))

HNF_Bac_A <- Bac_selec %>%
  inner_join(HNF_selec, by = c("Var2" = "Var2")) %>%
  inner_join(Bac_A, by = c("Var2" = "Site")) %>%
  inner_join(HNF_A, by = c("Var2" = "Site")) %>%
  inner_join(Vars, by = c("Var2" = "SampleID")) %>%
  filter(!is.na(NF_Biom)) %>%
  mutate(ln.Bac_SR = log(Bac_SR),
         ln.HNF_SR = log(HNF_SR),
         ln.Bac_Simpson = log(Bac_Simpson),
         ln.HNF_Simpson = log(HNF_Simpson),
         ln.Bac_Shannon = log(Bac_Shannon),
         ln.HNF_Shannon = log(HNF_Shannon),
         ln.Bac_Biom = log(Bac_Biom),
         ln.HNF_Biom = log(HNF_Biom),
         ln.Temp = log(Temp),
         ln.Sal = log(Sal),
         ln.PAR = log(PAR),
         ln.NO2 = log(NO2 + 0.0001),
         ln.NO3 = log(NO3 + 0.0001),
         ln.DIN = log(DIN + 0.0001),
         ln.PO3 = log(PO3 + 0.0001), 
         ln.Chla = log(Chla + 0.0001))
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
# names(HNF_Bac_A)
# p_Adiv_pairs <- HNF_Bac_A %>%
#   select(Bac_Shannon, HNF_Shannon, ln.Bac_Shannon, ln.HNF_Shannon, ln.Bac_Biom, ln.HNF_Biom, Bac_select, HNF_select) %>%
#     ggpairs(columns = c("Bac_Shannon", "HNF_Shannon", "ln.Bac_Shannon", "ln.HNF_Shannon", 
#                         "ln.Bac_Biom", "ln.HNF_Biom", "Bac_select", "HNF_select"),
#           columnLabels = c("Bacteria\nShannon diversity", "HNF\nShannon diversity", 
#                            "log(Bacteria\nShannon diversity)", "log(HNF\nShannon diversity)",
#                            "log(Bacteria\nbiomass)", "log(HNF\nbiomass)", "Bacteria\nselection", "HNF\nselection"),
#           #mapping = ggplot2::aes(colour = Cruise),
#           upper = list(continuous = cor_fun),
#           lower = list(continuous = fit_fun)) +
#   theme(strip.text.x = element_text(color = "black", size = 14),
#         strip.text.y = element_text(angle = 45, color = "black", size = 14))
# p_Adiv_pairs
# #ggsave(p_Adiv_pairs, file = "D:/Research/PdPy_Div_Results/p_ADiv_pairs.jpeg", dpi = 600, width = 34, height = 28, units = "cm")
p_Adiv_pairs <- HNF_Bac_A %>%
  ggpairs(columns = c("ln.Bac_SR", "ln.HNF_SR", "ln.Bac_Shannon", "ln.HNF_Shannon", "ln.Bac_Simpson", "ln.HNF_Simpson",  
                      "ln.Bac_Biom", "ln.HNF_Biom", "Bac_select", "HNF_select"),
          columnLabels = c("Bacteria\nspecies\nrichness", "HNF\nspecies\nrichness",
                           "Bacteria\nShannon\ndiversity", "HNF\nShannon\ndiversity", 
                           "Bacteria\nSimpson\ndiversity", "HNF\nSimpson\ndiversity", 
                           "log(Bacteria\nbiomass)", "log(HNF\nbiomass)", "Bacteria\nselection", "HNF\nselection"),
          #mapping = ggplot2::aes(colour = Cruise),
          upper = list(continuous = cor_fun),
          lower = list(continuous = fit_fun)) +
  theme(strip.text.x = element_text(color = "black", size = 14),
        strip.text.y = element_text(angle = 45, color = "black", size = 14))
p_Adiv_pairs
# ggsave(p_Adiv_pairs, file = "D:/Research/PdPy_Div_Results/p_ADiv_pairs_ln.jpeg", 
#        dpi = 600, width = 34, height = 28, units = "cm")

# zooming 
p_Adiv_pairs <- HNF_Bac_A %>%
  ggpairs(columns = c("ln.Bac_SR", "ln.HNF_SR", "ln.HNF_Shannon", "ln.HNF_Simpson",  
                      "ln.Bac_Biom", "ln.HNF_Biom", "Bac_select", "HNF_select"),
          columnLabels = c("Bacteria\nspecies\nrichness", 
                           "HNF\nspecies\nrichness", "HNF\nShannon\ndiversity", "HNF\nSimpson\ndiversity", 
                           "log(Bacteria\nbiomass)", "log(HNF\nbiomass)", "Bacteria\nselection", "HNF\nselection"),
          #mapping = ggplot2::aes(colour = Cruise),
          upper = list(continuous = cor_fun),
          lower = list(continuous = fit_fun)) +
  theme(strip.text.x = element_text(color = "black", size = 14),
        strip.text.y = element_text(angle = 45, color = "black", size = 14))
p_Adiv_pairs
# ggsave(p_Adiv_pairs, file = "D:/Research/PdPy_Div_Results/p_ADiv_pairsZoomIn_ln.jpeg", 
#        dpi = 600, width = 34, height = 28, units = "cm")

##### Pair-wise plot of bio-variables ##########

##### Path model analysis : Bac_q0 vs HNF_q1 ##########
### Step 1 
Bacq0_HNFq0_mod1.0 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon
    ln.Bac_Biom ~ ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_Biom + ln.Bac_SR + ln.HNF_Shannon
'
Bacq0_HNFq0_mod1.1 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.HNF_Biom
    ln.Bac_Biom ~ ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_Biom + ln.HNF_Shannon
'
Bacq0_HNFq0_mod1.2 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon
    ln.Bac_Biom ~ ln.Bac_SR
    ln.HNF_Biom ~ ln.Bac_Biom + ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_Biom
'
Bacq0_HNFq0_mod1.3 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.HNF_Biom
    ln.Bac_Biom ~ ln.Bac_SR
    ln.HNF_Biom ~ ln.Bac_Biom + ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_Biom
'
Bacq0_HNFq0_mod1.4 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon
    ln.Bac_Biom ~ ln.HNF_Biom + ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_SR + ln.HNF_Shannon
'
Bacq0_HNFq0_mod1.5 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.HNF_Biom
    ln.Bac_Biom ~ ln.HNF_Biom + ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Biom ~ ln.HNF_Shannon
'
Bacq0_HNFq0_mod1.6 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon
    ln.Bac_Biom ~ ln.HNF_Biom + ln.Bac_SR
    ln.HNF_Biom ~ ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_Biom
'
Bacq0_HNFq0_mod1.7 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.HNF_Biom
    ln.Bac_Biom ~ ln.HNF_Biom + ln.Bac_SR
    ln.HNF_Biom ~ ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_Biom
'
Bacq0_HNFq0_mod1.8 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.Bac_Biom
    ln.Bac_Biom ~ ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_Biom + ln.Bac_SR + ln.HNF_Shannon
'
Bacq0_HNFq0_mod1.9 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.Bac_Biom + ln.HNF_Biom
    ln.Bac_Biom ~ ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_Biom + ln.HNF_Shannon
'
Bacq0_HNFq0_mod1.10 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.Bac_Biom
    ln.HNF_Biom ~ ln.Bac_Biom + ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_Biom
'
Bacq0_HNFq0_mod1.11 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.Bac_Biom + ln.HNF_Biom
    ln.HNF_Biom ~ ln.Bac_Biom + ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_Biom
'
Bacq0_HNFq0_mod1.12 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.Bac_Biom
    ln.Bac_Biom ~ ln.HNF_Biom + ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_SR + ln.HNF_Shannon
'
Bacq0_HNFq0_mod1.13 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.Bac_Biom + ln.HNF_Biom
    ln.Bac_Biom ~ ln.HNF_Biom + ln.HNF_Shannon
    ln.HNF_Biom ~ ln.HNF_Shannon
'
Bacq0_HNFq0_mod1.14 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.Bac_Biom
    ln.Bac_Biom ~ ln.HNF_Biom
    ln.HNF_Biom ~ ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_Biom
'
Bacq0_HNFq0_mod1.15 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.Bac_Biom + ln.HNF_Biom
    ln.Bac_Biom ~ ln.HNF_Biom
    ln.HNF_Biom ~ ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_Biom
'
Bacq0_HNFq0_mod1.16 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon
    ln.Bac_Biom ~ ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_Biom + ln.Bac_SR
    ln.HNF_Shannon ~ ln.HNF_Biom
'
Bacq0_HNFq0_mod1.17 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.HNF_Biom
    ln.Bac_Biom ~ ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_Biom
    ln.HNF_Shannon ~ ln.HNF_Biom
'
Bacq0_HNFq0_mod1.18 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon
    ln.Bac_Biom ~ ln.Bac_SR
    ln.HNF_Biom ~ ln.Bac_Biom + ln.Bac_SR
    ln.HNF_Shannon ~ ln.Bac_Biom + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.19 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.HNF_Biom
    ln.Bac_Biom ~ ln.Bac_SR
    ln.HNF_Biom ~ ln.Bac_Biom
    ln.HNF_Shannon ~ ln.Bac_Biom + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.20 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon
    ln.Bac_Biom ~ ln.HNF_Biom + ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_SR
    ln.HNF_Shannon ~ ln.HNF_Biom
'
Bacq0_HNFq0_mod1.21 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.HNF_Biom
    ln.Bac_Biom ~ ln.HNF_Biom + ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.HNF_Biom
'
Bacq0_HNFq0_mod1.22 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon
    ln.Bac_Biom ~ ln.HNF_Biom + ln.Bac_SR
    ln.HNF_Biom ~ ln.Bac_SR
    ln.HNF_Shannon ~ ln.Bac_Biom + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.23 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.HNF_Biom
    ln.Bac_Biom ~ ln.HNF_Biom + ln.Bac_SR
    ln.HNF_Shannon ~ ln.Bac_Biom + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.24 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.Bac_Biom
    ln.Bac_Biom ~ ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_Biom + ln.Bac_SR
    ln.HNF_Shannon ~ ln.HNF_Biom
'
Bacq0_HNFq0_mod1.25 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.Bac_Biom + ln.HNF_Biom
    ln.Bac_Biom ~ ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_Biom
    ln.HNF_Shannon ~ ln.HNF_Biom
'
Bacq0_HNFq0_mod1.26 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.Bac_Biom
    ln.HNF_Biom ~ ln.Bac_Biom + ln.Bac_SR
    ln.HNF_Shannon ~ ln.Bac_Biom + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.27 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.Bac_Biom + ln.HNF_Biom
    ln.HNF_Biom ~ ln.Bac_Biom
    ln.HNF_Shannon ~ ln.Bac_Biom + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.28 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.Bac_Biom
    ln.Bac_Biom ~ ln.HNF_Biom + ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_SR
    ln.HNF_Shannon ~ ln.HNF_Biom
'
Bacq0_HNFq0_mod1.29 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.Bac_Biom + ln.HNF_Biom
    ln.Bac_Biom ~ ln.HNF_Biom + ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.HNF_Biom
'
Bacq0_HNFq0_mod1.30 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.Bac_Biom
    ln.Bac_Biom ~ ln.HNF_Biom
    ln.HNF_Biom ~ ln.Bac_SR
    ln.HNF_Shannon ~ ln.Bac_Biom + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.31 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Shannon + ln.Bac_Biom + ln.HNF_Biom
    ln.Bac_Biom ~ ln.HNF_Biom
    ln.HNF_Shannon ~ ln.Bac_Biom + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.32 <- '
  # regressions
    ln.Bac_Biom ~ ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_Biom + ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_SR
'
Bacq0_HNFq0_mod1.33 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Biom
    ln.Bac_Biom ~ ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_Biom + ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_SR
'
Bacq0_HNFq0_mod1.34 <- '
  # regressions
    ln.Bac_Biom ~ ln.Bac_SR
    ln.HNF_Biom ~ ln.Bac_Biom + ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_SR + ln.Bac_Biom
'
Bacq0_HNFq0_mod1.35 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Biom
    ln.Bac_Biom ~ ln.Bac_SR
    ln.HNF_Biom ~ ln.Bac_Biom + ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_SR + ln.Bac_Biom
'
Bacq0_HNFq0_mod1.36 <- '
  # regressions
    ln.Bac_Biom ~ ln.HNF_Biom + ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_SR
'
Bacq0_HNFq0_mod1.37 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Biom
    ln.Bac_Biom ~ ln.HNF_Biom + ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Biom ~ ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_SR
'
Bacq0_HNFq0_mod1.38 <- '
  # regressions
    ln.Bac_Biom ~ ln.HNF_Biom + ln.Bac_SR
    ln.HNF_Biom ~ ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_SR + ln.Bac_Biom
'
Bacq0_HNFq0_mod1.39 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Biom
    ln.Bac_Biom ~ ln.HNF_Biom + ln.Bac_SR
    ln.HNF_Biom ~ ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_SR + ln.Bac_Biom
'
Bacq0_HNFq0_mod1.40 <- '
  # regressions
    ln.Bac_SR ~ ln.Bac_Biom
    ln.Bac_Biom ~ ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_Biom + ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_SR
'
Bacq0_HNFq0_mod1.41 <- '
  # regressions
    ln.Bac_SR ~ ln.Bac_Biom + ln.HNF_Biom
    ln.Bac_Biom ~ ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_Biom + ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_SR
'
Bacq0_HNFq0_mod1.42 <- '
  # regressions
    ln.Bac_SR ~ ln.Bac_Biom
    ln.HNF_Biom ~ ln.Bac_Biom + ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_SR + ln.Bac_Biom
'
Bacq0_HNFq0_mod1.43 <- '
  # regressions
    ln.Bac_SR ~ ln.Bac_Biom + ln.HNF_Biom
    ln.HNF_Biom ~ ln.Bac_Biom + ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_SR + ln.Bac_Biom
'
Bacq0_HNFq0_mod1.44 <- '
  # regressions
    ln.Bac_SR ~ ln.Bac_Biom
    ln.Bac_Biom ~ ln.HNF_Biom + ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_SR
'
Bacq0_HNFq0_mod1.45 <- '
  # regressions
    ln.Bac_SR ~ ln.Bac_Biom + ln.HNF_Biom
    ln.Bac_Biom ~ ln.HNF_Biom + ln.HNF_Shannon
    ln.HNF_Biom ~ ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_SR
'
Bacq0_HNFq0_mod1.46 <- '
  # regressions
    ln.Bac_SR ~ ln.Bac_Biom
    ln.Bac_Biom ~ ln.HNF_Biom
    ln.HNF_Biom ~ ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_SR + ln.Bac_Biom
'
Bacq0_HNFq0_mod1.47 <- '
  # regressions
    ln.Bac_SR ~ ln.Bac_Biom + ln.HNF_Biom
    ln.Bac_Biom ~ ln.HNF_Biom
    ln.HNF_Biom ~ ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_SR + ln.Bac_Biom
'
Bacq0_HNFq0_mod1.48 <- '
  # regressions
    ln.Bac_Biom ~ ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_Biom + ln.Bac_SR
    ln.HNF_Shannon ~ ln.Bac_SR + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.49 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Biom
    ln.Bac_Biom ~ ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_Biom
    ln.HNF_Shannon ~ ln.Bac_SR + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.50 <- '
  # regressions
    ln.Bac_Biom ~ ln.Bac_SR
    ln.HNF_Biom ~ ln.Bac_Biom + ln.Bac_SR
    ln.HNF_Shannon ~ ln.Bac_SR + ln.Bac_Biom + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.51 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Biom
    ln.Bac_Biom ~ ln.Bac_SR
    ln.HNF_Biom ~ ln.Bac_Biom
    ln.HNF_Shannon ~ ln.Bac_SR + ln.Bac_Biom + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.52 <- '
  # regressions
    ln.Bac_Biom ~ ln.HNF_Biom + ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_SR
    ln.HNF_Shannon ~ ln.Bac_SR + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.53 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Biom
    ln.Bac_Biom ~ ln.HNF_Biom + ln.Bac_SR + ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_SR + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.54 <- '
  # regressions
    ln.Bac_Biom ~ ln.HNF_Biom + ln.Bac_SR
    ln.HNF_Biom ~ ln.Bac_SR
    ln.HNF_Shannon ~ ln.Bac_SR + ln.Bac_Biom + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.55 <- '
  # regressions
    ln.Bac_SR ~ ln.HNF_Biom
    ln.Bac_Biom ~ ln.HNF_Biom + ln.Bac_SR
    ln.HNF_Shannon ~ ln.Bac_SR + ln.Bac_Biom + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.56 <- '
  # regressions
    ln.Bac_SR ~ ln.Bac_Biom
    ln.Bac_Biom ~ ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_Biom + ln.Bac_SR
    ln.HNF_Shannon ~ ln.Bac_SR + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.57 <- '
  # regressions
    ln.Bac_SR ~ ln.Bac_Biom + ln.HNF_Biom
    ln.Bac_Biom ~ ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_Biom
    ln.HNF_Shannon ~ ln.Bac_SR + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.58 <- '
  # regressions
    ln.Bac_SR ~ ln.Bac_Biom
    ln.HNF_Biom ~ ln.Bac_Biom + ln.Bac_SR
    ln.HNF_Shannon ~ ln.Bac_SR + ln.Bac_Biom + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.59 <- '
  # regressions
    ln.Bac_SR ~ ln.Bac_Biom + ln.HNF_Biom
    ln.HNF_Biom ~ ln.Bac_Biom
    ln.HNF_Shannon ~ ln.Bac_SR + ln.Bac_Biom + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.60 <- '
  # regressions
    ln.Bac_SR ~ ln.Bac_Biom
    ln.Bac_Biom ~ ln.HNF_Biom + ln.HNF_Shannon
    ln.HNF_Biom ~ ln.Bac_SR
    ln.HNF_Shannon ~ ln.Bac_SR + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.61 <- '
  # regressions
    ln.Bac_SR ~ ln.Bac_Biom + ln.HNF_Biom
    ln.Bac_Biom ~ ln.HNF_Biom + ln.HNF_Shannon
    ln.HNF_Shannon ~ ln.Bac_SR + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.62 <- '
  # regressions
    ln.Bac_SR ~ ln.Bac_Biom
    ln.Bac_Biom ~ ln.HNF_Biom
    ln.HNF_Biom ~ ln.Bac_SR
    ln.HNF_Shannon ~ ln.Bac_SR + ln.Bac_Biom + ln.HNF_Biom
'
Bacq0_HNFq0_mod1.63 <- '
  # regressions
    ln.Bac_SR ~ ln.Bac_Biom + ln.HNF_Biom
    ln.Bac_Biom ~ ln.HNF_Biom
    ln.HNF_Shannon ~ ln.Bac_SR + ln.Bac_Biom + ln.HNF_Biom
'


Bacq0_HNFq0_lavaan1.0 <- sem(Bacq0_HNFq0_mod1.0, data = HNF_Bac_A)#, se = "bootstrap")
Bacq0_HNFq0_lavaan1.1 <- sem(Bacq0_HNFq0_mod1.1, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.2 <- sem(Bacq0_HNFq0_mod1.2, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.3 <- sem(Bacq0_HNFq0_mod1.3, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.4 <- sem(Bacq0_HNFq0_mod1.4, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.5 <- sem(Bacq0_HNFq0_mod1.5, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.6 <- sem(Bacq0_HNFq0_mod1.6, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.7 <- sem(Bacq0_HNFq0_mod1.7, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.8 <- sem(Bacq0_HNFq0_mod1.8, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.9 <- sem(Bacq0_HNFq0_mod1.9, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.10 <- sem(Bacq0_HNFq0_mod1.10, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.11 <- sem(Bacq0_HNFq0_mod1.11, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.12 <- sem(Bacq0_HNFq0_mod1.12, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.13 <- sem(Bacq0_HNFq0_mod1.13, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.14 <- sem(Bacq0_HNFq0_mod1.14, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.15 <- sem(Bacq0_HNFq0_mod1.15, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.16 <- sem(Bacq0_HNFq0_mod1.16, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.17 <- sem(Bacq0_HNFq0_mod1.17, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.18 <- sem(Bacq0_HNFq0_mod1.18, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.19 <- sem(Bacq0_HNFq0_mod1.19, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.20 <- sem(Bacq0_HNFq0_mod1.20, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.21 <- sem(Bacq0_HNFq0_mod1.21, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.22 <- sem(Bacq0_HNFq0_mod1.22, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.23 <- sem(Bacq0_HNFq0_mod1.23, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.24 <- sem(Bacq0_HNFq0_mod1.24, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.25 <- sem(Bacq0_HNFq0_mod1.25, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.26 <- sem(Bacq0_HNFq0_mod1.26, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.27 <- sem(Bacq0_HNFq0_mod1.27, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.28 <- sem(Bacq0_HNFq0_mod1.28, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.29 <- sem(Bacq0_HNFq0_mod1.29, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.30 <- sem(Bacq0_HNFq0_mod1.30, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.31 <- sem(Bacq0_HNFq0_mod1.31, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.32 <- sem(Bacq0_HNFq0_mod1.32, data = HNF_Bac_A)#, se = "bootstrap")
Bacq0_HNFq0_lavaan1.33 <- sem(Bacq0_HNFq0_mod1.33, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.34 <- sem(Bacq0_HNFq0_mod1.34, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.35 <- sem(Bacq0_HNFq0_mod1.35, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.36 <- sem(Bacq0_HNFq0_mod1.36, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.37 <- sem(Bacq0_HNFq0_mod1.37, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.38 <- sem(Bacq0_HNFq0_mod1.38, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.39 <- sem(Bacq0_HNFq0_mod1.39, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.40 <- sem(Bacq0_HNFq0_mod1.40, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.41 <- sem(Bacq0_HNFq0_mod1.41, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.42 <- sem(Bacq0_HNFq0_mod1.42, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.43 <- sem(Bacq0_HNFq0_mod1.43, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.44 <- sem(Bacq0_HNFq0_mod1.44, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.45 <- sem(Bacq0_HNFq0_mod1.45, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.46 <- sem(Bacq0_HNFq0_mod1.46, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.47 <- sem(Bacq0_HNFq0_mod1.47, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.48 <- sem(Bacq0_HNFq0_mod1.48, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.49 <- sem(Bacq0_HNFq0_mod1.49, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.50 <- sem(Bacq0_HNFq0_mod1.50, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.51 <- sem(Bacq0_HNFq0_mod1.51, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.52 <- sem(Bacq0_HNFq0_mod1.52, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.53 <- sem(Bacq0_HNFq0_mod1.53, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.54 <- sem(Bacq0_HNFq0_mod1.54, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.55 <- sem(Bacq0_HNFq0_mod1.55, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.56 <- sem(Bacq0_HNFq0_mod1.56, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.57 <- sem(Bacq0_HNFq0_mod1.57, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.58 <- sem(Bacq0_HNFq0_mod1.58, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.59 <- sem(Bacq0_HNFq0_mod1.59, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.60 <- sem(Bacq0_HNFq0_mod1.60, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.61 <- sem(Bacq0_HNFq0_mod1.61, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.62 <- sem(Bacq0_HNFq0_mod1.62, data = HNF_Bac_A)
Bacq0_HNFq0_lavaan1.63 <- sem(Bacq0_HNFq0_mod1.63, data = HNF_Bac_A)

AICstep1 <- AIC(Bacq0_HNFq0_lavaan1.0, Bacq0_HNFq0_lavaan1.1, Bacq0_HNFq0_lavaan1.2, Bacq0_HNFq0_lavaan1.3,
                Bacq0_HNFq0_lavaan1.4, Bacq0_HNFq0_lavaan1.5, Bacq0_HNFq0_lavaan1.6, Bacq0_HNFq0_lavaan1.7,
                Bacq0_HNFq0_lavaan1.8, Bacq0_HNFq0_lavaan1.9, Bacq0_HNFq0_lavaan1.10, Bacq0_HNFq0_lavaan1.11,
                Bacq0_HNFq0_lavaan1.12, Bacq0_HNFq0_lavaan1.13, Bacq0_HNFq0_lavaan1.14, Bacq0_HNFq0_lavaan1.15,
                Bacq0_HNFq0_lavaan1.16, Bacq0_HNFq0_lavaan1.17, Bacq0_HNFq0_lavaan1.18, Bacq0_HNFq0_lavaan1.19,
                Bacq0_HNFq0_lavaan1.20, Bacq0_HNFq0_lavaan1.21, Bacq0_HNFq0_lavaan1.22, Bacq0_HNFq0_lavaan1.23,
                Bacq0_HNFq0_lavaan1.24, Bacq0_HNFq0_lavaan1.25, Bacq0_HNFq0_lavaan1.26, Bacq0_HNFq0_lavaan1.27,
                Bacq0_HNFq0_lavaan1.28, Bacq0_HNFq0_lavaan1.29, Bacq0_HNFq0_lavaan1.30, Bacq0_HNFq0_lavaan1.31, 
                Bacq0_HNFq0_lavaan1.32, Bacq0_HNFq0_lavaan1.33, Bacq0_HNFq0_lavaan1.34, Bacq0_HNFq0_lavaan1.35,
                Bacq0_HNFq0_lavaan1.36, Bacq0_HNFq0_lavaan1.37, Bacq0_HNFq0_lavaan1.38, Bacq0_HNFq0_lavaan1.39,
                Bacq0_HNFq0_lavaan1.40, Bacq0_HNFq0_lavaan1.41, Bacq0_HNFq0_lavaan1.42, Bacq0_HNFq0_lavaan1.43,
                Bacq0_HNFq0_lavaan1.44, Bacq0_HNFq0_lavaan1.45, Bacq0_HNFq0_lavaan1.46, Bacq0_HNFq0_lavaan1.47,
                Bacq0_HNFq0_lavaan1.48, Bacq0_HNFq0_lavaan1.49, Bacq0_HNFq0_lavaan1.50, Bacq0_HNFq0_lavaan1.51,
                Bacq0_HNFq0_lavaan1.52, Bacq0_HNFq0_lavaan1.53, Bacq0_HNFq0_lavaan1.54, Bacq0_HNFq0_lavaan1.55,
                Bacq0_HNFq0_lavaan1.56, Bacq0_HNFq0_lavaan1.57, Bacq0_HNFq0_lavaan1.58, Bacq0_HNFq0_lavaan1.59,
                Bacq0_HNFq0_lavaan1.60, Bacq0_HNFq0_lavaan1.61, Bacq0_HNFq0_lavaan1.62, Bacq0_HNFq0_lavaan1.63)
AICstep1 <- AICstep1 %>% cbind(row.names(AICstep1)) %>%
  arrange(AIC)
# AICstep1
moreFitIndices(Bacq0_HNFq0_lavaan1.21, fit.measures = "all", nPrior = 1)
moreFitIndices(Bacq0_HNFq0_lavaan1.23, fit.measures = "all", nPrior = 1)
moreFitIndices(Bacq0_HNFq0_lavaan1.29, fit.measures = "all", nPrior = 1)
moreFitIndices(Bacq0_HNFq0_lavaan1.31, fit.measures = "all", nPrior = 1)
moreFitIndices(Bacq0_HNFq0_lavaan1.53, fit.measures = "all", nPrior = 1)
moreFitIndices(Bacq0_HNFq0_lavaan1.55, fit.measures = "all", nPrior = 1)
moreFitIndices(Bacq0_HNFq0_lavaan1.61, fit.measures = "all", nPrior = 1)
moreFitIndices(Bacq0_HNFq0_lavaan1.63, fit.measures = "all", nPrior = 1)

### Step 2 : include selection processes as the interaction terms and grouping variables (random effects)
HNF_Bac_A <- HNF_Bac_A %>%
  mutate(HNF_Bac_int = ln.HNF_Shannon * Bac_select,
         HNF_HNF_int = ln.HNF_Shannon * HNF_select)
HNF_Bac_A <- as.data.frame(HNF_Bac_A)

Bacq0_HNFq0_mod2.Cr <- psem(
  lme(ln.Bac_SR ~ ln.HNF_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Cruise, data = HNF_Bac_A),
  lme(ln.Bac_Biom ~ ln.HNF_Biom + ln.Bac_SR + ln.HNF_Shannon + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Cruise, data = HNF_Bac_A),
  lme(ln.HNF_Shannon ~ ln.Bac_SR + ln.HNF_Biom + HNF_Bac_int + HNF_HNF_int + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Cruise, data = HNF_Bac_A),
  lme(ln.HNF_Biom ~ ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Cruise, data = HNF_Bac_A),
  data = HNF_Bac_A
)
Bacq0_HNFq0_mod2.Season <- psem(
  lme(ln.Bac_SR ~ ln.HNF_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Season, data = HNF_Bac_A),
  lme(ln.Bac_Biom ~ ln.HNF_Biom + ln.Bac_SR + ln.HNF_Shannon + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Season, data = HNF_Bac_A),
  lme(ln.HNF_Shannon ~ ln.Bac_SR + ln.HNF_Biom + HNF_Bac_int + HNF_HNF_int + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Season, data = HNF_Bac_A),
  lme(ln.HNF_Biom ~ ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Season, data = HNF_Bac_A),
  data = HNF_Bac_A
)
Bacq0_HNFq0_mod2.St <- psem(
  lme(ln.Bac_SR ~ ln.HNF_Biom + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Station, data = HNF_Bac_A),
  lme(ln.Bac_Biom ~ ln.HNF_Biom + ln.Bac_SR + ln.HNF_Shannon + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Station, data = HNF_Bac_A),
  lme(ln. ~ ln.Bac_SR + ln.HNF_Biom + HNF_Bac_int + HNF_HNF_int + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Station, data = HNF_Bac_A),
  lme(ln.HNF_Biom ~ ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Season, data = HNF_Bac_A),
  data = HNF_Bac_A
)

AIC(Bacq0_HNFq0_mod2.Cr)
AIC(Bacq0_HNFq0_mod2.Season)
AIC(Bacq0_HNFq0_mod2.St)

summary(Bacq0_HNFq0_mod2.Cr)
summary(Bacq0_HNFq0_mod2.Season)
summary(Bacq0_HNFq0_mod2.St)
##### Path model analysis : Bac_q0 vs HNF_q0 ##########


##### Univariate HNF-Bac A diversity relationship ##########
### Plotting
p_ADiv_BacSelect <- HNF_Bac_A %>% 
  ggplot() + 
  geom_point(aes(x = Bac_Shannon, y = HNF_Shannon, color = Bac_select), size = 3) + 
  scale_colour_viridis(alpha = 0.7)
p_ADiv_BacSelect
#ggsave(p_ADiv_BacSelect, file = "D:/Research/PdPy_Div_Results/p_ADiv_BacSelect.jpeg")
p_ADiv_HNFSelect <- HNF_Bac_A %>% 
  ggplot() + 
  geom_point(aes(x = Bac_Shannon, y = HNF_Shannon, color = HNF_select), size = 3) + 
  scale_colour_viridis(alpha = 0.7)
p_ADiv_HNFSelect
#ggsave(p_ADiv_HNFSelect, file = "D:/Research/PdPy_Div_Results/p_ADiv_HNFSelect.jpeg")
### Plotting


##### Univariate HNF-Bac A diversity relationship ##########



##### Analyzing ##########
### Univariate relationships

# Selection vs alpha diversity
p_Bac_Selec <- HNF_Bac_A %>%
  select(ln.Bac_Shannon, Bac_select, HNF_select) %>%
  gather(Community, Select, -ln.Bac_Shannon) %>%
  ggplot(aes(x = Select, y = ln.Bac_Shannon)) + 
    geom_point(size = 3) +
    facet_grid(cols = vars(Community), scales = "free") + 
    geom_smooth(formula = y ~ x, method = "lm", se = TRUE) + 
    geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red") 
    #geom_smooth(formula = y ~ x, method = lm, se = TRUE, )
ggsave(p_Bac_Selec, file = "D:/Research/PdPy_Div_Results/p_Bac_Selec.jpeg")

gam0 <- gam(ln.Bac_Shannon ~ Bac_select + HNF_select + s(Bac_select) + s(HNF_select), data = HNF_Bac_A)
summary(gam0)

lm0_BacADiv_Shannon_Sea <- lme(Bac_Shannon ~ Bac_select + HNF_select, random = ~1 | Season, data = HNF_Bac_A)
lm0_BacADiv_Shannon_Cr <- lme(Bac_Shannon ~ Bac_select + HNF_select, random = ~1 | Cruise, data = HNF_Bac_A)
lm0_BacADiv_Shannon_St <- lme(Bac_Shannon ~ Bac_select + HNF_select, random = ~1 | Station, data = HNF_Bac_A)
AIC(lm0_BacADiv_Shannon_Sea, lm0_BacADiv_Shannon_Cr, lm0_BacADiv_Shannon_St)
summary(lm0_BacADiv_Shannon_Sea)



###############################################################################################
##### Alpha level analyses ####################################################################
###############################################################################################
