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

Bac_selec <- Bac_MNTD %>%
  mutate(Cr_V1 = substr(Var1, start = 1, stop = 9),
         Cr_V2 = substr(Var2, start = 1, stop = 9)) %>%
  filter(Cr_V1 == Cr_V2) %>%
  #filter(Bac_select_p < 0.05) %>%
  group_by(Var2) %>%
  # group_by(St, Cr) %>%
  summarize(Bac_select = mean(Bac_select_strength, na.rm = TRUE))

HNF_selec <- HNF_MNTD %>%
  mutate(Cr_V1 = substr(Var1, start = 1, stop = 9),
         Cr_V2 = substr(Var2, start = 1, stop = 9)) %>%
  filter(Cr_V1 == Cr_V2) %>%
  #filter(HNF_select_p < 0.05) %>%
  group_by(Var2) %>%
  summarize(HNF_select = mean(HNF_select_strength, na.rm = TRUE))

HNF_Bac_A <- Bac_selec %>%
  inner_join(HNF_selec, by = c("Var2" = "Var2")) %>%
  inner_join(Bac_A, by = c("Var2" = "Site")) %>%
  inner_join(HNF_A, by = c("Var2" = "Site")) %>%
  inner_join(Vars, by = c("Var2" = "SampleID")) %>%
  filter(!is.na(NF_Biom)) %>%
  mutate(ln.Bac_q0 = log(Bac_q0),
         ln.HNF_q0 = log(HNF_q0),
         ln.Bac_q1 = log(Bac_q1),
         ln.HNF_q1 = log(HNF_q1),
         ln.Bac_q2 = log(Bac_q2),
         ln.HNF_q2 = log(HNF_q2),
         ln.Bac_Biom = log(Bac_Biom),
         ln.HNF_Biom = log(HNF_Biom),
         ln.Bac_select = log(Bac_select - min(Bac_select) + 0.1),
         ln.HNF_select = log(HNF_select - min(HNF_select) + 0.1),
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
p_Adiv_pairs <- HNF_Bac_A %>%
  ggpairs(columns = c("ln.Bac_q0", "ln.HNF_q0", "ln.Bac_q1", "ln.HNF_q1", "ln.Bac_q2", "ln.HNF_q2",  
                      "ln.Bac_Biom", "ln.HNF_Biom", "Bac_select", "HNF_select"),
          columnLabels = c("Bacteria\nspecies\nrichness", "HNF\nspecies\nrichness",
                           "Bacteria\nShannon\ndiversity", "HNF\nShannon\ndiversity", 
                           "Bacteria\nSimpson\ndiversity", "HNF\nSimpson\ndiversity", 
                           "log(Bacteria\nbiomass)", "log(HNF\nbiomass)", "Bacteria\nbNTI", "HNF\nbNTI"),
          #mapping = ggplot2::aes(colour = Cruise),
          upper = list(continuous = cor_fun),
          lower = list(continuous = fit_fun)) +
  theme(strip.text.x = element_text(color = "black", size = 14),
        strip.text.y = element_text(angle = 45, color = "black", size = 14))
p_Adiv_pairs
# ggsave(p_Adiv_pairs, file = "D:/Research/PdPy_Div_Results/Figs/p_ADiv_pairs_ln.jpeg",
#        dpi = 600, width = 34, height = 28, units = "cm")

# zooming 
p_Adiv_pairs <- HNF_Bac_A %>%
  ggpairs(columns = c("ln.Bac_q0", "ln.HNF_q0", "ln.HNF_q1", "ln.HNF_q2",  
                      "ln.Bac_Biom", "ln.HNF_Biom", "Bac_select", "HNF_select"),
          columnLabels = c("Bacteria\nspecies\nrichness", 
                           "HNF\nspecies\nrichness", "HNF\nShannon\ndiversity", "HNF\nSimpson\ndiversity", 
                           "log(Bacteria\nbiomass)", "log(HNF\nbiomass)", "Bacteria\nbNTI", "HNF\nbNTI"),
          #mapping = ggplot2::aes(colour = Cruise),
          upper = list(continuous = cor_fun),
          lower = list(continuous = fit_fun)) +
  theme(strip.text.x = element_text(color = "black", size = 14),
        strip.text.y = element_text(angle = 45, color = "black", size = 14))
p_Adiv_pairs
# ggsave(p_Adiv_pairs, file = "D:/Research/PdPy_Div_Results/Figs/p_ADiv_pairsZoomIn_ln.jpeg",
#        dpi = 600, width = 34, height = 28, units = "cm")

##### Pair-wise plot of bio-variables ##########
###############################################################################################
##### Envi exploratory factor analyses and  pair-wise correlations ############################
###############################################################################################

###############################################################################################
##### Simple HNF and Bac alppha diversity relationship  #######################################
###############################################################################################
##### without selection #####
### linear model
# Statistical check on non-linear relationship between Bac and HNF alpha diversity
gam0 <- gam(ln.HNF_q0 ~ ln.Bac_q0, data = HNF_Bac_A)
gam1 <- gam(ln.HNF_q0 ~ s(ln.Bac_q0), data = HNF_Bac_A)
anova(gam0, gam1, test = "F")
### plotting (without selection)
p_HNFq0_Bacq0 <- HNF_Bac_A %>% 
  select(ln.Bac_q0, ln.HNF_q0, Bac_select, HNF_select) %>%
  ggplot(aes(x = ln.Bac_q0, y = ln.HNF_q0)) + 
  geom_point(size = 3) + 
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE) + 
  geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "dotted") + 
  scale_colour_viridis(alpha = 0.7) + 
  labs(x = expression("Log[ Bacteria species richness (Hill number = 0) ]"),
       y = expression("Log[ HNF species richness (Hill number = 0) ]")) + 
  theme(
    strip.text.x = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16)
  )
p_HNFq0_Bacq0
# ggsave(p_HNFq0_Bacq0, file = "D:/Research/PdPy_Div_Results/Figs/p_ADiv_Bacq0_HNFq0.jpeg",
#        dpi = 600, width = 34, height = 28, units = "cm")
##### without selection #####

##### add selection #####
### linear model
HNFq0_Bacq0_0 <- glm(ln.HNF_q0 ~ ln.Bac_q0 + HNF_select, data = HNF_Bac_A)
HNFq0_Bacq0_1 <- glm(ln.HNF_q0 ~ ln.Bac_q0 * HNF_select + HNF_select * ln.Bac_q0, data = HNF_Bac_A)
anova(HNFq0_Bacq0_0, HNFq0_Bacq0_1, test = "F")
summary(HNFq0_Bacq0_0)

# Plotting HNF and Bac alppha diversity relationship with selection as color code
community.labs <- c("Colored by selection on bacteria community", "Colored by selection on HNF community")
names(community.labs) <- c("Bac_select", "HNF_select")

p_ADiv_Select <- HNF_Bac_A %>% 
  select(ln.Bac_q0, ln.HNF_q0, Bac_select, HNF_select) %>%
  gather(key = "community", value = "Selection", -c(ln.Bac_q0, ln.HNF_q0)) %>%
  ggplot(aes(x = ln.Bac_q0, y = ln.HNF_q0, color = Selection)) + 
  geom_point(size = 3) + 
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE) + 
  geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "dotted") + 
  facet_grid(~ community, labeller = labeller(community = community.labs)) +  
  scale_colour_viridis(alpha = 0.7) + 
  labs(x = expression("Log[ Bacteria species richness (Hill number = 0) ]"),
       y = expression("Log[ HNF species richness (Hill number = 0) ]"),
       colour = expression(paste("\U03B2", "NTI"))) + 
  theme(
    strip.text.x = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16)
  )
p_ADiv_Select
# ggsave(p_ADiv_Select, file = "D:/Research/PdPy_Div_Results/Figs/p_ADiv_Bacq0_HNFq0_Select.jpeg",
#        dpi = 600, width = 34, height = 28, units = "cm")
##### add selection #####
###############################################################################################
##### Simple HNF and Bac alppha diversity relationship  #######################################
###############################################################################################

###############################################################################################
##### Testing hypothesis: HNFq0 -> Bac selection -> Bacq0 & Bacq0 -> Bac selection ->  HNFq0 ##
###############################################################################################
##### HNFq0 -> Bac selection -> Bacq0 ##########
### linear model testing
BacS_HNFq0.0 <- lm(Bac_select ~ ln.HNF_q0, data = HNF_Bac_A)
BacS_HNFq0.Cr <- lme(Bac_select ~ ln.HNF_q0, random = ~ 1 | Cruise, data = HNF_Bac_A)
BacS_HNFq0.Season <- lme(Bac_select ~ ln.HNF_q0, random = ~ 1 | Season, data = HNF_Bac_A)
BacS_HNFq0.St <- lme(Bac_select ~ ln.HNF_q0, random = ~ 1 | Station, data = HNF_Bac_A)
AIC(BacS_HNFq0.0, BacS_HNFq0.Cr, BacS_HNFq0.Season, BacS_HNFq0.St)
summary(BacS_HNFq0.Cr)
summary(gam(Bac_select ~ s(ln.HNF_q0, bs = "ts"), data = HNF_Bac_A))

Bacq0_BacS.0 <- lm(ln.Bac_q0 ~ Bac_select, data = HNF_Bac_A)
Bacq0_BacS.Cr <- lme(ln.Bac_q0 ~ Bac_select, random = ~ 1 | Cruise, data = HNF_Bac_A)
Bacq0_BacS.Season <- lme(ln.Bac_q0 ~ Bac_select, random = ~ 1 | Season, data = HNF_Bac_A)
Bacq0_BacS.St <- lme(ln.Bac_q0 ~ Bac_select, random = ~ 1 | Station, data = HNF_Bac_A)
AIC(Bacq0_BacS.0, Bacq0_BacS.Cr, Bacq0_BacS.Season, Bacq0_BacS.St)
summary(Bacq0_BacS.Cr)
summary(gam(ln.Bac_q0 ~ s(Bac_select, bs = "ts"), data = HNF_Bac_A))

### plotting
p_HNFq0_BacSelect <- HNF_Bac_A %>% 
  select(ln.Bac_q0, ln.HNF_q0, Bac_select, HNF_select, Cruise) %>%
  ggplot(aes(x = ln.HNF_q0, y = Bac_select)) + 
  geom_point(aes(color = Cruise), size = 3) + 
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "dotted") + 
  geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "dotted") + 
  scale_colour_viridis(alpha = 0.7, discrete=TRUE) + 
  labs(x = expression("Log[ HNF species richness (Hill number = 0) ]"),
       y = expression("Deterministic assembly processes ( \U03B2NTI) of Bacteria community ")) + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))

p_BacSelect_Bacq0 <- HNF_Bac_A %>% 
  select(ln.Bac_q0, ln.HNF_q0, Bac_select, HNF_select, Cruise) %>%
  ggplot(aes(x = Bac_select, y = ln.Bac_q0)) + 
  geom_point(aes(color = Cruise), size = 3) + 
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "dotted") + 
  geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "solid") + 
  scale_colour_viridis(alpha = 0.7, discrete=TRUE) + 
  labs(x = expression("Deterministic assembly processes ( \U03B2NTI) of Bacteria community "),
       y = expression("Log[ Bacteria species richness (Hill number = 0) ]")) + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))

legend <- get_legend(
  p_HNFq0_BacSelect + theme(legend.box.margin = margin(0, 0, 0, 12))
)
p_HNFq0_BacSelect_Bacq0 <- plot_grid(
  plot_grid(p_HNFq0_BacSelect + theme(legend.position = "none"),
            p_BacSelect_Bacq0 + theme(legend.position = "none")),
  legend, rel_widths = c(3, .4))
p_HNFq0_BacSelect_Bacq0
# ggsave(p_HNFq0_BacSelect_Bacq0, file = "D:/Research/PdPy_Div/Presentation/20200110_NTOU/p_HNFq0_BacSelect_Bacq0.png",
#        dpi = 600, width = 56, height = 28, units = "cm")
##### HNFq0 -> Bac selection -> Bacq0 ##########

##### Bacq0 -> Bac selection ->  HNFq0 ##########
### linear model testing
HNFS_Bacq0.0 <- lm(HNF_select ~ ln.Bac_q0, data = HNF_Bac_A)
HNFS_Bacq0.Cr <- lme(HNF_select ~ ln.Bac_q0, random = ~ 1 | Cruise, data = HNF_Bac_A)
HNFS_Bacq0.Season <- lme(HNF_select ~ ln.Bac_q0, random = ~ 1 | Season, data = HNF_Bac_A)
HNFS_Bacq0.St <- lme(HNF_select ~ ln.Bac_q0, random = ~ 1 | Station, data = HNF_Bac_A)
AIC(HNFS_Bacq0.0, HNFS_Bacq0.Cr, HNFS_Bacq0.Season, HNFS_Bacq0.St)
summary(HNFS_Bacq0.Cr)
summary(gam(HNF_select ~ s(ln.Bac_q0, bs = "ts"), data = HNF_Bac_A))

HNFq0_HNFS.0 <- lm(ln.HNF_q0 ~ HNF_select, data = HNF_Bac_A)
HNFq0_HNFS.Cr <- lme(ln.HNF_q0 ~ HNF_select, random = ~ 1 | Cruise, data = HNF_Bac_A)
HNFq0_HNFS.Season <- lme(ln.HNF_q0 ~ HNF_select, random = ~ 1 | Season, data = HNF_Bac_A)
HNFq0_HNFS.St <- lme(ln.HNF_q0 ~ HNF_select, random = ~ 1 | Station, data = HNF_Bac_A)
AIC(HNFq0_HNFS.0, HNFq0_HNFS.Cr, HNFq0_HNFS.Season, HNFq0_HNFS.St)
summary(HNFq0_HNFS.Cr)
summary(gam(ln.HNF_q0 ~ s(HNF_select, bs = "ts"), data = HNF_Bac_A))

### plotting
p_Bacq0_HNFSelect <- HNF_Bac_A %>% 
  select(ln.Bac_q0, ln.HNF_q0, Bac_select, HNF_select, Cruise) %>%
  ggplot(aes(x = ln.Bac_q0, y = HNF_select)) + 
  geom_point(aes(color = Cruise), size = 3) + 
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "dotted") + 
  geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "dotted") + 
  scale_colour_viridis(alpha = 0.7, discrete=TRUE) + 
  labs(x = expression("Log[ Bacteria species richness (Hill number = 0) ]"),
       y = expression("Deterministic assembly processes ( \U03B2NTI) of HNF community ")) + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))

p_HNFSelect_HNFq0 <- HNF_Bac_A %>% 
  select(ln.Bac_q0, ln.HNF_q0, Bac_select, HNF_select, Cruise) %>%
  ggplot(aes(x = HNF_select, y = ln.HNF_q0)) + 
  geom_point(aes(color = Cruise), size = 3) + 
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "dotted") + 
  geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "solid") + 
  scale_colour_viridis(alpha = 0.7, discrete=TRUE) + 
  labs(x = expression("Deterministic assembly processes ( \U03B2NTI) of HNF community "),
       y = expression("Log[ HNF species richness (Hill number = 0) ]")) + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))

legend <- get_legend(
  p_Bacq0_HNFSelect + theme(legend.box.margin = margin(0, 0, 0, 12))
)
p_Select_Bacq0_HNFq0 <- plot_grid(
  plot_grid(p_Bacq0_HNFSelect + theme(legend.position="none"),
            p_HNFSelect_HNFq0 + theme(legend.position="none")),
  legend, rel_widths = c(3, .4))
p_Select_Bacq0_HNFq0
# ggsave(p_Select_Bacq0_HNFq0, file = "D:/Research/PdPy_Div/Presentation/20200110_NTOU/p_Bacq0_HNFSelect_HNFq0.png",
#        dpi = 600, width = 56, height = 28, units = "cm")
##### Bacq0 -> Bac selection ->  HNFq0 ##########

##### Path model analysis : Bac_q0 vs HNF_q0 ##########
##### Step 1: no random effects #####
Bacq0_HNFq0_mod1.0 <- '
  # regressions
    ln.Bac_q0 ~ Bac_select
    ln.HNF_q0 ~ HNF_select
    Bac_select ~ ln.HNF_q0
    HNF_select ~ ln.Bac_q0
  
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q0
'

Bacq0_HNFq0_mod1.1 <- '
  # regressions
    ln.Bac_q0 ~ Bac_select
    ln.HNF_q0 ~ HNF_select
    Bac_select ~ ln.HNF_q0
    HNF_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0 + ln.HNF_q0
    ln.HNF_Biom ~ ln.HNF_q0 + ln.Bac_q0
  
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q0
    ln.Bac_Biom ~~ ln.HNF_Biom
'

Bacq0_HNFq0_mod1.2 <- '
  # regressions
    ln.Bac_q0 ~ Bac_select
    ln.HNF_q0 ~ HNF_select + ln.Bac_Biom
    Bac_select ~ ln.HNF_q0
    HNF_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0
    ln.HNF_Biom ~ ln.HNF_q0 + ln.Bac_q0
  
  # correlated error  
    ln.Bac_q0 ~~ ln.HNF_q0
    ln.Bac_Biom ~~ ln.HNF_Biom
'

Bacq0_HNFq0_mod1.3 <- '
  # regressions
    ln.Bac_q0 ~ Bac_select + ln.HNF_Biom
    ln.HNF_q0 ~ HNF_select + ln.Bac_Biom
    Bac_select ~ ln.HNF_q0
    HNF_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0
    ln.HNF_Biom ~ ln.HNF_q0
  
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q0
    ln.Bac_Biom ~~ ln.HNF_Biom
'

Bacq0_HNFq0_mod1.4 <- '
  # regressions
    ln.Bac_q0 ~ Bac_select + ln.HNF_Biom
    ln.HNF_q0 ~ HNF_select
    Bac_select ~ ln.HNF_q0
    HNF_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0 + ln.HNF_q0
    ln.HNF_Biom ~ ln.HNF_q0
  
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q0
    ln.Bac_Biom ~~ ln.HNF_Biom
'

Bacq0_HNFq0_mod1.5 <- '
  # regressions
    ln.Bac_q0 ~ Bac_select
    ln.HNF_q0 ~ HNF_select
    Bac_select ~ ln.HNF_q0
    HNF_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0 + ln.HNF_q0
    ln.HNF_Biom ~ ln.HNF_q0 + ln.Bac_q0 + ln.Bac_Biom
    
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q0
'

Bacq0_HNFq0_mod1.6 <- '
  # regressions
    ln.Bac_q0 ~ Bac_select
    ln.HNF_q0 ~ HNF_select + ln.Bac_Biom
    Bac_select ~ ln.HNF_q0
    HNF_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0
    ln.HNF_Biom ~ ln.HNF_q0 + ln.Bac_q0 + ln.Bac_Biom
    
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q0
'

Bacq0_HNFq0_mod1.7 <- '
  # regressions
    ln.Bac_q0 ~ Bac_select + ln.HNF_Biom
    ln.HNF_q0 ~ HNF_select + ln.Bac_Biom
    Bac_select ~ ln.HNF_q0
    HNF_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0
    ln.HNF_Biom ~ ln.HNF_q0 + ln.Bac_Biom
    
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q0
'

Bacq0_HNFq0_mod1.8 <- '
  # regressions
    ln.Bac_q0 ~ Bac_select + ln.HNF_Biom
    ln.HNF_q0 ~ HNF_select
    Bac_select ~ ln.HNF_q0
    HNF_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0 + ln.HNF_q0
    ln.HNF_Biom ~ ln.HNF_q0 + ln.Bac_Biom
    
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q0
'

Bacq0_HNFq0_mod1.9 <- '
  # regressions
    ln.Bac_q0 ~ Bac_select
    ln.HNF_q0 ~ HNF_select
    Bac_select ~ ln.HNF_q0
    HNF_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0 + ln.HNF_q0 + ln.HNF_Biom
    ln.HNF_Biom ~ ln.HNF_q0 + ln.Bac_q0
    
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q0
'

Bacq0_HNFq0_mod1.10 <- '
  # regressions
    ln.Bac_q0 ~ Bac_select
    ln.HNF_q0 ~ HNF_select + ln.Bac_Biom
    Bac_select ~ ln.HNF_q0
    HNF_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0 + ln.HNF_Biom
    ln.HNF_Biom ~ ln.HNF_q0 + ln.Bac_q0
    
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q0
'

Bacq0_HNFq0_mod1.11 <- '
  # regressions
    ln.Bac_q0 ~ Bac_select + ln.HNF_Biom
    ln.HNF_q0 ~ HNF_select + ln.Bac_Biom
    Bac_select ~ ln.HNF_q0
    HNF_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0 + ln.HNF_Biom
    ln.HNF_Biom ~ ln.HNF_q0
    
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q0
'

Bacq0_HNFq0_mod1.12 <- '
  # regressions
    ln.Bac_q0 ~ Bac_select + ln.HNF_Biom
    ln.HNF_q0 ~ HNF_select
    Bac_select ~ ln.HNF_q0
    HNF_select ~ ln.Bac_q0
    ln.Bac_Biom ~ ln.Bac_q0 + ln.HNF_q0 + ln.HNF_Biom
    ln.HNF_Biom ~ ln.HNF_q0
    
  # correlated error
    ln.Bac_q0 ~~ ln.HNF_q0
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

AICstep1 <- AIC(Bacq0_HNFq0_lavaan1.0, Bacq0_HNFq0_lavaan1.1, Bacq0_HNFq0_lavaan1.2, Bacq0_HNFq0_lavaan1.3,
                Bacq0_HNFq0_lavaan1.4, Bacq0_HNFq0_lavaan1.5, Bacq0_HNFq0_lavaan1.6, Bacq0_HNFq0_lavaan1.7,
                Bacq0_HNFq0_lavaan1.8, Bacq0_HNFq0_lavaan1.9, Bacq0_HNFq0_lavaan1.10, Bacq0_HNFq0_lavaan1.11,
                Bacq0_HNFq0_lavaan1.12)
AICstep1 <- AICstep1 %>% cbind(row.names(AICstep1)) %>%
  arrange(AIC)
AICstep1

Bacq0_HNFq0_mod1.11 <- '
  ln.Bac_q0 ~ Bac_select + ln.HNF_Biom
  ln.HNF_q0 ~ HNF_select + ln.Bac_Biom
  Bac_select ~ ln.HNF_q0
  HNF_select ~ ln.Bac_q0
  ln.Bac_Biom ~ ln.Bac_q0 + ln.HNF_Biom
  ln.HNF_Biom ~ ln.HNF_q0
  
  ln.Bac_q0 ~~ ln.HNF_q0
'
Bacq0_HNFq0_lavaan1.11 <- sem(Bacq0_HNFq0_mod1.11, data = HNF_Bac_A)
summary(Bacq0_HNFq0_lavaan1.11, fit.measures = TRUE)
##### Step 1: no random effects #####
##### Step 2 : include grouping variables (random effects) and environmental variables #####
Bacq0_HNFq0_psem2.0 <- psem(
  lm(ln.Bac_q0 ~ Bac_select + ln.HNF_Biom, # + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      data = HNF_Bac_A),
  lm(ln.HNF_q0 ~ HNF_select + ln.Bac_Biom, # + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      data = HNF_Bac_A),
  lm(Bac_select ~ ln.HNF_q0+ ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
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
  lme(ln.Bac_q0 ~ Bac_select + ln.HNF_Biom, # + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Cruise, data = HNF_Bac_A),
  lme(ln.HNF_q0 ~ HNF_select + ln.Bac_Biom, # + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
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
  lme(ln.Bac_q0 ~ Bac_select + ln.HNF_Biom, #  + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Season, data = HNF_Bac_A),
  lme(ln.HNF_q0 ~ HNF_select + ln.Bac_Biom, #  + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
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
  lme(ln.Bac_q0 ~ Bac_select + ln.HNF_Biom, #  + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
      random = ~ 1 | Station, data = HNF_Bac_A),
  lme(ln.HNF_q0 ~ HNF_select + ln.Bac_Biom, #  + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3 + ln.Chla, 
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

anova(Bacq0_HNFq0_psem2.0, Bacq0_HNFq0_psem2.Cr, Bacq0_HNFq0_psem2.Season, Bacq0_HNFq0_psem2.St)
anova(Bacq0_HNFq0_psem2.0, Bacq0_HNFq0_psem2.Cr)
summary(Bacq0_HNFq0_psem2.Cr)
##### Step 2 : include grouping variables (random effects) and environmental variables #####
##### Path model analysis : Bac_q0 vs HNF_q0 ##########
###############################################################################################
##### Testing hypothesis: HNFq0 -> Bac selection -> Bacq0 & Bacq0 -> Bac selection ->  HNFq0 ##
###############################################################################################
