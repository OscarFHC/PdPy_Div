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
##### Loading nulls and alpha diversity #######################################################
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

Bac_A <- iNEXT(t(Bac_comm), q = 0, datatype = "abundance", size = max(colSums(Bac_comm)) + 100000)$AsyEst %>% 
  select(Site, Diversity, Estimator) %>% 
  spread(Diversity, Estimator) %>%
  rename(Bac_SR = "Species richness", Bac_Shannon = "Shannon diversity", Bac_Simpson = "Simpson diversity") %>%
  mutate(Site = rownames(Bac_comm))

HNF_MNTD_null <- read.table(file = "D:/Research/PdPy_Div_Results/HNF_MNTD_null.csv", sep = ",", 
                            header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
HNF_MNTD <- HNF_MNTD_null %>% select(c(obs, Var1, Var2)) %>%
  mutate(MNTD_null_mean = apply(HNF_MNTD_null[, !names(HNF_MNTD_null) %in% c("obs", "Var1", "Var2")], 1, mean),
         MNTD_null_sd = apply(HNF_MNTD_null[, !names(HNF_MNTD_null) %in% c("obs", "Var1", "Var2")], 1, sd),
         HNF_select_strength = (obs - MNTD_null_mean) / MNTD_null_sd,
         HNF_select_p = pnorm(-abs(HNF_select_strength), 0, 1))
HNF_Chao_null <- read.table(file = "D:/Research/PdPy_Div_Results/HNF_Chao_null.csv", sep = ",", 
                            header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
HNF_BDiv_Chao <- HNF_Chao_null %>% 
  select(c(obs, Var1, Var2)) %>%
  mutate(Chao_null_mean = apply(HNF_Chao_null[, !names(HNF_Chao_null) %in% c("obs", "Var1", "Var2")], 1, mean),
         Chao_null_sd = apply(HNF_Chao_null[, !names(HNF_Chao_null) %in% c("obs", "Var1", "Var2")], 1, sd),
         HNF_disp_strength = (obs - Chao_null_mean) / Chao_null_sd,
         HNF_disp_p = pnorm(HNF_disp_strength, 0, 1))

HNF_A <- iNEXT(t(HNF_comm), q = 0, datatype = "abundance", size = max(colSums(HNF_comm)) + 100000)$AsyEst %>% 
  select(Site, Diversity, Estimator) %>% 
  spread(Diversity, Estimator) %>%
  rename(HNF_SR = "Species richness", HNF_Shannon = "Shannon diversity", HNF_Simpson = "Simpson diversity") %>%
  mutate(Site = rownames(HNF_comm))
###############################################################################################
##### Loading nulls and alpha diversity #######################################################
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
  filter(!is.na(NF_Biom))
head(HNF_Bac_A)
##### Preping data ##########

##### Exploratory Analyses ##########
### exploratory factor analyses on environmental data
p_Envi_pairs <- HNF_Bac_A %>%
  ggpairs(columns = c("Temp", "Sal", "PAR", "NO2", "NO3", "DIN", "PO3"),
          columnLabels = c("Temperature", "Salinity", "PAR", "Nitrite", "Nitrate", "TN", "TP"),
          upper = list(continuous = cor_fun),
          lower = list(continuous = fit_fun)) +
  theme(strip.text.x = element_text(color = "black", size = 14),
        strip.text.y = element_text(angle = 45, color = "black", size = 14))
p_Envi_pairs
#ggsave(p_Envi_pairs, file = "D:/Research/PdPy_Div_Results/p_Envi_pairs.jpeg", dpi = 600, width = 34, height = 28, units = "cm")

Envi <- HNF_Bac_A[, c("Temp", "Sal", "PAR", "NO2", "NO3", "PO3")] #, "DIN"
fa <- fa.parallel(Envi, fm = "ml", fa = 'fa')

fa1 <- fa(Envi, nfactors = 1, rotate = "varimax", fm = "ml")
print(fa1)
fa2 <- fa(Envi, nfactors = 2, rotate = "oblimin", fm = "ml")
print(fa2)
fa3 <- fa(Envi, nfactors = 3, rotate = "oblimin", fm = "ml")
print(fa3)
fa4 <- fa(Envi, nfactors = 4, rotate = "oblimin", fm = "ml")
print(fa4)
fa5 <- fa(Envi, nfactors = 5, rotate = "oblimin", fm = "ml")
print(fa5)
fa.diagram(fa10)
### plotting


### Bio variables 
p_Adiv_pairs <- HNF_Bac_A %>%
  mutate(ln_Bac_Biom = log(Bac_Biom),
         ln_HNF_Biom = log(HNF_Biom)) %>%
  ggpairs(columns = c("Bac_Shannon", "HNF_Shannon", "ln_Bac_Biom", "ln_HNF_Biom", "Bac_select", "HNF_select"),
          columnLabels = c("Bacteria\nShannon diversity", "HNF\nShannon diversity", "log(Bacteria\nbiomass)", "log(HNF\nbiomass)", 
                           "Bacteria\nselection", "HNF\nselection"),
          #mapping = ggplot2::aes(colour = Cruise),
          upper = list(continuous = cor_fun),
          lower = list(continuous = fit_fun)) +
  theme(strip.text.x = element_text(color = "black", size = 14),
        strip.text.y = element_text(angle = 45, color = "black", size = 14))
p_Adiv_pairs
ggsave(p_Adiv_pairs, file = "D:/Research/PdPy_Div_Results/p_ADiv_pairs.jpeg", dpi = 600, width = 34, height = 28, units = "cm")
##### Exploratory Analyses ##########

##### Plotting ##########
p_ADiv_BacSelect <- HNF_Bac_A %>% 
  ggplot() + 
    geom_point(aes(x = Bac_Shannon, y = HNF_Shannon, color = Bac_select), size = 3) + 
    scale_colour_viridis(alpha = 0.7)
p_ADiv_BacSelect
ggsave(p_ADiv_BacSelect, file = "D:/Research/PdPy_Div_Results/p_ADiv_BacSelect.jpeg")
p_ADiv_HNFSelect <- HNF_Bac_A %>% 
  ggplot() + 
    geom_point(aes(x = Bac_Shannon, y = HNF_Shannon, color = HNF_select), size = 3) + 
    scale_colour_viridis(alpha = 0.7)
p_ADiv_HNFSelect
ggsave(p_ADiv_HNFSelect, file = "D:/Research/PdPy_Div_Results/p_ADiv_HNFSelect.jpeg")
##### Plotting ##########

##### Analyzing ##########
### linear model
lm0_ADiv_Shannon <- lm(Bac_Shannon ~ HNF_Shannon, data = HNF_Bac_A)
lm0_ADiv_Shannon_Sea <- lme(Bac_Shannon ~ HNF_Shannon, random = ~1 | Season, data = HNF_Bac_A)
lm0_ADiv_Shannon_Cr <- lme(Bac_Shannon ~ HNF_Shannon, random = ~1 | Cruise, data = HNF_Bac_A)
lm0_ADiv_Shannon_St <- lme(Bac_Shannon ~ HNF_Shannon, random = ~1 | Station, data = HNF_Bac_A)
AIC(lm0_ADiv_Shannon, lm0_ADiv_Shannon_Sea, lm0_ADiv_Shannon_Cr, lm0_ADiv_Shannon_St)
summary(lm0_ADiv_Shannon_Cr)
lm_ADiv_Shannon_Cr

# lm21_ADiv_Shannon <- lmodel2(Bac_Shannon ~ HNF_Shannon,# + HNF_BDiv_Chao:Bac_select_strength + HNF_BDiv_Chao:HNF_select_strength, 
#                   data = HNF_Bac_A)
# lm21_ADiv_Shannon
lm1_ADiv_Shannon <- lm(Bac_Shannon ~ HNF_Shannon*Bac_select + HNF_Shannon:HNF_select, data = HNF_Bac_A)
lm1_ADiv_Shannon_Sea <- lme(Bac_Shannon ~ HNF_Shannon*Bac_select + HNF_Shannon:HNF_select, random = ~1 | Season, data = HNF_Bac_A)
lm1_ADiv_Shannon_Cr <- lme(Bac_Shannon ~ HNF_Shannon*Bac_select + HNF_Shannon:HNF_select, random = ~1 | Cruise, data = HNF_Bac_A)
lm1_ADiv_Shannon_St <- lme(Bac_Shannon ~ HNF_Shannon*Bac_select + HNF_Shannon:HNF_select, random = ~1 | Station, data = HNF_Bac_A)
AIC(lm1_ADiv_Shannon, lm1_ADiv_Shannon_Sea, lm1_ADiv_Shannon_Cr, lm1_ADiv_Shannon_St)
summary(lm1_ADiv_Shannon_Sea)
summary(lm1_ADiv_Shannon_Cr) # this is the best model
summary(lm1_ADiv_Shannon_St)

HNF_Bac_A <- as.data.frame(HNF_Bac_A)
### SEM
psem0_ADiv <- psem(
  # lme(Bac_Shannon ~ HNF_Shannon*Bac_select + HNF_Shannon:HNF_select, random = ~1 | as.factor(Cruise), data = HNF_Bac_A),
  lm(Bac_Shannon ~ HNF_Shannon + Bac_select + HNF_select, data = HNF_Bac_A),
  # glm(Bac_Biom ~ Bac_Shannon + HNF_Shannon, data = HNF_Bac_A),
  # lme(Bac_Shannon ~ HNF_Shannon*Bac_select + HNF_Shannon:HNF_select, random = ~1 | Season, data = HNF_Bac_A),
  # lme(Bac_Biom ~ Bac_Shannon + HNF_Shannon, random = ~1 | Season, data = HNF_Bac_A),
  # lme(HNF_Biom ~ Bac_Shannon + HNF_Shannon, random = ~1 | Season, data = HNF_Bac_A),
  data = HNF_Bac_A
)
summary(psem0_ADiv)

lavaan_reg <- '
  # regressions
    Bac_Shannon ~ HNF_Shannon + Bac_select + HNF_select
'

lavaan_latent <- '
  # measurement model
    Bac_Shannon =~ Bac_select + HNF_select
  # regressions
    Bac_Shannon ~ HNF_Shannon
'

fit_reg <- sem(lavaan_reg, data = HNF_Bac_A)
fit_lat <- sem(lavaan_latent, data = HNF_Bac_A)
inspect(fit_reg, "r2")
summary(fit_reg)
summary(fit_lat)
Bac_bf <- bf(Bac_Shannon ~ HNF_Shannon*Bac_select + HNF_Shannon:HNF_select + 1 | Season)
HNF_bf <- bf(Bac_Biom ~ Bac_Shannon + HNF_Shannon + 1 | Season)
Mod0 <- brm(formula = Bac_bf + HNF_bf, data = HNF_Bac_A, family = gaussian(), control = list(adapt_delta = 0.9))



##### Analyzing ##########
###############################################################################################
##### Alpha level analyses ####################################################################
###############################################################################################

###############################################################################################
##### Beta level analyses #####################################################################
###############################################################################################
##### Preping data ##########
HNF_Bac_B <- Bac_MNTD %>%
  inner_join(Bac_BDiv_Chao, by = c("Var2" = "Var2", "Var1" = "Var1")) %>%
  rename(Bac_BDiv_MNTD = "obs.x", Bac_BDiv_Chao = "obs.y") %>%
  inner_join(HNF_MNTD, by = c("Var2" = "Var2", "Var1" = "Var1")) %>%
  rename(HNF_BDiv_MNTD = "obs") %>%
  inner_join(HNF_BDiv_Chao, by = c("Var2" = "Var2", "Var1" = "Var1")) %>%
  rename(HNF_BDiv_Chao = "obs") %>%
  inner_join(Vars, by = c("Var2" = "SampleID")) %>%
  filter(Bac_BDiv_MNTD != 0 & HNF_BDiv_MNTD != 0) %>% 
  filter(!is.na(NF_Biom))

#apply(is.na(HNF_Bac_B), 2, which) 

head(HNF_Bac_B)
##### Preping data ##########

##### Exploratory Analyses ##########
p_Bdiv_pairs <- HNF_Bac_B %>%
  mutate(ln_Bac_Biom = log(Bac_Biom),
         ln_HNF_Biom = log(HNF_Biom)) %>%
  ggpairs(columns = c("Bac_BDiv_Chao", "HNF_BDiv_Chao", "ln_Bac_Biom", "ln_HNF_Biom", "Bac_select_strength", "HNF_select_strength"),
          columnLabels = c("Bacteria\nbeta diversity", "HNF\nbeta diversity", "log(Bacteria\nbiomass)", "log(HNF\nbiomass)", 
                           "Bacteria\nselection", "HNF\nselection"),
          #mapping = ggplot2::aes(colour = Cruise),
          upper = list(continuous = cor_fun),
          lower = list(continuous = fit_fun)) +
  theme(strip.text.x = element_text(color = "black", size = 14),
        strip.text.y = element_text(angle = 45, color = "black", size = 14))
p_Bdiv_pairs
ggsave(p_Bdiv_pairs, file = "D:/Research/PdPy_Div_Results/p_BDiv_pairs.jpeg", dpi = 600, width = 34, height = 28, units = "cm")
##### Exploratory Analyses ##########


##### Plotting ##########
p_selec <- HNF_Bac_B %>%
  filter(Bac_select_p <0.05 & HNF_select_p<0.05) %>%
  ggplot() + 
  geom_point(aes(x = Bac_select_strength, y = HNF_select_strength))
p_selec

p_MNTD_Bac_selec <- HNF_Bac_B %>%
  ggplot() + 
  geom_point(aes(x = Bac_BDiv_MNTD, y = HNF_BDiv_MNTD, colour = Bac_select_strength)) + 
  scale_colour_viridis(alpha = 0.7)
p_MNTD_Bac_selec

p_MNTD_HNF_selec <- HNF_Bac_B %>%
  ggplot() + 
  geom_point(aes(x = Bac_BDiv_MNTD, y = HNF_BDiv_MNTD, colour = Bac_select_strength)) + 
  scale_colour_viridis(alpha = 0.7)
p_MNTD_HNF_selec

p_Chao_Bac_selec <- HNF_Bac_B %>%
  ggplot() + 
  geom_point(aes(x = Bac_BDiv_Chao, y = HNF_BDiv_Chao, colour = Bac_select_strength)) + 
  scale_colour_viridis(alpha = 0.7)
p_Chao_Bac_selec
ggsave(p_Chao_Bac_selec, file = "D:/Research/PdPy_Div_Results/p_BDiv_BacSelect.jpeg")

p_Chao_HNF_selec <- HNF_Bac_B %>%
  ggplot() + 
  geom_point(aes(x = Bac_BDiv_Chao, y = HNF_BDiv_Chao, colour = HNF_select_strength)) + 
  scale_colour_viridis(alpha = 0.7)
p_Chao_HNF_selec
ggsave(p_Chao_HNF_selec, file = "D:/Research/PdPy_Div_Results/p_BDiv_HNFSelect.jpeg")
##### Plotting ##########

##### Analyzing ##########
### Major axis linear model #####
mod_MA <- lmodel2(Bac_BDiv_Chao ~ HNF_BDiv_Chao,# + HNF_BDiv_Chao:Bac_select_strength + HNF_BDiv_Chao:HNF_select_strength, 
                  data = HNF_Bac_B, nperm = 5000)
mod_MA
### Major axis linear model #####
### linear model #####
lm1_BDiv_Season <- lme(Bac_BDiv_Chao ~ HNF_BDiv_Chao + HNF_BDiv_Chao:Bac_select_strength + HNF_BDiv_Chao:HNF_select_strength,
                       random = ~ 1 | Season,
                       data = HNF_Bac_B)
lm1_BDiv_Cr <- lme(Bac_BDiv_Chao ~ HNF_BDiv_Chao + HNF_BDiv_Chao:Bac_select_strength + HNF_BDiv_Chao:HNF_select_strength,
                   random = ~ 1 | Cruise,
                   data = HNF_Bac_B)
lm1_BDiv_St <- lme(Bac_BDiv_Chao ~ HNF_BDiv_Chao + HNF_BDiv_Chao:Bac_select_strength + HNF_BDiv_Chao:HNF_select_strength,
                   random = ~ 1 | Station,
                   data = HNF_Bac_B)
AIC(lm1_BDiv_Season, lm1_BDiv_Cr, lm1_BDiv_St)
summary(lm1_BDiv_Season)

### linear model #####



mod_0_psem <- psem(
  #Bac_BDiv_Chao %~~% HNF_BDiv_Chao,
  lme(Bac_BDiv_Chao ~ HNF_BDiv_Chao + HNF_BDiv_Chao:Bac_select_strength + HNF_BDiv_Chao:HNF_select_strength,
      random = ~ 1 | Season,
      data = HNF_Bac_B),
  # lme(HNF_BDiv_Chao ~ Bac_BDiv_Chao*Bac_select_strength + Bac_BDiv_Chao:HNF_select_strength,
  #     random = ~ 1 | Season,
  #     data = HNF_Bac_B),
  data = HNF_Bac_B
)
summary(mod_0_psem)

mod_1_psem <- psem(
  lme(Bac_BDiv_Chao ~ HNF_BDiv_Chao + HNF_BDiv_Chao:Bac_select_strength + HNF_BDiv_Chao:HNF_select_strength,
      random = ~ 1 | Season,
      data = HNF_Bac_B),
  lme(log(Bac_Biom) ~ Bac_BDiv_Chao + HNF_BDiv_Chao,
      random = ~ 1 | Season,
      data = HNF_Bac_B),
  lme(log(HNF_Biom) ~ Bac_BDiv_Chao + HNF_BDiv_Chao,
      random = ~ 1 | Season,
      data = HNF_Bac_B),
  data = HNF_Bac_B
)
summary(mod_1_psem)
coef <- coefs(mod_1_psem, standardize = "scale")
coef[,10] <- coef$Estimate - coef$Std.Error * qnorm(0.975)
coef[,11] <- coef$Estimate + coef$Std.Error * qnorm(0.975)



summary(lme(Bac_BDiv_MNTD ~ HNF_BDiv_MNTD + 
                            Bac_select_strength + HNF_select_strength + Temp + Sal + PAR + DIN + PO3,
            random = ~ Bac_select_strength + HNF_select_strength + Temp + Sal + PAR + DIN + PO3 | Season,
            data = HNF_Bac_B))


dat <- data.frame(x1 = runif(50), y1 = runif(50), y2 = runif(50), y3 = runif(50))
model <- psem(lm(y1 ~ x1 + y2, dat), lm(y2 ~ x1, dat), lm(y3 ~ y1, dat))
summary(model)

Bac_HNF <- bf(Bac_BDiv_MNTD ~~ HNF_BDiv_MNTD)
Bac_bf <- bf(Bac_BDiv_MNTD ~ HNF_BDiv_MNTD + 
     Bac_select_strength + HNF_select_strength + 
     Temp + Sal + PAR + DIN + PO3)
HNF_bf <- bf(HNF_BDiv_MNTD ~ 
     Bac_select_strength + HNF_select_strength +
     Temp + Sal + PAR + DIN + PO3)
Mod0 <- brm(formula = Bac_bf + HNF_bf, data = HNF_Bac_B, family = gaussian(), control = list(adapt_delta = 0.9))

### Analyzing ##########
###############################################################################################
##### Beta level analyses #####################################################################
###############################################################################################





##### Beta level analyses #########

##### Beta level analyses #########




HNF_Bac %>% 
  #filter(HNF_select_p < 0.05 & Bac_select_p < 0.05 ) %>%
  ggplot() + 
  geom_point(aes(x = HNF_select, y = Bac_select))
    

HNF_Bac %>% 
  ggplot() + 
    geom_point(aes(x = Bac_SR, y = HNF_SR, colour = HNF_select)) + 
    scale_colour_viridis()

fn0 <- bf(Bac_SR ~ HNF_SR + HNF_SR:Bac_select + HNF_SR:HNF_select)
fn1_season <- bf(Bac_SR ~ HNF_SR + HNF_SR:Bac_select + (1 + HNF_SR + HNF_SR:Bac_select | Season))
fn2_station <- bf(Bac_SR ~ HNF_SR + HNF_SR:Bac_select + (1 + HNF_SR + HNF_SR:Bac_select | Station))
Mod0 <- brm(formula = fn0, data = HNF_Bac, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0 <- brm(formula = fn0, data = HNF_Bac, family = gaussian(), control = list(adapt_delta = 0.9), autocor = cor_s)
Mod1_season <- brm(formula = fn1_season, data = Bac, family = gaussian(), control = list(adapt_delta = 0.9), iter = 5000, cores = 4)
Mod2_station <- brm(formula = fn2_station, data = Bac, family = gaussian(), control = list(adapt_delta = 0.9), iter = 5000, cores = 4)
summary(Mod1)
marginal_effects(Mod1)


summary(lm(Bac_SR ~ HNF_SR, data = Bac))
summary(lm(Bac_SR ~ HNF_SR*Bac_select, data = Bac))
lme(Bac_SR ~ HNF_SR, random = ~ HNF_SR | Bac_select, data = Bac)