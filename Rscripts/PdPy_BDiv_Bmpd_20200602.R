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

if (!require(sjstats)) {
  install.packages("sjstats", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(sjstats)
}else{library(sjstats)}
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
  file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_4/sECS_Bac_seqXst_PR2_4.csv", 
  sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
Bac_ra_comm <- Bac_comm / rowSums(Bac_comm)
Bac_phylo<- read.tree(
  file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_4/sECS_Bac_treeNJ_PR2_4.tree")

HNF_comm <- as.data.frame(t(read.table(
  file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_4/sECS_HNF_seqXst_PR2_4.csv", 
  sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
HNF_ra_comm <- HNF_comm / rowSums(HNF_comm)
HNF_phylo<- read.tree(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_4/sECS_HNF_treeNJ_PR2_4.tree")

Vars <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_Vars.csv", sep = ",", 
                   header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
###############################################################################################
##### Loading data ############################################################################
###############################################################################################

###############################################################################################
##### Loading nulls ###########################################################################
###############################################################################################
### Bacteria
Bac_Bmpd_null <- read.table(file = "D:/Research/PdPy_Div_Results/nulls_PR2_4/Bac_Bmpd_null_4.csv", 
                            sep = ",", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
Bac_Bmpd <- Bac_Bmpd_null %>% 
  select(c(obs, Var1, Var2))  %>%
  mutate(Bmpd_null_mean = apply(Bac_Bmpd_null[, !names(Bac_Bmpd_null) %in% c("obs", "Var1", "Var2")], 1, mean),
         Bmpd_null_sd = apply(Bac_Bmpd_null[, !names(Bac_Bmpd_null) %in% c("obs", "Var1", "Var2")], 1, sd),
         Bac_Bmpti = (obs - Bmpd_null_mean) / Bmpd_null_sd,
         Bac_Bmpti_select_p = pnorm(-abs(Bac_Bmpti), 0, 1)) %>%
  mutate(ind1 = rep(seq(1, nrow(Bac_comm)), nrow(Bac_comm)),
         ind2 = rep(seq(1, nrow(Bac_comm)), each = nrow(Bac_comm))) %>%
  filter(ind1 > ind2) %>%
  select(-c(ind1, ind2))

### HNF
HNF_Bmpd_null <- read.table(file = "D:/Research/PdPy_Div_Results/nulls_PR2_4/HNF_Bmpd_null_4.csv", 
                            sep = ",", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
HNF_Bmpd <- HNF_Bmpd_null %>% 
  select(c(obs, Var1, Var2)) %>%
  mutate(Bmpd_null_mean = apply(HNF_Bmpd_null[, !names(HNF_Bmpd_null) %in% c("obs", "Var1", "Var2")], 1, mean),
         Bmpd_null_sd = apply(HNF_Bmpd_null[, !names(HNF_Bmpd_null) %in% c("obs", "Var1", "Var2")], 1, sd),
         HNF_Bmpti = (obs - Bmpd_null_mean) / Bmpd_null_sd,
         HNF_Bmpti_select_p = pnorm(-abs(HNF_Bmpti), 0, 1)) %>%
  mutate(ind1 = rep(seq(1, nrow(HNF_comm)), nrow(HNF_comm)),
         ind2 = rep(seq(1, nrow(HNF_comm)), each = nrow(HNF_comm))) %>%
  filter(ind1 > ind2) %>%
  select(-c(ind1, ind2))
###############################################################################################
##### Loading nulls ###########################################################################
###############################################################################################

###############################################################################################
##### Preping data ############################################################################
###############################################################################################
Bac_B <- expand.grid(row.names(Bac_comm), row.names(Bac_comm)) %>%
  mutate(ind1 = rep(seq(1, nrow(Bac_comm)), nrow(Bac_comm)),
         ind2 = rep(seq(1, nrow(Bac_comm)), each = nrow(Bac_comm))) %>%
  filter(ind1 > ind2) %>%
  mutate(Bac_chao = matrix(vegdist(Bac_comm, method = "chao")),
         Var1 = as.character(Var1),
         Var2 = as.character(Var2)) %>%
  select(-c(ind1, ind2))

HNF_B <- expand.grid(row.names(HNF_comm), row.names(HNF_comm)) %>%
  mutate(ind1 = rep(seq(1, nrow(HNF_comm)), nrow(HNF_comm)),
         ind2 = rep(seq(1, nrow(HNF_comm)), each = nrow(HNF_comm))) %>%
  filter(ind1 > ind2) %>%
  mutate(HNF_chao = matrix(vegdist(HNF_comm, method = "chao")),
         Var1 = as.character(Var1),
         Var2 = as.character(Var2)) %>%
  select(-c(ind1, ind2))
  
BDiv <- Bac_B %>%
  inner_join(Bac_Bmpd, by = c("Var1" = "Var1", "Var2" = "Var2")) %>%
  inner_join(HNF_B, by = c("Var1" = "Var1", "Var2" = "Var2")) %>%
  inner_join(HNF_Bmpd, by = c("Var1" = "Var1", "Var2" = "Var2")) %>%
  select(c(Var1, Var2, Bac_chao, HNF_chao, Bac_Bmpti, HNF_Bmpti)) %>%
  inner_join(Vars, by = c("Var2" = "SampleID")) %>%
  filter(!is.na(NF_Biom)) %>%
  mutate(CruiseY = substr(Var1, 4, 9))
BDiv <- as.data.frame(BDiv)

BDiv_mean <- BDiv %>%
  group_by(CruiseY, Station) %>%
  summarize(Bac_chao_mean = mean(Bac_chao),
            HNF_chao_mean = mean(HNF_chao),
            Bac_Bmpti_mean = mean(Bac_Bmpti),
            HNF_Bmpti_mean = mean(HNF_Bmpti)) %>%
  mutate(SampleID = paste0("Or2", CruiseY, Station)) %>%
  inner_join(Vars, by = c("SampleID" = "SampleID"))  %>%
  filter(!is.na(NF_Biom)) %>%
  mutate(ln.Bac_Biom = log(Bac_Biom),
         ln.HNF_Biom = log(HNF_Biom),
         ln.Temp = log(Temp),
         ln.Sal = log(Sal),
         ln.PAR = log(PAR),
         ln.NO2 = log(NO2 + 0.0001),
         ln.NO3 = log(NO3 + 0.0001),
         ln.DIN = log(DIN + 0.0001),
         ln.PO3 = log(PO3 + 0.0001), 
         ln.Chla = log(Chla + 0.00001)) %>%
  select(-c(Cruise, Station.y)) %>%
  rename(Cruise = CruiseY, Station = Station.x)

head(BDiv_mean)
head(BDiv)
###############################################################################################
##### Preping data ############################################################################
###############################################################################################

###############################################################################################
##### Simple HNF and Bac Chao dissimilarity relationship  #####################################
###############################################################################################

# Check non-linear association between HNF and Bac chao diversity
gam0 <- gam(HNF_chao_mean ~ Bac_chao_mean, data = BDiv_mean)
gam1 <- gam(HNF_chao_mean ~ s(Bac_chao_mean), data = BDiv_mean)
anova(gam0, gam1, test = "F")
summary(gam0)

HNFchao_Bacchao.0 <- lm(HNF_chao_mean ~ Bac_chao_mean, data = BDiv_mean)
HNFchao_Bacchao.St <- lme(HNF_chao_mean ~ Bac_chao_mean, random = ~ 1 | Station, data = BDiv_mean, method = "ML")
#HNFchao_Bacchao.Cr <- lme(HNF_chao_mean ~ Bac_chao_mean, random = ~ 1 | Cruise, data = BDiv_mean, method = "ML")
HNFchao_Bacchao.Season <- lme(HNF_chao_mean ~ Bac_chao_mean, random = ~ 1 | Season, data = BDiv_mean, method = "ML")
AIC(HNFchao_Bacchao.0, HNFchao_Bacchao.St, HNFchao_Bacchao.Season) #, HNFchao_Bacchao.Cr
summary(HNFchao_Bacchao.Season)
performance::r2(HNFchao_Bacchao.Season)

p_HNFchao_Bacchao <- BDiv_mean %>% 
  select(Bac_chao_mean, HNF_chao_mean, Bac_Bmpti_mean, HNF_Bmpti_mean, Cruise, Season) %>%
  ggplot(aes(x = Bac_chao_mean, y = HNF_chao_mean)) + 
    geom_point(size = 5, aes(color = Season)) + 
    scale_colour_viridis(alpha = 1, discrete=TRUE) +
    geom_smooth(formula = y ~ x, method = "lm", se = TRUE) + 
    geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "solid") + 
    #scale_x_continuous(expand = c(0, 0)) +
    #scale_y_continuous(expand = c(0, 0)) +
    labs(x = bquote("Bacteria Chao similarity"),
         y = bquote("HNF Chao similarity")) + 
    annotate("text", x = 0.18, y = 1, label = "paste( \"conditional \", italic(R) ^ 2, \" = 0.18\")", parse = TRUE, size = 6) + 
    annotate("text", x = 0.18, y = 0.96, label = "paste( \"marginal \", italic(R) ^ 2, \" = 0.12\")", parse = TRUE, size = 6) + 
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    )

p_HNFchao_Bacchao
ggsave(p_HNFchao_Bacchao, file = "D:/Research/PdPy_Div_Results/Figs/Lab422 meeting_20200602/BDivCor.png",
       dpi = 600, width = 34, height = 28, units = "cm")
###############################################################################################
##### Simple HNF and Bac Chao dissimilarity relationship  #####################################
###############################################################################################

###############################################################################################
##### Testing hypothesis: HNFChao -> Bac BMPTI -> BacChao & BacChao -> HNF BMPTI -> HNFChao ###
###############################################################################################
##### HNFChao -> Bac BMPTI -> BacChao ##########
### linear model testing
BacS_HNFchao.0 <- lm(Bac_Bmpti_mean ~ HNF_chao_mean + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, data = BDiv_mean)
BacS_HNFchao.St <- lme(Bac_Bmpti_mean ~ HNF_chao_mean + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, random = ~ 1 | Station, data = BDiv_mean, method = "ML")
#BacS_HNFchao.Cr <- lme(Bac_Bmpti_mean ~ HNF_chao_mean + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, random = ~ 1 | Cruise, data = BDiv_mean, method = "ML")
BacS_HNFchao.Season <- lme(Bac_Bmpti_mean ~ HNF_chao_mean + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, random = ~ 1 | Season, data = BDiv_mean, method = "ML")
AIC(BacS_HNFchao.0, BacS_HNFchao.St, BacS_HNFchao.Season) #, BacS_HNFchao.Cr
summary(BacS_HNFchao.Season)
summary(gam(Bac_Bmpti_mean ~ s(HNF_chao_mean, bs = "ts"), data = BDiv_mean))
performance::r2(BacS_HNFchao.Season)

Bacchao_BacS.0 <- lm(Bac_chao_mean ~ Bac_Bmpti_mean + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, data = BDiv_mean)
Bacchao_BacS.St <- lme(Bac_chao_mean ~ Bac_Bmpti_mean + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, random = ~ 1 | Station, data = BDiv_mean, method = "ML")
#Bacchao_BacS.Cr <- lme(Bac_chao_mean ~ Bac_Bmpti_mean + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, random = ~ 1 | Cruise, data = BDiv_mean, method = "ML")
Bacchao_BacS.Season <- lme(Bac_chao_mean ~ Bac_Bmpti_mean + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, random = ~ 1 | Season, data = BDiv_mean, method = "ML")
AIC(Bacchao_BacS.0, Bacchao_BacS.St, Bacchao_BacS.Season) #, BacS_HNFchao.Cr
summary(Bacchao_BacS.Season)
summary(gam(Bac_chao ~ s(Bac_Bmpti, bs = "ts"), data = BDiv))
performance::r2(Bacchao_BacS.Season)

### plotting
p_HNFchao_BacBmpti <- BDiv_mean %>% 
  select(Bac_chao_mean, HNF_chao_mean, Bac_Bmpti_mean, HNF_Bmpti_mean, Cruise, Season) %>%
  ggplot(aes(x = HNF_chao_mean, y = Bac_Bmpti_mean)) + 
    geom_point(aes(color = Season), size = 5) + 
    geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "solid") + 
    geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "solid") + 
    scale_colour_viridis(alpha = 1, discrete=TRUE) + 
    labs(x = bquote(atop("HNF Chao similarity", " ")),
         y = bquote(atop("Deterministic assembly processes (\U03B2MPTI)", "of Bacteria community"))) + 
    annotate("text", x = 0.3, y = -0.45, label = "paste( \"conditional \", italic(R) ^ 2, \" = 0.23\")", parse = TRUE, size = 6) + 
    annotate("text", x = 0.3, y = -0.52, label = "paste( \"marginal \", italic(R) ^ 2, \" = 0.23\")", parse = TRUE, size = 6) +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)
      )
p_BacBmpti_Bacchao <- BDiv_mean %>% 
  select(Bac_chao_mean, HNF_chao_mean, Bac_Bmpti_mean, HNF_Bmpti_mean, Cruise, Season) %>%
  ggplot(aes(x = Bac_Bmpti_mean, y = Bac_chao_mean)) + 
    geom_point(aes(color = Season), size = 5) + 
    geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "solid") + 
    geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "solid") + 
    scale_colour_viridis(alpha = 1, discrete=TRUE) + 
    labs(x = bquote(atop("Deterministic assembly processes (\U03B2MPTI)", "of Bacteria community")),
         y = bquote("Bacterial Chao similarity")) + 
    annotate("text", x = -1.98, y = 0.8, label = "paste( \"conditional \", italic(R) ^ 2, \" = 0.26\")", parse = TRUE, size = 6) + 
    annotate("text", x = -1.98, y = 0.77, label = "paste( \"marginal \", italic(R) ^ 2, \" = 0.26\")", parse = TRUE, size = 6) + 
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    )
legend <- get_legend(
  p_BacBmpti_Bacchao + theme(legend.box.margin = margin(0, 0, 0, 12))
)

p_HNFchao_BacBmpti_Bacchao <- plot_grid(
  plot_grid(p_HNFchao_BacBmpti + theme(legend.position = "none"),
            p_BacBmpti_Bacchao + theme(legend.position = "none"),
            ncol = 2, labels = "AUTO", hjust = 0),
  legend, rel_widths = c(3, .4))
p_HNFchao_BacBmpti_Bacchao
ggsave(p_HNFchao_BacBmpti_Bacchao, file = "D:/Research/PdPy_Div_Results/Figs/Lab422 meeting_20200602/p_HNFchao_BacBmpti_Bacchao.png",
       dpi = 600, width = 64, height = 32, units = "cm")
##### HNFChao -> Bac BMPTI -> BacChao ##########

##### BacChao -> HNF BMPTI -> HNFChao ##########
### linear model testing
HNFBmpti_Bacchao.0 <- lm(HNF_Bmpti_mean ~ Bac_chao_mean + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, data = BDiv_mean)
HNFBmpti_Bacchao.St <- lme(HNF_Bmpti_mean ~ Bac_chao_mean + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, random = ~ 1 | Station, data = BDiv_mean, method = "ML")
#HNFBmpti_Bacchao.Cr <- lme(HNF_Bmpti_mean ~ Bac_chao_mean + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, random = ~ 1 | Cruise, data = BDiv_mean, method = "ML")
HNFBmpti_Bacchao.Season <- lme(HNF_Bmpti_mean ~ Bac_chao_mean + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, random = ~ 1 | Season, data = BDiv_mean, method = "ML")
AIC(HNFBmpti_Bacchao.0, HNFBmpti_Bacchao.St, HNFBmpti_Bacchao.Season) #HNFBmpti_Bacchao.Cr
summary(HNFBmpti_Bacchao.Season)
summary(gam(HNF_Bmpti_mean ~ s(Bac_chao_mean, bs = "ts"), data = BDiv_mean))
performance::r2(HNFBmpti_Bacchao.Season)

HNFchao_HNFBmpti.0 <- lm(HNF_chao_mean ~ HNF_Bmpti_mean + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, data = BDiv_mean)
HNFchao_HNFBmpti.St <- lme(HNF_chao_mean ~ HNF_Bmpti_mean + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, random = ~ 1 | Station, data = BDiv_mean, method = "ML")
#HNFchao_HNFBmpti.Cr <- lme(HNF_chao_mean ~ HNF_Bmpti_mean + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, random = ~ 1 | Cruise, data = BDiv_mean, method = "ML")
HNFchao_HNFBmpti.Season <- lme(HNF_chao_mean ~ HNF_Bmpti_mean + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, random = ~ 1 | Season, data = BDiv_mean, method = "ML")
AIC(HNFchao_HNFBmpti.0, HNFchao_HNFBmpti.St, HNFchao_HNFBmpti.Season) #HNFchao_HNFBmpti.Cr
summary(HNFchao_HNFBmpti.Season)
summary(gam(HNF_chao_mean ~ s(HNF_Bmpti_mean, bs = "ts"), data = BDiv_mean))
performance::r2(HNFchao_HNFBmpti.Season)

### plotting
p_Bacchao_HNFBmpti <- BDiv_mean %>% 
  select(Bac_chao_mean, HNF_chao_mean, Bac_Bmpti_mean, HNF_Bmpti_mean, Cruise, Season) %>%
  ggplot(aes(x = Bac_chao_mean, y = HNF_Bmpti_mean)) + 
    geom_point(aes(color = Season), size = 5) + 
    geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "solid") + 
    geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "solid") + 
    scale_colour_viridis(alpha = 1, discrete=TRUE) + 
    labs(x = bquote(atop("Bacterial Chao similarity", " ")),
         y = bquote(atop("Deterministic assembly processes (\U03B2MPTI)", "of HNF community"))) + 
    annotate("text", x = 0.2, y = -0.6, label = "paste( \"conditional \", italic(R) ^ 2, \" = 0.49\")", parse = TRUE, size = 6) + 
    annotate("text", x = 0.2, y = -0.8, label = "paste( \"marginal \", italic(R) ^ 2, \" = 0.41\")", parse = TRUE, size = 6) +
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    )
p_HNFBmpti_HNFchao <- BDiv_mean %>% 
  select(Bac_chao_mean, HNF_chao_mean, Bac_Bmpti_mean, HNF_Bmpti_mean, Cruise, Season) %>%
  ggplot(aes(x = HNF_Bmpti_mean, y = HNF_chao_mean)) + 
    geom_point(aes(color = Season), size = 5) + 
    geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "solid") + 
    geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "solid") + 
    scale_colour_viridis(alpha = 1, discrete = TRUE) + 
    labs(x = bquote(atop("Deterministic assembly processes (\U03B2MPTI)", "of HNF community")),
         y = bquote("HNF Chao dissimilarity")) + 
    annotate("text", x = -4.65, y = 0.9, label = "paste( \"conditional \", italic(R) ^ 2, \" = 0.26\")", parse = TRUE, size = 6) + 
    annotate("text", x = -4.65, y = 0.87, label = "paste( \"marginal \", italic(R) ^ 2, \" = 0.26\")", parse = TRUE, size = 6) + 
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    )

legend <- get_legend(
  p_HNFBmpti_HNFchao + theme(legend.box.margin = margin(0, 0, 0, 12))
)
p_Bacchao_HNFBmpti_HNFchao <- plot_grid(
  plot_grid(p_Bacchao_HNFBmpti + theme(legend.position = "none"),
            p_HNFBmpti_HNFchao + theme(legend.position = "none"),
            ncol = 2, labels = "AUTO", hjust = 0),
  legend, rel_widths = c(3, .4))
p_Bacchao_HNFBmpti_HNFchao
ggsave(p_Bacchao_HNFBmpti_HNFchao, file = "D:/Research/PdPy_Div_Results/Figs/Lab422 meeting_20200602/p_Bacchao_HNFBmpti_HNFchao.png",
       dpi = 600, width = 64, height = 32, units = "cm")
##### BacChao -> HNF BMPTI -> HNFChao ##########
###############################################################################################
##### Testing hypothesis: HNFChao -> Bac BMPTI -> BacChao & BacChao -> HNF BMPTI -> HNFChao ###
###############################################################################################

###############################################################################################
##### Ampti-> HNFq1_Bacq1 association  ########################################################
###############################################################################################
### statistical test
BDiv_associ <- BDiv_mean %>%
  select(Cruise, Season, HNF_chao_mean, Bac_chao_mean, HNF_Bmpti_mean, Bac_Bmpti_mean) %>%
  group_by(Cruise) %>%
  mutate(
    R2 = summary(lm(HNF_chao_mean ~ Bac_chao_mean))$r.squared,
    coef = summary(lm(HNF_chao_mean ~ Bac_chao_mean))$coefficients[2, 1]
  )  %>%
  summarize(
    Bac_Bmpti_CrMean = mean(Bac_Bmpti_mean),
    HNF_Bmpti_CrMean = mean(HNF_Bmpti_mean),
    R2 = mean(R2),
    coef = mean(coef)
  )
for (i in 1: nrow(BDiv_associ)){
  BDiv_associ[i,"Season"] = BDiv_mean$Season[which(BDiv_associ$Cruise[i] == BDiv_mean$Cruise)][1]
}

### R2
Mod_R2_BacBmpti_CrMean <- lme(R2 ~ Bac_Bmpti_CrMean, random = ~ 1 | Season, data = BDiv_associ)
Mod_R2_HNFBmpti_CrMean <- lme(R2 ~ HNF_Bmpti_CrMean, random = ~ 1 | Season, data = BDiv_associ) #[which(BDiv_associ$coef < 1),]
summary(Mod_R2_BacBmpti_CrMean)
summary(Mod_R2_HNFBmpti_CrMean)

p_BDivR2_BacBmpti <- ggplot(data = BDiv_associ, aes(x = Bac_Bmpti_CrMean, y = R2)) + 
  geom_point(aes(color = Season), size = 5) +
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "dashed") + 
  scale_colour_viridis(alpha = 1, discrete = TRUE) + 
  scale_y_continuous(lim= c(-0.5, 0.9)) + 
  labs(x = expression(atop("Deterministic assembly processes (\U03B2MPTI)", "of Bacteria community")),
       y = bquote("R"^2~ "of diversity association")) + 
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )
p_BDivR2_HNFBmpti <- ggplot(data = BDiv_associ, aes(x = HNF_Bmpti_CrMean, y = R2)) + 
  geom_point(aes(color = Season), size = 5) +
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "dashed") + 
  scale_colour_viridis(alpha = 1, discrete = TRUE) + 
  scale_y_continuous(lim= c(-0.5, 0.9)) + 
  labs(x = expression(atop("Deterministic assembly processes (\U03B2MPTI)", "of HNF community")),
       y = bquote("R"^2~ "of diversity association")) + 
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )
legend <- get_legend(
  p_BDivR2_BacBmpti + theme(legend.box.margin = margin(0, 0, 0, 12))
)
p_BDivR2_Bmpti <- plot_grid(
  plot_grid(p_BDivR2_BacBmpti + theme(legend.position = "none"),
            p_BDivR2_HNFBmpti + theme(legend.position = "none"),
            ncol = 2, labels = "AUTO", hjust = 0),
  legend, rel_widths = c(3, .4))
p_BDivR2_Bmpti
ggsave(p_BDivR2_Bmpti, file = "D:/Research/PdPy_Div_Results/Figs/Lab422 meeting_20200602/p_BDivR2_Bmpti.png",
       dpi = 600, width = 64, height = 32, units = "cm")

### Coef
Mod_coef_BacBmpti_CrMean <- lme(coef ~ Bac_Bmpti_CrMean, random = ~ 1 | Season, data = BDiv_associ)
Mod_coef_HNFBmpti_CrMean <- lme(coef ~ HNF_Bmpti_CrMean, random = ~ 1 | Season, data = BDiv_associ) #[which(BDiv_associ$coef < 1),]
summary(Mod_coef_BacBmpti_CrMean)
summary(Mod_coef_HNFBmpti_CrMean)

p_BDivcoef_BacBmpti <- ggplot(data = BDiv_associ, aes(x = Bac_Bmpti_CrMean, y = coef)) + #[which(Adiv_cr$coef > -1), ]
  geom_point(aes(color = Season), size = 5) +
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "dashed") + 
  scale_colour_viridis(alpha = 1, discrete = TRUE) + 
  scale_y_continuous(lim= c(-0.8, 0.8)) + 
  labs(x = expression(atop("Deterministic assembly processes (\U03B2MPTI)", "of Bacteria community")),
       y = bquote(atop("Regression coefficient", "of diversity association"))) + 
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )
p_BDivcoef_HNFBmpti <- ggplot(data = BDiv_associ, aes(x = HNF_Bmpti_CrMean, y = coef)) + #[which(ADiv_Cr$coef > -1), ]
  geom_point(aes(color = Season), size = 5) +
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "dashed") + 
  scale_colour_viridis(alpha = 1, discrete = TRUE) + 
  scale_y_continuous(lim= c(-0.8, 0.8)) + 
  labs(x = expression(atop("Deterministic assembly processes (\U03B2MPTI)", "of HNF community")),
       y = bquote(atop("Regression coefficient", "of diversity association"))) + 
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )
legend <- get_legend(
  p_BDivcoef_BacBmpti + theme(legend.box.margin = margin(0, 0, 0, 12))
)
p_BDivcoef_Bmpti <- plot_grid(
  plot_grid(p_BDivcoef_BacBmpti + theme(legend.position = "none"),
            p_BDivcoef_HNFBmpti + theme(legend.position = "none"),
            ncol = 2, labels = "AUTO", hjust = 0),
  legend, rel_widths = c(3, .4))
p_BDivcoef_Bmpti
ggsave(p_BDivcoef_Bmpti, file = "D:/Research/PdPy_Div_Results/Figs/Lab422 meeting_20200602/p_BDivcoef_Bmpti.png",
       dpi = 600, width = 64, height = 32, units = "cm")

###############################################################################################
##### Ampti-> HNFq1_Bacq1 association  ########################################################
###############################################################################################

############################## Supplementary analyses ##################################################################
###############################################################################################
##### sampling maps ###########################################################################
###############################################################################################
St_sECS = read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/raw/St_Location.csv", 
                     header = T, fill = T, sep = ",")
#library(ggmap)
if (!require(maps)) {
  install.packages("maps", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(maps)
}else{library(maps)}

world <- map_data("world")
map <- ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group)) +
  geom_point(data = St_sECS[1:6,], aes(x = Lon, y = Lat), size = 5) + 
  coord_fixed(1, xlim = c(119.5, 122.75), ylim = c(24.5, 27.5)) + 
  annotate("text", x = 121.35, y = 24.75, label = "Taiwan", color = "white", size = 7) + 
  geom_text(data = St_sECS[1:5,], aes(x = Lon, y = Lat, label = Station), hjust = -0.3, vjust = -0.2, size = 7) + 
  geom_text(data = St_sECS[6,], aes(x = Lon, y = Lat, label = Station), hjust = 0.4, vjust = 1.7, size = 7) + 
  labs(x = expression("Longitude"),
       y = expression("Latitude")) + 
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16))
map
ggsave(map, file = "D:/Manuscripts/PdPy_Div_MS/ms_Figs/Fig2_map.png", dpi = 600) 
###############################################################################################
##### sampling maps ###########################################################################
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
pairs(BDiv[, c("Bac_chao", "HNF_chao", "Bac_Bmpti", "HNF_Bmpti")])
cor(BDiv[, c("Bac_chao", "HNF_chao", "Bac_Bmpti", "HNF_Bmpti")])
corr.test(BDiv[, c("Bac_chao", "HNF_chao", "Bac_Bmpti", "HNF_Bmpti")])
p_Bdiv_BAssemb_pairs <- BDiv %>%
  ggpairs(columns = c("Bac_chao", "HNF_chao", "Bac_Bmpti", "HNF_Bmpti"),
          columnLabels = c("Bacteria\nchao\ndiversity", "HNF\nchao\ndiversity",
                           "Bacteria\n\U03B2MPTI", "HNF\n\U03B2MPTI"),
          #mapping = ggplot2::aes(colour = Cruise),
          upper = list(continuous = cor_fun),
          lower = list(continuous = fit_fun)) +
  theme(strip.text.x = element_text(color = "black", size = 14),
        strip.text.y = element_text(angle = 45, color = "black", size = 14))
p_Bdiv_BAssemb_pairs
ggsave(p_Adiv_AAssemb_pairs, file = "D:/Research/PdPy_Div_Results/Figs/PR2_3/FigS1_Adiv_AMPTI_pairs.png",
       dpi = 600, width = 34, height = 28, units = "cm")
##### Pair-wise plot of bio-variables ##########
###############################################################################################
##### Envi exploratory factor analyses and  pair-wise correlations ############################
###############################################################################################

###############################################################################################
##### Phylogenetic signal #####################################################################
###############################################################################################
##### Bacteria ##########
Bac_Niche <- as.data.frame(matrix(0, ncol(Bac_comm), 6)) %>%
  #mutate(Bac_ID = names(Bac_comm)) %>%
  rename(Temp = V1, Sal = V2, PAR = V3, DIN = V4, PO3 = V5, Chla = V6) #, NO2 = V4, NO3 = V5

for(i in 1:ncol(Bac_comm)){
  Bac_Niche[i, "Temp"] = sum(Bac_comm[,names(Bac_comm)[i]] * Vars$Temp) / sum(Bac_comm[,names(Bac_comm)[i]])
  Bac_Niche[i, "Sal"] = sum(Bac_comm[,names(Bac_comm)[i]] * Vars$Sal) / sum(Bac_comm[,names(Bac_comm)[i]])
  Bac_Niche[i, "PAR"] = sum(Bac_comm[,names(Bac_comm)[i]] * Vars$PAR) / sum(Bac_comm[,names(Bac_comm)[i]])
  #Bac_Niche[i, "NO2"] = sum(Bac_comm[,names(Bac_comm)[i]] * Vars$NO2) / sum(Bac_comm[,names(Bac_comm)[i]])
  #Bac_Niche[i, "NO3"] = sum(Bac_comm[,names(Bac_comm)[i]] * Vars$NO3) / sum(Bac_comm[,names(Bac_comm)[i]])
  Bac_Niche[i, "DIN"] = sum(Bac_comm[,names(Bac_comm)[i]] * Vars$DIN) / sum(Bac_comm[,names(Bac_comm)[i]])
  Bac_Niche[i, "PO3"] = sum(Bac_comm[,names(Bac_comm)[i]] * Vars$PO3) / sum(Bac_comm[,names(Bac_comm)[i]])
  Bac_Niche[i, "Chla"] = sum(Bac_comm[,names(Bac_comm)[i]] * Vars$Chla) / sum(Bac_comm[,names(Bac_comm)[i]])
}

Bac_Niche <- as.data.frame(scale(Bac_Niche, center = TRUE, scale = TRUE))
Bac_Niche_Dist <- vegdist(Bac_Niche, method = "euclidean")
Bac_Phylo_Dist <- cophenetic(Bac_phylo)
Bac_PhySig <- mantel.correlog(Bac_Niche_Dist, Bac_Phylo_Dist, n.class = 50)
##### Bacteria ##########

##### HNF ##########
HNF_noZeroComm <- HNF_comm[, which(colSums(HNF_comm) != 0)]
HNF_Niche <- as.data.frame(matrix(0, ncol(HNF_noZeroComm), 6)) %>%
  # mutate(HNF_ID = names(HNF_comm)) %>%
  rename(Temp = V1, Sal = V2, PAR = V3, DIN = V4, PO3 = V5, Chla = V6)
Vars_HNF <- Vars[Vars$SampleID %in% row.names(HNF_noZeroComm),]

for(i in 1:ncol(HNF_noZeroComm)){
  HNF_Niche[i, "Temp"] = sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]] * Vars_HNF$Temp) / sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]])
  HNF_Niche[i, "Sal"] = sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]] * Vars_HNF$Sal) / sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]])
  HNF_Niche[i, "PAR"] = sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]] * Vars_HNF$PAR) / sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]])
  #HNF_Niche[i, "NO2"] = sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]] * Vars_HNF$NO2) / sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]])
  #HNF_Niche[i, "NO3"] = sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]] * Vars_HNF$NO3) / sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]])
  HNF_Niche[i, "DIN"] = sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]] * Vars_HNF$DIN) / sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]])
  HNF_Niche[i, "PO3"] = sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]] * Vars_HNF$PO3) / sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]])
  HNF_Niche[i, "Chla"] = sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]] * Vars_HNF$Chla) / sum(HNF_noZeroComm[,names(HNF_noZeroComm)[i]])
}

HNF_Niche <- as.data.frame(scale(HNF_Niche, center = TRUE, scale = TRUE))
HNF_Niche_Dist <- vegdist(HNF_Niche, method = "euclidean")
HNF_Phylo_Dist <- cophenetic(drop.tip(HNF_phylo, which(colSums(HNF_comm) == 0)))
HNF_PhySig <- mantel.correlog(HNF_Niche_Dist, HNF_Phylo_Dist, n.class = 50)
##### HNF ##########

##### Visualization #########
Bac_PhySig_p <- as.data.frame(Bac_PhySig$mantel.res) %>%
  rename(PrCorrected = "Pr(corrected)", PrMantel = "Pr(Mantel)") %>%
  filter(PrCorrected != "NA") %>%
  mutate(sig = ifelse(PrCorrected < 0.05, "significant", "non-signifant")) %>%
  ggplot() + 
  geom_point(aes(x = class.index, y = Mantel.cor, color = as.factor(sig)), size = 3) + 
  geom_hline(yintercept = 0, color = "red") + 
  scale_colour_manual( values = c("grey", "black")) + 
  labs(x = expression("Phylogenetic distance"),
       y = expression("Mantel correlation")
  ) + 
  theme_classic() + 
  theme(
    strip.text.x = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 20),
    legend.title = element_blank()
  )

HNF_PhySig_p <- as.data.frame(HNF_PhySig$mantel.res) %>%
  rename(PrCorrected = "Pr(corrected)", PrMantel = "Pr(Mantel)") %>%
  filter(PrCorrected != "NA") %>%
  mutate(sig = ifelse(PrCorrected < 0.1, "significant", "non-signifant")) %>%
  ggplot() + 
  geom_point(aes(x = class.index, y = Mantel.cor, color = as.factor(sig)), size = 3) + 
  geom_hline(yintercept = 0, color = "red") + 
  scale_colour_manual(values = c("grey", "grey32"), 
                      labels = c(expression(atop("non-", "significant")), expression(atop("marginal", "significant")))) + 
  labs(x = expression("Phylogenetic distance"),
       y = expression("Mantel correlation")
  ) + 
  theme_classic() + 
  theme(
    strip.text.x = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 20),
    legend.title = element_blank()
  )

legend <- get_legend(
  HNF_PhySig_p + theme(legend.box.margin = margin(0, 0, 0, 12))
)
p_PhySig <- plot_grid(
  plot_grid(Bac_PhySig_p + theme(legend.position = "none"),
            HNF_PhySig_p + theme(legend.position = "none"),
            ncol = 1, labels = "AUTO", hjust = 0),
  legend, rel_widths = c(3, .4))
p_PhySig

ggsave(p_PhySig, file = "D:/Manuscript/PdPy_Div_MS/ms_Figs/FigS2_PhySig.png",
       dpi = 600, width = 34, height = 28, units = "cm")
##### Visualization ##########
###############################################################################################
##### Phylogenetic signal #####################################################################
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

