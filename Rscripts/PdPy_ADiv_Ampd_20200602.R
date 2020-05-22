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
Bac_Ampd_null <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/Anulls_PR2_4/Bac_Ampd_null_4.csv", 
                            sep = ",", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
Bac_Ampd <- Bac_Ampd_null %>% 
  select(c(obs, Site)) %>%
  mutate(Ampd_null_mean = apply(Bac_Ampd_null[, !names(Bac_Ampd_null) %in% c("obs", "Site")], 1, mean),
         Ampd_null_sd = apply(Bac_Ampd_null[, !names(Bac_Ampd_null) %in% c("obs", "Site")], 1, sd),
         Bac_Ampti = (obs - Ampd_null_mean) / Ampd_null_sd,
         Bac_Ampti_p = pnorm(-abs(Bac_Ampti), 0, 1))

HNF_Ampd_null <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/Anulls_PR2_4/HNF_Ampd_null_4.csv", 
                            sep = ",", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
HNF_Ampd <- HNF_Ampd_null %>% 
  select(c(obs, Site)) %>%
  mutate(Ampd_null_mean = apply(HNF_Ampd_null[, !names(HNF_Ampd_null) %in% c("obs", "Site")], 1, mean),
         Ampd_null_sd = apply(HNF_Ampd_null[, !names(HNF_Ampd_null) %in% c("obs", "Site")], 1, sd),
         HNF_Ampti = (obs - Ampd_null_mean) / Ampd_null_sd,
         HNF_Ampti_p = pnorm(-abs(HNF_Ampti), 0, 1))
###############################################################################################
##### Loading nulls ###########################################################################
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

HNF_Bac_A <- Bac_A %>%
  inner_join(Bac_Ampd, by = c("Site" = "Site")) %>%
  inner_join(HNF_A, by = c("Site" = "Site")) %>%
  inner_join(HNF_Ampd, by = c("Site" = "Site")) %>%
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
##### Simple HNF and Bac alppha diversity relationship  #######################################
###############################################################################################

# Check non-linear correlation between HNF and Bac Shannon diversity
gam0 <- gam(ln.HNF_q1 ~ ln.Bac_q1, data = HNF_Bac_A)
gam1 <- gam(ln.HNF_q1 ~ s(ln.Bac_q1), data = HNF_Bac_A)
anova(gam0, gam1, test = "F")
summary(gam0)

# plotting correlation between HNF and Bac Shannon diversity
p_HNFq1_Bacq1 <- HNF_Bac_A %>% 
  select(ln.Bac_q1, ln.HNF_q1, Bac_Ampd_select, HNF_Ampd_select, Cruise) %>%
  ggplot(aes(x = ln.Bac_q1, y = ln.HNF_q1)) + 
    geom_point(size = 3, aes(color = Cruise)) + 
    geom_smooth(formula = y ~ x, method = "lm", se = TRUE) + 
    scale_colour_viridis(alpha = 0.7, discrete=TRUE) +
    # geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "dotted") + 
    labs(x = expression("Log[ Bacteria Shannon diversity ]"),
         y = expression("Log[ HNF Shannon diversity ]")) + 
    #theme_classic() + 
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.text.x = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    )
p_HNFq1_Bacq1
ggsave(p_HNFq1_Bacq1, file = "D:/Manuscripts/PdPy_Div_MS/ms_Figs/Fig3_HNFq1_Bacq1.png",
       dpi = 600, width = 34, height = 28, units = "cm")
###############################################################################################
##### Simple HNF and Bac alppha diversity relationship  #######################################
###############################################################################################

###############################################################################################
##### Testing hypothesis: HNFq0 -> Bac selection -> Bacq0 & Bacq0 -> Bac selection ->  HNFq0 ##
###############################################################################################
##### HNFq1 -> Bac selection -> Bacq1 ##########
### linear model testing
BacS_HNFq1.0 <- lm(Bac_Ampti ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, data = HNF_Bac_A)
BacS_HNFq1.Cr <- lme(Bac_Ampti ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, random = ~ 1 | Cruise, data = HNF_Bac_A, method = "ML")
BacS_HNFq1.Season <- lme(Bac_Ampti ~ ln.HNF_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, random = ~ 1 | Season, data = HNF_Bac_A, method = "ML")
AIC(BacS_HNFq1.0, BacS_HNFq1.Cr, BacS_HNFq1.Season)
summary(BacS_HNFq1.Cr)

Bacq1_BacS.0 <- lm(ln.Bac_q1 ~ Bac_Ampti + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, data = HNF_Bac_A)
Bacq1_BacS.Cr <- lme(ln.Bac_q1 ~ Bac_Ampti + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, random = ~ 1 | Cruise, data = HNF_Bac_A, method = "ML")
Bacq1_BacS.Season <- lme(ln.Bac_q1 ~ Bac_Ampti + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, random = ~ 1 | Season, data = HNF_Bac_A, method = "ML")
AIC(Bacq1_BacS.0, Bacq1_BacS.Cr, Bacq1_BacS.Season)
summary(Bacq1_BacS.Cr)

### plotting
p_BacSelect_Bacq1 <- HNF_Bac_A %>% 
  select(ln.Bac_q1, ln.HNF_q1, Bac_Ampti, HNF_Ampti, Cruise) %>%
  ggplot(aes(x = Bac_Ampti, y = ln.Bac_q1)) + 
    geom_point(aes(color = Cruise), size = 3) + 
    geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "solid") + 
    #geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "dotted") + 
    scale_colour_viridis(alpha = 0.7, discrete=TRUE) + 
    labs(x = expression("Deterministic assembly processes (\U03B1MPTI) of Bacteria community "),
         y = expression(atop("Log[ Bacterial Shannon diversity ]"))) + 
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.text.x = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    )

p_HNFq1_BacSelect <- HNF_Bac_A %>% 
  select(ln.Bac_q1, ln.HNF_q1, Bac_Ampti, HNF_Ampti, Cruise) %>%
  ggplot(aes(x = ln.HNF_q1, y = Bac_Ampti)) + 
    geom_point(aes(color = Cruise), size = 3) + 
    geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "solid") + 
    #geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "solid") + 
    scale_colour_viridis(alpha = 0.7, discrete=TRUE) + 
    labs(x = expression("Log[ HNF Shannon diversity ]"),
         y = expression(atop("Deterministic assembly processes (\U03B1MPTI)", "of Bacteria community"))) + 
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.text.x = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    )

legend <- get_legend(
  p_HNFq1_BacSelect + theme(legend.box.margin = margin(0, 0, 0, 12))
)
p_HNFq1_BacSelect_Bacq1 <- plot_grid(
  plot_grid(p_BacSelect_Bacq1 + theme(legend.position = "none"),
            p_HNFq1_BacSelect + theme(legend.position = "none"),
            ncol = 1, labels = "AUTO", hjust = 0),
  legend, rel_widths = c(3, .4))
p_HNFq1_BacSelect_Bacq1
ggsave(p_HNFq1_BacSelect_Bacq1, file = "D:/Manuscripts/PdPy_Div_MS/ms_Figs/Fig4_HNFq1_BacAmpd_Bacq1.png",
       dpi = 600, width = 36, height = 42, units = "cm")
##### HNFq1 -> Bac selection -> Bacq1 ##########

##### Bacq1 -> HNF selection -> HNFq1 ##########
### linear model testing
HNFS_Bacq1.0 <- lm(HNF_Ampti ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, data = HNF_Bac_A)
HNFS_Bacq1.Cr <- lme(HNF_Ampti ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, random = ~ 1 | Cruise, data = HNF_Bac_A, method = "ML")
HNFS_Bacq1.Season <- lme(HNF_Ampti ~ ln.Bac_q1 + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, random = ~ 1 | Season, data = HNF_Bac_A, method = "ML")
AIC(HNFS_Bacq1.0, HNFS_Bacq1.Cr, HNFS_Bacq1.Season)
summary(HNFS_Bacq1.Cr)

HNFq1_HNFS.0 <- lm(ln.HNF_q1 ~ HNF_Ampti + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, data = HNF_Bac_A)
HNFq1_HNFS.Cr <- lme(ln.HNF_q1 ~ HNF_Ampti + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, random = ~ 1 | Cruise, data = HNF_Bac_A, method = "ML")
HNFq1_HNFS.Season <- lme(ln.HNF_q1 ~ HNF_Ampti + ln.Temp + ln.Sal + ln.PAR + ln.DIN + ln.PO3, random = ~ 1 | Season, data = HNF_Bac_A, method = "ML")
AIC(HNFq1_HNFS.0, HNFq1_HNFS.Cr, HNFq1_HNFS.Season)
summary(HNFq1_HNFS.Cr)

### plotting
p_HNFSelect_HNFq1 <- HNF_Bac_A %>% 
  select(ln.Bac_q1, ln.HNF_q1, Bac_Ampti, HNF_Ampti, Cruise) %>%
  ggplot(aes(x = HNF_Ampti, y = ln.HNF_q1)) + 
    geom_point(aes(color = Cruise), size = 3) + 
    geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "dashed") + 
    #geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "dotted") + 
    scale_colour_viridis(alpha = 0.7, discrete=TRUE) + 
    labs(x = expression("Deterministic assembly processes (\U03B1MPTI) of HNF community "),
         y = expression(atop("Log[ HNF Shannon diversity ]"))) + 
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.text.x = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    )

p_Bacq1_HNFSelect <- HNF_Bac_A %>% 
  select(ln.Bac_q1, ln.HNF_q1, Bac_Ampti, HNF_Ampti, Cruise) %>%
  ggplot(aes(x = ln.Bac_q1, y = HNF_Ampti)) + 
  geom_point(aes(color = Cruise), size = 3) + 
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "dashed") + 
  #geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = TRUE, color = "red", linetype = "dotted") + 
  scale_colour_viridis(alpha = 0.7, discrete=TRUE) + 
  labs(x = expression("Log[ Bacterial Shannon diversity ]"),
       y = expression(atop("Deterministic assembly processes (\U03B1MPTI)", "of HNF community "))) + 
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.text.x = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )
legend <- get_legend(
  p_Bacq1_HNFSelect + theme(legend.box.margin = margin(0, 0, 0, 12))
)
p_Bacq1_HNFSelect_HNFq1 <- plot_grid(
  plot_grid(p_HNFSelect_HNFq1 + theme(legend.position = "none"),
            p_Bacq1_HNFSelect + theme(legend.position = "none"),
            ncol = 1, labels = "AUTO", hjust = 0),
  legend, rel_widths = c(3, .4))
p_Bacq1_HNFSelect_HNFq1
ggsave(p_Bacq1_HNFSelect_HNFq1, file = "D:/Manuscripts/PdPy_Div_MS/ms_Figs/Fig5_Bacq1_HNFAmpd_HNFq1.png",
       dpi = 600, width = 36, height = 42, units = "cm")
##### Bacq1 -> HNF selection -> HNFq1 ##########

###############################################################################################
##### Testing hypothesis: HNFq0 -> Bac selection -> Bacq0 & Bacq0 -> Bac selection ->  HNFq0 ##
###############################################################################################

###############################################################################################
##### Ampti-> HNFq1_Bacq1 association  ########################################################
###############################################################################################
### statistical test
HNF_Bac_A <- HNF_Bac_A %>%
  group_by(Cruise) %>%
  mutate(
    R2 = summary(lm(ln.HNF_q1 ~ ln.Bac_q1))$r.squared,
    coef = summary(lm(ln.HNF_q1 ~ ln.Bac_q1))$coefficients[2, 1]
  ) 
ADiv_Cr <- HNF_Bac_A %>%
  summarize(
    meanBac_Ampti = mean(Bac_Ampti),
    meanHNF_Ampti = mean(HNF_Ampti),
    R2 = mean(R2),
    coef = mean(coef)
  )

# Mod_R2_BacAmpti <- glm(R2 ~ Bac_Ampti, data = HNF_Bac_A)
# Mod_R2_HNFAmpti <- glm(R2 ~ HNF_Ampti, data = HNF_Bac_A)
# summary(Mod_R2_BacAmpti)
# summary(Mod_R2_HNFAmpti)

Mod_R2_meanBacAmpti <- glm(R2 ~ meanBac_Ampti, data = ADiv_Cr)
Mod_R2_meanHNFAmpti <- glm(R2 ~ meanHNF_Ampti, data = ADiv_Cr)
summary(Mod_R2_meanBacAmpti)
summary(Mod_R2_meanHNFAmpti)

# Mod_coef_BacAmpti <- glm(coef ~ Bac_Ampti, data = HNF_Bac_A[which(HNF_Bac_A$coef > -1), ])
# Mod_coef_HNFAmpti <- glm(coef ~ HNF_Ampti, data = HNF_Bac_A[which(HNF_Bac_A$coef > -1), ])
# summary(Mod_coef_BacAmpti)
# summary(Mod_coef_HNFAmpti)

Mod_coef_meanBacAmpti <- glm(coef ~ meanBac_Ampti, data = ADiv_Cr[which(ADiv_Cr$coef > -1), ])
Mod_coef_meanHNFAmpti <- glm(coef ~ meanHNF_Ampti, data = ADiv_Cr[which(ADiv_Cr$coef > -1), ])
summary(Mod_coef_meanBacAmpti)
summary(Mod_coef_meanHNFAmpti)

### visualization
p_ADivR2_BacAmpd <- ggplot(data = ADiv_Cr, aes(x = meanBac_Ampti, y = R2)) + 
    geom_point(size = 4) +    
    geom_point(data = HNF_Bac_A, aes(x = Bac_Ampti, y = R2), size = 2, color = "grey") +
    geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "dashed") + 
    labs(x = expression(atop("Deterministic assembly processes (\U03B1MPTI)", "of Bacteria community")),
         y = bquote("R"^2~ "of diversity association")) + 
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.text.x = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    )
p_ADivR2_HNFAmpd <- ggplot(data = ADiv_Cr, aes(x = meanHNF_Ampti, y = R2)) + 
    geom_point(size = 4) +  
    geom_point(data = HNF_Bac_A, aes(x = HNF_Ampti, y = R2), size = 2, color = "grey") + 
    geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "dashed") + 
    labs(x = expression(atop("Deterministic assembly processes (\U03B1MPTI)", "of HNF community")),
         y = bquote("R"^2~ "of diversity association")) + 
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.text.x = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    )
p_ADivR2_Ampd <- plot_grid(p_ADivR2_BacAmpd, p_ADivR2_HNFAmpd, ncol = 1, labels = "AUTO", hjust = 0)
p_ADivR2_Ampd
ggsave(p_ADivR2_Ampd, file = "D:/Research/PdPy_Div_Results/Figs/Lab422 meeting_20200602/p_ADivR2_Ampd.png",
       dpi = 600, width = 36, height = 42, units = "cm")

p_ADivcoef_BacAmpd <- ggplot(data = ADiv_Cr, aes(x = meanBac_Ampti, y = coef)) + #[which(Adiv_cr$coef > -1), ]
    geom_point(size = 4) + 
    geom_point(data = HNF_Bac_A, aes(x = Bac_Ampti, y = coef), size = 2, color = "grey") + 
    geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "dashed") + 
    labs(x = expression(atop("Deterministic assembly processes (\U03B1MPTI)", "of Bacteria community")),
         y = bquote(atop("Regression coefficient", "of diversity association"))) + 
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.text.x = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    )
p_ADivcoef_HNFAmpd <- ggplot(data = ADiv_Cr, aes(x = meanHNF_Ampti, y = coef)) + #[which(ADiv_Cr$coef > -1), ]
    geom_point(size = 4) +
    geom_point(data = HNF_Bac_A, aes(x = HNF_Ampti, y = coef), size = 2, color = "grey") + 
    geom_smooth(formula = y ~ x, method = "lm", se = TRUE, linetype = "dashed") + 
    labs(x = expression(atop("Deterministic assembly processes (\U03B1MPTI)", "of HNF community")),
         y = bquote(atop("Regression coefficient", "of diversity association"))) + 
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.text.x = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    )
p_ADivcoef_Ampd <- plot_grid(p_ADivcoef_BacAmpd, p_ADivcoef_HNFAmpd, ncol = 1, labels = "AUTO", hjust = 0)
p_ADivcoef_Ampd
ggsave(p_ADivcoef_Ampd, file = "D:/Research/PdPy_Div_Results/Figs/Lab422 meeting_20200602/p_ADivcoef_Ampd.png",
       dpi = 600, width = 36, height = 42, units = "cm")

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
p_Adiv_AAssemb_pairs <- HNF_Bac_A %>%
  ggpairs(columns = c("ln.Bac_q0", "ln.HNF_q0", "ln.Bac_q1", "ln.HNF_q1", "ln.Bac_q2", "ln.HNF_q2", 
                      "Bac_Ampd_select", "HNF_Ampd_select"),
          columnLabels = c("Bacteria\nspecies\nrichness", "HNF\nspecies\nrichness",
                           "Bacteria\nShannon\ndiversity", "HNF\nShannon\ndiversity", 
                           "Bacteria\nSimpson\ndiversity", "HNF\nSimpson\ndiversity",
                           "Bacteria\n\U03B1MPTI", "HNF\n\U03B1MPTI"),
          #mapping = ggplot2::aes(colour = Cruise),
          upper = list(continuous = cor_fun),
          lower = list(continuous = fit_fun)) +
  theme(strip.text.x = element_text(color = "black", size = 14),
        strip.text.y = element_text(angle = 45, color = "black", size = 14))
p_Adiv_AAssemb_pairs
ggsave(p_Adiv_AAssemb_pairs, file = "D:/Research/PdPy_Div_Results/Figs/PR2_3/FigS1_Adiv_AMPTI_pairs.png",
       dpi = 600, width = 34, height = 28, units = "cm")

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

##### Visualization ##########
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
