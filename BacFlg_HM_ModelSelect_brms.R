######## Loading necessary libraries ########################################################################################
if (!require(tidyverse)) {
  install.packages("tidyverse", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(tidyverse)
}else{library(tidyverse)}

if (!require(vegan)) {
  install.packages("vegan", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(vegan)
}else{library(vegan)}

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

if (!require(shinystan)) {
  install.packages("shinystan", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(shinystan)
}else{library(shinystan)}

# if (!require(rstan)) {
#   install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
#   library(rstan)
# }else{library(rstan)
#   rstan_options(auto_write = TRUE)
#   options(mc.cores = parallel::detectCores())}

library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

if (!require(loo)) {
  install.packages("loo", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(loo)
}else{library(loo)}

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

if (!require(ggpmisc)) {
  install.packages("ggpmisc", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(ggpmisc)
}else{library(ggpmisc)}
######## Loading necessary libraries ########################################################################################

##### Reading data from Github ##############################################################################################
Dat.raw <- read.table(file = "D:/Research/PdPy_Div/BacFlg.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
##### Reading data from Github ##############################################################################################
Dat <- Dat.raw %>%
  drop_na() %>%
  select("SampleID.Bac.", "Cruise", "Station", "Year", "Month", "Season", 
         "q0_obs_HNF_inext", "q1_obs_HNF_inext", "q2_obs_HNF_inext", 
         "q0_obs_Bac_inext", "q1_obs_Bac_inext", "q2_obs_Bac_inext", 
         "prey_biomass", "predator_biomass", 
         "ABB", "HBB", "PNFB", "HNFB", "APPBR", "HPPBR", "TPPBR",
         "temperature", "salinity", "PAR", "nitrite", "nitrate", "phosphate", "Chla", "Q4temp") %>%
  mutate(Season = factor(Season, levels = c("Spring", "Summer", "Autumn")),
         q0_obs_HNF_inext = as.numeric(q0_obs_HNF_inext),
         q0_obs_Bac_inext = as.numeric(q0_obs_Bac_inext),
         prey_biom = prey_biomass, 
         predator_biom = predator_biomass,
         ABac_biom = ABB,
         Hbac_biom = HBB, 
         PNfl_biom = PNFB,
         HNfl_biom = HNFB,
         A_PPBR = APPBR, 
         H_PPBR = HPPBR,
         T_PPBR = TPPBR,
         TN = nitrite + nitrate,
         TP = phosphate,
         ln_q0_HNF = log(q0_obs_HNF_inext),
         ln_q1_HNF = log(q1_obs_HNF_inext),
         ln_q2_HNF = log(q2_obs_HNF_inext),
         ln_q0_Bac = log(q0_obs_Bac_inext),
         ln_q1_Bac = log(q1_obs_Bac_inext),
         ln_q2_Bac = log(q2_obs_Bac_inext)) %>%
  select("SampleID.Bac.", "Cruise", "Station", "Year", "Month", "Season", 
         "q0_obs_HNF_inext", "q1_obs_HNF_inext", "q2_obs_HNF_inext", 
         "q0_obs_Bac_inext", "q1_obs_Bac_inext", "q2_obs_Bac_inext", 
         "ln_q0_HNF", "ln_q1_HNF", "ln_q2_HNF", "ln_q0_Bac", "ln_q1_Bac", "ln_q2_Bac",
         "prey_biom", "predator_biomass",
         "ABac_biom", "Hbac_biom", "PNfl_biom", "HNfl_biom",
         "A_PPBR", "H_PPBR", "T_PPBR", 
         "temperature", "salinity", "PAR", "TN", "TP", "Chla", "Q4temp") 
str(Dat)

ggpairs(Dat, columns = c("q0_obs_HNF_inext", "q1_obs_HNF_inext", "q2_obs_HNF_inext", 
                         "q0_obs_Bac_inext", "q1_obs_Bac_inext", "q2_obs_Bac_inext"),
        ggplot2::aes(colour = Season, alpha = 0.4))

Bac_labs <- c("q0 (richness)", "q1 (Shannon)", "q2 (inverse Simpson)")
names(Bac_labs) <- c("q0_obs_Bac_inext", "q1_obs_Bac_inext", "q2_obs_Bac_inext")

Dat %>%
  gather(key = "Bac_Div_q", value = "Bac_Div", 
         q0_obs_Bac_inext, q1_obs_Bac_inext, q2_obs_Bac_inext) %>%
  ggplot() +
    geom_point(aes(y = q0_obs_HNF_inext, x = Bac_Div)) + 
    geom_smooth(aes(y = q0_obs_HNF_inext, x = Bac_Div), method = "lm") +
    facet_grid(cols = vars(Bac_Div_q), rows = vars(Season), 
               scales = "free_x", labeller = labeller(Bac_Div_q = Bac_labs)) + 
    labs(x = "Bacteria diversity index", y = "")

q0q0 <- Dat %>%
  group_by(Season) %>%
  mutate(lm_pval = ifelse(summary(lm(q0_obs_HNF_inext ~ q0_obs_Bac_inext))$coefficients[2,4] < 0.05, 0, 1)) %>% 
  ggplot(aes(x = q0_obs_Bac_inext, y = q0_obs_HNF_inext)) +
    geom_point() + 
    facet_grid(rows = vars(Season)) + 
    geom_smooth(method = "lm", formula = y ~ x, aes(linetype = factor(lm_pval, levels = c(0, 1))), show.legend = FALSE) +
    scale_linetype_manual(values = c("solid", "dotted")) + 
    labs(x = "Bacteria q0 diversity (richness)", y = "Heterotrophic nanoflagellate q0 diversity (richness)") + 
    stat_fit_glance(aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")), 
                    method = 'lm', method.args = list(formula = y ~ x),
                    label.x= "right", label.y = 0.85, size = 4) + 
    stat_poly_eq(formula = y ~ x,#, label = paste(..rr.label..)
                 label.x.npc = "right", label.y.npc = 0.75, parse = TRUE, size = 4)

q1q0 <- Dat %>%
  group_by(Season) %>%
  mutate(lm_pval = ifelse(summary(lm(q0_obs_HNF_inext ~ q1_obs_Bac_inext))$coefficients[2,4] < 0.05, 0, 1)) %>% 
  ggplot(aes(x = q1_obs_Bac_inext, y = q0_obs_HNF_inext)) +
    geom_point() + 
    geom_smooth(method = "lm", formula = y ~ x, aes(linetype = factor(lm_pval, levels = c(0, 1))), show.legend = FALSE) +
    scale_linetype_manual(values = c("dotted")) + 
    facet_grid(rows = vars(Season)) + 
    labs(x = "Bacteria q1 diversity (Shannon)", y = "Heterotrophic nanoflagellate q0 diversity (richness)") + 
    stat_fit_glance(aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")), 
                    method = 'lm', method.args = list(formula = y ~ x),
                    label.x= "right", label.y = 0.85, size = 4) + 
    stat_poly_eq(formula = y ~ x,#, label = paste(..rr.label..)
                 label.x.npc = "right", label.y.npc = 0.75, parse = TRUE, size = 4)
q2q0 <- Dat %>%
  group_by(Season) %>%
  mutate(lm_pval = ifelse(summary(lm(q0_obs_HNF_inext ~ q2_obs_Bac_inext))$coefficients[2,4] < 0.05, 0, 1)) %>% 
  ggplot(aes(x = q2_obs_Bac_inext, y = q0_obs_HNF_inext)) +
    geom_point() + 
    geom_smooth(method = "lm", formula = y ~ x, aes(linetype = factor(lm_pval, levels = c(0, 1))), show.legend = FALSE) +
    scale_linetype_manual(values = c("dotted")) + 
    facet_grid(rows = vars(Season)) + 
    labs(x = "Bacteria q2 diversity (inverse Simpson)", y = "Heterotrophic nanoflagellate q0 diversity (richness)") + 
    stat_fit_glance(aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")), 
                    method = 'lm', method.args = list(formula = y ~ x),
                    label.x= "right", label.y = 0.85, size = 4) + 
    stat_poly_eq(formula = y ~ x,#, label = paste(..rr.label..)
                 label.x.npc = "right", label.y.npc = 0.75, parse = TRUE, size = 4)
plot_grid(q0q0, q1q0, q2q0, nrow = 1)

ggpairs(Dat, columns = c("q0_obs_HNF_inext", "q1_obs_HNF_inext", "q2_obs_HNF_inext", 
                         "q0_obs_Bac_inext", "q1_obs_Bac_inext", "q2_obs_Bac_inext"),
        ggplot2::aes(colour = Season, alpha = 0.4))


#############################################################################################################################
##### Finding the best model  ###############################################################################################
#############################################################################################################################
##### Bac ~ Flg #########################################################################################################
### Step 1
fn0.1 <- bf(log(q0_obs_HNF_inext) ~ log(q0_obs_HNF_inext) + (1 + q0_obs_HNF_inext | Cruise))
fn0.2 <- bf(q0_obs_HNF_inext ~ q0_obs_HNF_inext + (1 + q0_obs_HNF_inext | Station))
fn0.3 <- bf(q0_obs_HNF_inext ~ q0_obs_HNF_inext + (1 + q0_obs_HNF_inext | Month))
fn0.4 <- bf(q0_obs_HNF_inext ~ q0_obs_HNF_inext + (1 + q0_obs_HNF_inext | Season))
Mod0.1 <- brm(formula = fn0.1, data = Dat, family = gaussian(), control = list(adapt_delta = 0.9))










##### Reading data from Github ##############################################################################################
bio07_scale = read.csv(file = "https://raw.githubusercontent.com/OscarFHC/NLA_Data/master/Final_cleaned/bio07_scaled.csv", 
                 header = TRUE, stringsAsFactors = FALSE, fill = TRUE) 

bio12_scale = read.csv(file = "https://raw.githubusercontent.com/OscarFHC/NLA_Data/master/Final_cleaned/bio12_scaled.csv", 
                 header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
##### Reading data from Github ##############################################################################################

#############################################################################################################################
##### Finding the best model  ###############################################################################################
#############################################################################################################################
##### 2012 #########################################################################################################
### Step 1
fn0.1_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + 
                 (1 + ln_zpSRr*WOmni | ECO3))
fn0.2_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_phySRr + 
                 (1 + ln_zpSRr*WOmni + ln_phySRr | ECO3))
fn0.3_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_zpDen + 
                 (1 + ln_zpSRr*WOmni + ln_zpDen | ECO3))
fn0.4_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_NTL + 
                 (1 + ln_zpSRr*WOmni + ln_NTL | ECO3))
fn0.5_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + 
                 (1 + ln_zpSRr*WOmni + ln_PTL | ECO3))
fn0.6_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + 
                 (1 + ln_zpSRr*WOmni + Lat | ECO3))
fn0.7_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_DOC + 
                 (1 + ln_zpSRr*WOmni + ln_DOC | ECO3))
fn0.8_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_Silica + 
                 (1 + ln_zpSRr*WOmni + ln_Silica | ECO3))
fn0.9_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + pH + 
                 (1 + ln_zpSRr*WOmni + pH | ECO3))
fn0.10_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_Temp + 
                  (1 + ln_zpSRr*WOmni + ln_Temp | ECO3))
fn0.11_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_Sol + 
                  (1 + ln_zpSRr*WOmni + ln_Sol | ECO3))
fn0.12_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_Cond + 
                  (1 + ln_zpSRr*WOmni + ln_Cond | ECO3))
Mod0.1_12 <- brm(formula = fn0.1_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.2_12 <- brm(formula = fn0.2_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.3_12 <- brm(formula = fn0.3_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.4_12 <- brm(formula = fn0.4_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.5_12 <- brm(formula = fn0.5_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.6_12 <- brm(formula = fn0.6_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.7_12 <- brm(formula = fn0.7_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.8_12 <- brm(formula = fn0.8_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.9_12 <- brm(formula = fn0.9_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.10_12 <- brm(formula = fn0.10_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.11_12 <- brm(formula = fn0.11_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.12_12 <- brm(formula = fn0.12_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
loo0.1_12 <- loo(Mod0.1_12, reloo = TRUE)
loo0.2_12 <- loo(Mod0.2_12, reloo = TRUE)
loo0.3_12 <- loo(Mod0.3_12, reloo = TRUE)
loo0.4_12 <- loo(Mod0.4_12, reloo = TRUE)
loo0.5_12 <- loo(Mod0.5_12, reloo = TRUE)
loo0.6_12 <- loo(Mod0.6_12, reloo = TRUE)
loo0.7_12 <- loo(Mod0.7_12, reloo = TRUE)
loo0.8_12 <- loo(Mod0.8_12, reloo = TRUE)
loo0.9_12 <- loo(Mod0.9_12, reloo = TRUE)
loo0.10_12 <- loo(Mod0.10_12, reloo = TRUE)
loo0.11_12 <- loo(Mod0.11_12, reloo = TRUE)
loo0.12_12 <- loo(Mod0.12_12, reloo = TRUE)

looic <- rbind(loo0.1_12$estimates["elpd_loo","Estimate"], loo0.2_12$estimates["elpd_loo","Estimate"],
               loo0.3_12$estimates["elpd_loo","Estimate"], loo0.4_12$estimates["elpd_loo","Estimate"],
               loo0.5_12$estimates["elpd_loo","Estimate"], loo0.6_12$estimates["elpd_loo","Estimate"],
               loo0.7_12$estimates["elpd_loo","Estimate"], loo0.8_12$estimates["elpd_loo","Estimate"],
               loo0.9_12$estimates["elpd_loo","Estimate"], loo0.10_12$estimates["elpd_loo","Estimate"],
               loo0.11_12$estimates["elpd_loo","Estimate"], loo0.12_12$estimates["elpd_loo","Estimate"])
looic_SE <- rbind(loo0.1_12$estimates["elpd_loo","SE"], loo0.2_12$estimates["elpd_loo","SE"],
                  loo0.3_12$estimates["elpd_loo","SE"], loo0.4_12$estimates["elpd_loo","SE"],
                  loo0.5_12$estimates["elpd_loo","SE"], loo0.6_12$estimates["elpd_loo","SE"],
                  loo0.7_12$estimates["elpd_loo","SE"], loo0.8_12$estimates["elpd_loo","SE"],
                  loo0.9_12$estimates["elpd_loo","SE"], loo0.10_12$estimates["elpd_loo","SE"],
                  loo0.11_12$estimates["elpd_loo","SE"], loo0.12_12$estimates["elpd_loo","SE"])
stacking_wts <- stacking_weights(cbind(loo0.1_12$pointwise[, "elpd_loo"], loo0.2_12$pointwise[, "elpd_loo"],
                                       loo0.3_12$pointwise[, "elpd_loo"], loo0.4_12$pointwise[, "elpd_loo"],
                                       loo0.5_12$pointwise[, "elpd_loo"], loo0.6_12$pointwise[, "elpd_loo"],
                                       loo0.7_12$pointwise[, "elpd_loo"], loo0.8_12$pointwise[, "elpd_loo"],
                                       loo0.9_12$pointwise[, "elpd_loo"], loo0.10_12$pointwise[, "elpd_loo"],
                                       loo0.11_12$pointwise[, "elpd_loo"], loo0.12_12$pointwise[, "elpd_loo"]))
Step1_12 <- round(cbind(looic, looic_SE, stacking_wts), 3)
# write.table(Step1_12, file = "D:/Research/OmniBEF_NLA/NLA_HM/ModResults/brms_Omni/ModSelec_12_step1.csv", sep = ",",
#             row.names = TRUE, col.names = TRUE)


### Step 2
fn0.5_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + 
                 (1 + ln_zpSRr*WOmni + ln_PTL | ECO3))
fn1.1_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_phySRr +
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_phySRr| ECO3))
fn1.2_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_zpDen + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_zpDen | ECO3))
fn1.3_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL | ECO3))
fn1.4_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + Lat + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + Lat | ECO3))
fn1.5_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_DOC + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_DOC | ECO3))
fn1.6_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_Silica + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_Silica | ECO3))
fn1.7_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + pH + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + pH | ECO3))
fn1.8_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_Temp + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_Temp | ECO3))
fn1.9_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_Sol +
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_Sol | ECO3))
fn1.10_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_Cond + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_Cond | ECO3))
Mod0.5_12 <- brm(formula = fn0.5_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.1_12 <- brm(formula = fn1.1_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.2_12 <- brm(formula = fn1.2_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.3_12 <- brm(formula = fn1.3_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.4_12 <- brm(formula = fn1.4_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.5_12 <- brm(formula = fn1.5_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.6_12 <- brm(formula = fn1.6_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.7_12 <- brm(formula = fn1.7_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.8_12 <- brm(formula = fn1.8_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.9_12 <- brm(formula = fn1.9_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.10_12 <- brm(formula = fn1.10_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
loo0.5_12 <- loo(Mod0.5_12, reloo = TRUE)
loo1.1_12 <- loo(Mod1.1_12, reloo = TRUE)
loo1.2_12 <- loo(Mod1.2_12, reloo = TRUE)
loo1.3_12 <- loo(Mod1.3_12, reloo = TRUE)
loo1.4_12 <- loo(Mod1.4_12, reloo = TRUE)
loo1.5_12 <- loo(Mod1.5_12, reloo = TRUE)
loo1.6_12 <- loo(Mod1.6_12, reloo = TRUE)
loo1.7_12 <- loo(Mod1.7_12, reloo = TRUE)
loo1.8_12 <- loo(Mod1.8_12, reloo = TRUE)
loo1.9_12 <- loo(Mod1.9_12, reloo = TRUE)
loo1.10_12 <- loo(Mod1.10_12, reloo = TRUE)

looic <- rbind(loo0.5_12$estimates["elpd_loo","Estimate"],
               loo1.1_12$estimates["elpd_loo","Estimate"], loo1.2_12$estimates["elpd_loo","Estimate"],
               loo1.3_12$estimates["elpd_loo","Estimate"], loo1.4_12$estimates["elpd_loo","Estimate"],
               loo1.5_12$estimates["elpd_loo","Estimate"], loo1.6_12$estimates["elpd_loo","Estimate"],
               loo1.7_12$estimates["elpd_loo","Estimate"], loo1.8_12$estimates["elpd_loo","Estimate"],
               loo1.9_12$estimates["elpd_loo","Estimate"], loo1.10_12$estimates["elpd_loo","Estimate"])
looic_SE <- rbind(loo0.5_12$estimates["elpd_loo","SE"],
                  loo1.1_12$estimates["elpd_loo","SE"], loo1.2_12$estimates["elpd_loo","SE"],
                  loo1.3_12$estimates["elpd_loo","SE"], loo1.4_12$estimates["elpd_loo","SE"],
                  loo1.5_12$estimates["elpd_loo","SE"], loo1.6_12$estimates["elpd_loo","SE"],
                  loo1.7_12$estimates["elpd_loo","SE"], loo1.8_12$estimates["elpd_loo","SE"],
                  loo1.9_12$estimates["elpd_loo","SE"], loo1.10_12$estimates["elpd_loo","SE"])
stacking_wts <- stacking_weights(cbind(loo0.5_12$pointwise[, "elpd_loo"],
                                       loo1.1_12$pointwise[, "elpd_loo"], loo1.2_12$pointwise[, "elpd_loo"],
                                       loo1.3_12$pointwise[, "elpd_loo"], loo1.4_12$pointwise[, "elpd_loo"],
                                       loo1.5_12$pointwise[, "elpd_loo"], loo1.6_12$pointwise[, "elpd_loo"],
                                       loo1.7_12$pointwise[, "elpd_loo"], loo1.8_12$pointwise[, "elpd_loo"],
                                       loo1.9_12$pointwise[, "elpd_loo"], loo1.10_12$pointwise[, "elpd_loo"]))
Step2_12 <- round(cbind(looic, looic_SE, stacking_wts), 3)
# write.table(Step2_12, file = "D:/Research/OmniBEF_NLA/NLA_HM/ModResults/brms_Omni/ModSelec_12_step2.csv", sep = ",",
#             row.names = TRUE, col.names = TRUE)
#summary(Mod1.3_12)
### Step 3
fn1.3_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL | ECO3))
fn2.1_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr | ECO3))
fn2.2_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_zpDen + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_zpDen | ECO3))
fn2.3_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + Lat + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + Lat | ECO3))
fn2.4_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_DOC + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_DOC | ECO3))
fn2.5_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_Silica + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_Silica | ECO3))
fn2.6_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + pH + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + pH | ECO3))
fn2.7_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_Temp + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_Temp | ECO3))
fn2.8_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_Sol + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_Sol | ECO3))
fn2.9_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_Cond + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_Cond | ECO3))
Mod1.3_12 <- brm(formula = fn1.3_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.1_12 <- brm(formula = fn2.1_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.2_12 <- brm(formula = fn2.2_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.3_12 <- brm(formula = fn2.3_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.4_12 <- brm(formula = fn2.4_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.5_12 <- brm(formula = fn2.5_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.6_12 <- brm(formula = fn2.6_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.7_12 <- brm(formula = fn2.7_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.8_12 <- brm(formula = fn2.8_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.9_12 <- brm(formula = fn2.9_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
loo1.3_12 <- loo(Mod1.3_12, reloo = TRUE)
loo2.1_12 <- loo(Mod2.1_12, reloo = TRUE)
loo2.2_12 <- loo(Mod2.2_12, reloo = TRUE)
loo2.3_12 <- loo(Mod2.3_12, reloo = TRUE)
loo2.4_12 <- loo(Mod2.4_12, reloo = TRUE)
loo2.5_12 <- loo(Mod2.5_12, reloo = TRUE)
loo2.6_12 <- loo(Mod2.6_12, reloo = TRUE)
loo2.7_12 <- loo(Mod2.7_12, reloo = TRUE)
loo2.8_12 <- loo(Mod2.8_12, reloo = TRUE)
loo2.9_12 <- loo(Mod2.9_12, reloo = TRUE)

looic <- rbind(loo1.3_12$estimates["elpd_loo","Estimate"],
               loo2.1_12$estimates["elpd_loo","Estimate"], loo2.2_12$estimates["elpd_loo","Estimate"],
               loo2.3_12$estimates["elpd_loo","Estimate"], loo2.4_12$estimates["elpd_loo","Estimate"],
               loo2.5_12$estimates["elpd_loo","Estimate"], loo2.6_12$estimates["elpd_loo","Estimate"],
               loo2.7_12$estimates["elpd_loo","Estimate"], loo2.8_12$estimates["elpd_loo","Estimate"],
               loo2.9_12$estimates["elpd_loo","Estimate"])
looic_SE <- rbind(loo1.3_12$estimates["elpd_loo","SE"],
                  loo2.1_12$estimates["elpd_loo","SE"], loo2.2_12$estimates["elpd_loo","SE"],
                  loo2.3_12$estimates["elpd_loo","SE"], loo2.4_12$estimates["elpd_loo","SE"],
                  loo2.5_12$estimates["elpd_loo","SE"], loo2.6_12$estimates["elpd_loo","SE"],
                  loo2.7_12$estimates["elpd_loo","SE"], loo2.8_12$estimates["elpd_loo","SE"],
                  loo2.9_12$estimates["elpd_loo","SE"])
stacking_wts <- stacking_weights(cbind(loo1.3_12$pointwise[, "elpd_loo"],
                                       loo2.1_12$pointwise[, "elpd_loo"], loo2.2_12$pointwise[, "elpd_loo"],
                                       loo2.3_12$pointwise[, "elpd_loo"], loo2.4_12$pointwise[, "elpd_loo"],
                                       loo2.5_12$pointwise[, "elpd_loo"], loo2.6_12$pointwise[, "elpd_loo"],
                                       loo2.7_12$pointwise[, "elpd_loo"], loo2.8_12$pointwise[, "elpd_loo"],
                                       loo2.9_12$pointwise[, "elpd_loo"]))
Step3_12 <- round(cbind(looic, looic_SE, stacking_wts), 3)
# write.table(Step3_12, file = "D:/Research/OmniBEF_NLA/NLA_HM/ModResults/brms_Omni/ModSelec_12_step3.csv", sep = ",",
#             row.names = TRUE, col.names = TRUE)
#summary(Mod2.1_12)

### Step 4
fn2.1_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr | ECO3))
fn3.1_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_zpDen + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_zpDen | ECO3))
fn3.2_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + Lat + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + Lat | ECO3))
fn3.3_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC | ECO3))
fn3.4_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_Silica + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_Silica | ECO3))
fn3.5_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + pH + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + pH | ECO3))
fn3.6_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_Temp + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_Temp | ECO3))
fn3.7_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_Sol + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_Sol | ECO3))
fn3.8_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_Cond + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_Cond | ECO3))
Mod2.1_12 <- brm(formula = fn2.1_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.1_12 <- brm(formula = fn3.1_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.2_12 <- brm(formula = fn3.2_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.3_12 <- brm(formula = fn3.3_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.4_12 <- brm(formula = fn3.4_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.5_12 <- brm(formula = fn3.5_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.6_12 <- brm(formula = fn3.6_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.7_12 <- brm(formula = fn3.7_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.8_12 <- brm(formula = fn3.8_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
loo2.1_12 <- loo(Mod2.1_12, reloo = TRUE)
loo3.1_12 <- loo(Mod3.1_12, reloo = TRUE)
loo3.2_12 <- loo(Mod3.2_12, reloo = TRUE)
loo3.3_12 <- loo(Mod3.3_12, reloo = TRUE)
loo3.4_12 <- loo(Mod3.4_12, reloo = TRUE)
loo3.5_12 <- loo(Mod3.5_12, reloo = TRUE)
loo3.6_12 <- loo(Mod3.6_12, reloo = TRUE)
loo3.7_12 <- loo(Mod3.7_12, reloo = TRUE)
loo3.8_12 <- loo(Mod3.8_12, reloo = TRUE)

looic <- rbind(loo2.1_12$estimates["elpd_loo","Estimate"],
               loo3.1_12$estimates["elpd_loo","Estimate"], loo3.2_12$estimates["elpd_loo","Estimate"],
               loo3.3_12$estimates["elpd_loo","Estimate"], loo3.4_12$estimates["elpd_loo","Estimate"],
               loo3.5_12$estimates["elpd_loo","Estimate"], loo3.6_12$estimates["elpd_loo","Estimate"],
               loo3.7_12$estimates["elpd_loo","Estimate"], loo3.8_12$estimates["elpd_loo","Estimate"])
looic_SE <- rbind(loo2.1_12$estimates["elpd_loo","SE"],
                  loo3.1_12$estimates["elpd_loo","SE"], loo3.2_12$estimates["elpd_loo","SE"],
                  loo3.3_12$estimates["elpd_loo","SE"], loo3.4_12$estimates["elpd_loo","SE"],
                  loo3.5_12$estimates["elpd_loo","SE"], loo3.6_12$estimates["elpd_loo","SE"],
                  loo3.7_12$estimates["elpd_loo","SE"], loo3.8_12$estimates["elpd_loo","SE"])
stacking_wts <- stacking_weights(cbind(loo2.1_12$pointwise[, "elpd_loo"],
                                       loo3.1_12$pointwise[, "elpd_loo"], loo3.2_12$pointwise[, "elpd_loo"],
                                       loo3.3_12$pointwise[, "elpd_loo"], loo3.4_12$pointwise[, "elpd_loo"],
                                       loo3.5_12$pointwise[, "elpd_loo"], loo3.6_12$pointwise[, "elpd_loo"],
                                       loo3.7_12$pointwise[, "elpd_loo"], loo3.8_12$pointwise[, "elpd_loo"]))
Step4_12 <- round(cbind(looic, looic_SE, stacking_wts), 3)
# write.table(Step4_12, file = "D:/Research/OmniBEF_NLA/NLA_HM/ModResults/brms_Omni/ModSelec_12_step4.csv", sep = ",",
#             row.names = TRUE, col.names = TRUE)
summary(Mod3.3_12)
### Step 5
fn3.3_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC | ECO3))
fn4.1_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + ln_zpDen + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + ln_zpDen | ECO3))
fn4.2_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + Lat + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + Lat | ECO3))
fn4.3_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + ln_Silica + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + ln_Silica | ECO3))
fn4.4_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + pH + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + pH | ECO3))
fn4.5_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + ln_Temp + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + ln_Temp | ECO3))
fn4.6_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + ln_Sol + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + ln_Sol | ECO3))
fn4.7_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + ln_Cond + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + ln_Cond | ECO3))
Mod3.3_12 <- brm(formula = fn3.3_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod4.1_12 <- brm(formula = fn4.1_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod4.2_12 <- brm(formula = fn4.2_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod4.3_12 <- brm(formula = fn4.3_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod4.4_12 <- brm(formula = fn4.4_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod4.5_12 <- brm(formula = fn4.5_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod4.6_12 <- brm(formula = fn4.6_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod4.7_12 <- brm(formula = fn4.7_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
loo3.3_12 <- loo(Mod3.3_12, reloo = TRUE)
loo4.1_12 <- loo(Mod4.1_12, reloo = TRUE)
loo4.2_12 <- loo(Mod4.2_12, reloo = TRUE)
loo4.3_12 <- loo(Mod4.3_12, reloo = TRUE)
loo4.4_12 <- loo(Mod4.4_12, reloo = TRUE)
loo4.5_12 <- loo(Mod4.5_12, reloo = TRUE)
loo4.6_12 <- loo(Mod4.6_12, reloo = TRUE)
loo4.7_12 <- loo(Mod4.7_12, reloo = TRUE)

looic <- rbind(loo3.3_12$estimates["elpd_loo","Estimate"],
               loo4.1_12$estimates["elpd_loo","Estimate"], loo4.2_12$estimates["elpd_loo","Estimate"],
               loo4.3_12$estimates["elpd_loo","Estimate"], loo4.4_12$estimates["elpd_loo","Estimate"],
               loo4.5_12$estimates["elpd_loo","Estimate"], loo4.6_12$estimates["elpd_loo","Estimate"],
               loo4.7_12$estimates["elpd_loo","Estimate"])
looic_SE <- rbind(loo3.3_12$estimates["elpd_loo","SE"],
                  loo4.1_12$estimates["elpd_loo","SE"], loo4.2_12$estimates["elpd_loo","SE"],
                  loo4.3_12$estimates["elpd_loo","SE"], loo4.4_12$estimates["elpd_loo","SE"],
                  loo4.5_12$estimates["elpd_loo","SE"], loo4.6_12$estimates["elpd_loo","SE"],
                  loo4.7_12$estimates["elpd_loo","SE"])
stacking_wts <- stacking_weights(cbind(loo3.3_12$pointwise[, "elpd_loo"],
                                       loo4.1_12$pointwise[, "elpd_loo"], loo4.2_12$pointwise[, "elpd_loo"],
                                       loo4.3_12$pointwise[, "elpd_loo"], loo4.4_12$pointwise[, "elpd_loo"],
                                       loo4.5_12$pointwise[, "elpd_loo"], loo4.6_12$pointwise[, "elpd_loo"],
                                       loo4.7_12$pointwise[, "elpd_loo"]))
Step5_12 <- round(cbind(looic, looic_SE, stacking_wts), 3)
# write.table(Step5_12, file = "D:/Research/OmniBEF_NLA/NLA_HM/ModResults/brms_Omni/ModSelec_12_step5.csv", sep = ",",
#             row.names = TRUE, col.names = TRUE)
#summary(Mod4.4_12)

### Step 6
fn4.2_12 <- bf(ln_phyDen ~ ln_phySRr + ln_NTL + ln_PTL + ln_zpSRr + Lat)
fn5.1_12 <- bf(ln_phyDen ~ ln_phySRr + ln_NTL + ln_PTL + ln_zpSRr + Lat + ln_zpDen)
fn5.2_12 <- bf(ln_phyDen ~ ln_phySRr + ln_NTL + ln_PTL + ln_zpSRr + Lat + ln_DOC)
fn5.3_12 <- bf(ln_phyDen ~ ln_phySRr + ln_NTL + ln_PTL + ln_zpSRr + Lat + ln_Silica)
fn5.4_12 <- bf(ln_phyDen ~ ln_phySRr + ln_NTL + ln_PTL + ln_zpSRr + Lat + pH)
fn5.5_12 <- bf(ln_phyDen ~ ln_phySRr + ln_NTL + ln_PTL + ln_zpSRr + Lat + ln_Temp)
fn5.6_12 <- bf(ln_phyDen ~ ln_phySRr + ln_NTL + ln_PTL + ln_zpSRr + Lat + ln_Sol)
fn5.7_12 <- bf(ln_phyDen ~ ln_phySRr + ln_NTL + ln_PTL + ln_zpSRr + Lat + ln_Cond)
Mod4.2_12 <- brm(formula = fn4.2_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod5.1_12 <- brm(formula = fn5.1_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod5.2_12 <- brm(formula = fn5.2_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod5.3_12 <- brm(formula = fn5.3_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod5.4_12 <- brm(formula = fn5.4_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod5.5_12 <- brm(formula = fn5.5_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod5.6_12 <- brm(formula = fn5.6_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod5.7_12 <- brm(formula = fn5.7_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
loo4.2_12 <- loo(Mod4.2_12, reloo = TRUE)
loo5.1_12 <- loo(Mod5.1_12, reloo = TRUE)
loo5.2_12 <- loo(Mod5.2_12, reloo = TRUE)
loo5.3_12 <- loo(Mod5.3_12, reloo = TRUE)
loo5.4_12 <- loo(Mod5.4_12, reloo = TRUE)
loo5.5_12 <- loo(Mod5.5_12, reloo = TRUE)
loo5.6_12 <- loo(Mod5.6_12, reloo = TRUE)
loo5.7_12 <- loo(Mod5.7_12, reloo = TRUE)

looic <- rbind(loo4.2_12$estimates["elpd_loo","Estimate"],
               loo5.1_12$estimates["elpd_loo","Estimate"], loo5.2_12$estimates["elpd_loo","Estimate"],
               loo5.3_12$estimates["elpd_loo","Estimate"], loo5.4_12$estimates["elpd_loo","Estimate"],
               loo5.5_12$estimates["elpd_loo","Estimate"], loo5.6_12$estimates["elpd_loo","Estimate"],
               loo5.7_12$estimates["elpd_loo","Estimate"])
looic_SE <- rbind(loo4.2_12$estimates["elpd_loo","SE"],
                  loo5.1_12$estimates["elpd_loo","SE"], loo5.2_12$estimates["elpd_loo","SE"],
                  loo5.3_12$estimates["elpd_loo","SE"], loo5.4_12$estimates["elpd_loo","SE"],
                  loo5.5_12$estimates["elpd_loo","SE"], loo5.6_12$estimates["elpd_loo","SE"],
                  loo5.7_12$estimates["elpd_loo","SE"])
stacking_wts <- stacking_weights(cbind(loo4.2_12$pointwise[, "elpd_loo"],
                                       loo5.1_12$pointwise[, "elpd_loo"], loo5.2_12$pointwise[, "elpd_loo"],
                                       loo5.3_12$pointwise[, "elpd_loo"], loo5.4_12$pointwise[, "elpd_loo"],
                                       loo5.5_12$pointwise[, "elpd_loo"], loo5.6_12$pointwise[, "elpd_loo"],
                                       loo5.7_12$pointwise[, "elpd_loo"]))
Step6_12 <- round(cbind(looic, looic_SE, stacking_wts), 3)
# write.table(Step6_12, file = "D:/Research/OmniBEF_NLA/NLA_HM/ModResults/brms_Omni1/ModSelec_12_step6.csv", sep = ",",
#             row.names = TRUE, col.names = TRUE)
#summary(Mod5.4_12)

### Step 7
fn5.4_12 <- bf(ln_phyDen ~ ln_NTL + ln_PTL + ln_phySRr + ln_zpSRr + Lat + pH)
fn6.1_12 <- bf(ln_phyDen ~ ln_NTL + ln_PTL + ln_phySRr + ln_zpSRr + Lat + pH + ln_zpDen)
fn6.2_12 <- bf(ln_phyDen ~ ln_NTL + ln_PTL + ln_phySRr + ln_zpSRr + Lat + pH + ln_DOC)
fn6.3_12 <- bf(ln_phyDen ~ ln_NTL + ln_PTL + ln_phySRr + ln_zpSRr + Lat + pH + ln_Silica)
fn6.4_12 <- bf(ln_phyDen ~ ln_NTL + ln_PTL + ln_phySRr + ln_zpSRr + Lat + pH + ln_Temp)
fn6.5_12 <- bf(ln_phyDen ~ ln_NTL + ln_PTL + ln_phySRr + ln_zpSRr + Lat + pH + ln_Sol)
Mod5.4_12 <- brm(formula = fn5.4_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod6.1_12 <- brm(formula = fn6.1_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod6.2_12 <- brm(formula = fn6.2_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod6.3_12 <- brm(formula = fn6.3_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod6.4_12 <- brm(formula = fn6.4_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod6.5_12 <- brm(formula = fn6.5_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
loo5.4_12 <- loo(Mod5.4_12, reloo = TRUE)
loo6.1_12 <- loo(Mod6.1_12, reloo = TRUE)
loo6.2_12 <- loo(Mod6.2_12, reloo = TRUE)
loo6.3_12 <- loo(Mod6.3_12, reloo = TRUE)
loo6.4_12 <- loo(Mod6.4_12, reloo = TRUE)
loo6.5_12 <- loo(Mod6.5_12, reloo = TRUE)

looic <- rbind(loo5.4_12$estimates["elpd_loo","Estimate"],
               loo6.1_12$estimates["elpd_loo","Estimate"], loo6.2_12$estimates["elpd_loo","Estimate"],
               loo6.3_12$estimates["elpd_loo","Estimate"], loo6.4_12$estimates["elpd_loo","Estimate"],
               loo6.5_12$estimates["elpd_loo","Estimate"])
looic_SE <- rbind(loo5.4_12$estimates["elpd_loo","SE"],
                  loo6.1_12$estimates["elpd_loo","SE"], loo6.2_12$estimates["elpd_loo","SE"],
                  loo6.3_12$estimates["elpd_loo","SE"], loo6.4_12$estimates["elpd_loo","SE"],
                  loo6.5_12$estimates["elpd_loo","SE"])
stacking_wts <- stacking_weights(cbind(loo5.4_12$pointwise[, "elpd_loo"],
                                       loo6.1_12$pointwise[, "elpd_loo"], loo6.2_12$pointwise[, "elpd_loo"],
                                       loo6.3_12$pointwise[, "elpd_loo"], loo6.4_12$pointwise[, "elpd_loo"],
                                       loo6.5_12$pointwise[, "elpd_loo"]))
Step7_12 <- round(cbind(looic, looic_SE, stacking_wts), 3)
# write.table(Step7_12, file = "D:/Research/OmniBEF_NLA/NLA_HM/ModResults/brms_Omni1/ModSelec_12_step7.csv", sep = ",",
#             row.names = TRUE, col.names = TRUE)
# summary(Mod6.2_12)
# summary(Mod6.3_12)

### Step 8
fn6.3_12 <- bf(ln_phyDen ~ ln_NTL + ln_PTL + ln_phySRr + ln_zpSRr + Lat + pH + ln_Silica)
fn7.1_12 <- bf(ln_phyDen ~ ln_NTL + ln_PTL + ln_phySRr + ln_zpSRr + Lat + pH + ln_Silica + ln_zpDen)
fn7.2_12 <- bf(ln_phyDen ~ ln_NTL + ln_PTL + ln_phySRr + ln_zpSRr + Lat + pH + ln_Silica + ln_DOC)
fn7.3_12 <- bf(ln_phyDen ~ ln_NTL + ln_PTL + ln_phySRr + ln_zpSRr + Lat + pH + ln_Silica + ln_Temp)
fn7.4_12 <- bf(ln_phyDen ~ ln_NTL + ln_PTL + ln_phySRr + ln_zpSRr + Lat + pH + ln_Silica + ln_Sol)
Mod6.3_12 <- brm(formula = fn6.3_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod7.1_12 <- brm(formula = fn7.1_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod7.2_12 <- brm(formula = fn7.2_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod7.3_12 <- brm(formula = fn7.3_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod7.4_12 <- brm(formula = fn7.4_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
loo6.3_12 <- loo(Mod6.3_12, reloo = TRUE)
loo7.1_12 <- loo(Mod7.1_12, reloo = TRUE)
loo7.2_12 <- loo(Mod7.2_12, reloo = TRUE)
loo7.3_12 <- loo(Mod7.3_12, reloo = TRUE)
loo7.4_12 <- loo(Mod7.4_12, reloo = TRUE)

looic <- rbind(loo6.3_12$estimates["elpd_loo","Estimate"],
               loo7.1_12$estimates["elpd_loo","Estimate"], loo7.2_12$estimates["elpd_loo","Estimate"],
               loo7.3_12$estimates["elpd_loo","Estimate"], loo7.4_12$estimates["elpd_loo","Estimate"])
looic_SE <- rbind(loo6.3_12$estimates["elpd_loo","SE"],
                  loo7.1_12$estimates["elpd_loo","SE"], loo7.2_12$estimates["elpd_loo","SE"],
                  loo7.3_12$estimates["elpd_loo","SE"], loo7.4_12$estimates["elpd_loo","SE"])
stacking_wts <- stacking_weights(cbind(loo6.3_12$pointwise[, "elpd_loo"],
                                       loo7.1_12$pointwise[, "elpd_loo"], loo7.2_12$pointwise[, "elpd_loo"],
                                       loo7.3_12$pointwise[, "elpd_loo"], loo7.4_12$pointwise[, "elpd_loo"]))
Step8_12 <- round(cbind(looic, looic_SE, stacking_wts), 3)
# write.table(Step8_12, file = "D:/Research/OmniBEF_NLA/NLA_HM/ModResults/brms_Omni1/ModSelec_12_step8.csv", sep = ",",
#             row.names = TRUE, col.names = TRUE)
# summary(Mod6.3_12)
# summary(Mod7.2_12)
### Step 9
fn7.2_12 <- bf(ln_phyDen ~ ln_NTL + ln_PTL + ln_phySRr + ln_zpSRr + Lat + pH + ln_Silica + ln_DOC)
fn8.1_12 <- bf(ln_phyDen ~ ln_NTL + ln_PTL + ln_phySRr + ln_zpSRr + Lat + pH + ln_Silica + ln_DOC + ln_zpDen)
fn8.2_12 <- bf(ln_phyDen ~ ln_NTL + ln_PTL + ln_phySRr + ln_zpSRr + Lat + pH + ln_Silica + ln_DOC + ln_Temp)
fn8.3_12 <- bf(ln_phyDen ~ ln_NTL + ln_PTL + ln_phySRr + ln_zpSRr + Lat + pH + ln_Silica + ln_DOC + ln_Sol)
Mod7.2_12 <- brm(formula = fn7.2_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod8.1_12 <- brm(formula = fn8.1_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod8.2_12 <- brm(formula = fn8.2_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod8.3_12 <- brm(formula = fn8.3_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
loo7.2_12 <- loo(Mod7.2_12, reloo = TRUE)
loo8.1_12 <- loo(Mod8.1_12, reloo = TRUE)
loo8.2_12 <- loo(Mod8.2_12, reloo = TRUE)
loo8.3_12 <- loo(Mod8.3_12, reloo = TRUE)

looic <- rbind(loo7.2_12$estimates["elpd_loo","Estimate"],
               loo8.1_12$estimates["elpd_loo","Estimate"], loo8.2_12$estimates["elpd_loo","Estimate"],
               loo8.3_12$estimates["elpd_loo","Estimate"])
looic_SE <- rbind(loo7.2_12$estimates["elpd_loo","SE"],
                  loo8.1_12$estimates["elpd_loo","SE"], loo8.2_12$estimates["elpd_loo","SE"],
                  loo8.3_12$estimates["elpd_loo","SE"])
stacking_wts <- stacking_weights(cbind(loo7.2_12$pointwise[, "elpd_loo"],
                                       loo8.1_12$pointwise[, "elpd_loo"], loo8.2_12$pointwise[, "elpd_loo"],
                                       loo8.3_12$pointwise[, "elpd_loo"]))
Step9_12 <- round(cbind(looic, looic_SE, stacking_wts), 3)
# write.table(Step9_12, file = "D:/Research/OmniBEF_NLA/NLA_HM/ModResults/brms_Omni1/ModSelec_12_step9.csv", sep = ",",
#             row.names = TRUE, col.names = TRUE)

# The best model is fn6.3 <- bf(ln_phyDen ~ ln_NTL + ln_PTL + ln_phySRr + ln_zpSRr + Lat + pH + ln_Silica)

##### 2012 #########################################################################################################
##### Finding the best model  ###############################################################################################
#############################################################################################################################
