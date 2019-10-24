######## Loading necessary libraries ########################################################################################
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

##### Reading data ##############################################################################################
Dat.raw <- read.table(file = "D:/Research/PdPy_Div/data/16s_OTUtab.csv", 
                      sep = ",", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
##### Reading data ##############################################################################################

##### function for null Beta diversity ##########################################################################
### function using for loop which is much slower #####
# BDiv_null_func <- function(nsim, method, data){
#   
#   Dat <- data[, !(colnames(data) %in% c("CrSt"))]
#   #Dat <- Dat.raw[, !(colnames(Dat.raw) %in% c("CrSt", "size"))]
#   B_perm <- c()
#   Dat.perm <- Dat
#   for (i in 1:nsim){
#     for (j in 1:nrow(Dat)){
#       Dat.perm[j,] <- Dat[j, sample(c(1:ncol(Dat)), ncol(Dat))]
#     }
#     B_perm <- cbind(B_perm, colMeans(as.matrix(vegdist(Dat.perm, method = method, upper = TRUE))))
#   }
#   
#   null_mean <- c()
#   null_lo <- c()
#   null_hi <- c()
#   
#   for (i in 1:nrow(B_perm)){
#     null_mean <- c(null_mean, mean(B_perm[i,]))
#     null_lo <- c(null_lo, sort(B_perm[i,])[nsim*0.025])
#     null_hi <- c(null_hi, sort(B_perm[i,])[nsim*0.975])
#   }
#   null <- list(null_mean = null_mean, null_lo = null_lo, null_hi = null_hi)
#   null
# }
# Bac_BDiv_null <- BDiv_null_func(nsim = 20, method = "chao", data = Dat.raw)
### function using for loop which is much slower #####

##### calculating null Beta distribution with parallel computation #####
numCores <- detectCores()
numCores

cl <- makeCluster(numCores)

clusterEvalQ(cl, {
  library(vegan)
  library(tidyverse)
  library(ggplot2)
  library(stringr)
  Dat.raw <- read.table(file = "D:/Research/PdPy_Div/data/16s_OTUtab.csv", 
                        sep = ",", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
})

BDiv_null_fun = function(x){
  spXsite = Dat.raw
  plot_names_in_col1 = TRUE
  classic_metric = FALSE
  split_ties = TRUE
  set_all_species_equal = FALSE 
  as.distance.matrix = TRUE 
  report_similarity = FALSE
  
  ## expects a species by site matrix for spXsite, with row names for plots, or optionally plots named in column 1.  
  ## By default calculates a modification of the Raup-Crick metric (standardizing the metric to range from -1 to 1 instead of 0 to 1). 
  ## Specifying classic_metric=TRUE instead calculates the original Raup-Crick metric that ranges from 0 to 1. 
  ## The option split_ties (defaults to TRUE) adds half of the number of null observations that are equal to the observed number of shared species to the calculation- this is highly recommended.  
  ## The argument report_similarity defaults to FALSE so the function reports a dissimilarity (which is appropriate as a measure of beta diversity).  
  ## Setting report_similarity=TRUE returns a measure of similarity, as Raup and Crick originally specified.  
  ## If ties are split (as we recommend) the dissimilarity (default) and similarity (set report_similarity=TRUE) calculations can be flipped by multiplying by -1 (for our modification, which ranges from -1 to 1) or by subtracting the metric from 1 (for the classic metric which ranges from 0 to 1). 
  ## If ties are not split (and there are ties between the observed and expected shared number of species) this conversion will not work. 
  ## The argument reps specifies the number of randomizations (a minimum of 999 is recommended- default is 9999).  
  ## set_all_species_equal weights all species equally in the null model instead of weighting species by frequency of occupancy.  
  
  
  ##Note that the choice of how many plots (rows) to include has a real impact on the metric, as species and their occurrence frequencies across the set of plots is used to determine gamma and the frequency with which each species is drawn from the null model	
  
  
  ##this section moves plot names in column 1 (if specified as being present) into the row names of the matrix and drops the column of names
  if(plot_names_in_col1){
    row.names(spXsite) <- spXsite[,1]
    spXsite <- spXsite[,-1]
  }
  
  
  ## count number of sites and total species richness across all plots (gamma)
  n_sites <- nrow(spXsite)
  gamma <- ncol(spXsite)
  
  ##make the spXsite matrix into a pres/abs. (overwrites initial spXsite matrix):
  spXsite_p <- ceiling(spXsite/max(spXsite))
  
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur <- apply(spXsite, MARGIN = 2, FUN = sum)
  
  
  ##NOT recommended- this is a non-trivial change to the metric:
  ##sets all species to occur with equal frequency in the null model
  ##e.g.- discards any occupancy frequency information
  if(set_all_species_equal){
    occur <- rep(1,gamma)
  }
  
  
  ## determine how many unique species richness values are in the dataset
  ##this is used to limit the number of null communities that have to be calculated
  alpha_levels <- sort(apply(spXsite_p, MARGIN=1, FUN=sum))
  #alpha_levels <- apply(spXsite_p, MARGIN=1, FUN=sum)
  
  ##make_null:
  
  ##alpha_table is used as a lookup to help identify which null distribution to use for the tests later.  It contains one row for each combination of alpha richness levels. 
  
  alpha_table <- data.frame(c(NA), c(NA))
  names(alpha_table) <- c("smaller_alpha", "bigger_alpha")
  
  ## null_array will hold the actual null distribution values.  
  ## Each element of the array corresponds to a null distribution for each combination of alpha values.  
  ## The alpha_table is used to point to the correct null distribution- the row numbers of alpha_table correspond to the [[x]] indices of the null_array.  
  ## Later the function will find the row of alpha_table with the right combination of alpha values.  
  ## That row number is used to identify the element of null_array that contains the correct null distribution for that combination of alpha levels. 
  
  
  null_array <- c()
  
  ##looping over each combination of alpha levels:
  ##build a null distribution of the number of shared species for a pair of alpha values:

  for(a1 in 1:length(alpha_levels)){
    for(a2 in a1:length(alpha_levels)){
      ##two empty null communities of size gamma:
      com1<-rep(0,gamma)
      com2<-rep(0,gamma)
      
      row.names(spXsite) == names(alpha_levels[a1])
      
      ##add alpha1 number of species to com1, weighting by species occurrence frequencies:
      com1_index <- which(row.names(spXsite) == names(alpha_levels[a1])) #order(apply(spXsite_p, MARGIN = 1, FUN = sum))[a1]
      com1[sample(1:gamma, alpha_levels[a1], replace = FALSE, prob = occur)] <- 
        as.numeric(spXsite[com1_index, which(spXsite[com1_index, ] != 0)])
      ##same for com2:
      com2_index <- which(row.names(spXsite) == names(alpha_levels[a2]))#order(apply(spXsite_p, MARGIN = 1, FUN = sum))[a2]
      com2[sample(1:gamma, alpha_levels[a2], replace=FALSE, prob=occur)] <- 
        as.numeric(spXsite[com2_index, which(spXsite[com2_index,] != 0)])
      
      comU <- as.numeric(which(com1 + com2 != 0))
      ##how many species are shared in common?
      null_shared_spp <- vegdist(rbind(com1, com2)[,comU], method = "chao")
      #null_shared_spp[i]<-sum((com1+com2)>1)
      
      ##store null distribution, record values for alpha 1 and 2 in the alpha_table to help find the correct null distribution later:
      null_array <- c(null_array, null_shared_spp)
    }
  }
  null_array
}

nsim.list <- sapply(1:1000, list)
ini <- Sys.time()
test <- parLapply(cl, nsim.list, BDiv_null_fun)
Sys.time() - ini

stopCluster(cl)

##### calculating null Beta distribution with parallel computation #####

##### creating index table for null Beta distribution #####
spXsite <- Dat.raw
row.names(spXsite) <- spXsite[,1]
spXsite <- spXsite[,-1]
spXsite_p <- ceiling(spXsite/max(spXsite))
alpha_levels <- sort(apply(spXsite_p, MARGIN = 1, FUN = sum))
alpha_table <- data.frame(c(NA), c(NA), c(NA), c(NA))
names(alpha_table) <- c("smaller_alpha", "bigger_alpha", "matching", "BDiv_obs")
col_count <- 1
for(a1 in 1:length(alpha_levels)){
  for(a2 in a1:length(alpha_levels)){
    alpha_table[col_count, which(names(alpha_table) == "smaller_alpha")] <- names(alpha_levels[a1])
    alpha_table[col_count, which(names(alpha_table) == "bigger_alpha")] <- names( alpha_levels[a2])
    
    com1_index <- which(row.names(spXsite) == names(alpha_levels[a1]))
    com1 <- as.numeric(spXsite[com1_index, ])
    
    com2_index <- which(row.names(spXsite) == names(alpha_levels[a2]))
    com2 <- as.numeric(spXsite[com2_index, ])
    
    comU <- as.numeric(which(com1 + com2 != 0))
    alpha_table[col_count, "BDiv_obs"] <- vegdist(rbind(com1, com2)[,comU], method = "chao")
    col_count <- col_count + 1
  }
}
alpha_table$matching <- paste(alpha_table[,1], alpha_table[,2], sep = "_")
##### creating index table for null Beta distribution #####

Bac_BDiv_null <- data.frame(matrix(unlist(test), ncol = length(test), byrow = FALSE)) %>%
  cbind(alpha_table)
write.table(Bac_BDiv_null, file = "D:/Research/PdPy_Div/Bac_BDiv_null.csv", sep = ",", col.names = TRUE, row.names = FALSE)




### function using parallel calculation which is faster #####
##### function for null Beta diversity ##########################################################################
Dat.num <- t(as.matrix(Dat.raw[,-1])) # making the data into OUT (row) by site (col) format
colnames(Dat.num) <- Dat.raw[, 1]
#iNEXT(Dat.num[which(rowSums(Dat.num) != 0), ], q = 0, datatype = "abundance")

Bac_A <- iNEXT(Dat.num[which(rowSums(Dat.num) != 0), ], 
               q = 0, datatype = "abundance", size = max(colSums(Dat.num)) + 100000)$AsyEst

Bac_Div <-  Bac_A %>% select(Site, Diversity, Estimator) %>% spread(Diversity, Estimator) %>%
  rename(SR = "Species richness", Shannon = "Shannon diversity", Simpson = "Simpson diversity") %>%
  mutate(Beta_chao = colMeans(as.matrix(vegdist(Dat.raw[, !(colnames(Dat.raw) %in% c("CrSt"))], 
                                                method = "chao", upper = TRUE))),
         Beta_chao_null_mean <- Bac_BDiv_null$null_mean,
         Beta_chao_null_lo <- Bac_BDiv_null$null_lo,
         Beta_chao_null_hi <- Bac_BDiv_null$null_hi)
Bac_ADivSE <-  Bac_A %>% select(Site, Diversity, s.e.) %>% spread(Diversity, s.e.)







Dat.list <- list()
for(i in 1:nrow(Dat.raw)){
  Dat.list[[i]] <- sort(as.numeric(Dat.raw[i, -c(1, which(Dat.raw[i,] == 0))]), decreasing = TRUE)
}
names(Dat.list) <- Dat.raw[, 1]

head(Dat.num)
str(Dat.num)




str(Dat.raw[,-1])
for(i in 1:nrow(Dat.raw)){
  iNEXT(sort(as.numeric(Dat.raw[i, -c(1)]), decreasing = TRUE),
        q = 0, datatype = "abundance", size = max(rowSums(Dat.raw[,-1])))
}

i=34# is not working when q=0 and I have no idea what is the problem


Bac_A <- iNEXT(Dat.list, q = 0, datatype = "abundance", size = max(rowSums(Dat.raw[,-1])))$AsyEst
Bac_BDiv_null <- BDiv_null_func(nsim = 5000, method = "chao", data = Dat.raw)

Bac_Div <-  Bac_A %>% select(Site, Diversity, Estimator) %>% spread(Diversity, Estimator) %>%
  rename(SR = "Species richness", Shannon = "Shannon diversity", Simpson = "Simpson diversity") %>%
  mutate(Beta_chao = colMeans(as.matrix(vegdist(Dat.raw[, !(colnames(Dat.raw) %in% c("CrSt"))], 
                                                method = "chao", upper = TRUE))),
         Beta_chao_null_mean <- Bac_BDiv_null$null_mean,
         Beta_chao_null_lo <- Bac_BDiv_null$null_lo,
         Beta_chao_null_hi <- Bac_BDiv_null$null_hi)
Bac_ADivSE <-  Bac_A %>% select(Site, Diversity, s.e.) %>% spread(Diversity, s.e.)

Bac_Div %>% ggplot() + 
  geom_point(aes(x = SR, y = Beta_jac))






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
