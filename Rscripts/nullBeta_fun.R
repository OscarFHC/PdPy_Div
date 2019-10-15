BDiv_null_test = function(spXsite, 
                          plot_names_in_col1 = TRUE, 
                          classic_metric = FALSE, 
                          split_ties = TRUE, 
                          reps = 99, 
                          set_all_species_equal = FALSE, 
                          as.distance.matrix = TRUE, 
                          report_similarity = FALSE){
  
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
  n_sites<-nrow(spXsite)
  gamma<-ncol(spXsite)
  
  ##make the spXsite matrix into a pres/abs. (overwrites initial spXsite matrix):
  spXsite_p <- ceiling(spXsite/max(spXsite))
  
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur<-apply(spXsite, MARGIN=2, FUN=sum)
  
  
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
  col_count <- 1
  
  ## null_array will hold the actual null distribution values.  
  ## Each element of the array corresponds to a null distribution for each combination of alpha values.  
  ## The alpha_table is used to point to the correct null distribution- the row numbers of alpha_table correspond to the [[x]] indices of the null_array.  
  ## Later the function will find the row of alpha_table with the right combination of alpha values.  
  ## That row number is used to identify the element of null_array that contains the correct null distribution for that combination of alpha levels. 
  
  
  null_array<-list()
  
  ##looping over each combination of alpha levels:
  
  for(a1 in 1:length(alpha_levels)){
    for(a2 in a1:length(alpha_levels)){
      
      ##build a null distribution of the number of shared species for a pair of alpha values:
      null_shared_spp<-NULL
      for(i in 1:reps){
        
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
        null_shared_spp[i] <- vegdist(rbind(com1, com2)[,comU], method = "chao")
        #null_shared_spp[i]<-sum((com1+com2)>1)
        
      }
      
      ##store null distribution, record values for alpha 1 and 2 in the alpha_table to help find the correct null distribution later:
      null_array[[col_count]] <- null_shared_spp
      
      alpha_table[col_count, which(names(alpha_table) == "smaller_alpha")] <- names(alpha_levels[a1])
      alpha_table[col_count, which(names(alpha_table) == "bigger_alpha")] <- names( alpha_levels[a2])
      
      #increment the counter for the columns of the alpha table/ elements of the null array
      col_count <- col_count + 1
    }
    
  }
  
  ##create a new column with both alpha levels to match on:
  alpha_table$matching <- paste(alpha_table[,1], alpha_table[,2], sep = "_")
  
  
  #####################
  ##do the test:
  
  
  
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  pct_less <- matrix(data = NA, nrow = n_sites, ncol = n_sites, dimnames = list(row.names(spXsite), row.names(spXsite)))
  dist_obs <- matrix(data = NA, nrow = n_sites, ncol = n_sites, dimnames = list(row.names(spXsite), row.names(spXsite)))
  
  ##for each pair of sites (duplicates effort now to make a full matrix instead of a half one- but this part should be minimal time as compared to the null model building)
  for(i in 1:n_sites){
    for(j in 1:n_sites){
      
      ##how many species are shared between the two sites:
      comU <- as.numeric(which(spXsite[i, ] + spXsite[j, ] != 0))
      n_shared_obs <- vegdist(rbind(as.numeric(spXsite[i, ]), as.numeric(spXsite[j, ]))[, comU], method = "chao")
      #sum((spXsite[i, ] + spXsite[j,])>1)
      
      null_index <- ifelse(paste(row.names(spXsite[i,]), row.names(spXsite[j,]), sep = "_") %in% alpha_table$matching,
                           which(alpha_table$matching == paste(row.names(spXsite[i,]), row.names(spXsite[j,]), sep = "_")),
                           which(alpha_table$matching == paste(row.names(spXsite[j,]), row.names(spXsite[i,]), sep = "_")))
      
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null <- sum(null_array[[null_index]] == n_shared_obs)
      
      ##how many null values are bigger than the observed value?
      num_greater_in_null <- sum(null_array[[null_index]] < n_shared_obs)
      
      rc <- (num_greater_in_null)/reps
      
      if(split_ties){
        rc <- ((num_greater_in_null + (num_exact_matching_in_null)/2)/reps)
      }
      
      if(!classic_metric){
        ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
        rc <- (rc - 0.5) * 2
      }
      
      ## at this point rc represents an index of dissimilarity- multiply by -1 to convert to a similarity as specified in the original 1979 Raup Crick paper
      if(report_similarity & !classic_metric){
        rc <- rc * -1
      }
      
      ## the switch to similarity is done differently if the original 0 to 1 range of the metric is used:
      if(report_similarity & classic_metric){
        rc <- 1 - rc
      }
      
      ##store the metric in the results matrix:
      pct_less[i,j] <- rc
      dist_obs[i,j] <- n_shared_obs
      
    }
  }
  
  if(as.distance.matrix){
    pct_less <- as.dist(pct_less)
  }	
  
  return(list(dist_obs = dist_obs, pct_less = pct_less))
}

spXsite <- read.table(file = "D:/Research/PdPy_Div/data/16s_OTUtab.csv", 
                      sep = ",", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)

ini <- Sys.time()
BDiv_null_test(spXsite = spXsite, 
               plot_names_in_col1 = TRUE, 
               classic_metric = TRUE, 
               split_ties = TRUE, 
               reps = 1, 
               set_all_species_equal = FALSE, 
               as.distance.matrix = FALSE, 
               report_similarity = FALSE)
Sys.time() - ini


Bac_phylo<- read.tree(file = "D:/Research/PdPy_Div/data/treeNJ_16s.tree")
Bac_comm <- t(read.table(file = "D:/Research/PdPy_Div/data/16s_seqtab.csv", sep = ",", 
                         header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE))
Bac_ra_comm <- Bac_comm / rowSums(Bac_comm)

pd(Bac_comm, root(Bac_phylo, 1, resolve.root = TRUE))[,1]
UniFrac <- GUniFrac(Bac_comm, root(Bac_phylo, 1, resolve.root = TRUE), alpha = c(0, 1))
#UniFrac[, , 1]

str(UniFrac$unifracs)


#plot(Bac_phylo)
str(Bac_comm)
data(phylocom)
cophenetic(phylocom$phylo)

phydist <- cophenetic(Bac_phylo)
phyvcv <- vcv(Bac_phylo)
ses.mpd.result <- ses.mpd(Bac_comm, phydist, null.model = "taxa.labels",
                          abundance.weighted = FALSE, runs = 99)



samp <- read.table(file = "C:/Users/Oscar/Downloads/beta.example.sample.txt", sep = "\t", row.names = 1, header = TRUE)
ra.samp <- samp / rowSums(samp)
phylo <- read.tree(file = "C:/Users/Oscar/Downloads/beta.example.phylo.txt")

plot(phylo)
pd(as.matrix(samp[1,]), phylo)
unifrac(samp, phylo)
pd(as.matrix(phylocom$sample[1,]), phylocom$phylo, include.root = TRUE)

pd(Bac_comm, Bac_phylo)

Bac_phylo$edge
Bac_phylo$edge.length

node <- matrix(NA, nrow = length(Bac_phylo$edge.length), ncol = 6)
node[,1:2] <- Bac_phylo$edge
node[,3] <- Bac_phylo$edge.length

for(i in 1:dim(node)[1]){
  leave.node <- tips(Bac_phylo, node[i, 2])  
  node[i, 4] <- sum(Bac_ra_comm[1, leave.node])
  node[i, 5] <- sum(Bac_ra_comm[2, leave.node])
  node[i, 6] <- node[i, 3] * (abs(node[i, 4] - node[i, 5]))
}
head(node)
sum(node[, 6])

Bac_phylo_rt <- root(Bac_phylo, 1)
GUniFrac(Bac_comm, root(Bac_phylo, 2, resolve.root = TRUE))

data(bird.orders)
plot(bird.orders)
plot(root(bird.orders, 3, resolve.root = TRUE))
plot(root(bird.orders, 1:5))
tr <- root(bird.orders, 1)
plot(drop.tip(tr, 1))
is.rooted(drop.tip(tr, 1))


