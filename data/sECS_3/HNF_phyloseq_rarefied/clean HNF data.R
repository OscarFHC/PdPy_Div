HNF_raw <- as.data.frame(read.table(file = "D:/Research/PdPy_Div/data/sECS_3/HNF_phyloseq_rarefied/sECS_HNF_seqXst_PR2_3_raw.csv",
                                    sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE))
Bac_comm <- as.data.frame(t(read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_3/sECS_Bac_seqXst_PR2_3.csv",
                                       sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
HNF_index <- which(rownames(Bac_comm) %in% unique(substr(names(HNF_raw), 1, 13)))
HNF_comm <- as.data.frame(matrix(0, nrow(HNF_raw), length(HNF_index)))
for (i in 1: length(HNF_index)){
  pick <- which (substr(names(HNF_raw), 1, 13) == rownames(Bac_comm)[HNF_index[i]])
  print(c(i, substr(names(HNF_raw), 1, 13)[pick]))
  #temp <- ifelse(length(pick)>1, rowMeans(HNF_raw[,pick]), HNF_raw[,pick])
  ifelse(length(pick)>1, temp <- as.data.frame(rowMeans(HNF_raw[,pick])), temp <- as.data.frame(HNF_raw[,pick]))
  HNF_comm[, i] <- temp
}
colnames(HNF_comm) <- rownames(Bac_comm)[HNF_index]
rownames(HNF_comm) <- rownames(HNF_raw)
write.table(HNF_comm, file = "D:/Research/PdPy_Div/data/sECS_3/sECS_HNF_seqXst_PR2_3.csv", 
            sep = ",", col.names = TRUE, row.names = TRUE)
