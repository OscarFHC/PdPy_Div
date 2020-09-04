if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pasilla")

pasillaCountTable = read.table( datafile, header=TRUE, row.names=1 )

pasillaDesign = data.frame(
  row.names = colnames( pasillaCountTable ),
  condition = c( "untreated", "untreated", "untreated",
                 "untreated", "treated", "treated", "treated" ),
  libType = c( "single-end", "single-end", "paired-end",
               "paired-end", "single-end", "paired-end", "paired-end" ) )

pairedSamples = pasillaDesign$libType == "paired-end"
countTable = pasillaCountTable[ , pairedSamples ]
condition = pasillaDesign$condition[ pairedSamples ]
condition1 = rep("trmt", ncol(countTable))
cds = newCountDataSet( countTable, condition )
cds1 = newCountDataSet( countTable, condition )

cds1 <- estimateSizeFactors(cds1)
sizeFactors(cds1)

cds1 <- estimateDispersions(cds1)


