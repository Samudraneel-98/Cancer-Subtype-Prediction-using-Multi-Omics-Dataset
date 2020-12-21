for (package in c("TCGAbiolinks", "SummarizedExperiment","preprocessCore"))
{
  if (!require(package, character.only=T, quietly=F))
  {
    BiocManager::install(package)
    library(package, character.only=T)
  }
}
library(TCGAbiolinks)
library("preprocessCore")
setwd('C://Users//KIIT//Desktop//Bioinfo')
load("HNSC_Methylation.RData")

meth_mat = assay(hnsc.met)
dim(meth_mat)
head(meth_mat)[,1:3]
meth_mat <- na.omit(meth_mat)
#meth_mat = assay(hnsc.met)
#meth_mat = meth_mat[complete.cases(meth_mat), ]

meth_mat_norm = normalize.quantiles(meth_mat,copy=TRUE)

head(meth_mat_norm)[,1:3]
meth_mat_norm = as.data.frame(meth_mat_norm)
head(meth_mat_norm)[,1:3]
rownames(meth_mat_norm) = rownames(meth_mat)
colnames(meth_mat_norm) = colnames(meth_mat)
head(meth_mat_norm)[,1:3]

hnsc_met_snames = as.vector(colnames(assay(hnsc.met)))

hnsc_met_cnames = toupper(substr(hnsc_met_snames,1,16))
length(unique(hnsc_met_cnames))

colnames(meth_mat_norm) = hnsc_met_cnames
head(meth_mat_norm)[,1:3]

write.csv(meth_mat_norm, file = "hnsc_meth_mat_norm.csv", row.names = TRUE)

write.table(meth_mat_norm, file = "hnsc_meth_mat_norm.txt", sep = "\t", row.names = TRUE, col.names = NA)


rm(hnsc.met, query.met, meth_mat)
save.image("hnsc_meth_norm.RData")


