if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
for (package in c("TCGAbiolinks", "SummarizedExperiment"))
{
  if (!require(package, character.only=T, quietly=F)) 
  {
    BiocManager::install(package)
    library(package, character.only=T)
  }
}
library(TCGAbiolinks)
setwd('C:/Users/KIIT/Desktop')
load("hnsc_rnaseq_mme_counts.Rdata")

hnsc_mat = assay(hnsc.exp)

dim(hnsc_mat)
head(hnsc_mat)[,1:3]
hnsc_mat = hnsc_mat[complete.cases(hnsc_mat), ]
dim(hnsc_mat)
head(hnsc_mat)[,1:3]
if (!require("DESeq2", character.only=T, quietly=F)) {
  BiocManager::install("DESeq2")
  library("DESeq2")
}
library(DESeq2)
ls()
unique(hnsc.exp$definition)
dim(assay(hnsc.exp))
head(hnsc.exp)[,1:3]
assay(hnsc.exp)[1:3, 1:3]
colnames(colData(hnsc.exp))
unique(hnsc.exp$definition)
hnsc.exp$definition = ifelse(hnsc.exp$definition == "Primary solid Tumor","PrimarySolidTumor", ifelse(hnsc.exp$definition == "Solid Tissue Normal","SolidTissueNormal","Metastatic"))
unique(hnsc.exp$definition)
ddsSE <- DESeqDataSet(hnsc.exp, design = ~ definition)
head(assay(ddsSE))[,1:3]
dim(colData(hnsc.exp))
keep <- rowSums(counts(ddsSE)) >= 10
dim(counts(ddsSE))
dds <- ddsSE[keep,]
dim(counts(dds))
dds$definition <- relevel(dds$definition, ref = "SolidTissueNormal")
vsd <- vst(dds, blind = FALSE)
head(assay(vsd))[,1:3]
vst_hnsc_exp = assay(vsd)
psts_cid = which(hnsc.exp$definition == "PrimarySolidTumor")
length(psts_cid)
hnsc_exp_snames = as.vector(colnames(assay(hnsc.exp)))
hnsc_exp_snames[1:3]
length(unique(hnsc_exp_snames))
hnsc_exp_cnames = toupper(substr(hnsc_exp_snames,1, 16))
hnsc_exp_cnames[1:3]
length(unique(hnsc_exp_cnames))
table(hnsc.exp$definition)
head(vst_hnsc_exp)[,1:3]
dim(vst_hnsc_exp)
vst_hnsc_exp = as.data.frame(vst_hnsc_exp)
head(vst_hnsc_exp)
colnames(vst_hnsc_exp) = hnsc_exp_cnames
head(vst_hnsc_exp)
write.csv(vst_hnsc_exp, file = "hnsc_exp_norm.csv", row.names = TRUE)











