# Libraries
require(tidyverse)
require(readr)
require(data.table)
require(TCGAbiolinks)

# Cancer Type
typeCancer = "COAD"
typeTCGA = paste("TCGA", typeCancer, sep ="-")

# Ancestry
ancestryTCGA = read.delim(file = "./data/anc/TCGA_consensus_ancestry.tsv",
                          sep = "\t",
                          header = T) %>% as.data.table()


# Query
query.exp.TCGA.anc.AFR <- GDCquery(project = typeTCGA,
                                   legacy = F,
                                   data.category = "Transcriptome Profiling",
                                   data.type = "Gene Expression Quantification",
                                   workflow.type = "HTSeq - Counts",                  
                                   barcode = ancestryTCGA[consensus_ancestry %in% c("afr_admix","afr")],
                                   sample.type = c("Solid Tissue Normal", "Primary Tumor"))

query.exp.TCGA.anc.EUR <- GDCquery(project = typeTCGA,
                                   legacy = F,
                                   data.category = "Transcriptome Profiling",
                                   data.type = "Gene Expression Quantification",
                                   workflow.type = "HTSeq - Counts",                  
                                   barcode = ancestryTCGA[consensus_ancestry %in% c("eur_admix","eur")],
                                   sample.type = c("Solid Tissue Normal", "Primary Tumor"))

query.exp.TCGA.anc.AFR.fpkm <- GDCquery(project = typeTCGA,
                                   legacy = F,
                                   data.category = "Transcriptome Profiling",
                                   data.type = "Gene Expression Quantification",
                                   workflow.type = "HTSeq - FPKM",                  
                                   barcode = ancestryTCGA[consensus_ancestry %in% c("afr_admix","afr")],
                                   sample.type = c("Solid Tissue Normal", "Primary Tumor"))

query.exp.TCGA.anc.EUR.fpkm <- GDCquery(project = typeTCGA,
                                   legacy = F,
                                   data.category = "Transcriptome Profiling",
                                   data.type = "Gene Expression Quantification",
                                   workflow.type = "HTSeq - FPKM",                  
                                   barcode = ancestryTCGA[consensus_ancestry %in% c("eur_admix","eur")],
                                   sample.type = c("Solid Tissue Normal", "Primary Tumor"))


GDCdownload(query.exp.TCGA.anc.AFR.fpkm)
exp.TCGA.anc.AFR.fpkm <- GDCprepare(query.exp.TCGA.anc.AFR.fpkm)

GDCdownload(query.exp.TCGA.anc.EUR.fpkm)
exp.TCGA.anc.EUR.fpkm <- GDCprepare(query.exp.TCGA.anc.EUR.fpkm)

# why is the sample count the same? they aren't the same samples
colnames(exp.TCGA.anc.AFR.fpkm) %in% ancestryTCGA[consensus_ancestry %in% c("afr_admix","afr")]$patient %>% length()
colnames(exp.TCGA.anc.AFR.fpkm) %in% ancestryTCGA[consensus_ancestry %in% c("eur_admix","eur")]$patient %>% length()

# ensembl to hugo
library(biomaRt)
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")

exp.TCGA.anc.AFR.fpkm.ann <- getBM(attributes = "hgnc_symbol", "ensembl_gene_id", values = rownames(exp.TCGA.anc.AFR.fpkm), mart = mart)
exp.TCGA.anc.EUR.fpkm.ann <- getBM(attributes = "hgnc_symbol", "ensembl_gene_id", values = rownames(exp.TCGA.anc.EUR.fpkm), mart = mart)

rownames(exp.TCGA.anc.AFR.fpkm) = exp.TCGA.anc.AFR.fpkm.ann$hgnc_symbol
rownames(exp.TCGA.anc.EUR.fpkm) = exp.TCGA.anc.EUR.fpkm.ann$hgnc_symbol

write.table(x = exp.TCGA.anc.AFR.fpkm@assays@data@listData, file = "../scratch/Cancer-Health-Disparities-through-genomic-epidemiology-Lens/data/exp.TCGA.anc.AFR.fpkm.txt", sep = "\t", row.names = rownames(exp.TCGA.anc.AFR.fpkm), col.names = colnames(exp.TCGA.anc.AFR.fpkm), quote = F)
write.table(x = exp.TCGA.anc.EUR.fpkm@assays@data@listData, file = "../scratch/Cancer-Health-Disparities-through-genomic-epidemiology-Lens/data/exp.TCGA.anc.EUR.fpkm.txt", sep = "\t", row.names = rownames(exp.TCGA.anc.EUR.fpkm), col.names = colnames(exp.TCGA.anc.EUR.fpkm), quote = F)
