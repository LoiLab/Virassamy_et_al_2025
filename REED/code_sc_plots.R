rm(list = ls())
library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)
library(AUCell)
library(GSEABase)
library(viridis)
library(ggprism)
library(scran)
library(harmony)
library(ggpubr)
library(ggstats)


## This data processing was performed using an Institutional HPC with over 256GB of RAM memory. 
# The original h5ad was filtered for immune cell types and transformed to an R object (rds)

##link: https://cellxgene.cziscience.com/collections/48259aa8-f168-4bf5-b797-af8e88da6637

# library(zellkonverter)
# adata = zellkonverter::readH5AD('reed_integrated.h5ad')
# adata_Seurat <- as.Seurat(adata, counts = "X", data = NULL)
# ##Filter sorted and no parity cells
# adata_Seurat = subset(adata_Seurat,subset = parity!='unknown')
# adata_Seurat = subset(adata_Seurat,subset = FACS_status %in%c('live_sorted','not_sorted'))
# adata_Seurat = subset(adata_Seurat,subset = risk_status %in%c('AR','HR-Unk','HR-cUnk'))
# 
# saveRDS(adata_Seurat,'reed_demo_filtered.rds')

dat = readRDS("reed_demo_filtered.rds")

## Load signature set prepared with sets of signatures. 
signatures.set = readRDS('../assets/signature_set.rds')
signatures.set <- GeneSetCollection(signatures.set)
auc_res <- AUCell_run(as.SingleCellExperiment(dat), signatures.set)
AUC_table = getAUC(auc_res)%>%t()
colnames(AUC_table) = paste0('AUC_',colnames(AUC_table))
meta = dat@meta.data %>% cbind(AUC_table)
rownames(meta) = rownames(dat@meta.data)
dat@meta.data = meta

## Filter by original annotation (Kumar and Reed)
mixed_cells = c("CD8-Tem","CD8-Trm",'CD8_Trm','CD8_Tem')
cells = subset(dat,subset = original_celltype %in%mixed_cells)

### Normalise with Harmony
sce.filt = as.SingleCellExperiment(cells)
clusters <- quickCluster(sce.filt, use.ranks="FALSE", method="igraph",assay.type = "counts")

sce.filt <- computeSumFactors(sce.filt, clusters=clusters,assay.type = "counts")
sce.filt <- scater::logNormCounts(sce.filt,assay.type = "counts")

## Run PCA
sce.filt <- scater::runPCA(x = sce.filt)

## Harmony
harmony_embeddings_ds_donor = RunHarmony(sce.filt,group.by.vars = c('dataset','donor_id'),verbose=T)

## Dim plot
dim_plot = Seurat::DimPlot((as.Seurat(harmony_embeddings_ds_donor)), reduction="TSNE", group.by = "level3", cols = c('dodgerblue','red2')) + 
  ggtitle('Harmony Dataset, donor')

## Boxplot
b_plot = subset(as.Seurat(harmony_embeddings_ds_donor))@meta.data %>% 
  ggplot(aes(x = level3, y = AUC_TRM_Normal_breast_ALL)) +
  geom_boxplot(aes(fill=level3)) + ylab("P-TRM AUCell Enrichment") +
  scale_fill_manual(values = c(CD8_Tem = 'dodgerblue',CD8_Trm = 'red2')) + 
  theme_minimal() + 
  stat_compare_means()



