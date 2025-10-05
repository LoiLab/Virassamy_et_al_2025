rm(list = ls())
library(dplyr)
library(data.table)
library(genefu)
library(ggplot2)
library(ggpubr)
library(survminer)
library(ggh4x)
## Human homologs from Table S5, some symbols have changed in the metabric microarray annotation
signature_genes = fread("../assets/TRMvsAll_signature.csv")
## Combined expression of signature genes and clinical data from CBioportal
tcga_tbl = fread('../assets/tcga_data.csv') %>% data.frame(check.names = F)

# Only colour strips in x-direction
strip <- strip_themed(background_x = elem_list_rect(fill = c('red','deeppink','blue','cyan')))
## Calculate signature scores
x = data.frame(probe = signature_genes$Symbol,EntrezGene.ID=signature_genes$Symbol,
               coefficient=signature_genes$logFC)
annot = data.frame(EntrezGene.ID=signature_genes$Symbol)
scores = sig.score(x = x, data = tcga_tbl[,signature_genes$Symbol], annot = annot,signed = F)
tcga_tbl$signature_score = scores$score

ggscatter(data = tcga_tbl%>%filter(PAM50=='Basal'),x='TILS.avg',y = 'signature_score',
              add = "reg.line",  # Add regressin line
              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
              conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 20, label.y = 7) + facet_wrap2(~ PAM50, strip = strip) +
  theme(
    strip.text.x = element_text(
      size = 14, color = "black", face = "bold.italic"
    )) + ylab('P-TRM Enrichment')


