rm(list = ls())
library(dplyr)
library(data.table)
library(genefu)
library(ggplot2)
library(ggpubr)
library(survminer)

## Human homologs from Table S5, some symbols have changed in the metabric microarray annotation
signature_genes = fread("../assets/TRMvsAll_signature.csv")
## Combined expression of signature genes and clinical data from CBioportal
metabric_tbl = fread('../assets/metabric_data.csv') %>% data.frame(check.names = F)


## Calculate signature scores
x = data.frame(probe = signature_genes$Symbol,EntrezGene.ID=signature_genes$Symbol,
               coefficient=signature_genes$logFC)
annot = data.frame(EntrezGene.ID=signature_genes$Symbol)
scores = sig.score(x = x, data = metabric_tbl[,signature_genes$Symbol], annot = annot,signed = F)
metabric_tbl$signature_score = scores$score

## Boxplots


metabric_tbl$Pam50 = as.factor(metabric_tbl$Pam50Subtype)

metabric_tbl %>% 
  filter(Pam50%in%c('Basal','LumA',"LumB","Her2")) %>%
  ggplot(aes(x = Pam50, y = signature_score, fill = Pam50)) +
  geom_boxplot() + stat_compare_means()


# Filter Basal samples
metabric_tbl = metabric_tbl %>% filter(Pam50Subtype=='Basal')
metabric_tbl$signature_score_d = 
  ifelse(metabric_tbl$signature_score>median(metabric_tbl$signature_score),
         'high','low')

metabric_tbl$e.dss = ifelse(metabric_tbl$t.dss>10, 0, metabric_tbl$e.dss)
metabric_tbl$e.os = ifelse(metabric_tbl$t.os>10, 0, metabric_tbl$e.os)

## Survival Plots
## DSS
sfit <- survfit(Surv(t.dss, e.dss) ~ signature_score_d, data=metabric_tbl)
coxregre <- coxph(Surv(t.dss,e.dss)~signature_score,data = metabric_tbl)
hr = signif(summary(coxregre)$coef[1,"exp(coef)"],digits=2)
ci.up= signif(summary(coxregre)$conf.int[1,4],digits=2)
ci.down = signif(summary(coxregre)$conf.int[1,3],digits=2)
pval = summary(coxregre)$coefficients[1,5]
pval.txt = ifelse(pval<0.0001,"p ≈ 0",paste("p =",format(signif(pval,digits=1),scientific=F)))
hr.txt = c('HR: 95% CI')
hr.txt = paste(hr.txt,paste(hr,": ",ci.down,"-",ci.up,sep = ""))

ggsurvplot(data = metabric_tbl,
           sfit,                     # survfit object with calculated statistics.
           conf.int = F,
           xlab = "Time (Years)",
           font.x = c(8, "bold", "black"),
           ylab = "Disease Specific Survival",
           linetype = c(1,2),    # customize X axis label.
           break.time.by = 1,     # break X axis in time intervals by 200.
           ggtheme = theme_light(), # customize plot and risk table with a theme.
           tables.theme = theme_cleantable(),
           fontsize = 3,
           risk.table.y.text.col = T,# colour risk table text annotations.
           risk.table.y.text = T,# show bars instead of names in text annotations
           title = "Metabric: DSS Basal P-Trm signature",
           subtitle = paste(hr.txt,pval.txt,sep = ' '),
           font.subtitle = c(8, "bold.italic", "purple"),
           font.title = c(10, "bold", "blue"),
           surv.median.line = "none",  # add the median survival pointer.
           xlim = c(0, 10),
           pval = T


)

## OS
sfit <- survfit(Surv(t.os, e.os) ~ signature_score_d, data=metabric_tbl)
coxregre <- coxph(Surv(t.os,e.os)~signature_score,data = metabric_tbl)
hr = signif(summary(coxregre)$coef[1,"exp(coef)"],digits=2)
ci.up= signif(summary(coxregre)$conf.int[1,4],digits=2)
ci.down = signif(summary(coxregre)$conf.int[1,3],digits=2)
pval = summary(coxregre)$coefficients[1,5]
pval.txt = ifelse(pval<0.0001,"p ≈ 0",paste("p =",format(signif(pval,digits=1),scientific=F)))

hr.txt = c('HR: 95% CI')
hr.txt = paste(hr.txt,paste(hr,": ",ci.down,"-",ci.up,sep = ""))

ggsurvplot(data = metabric_tbl,
           sfit,                     # survfit object with calculated statistics.
           conf.int = F,
           xlab = "Time (Years)",
           font.x = c(8, "bold", "black"),
           ylab = "Overall Survival",
           linetype = c(1,2),    # customize X axis label.
           break.time.by = 1,     # break X axis in time intervals by 200.
           ggtheme = theme_light(), # customize plot and risk table with a theme.
           tables.theme = theme_cleantable(),
           fontsize = 3,
           risk.table.y.text.col = T,# colour risk table text annotations.
           risk.table.y.text = FALSE,# show bars instead of names in text annotations
           title = "Metabric: OS Basal P_Trm signature",
           subtitle = paste(hr.txt,pval.txt,sep = ' '),
           font.subtitle = c(8, "bold.italic", "purple"),
           font.title = c(10, "bold", "blue"),
           surv.median.line = "none",  # add the median survival pointer.
           xlim = c(0, 10),
           pval = T


)
