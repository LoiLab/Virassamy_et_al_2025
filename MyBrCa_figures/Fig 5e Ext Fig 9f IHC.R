library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(stringr)
library(betareg)
library(ggsignif)
library(emmeans)
library(broom)

# Session Info ----

"
R version 4.4.1 (2024-06-14)
Platform: aarch64-apple-darwin20
Running under: macOS 15.5

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] emmeans_1.11.0 ggsignif_0.6.4 betareg_3.2-1  dplyr_1.1.4  
 [5] stringr_1.5.1  ggplot2_3.5.1  magrittr_2.0.3 readr_2.1.5   
 [9] tidyr_1.3.1   
"



# Source data for these figures is protected. See data availability statement in main paper regarding how to access it.
# md5 hash of source data file 'Fig 5e Ext Fig 9 IHC MyBrCa protected source data.csv' confirmed to reproduce analysis and figures with code below: d32d2cb58b1ffe7f3cf8bc2a31fba5cd
# Do not change order of columns if loading data with the below readr expression
IHC.data.df <- read_csv("Fig 5e Ext Fig 9 IHC protected source data.csv",
                        col_types = "cccnncnn") %>%
  mutate(parity_bf_6mo_cat = factor(parity_bf_6mo_cat, levels = c("NP", "BF<6mo", "BF>6mo")))

# prepare breast feeding data
IHC.data_bf.df <- IHC.data.df %>%
  filter(!is.na(parity_bf_6mo_cat))


# Functions ----

model_betareg_funcdf <- function(subtype, idata, response, covars, pvar) {
  
  lm.df <- idata %>%
    filter(PAM50 == subtype) %>%
    dplyr::select(RNAseq_ID, PAM50, all_of(c(pvar, covars, response))) %>%
    na.omit
  
  # lm.res <- betareg(as.formula(paste(response, "~", paste(c(pvar, covars), collapse = "+"))), data = lm.df)
  lm.res <- betareg(reformulate(termlabels = c(pvar, covars), response = response), data = lm.df)
  
  emm_av <- emmeans(lm.res, specs = reformulate(pvar), type = "response")
  
  emm_ci <- confint(pairs(emm_av))
  
  return(list(data = lm.df,
              lmres = lm.res,
              tidylmres = broom::tidy(lm.res, conf.int = T),
              lmci = confint(lm.res),
              emm_av = emm_av,
              emm_res = list(mean_diff = -emm_ci$estimate, CI95 = c(-emm_ci$asymp.UCL, -emm_ci$asymp.LCL))))
}


sig_text_brgen <- function(model, term = c("live_birthP"), p_only = FALSE, addci = FALSE) {
  
  x <- model$tidylmres
  
  tindex <- which(x$term == term) - 1
  
  ci <- model$emm_res$CI95[(tindex + (tindex-1)):(2*tindex)]
  

  if(p_only) {
    
    paste("p=", signif(x[x$term == term,]$p.value, 2), sep = "")
    
  } else {
  
    if(addci) {
  paste("p=", signif(x[x$term == term,]$p.value, 2), ", AME=", 100*signif(model$emm_res$mean_diff[tindex], 2), "%(95%CI", paste(round(ci, 3)*100, collapse = "-"), ")", sep = "") 
    
      } else {
      
      paste("p=", signif(x[x$term == term,]$p.value, 2), ", AME=", 100*signif(model$emm_res$mean_diff[tindex], 2), "%", sep = "")
      
    }
    }
}


# CD8 density by parity status ------------------------------------------------------

## beta regression ----

# calculate N
IHC.data_CD8pop_N.df <- IHC.data.df %>% 
  dplyr::select(RNAseq_ID, live_birth, CD8_pos_area_pct, PAM50) %>%
  group_by(live_birth, PAM50) %>%
  tally %>%
  spread(PAM50, n)


basal_br_cd8pp.l <- model_betareg_funcdf(subtype = "Basal", idata = IHC.data.df, response = "CD8_pos_area_pct", pvar = "live_birth",  covars = c("AgeDiag1", "Grade1"))

her2_br_cd8pp.l <- model_betareg_funcdf(subtype = "Her2", idata = IHC.data.df, response = "CD8_pos_area_pct", pvar = "live_birth",  covars = c("AgeDiag1", "Grade1"))

lumA_br_cd8pp.l <- model_betareg_funcdf(subtype = "LumA", idata = IHC.data.df, response = "CD8_pos_area_pct", pvar = "live_birth",  covars = c("AgeDiag1", "Grade1"))

lumB_br_cd8pp.l <- model_betareg_funcdf(subtype = "LumB", idata = IHC.data.df, response = "CD8_pos_area_pct", pvar = "live_birth",  covars = c("AgeDiag1", "Grade1"))

PAM50_br_cd8pp.df <- rbind(basal_br_cd8pp.l$data, her2_br_cd8pp.l$data, lumA_br_cd8pp.l$data, lumB_br_cd8pp.l$data)


## prepare annotations ----

sig_text_alphabetical_cd8pp <- c(Basal = sig_text_brgen(basal_br_cd8pp.l), Her2 = sig_text_brgen(her2_br_cd8pp.l), LumA = sig_text_brgen(lumA_br_cd8pp.l), LumB = sig_text_brgen(lumB_br_cd8pp.l)) 


PAM50_br_cd8pp_N.df <- PAM50_br_cd8pp.df %>%
  group_by(live_birth, PAM50) %>%
  tally %>%
  mutate(live_birth_n = paste(live_birth, "\n(n = ", n, ")", sep = ""))

annotation_cd8pp.df <- data.frame(
  PAM50 = c("Basal", "Her2","LumA", "LumB"),
  start = PAM50_br_cd8pp_N.df$live_birth_n[1:4],
  end = PAM50_br_cd8pp_N.df$live_birth_n[5:8], 
  y = rep(0.3, 4),
  # y = c(0.25, 0.2, 0.6, 0.3),
  label = sig_text_alphabetical_cd8pp)


## Plot Fig 5e middle - CD8 density Basal subtype ----

annotation_cd8pp_basal.df <- annotation_cd8pp.df %>%
  filter(PAM50 == "Basal") %>%
  mutate(y = c(0.2)) 

cd8pp_basal.p <- IHC.data.df %>%
  left_join(PAM50_br_cd8pp.df,.) %>%
  left_join(PAM50_br_cd8pp_N.df) %>%
  filter(PAM50 == "Basal") %>%
  ggplot(aes(x = factor(live_birth_n), y = CD8_pos_area_pct)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = factor(live_birth)), width = 0.15, height = 0, shape = 21, size = 2, show.legend = F) +
  scale_fill_manual(values = c("#1C78B4", "#F57F1F")) +
  scale_y_continuous(breaks = c(0,0.1,0.2, 0.3), limits = c(0, 0.28), labels = scales::label_percent(suffix = "")) +
  geom_signif(comparisons = list(c(annotation_cd8pp_basal.df$start, annotation_cd8pp_basal.df$end)), textsize = 5, vjust = -0.2, annotations = annotation_cd8pp_basal.df$label) +
  labs(x = "", y = "Intratumoral CD8 density (% pixels)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title = element_text(size = 16),
    strip.background = element_rect(fill = "grey90", linetype = 0),
    strip.text = element_text(size = 11, face = "bold"),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 14),
    panel.grid = element_blank()
  )


cd8pp_basal.p + facet_wrap(~PAM50) +
  theme(strip.background = element_rect(fill = "#EE3323", color = "black"),
        strip.text = element_text(colour = "white"))



# CD8 density by breastfeeding status ----------------------------------------------

## beta regression ----

# calculate N
IHC.data_CD8pop_BFN.df <- IHC.data_bf.df %>% 
  dplyr::select(RNAseq_ID, parity_bf_6mo_cat, CD8_pos_area_pct, PAM50) %>%
  group_by(parity_bf_6mo_cat, PAM50) %>%
  tally %>%
  spread(PAM50, n)


basal_br_cd8pp_bf.l <- model_betareg_funcdf(subtype = "Basal", idata = IHC.data_bf.df, response = "CD8_pos_area_pct", pvar = "parity_bf_6mo_cat", covars = c("AgeDiag1", "Grade1"))

her2_br_cd8pp_bf.l <- model_betareg_funcdf(subtype = "Her2", idata = IHC.data_bf.df, response = "CD8_pos_area_pct", pvar = "parity_bf_6mo_cat", covars = c("AgeDiag1", "Grade1"))

lumA_br_cd8pp_bf.l <- model_betareg_funcdf(subtype = "LumA", idata = IHC.data_bf.df, response = "CD8_pos_area_pct", pvar = "parity_bf_6mo_cat", covars = c("AgeDiag1", "Grade1"))

lumB_br_cd8pp_bf.l <- model_betareg_funcdf(subtype = "LumB", idata = IHC.data_bf.df, response = "CD8_pos_area_pct", pvar = "parity_bf_6mo_cat", covars = c("AgeDiag1", "Grade1"))

PAM50_br_cd8pp_bf.df <- rbind(basal_br_cd8pp_bf.l$data, her2_br_cd8pp_bf.l$data, lumA_br_cd8pp_bf.l$data, lumB_br_cd8pp_bf.l$data)


## prepare annotations ----

sig_text_alphabetical_cd8pp_bfgt6m <- c(Basal = sig_text_brgen(basal_br_cd8pp_bf.l, term = "parity_bf_6mo_catBF>6mo"), Her2 = sig_text_brgen(her2_br_cd8pp_bf.l, term = "parity_bf_6mo_catBF>6mo"), LumA = sig_text_brgen(lumA_br_cd8pp_bf.l, term = "parity_bf_6mo_catBF>6mo"), LumB = sig_text_brgen(lumB_br_cd8pp_bf.l, term = "parity_bf_6mo_catBF>6mo")) 

sig_text_alphabetical_cd8pp_bflt6m <- c(Basal = sig_text_brgen(basal_br_cd8pp_bf.l, term = "parity_bf_6mo_catBF<6mo"), Her2 = sig_text_brgen(her2_br_cd8pp_bf.l, term = "parity_bf_6mo_catBF<6mo"), LumA = sig_text_brgen(lumA_br_cd8pp_bf.l, term = "parity_bf_6mo_catBF<6mo"), LumB = sig_text_brgen(lumB_br_cd8pp_bf.l, term = "parity_bf_6mo_catBF<6mo")) 


PAM50_br_cd8pp_bfN.df <- PAM50_br_cd8pp_bf.df %>%
  group_by(parity_bf_6mo_cat, PAM50) %>%
  tally %>%
  mutate(parity_bf_6mo_cat_n = paste(parity_bf_6mo_cat, "\n(n = ", n, ")", sep = "")) %>%
  mutate(parity_bf_6mo_cat_n = factor(parity_bf_6mo_cat_n, levels = parity_bf_6mo_cat_n))


annotation_cd8pp_bf1.df <- data.frame(
  PAM50 = c("Basal", "Her2","LumA", "LumB"),
  start = PAM50_br_cd8pp_bfN.df$parity_bf_6mo_cat_n[1:4],
  end = PAM50_br_cd8pp_bfN.df$parity_bf_6mo_cat_n[5:8], 
  y = rep(0.3, 4),
  label = sig_text_alphabetical_cd8pp_bflt6m)

annotation_cd8pp_bf2.df <- data.frame(
  PAM50 = c("Basal", "Her2","LumA", "LumB"),
  start = PAM50_br_cd8pp_bfN.df$parity_bf_6mo_cat_n[1:4],
  end = PAM50_br_cd8pp_bfN.df$parity_bf_6mo_cat_n[9:12], 
  y = rep(0.39, 4),
  label = sig_text_alphabetical_cd8pp_bfgt6m)

annotation_cd8pp_bf.df <- bind_rows(annotation_cd8pp_bf1.df, annotation_cd8pp_bf2.df)



## Plot Fig 5e right - CD8 density in breastfeeding Basal subtype ----

annotation_cd8pp_bf_basal.df <- annotation_cd8pp_bf.df %>%
  filter(PAM50 == "Basal") %>%
  mutate(y = c(0.2, 0.25)) 
  
cd8pp_bf_basal.p <- IHC.data_bf.df %>%
  left_join(PAM50_br_cd8pp_bfN.df) %>%
  filter(PAM50 == "Basal") %>%
  ggplot(aes(x = factor(parity_bf_6mo_cat_n), y = CD8_pos_area_pct)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = factor(parity_bf_6mo_cat)), width = 0.15, height = 0, shape = 21, size = 2, show.legend = F) +
  scale_fill_manual(values = c("#1C78B4", "#F69491", "#F57F1F")) +
  scale_y_continuous(breaks = c(0,0.1,0.2, 0.3, 0.4), limits = c(0, 0.28), labels = scales::label_percent(suffix = "")) +
  geom_signif(data = annotation_cd8pp_bf_basal.df, aes(xmin = start, xmax = end, annotations = label, y_position = y), textsize = 5, vjust = -0.5, manual = TRUE) +
  labs(x = "", y = "Intratumoral CD8 density (% pixels)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title = element_text(size = 16),
    strip.background = element_rect(fill = "grey90", linetype = 0),
    strip.text = element_text(size = 11, face = "bold"),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 14),
    panel.grid = element_blank()
  )

cd8pp_bf_basal.p + facet_wrap(~PAM50) +
  theme(strip.background = element_rect(fill = "#EE3323", color = "black"),
        strip.text = element_text(colour = "white"))


## Plot Ext Fig 9f - CD8 density in breastfeeding non-Basal subtypes ----

annotation_cd8pp_bf_nonbasal.df <- annotation_cd8pp_bf.df %>%
  filter(PAM50 != "Basal") 


cd8pp_bf_nonbasal.p <- IHC.data_bf.df %>%
  left_join(PAM50_br_cd8pp_bfN.df) %>%
  filter(PAM50 != "Basal") %>%
  ggplot(aes(x = factor(parity_bf_6mo_cat_n), y = CD8_pos_area_pct)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = factor(parity_bf_6mo_cat)), width = 0.15, height = 0, shape = 21, size = 2, show.legend = F) +
  scale_fill_manual(values = c("#1C78B4", "#F69491", "#F57F1F")) +
  scale_y_continuous(breaks = c(0,0.1,0.2, 0.3, 0.4), limits = c(0, 0.45), labels = scales::label_percent(suffix = "")) +
  facet_wrap(~PAM50, nrow = 1, scales = "free_x") +
  geom_signif(data = annotation_cd8pp_bf_nonbasal.df, aes(xmin = start, xmax = end, annotations = label, y_position = y), textsize = 5, vjust = -0.5, manual = TRUE) +
  labs(x = "", y = "Intratumoral CD8 density (% pixels)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title = element_text(size = 16),
    strip.background = element_rect(fill = "grey90", linetype = 0),
    strip.text = element_text(size = 11, face = "bold", colour = "white"),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 14),
    panel.grid = element_blank()
  )

cd8pp_bf_nonbasal.pg <- ggplot_gtable(ggplot_build(cd8pp_bf_nonbasal.p))
stripr <- which(grepl('strip-t', cd8pp_bf_nonbasal.pg$layout$name))
fills <- c("#F499C1","#1F449C", "#56C9EF")
alphas <- c(1, 1, 1)
k <- 1
for (i in stripr[1:3]) {
  j1 <- which(grepl('rect', cd8pp_bf_nonbasal.pg$grobs[[i]]$grobs[[1]]$childrenOrder))
  j2 <- which(grepl('rect', cd8pp_bf_nonbasal.pg$grobs[[i]]$grobs[[1]]$childrenOrder))
  j3 <- which(grepl('rect', cd8pp_bf_nonbasal.pg$grobs[[i]]$grobs[[1]]$childrenOrder))
  cd8pp_bf_nonbasal.pg$grobs[[i]]$grobs[[1]]$children[[j1]]$gp$fill <- fills[k]
  cd8pp_bf_nonbasal.pg$grobs[[i]]$grobs[[1]]$children[[j2]]$gp$fill <- fills[k]
  cd8pp_bf_nonbasal.pg$grobs[[i]]$grobs[[1]]$children[[j3]]$gp$fill <- fills[k]
  cd8pp_bf_nonbasal.pg$grobs[[i]]$grobs[[1]]$children[[j1]]$gp$alpha <- alphas[k]
  cd8pp_bf_nonbasal.pg$grobs[[i]]$grobs[[1]]$children[[j2]]$gp$alpha <- alphas[k]
  cd8pp_bf_nonbasal.pg$grobs[[i]]$grobs[[1]]$children[[j3]]$gp$alpha <- alphas[k]
  k <- k+1
}


plot(cd8pp_bf_nonbasal.pg)



# CD3 density Parity ------------------------------------------------------

## beta regression ----

# calculate N
IHC.data_CD3pop_N.df <- IHC.data.df %>% 
  dplyr::select(RNAseq_ID, live_birth, CD3_pos_area_pct, PAM50) %>%
  group_by(live_birth, PAM50) %>%
  tally %>%
  spread(PAM50, n)


basal_br_cd3pp.l <- model_betareg_funcdf(subtype = "Basal", idata = IHC.data.df, response = "CD3_pos_area_pct", pvar = "live_birth",  covars = c("AgeDiag1", "Grade1"))

her2_br_cd3pp.l <- model_betareg_funcdf(subtype = "Her2", idata = IHC.data.df, response = "CD3_pos_area_pct", pvar = "live_birth",  covars = c("AgeDiag1", "Grade1"))

lumA_br_cd3pp.l <- model_betareg_funcdf(subtype = "LumA", idata = IHC.data.df, response = "CD3_pos_area_pct", pvar = "live_birth",  covars = c("AgeDiag1", "Grade1"))

lumB_br_cd3pp.l <- model_betareg_funcdf(subtype = "LumB", idata = IHC.data.df, response = "CD3_pos_area_pct", pvar = "live_birth",  covars = c("AgeDiag1", "Grade1"))

PAM50_br_cd3pp.df <- rbind(basal_br_cd3pp.l$data, her2_br_cd3pp.l$data, lumA_br_cd3pp.l$data, lumB_br_cd3pp.l$data)

## prepare annotations ----

sig_text_alphabetical_cd3pp <- c(Basal = sig_text_brgen(basal_br_cd3pp.l), Her2 = sig_text_brgen(her2_br_cd3pp.l), LumA = sig_text_brgen(lumA_br_cd3pp.l), LumB = sig_text_brgen(lumB_br_cd3pp.l)) 


PAM50_br_cd3pp_N.df <- PAM50_br_cd3pp.df %>%
  group_by(live_birth, PAM50) %>%
  tally %>%
  mutate(live_birth_n = paste(live_birth, "\n(n = ", n, ")", sep = ""))


annotation_cd3pp.df <- data.frame(
  PAM50 = c("Basal", "Her2","LumA", "LumB"),
  start = PAM50_br_cd3pp_N.df$live_birth_n[1:4],
  end = PAM50_br_cd3pp_N.df$live_birth_n[5:8], 
  y = rep(0.45, 4),
  label = sig_text_alphabetical_cd3pp)

annotation_cd3pp_basal.df <- annotation_cd3pp.df %>%
  filter(PAM50 == "Basal") %>%
  mutate(y = c(0.45)) 


## Plot Fig 5e left - CD3 density Basal subtype ----

cd3pp_basal.p <- IHC.data.df %>%
  left_join(PAM50_br_cd3pp_N.df) %>%
  filter(PAM50 == "Basal") %>%
  ggplot(aes(x = factor(live_birth_n), y = CD3_pos_area_pct)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = factor(live_birth)), width = 0.15, height = 0, shape = 21, size = 2, show.legend = F) +
  scale_fill_manual(values = c("#1C78B4", "#F57F1F")) +
  scale_y_continuous(breaks = c(0,0.1,0.2, 0.3, 0.4, 0.5), limits = c(0, 0.5), labels = scales::label_percent(suffix = "")) +
  geom_signif(data = annotation_cd3pp_basal.df, aes(xmin = start, xmax = end, annotations = label, y_position = y), textsize = 5, vjust = -0.5, manual = TRUE) +
  labs(x = "", y = "Intratumoral CD3 density (% pixels)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title = element_text(size = 16),
    strip.background = element_rect(fill = "grey90", linetype = 0),
    strip.text = element_text(size = 11, face = "bold"),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 14),
    panel.grid = element_blank()
  )

cd3pp_basal.p + facet_wrap(~PAM50) +
  theme(strip.background = element_rect(fill = "#EE3323", color = "black"),
        strip.text = element_text(colour = "white"))



