library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggsignif)
library(broom)
library(ggsignif)

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
[1] stats     graphics  grDevices utils     datasets  methods  
[7] base     

other attached packages:
[1] broom_1.0.7    ggsignif_0.6.4 ggplot2_3.5.1 
[5] readr_2.1.5    tidyr_1.3.1    dplyr_1.1.4   
"


# Source data for this figure is protected. See data availability statement in main paper regarding how to access it.
# md5 hash of source data file 'Fig 5ab MyBrCa protected source data.csv' confirmed to reproduce analysis and figures with code below: 2ba477075bfc560f1b855cf777a9bf0b
# Do not change order of columns if loading data with the below readr expression
mb2lm.df <- read_csv("Fig 5ab MyBrCa protected source data.csv", 
                     col_types = "cncccc") %>%
  mutate(breastfed2 = factor(breastfed2, levels = c("NP","NoBF","BF"))) %>%
  mutate(sequencing_batch = as.character(sequencing_batch))


# Linear model Fig 5a ----

model_funcdf <- function(subtype, idata, response, vars) {
  
  lm.df <- idata %>%
    filter(PAM50 == subtype) %>%
    dplyr::select(RNAseq_ID, PAM50, all_of(c(vars, response))) %>%
    na.omit
  
  lm.res <- lm(as.formula(paste(response, "~", paste(vars, collapse = "+"))), data = lm.df)
  
  return(list(data = lm.df,
              lmres = lm.res,
              tidylmres = broom::tidy(lm.res, conf.int = T),
              lmci = confint(lm.res)))
  }


#### models ----
basal_mb2lm.l <- model_funcdf(subtype = "Basal", idata = mb2lm.df, response = "ESTIMATE", vars = c("live_birth", "sequencing_batch"))

lumA_mb2lm.l <- model_funcdf(subtype = "LumA", idata = mb2lm.df, response = "ESTIMATE", vars = c("live_birth", "sequencing_batch"))

lumB_mb2lm.l <- model_funcdf(subtype = "LumB", idata = mb2lm.df, response = "ESTIMATE", vars = c("live_birth", "sequencing_batch"))

her2_mb2lm.l <- model_funcdf(subtype = "Her2", idata = mb2lm.df, response = "ESTIMATE", vars = c("live_birth", "sequencing_batch"))


# ANOVA Fig 5b ----

anova_funcdf <- function(subtype, idata, group_var, covars, response) {
    
  aov.df <- idata %>%
    filter(PAM50 == subtype) %>%
    dplyr::select(RNAseq_ID, PAM50, all_of(c(group_var, response, covars))) %>%
    na.omit
    
  return(aov.df)
}


anova_res_detail <- function(r, data, response, group_var) {

  tukey.res <- TukeyHSD(r)

  shapiro.res <- shapiro.test(data[[response]])

  return(list(tukey = tukey.res,
              shap = shapiro.res
                ))

  }

#### models ----

# Basal
basal_bf_anova.df <- anova_funcdf(subtype = "Basal", idata = mb2lm.df, group_var = "breastfed2", response = "ESTIMATE", covars = c("sequencing_batch"))

basal_aov.res <- aov(ESTIMATE ~ breastfed2 + sequencing_batch, data = basal_bf_anova.df)

basal_bf2_aov.res.l <- list(anova.res = basal_aov.res,
        data = basal_bf_anova.df,
     tukey = TukeyHSD(basal_aov.res),
     bart = bartlett.res <- bartlett.test(ESTIMATE ~ breastfed2, data = basal_bf_anova.df),
     shap = shapiro.test(basal_bf_anova.df$ESTIMATE),
     shap_res = shapiro.test(residuals(basal_aov.res)))


# Luminal A
lumA_bf_anova.df <- anova_funcdf(subtype = "LumA", idata = mb2lm.df, group_var = "breastfed2", response = "ESTIMATE", covars = c("sequencing_batch"))

lumA_aov.res <- aov(ESTIMATE ~ breastfed2 + sequencing_batch, data = lumA_bf_anova.df)

lumA_bf2_aov.res.l <- list(anova.res = lumA_aov.res,
                            data = lumA_bf_anova.df,
                            tukey = TukeyHSD(lumA_aov.res),
                            bart = bartlett.res <- bartlett.test(ESTIMATE ~ breastfed2, data = lumA_bf_anova.df),
                            shap = shapiro.test(lumA_bf_anova.df$ESTIMATE),
                           shap_res = shapiro.test(residuals(lumA_aov.res)))


# Luminal B
lumB_bf_anova.df <- anova_funcdf(subtype = "LumB", idata = mb2lm.df, group_var = "breastfed2", response = "ESTIMATE", covars = c("sequencing_batch"))

lumB_aov.res <- aov(ESTIMATE ~ breastfed2 + sequencing_batch, data = lumB_bf_anova.df)

lumB_bf2_aov.res.l <- list(anova.res = lumB_aov.res,
                            data = lumB_bf_anova.df,
                            tukey = TukeyHSD(lumB_aov.res),
                            bart = bartlett.res <- bartlett.test(ESTIMATE ~ breastfed2, data = lumB_bf_anova.df),
                            shap = shapiro.test(lumB_bf_anova.df$ESTIMATE),
                           shap_res = shapiro.test(residuals(lumB_aov.res)))


# HER2
her2_bf_anova.df <- anova_funcdf(subtype = "Her2", idata = mb2lm.df, group_var = "breastfed2", response = "ESTIMATE", covars = c("sequencing_batch"))

her2_aov.res <- aov(ESTIMATE ~ breastfed2 + sequencing_batch, data = her2_bf_anova.df)

her2_bf2_aov.res.l <- list(anova.res = her2_aov.res,
                            data = her2_bf_anova.df,
                            tukey = TukeyHSD(her2_aov.res),
                            bart = bartlett.res <- bartlett.test(ESTIMATE ~ breastfed2, data = her2_bf_anova.df),
                            shap = shapiro.test(her2_bf_anova.df$ESTIMATE),
                           shap_res = shapiro.test(residuals(her2_aov.res)))



PAM50_bf_aov.l <- list(basal = basal_bf2_aov.res.l,
                       lumA = lumA_bf2_aov.res.l,
                       lumB = lumB_bf2_aov.res.l,
                       her2 = her2_bf2_aov.res.l)

# Plots ----
### All subtypes Fig 5a ----

# merge subtypes
PAM50_lm.df <- rbind(basal_mb2lm.l$data, lumA_mb2lm.l$data, lumB_mb2lm.l$data, her2_mb2lm.l$data )

# annotation text function
sig_text_gen <- function(model) {
  
  x <- model$tidylmres
  
  ci <- model$lmci

  paste("p = ", signif(x[x$term == "live_birthY",]$p.value, 2), ", β = ", signif(x[x$term == "live_birthY",]$estimate, 3), " (95%CI ", paste(round(ci["live_birthY",], 0), collapse = " - "), ")", sep = "") 
  
}

# annotation text
sig_text_alphabetical <- c(Basal = sig_text_gen(basal_mb2lm.l), Her2 = sig_text_gen(her2_mb2lm.l), LumA = sig_text_gen(lumA_mb2lm.l), LumB = sig_text_gen(lumB_mb2lm.l)) 

# count cases
PAM50_lm_N.df <- PAM50_lm.df %>%
  group_by(live_birth, PAM50) %>%
  tally %>%
  mutate(live_birth = ifelse(live_birth == "Y", "P", "NP")) %>%
  mutate(live_birth_n = paste(live_birth, "\n(n = ", n, ")", sep = ""))

# annotation data frame
annotation_df <- data.frame(
  PAM50 = c("Basal", "Her2","LumA", "LumB"),
  start = PAM50_lm_N.df$live_birth_n[1:4],
  end = PAM50_lm_N.df$live_birth_n[5:8], 
  y = rep(8500, 4),
  label = sig_text_alphabetical)

# generate base plot
p0 <- PAM50_lm.df %>%
  mutate(live_birth = ifelse(live_birth == "Y", "P", "NP")) %>%
  left_join(PAM50_lm_N.df) %>%
  ggplot(aes(x = live_birth_n, y = ESTIMATE)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = live_birth), width = 0.1, height = 0.1, shape = 21, size = 1.5, show.legend = F) +
  scale_fill_manual(values = c("#1C78B4", "#F57F1F")) +
  scale_y_continuous(limits = c(0,9000), breaks = c(0,2000,4000,6000,8000)) +
  facet_wrap(~PAM50, scales = "free_x", nrow = 1) +
  ggsignif::geom_signif(data = annotation_df, aes(xmin = start, xmax = end, annotations = label, y_position = y), textsize = 3, vjust = -0.5, manual = TRUE) +
  labs(x = "", y = "ESTIMATE\nImmune Infiltration Score", title = "Parity and immune infiltration in primary breast cancers") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0),
    axis.title = element_text(size = 16),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    panel.grid = element_blank(),
    strip.text = element_text(size = 14, colour = "white"),
    strip.background = element_rect(colour = "black", linewidth = 1)
  )

# modify facet label colours
p1 <- ggplot_gtable(ggplot_build(p0))
stripr <- which(grepl('strip-t', p1$layout$name))
fills <- c("#EE3323", "#F499C1","#1F449C", "#56C9EF")
alphas <- c(1, 1, 1, 1)
k <- 1
for (i in stripr[1:4]) {
  j1 <- which(grepl('rect', p1$grobs[[i]]$grobs[[1]]$childrenOrder))
  j2 <- which(grepl('rect', p1$grobs[[i]]$grobs[[1]]$childrenOrder))
  j3 <- which(grepl('rect', p1$grobs[[i]]$grobs[[1]]$childrenOrder))
  j4 <- which(grepl('rect', p1$grobs[[i]]$grobs[[1]]$childrenOrder))
  p1$grobs[[i]]$grobs[[1]]$children[[j1]]$gp$fill <- fills[k]
  p1$grobs[[i]]$grobs[[1]]$children[[j2]]$gp$fill <- fills[k]
  p1$grobs[[i]]$grobs[[1]]$children[[j3]]$gp$fill <- fills[k]
  p1$grobs[[i]]$grobs[[1]]$children[[j4]]$gp$fill <- fills[k]
  p1$grobs[[i]]$grobs[[1]]$children[[j1]]$gp$alpha <- alphas[k]
  p1$grobs[[i]]$grobs[[1]]$children[[j2]]$gp$alpha <- alphas[k]
  p1$grobs[[i]]$grobs[[1]]$children[[j3]]$gp$alpha <- alphas[k]
  p1$grobs[[i]]$grobs[[1]]$children[[j4]]$gp$alpha <- alphas[k]
  k <- k+1
}

# final plot
plot(p1)



### Breastfeeding Fig 5b ----

# merge subtypes
PAM50_bf_aov.df <- rbind(PAM50_bf_aov.l$basal$data, PAM50_bf_aov.l$lumA$data, PAM50_bf_aov.l$lumB$data, PAM50_bf_aov.l$her2$data)

# annotation text function
sig_text_aov_bf_func <- function(model) {
  
  x <- signif(model$tukey$breastfed2["BF-NP",],2)
  
  paste("p = ",  x["p adj"], ", Mean Diff = ", x["diff"], "\n(95%CI ", x["lwr"], " - ", x["upr"], ")", sep = "") 
  
}

# annotation text
sig_text_bf_aov <- c(Basal = sig_text_aov_bf_func(PAM50_bf_aov.l$basal), Her2 = sig_text_aov_bf_func(PAM50_bf_aov.l$her2), LumA = sig_text_aov_bf_func(PAM50_bf_aov.l$lumA), LumB = sig_text_aov_bf_func(PAM50_bf_aov.l$lumB)) 

# count cases
PAM50_bf_aov_N.df <- PAM50_bf_aov.df %>%
  group_by(breastfed2, PAM50) %>%
  tally %>%
  mutate(breastfed2_n = paste(breastfed2, "\n(n = ", n, ")", sep = "")) %>%
  mutate(breastfed2_n = factor(breastfed2_n, levels = unique(breastfed2_n)))

# annotation data frame
annotation_bf2_aov.df <- data.frame(
  PAM50 = c("Basal", "Her2","LumA", "LumB"),
  start = PAM50_bf_aov_N.df$breastfed2_n[1:4],
  end = PAM50_bf_aov_N.df$breastfed2_n[9:12], 
  y = rep(8500, 4),
  label = sig_text_bf_aov)

# generate base plot
pbf0 <- PAM50_bf_aov.df %>%
  left_join(PAM50_bf_aov_N.df) %>%
  ggplot(aes(x = breastfed2_n, y = ESTIMATE)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = breastfed2), width = 0.1, height = 0.1, shape = 21, size = 1.5, show.legend = F) +
  scale_fill_manual(values = c("#1C78B4", "#F69491", "#F57F1F")) +
  scale_y_continuous(limits = c(0,9000), breaks = c(0,2000,4000,6000,8000)) +
  facet_wrap(~PAM50, scales = "free_x", nrow = 1) +
  ggsignif::geom_signif(data = annotation_bf2_aov.df, aes(xmin = start, xmax = end, annotations = label, y_position = y), textsize = 3, vjust = 0, manual = TRUE) +
  labs(x = "", y = "ESTIMATE\nImmune Infiltration Score", title = "Breastfeeding and immune infiltration in primary breast cancers") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0),
    axis.title = element_text(size = 16),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    panel.grid = element_blank(),
    strip.text = element_text(size = 14, colour = "white"),
    strip.background = element_rect(colour = "black", linewidth = 1)
  )

# modify facet label colours
pbf1 <- ggplot_gtable(ggplot_build(pbf0))
stripr <- which(grepl('strip-t', pbf1$layout$name))
fills <- c("#EE3323", "#F499C1","#1F449C", "#56C9EF")
alphas <- c(1, 1, 1, 1)
k <- 1
for (i in stripr[1:4]) {
  j1 <- which(grepl('rect', pbf1$grobs[[i]]$grobs[[1]]$childrenOrder))
  j2 <- which(grepl('rect', pbf1$grobs[[i]]$grobs[[1]]$childrenOrder))
  j3 <- which(grepl('rect', pbf1$grobs[[i]]$grobs[[1]]$childrenOrder))
  j4 <- which(grepl('rect', pbf1$grobs[[i]]$grobs[[1]]$childrenOrder))
  pbf1$grobs[[i]]$grobs[[1]]$children[[j1]]$gp$fill <- fills[k]
  pbf1$grobs[[i]]$grobs[[1]]$children[[j2]]$gp$fill <- fills[k]
  pbf1$grobs[[i]]$grobs[[1]]$children[[j3]]$gp$fill <- fills[k]
  pbf1$grobs[[i]]$grobs[[1]]$children[[j4]]$gp$fill <- fills[k]
  pbf1$grobs[[i]]$grobs[[1]]$children[[j1]]$gp$alpha <- alphas[k]
  pbf1$grobs[[i]]$grobs[[1]]$children[[j2]]$gp$alpha <- alphas[k]
  pbf1$grobs[[i]]$grobs[[1]]$children[[j3]]$gp$alpha <- alphas[k]
  pbf1$grobs[[i]]$grobs[[1]]$children[[j4]]$gp$alpha <- alphas[k]
  k <- k+1
}

plot(pbf1)
