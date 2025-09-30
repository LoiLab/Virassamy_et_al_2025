library(dplyr)
library(tidyr)
library(readr)
library(glue)
library(ggplot2)
library(readxl)
library(betareg)
library(margins)
library(ggeffects)
library(ggh4x)
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
 [1] ggsignif_0.6.4  ggh4x_0.3.0     ggeffects_2.2.1 margins_0.3.28 
 [5] betareg_3.2-1   readxl_1.4.3    ggplot2_3.5.1   glue_1.8.0     
 [9] readr_2.1.5     tidyr_1.3.1     dplyr_1.1.4    
"



# Fig 1b ------------------------------------------------------------------

fig1b.df <- read_excel("Source Data Fig. 1.xlsx", sheet = "Fig. 1b", skip = 1)

# function for formatting p values
p_label_thresh <- function(p, threshold = 1e-4) {
  if (!is.finite(p) || is.na(p)) return("NA")
  if (p <= 0) return('"<" * 10^{-308}')  # renders as < 10^-308
  
  if (p < threshold) {
    e <- floor(log10(p))
    m <- signif(p / 10^e, 2)  # 2 significant digits in mantissa
    m_str <- sub("\\.?0+$", "", format(m, trim = TRUE, scientific = FALSE))
    sprintf("%s %%*%% 10^{%d}", m_str, e)  # e.g., 1.2 × 10^{-5}
  } else {
    signif(p, 2)        # strip trailing dot
  }
}

# plot expression
fig1b.df %>%
  gather(key = cell, value = pct, -case_id, -parity) %>%
  mutate(cell = factor(cell, levels = c("CD45+","CD3+","CD3+CD8+","CD69+CD103-","CD69+CD103+"))) %>%
  ggplot(aes(x = parity, y = pct)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = factor(parity)), width = 0.15, height = 0.15, shape = 21, size = 3, show.legend = F) +
  scale_fill_manual(values = c("#1C78B4", "#F57F1F")) +
  labs(y = "cell type freuqency (%)") +
  ggsignif::geom_signif(comparisons = list(c("N","P")), textsize = 6, margin_top = 0.1,  map_signif_level = p_label_thresh, parse = T) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.title.x = element_blank(),
    strip.background = element_rect(fill = "#D0D2D3", linetype = 0),
    strip.text = element_text(size = 11, face = "bold"),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 14),
    panel.grid = element_blank()
  ) +
  ggh4x::facet_wrap2(~cell, nrow = 1, scales = "free_y", strip = ggh4x::strip_themed(
    background_x = ggh4x::elem_list_rect(fill = c("#D0D2D3", "#D0D2D3", "#D0D2D3", alpha(c("#75AFDB","#F281A5"), 1))))) +
  ggh4x::facetted_pos_scales(
    y = list(
      scale_y_continuous(limits = c(0, 40), breaks = c(0,10,20,30,40)),
      scale_y_continuous(limits = c(0, 18), breaks = c(0,5,10,15)),
      scale_y_continuous(limits = c(0, 10), breaks = c(0,2.5,5,7.5)),
      scale_y_continuous(limits = c(0, 110), breaks = c(0,25,50,75,100)),
      scale_y_continuous(limits = c(0, 110), breaks = c(0,25,50,75,100))
    )
  )



# Fig 1d ------------------------------------------------------------------

fig1d.df <- read_excel("Source Data Fig. 1.xlsx", sheet = "Fig. 1d", skip = 1)


## Beta reg fits

### Prep data frame
facs_parity_lm.df <- fig1d.df %>%
  gather(key = cell, value = pct, `CD69+CD103-`, `CD69+CD103+`) %>%
  mutate(age50 = ifelse(age_tissue_donation < 50, "<50", "≥50")) %>%
  mutate(cell = factor(cell, levels = c("CD69+CD103-","CD69+CD103+"))) %>%
  mutate(frequency = pct/100)

### beta regression model for each T cell subset with age at tissue donation as co-variate 
birth2tissue_CD103neg_br.m <- facs_parity_lm.df %>%
  filter(cell == "CD69+CD103-") %>%
  betareg(frequency ~ last_birth_to_tissue_donation_years + age_tissue_donation, data = .)

birth2tissue_CD103pos_br.m <- facs_parity_lm.df %>%
  filter(cell == "CD69+CD103+") %>%
  betareg(frequency ~ last_birth_to_tissue_donation_years + age_tissue_donation, data = .)


### Margins package for Average Marginal Effect (AME)
#### CD103neg
birth2tissue_CD103neg_br_ames <- margins::margins(birth2tissue_CD103neg_br.m, type = "response", data = birth2tissue_CD103neg_br.m$model)

birth2tissue_CD103neg_br_ames.s <- summary(birth2tissue_CD103neg_br_ames)

ame_value_CD103neg <- birth2tissue_CD103neg_br_ames.s$AME["last_birth_to_tissue_donation_years"]
p_value_CD103neg <- birth2tissue_CD103neg_br_ames.s$p[2]

##### Create annotation
text_label_CD103neg <- glue::glue("AME = {round(100*ame_value_CD103neg, 2)}, p = {signif(p_value_CD103neg, 2)}")


#### CD103pos
birth2tissue_CD103pos_br_ames <- margins::margins(birth2tissue_CD103pos_br.m, type = "response", data = birth2tissue_CD103pos_br.m$model)

birth2tissue_CD103pos_br_ames.s <- summary(birth2tissue_CD103pos_br_ames)

ame_value_CD103pos <- birth2tissue_CD103pos_br_ames.s$AME["last_birth_to_tissue_donation_years"]
p_value_CD103pos <- birth2tissue_CD103pos_br_ames.s$p[2]

##### Create annotation
text_label_CD103pos <- glue::glue("AME = {round(100*ame_value_CD103pos, 2)}, p = {round(p_value_CD103pos, 2)}")


### Generate data for trend line and CIs with ggeffects           
birth2tissue_CD103neg_br.m.gge <- ggpredict(birth2tissue_CD103neg_br.m, terms = "last_birth_to_tissue_donation_years")
birth2tissue_CD103pos_br.m.gge <- ggpredict(birth2tissue_CD103pos_br.m, terms = "last_birth_to_tissue_donation_years")

#### Align to existing facets
gge_facet_bind.df <- facs_parity_lm.df %>%
  dplyr::select(cell) %>%
  unique


br_gge.df <- bind_rows(as.data.frame(birth2tissue_CD103neg_br.m.gge), as.data.frame(birth2tissue_CD103pos_br.m.gge)) %>%
  mutate(cell = c(rep.int("CD69+CD103-", 11),rep.int("CD69+CD103+", 11))) %>%
  left_join(gge_facet_bind.df)


# significance annotation
annotate_lb2t.df <- data.frame(label = c(text_label_CD103neg, text_label_CD103pos),
                               xpos = c(5,5),
                               ypos = c(90, 5)) %>%
  bind_cols(gge_facet_bind.df)


# plot expression
facs_parity_lm.df %>%
  ggplot(aes(x = last_birth_to_tissue_donation_years, y = pct, colour = age50, fill = age50)) +
  geom_point(shape = 21, show.legend = T, size = 2, colour = "black") +
  scale_colour_manual(values = c("#479A8B", "#AE59A1")) +
  scale_fill_manual(values = c("#479A8B", "#AE59A1")) +
  scale_y_continuous(limits = c(0,100)) +
  geom_line(data = br_gge.df, aes(x = x, y = 100*predicted), color = "blue", inherit.aes = F) +
  geom_ribbon(data = br_gge.df, aes(x = x, ymin = 100*conf.low, ymax = 100*conf.high), alpha = 0.1, fill = "blue", inherit.aes = F) +
  labs(x = "Years from last birth to tissue donation", y = "% of CD3+CD8+ T cells", fill = "Age") +
  geom_text(data = annotate_lb2t.df, aes(x = xpos, y = ypos, label = label), hjust = 0, inherit.aes = F) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", linetype = 0),
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)) +
  ggh4x::facet_wrap2(~cell, 
                     nrow = 1,
                     strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = alpha(c("#75AFDB","#F281A5"), 1))))

