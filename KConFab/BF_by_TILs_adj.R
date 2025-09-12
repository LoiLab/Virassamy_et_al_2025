library(MASS)
library(dplyr)
library(rstatix)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(emmeans)
library(readxl)
library(ggpubr)

file.p = "/home/llara/PPBC_Court/Parity Fig 5 BF data source-de-identifed.xlsx"
save.p = "/home/llara/PPBC_Court/cleaned_data/plots/"
df = read_xlsx(path = file.p, sheet = 2)#read.table(file = file.p, header = T,sep = "\t" )
df$BF_group_cat2 = df$BF.cat.duration
table(df$BF_group_cat2)
df = df[!df$BF_group_cat2 %in% c("Unknown","Nulliparous", "NA"),]
#df = df[!is.na(df$BF_group_cat2),]
df[df$BF_group_cat2 %in% "6-12",]$BF_group_cat2 = "6-12m"
table(df$BF_group_cat2)

df$BF_group_cat2 = factor(df$BF_group_cat2, levels = c("0","<6","6-12m",">12"))
df$StromalTILS = as.numeric( df$StromalTILS )
df = df[!is.na(df$StromalTILS),]

# Parametric adjust but non parametric testing
model = lm(StromalTILS ~ Age.dx, data = df)
df$TIL.adj = resid(model) + coef(model)[1]
# Re-scale from 0 - 100
df$TIL.adj.scaled <- (df$TIL.adj - min(df$TIL.adj)) / (max(df$TIL.adj) - min(df$TIL.adj)) * 100
KW.stats = df %>% kruskal_test(TIL.adj.scaled ~ BF_group_cat2)
# then
stat.val = df %>% dunn_test(TIL.adj.scaled ~ BF_group_cat2, p.adjust.method = "fdr") %>% 
  add_xy_position(x = "BF_group_cat2", dodge = 0.8, step.increase = 0.15)
stat.val$p.adj = round(stat.val$p.adj,digits = 3)
kw_label <- paste("Kruskal-Wallis, p = ",KW.stats$p) #paste0("Kruskal-Wallis, p = ", round(KW.stats$p,3))
x_offset <- min(df$Age.dx, na.rm = TRUE) + 0.05 * (max(df$Age.dx, na.rm = TRUE) - min(df$Age.dx, na.rm = TRUE))

# Create the plot
p = ggplot(df, aes(x = BF_group_cat2, y = TIL.adj.scaled, color=BF_group_cat2)) +
  geom_jitter(width = 0.07, alpha = 0.3) +
  stat_summary(
    fun.data = function(x) {
      data.frame(
        y = median(x),
        ymin = quantile(x, 0.25),
        ymax = quantile(x, 0.75)
      )
    },
    geom = "errorbar",
    width = 0.2
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.4,
    fatten = 0,
    fun.min = median,
    fun.max = median
  ) +
  labs(
    legend.position = "top",  
    axis.text = element_text(color="black"),
    axis.title.x=element_blank(),
    axis.title = element_text(face="bold"),
    plot.title = element_text(face="bold", hjust = 0.5),
    x = NULL,
    y = "TIL (Adjusted %)",
    color = "Breast Feeding"
  ) +
  scale_color_npg() +
  scale_x_discrete(labels = c("0\n(n=12)",
                              "<6\n(n=32)",
                              "6-12m\n(n=29)",
                              ">12\n(n=63)")) +
  annotate("text", x = 0.5, y = Inf, label = kw_label, 
           hjust = 0, vjust = 1, size = 5, fontface = "italic") +
  theme_pubr()

ggsave(filename = "BF_by_TILs_IQR.pdf", 
       plot = p,
       device = "pdf", 
       path = save.p, 
       dpi=320, height = 4, width = 5)
