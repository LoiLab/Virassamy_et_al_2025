library(survival)
library(survminer)
library(readxl)
library(survPen)
library(forestmodel)
library(forestplot)
library(dplyr)
library(ggsci)

#### 
file.p = "/home/llara/PPBC_Court/Parity Fig 5 BF data source-de-identifed.xlsx"
save.p = "/home/llara/PPBC_Court/cleaned_data/plots/"
df = read_excel(path = file.p, sheet = 1)#read.table(file = file.p,header = T, sep="\t" )
####
df$time.last.cat = df$`time.last.cat. Last live birth to BC diagnosis`

# Set Levels
df$BF = factor(df$BF, levels = c("No", "Yes"))
df$time.last.cat = factor(df$time.last.cat , levels = c("<10", ">10"))

univ = coxph(Surv(OS_yrs, Status_15) ~  
               BF,
             data = df  ) 

summary(univ)
forest_model(univ)
cox.zph( univ )

summary(univ)
s <- summary(univ)

HR_BF     <- round(s$conf.int     ["BFYes","exp(coef)"],2)
CI_lower  <- round(s$conf.int     ["BFYes","lower .95"],2)
CI_upper  <- round(s$conf.int     ["BFYes","upper .95"],2)
pvalue_BF <- round(s$coefficients["BFYes","Pr(>|z|)"],3)

p = ggsurvplot(survfit(Surv(OS_yrs, Status_15) ~ BF, data = df), data = df, # BF
               pval = paste0("HR:\nYes vs No: ",HR_BF," (",CI_lower,",",CI_upper,")","\nP=",pvalue_BF), 
               pval.method = F,    # Add p-value &  method name
               #surv.median.line = "hv",            # Add median survival lines
               legend.title = "",               # Change legend titles
               xlab="Time in years",
               ylab="OS",
               risk.table = TRUE,                  # Add No at risk table
               cumevents = TRUE,                   # Add cumulative No of events table
               tables.height = 0.15,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE,               # Hide tables y axis text
               #title = "Kaplan-Meier Unadjusted",
               censor = F
)

# add method to grid.draw
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}
ggsave(filename = "KM_BF_Binary.pdf",
       plot = p,
       device = "pdf", 
       dpi = 320,
       path = save.p,
       height = 6,
       width =8  )

multi.un = coxph(Surv(OS_yrs, Status_15) ~  
                   BF + 
                   Age.dx +
                   time.last.cat +
                   Treatment_chemo + 
                   Treatment_mastectomy,
                 data = df, 
                 robust = T  ) 

cox.zph( multi.un )
pf.mv = forest_model(multi.un)
car::vif(multi.un)

ggsave(filename = "Multivariate.pdf",
       plot = pf.mv,
       device = cairo_pdf, 
       dpi = 320,
       path = save.p,
       height = 5.5,width =7.5  )
