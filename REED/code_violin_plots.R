rm(list = ls())
gc(full = T)
library(dplyr)
library(ggpubr)
library(data.table)
library(rstatix)

# Plots needed:
## Figures
# reed_cd45_epi_ratio
# reed_cd3_epi_ratio
# reed_cd8_epi_ratio
# reed_trm_epi_ratio
## Ext Figures
# reed_macro_epi_ratio
# M1_cells_orig_epi_ratio
# M2_cells_orig_epi_ratio
# monocytes_cells_orig_epi_ratio
# all_dc_cells_orig_epi_ratio
# b_cells_cells_orig_epi_ratio
# orig_fibro_epi_ratio
# kumarfibroMatrix_epi_ratio
# kumarfibroSFRP4_epi_ratio



reed_cd45_cells = c("B_mem_switched","B_mem_unswitched","B_naive","CD4_Th","CD4_naive","CD8_Tc1","CD8_Tem","CD8_Trm","DC", "Macro","Macro-lipo","NK","NKT","Plasma_cell")
reed_cd3_cells = c("CD4_Th","CD4_naive","CD8_Tc1","CD8_Tem","CD8_Trm","NKT")
reed_cd8_cells = c("CD8_Tc1","CD8_Tem","CD8_Trm","NKT")

reed_epi_cells = c("BMYO1","BMYO2","LASP1","LASP2","LASP3","LASP4","LASP5","LHS1","LHS2","LHS3")	
orig_epi_cells = c("BMYO1","BMYO2","Basal","Basal-Myoepithelial","Epithelial","HRpos_Luminal","Luminal1-ALDH1A3","Luminal1-LTF","Luminal2-AREG","Luminal2-MUCL1","Luminal Progenitor","LummHR-SCGB","LummHR-active","LummHR-major","Lumsec-HLA","Lumsec-KIT","Lumsec-basal","Lumsec-lac","Lumsec-major","Lumsec-myo","Lumsec-prol","Mature Luminal","Stroma","LASP1","LASP2","LASP3","LASP4","LASP5","LHS1","LHS2","LHS3")

orig_fibro_cells = c("F1","F2","FB","FB1","FB2","FB3","FB4","Fibro-SFRP4","Fibro-major","Fibro-matrix","Fibro-prematrix","Fibroblast")
reed_macro_cells = c('Macro','Macro-lipo')

reed_fibro_cells = c("FB1","FB2", "FB3","FB4")
orig_M1 = c("Macro-m1","Macro-m1-CCL")
orig_M2 = c("Macro-m2","Macro-m2-CXCL")
orig_monocytes = c('Mono-classical','Mono-non-classical')
orig_all_dc = c('DC','cDC1','cDC2','mDC','pDC')
orig_b_cells = c("B_mem_switched","B_mem_unswitched","B_naive","I4_Bcell","I5_PlasmaCell","Plasma_cell","b_naive","plasma_IgA","plasma_IgG","bmem_switched","bmem_unswitched")
########################################################  

meta.data = readRDS('../assets/REED_metadata.rds')
meta.data = meta.data %>% filter(FACS_status %in% c('live_sorted','not_sorted') & !level2%in%c('Doublet','stripped_nuclei'))

cd45_cells_reed = meta.data %>% filter(level2%in%reed_cd45_cells) %>% group_by(donor_id) %>% summarise(reed_cd45_cells = n()) 
cd3_cells_reed = meta.data %>% filter(level2%in%reed_cd3_cells) %>% group_by(donor_id) %>% summarise(reed_cd3_cells = n()) 
cd8_cells_reed = meta.data %>% filter(level2%in%reed_cd8_cells) %>% group_by(donor_id) %>% summarise(reed_cd8_cells = n()) 
epi_cells_reed = meta.data %>% filter(level2%in%reed_epi_cells) %>% group_by(donor_id) %>% summarise(reed_epi_cells = n()) 
epi_cells_orig = meta.data %>% filter(original_celltype%in%orig_epi_cells) %>% group_by(donor_id) %>% summarise(orig_epi_cells = n()) 
trm_cells_reed = meta.data %>% filter(level2%in%c("CD8_Trm")) %>% group_by(donor_id) %>% summarise(reed_trm_cells = n())

fibro_cells_orig =  meta.data %>% filter(original_celltype%in%orig_fibro_cells) %>% group_by(donor_id) %>% summarise(orig_fibro_cells = n()) 
macro_cells_reed =  meta.data %>% filter(level2%in%reed_macro_cells) %>% group_by(donor_id) %>% summarise(reed_macro_cells = n()) 
fibroSFRP4_cells =  meta.data %>% filter(original_celltype%in%c('Fibro-SFRP4')) %>% group_by(donor_id) %>% summarise(fibroSFRP4_cells = n()) 
fibroMatrix_cells =  meta.data %>% filter(original_celltype%in%c('Fibro-matrix')) %>% group_by(donor_id) %>% summarise(fibroMatrix_cells = n()) 
M1_cells_orig =  meta.data %>% filter(original_celltype%in%orig_M1) %>% group_by(donor_id) %>% summarise(M1_cells_orig = n()) 
M2_cells_orig =  meta.data %>% filter(original_celltype%in%orig_M2) %>% group_by(donor_id) %>% summarise(M2_cells_orig = n()) 
monocytes_cells_orig =  meta.data %>% filter(original_celltype%in%orig_monocytes) %>% group_by(donor_id) %>% summarise(monocytes_cells_orig = n()) 
all_dc_cells_orig =  meta.data %>% filter(original_celltype%in%orig_all_dc) %>% group_by(donor_id) %>% summarise(all_dc_cells_orig = n()) 
b_cells_cells_orig =  meta.data %>% filter(original_celltype%in%orig_b_cells) %>% group_by(donor_id) %>% summarise(b_cells_cells_orig = n()) 

meta.data$parity = as.vector(meta.data$parity)
meta.data = meta.data %>% filter(parity!='unknown')
meta.data$parous = ifelse(as.numeric(meta.data$parity)>0,'Parous','Nulliparous')
total_cells = meta.data %>% group_by(donor_id) %>% summarise(total_cells = n()) 


demo = total_cells %>%
  left_join(trm_cells_reed) %>%
  left_join(cd45_cells_reed) %>%
  left_join(cd3_cells_reed) %>%
  left_join(cd8_cells_reed) %>%
  left_join(epi_cells_reed) %>%
  left_join(epi_cells_orig) %>%
  left_join(fibro_cells_orig) %>%
  left_join(macro_cells_reed) %>%
  left_join(fibroSFRP4_cells) %>%
  left_join(fibroMatrix_cells) %>%
  left_join(M1_cells_orig) %>%
  left_join(M2_cells_orig) %>%
  left_join(monocytes_cells_orig) %>%
  left_join(all_dc_cells_orig) %>%
  left_join(b_cells_cells_orig) %>%
  left_join(meta.data %>% dplyr::select(donor_id,parous,risk_status, self_reported_ethnicity, dataset) %>% distinct())


## Proportions
## Figures
demo$reed_cd45_epi_ratio = demo$reed_cd45_cells / demo$reed_epi_cells
demo$reed_cd3_epi_ratio = demo$reed_cd3_cells / demo$reed_epi_cells
demo$reed_cd8_epi_ratio = demo$reed_cd8_cells / demo$reed_epi_cells
demo$reed_trm_epi_ratio = demo$reed_trm_cells / demo$reed_epi_cells
## Extended figures
demo$reed_macro_epi_ratio = demo$reed_macro_cells / demo$reed_epi_cells
demo$orig_fibro_epi_ratio = demo$orig_fibro_cells / demo$orig_epi_cells
demo$kumarfibroSFRP4_epi_ratio = demo$fibroSFRP4_cells / demo$orig_epi_cells
demo$kumarfibroMatrix_epi_ratio = demo$fibroMatrix_cells / demo$orig_epi_cells
demo$M1_cells_orig_epi_ratio = demo$M1_cells_orig / demo$orig_epi_cells
demo$M2_cells_orig_epi_ratio = demo$M2_cells_orig / demo$orig_epi_cells
demo$monocytes_cells_orig_epi_ratio = demo$monocytes_cells_orig / demo$orig_epi_cells
demo$all_dc_cells_orig_epi_ratio = demo$all_dc_cells_orig / demo$orig_epi_cells
demo$b_cells_cells_orig_epi_ratio = demo$b_cells_cells_orig / demo$orig_epi_cells



demo = demo%>%filter(risk_status %in% c('AR','HR-Unk','HR-cUnk'))
ratios_to_plot = colnames(demo)[grepl( "ratio" , colnames(demo) )]

# sample size
sample_size = demo %>% group_by(parous) %>% summarize(num=n())
orange = "#FF8000"
blue = "#5757F9"
for (r in ratios_to_plot){
  
  p1 = ggplot(data = demo%>%dplyr::filter(!!as.symbol(r) < 1 & !!as.symbol(r) > 0) %>% left_join(sample_size)%>% mutate(Parity = paste0(parous, "\n", "n=", num)) , aes_string(x = "Parity", y = r)) +   ggtitle(r) +
    scale_fill_manual(values=c( blue,orange)) + theme_minimal() +
    geom_violin(aes(fill = parous))  + ylim(c(0,1)) +
    stat_compare_means(label.x = 1.4,label.y = 0.7, size = 4,aes(label = sprintf("p = %5.4f", as.numeric(..p.format..)))) +
    geom_boxplot(width = 0.05) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 16))
  
  
  print(p1)
  
  
}
