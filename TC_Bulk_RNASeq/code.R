rm(list = ls())
library(fgsea)
library(msigdbr)
library(dplyr)
library(openxlsx)
library(data.table)
library(genefu)
library(knitr)
library(kableExtra)
library(ruv)
library(edgeR)
library(ggplot2)
library(ggrepel)

PrintVolcano <- function(dat, x.col, y.col, name.col, plot.title = "Volcano ", labels = NULL, print.labels = F,
                         x.title = "log Fold Change", y.title = "-log10(p.value)", y.cut = 0.05, x.cut = 1, font.size = .5,
                         high_lbl = 'higher LogFC', low_lbl = 'lower LogFC', show_legend = F) {
  x.cut.txt = format(round(x.cut, 4), nsmall = 2)
  
  options(ggrepel.max.overlaps = Inf)
  blue = "#004C99"
  orange = "#FF8000"
  red = "#990000"
  black = "#000000"
  gray = "#808080"
  white = "#FFFFFF"
  x.cut <- log2(x.cut)
  up_label = paste("Up:" , nrow(dat[which(dat[, y.col] <= y.cut &dat[, x.col] >= x.cut), ]))
  down_label = paste("Down:" , nrow(dat[which(dat[, y.col] <= y.cut &dat[, x.col] < -x.cut),]))
  dat[, "col"] <- "ns"
  dat[which(dat[, y.col] <= y.cut &dat[, x.col] >= x.cut), "col"] <- up_label
  dat[which(dat[, y.col] <= y.cut &dat[, x.col] < -x.cut), "col"] <- down_label
  x_limit = c(-max(abs(dat[,x.col])),max(abs(dat[,x.col])))
  
  p <- ggplot(dat, aes(dat[, x.col], -log10(dat[, y.col]))) +
    ggtitle(label = plot.title) +
    theme(legend.position="top",
          legend.text=element_text(size=20),
          axis.title=element_text(size=18,face="bold"),
          panel.background = element_rect(fill = "white"),
          # panel.grid.major = element_line(size = 0.15, linetype = 'solid',
          #                                 colour = "gray"),
          # panel.grid.minor = element_line(size = 0.15, linetype = 'solid',
          #                                 colour = "gray")) +
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    
    xlim(x_limit)+ #
    #ylim(c(0, max(-log10(dat[, y.col]))))
    ylim(c(-1,max(-log10(dat[, y.col])))+1)
  p <- p + geom_point(aes(dat[, x.col], -log10(dat[, y.col]), color = dat[, "col"]), show.legend = show_legend,size=0.5) +
    # coord_cartesian(xlim=c(-1,1)) +
    labs(x = x.title, y = y.title, color = dat[,"col"]) +#$ ylim(0,max(-log10(dat[, y.col]))+10) +
    scale_color_manual(name = "Filter", labels = c( down_label,
                                                    up_label,
                                                    'ns','ns'),
                       breaks=c(down_label,up_label,'ns','white'),
                       values = c(blue,orange,gray,'white'), limits= c(down_label,up_label),
                       guide= guide_legend(title = "",override.aes = list(size=3)))
  
  
  if (!is.null(labels)) {
    #
    #dat[which(dat[, y.col] >= y.cut | abs(dat[, x.col]) < x.cut), name.col] <- ""
    dat[which(!dat[, name.col] %in% labels), name.col] <- ""
  } else {
    dat[which(dat[, y.col] >= y.cut | abs(dat[, x.col]) < x.cut), name.col] <- ""
  }
  
  if (print.labels) {
    p <- p + geom_label_repel(aes(dat[, x.col], -log10(dat[, y.col]), label = dat[, name.col], fill =  dat[,"col"]),
                              size = font.size, fontface = "bold",color='white',box.padding = 0.5,
                              segment.colour = "grey50"
    )+
      scale_fill_manual(name = "Filter", breaks=c(up_label,down_label,'ns'),
                        values = c(orange,blue,gray),guide = F)
  }
  
  # Higher in label
  # print(
  #   p + annotate("label", x=0,
  #                y=max(-log10(dat[, y.col])+5),
  #                label= paste0(high_lbl), fill = red)
  # )
  return(p)
  
}


NormaliseVoom <- function(raw.data = data, design = NULL, quanseq = F, pop.pctg = 0.60, with.weigths = F, plot = T) {
  if (is.null(raw.data)) {
    return()
  }
  ### Filter low read genes####
  L <- min(colSums(raw.data)) # min library size
  P <- round(ncol(raw.data) * pop.pctg) # % population
  # print(paste('Keep gene if low in less than',P,'samples'))
  keep <- rowSums(cpm(raw.data) > 5 / L * 1e6) >= P
  
  if (quanseq == T) {
    keep <- rowSums(cpm(raw.data) > 5) > P
  }
  dge <- DGEList(raw.data[keep, ])
  ## Normalise
  
  #dge = cpm(dge,method='TMM',log=T)
  # print('Normalising')
  if (with.weigths == T) {
    v <- voomWithQualityWeights(dge, design = design, plot = plot, normalize.method = 'none')
  } else {
    v <- voom(dge, design = design, plot = plot, normalize.method = 'none')
  }
  return(v)
}


PrintDiffExpTable <- function(fit, dir=NULL, contrast.name, gene.annot = NULL, file.sufix = "", print_to_file=F) {
  top.t <- topTable(fit, coef = contrast.name, sort.by = "p", number = Inf, p.value = 1, adjust.method = "fdr")
  top.t$Symbol <- rownames(top.t)
  if (!is.null(gene.annot)) {
    top.t$ENTREZ.ID <- as.integer(rownames(top.t))
    top.t$GENE.SYMBOL <- gene.annot[match(rownames(top.t), gene.annot[, "Entrez.Gene.ID"]), "Approved.Symbol"]
    top.t$CHROMOSOME <- gene.annot[match(rownames(top.t), gene.annot[, "Entrez.Gene.ID"]), "Chrm"]
  }
  if(print_to_file==T){
    write.csv(top.t, paste(dir, contrast.name, file.sufix, ".csv", sep = ""), row.names = F)
  }
  
  return(top.t)
}


demo.file = '../assets/demo.csv'
raw.data.file = '../assets/rawCounts.txt'

## Load data
demo = read.csv(demo.file, stringsAsFactors = F)
raw.data = fread(raw.data.file)
raw.data = data.frame(raw.data,check.names = F)
row.names(raw.data) = raw.data$Symbol
raw.data$Symbol=NULL
###Create dge counts object#####
rownames(demo) = demo$Sample.ID
raw.data = raw.data[,demo$Sample.ID]


## Normalisation (RUVIII, Voom)
demo$Sample.Design = demo$Type

# Remove low count genes
lvls = levels(as.factor(demo$Sample.Design))
min_count <- 5
keep_s = T
keep = F
for(i in lvls){
  samps = filter(demo,i==Sample.Design)%>% pull(Sample.ID)
  for(s in samps){
    keep_s = keep_s & raw.data[,s]>min_count
  }
  keep = keep | keep_s
  
}

raw_data_filt <- raw.data[keep,]
design <- model.matrix(~0+demo$Type+demo$Patient)
colnames(design) <- c(levels(as.factor(as.vector(demo$Type))),'MP',"SS")
rownames(design) <- demo$Sample.ID
#Normalisation for DE
v = NormaliseVoom(raw.data = raw_data_filt,design = design,quanseq = F,with.weigths = T)

## DE Analysis with Limma

g <- factor(demo$Sample.Design)
rep_matrix = replicate.matrix(g)
rownames(rep_matrix) = demo$Sample.ID
fit <- lmFit(object = v,design = rep_matrix)
cm <- makeContrasts(
  
  
  TRMvsAll = TRM - ((Memory_T+Naive_T)/2),
  levels = rep_matrix
)

fit2 = contrasts.fit(fit,cm)
fit2 = eBayes(fit2)

top.t = PrintDiffExpTable(contrast = 'TRMvsAll', fit = fit2, dir = output.dir)
top.t$Symbol = rownames(top.t)

## P_TRM Signature
p_trm_signature = top.t %>% filter(logFC>1,adj.P.Val<0.01) %>% arrange(-logFC)
top.genes = p_trm_signature

## Volcano plot
top.genes$a_logFC = abs(top.genes$logFC)
top10.genes=top.genes %>% top_n(10,wt = (a_logFC))

p = PrintVolcano(dat = top.t,x.col = 'logFC',y.col = 'P.Value', 
                 print.labels = T, labels= top10.genes$Symbol,
                 name.col = 'Symbol',plot.title = 'Volcano FDR<0.05, label top 10',
                 x.title = 'log Fold Change',y.title = 'P.Value',
                 x.cut = 1,y.cut =  0.05, font.size = 2)




