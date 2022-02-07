library(pheatmap)
library(org.Mm.eg.db)
library(grid)

setwd("C:/Users/asbarros/Desktop/Bioinfo/RNASeq_Set17/R_161118")
setwd("04_01_results/")


genes_interest<-c("Tnf","Il12b","Il6","Il1b","Il10","Cxcl10","Ccl2","Ccl5","Nos2","Ifnb1")

FC_symbol<-mapIds(org.Mm.eg.db,keys=rownames(FC),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")

pheatmap<-pheatmap(FC,color=colors2,breaks=bks1,cluster_cols=F,cluster_rows=F,legend=T,
                   labels_col=c("LPS v Cells","Epi v Ctrl", "Acla v Ctrl", "Epi + LPS v LPS","Acla + LPS v LPS"), 
                   labels_row=FC_symbol,
                   main="Overall Gene Expression",show_rownames = T,show_colnames=T)

heatmap <- pheatmap$gtable

new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 

new.label$label <- ifelse(new.label$label %in% genes_interest, 
                          new.label$label, "")


repel.degree<-0.01

repelled.y <- function(d, d.select, k = repel.degree){
  # d = vector of distances for labels
  # d.select = vector of T/F for which labels are significant
  
  # recursive function to get current label positions
  # (note the unit is "npc" for all components of each distance)
  strip.npc <- function(dd){
    if(!"unit.arithmetic" %in% class(dd)) {
      return(as.numeric(dd))
    }
    
    d1 <- strip.npc(dd$arg1)
    d2 <- strip.npc(dd$arg2)
    fn <- dd$fname
    return(lazyeval::lazy_eval(paste(d1, fn, d2)))
  }
  
  full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
  selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
  
  return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                  to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                  length.out = sum(d.select)),"npc"))
}
new.y.positions <- repelled.y(new.label$y,
                              d.select = new.label$label != "")
new.flag <- segmentsGrob(x0 = new.label$x,
                         x1 = new.label$x + unit(0.15, "npc"),
                         y0 = new.label$y[new.label$label != ""],
                         y1 = new.y.positions)

# shift position for selected labels
new.label$x <- new.label$x + unit(0.2, "npc")
new.label$y[new.label$label != ""] <- new.y.positions

# add flag to heatmap
heatmap <- gtable::gtable_add_grob(x = heatmap,
                                   grobs = new.flag,
                                   t = 4, 
                                   l = 4
)

# replace label positions in heatmap
heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label


tiff("heatmap_leg_genes.tiff", width = 1200, height = 1500,res = 200)
grid.draw(heatmap)
dev.off()
