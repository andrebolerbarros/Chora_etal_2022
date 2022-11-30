library(pheatmap)
library(org.Mm.eg.db)
library(grid)

rm(list = ls())

source("DESeq_Anal.R")

setwd("04_01_results/")
writeLines(capture.output(sessionInfo()),paste0("environments/04_01_heatmaps_SessionInfo_",Sys.Date(),".txt",sep=""))

#######################################################################################
res1<-results(dm,contrast=c("group","ctrlLPS","ctrlno_LPS"),alpha=0.05)
res2<-results(dm,contrast=c("group","epino_LPS","ctrlno_LPS"),alpha=0.05)
res3<-results(dm,contrast=c("group","aclano_LPS","ctrlno_LPS"),alpha=0.05)
res4<-results(dm,contrast=c("group","epiLPS","ctrlLPS"),alpha=0.05)
res5<-results(dm,contrast=c("group","aclaLPS","ctrlLPS"),alpha=0.05)
              
ashr_res1<-lfcShrink(dm,contrast=c("group","ctrlLPS","ctrlno_LPS"),type="ashr",res=res1)
ashr_res2<-lfcShrink(dm,contrast=c("group","epino_LPS","ctrlno_LPS"),type="ashr",res=res2)
ashr_res3<-lfcShrink(dm,contrast=c("group","aclano_LPS","ctrlno_LPS"),type="ashr",res=res3)
ashr_res4<-lfcShrink(dm,contrast=c("group","epinLPS","ctrlLPS"),type="ashr",res=res4)
ashr_res5<-lfcShrink(dm,contrast=c("group","aclaLPS","ctrlLPS"),type="ashr",res=res5)

res1<-res1[,-4]
res2<-res2[,-4]
res3<-res3[,-4]
res4<-res4[,-4]
res5<-res5[,-4]

check1<-vector()
check2<-vector()
check3<-vector()
check4<-vector()
check5<-vector()

for (i in 1:ncol(ashr_res1)) {
  check1[i]<-all(res1[,i] == ashr_res1[,i],na.rm=T)
  check2[i]<-all(res2[,i] == ashr_res2[,i],na.rm=T)
  check3[i]<-all(res3[,i] == ashr_res3[,i],na.rm=T)
  check4[i]<-all(res4[,i] == ashr_res4[,i],na.rm=T)
  check5[i]<-all(res5[,i] == ashr_res5[,i],na.rm=T)
  
}

check1
check2
check3
check4
check5

res1$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res1$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")

sum<-data.frame(rownames(res1),res1$symbol,res1$geneName,
                ashr_res1$log2FoldChange,ashr_res1$padj,
                ashr_res2$log2FoldChange,ashr_res2$padj,
                ashr_res3$log2FoldChange,ashr_res3$padj,
                ashr_res4$log2FoldChange,ashr_res4$padj,
                ashr_res5$log2FoldChange,ashr_res5$padj,row.names = 1)

colnames(sum)<-c("Gene Symbol","Gene Name",
                 "log2FC_ashr Ctrl LPS v Ctrl","Adjusted_p-value Ctrl LPS v Ctrl",
                 "log2FC_ashr Epi v Ctrl (No LPS)","Adjusted_p-value Epi v Ctrl (No LPS)",
                 "log2FC_ashr Acla v Ctrl (No LPS)","Adjusted_p-value Acla v Ctrl (No LPS)",                 
                 "log2FC_ashr Epi v Ctrl (LPS)","Adjusted_p-value Epi v Ctrl (LPS)",
                 "log2FC_ashr Acla v Ctrl (LPS)","Adjusted_p-value Acla v Ctrl (LPS)")

saveRDS(sum, paste0("variables/sum_full_",Sys.Date(),".rds",sep=""))

####################################################################

sum<-sum[!is.na(sum$`Adjusted_p-value Ctrl LPS v Ctrl`),]
sum<-sum[sum$`Adjusted_p-value Ctrl LPS v Ctrl` < 0.05,]

saveRDS(sum, paste0("variables/sum_filt_pval_",Sys.Date(),".rds",sep=""))

sum_FC5<-sum[sum$`log2FC_ashr Ctrl LPS v Ctrl` > 5,]


#############################################################
FC<-data.frame(sum_FC5$`log2FC_ashr Ctrl LPS v Ctrl`, 
               sum_FC5$`log2FC_ashr Epi v Ctrl (No LPS)`,
               sum_FC5$`log2FC_ashr Acla v Ctrl (No LPS)`,
               sum_FC5$`log2FC_ashr Epi v Ctrl (LPS)`,
               sum_FC5$`log2FC_ashr Acla v Ctrl (LPS)`)

rownames(FC)<-rownames(sum_FC5)

FC$symbol<-mapIds(org.Mm.eg.db,keys=rownames(FC),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
FC$geneName <- mapIds(org.Mm.eg.db, keys=rownames(FC), column="GENENAME", keytype="ENSEMBL", multiVals="first")


#FC<-FC[order(FC$sum_FC5..log2FC_ashr.Ctrl.LPS.v.Ctrl.,decreasing = T),]

bks1<-c(seq(-10,-0.0005,length=25),0,seq(0.0005,12,length=25))
col1<-colorRampPalette(c("red","bisque"))(25)
col2<-c("lightgrey")
col3<-colorRampPalette(c("aliceblue","blue"))(25)
colors2<-c(col1,col2,col3)

png("heatmap_leg.png", width = 550, height = 700)

pheatmap(FC,color=colors2,breaks=bks1,cluster_cols=F,cluster_rows=F,legend=T,labels_col=c("LPS v Cells","Epi v Ctrl", "Acla v Ctrl", "Epi + LPS v LPS","Acla + LPS v LPS"), main="Overall Gene Expression",show_rownames = F,show_colnames=T)

dev.off()

png("heatmap_noleg.png", width = 550, height = 700)

pheatmap(FC,color=colors2,breaks=bks1,cluster_cols=F,cluster_rows=F,legend=T,labels_col=c("LPS v Cells","Epi v Ctrl", "Acla v Ctrl", "Epi + LPS v LPS","Acla + LPS v LPS"), main="Overall Gene Expression",show_rownames = F,show_colnames=F)

dev.off()

save.image(file = paste0("environments/04_01_heatmaps_",Sys.Date(),".RData",sep=""))
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
