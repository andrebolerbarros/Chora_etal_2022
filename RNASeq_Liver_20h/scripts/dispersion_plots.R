library("org.Mm.eg.db")
library("xlsx")

setwd("~/Colaco_etal_2019/RNASeq_Liver_20h/")
options(java.parameters = "- Xmx2048m")
source("scripts/DESeq_Model.R",echo=T)
source("scripts/theme_geometry_function.R")

cbPalette <- c("#999999","#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

##########################################################################################################################################################
dir.create(path = "DispersionPlots/", showWarnings = FALSE)

res2<-results(dm,contrast = c("Drug_Inf", "Doxy_Ctrl","PBS_Ctrl"),alpha=0.05)
res3<-results(dm,contrast = c("Drug_Inf", "Doxy_Inf","PBS_Inf"),alpha = 0.05)


res2$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res2),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res2$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res2), column="GENENAME", keytype="ENSEMBL", multiVals="first")

res3$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res3),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res3$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res3), column="GENENAME", keytype="ENSEMBL", multiVals="first")

summary(res2)
summary(res3)

ashr_de2<-lfcShrink(dm, contrast = c("Drug_Inf", "Doxy_Ctrl","PBS_Ctrl"),type="ashr",res = res2)
ashr_de3<-lfcShrink(dm, contrast = c("Drug_Inf", "Doxy_Inf","PBS_Inf"),type="ashr",res = res3)


res2_sub<-res2[,-4]
res3_sub<-res3[,-4]

check2_ashr<-vector()
check3_ashr<-vector()

for (i in 1:ncol(ashr_de2)) {
  check2_ashr[i]<-all(as.data.frame(res2_sub)[,i] == as.data.frame(ashr_de2)[,i],na.rm=T)
  check3_ashr[i]<-all(as.data.frame(res3_sub)[,i] == as.data.frame(ashr_de3)[,i],na.rm=T)
}

print(check2_ashr)
print(check3_ashr)

###############################################################################

Joint<-data.frame(GeneSymbol=res2$symbol,GeneName=res2$geneName,
                  pAdj_NI=ashr_de2$padj, log2FC_NI=ashr_de2$log2FoldChange,
                  pAdj_Inf=ashr_de3$padj,log2FC_Inf=ashr_de3$log2FoldChange)


Joint$Sign<-"NS"
Joint$Sign[Joint$pAdj_NI < 0.05 & Joint$pAdj_Inf < 0.05]<-"Significant for Both"
Joint$Sign[Joint$pAdj_NI < 0.05 & Joint$pAdj_Inf >= 0.05]<-"Only Non-Infected"
Joint$Sign[Joint$pAdj_NI >= 0.05 & Joint$pAdj_Inf < 0.05]<-"Only Infected"
Joint$line<-seq(-20,20,length.out = nrow(Joint))

tiff("DispersionPlots/DP_Doxy_vs_PBS.tiff", 
     width = 1000, height = 1300,res = 300)

ggplot(Joint) +
  theme_geometry(Joint$log2FC_NI, Joint$log2FC_Inf,interval = 5)+
  annotate("text", x = -0.1, y = -0.1, label = "0", size = 3)+
  annotate("text", x = 26, y = 0, label = "Log2FC in Non-Infection", fontface = 'italic',size = 4,angle = 90,)+
  annotate("text", x = 0, y = 26, label = "Log2FC in Infection", size = 4,fontface = 'italic')+
  geom_segment(aes(x = -20, y = 0, xend=24,yend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  geom_segment(aes(x = 0, y = -20, yend=24,xend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  xlim(-28, 28)+ 
  geom_line(aes(x=line,y=line),linetype = "dashed",color="darkgrey")+
  ylim(-28, 28)+
  ggtitle("Dispersion Plot")+
  theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank(),legend.position = "bottom")+
  geom_point(aes(x = log2FC_NI, y = log2FC_Inf,color=Sign), alpha=0.7, size=3)+
  scale_colour_manual(values=cbPalette)

dev.off()
##########################################################################################################################################################
res2<-results(dm,contrast = c("Drug_Inf", "Phen_Ctrl","PBS_Ctrl"),alpha=0.05)
res3<-results(dm,contrast = c("Drug_Inf", "Phen_Inf","PBS_Inf"),alpha = 0.05)


res2$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res2),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res2$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res2), column="GENENAME", keytype="ENSEMBL", multiVals="first")

res3$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res3),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res3$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res3), column="GENENAME", keytype="ENSEMBL", multiVals="first")

summary(res2)
summary(res3)

ashr_de2<-lfcShrink(dm, contrast = c("Drug_Inf", "Phen_Ctrl","PBS_Ctrl"),type="ashr",res = res2)
ashr_de3<-lfcShrink(dm, contrast = c("Drug_Inf", "Phen_Inf","PBS_Inf"),type="ashr",res = res3)


res2_sub<-res2[,-4]
res3_sub<-res3[,-4]

check2_ashr<-vector()
check3_ashr<-vector()

for (i in 1:ncol(ashr_de2)) {
  check2_ashr[i]<-all(as.data.frame(res2_sub)[,i] == as.data.frame(ashr_de2)[,i],na.rm=T)
  check3_ashr[i]<-all(as.data.frame(res3_sub)[,i] == as.data.frame(ashr_de3)[,i],na.rm=T)
}

print(check2_ashr)
print(check3_ashr)

###############################################################################

Joint<-data.frame(GeneSymbol=res2$symbol,GeneName=res2$geneName,
                  pAdj_NI=ashr_de2$padj, log2FC_NI=ashr_de2$log2FoldChange,
                  pAdj_Inf=ashr_de3$padj,log2FC_Inf=ashr_de3$log2FoldChange)


Joint$Sign<-"NS"
Joint$Sign[Joint$pAdj_NI < 0.05 & Joint$pAdj_Inf < 0.05]<-"Significant for Both"
Joint$Sign[Joint$pAdj_NI < 0.05 & Joint$pAdj_Inf >= 0.05]<-"Only Non-Infected"
Joint$Sign[Joint$pAdj_NI >= 0.05 & Joint$pAdj_Inf < 0.05]<-"Only Infected"
Joint$line<-seq(-25,25,length.out = nrow(Joint))

tiff("DispersionPlots/DP_Phen_vs_PBS.tiff", 
     width = 1000, height = 1300,res = 300)

ggplot(Joint) +
  theme_geometry(Joint$log2FC_NI, Joint$log2FC_Inf,interval = 5,
                 xlimit = 25,ylimit = 25)+
  annotate("text", x = -0.1, y = -0.1, label = "0", size = 3)+
  annotate("text", x = 31, y = 0, label = "Log2FC in Non-Infection", fontface = 'italic',size = 4,angle = 90,)+
  annotate("text", x = 0, y = 31, label = "Log2FC in Infection", size = 4,fontface = 'italic')+
  geom_segment(aes(x = -25, y = 0, xend=28,yend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  geom_segment(aes(x = 0, y = -25, yend=28,xend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  xlim(-33, 33)+ 
  geom_line(aes(x=line,y=line),linetype = "dashed",color="darkgrey")+
  ylim(-33, 33)+
  ggtitle("Dispersion Plot")+
  guides(col=guide_legend(nrow=2,byrow=TRUE))+
  theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank(),legend.position = "bottom",legend.text=element_text(size=10))+
  geom_point(aes(x = log2FC_NI, y = log2FC_Inf,color=Sign), alpha=0.7, size=3)+
  scale_colour_manual(values=cbPalette)

dev.off()

##########################################################################################################################################################
res2<-results(dm,contrast = c("Drug_Inf", "Epi_Ctrl","PBS_Ctrl"),alpha=0.05)
res3<-results(dm,contrast = c("Drug_Inf", "Epi_Inf","PBS_Inf"),alpha = 0.05)


res2$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res2),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res2$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res2), column="GENENAME", keytype="ENSEMBL", multiVals="first")

res3$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res3),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res3$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res3), column="GENENAME", keytype="ENSEMBL", multiVals="first")

summary(res2)
summary(res3)

ashr_de2<-lfcShrink(dm, contrast = c("Drug_Inf", "Epi_Ctrl","PBS_Ctrl"),type="ashr",res = res2)
ashr_de3<-lfcShrink(dm, contrast = c("Drug_Inf", "Epi_Inf","PBS_Inf"),type="ashr",res = res3)


res2_sub<-res2[,-4]
res3_sub<-res3[,-4]

check2_ashr<-vector()
check3_ashr<-vector()

for (i in 1:ncol(ashr_de2)) {
  check2_ashr[i]<-all(as.data.frame(res2_sub)[,i] == as.data.frame(ashr_de2)[,i],na.rm=T)
  check3_ashr[i]<-all(as.data.frame(res3_sub)[,i] == as.data.frame(ashr_de3)[,i],na.rm=T)
}

print(check2_ashr)
print(check3_ashr)

###############################################################################

Joint<-data.frame(GeneSymbol=res2$symbol,GeneName=res2$geneName,
                  pAdj_NI=ashr_de2$padj, log2FC_NI=ashr_de2$log2FoldChange,
                  pAdj_Inf=ashr_de3$padj,log2FC_Inf=ashr_de3$log2FoldChange)


Joint$Sign<-"NS"
Joint$Sign[Joint$pAdj_NI < 0.05 & Joint$pAdj_Inf < 0.05]<-"Significant for Both"
Joint$Sign[Joint$pAdj_NI < 0.05 & Joint$pAdj_Inf >= 0.05]<-"Only Non-Infected"
Joint$Sign[Joint$pAdj_NI >= 0.05 & Joint$pAdj_Inf < 0.05]<-"Only Infected"
Joint$line<-seq(-20,20,length.out = nrow(Joint))

tiff("DispersionPlots/DP_Epi_vs_PBS.tiff", 
     width = 1000, height = 1300,res = 300)

ggplot(Joint) +
  theme_geometry(Joint$log2FC_NI, Joint$log2FC_Inf,interval = 5,
                 xlimit = 20, ylimit = 20)+
  annotate("text", x = -0.1, y = -0.1, label = "0", size = 3)+
  annotate("text", x = 24, y = 0, label = "Log2FC in Non-Infection", fontface = 'italic',size = 4,angle = 90,)+
  annotate("text", x = 0, y = 24, label = "Log2FC in Infection", size = 4,fontface = 'italic')+
  geom_segment(aes(x = -20, y = 0, xend=22,yend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  geom_segment(aes(x = 0, y = -20, yend=22,xend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  xlim(-26, 26)+ 
  geom_line(aes(x=line,y=line),linetype = "dashed",color="darkgrey")+
  ylim(-26, 26)+
  ggtitle("Dispersion Plot")+
  guides(col=guide_legend(nrow=2,byrow=TRUE))+
  theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank(),legend.position = "bottom",legend.text=element_text(size=10))+
  geom_point(aes(x = log2FC_NI, y = log2FC_Inf,color=Sign), alpha=0.7, size=3)+
  scale_colour_manual(values=cbPalette)

dev.off()
###################################################################################
save.image(file = paste0("environments/DispersonPlots_",Sys.Date(),".RData",sep=""))
