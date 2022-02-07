library("org.Mm.eg.db")
library("xlsx")
library("data.table")
library("ggrepel")

setwd("~/Colaco_etal_2019/RNASeq_Liver_20h/")
options(java.parameters = "- Xmx2048m")
source("scripts/DESeq_Model.R",echo=T)
source("scripts/theme_geometry_function_ZoomIn.R")

cbPalette <- c("#999999","#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

##########################################################################################################################################################
dir.create(path = "DispersionPlots/", showWarnings = FALSE)

res2<-results(dm,contrast = c("Drug_Inf", "Doxy_Inf","PBS_Inf"),alpha=0.05)
res3<-results(dm,contrast = c("Drug_Inf", "Phen_Inf","PBS_Inf"),alpha = 0.05)


res2$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res2),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res2$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res2), column="GENENAME", keytype="ENSEMBL", multiVals="first")

res3$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res3),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res3$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res3), column="GENENAME", keytype="ENSEMBL", multiVals="first")

summary(res2)
summary(res3)

ashr_de2<-lfcShrink(dm, contrast = c("Drug_Inf", "Doxy_Inf","PBS_Inf"),type="ashr",res = res2)
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
        pAdj_Doxy=ashr_de2$padj, log2FC_Doxy=ashr_de2$log2FoldChange,
        pAdj_Phen=ashr_de3$padj, log2FC_Phen=ashr_de3$log2FoldChange)


Joint$Sign<-"NS"
Joint$Sign[Joint$pAdj_Doxy < 0.05 & Joint$pAdj_Phen < 0.05]<-"Significant for Both"
Joint$Sign[Joint$pAdj_Doxy < 0.05 & Joint$pAdj_Phen >= 0.05]<-"Only Doxy"
Joint$Sign[Joint$pAdj_Doxy >= 0.05 & Joint$pAdj_Phen < 0.05]<-"Only Phen"

Joint$line<-seq(-24,24,length.out = nrow(Joint))

tiff("DispersionPlots/DP_Doxy_vs_Phen_Inf.tiff", 
     width = 1000, height = 1300,res = 300)

ggplot(Joint) +
  theme_geometry(Joint$log2FC_Doxy, Joint$log2FC_Phen,
                 xlimit = 24,ylimit = 24,epsilon=1.5,interval=12)+
  annotate("text", x = -0.1, y = -0.1, label = "0", size = 3)+
  annotate("text", x = 32, y = 0, label = "Log2FC for Doxy", fontface = 'italic',size = 4,angle = 90,)+
  annotate("text", x = 0, y = 32, label = "Log2FC for Phen", size = 4,fontface = 'italic')+
  geom_segment(aes(x = -24, y = 0, xend=29,yend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  geom_segment(aes(x = 0, y = -24, yend=29,xend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  xlim(-35, 35)+ 
  geom_line(aes(x=line,y=line),linetype = "dashed",color="darkgrey")+
  ylim(-35, 35)+
  ggtitle("Dispersion Plot")+
  guides(col=guide_legend(nrow=2,byrow=TRUE))+
  theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank(),legend.position = "bottom")+
  geom_point(aes(x = log2FC_Doxy, y = log2FC_Phen,color=Sign), alpha=0.7, size=3)+
  geom_point(aes(x = log2FC_Doxy, y = log2FC_Phen), 
             color="red", size=3,
    data=Joint[Joint$GeneSymbol %in% c("Acnat1","Acnat2"),])+
  geom_point(aes(x = log2FC_Doxy, y = log2FC_Phen), 
             color="blue", size=3,
             data=Joint[Joint$GeneSymbol %like% "Adra" |
                        Joint$GeneSymbol %like% "Adrb",])+
  scale_colour_manual(values=cbPalette)

dev.off()


setEPS()
postscript("DispersionPlots/DP_Doxy_vs_Phen_Inf.eps")

ggplot(Joint) +
  theme_geometry(Joint$log2FC_Doxy, Joint$log2FC_Phen,
                 xlimit = 24,ylimit = 24,epsilon=0.5,interval=12)+
  annotate("text", x = -0.1, y = -0.1, label = "0", size = 3)+
  annotate("text", x = 32, y = 0, label = "Log2FC for Doxy", fontface = 'italic',size = 4,angle = 90,)+
  annotate("text", x = 0, y = 32, label = "Log2FC for Phen", size = 4,fontface = 'italic')+
  geom_segment(aes(x = -24, y = 0, xend=29,yend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  geom_segment(aes(x = 0, y = -24, yend=29,xend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  xlim(-35, 35)+ 
  geom_line(aes(x=line,y=line),linetype = "dashed",color="darkgrey")+
  ylim(-35, 35)+
  ggtitle("Dispersion Plot")+
  guides(col=guide_legend(nrow=2,byrow=TRUE))+
  theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank(),legend.position = "bottom")+
  geom_point(aes(x = log2FC_Doxy, y = log2FC_Phen,color=Sign), size=3)+
  geom_point(aes(x = log2FC_Doxy, y = log2FC_Phen), 
             color="red", size=3,
             data=Joint[Joint$GeneSymbol %in% c("Acnat1","Acnat2"),])+
  geom_point(aes(x = log2FC_Doxy, y = log2FC_Phen), 
             color="blue", size=3,
             data=Joint[Joint$GeneSymbol %like% "Adra" |
                          Joint$GeneSymbol %like% "Adrb",])+
  scale_colour_manual(values=cbPalette)

dev.off()

###################################################################################
Joint$line<-seq(-6,6,length.out = nrow(Joint))

tiff("DispersionPlots/DP_Doxy_vs_Phen_Inf_ZoomIn.tiff", 
     width = 1700, height = 2000,res = 600)

ggplot(Joint) +
  theme_geometry(Joint$log2FC_Doxy, Joint$log2FC_Phen,
                 xlimit = 6,ylimit = 6,epsilon=0.3,interval=3)+
  annotate("text", x = -0.1, y = -0.1, label = "0", size = 3)+
  annotate("text", x = 8, y = 0, label = "Log2FC for Doxy", fontface = 'italic',size = 4,angle = 90,)+
  annotate("text", x = 0, y = 8, label = "Log2FC for Phen", size = 4,fontface = 'italic')+
  geom_segment(aes(x = -6, y = 0, xend=7,yend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  geom_segment(aes(x = 0, y = -6, yend=7,xend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  xlim(-9, 9)+ 
  geom_line(aes(x=line,y=line),linetype = "dashed",color="darkgrey")+
  ylim(-9, 9)+
  ggtitle("Dispersion Plot")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.title=element_blank(),legend.position = "bottom")+
  geom_point(aes(x = log2FC_Doxy, y = log2FC_Phen,color=Sign), 
             alpha=0.7, size=0.5)+
  geom_point(aes(x = log2FC_Doxy, y = log2FC_Phen), 
             color="red", size=0.5,
             data=Joint[Joint$GeneSymbol %in% c("Acnat1","Acnat2"),])+
  geom_text(aes(x = log2FC_Doxy, y = log2FC_Phen, label=GeneSymbol), 
             color="red", vjust = 0, nudge_x = -1,size=2,
             data=Joint[Joint$GeneSymbol %in% c("Acnat1","Acnat2"),])+
  
  geom_point(aes(x = log2FC_Doxy, y = log2FC_Phen), 
             color="blue", size=0.5,
             data=Joint[Joint$GeneSymbol %like% "Adra" |
                          Joint$GeneSymbol %like% "Adrb",])+
  geom_text_repel(aes(x = log2FC_Doxy, y = log2FC_Phen,label=GeneSymbol), 
             color="blue", size=2,
             data=Joint[Joint$GeneSymbol %like% "Adra" |
                          Joint$GeneSymbol %like% "Adrb",])+
  guides(col=guide_legend(nrow=2,byrow=TRUE,
                          override.aes = list(size=3)))+
  scale_colour_manual(values=cbPalette)

dev.off()


setEPS()
postscript("DispersionPlots/DP_Doxy_vs_Phen_Inf_ZoomIn.eps")

ggplot(Joint) +
  theme_geometry(Joint$log2FC_Doxy, Joint$log2FC_Phen,
                 xlimit = 6,ylimit = 6,epsilon=0.3,interval=3)+
  annotate("text", x = -0.1, y = -0.1, label = "0", size = 3)+
  annotate("text", x = 8, y = 0, label = "Log2FC for Doxy", fontface = 'italic',size = 4,angle = 90,)+
  annotate("text", x = 0, y = 8, label = "Log2FC for Phen", size = 4,fontface = 'italic')+
  geom_segment(aes(x = -6, y = 0, xend=7,yend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  geom_segment(aes(x = 0, y = -6, yend=7,xend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  xlim(-9, 9)+ 
  geom_line(aes(x=line,y=line),linetype = "dashed",color="darkgrey")+
  ylim(-9, 9)+
  ggtitle("Dispersion Plot")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.title=element_blank(),legend.position = "bottom")+
  geom_point(aes(x = log2FC_Doxy, y = log2FC_Phen,color=Sign), 
             alpha=0.7, size=0.5)+
  geom_point(aes(x = log2FC_Doxy, y = log2FC_Phen), 
             color="red", size=0.5,
             data=Joint[Joint$GeneSymbol %in% c("Acnat1","Acnat2"),])+
  geom_text(aes(x = log2FC_Doxy, y = log2FC_Phen, label=GeneSymbol), 
            color="red", vjust = 0, nudge_x = -1,size=2,
            data=Joint[Joint$GeneSymbol %in% c("Acnat1","Acnat2"),])+
  
  geom_point(aes(x = log2FC_Doxy, y = log2FC_Phen), 
             color="blue", size=0.5,
             data=Joint[Joint$GeneSymbol %like% "Adra" |
                          Joint$GeneSymbol %like% "Adrb",])+
  geom_text_repel(aes(x = log2FC_Doxy, y = log2FC_Phen,label=GeneSymbol), 
                  color="blue", size=2,
                  data=Joint[Joint$GeneSymbol %like% "Adra" |
                               Joint$GeneSymbol %like% "Adrb",])+
  guides(col=guide_legend(nrow=2,byrow=TRUE,
                          override.aes = list(size=3)))+
  scale_colour_manual(values=cbPalette)

dev.off()

tiff("DispersionPlots/DP_Doxy_vs_Phen_Inf_ZoomIn_020720.tiff", 
     width = 1700, height = 2000,res = 600)

ggplot(Joint) +
  theme_geometry(Joint$log2FC_Doxy, Joint$log2FC_Phen,
                 xlimit = 6,ylimit = 6,epsilon=0.3,interval=3)+
  annotate("text", x = -0.1, y = -0.1, label = "0", size = 3)+
  annotate("text", x = 8, y = 0, label = "Log2FC for Doxy", fontface = 'italic',size = 4,angle = 90,)+
  annotate("text", x = 0, y = 8, label = "Log2FC for Phen", size = 4,fontface = 'italic')+
  geom_segment(aes(x = -6, y = 0, xend=7,yend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  geom_segment(aes(x = 0, y = -6, yend=7,xend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  xlim(-9, 9)+ 
  geom_line(aes(x=line,y=line),linetype = "dashed",color="darkgrey")+
  ylim(-9, 9)+
  ggtitle("Dispersion Plot")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.title=element_blank(),legend.position = "bottom")+
  geom_point(aes(x = log2FC_Doxy, y = log2FC_Phen,color=Sign), 
             alpha=0.7, size=0.5)+
  geom_point(aes(x = log2FC_Doxy, y = log2FC_Phen), 
             color="red", size=0.5,
             data=Joint[Joint$GeneSymbol %in% c("Acnat1","Acnat2"),])+
  geom_text(aes(x = log2FC_Doxy, y = log2FC_Phen, label=GeneSymbol), 
            color="red", vjust = 0, nudge_x = -1,size=2,
            data=Joint[Joint$GeneSymbol %in% c("Acnat1","Acnat2"),])+
  
  geom_point(aes(x = log2FC_Doxy, y = log2FC_Phen), 
             color="blue", size=0.5,
             data=Joint[Joint$GeneSymbol %like% "Adra" |
                          Joint$GeneSymbol %like% "Adrb",])+
  geom_text_repel(aes(x = log2FC_Doxy, y = log2FC_Phen,label=GeneSymbol), 
                  color="blue", size=2,
                  data=Joint[Joint$GeneSymbol %like% "Adra" |
                               Joint$GeneSymbol %like% "Adrb",])+
  geom_point(aes(x = log2FC_Doxy, y = log2FC_Phen), 
             color="black", size=0.5,
             data=Joint[Joint$GeneSymbol %in% c("Olfr1033","Olfr1034"),])+
  geom_text(aes(x = log2FC_Doxy, y = log2FC_Phen,label=GeneSymbol), 
                  color="black", size=2,hjust=-.1,
                  data=Joint[Joint$GeneSymbol %in% c("Olfr1033","Olfr1034"),])+
  guides(col=guide_legend(nrow=2,byrow=TRUE,
                          override.aes = list(size=3)))+
  scale_colour_manual(values=cbPalette)

dev.off()

tiff("DispersionPlots/DP_Doxy_vs_Phen_Inf_ZoomIn_050820.tiff", 
     width = 1700, height = 2000,res = 600)

ggplot(Joint) +
  theme_geometry(Joint$log2FC_Doxy, Joint$log2FC_Phen,
                 xlimit = 6,ylimit = 6,epsilon=0.3,interval=3)+
  annotate("text", x = -0.1, y = -0.1, label = "0", size = 3)+
  annotate("text", x = 8, y = 0, label = "Log2FC for Doxy", fontface = 'italic',size = 4,angle = 90,)+
  annotate("text", x = 0, y = 8, label = "Log2FC for Phen", size = 4,fontface = 'italic')+
  geom_segment(aes(x = -6, y = 0, xend=7,yend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  geom_segment(aes(x = 0, y = -6, yend=7,xend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  xlim(-9, 9)+ 
  geom_line(aes(x=line,y=line),linetype = "dashed",color="darkgrey")+
  ylim(-9, 9)+
  ggtitle("Dispersion Plot")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.title=element_blank(),legend.position = "bottom")+
  geom_point(aes(x = log2FC_Doxy, y = log2FC_Phen,color=Sign), 
             alpha=0.7, size=0.5)+
  geom_point(aes(x = log2FC_Doxy, y = log2FC_Phen), 
             color="red", size=0.5,
             data=Joint[Joint$GeneSymbol %in% c("Mavs","Nlrp6","Ccl20"),])+
  geom_text_repel(aes(x = log2FC_Doxy, y = log2FC_Phen, label=GeneSymbol), 
            color="red", vjust = 0, nudge_x = -1,nudge_y = 1.5,size=2,
             data=Joint[Joint$GeneSymbol %in% c("Mavs","Nlrp6"),])+  
  geom_text_repel(aes(x = log2FC_Doxy, y = log2FC_Phen, label=GeneSymbol), 
                  color="red", vjust = 0, nudge_x = .5,nudge_y = -2.2,size=2,
                  data=Joint[Joint$GeneSymbol %in% c("Ccl20"),])+
  geom_point(aes(x = log2FC_Doxy, y = log2FC_Phen), 
             color="blue", size=0.5,
             data=Joint[Joint$GeneSymbol %like% "Adra" |
                          Joint$GeneSymbol %like% "Adrb",])+
  geom_text_repel(aes(x = log2FC_Doxy, y = log2FC_Phen,label=GeneSymbol), 
                  color="blue", size=2,
                  data=Joint[Joint$GeneSymbol %like% "Adra" |
                               Joint$GeneSymbol %like% "Adrb",])+
  geom_point(aes(x = log2FC_Doxy, y = log2FC_Phen), 
             color="black", size=0.5,
             data=Joint[Joint$GeneSymbol %in% c("Olfr1033","Olfr1034"),])+
  geom_text(aes(x = log2FC_Doxy, y = log2FC_Phen,label=GeneSymbol), 
                  color="black", size=2,hjust=-.1,
                  data=Joint[Joint$GeneSymbol %in% c("Olfr1033","Olfr1034"),])+
  guides(col=guide_legend(nrow=2,byrow=TRUE,
                          override.aes = list(size=3)))+
  scale_colour_manual(values=cbPalette)

dev.off()


###################################################################################
save.image(file = paste0("environments/DispersonPlots_Doxy&Phen_",Sys.Date(),".RData",sep=""))
