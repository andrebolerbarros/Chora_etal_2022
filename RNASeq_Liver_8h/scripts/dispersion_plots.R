#Load the required libraries
library("org.Mm.eg.db")
library("xlsx")

#Determine the working directory
setwd("~/Colaco_etal_2019/RNASeq_Liver_8h/")
options(java.parameters = "- Xmx2048m")

##### Re-running DESeq Model #### 
source("scripts/DESeq_Model.R",echo = T)

#Importing a necessary function for the plot (function adapted from the one in: https://stackoverflow.com/questions/17753101/center-x-and-y-axis-with-ggplot2)
source("scripts/theme_geometry_function.R")

#Determining a colorblind-friendly pallete
cbPalette <- c("#999999","#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

##### Differential Expression Analysis ########
dir.create(path = "DispersionPlots/", showWarnings = FALSE)

#Perform the pairwise comparisons
res1<-results(dm,contrast = c("Group", "PBS_Inf","PBS_Ctrl"),alpha = 0.05)
res2<-results(dm,contrast = c("Group", "Doxy_Ctrl","PBS_Ctrl"),alpha=0.05)
res3<-results(dm,contrast = c("Group", "Doxy_Inf","PBS_Inf"),alpha = 0.05)

#Create new columns for the results to include gene symbols & gene names
res1$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res1$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")

res2$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res2),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res2$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res2), column="GENENAME", keytype="ENSEMBL", multiVals="first")

res3$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res3),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res3$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res3), column="GENENAME", keytype="ENSEMBL", multiVals="first")

#Perform the shrinkage for the log2FC values, using the ashr correction
ashr_de1<-lfcShrink(dm, contrast = c("Group", "PBS_Inf","PBS_Ctrl"),type="ashr",res = res1)
ashr_de2<-lfcShrink(dm, contrast = c("Group", "Doxy_Ctrl","PBS_Ctrl"),type="ashr",res = res2)
ashr_de3<-lfcShrink(dm, contrast = c("Group", "Doxy_Inf","PBS_Inf"),type="ashr",res = res3)

#Perform an internal check to assess if the non-shrinkage results match up with the shrinkage ones. In this case, only the columns corresponding to the log2FC & lfcSE should appear as FALSE

#Perform an internal check to assess if the non-shrinkage results match up with the shrinkage ones. In this case, only the columns corresponding to the log2FC & lfcSE should appear as FALSE
res1_sub<-res1[,-4]
res2_sub<-res2[,-4]
res3_sub<-res3[,-4]

check1_ashr<-vector()
check2_ashr<-vector()
check3_ashr<-vector()


for (i in 1:ncol(ashr_de1)) {
  check1_ashr[i]<-all(as.data.frame(res1_sub)[,i] == as.data.frame(ashr_de1)[,i],na.rm=T)
  check2_ashr[i]<-all(as.data.frame(res2_sub)[,i] == as.data.frame(ashr_de2)[,i],na.rm=T)
  check3_ashr[i]<-all(as.data.frame(res3_sub)[,i] == as.data.frame(ashr_de3)[,i],na.rm=T)
}

print(check1_ashr)
print(check2_ashr)
print(check3_ashr)

########### Creating the plot #######################################

#Create a new data-frame
Joint<-data.frame(GeneSymbol=res2$symbol,GeneName=res2$geneName,
                  pAdj_NI=ashr_de2$padj, log2FC_NI=ashr_de2$log2FoldChange,
                  pAdj_Inf=ashr_de3$padj,log2FC_Inf=ashr_de3$log2FoldChange)

#Create a new variable to determine the overlap of significance between the two groups
Joint$Sign<-"NS"
Joint$Sign[Joint$pAdj_NI < 0.05 & Joint$pAdj_Inf < 0.05]<-"Significant for Both"
Joint$Sign[Joint$pAdj_NI < 0.05 & Joint$pAdj_Inf >= 0.05]<-"Only Non-Infected"
Joint$Sign[Joint$pAdj_NI >= 0.05 & Joint$pAdj_Inf < 0.05]<-"Only Infected"

#Create a new variable for the lines in the plot
Joint$line<-seq(-30,30,length.out = nrow(Joint))

#Create the plot & saving it as tiff
tiff("DispersionPlots/DP_Doxy_vs_PBS.tiff", 
     width = 1000, height = 1300,res = 300)

ggplot(Joint) +
  theme_geometry(Joint$log2FC_NI, Joint$log2FC_Inf,interval = 8,
                 epsilon = 1.5)+
  annotate("text", x = -0.1, y = -0.1, label = "0", size = 3)+
  annotate("text", x = 40, y = 0, label = "Log2FC in Non-Infection", fontface = 'italic',size = 4,angle = 90,)+
  annotate("text", x = 0, y = 40, label = "Log2FC in Infection", size = 4,fontface = 'italic')+
  geom_segment(aes(x = -32, y = 0, xend=36,yend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  geom_segment(aes(x = 0, y = -32, yend=36,xend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  xlim(-42, 42)+ 
  geom_line(aes(x=line,y=line),linetype = "dashed",color="darkgrey")+
  ylim(-42, 42)+
  ggtitle("Dispersion Plot")+
  theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank(),legend.position = "bottom",legend.text=element_text(size=10))+
  geom_point(aes(x = log2FC_NI, y = log2FC_Inf,color=Sign), alpha=0.7, size=3)+
  scale_colour_manual(values=cbPalette)

dev.off()

setEPS()
postscript("DispersionPlots/DP_Doxy_vs_PBS.eps")

ggplot(Joint) +
  theme_geometry(Joint$log2FC_NI, Joint$log2FC_Inf,interval = 8, 
                 epsilon = 1.5)+
  annotate("text", x = -0.1, y = -0.1, label = "0", size = 3)+
  annotate("text", x = 40, y = 0, label = "Log2FC in Non-Infection", fontface = 'italic',size = 4,angle = 90,)+
  annotate("text", x = 0, y = 40, label = "Log2FC in Infection", size = 4,fontface = 'italic')+
  geom_segment(aes(x = -32, y = 0, xend=36,yend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  geom_segment(aes(x = 0, y = -32, yend=36,xend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  xlim(-42, 42)+ 
  geom_line(aes(x=line,y=line),linetype = "dashed",color="darkgrey")+
  ylim(-42, 42)+
  ggtitle("Dispersion Plot")+
  theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank(),legend.position = "bottom",legend.text=element_text(size=10))+
  geom_point(aes(x = log2FC_NI, y = log2FC_Inf,color=Sign), #alpha=0.7,
             size=3)+
  scale_colour_manual(values=cbPalette)

dev.off()


####### Save the R environment
save.image(file = paste0("environments/DispersonPlots_",Sys.Date(),".RData",sep=""))
