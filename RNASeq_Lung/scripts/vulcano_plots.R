#Determine the working directory
setwd("~/Colaco_etal_2019/RNASeq_Lung/")

##### Re-running DESeq Model #### 
source("scripts/DESeq_Model.R",echo = T)

#Create the required folder to save the plots
dir.create(path = "Vulcano_Plots/", showWarnings = FALSE)

######### Vulcano Plot PBS Inf vs PBS Non Inf#####################

#Perform the pairwise comparison
res1<-results(dm,contrast=c("group","PBS_Inf","PBS_Ctrl"),alpha=0.05)

#Perform the shrinkage for the log2FC values, using the ashr correction
ashr_de1<-lfcShrink(dm,contrast = c("group","PBS_Inf","PBS_Ctrl"),res=res1,type="ashr")

#Perform an internal check to assess if the non-shrinkage results match up with the shrinkage ones. In this case, only the columns corresponding to the log2FC & lfcSE should appear as FALSE

#The shrinkage results do not contain the 'stat' column so, for this step, we need to remove that specific column
res1_sub<-res1[,-4]

check1<-vector()

for (j in 1:ncol(ashr_de1)) {
  check1[j]<-all(as.data.frame(res1_sub)[,j] == as.data.frame(ashr_de1)[,j],na.rm=T)
}

print(check1)

#Create new columns for the results to include gene symbols & gene names
ashr_de1$symbol<-mapIds(org.Mm.eg.db,keys=rownames(ashr_de1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
ashr_de1$geneName <- mapIds(org.Mm.eg.db, keys=rownames(ashr_de1), column="GENENAME", keytype="ENSEMBL", multiVals="first")

#Replace the NA's with value 1
ashr_de1$padj[is.na(ashr_de1$padj)]<-1

#Determine a new column to be able to color the entries accordingly with their cateogry
ashr_de1$threshold<-"NoVariance"
ashr_de1$threshold[ashr_de1$padj < 0.05]<-c("-5 < FC < 5")
ashr_de1$threshold[ashr_de1$padj < 0.05 & ashr_de1$log2FoldChange > 5] <- c("FC > 5")
ashr_de1$threshold[ashr_de1$padj < 0.05 & ashr_de1$log2FoldChange < -5] <- c("FC < -5")

#Transform the results into a data-frame
ashr_de1<-as.data.frame(ashr_de1)

#Save the new variable as factor & relevel it 
ashr_de1$threshold<-factor(ashr_de1$threshold)
ashr_de1$threshold<-factor(ashr_de1$threshold, 
                           levels(ashr_de1$threshold)[c(4,2,1,3)])

#Produce the plot & save it as tiff
tiff("Vulcano_Plots/VP_PBSInf_vs_PBSNonInf.tiff",
     width = 1400, height = 1000,res = 300)

ggplot(data=ashr_de1, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlim(c(-15,15))+
  #xlim(c(-10, 20)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("darkgrey","red","lightgrey", "blue")) +
  geom_vline(xintercept = 5,colour="blue", linetype = "longdash") +
  geom_vline(xintercept = -5,colour="red", linetype = "longdash") +
  geom_hline(yintercept = -log10(0.05),colour="blue", linetype = "longdash") +
  labs(title = "Infected v Non-Infected (PBS)") + 
  labs(colour="Expression Level") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  annotate("text", x = 12, y =175, label = length(ashr_de1$threshold[ashr_de1$threshold=="FC > 5"]))+
  annotate("text", x = -12, y = 175, label = length(ashr_de1$threshold[ashr_de1$threshold=="FC < -5"]))+ 
  theme(panel.grid.major = element_line(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        plot.title=element_text(hjust=0.5), 
        legend.position = "none", legend.key = element_blank())

dev.off()

#Produce the plot & save it as tiff
setEPS()
postscript("Vulcano_Plots/VP_PBSInf_vs_PBSNonInf.eps")

ggplot(data=ashr_de1, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(#alpha=0.4, 
             size=1.75) +
  xlim(c(-15,15))+
  #xlim(c(-10, 20)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("darkgrey","red","lightgrey", "blue")) +
  geom_vline(xintercept = 5,colour="blue", linetype = "longdash") +
  geom_vline(xintercept = -5,colour="red", linetype = "longdash") +
  geom_hline(yintercept = -log10(0.05),colour="blue", linetype = "longdash") +
  labs(title = "Infected v Non-Infected (PBS)") + 
  labs(colour="Expression Level") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  annotate("text", x = 12, y =175, label = length(ashr_de1$threshold[ashr_de1$threshold=="FC > 5"]))+
  annotate("text", x = -12, y = 175, label = length(ashr_de1$threshold[ashr_de1$threshold=="FC < -5"]))+ 
  theme(panel.grid.major = element_line(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        plot.title=element_text(hjust=0.5), 
        legend.position = "none", legend.key = element_blank())

dev.off()

######### Vulcano Plot Doxy Non Inf vs PBS Non Inf#####################
#Perform the pairwise comparison
res2<-results(dm,contrast=c("group","Doxy_Ctrl","PBS_Ctrl"),alpha=0.05)

#Perform the shrinkage for the log2FC values, using the ashr correction
ashr_de2<-lfcShrink(dm,contrast=c("group","Doxy_Ctrl","PBS_Ctrl"),res=res2,type="ashr")

#Perform an internal check to assess if the non-shrinkage results match up with the shrinkage ones. In this case, only the columns corresponding to the log2FC & lfcSE should appear as FALSE

#The shrinkage results do not contain the 'stat' column so, for this step, we need to remove that specific column
res2_sub<-res2[,-4]

check2<-vector()

for (j in 1:ncol(ashr_de2)) {
  check2[j]<-all(as.data.frame(res2_sub)[,j] == as.data.frame(ashr_de2)[,j],na.rm=T)
}

print(check2)

#Create new columns for the results to include gene symbols & gene names
ashr_de2$symbol<-mapIds(org.Mm.eg.db,keys=rownames(ashr_de2),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
ashr_de2$geneName <- mapIds(org.Mm.eg.db, keys=rownames(ashr_de2), column="GENENAME", keytype="ENSEMBL", multiVals="first")

#Replace the NA's with value 1
ashr_de2$padj[is.na(ashr_de2$padj)]<-1

#Determine a new column to be able to color the entries accordingly with their cateogry
ashr_de2$threshold<-"NoVariance"
ashr_de2$threshold[ashr_de2$padj < 0.05]<-c("-5 < FC < 5")
ashr_de2$threshold[ashr_de2$padj < 0.05 & ashr_de2$log2FoldChange > 5] <- c("FC > 5")
ashr_de2$threshold[ashr_de2$padj < 0.05 & ashr_de2$log2FoldChange < -5] <- c("FC < -5")

#Transform the results into a data-frame
ashr_de2<-as.data.frame(ashr_de2)

#Save the new variable as factor & relevel it 
ashr_de2$threshold<-factor(ashr_de2$threshold)
ashr_de2$threshold<-factor(ashr_de2$threshold, 
                          levels(ashr_de2$threshold)[c(3,2,1)])

#Produce the plot & save it as tiff
tiff("Vulcano_Plots/VP_DoxyNonInf_vs_PBSNonInf.tiff", 
     width = 1400, height = 1000,res = 300)

ggplot(data=ashr_de2, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlim(c(-30, 30)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("darkgrey","red","lightgrey", "blue")) +
  geom_vline(xintercept = 5,colour="blue", linetype = "longdash") +
  geom_vline(xintercept = -5,colour="red", linetype = "longdash") +
  geom_hline(yintercept = -log10(0.05),colour="blue", linetype = "longdash") +
  labs(title = "Doxy vs PBS (No Infection)") + 
  labs(colour="Expression Level") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  annotate("text", x = 10, y =10, label = length(ashr_de2$threshold[ashr_de2$threshold=="FC > 5"]))+
  annotate("text", x = -9, y = 10, label = length(ashr_de2$threshold[ashr_de2$threshold=="FC < -5"]))+ 
  theme(panel.grid.major = element_line(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        plot.title=element_text(hjust=0.5), 
        legend.position = "none", legend.key = element_blank())
dev.off()


setEPS()
postscript("Vulcano_Plots/VP_DoxyNonInf_vs_PBSNonInf.eps")

ggplot(data=ashr_de2, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(#alpha=0.4, 
             size=1.75) +
  xlim(c(-30, 30)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("darkgrey","red","lightgrey", "blue")) +
  geom_vline(xintercept = 5,colour="blue", linetype = "longdash") +
  geom_vline(xintercept = -5,colour="red", linetype = "longdash") +
  geom_hline(yintercept = -log10(0.05),colour="blue", linetype = "longdash") +
  labs(title = "Doxy vs PBS (No Infection)") + 
  labs(colour="Expression Level") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  annotate("text", x = 10, y =10, label = length(ashr_de2$threshold[ashr_de2$threshold=="FC > 5"]))+
  annotate("text", x = -9, y = 10, label = length(ashr_de2$threshold[ashr_de2$threshold=="FC < -5"]))+ 
  theme(panel.grid.major = element_line(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        plot.title=element_text(hjust=0.5), 
        legend.position = "none", legend.key = element_blank())
dev.off()


######### Vulcano Plot Doxy Inf vs PBS Inf#####################

#Perform the pairwise comparison
res3<-results(dm,contrast=c("group","Doxy_Inf","PBS_Inf"),alpha=0.05)

#Perform the shrinkage for the log2FC values, using the ashr correction
ashr_de3<-lfcShrink(dm,contrast=c("group","Doxy_Inf","PBS_Inf"),res=res3,type="ashr")

#Perform an internal check to assess if the non-shrinkage results match up with the shrinkage ones. In this case, only the columns corresponding to the log2FC & lfcSE should appear as FALSE

#The shrinkage results do not contain the 'stat' column so, for this step, we need to remove that specific column
res3_sub<-res3[,-4]
check3<-vector()

for (j in 1:ncol(ashr_de3)) {
  check3[j]<-all(as.data.frame(res3_sub)[,j] == as.data.frame(ashr_de3)[,j],na.rm=T)
}

print(check3)

#Create new columns for the results to include gene symbols & gene names
ashr_de3$symbol<-mapIds(org.Mm.eg.db,keys=rownames(ashr_de3),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
ashr_de3$geneName <- mapIds(org.Mm.eg.db, keys=rownames(ashr_de3), column="GENENAME", keytype="ENSEMBL", multiVals="first")

#Replace the NA's with value 1
ashr_de3$padj[is.na(ashr_de3$padj)]<-1

#Determine a new column to be able to color the entries accordingly with their cateogryashr_de3$threshold<-"NoVariance"
ashr_de3$threshold<-"NoVariance"
ashr_de3$threshold[ashr_de3$padj < 0.05]<-c("-5 < FC < 5")
ashr_de3$threshold[ashr_de3$padj < 0.05 & ashr_de3$log2FoldChange > 5] <- c("FC > 5")
ashr_de3$threshold[ashr_de3$padj < 0.05 & ashr_de3$log2FoldChange < -5] <- c("FC < -5")

#Transform the results into a data-frame
ashr_de3<-as.data.frame(ashr_de3)

#Save the new variable as factor & relevel it 
ashr_de3$threshold<-factor(ashr_de3$threshold)
ashr_de3$threshold<-factor(ashr_de3$threshold, 
                            levels(ashr_de3$threshold)[c(3,2,1)])

#Produce the plot & save it as tiff
tiff("Vulcano_Plots/VP_DoxyInf_vs_PBSInf.tiff",
     width = 1400, height = 1000,res = 300)

ggplot(data=ashr_de3, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlim(c(-30, 30)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("darkgrey","red","lightgrey", "blue")) +
  geom_vline(xintercept = 5,colour="blue", linetype = "longdash") +
  geom_vline(xintercept = -5,colour="red", linetype = "longdash") +
  geom_hline(yintercept = -log10(0.05),colour="blue", linetype = "longdash") +
  labs(title = "Doxy vs PBS (with Infection)") + 
  labs(colour="Expression Level") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  annotate("text", x = 10, y =10, label = length(ashr_de3$threshold[ashr_de3$threshold=="FC > 5"]))+
  annotate("text", x = -9, y = 10, label = length(ashr_de3$threshold[ashr_de3$threshold=="FC < -5"]))+ 
  theme(panel.grid.major = element_line(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        plot.title=element_text(hjust=0.5), 
        legend.position = "none", legend.key = element_blank())
dev.off()

setEPS()
postscript("Vulcano_Plots/VP_DoxyInf_vs_PBSInf.eps")

ggplot(data=ashr_de3, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(#alpha=0.4, 
             size=1.75) +
  xlim(c(-30, 30)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("darkgrey","red","lightgrey", "blue")) +
  geom_vline(xintercept = 5,colour="blue", linetype = "longdash") +
  geom_vline(xintercept = -5,colour="red", linetype = "longdash") +
  geom_hline(yintercept = -log10(0.05),colour="blue", linetype = "longdash") +
  labs(title = "Doxy vs PBS (with Infection)") + 
  labs(colour="Expression Level") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  annotate("text", x = 10, y =10, label = length(ashr_de3$threshold[ashr_de3$threshold=="FC > 5"]))+
  annotate("text", x = -9, y = 10, label = length(ashr_de3$threshold[ashr_de3$threshold=="FC < -5"]))+ 
  theme(panel.grid.major = element_line(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        plot.title=element_text(hjust=0.5), 
        legend.position = "none", legend.key = element_blank())
dev.off()

####### Save the R environment
save.image(file = paste0("environments/Vulcano_Plots_",Sys.Date(),".RData",sep=""))
