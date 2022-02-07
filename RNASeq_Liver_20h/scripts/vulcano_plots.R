setwd("~/Colaco_etal_2019/RNASeq_Liver_20h/")
options(java.parameters = "- Xmx2048m")
source("scripts/DESeq_Model.R",echo=T)
#############################################################################
dir.create(path = "Vulcano_Plots/", showWarnings = FALSE)

res1<-results(dm,contrast = c("Drug_Inf", "PBS_Inf","PBS_Ctrl"),alpha = 0.05)
ashr_de1<-lfcShrink(dm, contrast = c("Drug_Inf", "PBS_Inf","PBS_Ctrl"),type="ashr",res = res1)

res1_sub<-res1[,-4]
check1<-vector()

for (j in 1:ncol(ashr_de1)) {
  check1[j]<-all(as.data.frame(res1_sub)[,j] == as.data.frame(ashr_de1)[,j],na.rm=T)
}

print(check1)

ashr_de1$symbol<-mapIds(org.Mm.eg.db,keys=rownames(ashr_de1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
ashr_de1$geneName <- mapIds(org.Mm.eg.db, keys=rownames(ashr_de1), column="GENENAME", keytype="ENSEMBL", multiVals="first")

ashr_de1$padj[is.na(ashr_de1$padj)]<-1
ashr_de1$threshold<-"NoVariance"
ashr_de1$threshold[ashr_de1$padj < 0.05]<-c("-5 < FC < 5")
ashr_de1$threshold[ashr_de1$padj < 0.05 & ashr_de1$log2FoldChange > 5] <- c("FC > 5")
ashr_de1$threshold[ashr_de1$padj < 0.05 & ashr_de1$log2FoldChange < -5] <- c("FC < -5")

ashr_de1<-as.data.frame(ashr_de1)
ashr_de1$threshold<-factor(ashr_de1$threshold)
ashr_de1$threshold<-factor(ashr_de1$threshold, levels(ashr_de1$threshold)
                       [c(4,2,1,3)])

tiff("Vulcano_Plots/VP_PBSInf_vs_PBSNonInf.tiff", 
     width = 1400, height = 1000,res = 300)

ggplot(data=ashr_de1, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlim(c(-25,25))+
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

#############################################################################
res2<-results(dm,contrast = c("Drug_Inf", "Doxy_Ctrl","PBS_Ctrl"),alpha=0.05)
ashr_de2<-lfcShrink(dm, contrast = c("Drug_Inf", "Doxy_Ctrl","PBS_Ctrl"),type="ashr",res = res2)

res2_sub<-res2[,-4]
check2<-vector()

for (j in 1:ncol(ashr_de2)) {
  check2[j]<-all(as.data.frame(res2_sub)[,j] == as.data.frame(ashr_de2)[,j],na.rm=T)
}

print(check2)


ashr_de2$symbol<-mapIds(org.Mm.eg.db,keys=rownames(ashr_de2),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
ashr_de2$geneName <- mapIds(org.Mm.eg.db, keys=rownames(ashr_de2), column="GENENAME", keytype="ENSEMBL", multiVals="first")

ashr_de2$padj[is.na(ashr_de2$padj)]<-1
ashr_de2$threshold<-"NoVariance"
ashr_de2$threshold[ashr_de2$padj < 0.05]<-c("-5 < FC < 5")
ashr_de2$threshold[ashr_de2$padj < 0.05 & ashr_de2$log2FoldChange > 5] <- c("FC > 5")
ashr_de2$threshold[ashr_de2$padj < 0.05 & ashr_de2$log2FoldChange < -5] <- c("FC < -5")

ashr_de2<-as.data.frame(ashr_de2)
ashr_de2$threshold<-factor(ashr_de2$threshold)
ashr_de2$threshold<-factor(ashr_de2$threshold, levels(ashr_de2$threshold)
                       [c(3,1,2)])

tiff("Vulcano_Plots/VP_DoxyNonInf_vs_PBSNonInf.tiff", width = 1400, height = 1000,res = 300)

ggplot(data=ashr_de2, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlim(c(-30, 30)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("darkgrey","lightgrey", "blue"))+
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

#############################################################################
res3<-results(dm,contrast = c("Drug_Inf", "Doxy_Inf","PBS_Inf"),alpha = 0.05)
ashr_de3<-lfcShrink(dm, contrast = c("Drug_Inf", "Doxy_Inf","PBS_Inf"),type="ashr",res = res3)

res3_sub<-res3[,-4]
check3<-vector()

for (j in 1:ncol(ashr_de3)) {
  check3[j]<-all(as.data.frame(res3_sub)[,j] == as.data.frame(ashr_de3)[,j],na.rm=T)
}

print(check3)


ashr_de3$symbol<-mapIds(org.Mm.eg.db,keys=rownames(ashr_de3),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
ashr_de3$geneName <- mapIds(org.Mm.eg.db, keys=rownames(ashr_de3), column="GENENAME", keytype="ENSEMBL", multiVals="first")

ashr_de3$padj[is.na(ashr_de3$padj)]<-1
ashr_de3$threshold<-"NoVariance"
ashr_de3$threshold[ashr_de3$padj < 0.05]<-c("-5 < FC < 5")
ashr_de3$threshold[ashr_de3$padj < 0.05 & ashr_de3$log2FoldChange > 5] <- c("FC > 5")
ashr_de3$threshold[ashr_de3$padj < 0.05 & ashr_de3$log2FoldChange < -5] <- c("FC < -5")

ashr_de3<-as.data.frame(ashr_de3)
ashr_de3$threshold<-factor(ashr_de3$threshold)
ashr_de3$threshold<-factor(ashr_de3$threshold, levels(ashr_de3$threshold)
                           [c(2,1)])

tiff("Vulcano_Plots/VP_DoxyInf_vs_PBSInf.tiff", width = 1400, height = 1000,res = 300)

ggplot(data=ashr_de3, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlim(c(-31, 31)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("darkgrey","lightgrey")) +
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

#############################################################################
res2<-results(dm,contrast = c("Drug_Inf", "Phen_Ctrl","PBS_Ctrl"),alpha=0.05)
ashr_de2<-lfcShrink(dm, contrast = c("Drug_Inf", "Phen_Ctrl","PBS_Ctrl"),type="ashr",res = res2)

res2_sub<-res2[,-4]
check2<-vector()

for (j in 1:ncol(ashr_de2)) {
  check2[j]<-all(as.data.frame(res2_sub)[,j] == as.data.frame(ashr_de2)[,j],na.rm=T)
}

print(check2)


ashr_de2$symbol<-mapIds(org.Mm.eg.db,keys=rownames(ashr_de2),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
ashr_de2$geneName <- mapIds(org.Mm.eg.db, keys=rownames(ashr_de2), column="GENENAME", keytype="ENSEMBL", multiVals="first")

ashr_de2$padj[is.na(ashr_de2$padj)]<-1
ashr_de2$threshold<-"NoVariance"
ashr_de2$threshold[ashr_de2$padj < 0.05]<-c("-5 < FC < 5")
ashr_de2$threshold[ashr_de2$padj < 0.05 & ashr_de2$log2FoldChange > 5] <- c("FC > 5")
ashr_de2$threshold[ashr_de2$padj < 0.05 & ashr_de2$log2FoldChange < -5] <- c("FC < -5")

ashr_de2<-as.data.frame(ashr_de2)
ashr_de2$threshold<-factor(ashr_de2$threshold)
ashr_de2$threshold<-factor(ashr_de2$threshold, 
    levels(ashr_de2$threshold)[c(3,2,1)])

tiff("Vulcano_Plots/VP_PhenNonInf_vs_PBSNonInf.tiff", width = 1400, height = 1000,res = 300)

ggplot(data=ashr_de2, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlim(c(-30, 30)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("darkgrey","red","lightgrey"))+
  geom_vline(xintercept = 5,colour="blue", linetype = "longdash") +
  geom_vline(xintercept = -5,colour="red", linetype = "longdash") +
  geom_hline(yintercept = -log10(0.05),colour="blue", linetype = "longdash") +
  labs(title = "Phen vs PBS (No Infection)") + 
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

#############################################################################
res3<-results(dm,contrast = c("Drug_Inf", "Phen_Inf","PBS_Inf"),alpha = 0.05)
ashr_de3<-lfcShrink(dm, contrast = c("Drug_Inf", "Phen_Inf","PBS_Inf"),type="ashr",res = res3)

res3_sub<-res3[,-4]
check3<-vector()

for (j in 1:ncol(ashr_de3)) {
  check3[j]<-all(as.data.frame(res3_sub)[,j] == as.data.frame(ashr_de3)[,j],na.rm=T)
}

print(check3)


ashr_de3$symbol<-mapIds(org.Mm.eg.db,keys=rownames(ashr_de3),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
ashr_de3$geneName <- mapIds(org.Mm.eg.db, keys=rownames(ashr_de3), column="GENENAME", keytype="ENSEMBL", multiVals="first")

ashr_de3$padj[is.na(ashr_de3$padj)]<-1
ashr_de3$threshold<-"NoVariance"
ashr_de3$threshold[ashr_de3$padj < 0.05]<-c("-5 < FC < 5")
ashr_de3$threshold[ashr_de3$padj < 0.05 & ashr_de3$log2FoldChange > 5] <- c("FC > 5")
ashr_de3$threshold[ashr_de3$padj < 0.05 & ashr_de3$log2FoldChange < -5] <- c("FC < -5")

ashr_de3<-as.data.frame(ashr_de3)
ashr_de3$threshold<-factor(ashr_de3$threshold)
ashr_de3$threshold<-factor(ashr_de3$threshold, levels(ashr_de3$threshold)
                           [c(2,1)])

tiff("Vulcano_Plots/VP_PhenInf_vs_PBSInf.tiff", width = 1400, height = 1000,res = 300)

ggplot(data=ashr_de3, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlim(c(-31, 31)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("darkgrey","lightgrey")) +
  geom_vline(xintercept = 5,colour="blue", linetype = "longdash") +
  geom_vline(xintercept = -5,colour="red", linetype = "longdash") +
  geom_hline(yintercept = -log10(0.05),colour="blue", linetype = "longdash") +
  labs(title = "Phen vs PBS (with Infection)") + 
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

#############################################################################
res2<-results(dm,contrast = c("Drug_Inf", "Epi_Ctrl","PBS_Ctrl"),alpha=0.05)
ashr_de2<-lfcShrink(dm, contrast = c("Drug_Inf", "Epi_Ctrl","PBS_Ctrl"),type="ashr",res = res2)

res2_sub<-res2[,-4]
check2<-vector()

for (j in 1:ncol(ashr_de2)) {
  check2[j]<-all(as.data.frame(res2_sub)[,j] == as.data.frame(ashr_de2)[,j],na.rm=T)
}

print(check2)


ashr_de2$symbol<-mapIds(org.Mm.eg.db,keys=rownames(ashr_de2),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
ashr_de2$geneName <- mapIds(org.Mm.eg.db, keys=rownames(ashr_de2), column="GENENAME", keytype="ENSEMBL", multiVals="first")

ashr_de2$padj[is.na(ashr_de2$padj)]<-1
ashr_de2$threshold<-"NoVariance"
ashr_de2$threshold[ashr_de2$padj < 0.05]<-c("-5 < FC < 5")
ashr_de2$threshold[ashr_de2$padj < 0.05 & ashr_de2$log2FoldChange > 5] <- c("FC > 5")
ashr_de2$threshold[ashr_de2$padj < 0.05 & ashr_de2$log2FoldChange < -5] <- c("FC < -5")

ashr_de2<-as.data.frame(ashr_de2)
ashr_de2$threshold<-factor(ashr_de2$threshold)
ashr_de2$threshold<-factor(ashr_de2$threshold, 
            levels(ashr_de2$threshold)[c(2,1)])

tiff("Vulcano_Plots/VP_EpiNonInf_vs_PBSNonInf.tiff", width = 1400, height = 1000,res = 300)

ggplot(data=ashr_de2, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlim(c(-30, 30)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("darkgrey","lightgrey")) +
  geom_vline(xintercept = 5,colour="blue", linetype = "longdash") +
  geom_vline(xintercept = -5,colour="red", linetype = "longdash") +
  geom_hline(yintercept = -log10(0.05),colour="blue", linetype = "longdash") +
  labs(title = "Epi vs PBS (No Infection)") + 
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

#############################################################################
res3<-results(dm,contrast = c("Drug_Inf", "Epi_Inf","PBS_Inf"),alpha = 0.05)
ashr_de3<-lfcShrink(dm, contrast = c("Drug_Inf", "Epi_Inf","PBS_Inf"),type="ashr",res = res3)

res3_sub<-res3[,-4]
check3<-vector()

for (j in 1:ncol(ashr_de3)) {
  check3[j]<-all(as.data.frame(res3_sub)[,j] == as.data.frame(ashr_de3)[,j],na.rm=T)
}

print(check3)


ashr_de3$symbol<-mapIds(org.Mm.eg.db,keys=rownames(ashr_de3),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
ashr_de3$geneName <- mapIds(org.Mm.eg.db, keys=rownames(ashr_de3), column="GENENAME", keytype="ENSEMBL", multiVals="first")

ashr_de3$padj[is.na(ashr_de3$padj)]<-1
ashr_de3$threshold<-"NoVariance"
ashr_de3$threshold[ashr_de3$padj < 0.05]<-c("-5 < FC < 5")
ashr_de3$threshold[ashr_de3$padj < 0.05 & ashr_de3$log2FoldChange > 5] <- c("FC > 5")
ashr_de3$threshold[ashr_de3$padj < 0.05 & ashr_de3$log2FoldChange < -5] <- c("FC < -5")

ashr_de3<-as.data.frame(ashr_de3)
ashr_de3$threshold<-factor(ashr_de3$threshold)
ashr_de3$threshold<-factor(ashr_de3$threshold, levels(ashr_de3$threshold)
                           [c(2,1)])

tiff("Vulcano_Plots/VP_EpiInf_vs_PBSInf.tiff", width = 1400, height = 1000,res = 300)

ggplot(data=ashr_de3, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlim(c(-31, 31)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("darkgrey","lightgrey")) +
  geom_vline(xintercept = 5,colour="blue", linetype = "longdash") +
  geom_vline(xintercept = -5,colour="red", linetype = "longdash") +
  geom_hline(yintercept = -log10(0.05),colour="blue", linetype = "longdash") +
  labs(title = "Epi vs PBS (with Infection)") + 
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

###############################################################################


save.image(file = paste0("environments/Vulcano_Plots",Sys.Date(),".RData",sep=""))
