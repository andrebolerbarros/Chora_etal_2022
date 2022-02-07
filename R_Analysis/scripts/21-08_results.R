rm(list = ls())

setwd("C:/Users/asbarros/Desktop/Bioinfo/RNASeq_Set17/R_161118")
source("DESeq_Anal.R")

setwd("21_08_results/")
writeLines(capture.output(sessionInfo()),paste0("environments/21_08_NoLPS_SessionInfo_",Sys.Date(),".txt",sep=""))

################################### No LPS Conditions #######################################
res1<-results(dm,contrast=c("group","epino_LPS","ctrlno_LPS"),alpha = 0.05)
res2<-results(dm,contrast=c("group","aclano_LPS","ctrlno_LPS"),alpha = 0.05)
res3<-results(dm,contrast=c("group","epino_LPS","aclano_LPS"),alpha = 0.05)

ashr_res1<-lfcShrink(dm,contrast=c("group","epino_LPS","ctrlno_LPS"),type="ashr",res = res1)
ashr_res2<-lfcShrink(dm,contrast=c("group","aclano_LPS","ctrlno_LPS"),type="ashr",res = res2)
ashr_res3<-lfcShrink(dm,contrast=c("group","epino_LPS","aclano_LPS"),type="ashr", res=res3)

res1<-res1[,-4]
res2<-res2[,-4]
res3<-res3[,-4]

check1<-vector()
check2<-vector()
check3<-vector()
for (i in 1:ncol(ashr_res1)) {
  check1[i]<-all(res1[,i] == ashr_res1[,i],na.rm=T)
  check2[i]<-all(res2[,i] == ashr_res2[,i],na.rm=T)
  check3[i]<-all(res3[,i] == ashr_res3[,i],na.rm=T)
}

check1
check2
check3

res1$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res1$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")

res2$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res2$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")

res3$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res3$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")


de1<-data.frame(rownames(res1),res1$symbol,res1$geneName,res1$log2FoldChange,ashr_res1$log2FoldChange,res1$padj,row.names = 1)
colnames(de1)<-c("Gene Symbol","Gene Name", "log2FC_nocorrection","log2FC_ashr","Adjusted p-value")
de1<-de1[!is.na(de1$`Adjusted p-value`),]
de1<-de1[de1$`Adjusted p-value` < 0.05, ]

qplot(de1$log2FC_nocorrection,de1$log2FC_ashr,main="Epi v Ctrl (No LPS)", xlab = "Log2FC No Correction",
      ylab="Log2FC ashr Correction",xlim=c(-10,10),ylim=c(-10,10)) +
  geom_abline(intercept = 0)+
  theme(plot.title = element_text(hjust = 0.5))


de2<-data.frame(rownames(res2),res2$symbol,res2$geneName,res2$log2FoldChange,ashr_res2$log2FoldChange,res2$padj,row.names = 1)
colnames(de2)<-c("Gene Symbol","Gene Name", "log2FC_nocorrection","log2FC_ashr","Adjusted p-value")
de2<-de2[!is.na(de2$`Adjusted p-value`),]
de2<-de2[de2$`Adjusted p-value` < 0.05, ]

qplot(de2$log2FC_nocorrection,de2$log2FC_ashr,main="Acla v Ctrl (No LPS)", xlab = "Log2FC No Correction",
      ylab="Log2FC ashr Correction",xlim=c(-10,10),ylim=c(-10,10)) +
  geom_abline(intercept = 0)+
  theme(plot.title = element_text(hjust = 0.5))


de3<-data.frame(rownames(res3),res3$symbol,res3$geneName,res3$log2FoldChange,ashr_res3$log2FoldChange,res3$padj,row.names = 1)
colnames(de3)<-c("Gene Symbol","Gene Name", "log2FC_nocorrection","log2FC_ashr","Adjusted p-value")
de3<-de3[!is.na(de3$`Adjusted p-value`),]
de3<-de3[de3$`Adjusted p-value` < 0.05, ]

qplot(de3$log2FC_nocorrection,de3$log2FC_ashr,main="Epi v Acla (No LPS)", xlab = "Log2FC No Correction",
      ylab="Log2FC ashr Correction",xlim=c(-10,10),ylim=c(-10,10)) +
  geom_abline(intercept = 0)+
  theme(plot.title = element_text(hjust = 0.5))


write.xlsx(x = de1,file = "Diffs_NoLPS.xls",sheetName = "Epi vs Ctrl", append = F)
write.xlsx(x = de2,file = "Diffs_NoLPS.xls",sheetName = "Acla vs Ctrl", append = T)
write.xlsx(x = de3,file = "Diffs_NoLPS.xls",sheetName = "Epi vs Acla", append = T)


saveRDS(de1, paste0("variables/Epi_v_Ctrl_NoLPs_",Sys.Date(),".rds",sep=""))
saveRDS(de2, paste0("variables/Acla_v_Ctrl_NoLPs_",Sys.Date(),".rds",sep=""))
saveRDS(de3, paste0("variables/Epi_V_Acla_NoLPs_",Sys.Date(),".rds",sep=""))
save.image(file = paste0("environments/21_08_NoLPS_",Sys.Date(),".RData",sep=""))

############################################################################################
############################ LPS Conditions ################################################
rm(list = ls())

setwd("C:/Users/asbarros/Desktop/Bioinfo/RNASeq_Set17/R_161118")
source("DESeq_Anal.R")

setwd("21_08_results/")
writeLines(capture.output(sessionInfo()),paste0("environments/21_08_LPS_SessionInfo_",Sys.Date(),".txt",sep=""))


res1<-results(dm,contrast=c("group","epiLPS","ctrlLPS"),alpha=0.05)
res2<-results(dm,contrast=c("group","aclaLPS","ctrlLPS"),alpha=0.05)
res3<-results(dm,contrast=c("group","epiLPS","aclaLPS"),alpha=0.05)

ashr_res1<-lfcShrink(dm,contrast=c("group","epiLPS","ctrlLPS"),type="ashr",res=res1)
ashr_res2<-lfcShrink(dm,contrast=c("group","aclaLPS","ctrlLPS"),type="ashr",res=res2)
ashr_res3<-lfcShrink(dm,contrast=c("group","epiLPS","aclaLPS"),type="ashr",res=res3)


res1<-res1[,-4]
res2<-res2[,-4]
res3<-res3[,-4]

check1<-vector()
check2<-vector()
check3<-vector()
for (i in 1:ncol(ashr_res1)) {
  check1[i]<-all(res1[,i] == ashr_res1[,i],na.rm=T)
  check2[i]<-all(res2[,i] == ashr_res2[,i],na.rm=T)
  check3[i]<-all(res3[,i] == ashr_res3[,i],na.rm=T)
}

check1
check2
check3

res1$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res1$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")

res2$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res2$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")

res3$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res3$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")


de1<-data.frame(rownames(res1),res1$symbol,res1$geneName,res1$log2FoldChange,ashr_res1$log2FoldChange,res1$padj,row.names = 1)
colnames(de1)<-c("Gene Symbol","Gene Name", "log2FC_nocorrection","log2FC_ashr","Adjusted p-value")
de1<-de1[!is.na(de1$`Adjusted p-value`),]
de1<-de1[de1$`Adjusted p-value` < 0.05, ]

qplot(de1$log2FC_nocorrection,de1$log2FC_ashr,main="Epi v Ctrl (LPS)", xlab = "Log2FC No Correction",
      ylab="Log2FC ashr Correction",xlim=c(-10,10),ylim=c(-10,10)) +
      geom_abline(intercept = 0)+
      theme(plot.title = element_text(hjust = 0.5))
  
de2<-data.frame(rownames(res2),res2$symbol,res2$geneName,res2$log2FoldChange,ashr_res2$log2FoldChange,res2$padj,row.names = 1)
colnames(de2)<-c("Gene Symbol","Gene Name", "log2FC_nocorrection","log2FC_ashr","Adjusted p-value")
de2<-de2[!is.na(de2$`Adjusted p-value`),]
de2<-de2[de2$`Adjusted p-value` < 0.05, ]

qplot(de2$log2FC_nocorrection,de2$log2FC_ashr,main="Acla v Ctrl (LPS)", xlab = "Log2FC No Correction",
      ylab="Log2FC ashr Correction",xlim=c(-10,10),ylim=c(-10,10)) +
  geom_abline(intercept = 0)+
  theme(plot.title = element_text(hjust = 0.5))


de3<-data.frame(rownames(res3),res3$symbol,res3$geneName,res3$log2FoldChange,ashr_res3$log2FoldChange,res3$padj,row.names = 1)
colnames(de3)<-c("Gene Symbol","Gene Name", "log2FC_nocorrection","log2FC_ashr","Adjusted p-value")
de3<-de3[!is.na(de3$`Adjusted p-value`),]
de3<-de3[de3$`Adjusted p-value` < 0.05, ]

qplot(de3$log2FC_nocorrection,de3$log2FC_ashr,main="Epi v Acla (LPS)", xlab = "Log2FC No Correction",
      ylab="Log2FC ashr Correction",xlim=c(-10,10),ylim=c(-10,10)) +
  geom_abline(intercept = 0)+
  theme(plot.title = element_text(hjust = 0.5))

write.xlsx(x = de1,file = "Diffs_LPS.xls",sheetName = "Epi vs Ctrl", append = F)
write.xlsx(x = de2,file = "Diffs_LPS.xls",sheetName = "Acla vs Ctrl", append = T)
write.xlsx(x = de3,file = "Diffs_LPS.xls",sheetName = "Epi vs Acla", append = T)

saveRDS(de1, paste0("variables/Epi_v_Ctrl_LPs_",Sys.Date(),".rds",sep=""))
saveRDS(de2, paste0("variables/Acla_v_Ctrl_LPs_",Sys.Date(),".rds",sep=""))
saveRDS(de3, paste0("variables/Epi_V_Acla_LPs_",Sys.Date(),".rds",sep=""))
save.image(file = paste0("environments/21_08_LPS_",Sys.Date(),".RData",sep=""))