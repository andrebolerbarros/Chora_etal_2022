rm(list=ls())

library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("org.Mm.eg.db")
library("ggplot2")
library("genefilter")
library("xlsx")
library("stringr")
library("vsn")

setwd("~/Colaco_etal_2019/RNASeq_Liver_20h/")
options(java.parameters = "- Xmx50000m")
source("scripts/DESeq_Model.R",echo=T)
######################################################################################
dir.create(path = "Doxy_&_Phenp_Overlap_Results/", showWarnings = FALSE)

res1<-results(dm,contrast = c("Drug_Inf", "Doxy_Inf","PBS_Inf"),alpha=0.05)
res2<-results(dm,contrast = c("Drug_Inf", "Phen_Inf","PBS_Inf"),alpha = 0.05)

res1$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res1$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")


res2$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res2),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res2$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res2), column="GENENAME", keytype="ENSEMBL", multiVals="first")

summary(res1)
summary(res2)

ashr_de1<-lfcShrink(dm,contrast = c("Drug_Inf", "Doxy_Ctrl","PBS_Ctrl"),type="ashr",res = res1)
ashr_de2<-lfcShrink(dm,contrast = c("Drug_Inf", "Doxy_Inf","PBS_Inf"),type="ashr",res = res2)

res1_sub<-res1[,-4]
res2_sub<-res2[,-4]

check1<-vector()
check2<-vector()

for (j in 1:ncol(ashr_de1)) {
  check1[j]<-all(as.data.frame(res1_sub)[,j] == as.data.frame(ashr_de1)[,j],na.rm=T)
  check2[j]<-all(as.data.frame(res2_sub)[,j] == as.data.frame(ashr_de2)[,j],na.rm=T)
}

print(check1)
print(check2)

sum<-data.frame(Gene_Symbol = res1$symbol, Gene_Name=res1$geneName,
                Doxy_PAdj=res1$padj, Doxy_Log2FC=ashr_de1$log2FoldChange,
                Phen_PAdj=res2$padj, Phen_Lof2FC=ashr_de2$log2FoldChange)

attach(sum)

############### Up Regulated ################
up_doxy<-sum[!is.na(Doxy_PAdj) & Doxy_PAdj < 0.05 & Doxy_Log2FC > 1,]
dim(up_doxy)

up_phen<-sum[!is.na(Phen_PAdj) & Phen_PAdj < 0.05 & Phen_Lof2FC > 1,]
dim(up_phen)

up_joint<-sum[rownames(sum) %in% rownames(up_doxy) &
            rownames(sum) %in% rownames(up_phen),]
dim(up_joint)

up_only_Doxy<-sum[rownames(sum) %in% rownames(up_doxy) &
                 !(rownames(sum) %in% rownames(up_phen)),]

dim(up_only_Doxy)
up_only_Phen<-sum[!(rownames(sum) %in% rownames(up_doxy)) &
                   rownames(sum) %in% rownames(up_phen),]
dim(up_only_Phen)

write.xlsx(up_only_Doxy,file = paste0("Doxy_&_Phenp_Overlap_Results/Doxy_&_Phenp_Overlap_Up_Regulated_",Sys.Date(),".xlsx"),sheetName = "Only Doxy",append = F)

write.xlsx(up_joint,file = paste0("Doxy_&_Phenp_Overlap_Results/Doxy_&_Phenp_Overlap_Up_Regulated_",Sys.Date(),".xlsx"),sheetName = "Doxy & Phen",append = T)

write.xlsx(up_only_Phen,file = paste0("Doxy_&_Phenp_Overlap_Results/Doxy_&_Phenp_Overlap_Up_Regulated_",Sys.Date(),".xlsx"),sheetName = "Only Phen",append = T)

##################### Down regulated ##########
down_doxy<-sum[!is.na(Doxy_PAdj) & Doxy_PAdj < 0.05 & Doxy_Log2FC < -1,]
dim(down_doxy)

down_phen<-sum[!is.na(Phen_PAdj) & Phen_PAdj < 0.05 & Phen_Lof2FC < -1,]
dim(down_phen)

down_joint<-sum[rownames(sum) %in% rownames(down_doxy) &
                rownames(sum) %in% rownames(down_phen),]
dim(down_joint)

down_only_Doxy<-sum[rownames(sum) %in% rownames(down_doxy) &
                    !(rownames(sum) %in% rownames(down_phen)),]

dim(down_only_Doxy)
down_only_Phen<-sum[!(rownames(sum) %in% rownames(down_doxy)) &
                    rownames(sum) %in% rownames(down_phen),]
dim(down_only_Phen)

write.xlsx(down_only_Doxy,file = paste0("Doxy_&_Phenp_Overlap_Results/Doxy_&_Phenp_Overlap_Down_Regulated_",Sys.Date(),".xlsx"),
           sheetName = "Only Doxy",append = F)

write.xlsx(down_joint,file = paste0("Doxy_&_Phenp_Overlap_Results/Doxy_&_Phenp_Overlap_Down_Regulated_",Sys.Date(),".xlsx"),sheetName = "Doxy & Phen",append = T)

write.xlsx(down_only_Phen,file = paste0("Doxy_&_Phenp_Overlap_Results/Doxy_&_Phenp_Overlap_Down_Regulated_",Sys.Date(),".xlsx"),sheetName = "Only Phen",append = T)
