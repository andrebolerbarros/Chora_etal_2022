rm(list=ls())

#Load the required libraries
library("org.Mm.eg.db")
library("xlsx")

#Determine the working directory
setwd("~/Colaco_etal_2019/RNASeq_Liver_20h/")

options(java.parameters = "- Xmx50000m")

##### Re-running DESeq Model #### 
source("scripts/DESeq_Model.R",echo=T)

##### Differential Expression Analysis ########

#Create the folders where we are going to save the results
dir.create(path = "DiffExp_Results/", showWarnings = FALSE)

#### PBS Inf vs NI ###
#Perform the pairwise comparisons
res1<-results(dm,contrast = c("Drug_Inf", "PBS_Inf","PBS_Ctrl"),alpha=0.05)

#Create new columns for the results to include gene symbols & gene names
res1$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res1$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")

#Print the summaries for the comparisons
summary(res1)

#Perform the shrinkage for the log2FC values, using the ashr correction
ashr_de1<-lfcShrink(dm,contrast = c("Drug_Inf", "PBS_Inf","PBS_Ctrl"),type="ashr",res = res1)

#Perform the shrinkage for the log2FC values, using the normal (beta a-priori distribution) correction
normal_de1<-lfcShrink(dm,contrast = c("Drug_Inf", "PBS_Inf","PBS_Ctrl"),type="normal",res = res1)

#Perform an internal check to assess if the non-shrinkage results match up with the shrinkage ones. In this case, only the columns corresponding to the log2FC & lfcSE should appear as FALSE

#The shrinkage results do not contain the 'stat' column so, for this step, we need to remove that specific column
res1_sub<-res1[,-4]

check1_ashr<-vector()

check1_normal<-vector()

for (i in 1:ncol(normal_de1)) {
  check1_normal[i]<-all(as.data.frame(res1)[,i] == as.data.frame(normal_de1)[,i],na.rm=T)
  
}

for (i in 1:ncol(ashr_de1)) {
  check1_ashr[i]<-all(as.data.frame(res1_sub)[,i] == as.data.frame(ashr_de1)[,i],na.rm=T)
  }

print(check1_ashr)

print(check1_normal)

### For PBS Inf vs NI ###
#Create a summary data-frame, with information regarding the gene symbols & names, Adjusted p-value and Log2FC pre-shrinkage and with both corrections (normal & ashr)
sum_res1<-data.frame(Gene_Symbol = res1$symbol, Gene_Name=res1$geneName,
                Adjusted_Pvalue=res1$padj,
                log2FC_prefilt=res1$log2FoldChange,
                log2FC_ashr=ashr_de1$log2FoldChange,
                log2FC_normal=normal_de1$log2FoldChange,
                row.names = rownames(res1))

#Remove the entries with adjusted p-value set as NA
sum_res1<-sum_res1[!is.na(sum_res1$Adjusted_Pvalue),]

#Remove the entries with adjusted p-value higher than 0.05
sum_res1<-sum_res1[sum_res1$Adjusted_Pvalue <0.05,]

#Create new columns to assess the ratio of log2FC corrected vs pre-correction for each correction
sum_res1$ashr_correction<-sum_res1$log2FC_ashr/sum_res1$log2FC_prefilt
sum_res1$normal_correction<-sum_res1$log2FC_normal/sum_res1$log2FC_prefilt

write.xlsx(sum_res1,file = paste0("DiffExp_Results/DE_PBS_",Sys.Date(),".xlsx"),sheetName = "PBS Inf vs NI",
  append = F)

#### Doxy vs PBS ###
#Perform the pairwise comparisons
res1<-results(dm,contrast = c("Drug_Inf", "Doxy_Ctrl","PBS_Ctrl"),alpha=0.05)
res2<-results(dm,contrast = c("Drug_Inf", "Doxy_Inf","PBS_Inf"),alpha = 0.05)

#Create new columns for the results to include gene symbols & gene names
res1$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res1$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")

res2$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res2),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res2$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res2), column="GENENAME", keytype="ENSEMBL", multiVals="first")

#Print the summaries for the comparisons
summary(res1)
summary(res2)

#Perform the shrinkage for the log2FC values, using the ashr correction
ashr_de1<-lfcShrink(dm,contrast = c("Drug_Inf", "Doxy_Ctrl","PBS_Ctrl"),type="ashr",res = res1)
ashr_de2<-lfcShrink(dm,contrast = c("Drug_Inf", "Doxy_Inf","PBS_Inf"),type="ashr",res = res2)

#Perform the shrinkage for the log2FC values, using the normal (beta a-priori distribution) correction
normal_de1<-lfcShrink(dm,contrast = c("Drug_Inf", "Doxy_Ctrl","PBS_Ctrl"),type="normal",res = res1)
normal_de2<-lfcShrink(dm,contrast = c("Drug_Inf", "Doxy_Inf","PBS_Inf"),type="normal",res = res2)

#Perform an internal check to assess if the non-shrinkage results match up with the shrinkage ones. In this case, only the columns corresponding to the log2FC & lfcSE should appear as FALSE

#The shrinkage results do not contain the 'stat' column so, for this step, we need to remove that specific column
res1_sub<-res1[,-4]
res2_sub<-res2[,-4]

check1_ashr<-vector()
check2_ashr<-vector()

check1_normal<-vector()
check2_normal<-vector()

for (i in 1:ncol(normal_de1)) {
  check1_normal[i]<-all(as.data.frame(res1)[,i] == as.data.frame(normal_de1)[,i],na.rm=T)
  check2_normal[i]<-all(as.data.frame(res2)[,i] == as.data.frame(normal_de2)[,i],na.rm=T)
  
}

for (i in 1:ncol(ashr_de1)) {
  check1_ashr[i]<-all(as.data.frame(res1_sub)[,i] == as.data.frame(ashr_de1)[,i],na.rm=T)
  check2_ashr[i]<-all(as.data.frame(res2_sub)[,i] == as.data.frame(ashr_de2)[,i],na.rm=T)
  }

print(check1_ashr)
print(check2_ashr)

print(check1_normal)
print(check2_normal)

### For Doxy vs PBS, Non-Infected ###
#Create a summary data-frame, with information regarding the gene symbols & names, Adjusted p-value and Log2FC pre-shrinkage and with both corrections (normal & ashr)
sum_res1<-data.frame(Gene_Symbol = res1$symbol, Gene_Name=res1$geneName,
                Adjusted_Pvalue=res1$padj,
                log2FC_prefilt=res1$log2FoldChange,
                log2FC_ashr=ashr_de1$log2FoldChange,
                log2FC_normal=normal_de1$log2FoldChange,
                row.names = rownames(res1))

#Remove the entries with adjusted p-value set as NA
sum_res1<-sum_res1[!is.na(sum_res1$Adjusted_Pvalue),]

#Remove the entries with adjusted p-value higher than 0.05
sum_res1<-sum_res1[sum_res1$Adjusted_Pvalue <0.05,]

#Create new columns to assess the ratio of log2FC corrected vs pre-correction for each correction
sum_res1$ashr_correction<-sum_res1$log2FC_ashr/sum_res1$log2FC_prefilt
sum_res1$normal_correction<-sum_res1$log2FC_normal/sum_res1$log2FC_prefilt

### For Doxy vs PBS, Non-Infected ###
#Create a summary data-frame, with information regarding the gene symbols & names, Adjusted p-value and Log2FC pre-shrinkage and with both corrections (normal & ashr)
sum_res2<-data.frame(Gene_Symbol = res2$symbol, Gene_Name=res2$geneName,
                Adjusted_Pvalue=res2$padj,
                log2FC_prefilt=res2$log2FoldChange,log2FC_ashr=ashr_de2$log2FoldChange,
                log2FC_normal=normal_de2$log2FoldChange,row.names = rownames(res2))

#Remove the entries with adjusted p-value set as NA
sum_res2<-sum_res2[!is.na(sum_res2$Adjusted_Pvalue),]

#Remove the entries with adjusted p-value higher than 0.05
sum_res2<-sum_res2[sum_res2$Adjusted_Pvalue <0.05,]

#Create new columns to assess the ratio of log2FC corrected vs pre-correction for each correction
sum_res2$ashr_correction<-sum_res2$log2FC_ashr/sum_res2$log2FC_prefilt
sum_res2$normal_correction<-sum_res2$log2FC_normal/sum_res2$log2FC_prefilt

write.xlsx(sum_res1,file = paste0("DiffExp_Results/DE_Doxy_",Sys.Date(),".xlsx"),sheetName = "Doxy_NI vs PBS_NI",
  append = F)
write.xlsx(sum_res2,file = paste0("DiffExp_Results/DE_Doxy_",Sys.Date(),".xlsx"),sheetName = "Doxy_Inf vs PBS_Inf",
  append = T)

#############################################################
#### Phen vs PBS ###
#Perform the pairwise comparisons
res1<-results(dm,contrast = c("Drug_Inf", "Phen_Ctrl","PBS_Ctrl"),alpha=0.05)
res2<-results(dm,contrast = c("Drug_Inf", "Phen_Inf","PBS_Inf"),alpha = 0.05)

#Create new columns for the results to include gene symbols & gene names
res1$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res1$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")

res2$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res2),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res2$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res2), column="GENENAME", keytype="ENSEMBL", multiVals="first")

#Print the summaries for the comparisons
summary(res1)
summary(res2)

#Perform the shrinkage for the log2FC values, using the ashr correction
ashr_de1<-lfcShrink(dm,contrast = c("Drug_Inf", "Phen_Ctrl","PBS_Ctrl"),type="ashr",res = res1)
ashr_de2<-lfcShrink(dm,contrast = c("Drug_Inf", "Phen_Inf","PBS_Inf"),type="ashr",res = res2)

#Perform the shrinkage for the log2FC values, using the normal (beta a-priori distribution) correction
normal_de1<-lfcShrink(dm,contrast = c("Drug_Inf", "Phen_Ctrl","PBS_Ctrl"),type="normal",res = res1)
normal_de2<-lfcShrink(dm,contrast = c("Drug_Inf", "Phen_Inf","PBS_Inf"),type="normal",res = res2)

#Perform an internal check to assess if the non-shrinkage results match up with the shrinkage ones. In this case, only the columns corresponding to the log2FC & lfcSE should appear as FALSE

#The shrinkage results do not contain the 'stat' column so, for this step, we need to remove that specific column
res1_sub<-res1[,-4]
res2_sub<-res2[,-4]

check1_ashr<-vector()
check2_ashr<-vector()

check1_normal<-vector()
check2_normal<-vector()

for (i in 1:ncol(normal_de1)) {
  check1_normal[i]<-all(as.data.frame(res1)[,i] == as.data.frame(normal_de1)[,i],na.rm=T)
  check2_normal[i]<-all(as.data.frame(res2)[,i] == as.data.frame(normal_de2)[,i],na.rm=T)
  
}

for (i in 1:ncol(ashr_de1)) {
  check1_ashr[i]<-all(as.data.frame(res1_sub)[,i] == as.data.frame(ashr_de1)[,i],na.rm=T)
  check2_ashr[i]<-all(as.data.frame(res2_sub)[,i] == as.data.frame(ashr_de2)[,i],na.rm=T)
  }

print(check1_ashr)
print(check2_ashr)

print(check1_normal)
print(check2_normal)

### For Phen vs PBS, Non-Infected ###
#Create a summary data-frame, with information regarding the gene symbols & names, Adjusted p-value and Log2FC pre-shrinkage and with both corrections (normal & ashr)
sum_res1<-data.frame(Gene_Symbol = res1$symbol, Gene_Name=res1$geneName,
                Adjusted_Pvalue=res1$padj,
                log2FC_prefilt=res1$log2FoldChange,
                log2FC_ashr=ashr_de1$log2FoldChange,
                log2FC_normal=normal_de1$log2FoldChange,
                row.names = rownames(res1))

#Remove the entries with adjusted p-value set as NA
sum_res1<-sum_res1[!is.na(sum_res1$Adjusted_Pvalue),]

#Remove the entries with adjusted p-value higher than 0.05
sum_res1<-sum_res1[sum_res1$Adjusted_Pvalue <0.05,]

#Create new columns to assess the ratio of log2FC corrected vs pre-correction for each correction
sum_res1$ashr_correction<-sum_res1$log2FC_ashr/sum_res1$log2FC_prefilt
sum_res1$normal_correction<-sum_res1$log2FC_normal/sum_res1$log2FC_prefilt


### For Phen vs PBS, Non-Infected ###
#Create a summary data-frame, with information regarding the gene symbols & names, Adjusted p-value and Log2FC pre-shrinkage and with both corrections (normal & ashr)
sum_res2<-data.frame(Gene_Symbol = res2$symbol, Gene_Name=res2$geneName,
                Adjusted_Pvalue=res2$padj,
                log2FC_prefilt=res2$log2FoldChange,log2FC_ashr=ashr_de2$log2FoldChange,
                log2FC_normal=normal_de2$log2FoldChange,row.names = rownames(res2))

#Remove the entries with adjusted p-value set as NA
sum_res2<-sum_res2[!is.na(sum_res2$Adjusted_Pvalue),]

#Remove the entries with adjusted p-value higher than 0.05
sum_res2<-sum_res2[sum_res2$Adjusted_Pvalue <0.05,]

#Create new columns to assess the ratio of log2FC corrected vs pre-correction for each correction
sum_res2$ashr_correction<-sum_res2$log2FC_ashr/sum_res2$log2FC_prefilt
sum_res2$normal_correction<-sum_res2$log2FC_normal/sum_res2$log2FC_prefilt

write.xlsx(sum_res1,file = paste0("DiffExp_Results/DE_Phen_",Sys.Date(),".xlsx"),sheetName = "Phen_NI vs PBS_NI",
  append = F)
write.xlsx(sum_res2,file = paste0("DiffExp_Results/DE_Phen_",Sys.Date(),".xlsx"),sheetName = "Phen_Inf vs PBS_Inf",
  append = T)
#############################################################
#### Epi vs PBS ###
#Perform the pairwise comparisons
res1<-results(dm,contrast = c("Drug_Inf", "Epi_Ctrl","PBS_Ctrl"),alpha=0.05)
res2<-results(dm,contrast = c("Drug_Inf", "Epi_Inf","PBS_Inf"),alpha = 0.05)

#Create new columns for the results to include gene symbols & gene names
res1$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res1$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")

res2$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res2),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res2$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res2), column="GENENAME", keytype="ENSEMBL", multiVals="first")

#Print the summaries for the comparisons
summary(res1)
summary(res2)

#Perform the shrinkage for the log2FC values, using the ashr correction
ashr_de1<-lfcShrink(dm,contrast = c("Drug_Inf", "Epi_Ctrl","PBS_Ctrl"),type="ashr",res = res1)
ashr_de2<-lfcShrink(dm,contrast = c("Drug_Inf", "Epi_Inf","PBS_Inf"),type="ashr",res = res2)

#Perform the shrinkage for the log2FC values, using the normal (beta a-priori distribution) correction
normal_de1<-lfcShrink(dm,contrast = c("Drug_Inf", "Epi_Ctrl","PBS_Ctrl"),type="normal",res = res1)
normal_de2<-lfcShrink(dm,contrast = c("Drug_Inf", "Epi_Inf","PBS_Inf"),type="normal",res = res2)

#Perform an internal check to assess if the non-shrinkage results match up with the shrinkage ones. In this case, only the columns corresponding to the log2FC & lfcSE should appear as FALSE

#The shrinkage results do not contain the 'stat' column so, for this step, we need to remove that specific column
res1_sub<-res1[,-4]
res2_sub<-res2[,-4]

check1_ashr<-vector()
check2_ashr<-vector()

check1_normal<-vector()
check2_normal<-vector()

for (i in 1:ncol(normal_de1)) {
  check1_normal[i]<-all(as.data.frame(res1)[,i] == as.data.frame(normal_de1)[,i],na.rm=T)
  check2_normal[i]<-all(as.data.frame(res2)[,i] == as.data.frame(normal_de2)[,i],na.rm=T)
  
}

for (i in 1:ncol(ashr_de1)) {
  check1_ashr[i]<-all(as.data.frame(res1_sub)[,i] == as.data.frame(ashr_de1)[,i],na.rm=T)
  check2_ashr[i]<-all(as.data.frame(res2_sub)[,i] == as.data.frame(ashr_de2)[,i],na.rm=T)
  }

print(check1_ashr)
print(check2_ashr)

print(check1_normal)
print(check2_normal)

### For Epi vs PBS, Non-Infected ###
#Create a summary data-frame, with information regarding the gene symbols & names, Adjusted p-value and Log2FC pre-shrinkage and with both corrections (normal & ashr)
sum_res1<-data.frame(Gene_Symbol = res1$symbol, Gene_Name=res1$geneName,
                Adjusted_Pvalue=res1$padj,
                log2FC_prefilt=res1$log2FoldChange,
                log2FC_ashr=ashr_de1$log2FoldChange,
                log2FC_normal=normal_de1$log2FoldChange,
                row.names = rownames(res1))

#Remove the entries with adjusted p-value set as NA
sum_res1<-sum_res1[!is.na(sum_res1$Adjusted_Pvalue),]

#Remove the entries with adjusted p-value higher than 0.05
sum_res1<-sum_res1[sum_res1$Adjusted_Pvalue <0.05,]

#Create new columns to assess the ratio of log2FC corrected vs pre-correction for each correction
sum_res1$ashr_correction<-sum_res1$log2FC_ashr/sum_res1$log2FC_prefilt
sum_res1$normal_correction<-sum_res1$log2FC_normal/sum_res1$log2FC_prefilt

### For Epi vs PBS, Non-Infected ###
#Create a summary data-frame, with information regarding the gene symbols & names, Adjusted p-value and Log2FC pre-shrinkage and with both corrections (normal & ashr)
sum_res2<-data.frame(Gene_Symbol = res2$symbol, Gene_Name=res2$geneName,
                Adjusted_Pvalue=res2$padj,
                log2FC_prefilt=res2$log2FoldChange,log2FC_ashr=ashr_de2$log2FoldChange,
                log2FC_normal=normal_de2$log2FoldChange,row.names = rownames(res2))

#Remove the entries with adjusted p-value set as NA
sum_res2<-sum_res2[!is.na(sum_res2$Adjusted_Pvalue),]

#Remove the entries with adjusted p-value higher than 0.05
sum_res2<-sum_res2[sum_res2$Adjusted_Pvalue <0.05,]

#Create new columns to assess the ratio of log2FC corrected vs pre-correction for each correction
sum_res2$ashr_correction<-sum_res2$log2FC_ashr/sum_res2$log2FC_prefilt
sum_res2$normal_correction<-sum_res2$log2FC_normal/sum_res2$log2FC_prefilt

write.xlsx(sum_res1,file = paste0("DiffExp_Results/DE_Epi_",Sys.Date(),".xlsx"),sheetName = "Epi_NI vs PBS_NI",
  append = F)
write.xlsx(sum_res2,file = paste0("DiffExp_Results/DE_Epi_",Sys.Date(),".xlsx"),sheetName = "Epi_Inf vs PBS_Inf",
  append = T)

####### Save the R environment
save.image(file = paste0("environments/DE_Results_",Sys.Date(),".RData",sep=""))
