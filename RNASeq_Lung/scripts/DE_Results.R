#Load the required libraries
library("org.Mm.eg.db")
library("xlsx")

#Determine the working directory
setwd("~/Colaco_etal_2019/RNASeq_Lung/")

##### Re-running DESeq Model #### 
source("scripts/DESeq_Model.R",echo = T)

##### Differential Expression Analysis ########

#Create the folders where we are going to save the results
dir.create(path = "DiffExp_Results/", showWarnings = FALSE)

#Perform the pairwise comparisons
res1<-results(dm,contrast=c("group","PBS_Inf","PBS_NonInf"),alpha = 0.05)
res2<-results(dm,contrast=c("group","Doxy_NonInf","PBS_NonInf"),alpha = 0.05)
res3<-results(dm,contrast=c("group","Doxy_Inf","PBS_Inf"),alpha = 0.05)

#Create new columns for the results to include gene symbols & gene names
res1$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res1$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")

res2$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res2),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res2$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res2), column="GENENAME", keytype="ENSEMBL", multiVals="first")

res3$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res3),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res3$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res3), column="GENENAME", keytype="ENSEMBL", multiVals="first")

#Perform the shrinkage for the log2FC values, using the ashr correction
ashr_de1<-lfcShrink(dm,contrast = c("group","PBS_Inf","PBS_NonInf"),res=res1,type="ashr")
ashr_de2<-lfcShrink(dm,contrast=c("group","Doxy_NonInf","PBS_NonInf"),res=res2,type="ashr")
ashr_de3<-lfcShrink(dm,contrast=c("group","Doxy_Inf","PBS_Inf"),res=res3,type="ashr")

#Perform the shrinkage for the log2FC values, using the normal (beta a-priori distribution) correction
normal_de1<-lfcShrink(dm,contrast = c("group","PBS_Inf","PBS_NonInf"),res=res1,type="normal")
normal_de2<-lfcShrink(dm,contrast=c("group","Doxy_NonInf","PBS_NonInf"),res=res2,type="normal")
normal_de3<-lfcShrink(dm,contrast=c("group","Doxy_Inf","PBS_Inf"),res=res3,type="normal")

#Perform an internal check to assess if the non-shrinkage results match up with the shrinkage ones. In this case, only the columns corresponding to the log2FC & lfcSE should appear as FALSE

check1_ashr<-vector()
check2_ashr<-vector()
check3_ashr<-vector()

check1_normal<-vector()
check2_normal<-vector()
check3_normal<-vector()

#The shrinkage results do not contain the 'stat' column so, for this step, we need to remove that specific column
res1_sub<-res1[,-4]
res2_sub<-res2[,-4]
res3_sub<-res3[,-4]

for (i in 1:ncol(normal_de1)) {
  check1_normal[i]<-all(as.data.frame(res1)[,i] == as.data.frame(normal_de1)[,i],na.rm=T)
  check2_normal[i]<-all(as.data.frame(res2)[,i] == as.data.frame(normal_de2)[,i],na.rm=T)
  check3_normal[i]<-all(as.data.frame(res3)[,i] == as.data.frame(normal_de3)[,i],na.rm=T)
}

for (i in 1:ncol(ashr_de1)) {
  check1_ashr[i]<-all(as.data.frame(res1_sub)[,i] == as.data.frame(ashr_de1)[,i],na.rm=T)
  check2_ashr[i]<-all(as.data.frame(res2_sub)[,i] == as.data.frame(ashr_de2)[,i],na.rm=T)
  check3_ashr[i]<-all(as.data.frame(res3_sub)[,i] == as.data.frame(ashr_de3)[,i],na.rm=T)
  }

print(check1_ashr)
print(check2_ashr)
print(check3_ashr)
print(check1_normal)
print(check2_normal)
print(check3_normal)


############ Data-frame for Res 1 (PBS Inf vs PBS Non Inf) #############

#Create a summary data-frame, with information regarding the gene symbols & names, Adjusted p-value and Log2FC pre-shrinkage and with both corrections (normal & ashr)

sum_res1<-data.frame(res1$symbol,res1$geneName,res1$padj,
              res1$log2FoldChange,ashr_de1$log2FoldChange,
              normal_de1$log2FoldChange,row.names = rownames(res1))

#Rename the columns, to be more user friendly
colnames(sum_res1)<-c("symbol","geneName","Adjusted P-value","log2FC_prefilt","log2FC_ashr","log2FC_normal")

#Remove the entries with adjusted p-value set as NA
sum_res1<-sum_res1[!is.na(sum_res1$`Adjusted P-value`),]

#Remove the entries with adjusted p-value higher than 0.05
sum_res1<-sum_res1[sum_res1$`Adjusted P-value` <0.05,]

#Create new columns to assess the ratio of log2FC corrected vs pre-correction for each correction
sum_res1$ashr_correction<-sum_res1$log2FC_ashr/sum_res1$log2FC_prefilt
sum_res1$normal_correction<-sum_res1$log2FC_normal/sum_res1$log2FC_prefilt


############ Data-frame for Res 2 (Doxy Non Inf vs PBS Non Inf) ###########

#Create a summary data-frame, with information regarding the gene symbols & names, Adjusted p-value and Log2FC pre-shrinkage and with both corrections (normal & ashr)

sum_res2<-data.frame(res2$symbol,res2$geneName,res2$padj,
                res2$log2FoldChange,ashr_de2$log2FoldChange,
                normal_de2$log2FoldChange,row.names = rownames(res2))

#Rename the columns, to be more user friendly
colnames(sum_res2)<-c("symbol","geneName","Adjusted P-value","log2FC_prefilt","log2FC_ashr","log2FC_normal")

#Remove the entries with adjusted p-value set as NA
sum_res2<-sum_res2[!is.na(sum_res2$`Adjusted P-value`),]

#Remove the entries with adjusted p-value higher than 0.05
sum_res2<-sum_res2[sum_res2$`Adjusted P-value` <0.05,]

#Create new columns to assess the ratio of log2FC corrected vs pre-correction for each correction
sum_res2$ashr_correction<-sum_res2$log2FC_ashr/sum_res2$log2FC_prefilt
sum_res2$normal_correction<-sum_res2$log2FC_normal/sum_res2$log2FC_prefilt

############ Data-frame for Res 3 (Doxy Inf vs PBS Inf) ###########

#Create a summary data-frame, with information regarding the gene symbols & names, Adjusted p-value and Log2FC pre-shrinkage and with both corrections (normal & ashr)
sum_res3<-data.frame(res3$symbol,res3$geneName,res3$padj,
                res3$log2FoldChange,ashr_de3$log2FoldChange,
                normal_de3$log2FoldChange,row.names = rownames(res3))

#Rename the columns, to be more user friendly
colnames(sum_res3)<-c("symbol","geneName","Adjusted P-value","log2FC_prefilt","log2FC_ashr","log2FC_normal")

#Remove the entries with adjusted p-value set as NA
sum_res3<-sum_res3[!is.na(sum_res3$`Adjusted P-value`),]

#Remove the entries with adjusted p-value higher than 0.05
sum_res3<-sum_res3[sum_res3$`Adjusted P-value` <0.05,]

#Create new columns to assess the ratio of log2FC corrected vs pre-correction for each correction
sum_res3$ashr_correction<-sum_res3$log2FC_ashr/sum_res3$log2FC_prefilt
sum_res3$normal_correction<-sum_res3$log2FC_normal/sum_res3$log2FC_prefilt


####### Saving the data-frames of the results ##########
#To be easier to go through the data, the data-frames are saved in one Excel file
write.xlsx(sum_res1,file = paste0("DiffExp_Results/DE_",Sys.Date(),".xlsx"),sheetName = "PBS Inf vs NonInf",append = F)
write.xlsx(sum_res2,file = paste0("DiffExp_Results/DE_",Sys.Date(),".xlsx"),sheetName = "Doxy vs PBS NonInf",append = T)
write.xlsx(sum_res3,file = paste0("DiffExp_Results/DE_",Sys.Date(),".xlsx"),sheetName = "Doxy vs PBS Inf",append = T)

####### Save the R environment
save.image(file = paste0("environments/DE_Results_",Sys.Date(),".RData",sep=""))
