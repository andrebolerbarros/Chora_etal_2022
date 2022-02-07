library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("org.Mm.eg.db")
library("ggplot2")
library("genefilter")
library("xlsx")
library("stringr")

theme_set(theme_bw())

######################################################################################
setwd("C:/Users/asbarros/Desktop/Bioinfo/RNASeq_Set17/R_161118")

writeLines(capture.output(sessionInfo()),paste0("SessionInfo_",Sys.Date(),".txt",sep=""))

#Import the files
gc_tab<-"gene_count/"
file.names<-paste(gc_tab,dir(gc_tab,pattern="tab"),sep="")

outfile<-list()
for (i in file.names) {
  outfile[[i]]<-read.table(i,row.names=1)
}

#Remove the first four lines of the files, that contain unmappeds, dubious, and other artifacts
for (i in 1:length(outfile)) {
  outfile[[i]]<-outfile[[i]][-c(1:4),]
}

#Select one of the columns to use as the expression data (normally, select the one w/ most reads)
tab<-data.frame(outfile[[1]][,3])
rownames(tab)<-rownames(outfile[[1]])

for (i in 2:length(outfile)) {
  tab[,i]<-outfile[[i]][,3]
  rownames(tab)<-rownames(outfile[[1]])
}

#define the column names
colnames(tab)<-gsub("^.*(Sample_[0-9]+).*","\\1",names(outfile))

#import the design matrix
design<-read.csv2("metadata/sampleTable.csv",header=T,row.names = 1,sep=",")
design$group<-factor(paste0(design$Drug,design$LPS)) #grouping/combinating variables

all(rownames(design)==colnames(tab)) #check if the rownames of design match exactly w/ the column names

#DESeqResults
#########################################################################################

dm <- DESeqDataSetFromMatrix(countData = tab, colData = design, design = ~ group)
dm<-dm[rowSums(counts(dm)) > 0 , ] #remove the genes w/ 0 reads
dm<-estimateSizeFactors(dm) #each sample has a Size Factor, which will help further on in the normalization step

dlog<- rlog(dm,blind=F) #DESeq counts normalization (regular log transgormation)

dm<-estimateDispersions(dm)
plotDispEsts(dm,ylim=c(1e-6, 1e2))

###############################################################################
ddist<-dist(t(assay(dlog)))
heatmap1<-pheatmap(ddist,cluster_rows=T,cluster_cols=T,annotation_col = design[,c("Drug","LPS")])

data1<-plotPCA(dlog, intgroup = c("Drug", "LPS"))+
  ggtitle("Principal Component Analysis based on Overall Genetic Expression across all samples")+
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")

pca1<-data1

design(dm)<- ~ group
dm<-DESeq(dm)

saveRDS(dm, paste0("variables/dm_",Sys.Date(),".rds",sep="")) #remember that you need to use readRDS function to read this file
save.image(file = paste0("environments/DESeq_Anal_",Sys.Date(),".RData",sep=""))
