#!/usr/bin R

rm(list = ls()) 
options(stringsAsFactors = F)
setwd("E:/RNA-Seq/YH/hiasat2_stringtie/CPK")
library(DESeq2)
countData_raw <- read.csv("CPK_transcript_count_matrix.csv")
countData_gene <- countData_raw[!grepl("MSTRG", countData_raw$transcript_id),] #删除包含MSTRG的行
write.table(countData_gene, file="CPK_transcript_count_matrix_gene.csv", quote=F, sep=",",row.names=F, col.names=T)
countData <- as.matrix(read.csv("CPK_transcript_count_matrix_gene.csv",row.names=1))
countData <- countData[rowSums(countData)>2,]
colnames(countData) <- c("CPK_0_1","CPK_0_2","CPK_0_3","CPK_24_1","CPK_24_2","CPK_24_3","CPK_3_1","CPK_3_2","CPK_3_3","col_0_1","col_0_2","col_0_3","col_24_1","col_24_2","col_24_3","col_3_1","col_3_2","col_3_3")
condition <- factor(c(rep("CPK_0",3),rep("CPK_24",3),rep("CPK_3",3),rep("col_0",3),rep("col_24",3),rep("col_3",3)), levels = c("CPK_0","CPK_24","CPK_3","col_0","col_24","col_3"))
colData <- data.frame(row.names=colnames(countData), condition) 
all(rownames(colData) %in% colnames(countData)) 
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))
dds <- DESeqDataSetFromMatrix(countData, colData, design= ~ condition)
dds <- DESeq(dds)
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))
head(normalized_counts)
write.table(normalized_counts, file="CPK_DESeq2_normalized.xls", quote=F, sep="\t", row.names=T, col.names=T)

Get_DEG <- function(sampleA,sampleB){
  
  contrastV <- c("condition", sampleA, sampleB)
  res <- results(dds, contrast=contrastV) # 得到sampleA,B分析比较结果
  
  # 在 res 中添加 sampleA 定量信息
  baseA <- counts(dds, normalized=TRUE)[, colData(dds)$condition == sampleA]    
  if (is.vector(baseA)){
    baseMeanA <- as.data.frame(baseA)
  } else {
    baseMeanA <- as.data.frame(rowMeans(baseA))
  }
  colnames(baseMeanA) <- sampleA
  head(baseMeanA)
  
  # 添加 sampleB 定量信息
  baseB <- counts(dds, normalized=TRUE)[, colData(dds)$condition == sampleB]    
  if (is.vector(baseB)){
    baseMeanB <- as.data.frame(baseB)
  } else {
    baseMeanB <- as.data.frame(rowMeans(baseB))
  }
  colnames(baseMeanB) <- sampleB
  head(baseMeanB)
  
  res <- cbind(baseMeanA, baseMeanB, as.data.frame(res)) 
  res <- cbind(ID=rownames(res), as.data.frame(res)) # 添加基因名
  res$baseMean <- rowMeans(cbind(baseA, baseB))
  res$padj[is.na(res$padj)] <- 1 # 把 padj 为 NA 的代替为 1
  res$pvalue[is.na(res$pvalue)] <- 1
  
  comp314 <- paste(sampleA, "vs", sampleB)
  
  # 生成文件名
  file_base <- paste("DESeq2_normalized", comp314, sep="_")
  file_base1 <- paste(file_base, "results.xls", sep="_")
  write.table(as.data.frame(res), file=file_base1, sep="\t", quote=F, row.names=F)
  
  # 差异基因筛选，pvalue < 0.05
  res_de <- subset(res, res$pvalue<0.05, select=c('ID', sampleA, sampleB, 'log2FoldChange', 'pvalue','padj'))
  # foldchange > 2
  # 因为 sampleA 是未处理，所以是下调基因
  res_de_up <- subset(res_de, res_de$log2FoldChange>=1)
  file <- paste("DEG",sampleA,"higherThan",sampleB, 'down.xls', sep="_") 
  write.table(as.data.frame(res_de_up), file=file, sep="\t", quote=F, row.names=F)
  # 上调基因
  res_de_dw <- subset(res_de, res_de$log2FoldChange<=(-1)*1)
  file <- paste("DEG",sampleA, "lowerThan", sampleB, 'up.xls', sep="_") 
  write.table(as.data.frame(res_de_dw), file=file, sep="\t", quote=F, row.names=F)
}

sampleA <- c("col_0","col_0","col_3","col_0","col_3","col_24","CPK_0","CPK_0","CPK_3") 
sampleB <- c("col_3","col_24","col_24","CPK_0","CPK_3","CPK_24","CPK_3","CPK_24","CPK_24") 
data <- data.frame(sampleA,sampleB)


for(i in 1:nrow(data)){       
  A <- data$sampleA[i]
  B <- data$sampleB[i]
  Get_DEG(sampleA = A,sampleB = B)
}

