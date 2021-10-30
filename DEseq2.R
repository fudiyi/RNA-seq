#!/usr/bin R
# by fdy 20210825

rm(list = ls())
options(stringsAsFactors = F)

################  1.载入R包  ######### 

library(DESeq2)

################  2.读入 featurecounts 结果文件  #################

setwd("E:/转录组/RNA-Seq/DEseq2_test") #设置文件所放目录
countData_raw <- read.table("all_counts.txt",row.names = 1,header = T)
countData1 = countData_raw[,-c(1,2,3,4)] #去除不需要的列

################  3.去除掉样本总和 count<2 的值  ####################

countData1 <- countData1[rowSums(countData1[,2:ncol(countData1)])>2,] #去除掉样本总和 count<2 的值

###############  4.构建样本信息矩阵  #############
colnames(countData1) #查看列名
condition <- factor(c(rep("WT",3), rep("treat1",3),rep("treat2",3))) #设置样本对应的条件
colData <- data.frame(row.names=colnames(countData1[,2:ncol(countData1)]), condition)
all(rownames(colData) %in% colnames(countData1)) #判断名字是否一致

# 注：如果 countData_raw 文件中列名为非固定形式，使用以下脚本重命名为对应的样本名即可（将23行代码修改成以下两句命令）
# colnames(countData1) <- c("Length","CK-1-0","CK-1-24","CK-1-9","CK-2-0","CK-2-24","CK-2-9","CK-3-0","CK-3-24","CK-3-9")
# condition <- factor(c("CK-0","CK-24","CK-9","CK-0","CK-24","CK-9","CK-0","CK-24","CK-9"))

###############  5.获取表达矩阵  ###############
countData2 <- countData1[, rownames(colData)] #表达矩阵不需要Length此列
all(rownames(colData) == colnames(countData2))

###############  6.构建dds矩阵  ####################

dds <- DESeqDataSetFromMatrix(countData2, colData, design= ~ condition) 

############### 7.PCA分析 ####################

dds_pca <- estimateSizeFactors(dds) #计算每个样本的归一化系数
raw <- SummarizedExperiment(counts(dds_pca, normalized=FALSE),
                            colData=colData(dds_pca))
nor <- SummarizedExperiment(counts(dds_pca, normalized=TRUE),
                            colData=colData(dds_pca))
vsd <- vst(dds_pca)
rld <- rlog(dds_pca)
pdf("PCA_result.pdf")
plotPCA( DESeqTransform(raw), intgroup=c("condition") )
plotPCA( DESeqTransform(nor), intgroup=c("condition") )
plotPCA(vsd, intgroup=c("condition"))
plotPCA(rld, intgroup=c("condition"))
dev.off()


################ 8.利用DESeq函数标准化dds矩阵 ################

dds_DEG <- DESeq(dds)

################   9.获取标准化的 counts 及 TPM  ################

#获取 normalized_counts
normalized_counts <- as.data.frame(counts(dds_DEG, normalized=TRUE))
write.csv(normalized_counts, file="normalized.csv")

#获取 TPM 
kb <- countData1$Length / 1000
counts_data <- countData1[,2:ncol(countData1)]
rpk <- counts_data / kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
write.csv(tpm,file = "tpm_results.csv", quote=F, row.names=T)

#################  10.差异分析  #############

Get_DEG <- function(untreated,treated){
  
  contrastV <- c("condition", treated, untreated)
  res <- results(dds_DEG, contrast=contrastV)
  
  baseA <- counts(dds_DEG, normalized=TRUE)[, colData(dds_DEG)$condition == untreated]    
  if (is.vector(baseA)){
    baseMeanA <- as.data.frame(baseA)
  } else {
    baseMeanA <- as.data.frame(rowMeans(baseA))
  }
  colnames(baseMeanA) <-  untreated
  
  
  baseB <- counts(dds_DEG, normalized=TRUE)[, colData(dds_DEG)$condition == treated]    
  if (is.vector(baseB)){
    baseMeanB <- as.data.frame(baseB)
  } else {
    baseMeanB <- as.data.frame(rowMeans(baseB))
  }
  colnames(baseMeanB) <-  treated
  
  
  res <- cbind(baseMeanA, baseMeanB, as.data.frame(res)) 
  res <- cbind(ID=rownames(res), as.data.frame(res))
  res$baseMean <- rowMeans(cbind(baseA, baseB))
  res$padj[is.na(res$padj)] <- 1
  res$pvalue[is.na(res$pvalue)] <- 1
  
  res_de <- subset(res, res$pvalue < 0.05, select=c('ID',untreated, treated, 'log2FoldChange','pvalue','padj'))

  up_DEG <- subset(res_de, res_de$log2FoldChange >= 1)
  file <- paste("up_",treated,"_vs_",untreated,".csv") 
  write.csv(up_DEG,file = file,row.names=F)
  
  down_DEG <- subset(res_de, res_de$log2FoldChange <= -1)
  file <- paste("down_",treated,"_vs_",untreated,".csv") 
  write.csv(down_DEG,file = file,row.names=F)

}

untreated <- c("WT","WT")
treated <- c("treat1","treat2")

data <- data.frame(untreated,treated)

for(i in 1:nrow(data)){       
  A <- data$untreated[i]
  B <- data$treated[i]
  Get_DEG(untreated = A,treated = B)
  print(paste0(i," analysis finished"))
}

