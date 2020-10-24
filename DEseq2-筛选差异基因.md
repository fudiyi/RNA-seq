# DEseq2-筛选差异基因

#### 此脚本可一步提取所有差异基因，提供原始表达矩阵 countData 以及 colData 即可

```R
#!/usr/bin R

rm(list = ls()) # 清空环境
options(stringsAsFactors = F)
setwd("your path")
library(DESeq2)

################################ 读取 DEseq2 所需文件 ########################

countData <- as.matrix(read.csv("N7_transcript_count_matrix.csv", row.names=1)) # 读入各个样本的表达矩阵
> countData[1:4,1:4]
            N7_0_1 N7_0_2 N7_0_3 N7_24_1
AT4G04480.1      5      2      8       0
AT1G31380.1      0      0      0      12
AT1G38430.1      0      0      0       0
AT1G03340.1    224    319    235     205

countData <- countData[rowSums(countData)>2,] # 撇掉在多于两个样本中 count<1 的值
colData <- read.table(file= "compare_data.csv", sep=",", row.names=1,header=TRUE) # 需要进行差异分析的表
> colData
        condition
col_0_1      col_0
col_0_2      col_0
col_0_3      col_0
col_24_1    col_24
col_24_2    col_24
col_24_3    col_24
col_3_1      col_3
col_3_2      col_3
col_3_3      col_3

# 若不自己修改原始矩阵可通过如下方式在 R 中修改，若已自定义文件，请用 # 注释此部分

colnames(countData) <- c("N7_0_1","N7_0_2","N7_0_3","N7_24_1","N7_24_2","N7_24_3","N7_3_1","N7_3_2","N7_3_3","col_0_1","col_0_2","col_0_3","col_24_1","col_24_2","col_24_3","col_3_1","col_3_2","col_3_3")  # 修改列名方便结合 condition 设置分组信息
condition <- factor(c(rep("N7_0",3),rep("N7_24",3),rep("N7_3",3),rep("col_0",3),rep("col_24",3),rep("col_3",3)), levels = c("N7_0","N7_24","N7_3","col_0","col_24","col_3")) # 构建列名对应的实验条件，需和 colnames 保持一致
colData <- data.frame(row.names=colnames(countData), condition) 

# 判断countData与colData中的样本名是否一致

all(rownames(colData) %in% colnames(countData)) 
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

# 构建 dds

dds <- DESeqDataSetFromMatrix(countData, colData, design= ~ condition)
dds <- DESeq(dds)
dds
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE)) # 获取标准化后的数据
head(normalized_counts)
write.table(normalized_counts, file="N7_DESeq2_normalized.xls", quote=F, sep="\t", row.names=T, col.names=T)



################################## 自定义差异分析函数 #############################

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
res$padj[is.na(res$pvalue)] <- 1
    
comp314 <- paste(sampleA, "vs", sampleB)

# 生成文件名
file_base <- paste("DESeq2_normalized", comp314, sep="_")
file_base1 <- paste(file_base, "results.xls", sep="_")
write.table(as.data.frame(res), file=file_base1, sep="\t", quote=F, row.names=F)

# 差异基因筛选，pvalue < 0.05
res_de <- subset(res, res$pvalue<0.05, select=c('ID', sampleA, sampleB, 'log2FoldChange','pvalue','padj'))
# foldchange > 2
# 因为 sampleA 是未处理，所以是下调基因
res_de_up <- subset(res_de, res_de$log2FoldChange>=1)
file <- paste("DESeq2_normalized",sampleA,"higherThan",sampleB, 'down.xls', sep="_") 
write.table(as.data.frame(res_de_up), file=file, sep="\t", quote=F, row.names=F)
# 上调基因
res_de_dw <- subset(res_de, res_de$log2FoldChange<=(-1)*1)
file <- paste("DESeq2_normalized",sampleA, "lowerThan", sampleB, 'up.xls', sep="_") 
write.table(as.data.frame(res_de_dw), file=file, sep="\t", quote=F, row.names=F)
}

#################################### 自定义分组 ##############################

sampleA <- c("col_0","col_3","col_24","col_0") # untreat
sampleB <- c("N7_0","N7_3","N7_24","col_3") # treat
data <- data.frame(sampleA,sampleB)
> data
  sampleA sampleB
1   col_0    N7_0
2   col_3    N7_3
3  col_24   N7_24
4   col_0   col_3

# 使用 loop 进行差异分析

for(i in 1:nrow(data)){       # length(data[,1]),nrow(data) 统计行数
	A <- data$sampleA[i]
	B <- data$sampleB[i]
	Get_DEG(sampleA = A,sampleB = B)
}

```

