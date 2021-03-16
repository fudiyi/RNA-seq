#####
rm(list=ls())
setwd('E:/×ªÂ¼×é/RNA-Seq/YH/hisat2_featurecounts_DEseq2/CPK28/DEseq2')
mycounts <- read.table("CPK28_to_TPM.txt",header = T, row.names = 1)

####### TPM ############
kb <- mycounts$Length / 1000
countdata <- mycounts[,2:37]
rpk <- countdata / kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
write.table(tpm,file="CPK28_TPM_result.txt",sep="\t",quote = F)

