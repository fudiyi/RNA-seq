

# RNA-seq	procedure	for	Yanglab

## 如果没有相关RNA-seq分析基础，请先参考以下两篇经典文章

1. tophat + cufflink: https://www.ncbi.nlm.nih.gov/pubmed/22383036
2. hisat2 + stringtie: https://www.ncbi.nlm.nih.gov/pubmed/?term=Transcript-level+expression+analysis+of+RNA-seq+experiments+with+HISAT%2C+StringTie+and+Ballgown

**快速了解 RNA-seq是什么！**

**1. RNA-seq 小故事**

https://www.jianshu.com/p/d09e624efcab?utm_campaign=hugo&utm_medium=reader_share&utm_content=note&utm_source=weixin-friends

**2. 测序原理**

https://zhuanlan.zhihu.com/p/20702684

注：所有本文用到的软件均在官网有详细说明

### 根据自己需求：从以下方法选择一种组合即可

**A. 普通 RNA-Seq 分析**：步骤 1(质控) + 2(比对) + 7(定量) + 8(差异分析)

**B. 可变剪接分析**：步骤 1(质控) + 2(比对) + 10(可变剪接)

**C. 预测新的转录本**：步骤 1(质控) + 2(比对) + 5(拼接转录本) + 6(合并转录本) + 7(定量) + 8(差异分析) + 9(预测转录本)


### 数据样本：M,S,Col 此数据基本包括了所有转录本分析所需内容

### 数据类型：paired-end, 150bp, 10×， fr-firststrand(链特异性建库)

**什么是单端测序和双端测序？**

https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/paired-end-vs-single-read.html


**什么是链特异性建库？**

https://www.jianshu.com/p/a63595a41bed

此套数据目的之一是为了**预测拟南芥基因组上lncRNA**，lncRNA大多数处于基因的反义链上，所以在建库的时候使用了**ssRNA-Seq**，若无此需求使用普通建库即可

**什么是lncRNA**：https://en.wikipedia.org/wiki/Long_non-coding_RNA

`rawdata: Col-1-0_368368_all.R1.fastq.gz;Col-1-0_368368_all.R2.fastq.gz 
 #双端测序是在两端设计引物进行测序，因此有R1,R2两个fq文件`

**什么是fq文件**：https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html

## 1. 质控
在进行数据分析之前需要对下机数据进行质检，目的是为了判断数据是否达标，大部分公司返回的测序数据为**Cleandata（已去接头）**，质量均不错

**主要判断依据**：Perbase sequence quality  **Q20过滤法**：箱型图10%的线大于 Q=20  Q=-10*lg10(error P)

注：Per base sequence content 前几个碱基测序时候因为状态调整会导致测序略有偏差

**质控结果怎么看**：https://zhuanlan.zhihu.com/p/20731723

### fastqc（质控软件）

**reference：**

**fastqc**： http://darlinglab.org/tutorials/fastqc/
**mulitiqc**： https://multiqc.info/

```shell
fastqc /data/FDY_analysis/RNA_seq/FDY/mac3ab/rawdata/*.fastq.gz -o fastqc # *表示通配符，可对目录下所有.gz文件批量质控
或 ls *.gz | while read id ;do fastqc $id ;done # while read 一次读取一行
multiqc *fastqc.zip --ignore *.html # 整合质控结果
```


## 2. 序列比对

### hisat2 or STAR

**比对工具的选择**：https://www.jianshu.com/p/681e02e7f9af

**hisat2 怎么用**： https://daehwankimlab.github.io/hisat2/manual/

**STAR**： https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html

**添加环境变量：（在首次安装软件之前需配置环境，包括fastqc）**

**什么是环境变量：**

https://www.jianshu.com/p/9c2bf27c3921

http://blog.sciencenet.cn/home.php?mod=space&uid=118204&do=blog&id=1226040

```shell
vi ~/.bashrc #永久修改环境变量，可直接调用 hisat2
export PATH=/data/sly/tools/hisat2-2.0.4/hisat2:$PATH # 在 bashrc 中加入此命令
source ~/.bashrc #使修改生效
```
### 2.1 hisat2

**2.1.1 构建索引：**

**why index**：高通量测序遇到的第一个问题就是，成千上万甚至上几亿条read如果在合理的时间内比对到参考基因组上，并且保证错误率在接受范围内。为了提高比对速度，就需要根据参考基因组序列，经过BWT算法转换成index，而我们比对的序列其实是index的一个子集。当然转录组比对还要考虑到可变剪切的情况，所以更加复杂。
因此我们不是直接把read回贴到基因组上，而是把read和index进行比较

**reference**： http://www.biotrainee.com/thread-26-1-1.html

```shell
hisat2-build -p 10 Zea_mays.AGPv4.dna.toplevel.fa genome
```
**什么是fa文件**：https://zh.wikipedia.org/wiki/FASTA%E6%A0%BC%E5%BC%8F

**2.2.2 比对：**

初次进行比对时先尝试使用一组样本，再尝试批量比对

```shell
wkpath1=/data/FDY_analysis/RNA_seq/FDY/mac3ab/rawdata/cleandata #设置工作路径
wkpath2=/data/FDY_analysis/RNA_seq/FDY/mac3ab/rawdata/hisat2_results_for_cufflink
for i in $(ls ${wkpath1}/*.R1.fastq.gz) # $(your command) 用法等于 `your command`
do
    sample_name=`basename $i|sed s/.R1.*//g` #取文件名（sed 用于将 /.R1.*/ 前后的字符替换为空格）
    或 sample_name=`basename $i .R1.fastq.gz`
    hisat2 -p 12 \ #线程
    --dta-cufflinks \ #设置此参数用于后续分析结合cfflinks,若在拼接转录本时使用stringtie，则只加 --dta
    --rna-strandness RF \ #链特异性，若非此建库方式，则不使用此命令！！！
    -x /data/FDY_analysis/Arabidposis_index_hisat2/genome \ #hisat2-build 构建的索引文件
    -1 ${wkpath1}/${sample_name}.R1.fastq.gz \
    -2 ${wkpath1}/${sample_name}.R2.fastq.gz \
    -S ${wkpath2}/results_bam/${sample_name}.hisat2.sam #输出为sam文件格式
    samtools sort -@ 12 -o ${wkpath2}/results_bam/${sample_name}.hisat2.bam ${wkpath2}/results_bam/${sample_name}.hisat2.sam #转sam为bam文件并进行排序
done
```


**查看bam文件**

```shell
samtools view *.bam|less
```

**什么是bam文件**：https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/04_alignment_quality.html

```shell
ls *bam | while read id ;do (samtools flagstat -@ 10 $id > $(basename $id '.bam').flagstat) ;done  # flagstat 统计比对率
```
注：比对率高低只能说明样本纯度高低，若比对率不高（也有可能是因为 Ref 的问题，如玉米自交系很多但参考基因组很少）也不一定影响后续分析，有**足够数据量**就行

### 2.2 STAR


## 3. 使用IGV查看bam文件

**IGV 怎么用：** 

1. http://software.broadinstitute.org/software/igv/
2. https://www.jianshu.com/p/e5338858dd82

bam文件在导入IGV前需进行**排序及构建索引**

```shell
workpath=/data/FDY_analysis/RNA_seq/FDY/mac3ab/rawdata/hisat2_results
for i in $(ls ${workpath}/*_all.hisat2.bam)
do
    sample_name=`basename $i|sed s/_all.hisat2.bam//g`
    samtools sort -o ${sample_name}_sorted.bam ${sample_name}_all.hisat2.bam
    samtools index ${sample_name}_sorted.bam
done
```



## 4. 区分bam文件的正反链

数据为**链特异性数据**，因此在区分bam文件时需要区分正反链

**不同样本间比较需对wig文件进行标准化**：--normalizeUsing (Possible choices: RPKM, CPM, **BPM**, RPGC, None)

**reference**: https://mp.weixin.qq.com/s?__biz=MzIwODA1MzI4Mg==&mid=2456011255&idx=1&sn=c2d6c19150f2896a8d89fbf07e07522e&chksm=809f94bab7e81dacf31b62eb7159487b112dab496b97619686f01b54f1ce5133733dfa706948&scene=21#wechat_redirect

### bamCoverage

```txt
From: https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html

1. Versions after 2.2
--filterRNAstrand 此参数即可区分正反链

2. Versions before 2.2
For a stranded `paired-end` library
To get the file for transcripts that originated from the forward strand:

# include reads that are 2nd in a pair (128);
# exclude reads that are mapped to the reverse strand (16)
$ samtools view -b -f 128 -F 16 a.bam > a.fwd1.bam

# exclude reads that are mapped to the reverse strand (16) and
# first in a pair (64): 64 + 16 = 80
$ samtools view -b -f 80 a.bam > a.fwd2.bam

# combine the temporary files
$ samtools merge -f fwd.bam a.fwd1.bam a.fwd2.bam

# index the filtered BAM file
$ samtools index fwd.bam

# run bamCoverage
$ bamCoverage -b fwd.bam -o a.fwd.bigWig

# remove the temporary files
$ rm a.fwd*.bam

To get the file for transcripts that originated from the reverse strand:

# include reads that map to the reverse strand (128)
# and are second in a pair (16): 128 + 16 = 144
$ samtools view -b -f 144 a.bam > a.rev1.bam

# include reads that are first in a pair (64), but
# exclude those ones that map to the reverse strand (16)
$ samtools view -b -f 64 -F 16 a.bam > a.rev2.bam

# merge the temporary files
$ samtools merge -f rev.bam a.rev1.bam a.rev2.bam

# index the merged, filtered BAM file
$ samtools index rev.bam

# run bamCoverage
$ bamCoverage -b rev.bam -o a.rev.bw

# remove temporary files
$ rm a.rev*.bam
```

**Versions before 2.2**

```shell
workpath=/data/FDY_analysis/RNA_seq/FDY/mac3ab/rawdata/hisat2_results
env=/home/dell/anaconda2/bin/bamCoverage
for j in $(ls ${workpath}/*.bam) 
do
	i=`basename $j|sed s/.bam//g`
	samtools view -b -f 128 -F 16 ${workpath}/${i}.bam > 	${workpath}/strand_bam/${i}.fwd1.bam;
	samtools view -b -f 80 ${workpath}/${i}.bam > ${workpath}/strand_bam/${i}.fwd2.bam;
	samtools merge -f ${workpath}/strand_bam/${i}.fwd.bam ${workpath}/strand_bam/${i}.fwd1.bam ${workpath}/strand_bam/${i}.fwd2.bam;
	samtools index ${workpath}/strand_bam/${i}.fwd.bam;
	$env -b ${workpath}/strand_bam/${i}.fwd.bam -o ${workpath}/wig/${i}.fwd.bw;
	samtools view -b -f 144 ${workpath}/${i}.bam > ${workpath}/strand_bam/${i}.rev1.bam;
	samtools view -b -f 64 -F 16 ${workpath}/${i}.bam > ${workpath}/strand_bam/${i}.rev2.bam;
	samtools merge -f ${workpath}/strand_bam/${i}.rev.bam ${workpath}/strand_bam/${i}.rev1.bam ${workpath}/strand_bam/${i}.rev2.bam;
	samtools index ${workpath}/strand_bam/${i}.rev.bam;
	$env -b ${workpath}/strand_bam/${i}.rev.bam -o ${workpath}/wig/${i}.rev.bw;
done
```

**Versions after 2.2**


```shell
for i in $(ls ${wkpath_rawdata}/*.R1.fastq.gz)
do
    sample_name=`basename $i .R1.fastq.gz`
    bamCoverage -b ${wkpath_results}/results_bam/${sample_name}.hisat2.bam -o ${wkpath_results}/results_bam/wig/${sample_name}_fw.bw -p 20 -bs 10 --effectiveGenomeSize 119481543 --normalizeUsing BPM --filterRNAstrand forward
    bamCoverage -b ${wkpath_results}/results_bam/${sample_name}.hisat2.bam -o ${wkpath_results}/results_bam/wig/${sample_name}_re.bw -p 20 -bs 10 --effectiveGenomeSize 119481543 --normalizeUsing BPM --filterRNAstrand reverse
done
echo finished

--effectiveGenomeSize 119481543 # 拟南芥基因组碱基除 N 之外大小
--normalizeUsing BPM # BPM == TPM 
```


## 5. 拼接转录本（有参）

在**预测新的lncRNA时**常常需要进行转录本拼接，因为在参考 TAIR10 gff文件中并没有关于lncRNA的全部注释

### cufflink or stringtie or TACO

**cufflink**： http://cole-trapnell-lab.github.io/cufflinks/

**stringtie**： http://ccb.jhu.edu/software/stringtie/

#### 用于预测新的转录本

**5.1 cufflink**
```shell
for i in `ls *.bam`
do
    sample_name=${i%%_*} #拿掉变量i的第一个 _ 及其右边的字符串
    cufflinks -p 12 -u --library-type fr-firststrand 
    -b /data/FDY_analysis/Arabidposis_index_hisat2/Arabidopsis_TAIR10_gene_JYX.fa \
    -g /data/FDY_analysis/Ara_gff_file/TAIR10.GFF3.genes.gff \
    -o ./cufflink_all/${sample_name} $i
done
```

**results**

`cuffcom_split.transcripts.gtf.refmap;cuffcom_split.transcripts.gtf.tmap;genes.fpkm_tracking;isoforms.fpkm_tracking;transcripts.gtf`

**5.2 stringtie**
```shell
for i in `ls ${wkpath_N7_results}/results_bam/*.bam`
do
    sample_name=`basename $i .bam`
    echo ------- start stringtie ${sample_name}
    stringtie ${wkpath_N7_results}/results_bam/${sample_name}.bam -p 20 -G /data/FDY_analysis/Ara_gff_file/TAIR10.GFF3.genes.gff -o ${wkpath_N7_results}/results_bam/stringtie_gtf/${sample_name}_stringtie.gtf
    echo ------ finish ${sample_name} pinjie
done
```

**5.3 TACO**

**reference**： https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5199618/  

注：在组装样本大于10个的时候，stringtie 和 cufflink 更容易将染色体上距离相近的基因组装为一个，推荐 TACO

## 6. 合并转录本

### cufflinks-cuffmerge 

```shell
cuffmerge -p 12 
-s /data/FDY_analysis/Arabidposis_index_hisat2/Arabidopsis_TAIR10_gene_JYX.fa \
-g /data/FDY_analysis/Ara_gff_file/TAIR10.GFF3.genes.gtf  gtf_all_list.txt
```

**gtf_all_list.txt** 包含了所有拼接的 **transcripts.gtf** 的路径

### results

`results: merged_cufflinks.gff`

### stringtie --merge

```shell
echo ------ create a gtf list
ls ${wkpath_N7_results}/results_bam/stringtie_gtf/*.gtf > ${wkpath_N7_results}/results_bam/stringtie_gtf/mergelist.txt
echo ------ begin merge
stringtie --merge -p 20 -G /data/FDY_analysis/Ara_gff_file/TAIR10.GFF3.genes.gff -o ${wkpath_N7_results}/results_bam/stringtie_gtf/stringtie_merged.gtf ${wkpath_N7_results}/results_bam/stringtie_gtf/mergelist.txt
echo ------ merge done
```

## 7. 定量

### cuffdiff or stringtie or featurecounts

**千万注意：若是普通RNA-Seq分析/可变剪接，定量时使用的 gff 文件为原始 gff；若是预测新的转录本，此处是唯一一次需要修改 gff 为 merged.gff 的步骤**

**7.1 cuffdiff 定量及差异分析,包含 FPKM 结果**
```shell
cuffdiff -o diffout_all -p 12 
-b /data/FDY_analysis/Arabidposis_index_hisat2/Arabidopsis_TAIR10_gene_JYX.fa \
--library-type fr-firststrand \
-L col-0,col-3,col-24,mac-0,mac-3,mac-24,skip-0,skip-3,skip-24 \
-u ./merged_asm/merged.gtf \
Col-1-0_368368_all.hisat2.bam,Col-2-0_371371_all.hisat2.bam 
Col-1-3_369369_all.hisat2.bam,Col-2-3_372372_all.hisat2.bam
Col-1-24_370370_all.hisat2.bam,Col-2-24_373373_all.hisat2.bam 
M-1-0_374374_all.hisat2.bam,M-2-0_377377_all.hisat2.bam 
M-1-3_375375_all.hisat2.bam,M-2-3_378378_all.hisat2.bam 
M-1-24_376376_all.hisat2.bam,M-2-24_379379_all.hisat2.bam 
S-1-0_380380_all.hisat2.bam,S-2-0_383383_all.hisat2.bam 
S-1-3_381381_all.hisat2.bam,S-2-3_384384_all.hisat2.bam 
S-1-24_382382_all.hisat2.bam,S-2-24_385385_all.hisat2.bam
```

**7.2 stringtie 定量，包含 coverage,FPKM,TPM 结果，但未进行差异分析**
```shell
echo ------ begin quantify
for i in `ls ${wkpath_N7_results}/results_bam/*.bam`
do
    sample_name=`basename $i .bam`
    stringtie -e -B -p 20 -G ${wkpath_N7_results}/results_bam/stringtie_gtf/stringtie_merged.gtf -o ${wkpath_N7_results}/results_bam/stringtie_gtf/ballgown/${sample_name}/${sample_name}.gtf ${wkpath_N7_results}/results_bam/${sample_name}.bam
done
echo ------ quantify finished


```
结果说明：输出的结果在 bllgown 文件夹下，接下来可使用官方提供的 **ballgown** 进行差异；也可使用 **prepDE.py** 文件提取**原始 counts** 结合 DEseq2 进行差异分析 （推荐后者）

python prepDE.py 会产生两个文件：gene/transcript_count_matrix.csv



**7.3 featurecounts 定量**

reference： http://bioinf.wehi.edu.au/featureCounts/

```shell

/data/software/subread-2.0.0-Linux-x86_64/bin/featureCounts -T 16 -p -s 1 -t exon 
-g transcript_id \
-a ./merged_asm/merged_cufflinks.gtf \
-o all_counts_exon_strand_cufflink_gtf.txt *.bam
```



## 8. 差异分析

差异分析前必看**标准化**问题

**reference**：

1. https://www.jianshu.com/p/cd2888fec66b

2. https://www.jianshu.com/p/248228be3cf0

3. https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html （推荐看这个）

### 8.1	挑选差异基因

若是cuffdiff结果，直接用 **gene.exp.diff** 文件筛选差异基因；若是stringtie/featurecounts结果，使用 **DESeq2**

**阈值**：q/p < 0.05 , **|FC|** >= 1.5/2

```R
library(DESeq2)
countData <- as.matrix(read.csv("transcript_count_matrix.csv", row.names=1)) #counts文件
colData <- read.table(file= "pheno_data.csv", sep=",", row.names=1,header=TRUE) 
#样本对应处理条件
all(rownames(colData) %in% colnames(countData)) #判断countData与colData中的样本名是否一致
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))
dds <- DESeqDataSetFromMatrix(countData, colData, design= ~ treatment)
keep <- rowSums(counts(dds)) >= 10  
#在countData文件里，将每个基因所有的表达量加起来大于10的留下来，并输入到dds文件中
dds <- DESeq(dds)
dds
DESeq2_DEG=na.omit(DEG) #删除差异分析中缺少值的结果
res_n0_vs_c0 = results(dds, contrast=c("treatment", "treatment0", "control0"))
res_n0_vs_c0 = res_n0_vs_c0[order(res_n0_vs_c0$pvalue),]
head(res_n0_vs_c0)
summary(res_n0_vs_c0)
diff_gene_deseq2_res_n0_vs_c0 <-subset(res_n0_vs_c0, padj < 0.05 & abs(log2FoldChange) > 1) #阈值设定
dim(diff_gene_deseq2_res_n0_vs_c0)
head(diff_gene_deseq2_res_n0_vs_c0)
write.csv(diff_gene_deseq2_res_n0_vs_c0,file= "DEG_n0_vs_c0.csv")
```



**8.1.1	差异基因可视化**

```R
###################################### 热图 #######################################

setwd("E:/RNA-Seq/XSW") #设置文件路径
library(pheatmap)
data <- read.table("349_gene.txt",header=T,row.names=1,sep="\t")
pdf(file="file.pdf",width=4,height=6) #输出文件
pheatmap(data,scale="row",cluster_col=TRUE/FALSE,cluster_row=TRUE/FALSE,show_row/colnames=T/F,fontsize=x,breaks= bk,color=colorRampPalette(c("blue","yellow","red"))(100),border_color=FALSE)
dev.off()

参数说明：
	scale="row" 对数据均一化
	cluster_col/row 是否聚类
	show_row/colnames 是否显示行列名
	fontsize 字号
	breaks= bk 设置色度条

# 构建列注释信息 
annotation_col = data.frame( CellType = factor(rep(c("CT1", "CT2"), 5)), Time = 1:5 ) 
rownames(annotation_col) = paste("Test", 1:10, sep = "") 
head(annotation_col)

# 设置色度条在0-5
bk = unique(c(seq(0,5, length=100)))
```


```R
##################################### 火山图 ##################################

rm(list = ls())
data <- read.table("all_gene.txt")
library(ggplot2)
nrDEG[nrDEG$pvalue <0.05 & nrDEG$log2FoldChange >1,ncol(nrDEG)+1]="Up"
nrDEG[nrDEG$pvalue <0.05 & nrDEG$log2FoldChange < -1,ncol(nrDEG)]="Down"
nrDEG[nrDEG$pvalue>=0.05 | 1 > abs(nrDEG$log2FoldChange),ncol(nrDEG)]="Normal"
colnames(nrDEG)[ncol(nrDEG)]="Regulate"
nrDEG$Regulate=factor(nrDEG$Regulate,levels = c("Up","Down","Normal"),order=T)
head(nrDEG)
col=c("red","lightseagreen", "black")
p_volcano = ggplot(nrDEG,aes(x=log2FoldChange,y=-log10(pvalue)))+
  geom_point(aes(color=nrDEG$Regulate),alpha=0.5)+scale_color_manual(values =col)+
  geom_hline(yintercept=c(-log10(0.05)),linetype=4)+
  geom_vline(xintercept=c(-log2(2),log2(2)),linetype=4)+
  theme(
    legend.text = element_text(size = 10, face = "bold"),
    legend.position = 'right',
    legend.key.size=unit(0.5,'cm'),
    axis.text.x=element_text(size = 10,face = "bold", vjust = 0.5, hjust = 0.5),
    axis.text.y=element_text(size = 10,face = "bold", vjust = 0.5, hjust = 0.5),
    axis.title.x = element_text(size = 10,face = "bold", vjust = 0.5, hjust = 0.5),
    axis.title.y = element_text(size = 10,face = "bold", vjust = 0.5, hjust = 0.5),
    panel.background = element_rect(fill = "transparent",colour = "black"), 
    panel.grid.minor = element_line(color="lightgrey",size=0.1),
    panel.grid.major = element_line(color="lightgrey",size=0.1),
    plot.background = element_rect(fill = "transparent",colour = "black")
  )
plot(p_volcano)
```


### 8.2	GO-AgriGO

URL: http://systemsbiology.cau.edu.cn/agriGOv2/

**8.2.1	GO结果可视化**

```R
library(ggpubr) #载入包
data_bar <- read.table("col_up_mac.txt",header = T)	#读入数据，有行名
data$Pathway <- factor(data$Pathway,levels =c("Molecular_Function","Biological_Process"))
#改变横坐标因子顺序
ggbarplot(data_bar,x="Name",y="Score",
          fill="Pathway",rotate=T,
          palette = c("#ADC9D7","#E39E3E"),
          xlab = "",ylab="-log10(FDR)",
          sort.val = "asc",
          group="Pahway",sort.by.groups = TRUE)+
theme_base()
```





## 9. 预测 lncRNA

**鉴定全新的lncRNA**： https://www.jianshu.com/p/5b104830751b

### 9.1 预测新的lncRNA：cuffcompare or gffcompare

**cuffcompare**：http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/

**gffcompare**： http://ccb.jhu.edu/software/stringtie/gffcompare.shtml

**gtf文件与gff文件格式转换**：gffread

```shell
# gff to gtf
gffread my.gff3 -T -o my.gtf
# gtf to gff
gffread merged.gtf -o- > merged.gff3
```

**筛选 lncRNA**

**cuffcompare**
```shell
cuffcompare -r /data/FDY_analysis/Ara_gff_file/TAIR10.GFF3.genes.gtf -p 12 
-o ./cuffcom_split -i ../gtf_all_list.txt
awk '{if($7 >= 0.5 && $10 >1 && $11 > 200){print $0}}' cuffcom.merged_cufflinks.gtf.tmap > filter1.txt
awk '{if($3 == "u" || $3 == "i" || $3 == "x"){print $0}}' filter1.txt > filter2.txt
```

**gffcompare**
```shell
gffcompare -r /data/FDY_analysis/Ara_gff_file/TAIR10.GFF3.genes.gtf -o ./gffcompare/strtcmp ./stringtie_merged_skip.gtf
awk '{if($3 == "u" || $3 == "i" || $3 == "x"){print $0}}' strtcmp.stringtie_merged_skip.gtf.tmap > filter1.txt 
awk '{if($6 >1 && $10 > 200){print $0}}' filter1.txt > filter2.txt # 筛选exon大于1，length大于200的基因
```

### 9.2 **预测lncRNA是否编码**: CPC or CNCI

**lncRNA真真假假**： https://www.jianshu.com/p/0b355662c013

For CPC analysis:

URL: http://cpc.gao-lab.org/

```shell
# 制作 .bed 文件
format：
chr	start	end	ID	strand
1	4662573	4663562	TCONS_00011022	1	-
1	28970700	28974968	TCONS_00018341	1	-
2	11277948	11279070	TCONS_00021738	1	+
2	16737465	16831062	TCONS_00029400	1	-

# 获得 lncRNA 序列
bedtools getfasta -name -s -fi <fasta> -bed <bed/gff/vcf> 
```

For CNCI: (linux)

```shell
python /data/FDY_analysis/tools/CNCI-master/CNCI.py \
-f 346_lncRNA_strand.fa \
-o CNCI \
-m pl \
-p 8 \
```

注：CPC 也可在 linux 下运行，见官方文档说明



## 10. 可变剪接分析

**什么是可变剪接**：https://www.jianshu.com/p/759a5a714aa3

### rMATS（需要有重复）

reference： http://rnaseq-mats.sourceforge.net/

```shell
for i in {0,3,24}
do
	python2.7 /data/FDY_analysis/tools/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py --b1 col_${i}.txt --b2 mac_${i}.txt --gtf /data/FDY_analys
is/RNA_seq/FDY/mac3ab/rawdata/hisat2_results_for_cufflink/results_bam/merged_asm/merged.gtf --od rMATs/col_mac_${i} -t paired --readLength 150 --cstat 0.01 --nthread 8 --libType fr-firststrand
	python2.7 /data/FDY_analysis/tools/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py --b1 col_${i}.txt --b2 skip_${i}.txt --gtf /data/FDY_analy
sis/RNA_seq/FDY/mac3ab/rawdata/hisat2_results_for_cufflink/results_bam/merged_asm/merged.gtf --od rMATs/col_skip_${i} -t paired --readLength 150 --cstat 0.01 --nthread 8 --libType fr-firststrand
done

for j in {col,mac,skip}
do
	python2.7 /data/FDY_analysis/tools/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py --b1 ${j}_0.txt --b2 ${j}_3.txt --gtf /data/FDY_analysis/R
NA_seq/FDY/mac3ab/rawdata/hisat2_results_for_cufflink/results_bam/merged_asm/merged.gtf --od rMATs/${j}_0_3 -t paired --readLength 150 --cstat 0.01 --nthread 8 --libType fr-firststrand
	python2.7 /data/FDY_analysis/tools/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py --b1 ${j}_0.txt --b2 ${j}_24.txt --gtf /data/FDY_analysis/
RNA_seq/FDY/mac3ab/rawdata/hisat2_results_for_cufflink/results_bam/merged_asm/merged.gtf --od rMATs/${j}_0_24 -t paired --readLength 150 --cstat 0.01 --nthread 8 --libType fr-firststrand
	python2.7 /data/FDY_analysis/tools/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py --b1 ${j}_3.txt --b2 ${j}_24.txt --gtf /data/FDY_analysis/
RNA_seq/FDY/mac3ab/rawdata/hisat2_results_for_cufflink/results_bam/merged_asm/merged.gtf --od rMATs/${j}_3_24 -t paired --readLength 150 --cstat 0.01 --nthread 8 --libType fr-firststrand
done
```

**对可变剪接结果可视化：**

```R
library(upsetR)
data <- read.table("data_all.txt")
between <- function(row, min, max){
     newData <- (row["RI_in_mac"] < max) & (row["RI_in_mac"] > min)
}

	upset(data_all_RI,sets=c("RI_0_3","SE_0_3","A3SS_0_3","A5SS_0_3","MXE_0_3","RI_3_24","SE_3_24","A3SS_3_24","A5SS_3_24","MXE_3_24","RI_0_24","SE_0_24","A3SS_0_24","A5SS_0_24","MXE_0_24"),order.by = "freq",keep.order = TRUE,mainbar.y.label = "Gene Intersections", sets.x.label = "Splicing Form", mb.ratio = c(0.6, 0.4), 
		queries = list(list(query = intersects, params = list("RI_0_3")),list(query = between, params=list(1994,1996), color="red", active=T),
		list(query = intersects, params = list( "RI_0_24")),list(query = between, params=list(1995,1997), color="red", active=T),
		list(query = intersects, params = list("RI_3_24")),list(query = between, params=list(1996,1998), color="red", active=T),
		list(query = intersects, params = list("RI_0_3", "RI_0_24")),list(query = between, params=list(1997,1999), color="red", active=T),
		list(query = intersects, params = list("RI_0_3", "RI_3_24")),list(query = between, params=list(1998,2000), color="red", active=T),
		list(query = intersects, params = list("RI_0_24", "RI_3_24")),list(query = between, params=list(1999,2001), color="red", active=T),
		list(query = intersects, params = list("RI_0_3", "RI_0_24","RI_3_24")),list(query = between, params=list(2000,2002), color="red", active=T)))


##upsetR##
upset(movies, queries = list(list(query = intersects, params = list("Drama", 
    "Comedy", "Action"), color = "orange", active = T), list(query = intersects, 
    params = list("Drama"), color = "red", active = F), list(query = intersects, 
    params = list("Action", "Drama"), active = T)))
```

