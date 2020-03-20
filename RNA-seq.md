

# RNA-seq	procedure	for	Yanglab



## 数据：MAC3A/B,SKIP

### 类型：paired, 150x, fr-firststrand

`rawdata: Col-1-0_368368_all.R1.fastq.gz;Col-1-0_368368_all.R2.fastq.gz`



## 1. 质控

### fastqc

```shell
fastqc /data/FDY_analysis/RNA_seq/FDY/mac3ab/rawdata/*.fastq.gz -o fastqc
multiqc *fastqc.zip --ignore *.html
```

![1578378576(1)](D:\RNA-Seq-fdy\fastqc\1578378576(1).png)![1578378509(1)](D:\RNA-Seq-fdy\fastqc\1578378509(1).png)

## 2. 比对

### hisat2(√) or STAR

2.1	添加环境变量：

```shell
vi ~/.bashrc #永久修改环境变量，可直接调用 hisat2
export PATH=/data/sly/tools/hisat2-2.0.4/hisat2:$PATH
source ~/.bashrc #使修改生效
```

2.2	构建索引：

```shell
hisat2-build -p 10 Zea_mays.AGPv4.dna.toplevel.fa genome
```

2.3	比对：

```shell
wkpath1=/data/FDY_analysis/RNA_seq/FDY/mac3ab/rawdata/cleandata #设置工作路径
wkpath2=/data/FDY_analysis/RNA_seq/FDY/mac3ab/rawdata/hisat2_results_for_cufflink
for i in $(ls ${wkpath1}/*.R1.fastq.gz) 
do
    sample_name=`basename $i|sed s/.R1.*//g` #取文件名（sed 用于将 /.R1.*/ 前后的字符替换为空格）
    hisat2 -p 12 \ #线程
    --dta-cufflinks \ #设置此参数用于后续分析结合cfflinks
    --rna-strandness RF \ #链特异性
    -x /data/FDY_analysis/Arabidposis_index_hisat2/genome \ #hisat2-build 构建的索引文件
    -1 ${wkpath1}/${sample_name}.R1.fastq.gz \
    -2 ${wkpath1}/${sample_name}.R2.fastq.gz \
    -S ${wkpath2}/results_bam/${sample_name}.hisat2.sam #输出为sam文件格式
    samtools sort -@ 12 -o ${wkpath2}/results_bam/${sample_name}.hisat2.bam ${wkpath2}/results_bam/${sample_name}.hisat2.sam #转sam为bam文件并进行排序
done
```

### results

`results: Col-1-0_368368_all.hisat2.bam`

查看bam文件

```shell
samtools view *.bam|less
```

![bam](D:\RNA-Seq-fdy\markdown_file\bam.png)



## 3. 使用IGV查看bam文件

bam文件在导入IGV前需进行排序及构建索引

```shell
workpath=/data/FDY_analysis/RNA_seq/FDY/mac3ab/rawdata/hisat2_results
for i in $(ls ${workpath}/*_all.hisat2.bam)
do
    sample_name=`basename $i|sed s/_all.hisat2.bam//g`
    samtools sort -o ${sample_name}_sorted.bam ${sample_name}_all.hisat2.bam
    samtools index ${sample_name}_sorted.bam
done
```

IGV结果：

![IGV](D:\RNA-Seq-fdy\markdown_file\IGV.jpg)



## 4. 区分bam文件的正反链

### bamCoverage

```txt
From: https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html

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



## 5. 拼接转录本（有参）

### cufflink(√) or stringtie

#### 用于预测新的转录本

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

### results

`cuffcom_split.transcripts.gtf.refmap;cuffcom_split.transcripts.gtf.tmap;genes.fpkm_tracking;isoforms.fpkm_tracking;transcripts.gtf`



## 6. 合并转录本

### cuffmerge

```shell
cuffmerge -p 12 
-s /data/FDY_analysis/Arabidposis_index_hisat2/Arabidopsis_TAIR10_gene_JYX.fa \
-g /data/FDY_analysis/Ara_gff_file/TAIR10.GFF3.genes.gtf  gtf_all_list.txt
```

**gtf_all_list.txt** 包含了所有拼接的 **transcripts.gtf** 的路径

### results

`results: merged_cufflinks.gff`



## 7. 定量

### cuffdiff(√) or featurecounts

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

or

```shell
/data/software/subread-2.0.0-Linux-x86_64/bin/featureCounts -T 16 -p -s 1 -t exon 
-g transcript_id \
-a ./merged_asm/merged_cufflinks.gtf \
-o all_counts_exon_strand_cufflink_gtf.txt *.bam
```



## 8. 差异分析

### 8.1	挑选差异基因

若是cuffdiff结果，直接用 **gene.exp.diff** 文件筛选差异基因；若是featurecounts结果，使用DESeq2

阈值：q/p < 0.05 , |FC| >= 1.5/2

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



8.1.1	差异基因可视化

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

<img src="D:\RNA-Seq-fdy\markdown_file\1578382811(1).png" alt="1578382811(1)" style="zoom: 67%;" />

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

![1578383812(1)](D:\RNA-Seq-fdy\markdown_file\1578383812(1).png)

### 8.2	GO-AgriGO

URL: http://systemsbiology.cau.edu.cn/agriGOv2/

8.2.1	GO结果可视化

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



![1578382320(1)](D:\RNA-Seq-fdy\markdown_file\1578382320(1).png)



## 9. 预测 lncRNA

### cuffcompare

gtf文件与gff文件格式转换

```shell
# gff to gtf
gffread my.gff3 -T -o my.gtf
# gtf to gff
gffread merged.gtf -o- > merged.gff3
```

筛选 lncRNA

```shell
cuffcompare -r /data/FDY_analysis/Ara_gff_file/TAIR10.GFF3.genes.gtf -p 12 
-o ./cuffcom_split -i ../gtf_all_list.txt
awk '{if($7 >= 0.5 && $10 >1 && $11 > 200){print $0}}' cuffcom.merged_cufflinks.gtf.tmap > filter1.txt
awk '{if($3 == "u" || $3 == "i" || $3 == "u"){print $0}}' filter1.txt > filter2.txt
```

预测lncRNA是否编码: CPC or CNCI

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

### rMATS（需要有重复）

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

对可变剪接结果可视化：

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

