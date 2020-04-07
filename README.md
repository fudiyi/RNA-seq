# RNA-seq for yanglab

此文档包含了两种方法：hisat2 + cufflinks + cuffdiff； hisat2 + stringie + DEseq2 对高通量测序结果进行转录组分析，
对两种结果得到差异基因比较后发现结果相差较大（50% - 80%）。原因是两种差异分析方法的归一化原理不同（详细请自行查阅官方文档），
因此建议根据自己所需结果选择合适的方法。
