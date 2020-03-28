# python-筛选差异基因

#### 此脚本用于对RNA-Seq的cuffdiff结果进行批量差异基因筛选

```python
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#fdy 2020-3-28

#自定义一个函数，因为差异分析时经常为对照（Col）和实验组（mut）进行比较，因此定义函数时设定两个参数

def getcom(x,y):
    sample = str(x) + '_' + str(y) + '.txt'
    sample_DEG = str(x) + '_' + str(y) + '.DEG.txt'
    fl = open(sample, 'w') #用于存放结果
    fl_DEG = open(sample_DEG, 'w')
    fl.write('gene' + '\t' + 'sample_1' + '\t' + 'sample_2' + '\t' + 'value1' + '\t' + 'value2' + '\t' + 'log2(foldchange)' + '\t' + 'p' + '\t' + 'q' + '\n') #写入行名
    fl_DEG.write('gene' + '\t' + 'sample_1' + '\t' + 'sample_2' + '\t' + 'value1' + '\t' + 'value2' + '\t' + 'log2(foldchange)' + '\t' + 'p' + '\t' + 'q' + '\n')
    with open(r'D:\yanglab_data\gene_exp.diff','r') as a:
        for i in a:
            ii = i.strip().split('\t')
            if ii[4] == x and ii[5] == y: #判断cuffdiff的结果文件是否包含对应Col和mut
                fl.write(str(ii[2]) + '\t' + str(ii[4]) + '\t' + str(ii[5]) + '\t' + str(ii[7]) + '\t' + str(ii[8]) + '\t' + str(ii[9]) + '\t' + str(ii[11]) + '\t' + str(ii[12]) + '\n')
                #fl.close() #不能设置此参数，否则在后面调用函数时无法读取到文件
            if ii[4] == x and ii[5] == y and float(ii[11]) < float(0.05) and abs(float(ii[9])) >= int(1): #筛选差异基因（p<0.05,|log2FC|>=1）
                fl_DEG.write(str(ii[2]) + '\t' + str(ii[4]) + '\t' + str(ii[5]) + '\t' + str(ii[7]) + '\t' + str(ii[8]) + '\t' + str(ii[9]) + '\t' + str(ii[11]) + '\t' + str(ii[12]) + '\n')
    print('Get' + '\t' + sample)
    print('Get' + '\t' +sample_DEG)

com = [('q1','q2'),('q1','q3'),('q2','q3')] #设置所需比较，如('col-0','col-3')
for i in com:
    getcom(*i) # *i代表对com列表中元组的参数进行分配，即'q1'分配给x，'q2'分配给y

print('done')
```

