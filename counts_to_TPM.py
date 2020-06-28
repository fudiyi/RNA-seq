
import pandas as pd
import numpy as np


data = pd.read_table(r'D:\yanglab_data\RNA-Seq\FDY\reads.txt')



def get_RPKM():

    counts_sum_all = {}
    sample_rpkm = {}

    for i in data.keys()[2:]:

        counts_sum_all['%s' % i] = data[i].sum() # 对各个样本的counts求和

        sample_rpkm['%s' % i] = pd.Series((data[i]/(counts_sum_all[i]/1000000))/(data['len']/1000)) # 将样本的各个基因count除以10^6得到RPM,再将RPM除以基因长度

    index = pd.Series(data['gene'])

    data_rpkm = pd.concat([index,sample_rpkm['sampleA'],sample_rpkm['sampleB'],sample_rpkm['sampleC'],sample_rpkm['sampleD']],axis=1)

    data_rpkm.rename({0:'sampleA',1:'sampleB',2:'sampleC',3:'sampleD'},axis='columns',inplace=True)

    data_rpkm.loc['rpkm_sum'] = data_rpkm.apply(lambda x: x.sum())

    print(data_rpkm)

get_RPKM()


def get_TPM():

    sample_tpm = {}

    for i in data.keys()[2:]:

        gene_scale = data[i]/(data['len']/1000)

        sample_tpm['%s' % i] = gene_scale/(gene_scale.sum()/10000000)

    index = pd.Series(data['gene'])

    data_tpm = pd.concat([index,sample_tpm['sampleA'],sample_tpm['sampleB'],sample_tpm['sampleC'],sample_tpm['sampleD']],axis=1)

    data_tpm.rename({0:'sampleA',1:'sampleB',2:'sampleC',3:'sampleD'},axis='columns',inplace=True)

    data_tpm.loc['tpm_sum'] = data_tpm.apply(lambda x: x.sum())

    print(data_tpm)

get_TPM()



