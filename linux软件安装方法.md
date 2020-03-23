# conda安装与环境配置

### 在linux下载安装软件目前有两种方法

```shell
# 1.Installation with pip

pip install /PATH/filename
```

#### 此方法安装需注意软件所需安装环境且安装后需自己配置环境变量

```shell
#永久修改为全局环境变量，可直接调用

vi ~/.bashrc 
export PATH=/home/dell/miniconda3/bin/conda:$PATH #找到安装目录，添加环境变量
source ~/.bashrc #使修改生效
```



```shell
# 2. Installation with conda

##################################安装miniconda3############################
wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod 777 Miniconda3-latest-Linux-x86_64.sh #修改权限
sh Miniconda3-latest-Linux-x86_64.sh
vi ~/.bashrc #添加conda环境
source ~/.bashrc

#################################激活conda###################################
source ./activate or source /home/dell/miniconda3/bin/activate base #激活conda安装目录bin下的activate，linux命令行前会出现base，则表示已进入conda默认环境
conda env list #查看已有环境
conda create -n python2 python=2 #创建名为python2的环境
conda activate python2 #激活该环境
conda deactivate #关闭环境
vim ~/.condarc #有时候下载速度会很慢，更需改镜像

#进行转录组分析时可设置一个专门的RNA-Seq环境变量
conda install -n RNA-Seq
conda install -y hisat2 cutadapt star #在RNA-Seq环境下安装一系列分析所需软件

conda install -n python2 -y cufflinks #当python环境不对时也可使用-n指定特定环境进行安装
```



#### 推荐先尝试使用pip安装，可以学到很多环境变量配置的方法，但是很多问题确实无法单纯通过pip来解决，毕竟Linux的环境比较复杂且有时需要root权限，此时就可以用conda来解决