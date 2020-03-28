# linux实用命令

```shell
# 统计当前目录下文件/目录的数目

ls -l | grep "^-" | wc -l #统计文件
ls -l | grep "^d" | wc -l #统计目录

# 批量将目录下的所有 .R文件名写入到 .sh

for i in `ls *.R`;do echo "Rscript $i &";done > Rscript.sh
```

#### #，% 修改变量

```shell
先赋值一个变量为一个路径，如下：
file=/dir1/dir2/dir3/my.file.txt

#拿掉第一条 / 及其 左边 的字符串
echo ${file#*/}
dir1/dir2/dir3/my.file.txt

#拿掉最后一条 / 及其左边的字符串
echo ${file##*/}
my.file.txt

#拿掉第一个 . 及其左边的字符串
echo ${file#*.}
file.txt

#拿掉最后一个 . 及其左边的字符串
echo ${file##*.}
txt

#拿掉最后一条 / 及其右边的字符串
echo ${file%/*}
/dir1/dir2/dir3

#拿掉第一条 / 及其右边的字符串
echo ${file%%/*}
(空值)

#拿掉最后一个 . 及其右边的字符串
echo ${file%.*}
/dir1/dir2/dir3/my.file

#拿掉第一个 . 及其右边的字符串
echo ${file%%.*}
/dir1/dir2/dir3/my

注：
# 是去掉左边
% 是去掉右边
*是用来匹配不要的字符，也就是想要去掉的那部分
```

#### awk,sed,grep

```shell
# 同时匹配几个条件，如搜索txt文件中同时包含col-0，col-3的行

grep col-0 file.txt | grep col-3
grep -E 'col-0.*col-3' file.txt
sed -n '/col-0/{/col-3/p}' file.txt
awk '/col-0/&&/col-3/{ print $0 }' file.txt

# 匹配ABC 或 abc

sed -n '/\(ABC\|abc\)/p'
awk '/ABC/||/abc/{ print $0 }'
grep -E '(ABC|abc)' 或 egrep 'ABC|abc'
```


