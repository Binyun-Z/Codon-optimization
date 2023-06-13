# Codon-optimization
The software is in a specific host chassis proteins perform codon optimization  
该程序参考论文：[Rational design and construction of multi-copy biomanufacturing islands in mammalian cells](https://doi.org/10.1093/nar/gkab1214)，在其提供的源代码基础上提供了一个python3版本的程序  
该软件为以特定宿主底盘的蛋白质执行密码子优化
密码子使用表可从[此处](http://www.kazusa.or.jp/codon/)或者[金斯瑞官网](https://www.genscript.com/tools/codon-frequency-table)获得
## 使用方法：
```shell
python GD.py -p protein.fasta -n 10
```
## 参数说明：
`-p` : 要执行密码子优化的蛋白质氨基酸序列的fasta文件位置，默认为`test.fasta` ,可以指定该文件所在的位置，例如`./input/test.fasta`  
`-n` : 需要生成的不同基因序列数量

## 结果说明
 正确执行软件之后会输出名为：`codingSequenceVariants.csv`的结果文件，`sequence`列为优化后的基因序列，`CHI`为基因的密码子适应性指数 ，CAI 计算为每个密码子在基因序列长度 L 上的相对适应性的几何平均值，`GC`为基因序列中碱基`G`和碱基`C`的含量。  
![image](https://user-images.githubusercontent.com/35862583/234791351-31431315-7428-407f-959b-866d93601109.png)  
![image](https://user-images.githubusercontent.com/35862583/234791493-08b4ff42-4b1f-4d3f-aa05-83a8e2f5e198.png)    
同时，软件会输出生成的基因序列组中的平均、最大、最小汉明距离，以及序列组中两个基因序列的最远同源性。
（Mean, minimum and maximum hamming distances in the sequence）
（Longest stretch of homology between any two sequences (in bp)）
