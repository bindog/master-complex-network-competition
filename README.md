# 复杂网络中的关键节点发现（数据城堡复杂网络大师赛第四名代码）

## 大赛链接和数据资源

- [大赛链接](http://www.dcjingsai.com/common/cmpt/%E5%A4%A7%E5%B8%88%E8%B5%9B_%E7%AB%9E%E8%B5%9B%E4%BF%A1%E6%81%AF.html)
8个网络下载地址：[百度网盘链接](https://pan.baidu.com/s/1o8abtto)

## 评价算法

`judge.c`是用C实现的评价算法，实现思路与[官方版本](http://share.pkbigdata.com/ID.4407/Master_algorithm)一致

编译

```shell
g++ -std=c++11 judge.c -o judge
```

运行

```shell
./judge result.csv
```

注意`result.csv`的格式，因为是自己使用的，不用考虑节点数量错误之类的问题，所以没有加异常处理。

`result.csv`格式示例如下：

```shell
model1,34,21,...,8997,98
model2,394,340,...,2231,3
model3,7785,774,...,8821,778
model4,1123,6653,...,2213,6
real1,33214,776,...,89,99
real2,542,112,...,1123,82
real3,1277,123,...,3345,87
real4,1233,3567,...,90,900
```

## 说明

把数据下载到本地，新建一个`networks`文件夹，将8个网络的csv文件放入该文件夹，代码必须和该文件夹在同一级目录下。

`greedyPageRank.py`是基于贪心策略的`PageRank`算法，依赖`python-igraph`库，使用方法详见帮助选项。

`seqOptimize.c`是序列优化算法

编译

```shell
g++ -std=c++11 -fopenmp seqOptimize.c -o seqOptimize
```

使用方法详见帮助选项。

注意：`greedyPageRank.py`的输出文件可以作为`seqOptimize.c`的输入文件，得到所有网络的节点重要性序列后，需要合并到一个文件中，再用`judge.c`计算结果。