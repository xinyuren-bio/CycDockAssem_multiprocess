# CycDockAssem_overwrite

源项目 [https://github.com/victorPKU/CycAssem]

本仓库是一个对源项目进行部分重写及改善的版本，旨在通过并行化提升性能并引入多个自动化脚本提升效率。

## 重写
1. utility/targetframe.c:更改为python文件，替代之前导入frame.pdb的方式，只需要输入定义box1，box2的三个三维坐标即可。

2. fraglink：对fraglinl.c进行了重写，改善了数据量过大时报错的问题。
(1) 添加多进程并行，极大提升效率
(2) 添加参数修改，支持个性化处理


## 完善
1. SDOCK2.0-restrict/sdock：对接部分代码进行完善，支持并行，极大提升对接效率。

2. SDOCK2.0-restrict/watmap: 添加多进程，提升效率。

3. fragtermgeo/calTerminalgeo.c:之前代码计算后写入单独的geo文件，IO频繁，现改为先存为csv文件，并添加多线程，极大提升效率。

## 新增

1. 之前没有下载pdb的爬虫脚本，现已添加。
2. 没有将pdb转换为fragments的脚本，现已添加。