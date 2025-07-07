# CycDockAssem_multiprocess

源项目 [https://github.com/victorPKU/CycAssem]

该仓库是对源项目进行了部分改善的版本，旨在通过并行化提升性能并引入多个自动化脚本提升效率。注意主要流程请仍然使用源项目进行操作，当遇到耗时较长的步骤时，可以使用本仓库进行操作。

ps：大家如果有任何问题，随时在issues里提问QaQ.

## 重写
1. utility/targetframe.c:更改为python文件，替代之前导入frame.pdb的方式，只需要输入定义box1，box2的三个三维坐标即可对靶蛋白进行旋转。

2. fraglink：对fraglinl.c进行了重写，改善了数据量过大时报错的问题。
(1) 添加多进程并行，极大提升效率
(2) 添加参数修改，支持个性化处理


## 完善
1. SDOCK2.0-restrict/sdock：对接部分代码进行完善，支持并行，极大提升对接效率。

2. SDOCK2.0-restrict/watmap: 添加多进程，提升效率。

3. fragtermgeo/calTerminalgeo.c:之前代码计算后写入单独的geo文件，IO频繁，现改为先存为csv文件，并添加多线程，极大提升效率。

## 新增

1. 添加下载pdb文件的爬虫脚本。
2. 添加将pdb转换为fragments的脚本。