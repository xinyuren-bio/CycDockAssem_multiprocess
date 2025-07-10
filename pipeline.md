
## 对接前处理
### 1.下载pdb文件
```bash
cd ./Tools
python dowmload_pdb_from_rcsb.py --input_file $input_files output_dir $out_dir
```
### 2.将pdb文件切割转换为fragments
```bash
python pdb2fragments.py
```
### 3.对靶蛋白进行旋转
```bash
python rotate.py --input_pdb ../ALK1/ALK1orig.pdb --output ../ALK1/rotated.pdb --ref1 25.072,-20.981,-45.006 --ref2 26.114,-32.317,-41.141 --ref3 31.015,-30.723,-40.260
```
***输出：***
蛋白质已旋转并保存到: ../ALK1/rotated.pdb
旋转后参考原子的最终坐标:
box1         -6.011   0.000   0.000
box2          6.011   0.000  -0.000
3atm          5.216   5.168   0.000
box2-box1    12.022  -0.000  -0.000

### 4.生成靶蛋白的box文件以及watmap文件
```bash
../SDOCK2.0-restrict/preprocess ../ALK1/rotated.pdb -o ../ALK1/box1.pdb -a ../SDOCK2.0-restrict/ATM -m -6.011,0,0
../SDOCK2.0-restrict/preprocess ../ALK1/rotated.pdb -o ../ALK1/box2.pdb -a ../SDOCK2.0-restrict/ATM -m 6.011,0,0

../SDOCK2.0-restrict/watmap ../ALK1/box1.pdb ../ALK1/box1wat.pdb
../SDOCK2.0-restrict/watmap ../ALK1/box2.pdb ../ALK1/box2wat.pdb
```

### 5.处理对接片段
```bash
cd ../ALK1/fraglib
./fragpreprocess.sh
./fragwater.sh
```
***储存路径***：preprocess：node5:/data1/home/renxinyu/data/pre_fragments_all; node6:/data1/home/renxinyu/data/sdock/pre_fragments_all
***储存路径***：water：node5:/data1/home/renxinyu/data/w_all; node6:/data1/home/renxinyu/data/sdock/w_all

### 6.生成fragstruct,fragwmap以及recordf文件
```bash
python ./fragstruct.py --input_dir ../ALK1/fraglib/ --output_txt ../ALK1/fragstruct --tail ".pdb"
# 因为要保证fragstruct和fragwmap的顺序一致，所以复制一份fragstruct，把里面pre_fragments_all/pre替换成w_all/w
python ./struct2wmap.py --input_file ../ALK1/fragstruct --output_file ../ALK1/fragwmap --old_string preprocessed/pre --new_string watmap/w 
python ./recordf.py --input_file ../ALK1/fragstruct --output_file ../ALK1/box1recordf --prefix_string ../ALK1/dockresult/ALK1box1_
python ./recordf.py --input_file ../ALK1/fragstruct --output_file ../ALK1/box2recordf --prefix_string ../ALK1/dockresult/ALK1box2_
mkdir ../ALK1/dockresult
```

### 7.进行对接
```bash
mkdir ../ALK1/fragstruct_
python ./splitfragstruct.py --input_file ../ALK1/fragstruct --output_prefix ../ALK1/fragstruct_/fragstruct --n_parts 2
mkdir ../ALK1/fragwmap_
python ./splitfragstruct.py --input_file ../ALK1/fragwmap --output_prefix ../ALK1/fragwmap_/fragwmap --n_parts 2
mkdir ../ALK1/recordf_
python ./splitfragstruct.py --input_file ../ALK1/box1recordf --output_prefix ../ALK1/recordf_/box1recordf --n_parts 2
python ./splitfragstruct.py --input_file ../ALK1/box2recordf --output_prefix ../ALK1/recordf_/box2recordf --n_parts 2
# box1
python ./dock.py --num_tasks 2 --fixed_pdb ../ALK1/box1.pdb --fragstruct_dir ../ALK1/fragstruct_/ --fixed_wat_pdb ../ALK1/box1wat.pdb --fragwatmap_dir ../ALK1/fragwmap_/ --output_record_dir ../ALK1/recordf_/ --output_prefix box1recordf
# box2
python ./dock.py --num_tasks 2 --fixed_pdb ../ALK1/box2.pdb --fragstruct_dir ../ALK1/fragstruct_/ --fixed_wat_pdb ../ALK1/box2wat.pdb --fragwatmap_dir ../ALK1/fragwmap_/ --output_record_dir ../ALK1/recordf_/ --output_prefix box2recordf
```

### 8.每个box选取前20
```bash
python ./getbestdocking.py --recordf ../ALK1/box2recordf --output_path ../ALK1/box2bestdocking --top_n 20
python ./getbestdocking.py --recordf ../ALK1/box1recordf --output_path ../ALK1/box1bestdocking --top_n 20
```

### 9.生成docking结果
```bash
../utility/genBuildcommand ../ALK1/box1bestdocking ../ALK1/genbox1frag.sh ../SDOCK2.0-restrict/build ../ALK1/box1.pdb ../ALK1/dockmodel ../SDOCK2.0-restrict/so3layer_648.qua 0,0,0
../utility/genBuildcommand ../ALK1/box2bestdocking ../ALK1/genbox2frag.sh ../SDOCK2.0-restrict/build ../ALK1/box2.pdb ../ALK1/dockmodel ../SDOCK2.0-restrict/so3layer_648.qua -12.022,0,0
mkdir ../ALK1/dockmodel
bash ../ALK1/genbox1frag.sh
bash ../ALK1/genbox2frag.sh
ls ../ALK1/dockmodel/SDOCK_ALK1box1* > ../ALK1/box1dockfrag
ls ../ALK1/dockmodel/SDOCK_ALK1box2* > ../ALK1/box2dockfrag
python ./getfragpair.py ../ALK1/box1dockfrag ../ALK1/box2dockfrag -56 8 ../ALK1/ALK1dockfragpair
```

### 10.find linker
```bash
mkdir ../ALK1/fraglinking
python fraglink.py
mkdir ../ALK1/AssembledCyc
# 对c源码进行了修改
python ./assemblecyc.py --fraglinking_dir ../ALK1/fraglinking/ --assembled_cyc_dir ../ALK1/AssembledCyc/
../utility/SelectBuildcomplex ../ALK1/AssembledCyc/ALK1Cyc_4.pdb -55.0 2 -55.0 3 ../ALK
1/complexmodel/ ../ALK1/box1.pdb


