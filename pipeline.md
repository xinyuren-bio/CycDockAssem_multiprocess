## 编译：
### 编译sdock
```bash
cd ~/fragdesign
mkdir -p compiled_libs/fftw_install
mkdir -p compiled_sources

wget http://www.fftw.org/fftw-3.3.10.tar.gz -P compiled_sources/
tar -xzf compiled_sources/fftw-3.3.10.tar.gz -C compiled_sources/
cd compiled_sources/fftw-3.3.10/
./configure --prefix=~/fragdesign/compiled_libs/fftw_install --enable-shared --enable-static
make
make install

cd ~/fragdesign/SDOCK2.0-restrict
export FFTW_HOME=/data1/home/renxinyu/fragdesign/compiled_libs/fftw_install
export LD_LIBRARY_PATH="${FFTW_HOME}/lib" 
进入SDOCK2.0-restrict文件夹，找到MAKEFILE文件，将FFTW_HOME修改为上面配置的FFTW_HOME。
eg：FFTW_DIR        = /home/changsheng/Downloads/Softwares/MDL/SDOCK1p0/fftw-3.3.8 -> FFTW_DIR        = /data1/home/renxinyu/fragdesign/compiled_libs/fftw_install
make
cd ~/fragdesign/
```
### 编译其他工具：
```bash
cd ./utility
bash ./readme.sh 
```


# 运行
## 一键运行版本：
```bash
python ./pipeline.py ./config.yml
```
需要在ALK1中上传一个必要文件，即只包含靶蛋白的pdb文件。
在config.yml修改参数：
1. project_base_dir，修改为当前项目所在的目录。eg: "/data1/home/renxinyu/CycDockAssem_overwrite/"
2. target_protein_name，文件夹的名称，用来后续生成其他必要文件是作为prefix。eg: "ALK1"
3. sdock_base_dir，存储sdockc代码的文件夹。eg: "/data1/home/renxinyu/sdock/SDOCK2.0-restrict"
4. preprocessed_fragments_path，存储预处理片段的文件夹。eg: node5:"/data1/home/renxinyu/data/pre_fragments_all";node6:"/data1/home/renxinyu/data/sdock/pre_fragments_all" 注意：pre_fragments_all这个名称与后面代码存在耦合，不能更改。
5. water_fragments_path，存储水片段的文件夹。eg: node5:"/data1/home/renxinyu/data/water_fragments_all";node6:"/data1/home/renxinyu/data/sdock/water_fragments_all" 注意：water_fragments_all这个名称与后面代码存在耦合，不能更改。
6. ref1: "25.072,-20.981,-45.006"
   ref2: "26.114,-32.317,-41.141"
   ref3: "31.015,-30.723,-40.260"，要定义config.yml中step3三个原子的坐标，ref1-box1中心，ref2-box2中心，ref3-定义y轴。
7. 其他参数都是和dock，link相关的参数，如果没有特殊需求，使用默认值就可以。


## 逐步运行版本：
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
```

### 11.assemble cyc
```bash
mkdir ../ALK1/AssembledCyc
# 对c源码进行了修改
python ./assemblecyc.py --fraglinking_dir ../ALK1/fraglinking/ --assembled_cyc_dir ../ALK1/AssembledCyc/
../utility/SelectBuildcomplex ../ALK1/AssembledCyc/ALK1Cyc_4.pdb -55.0 2 -55.0 3 ../ALK1/complexmodel/ ../ALK1/box1.pdb
```

### 12.sequence design
```bash
# 使用MPNN对上面产生的pdb文件进行设计，得到设计好的.fa文件夹
# 在gene_resfile.py中修改参数，将.fa文件夹路径修改为上面生成的.fa文件夹路径，并运行，得到Rosetta所需要的resfile文件
python ./gene_resfile.py
# 执行完上面的命令后，我们就可获得resfile文件，有了resfile文件，就可以开始进行序列设计了，可以执行这个bash脚本，里面有一些参数需要修改，可以根据脚本内部的参数注解进行相应的修改（参数都在脚本上方），修改完成后就可以直接运行。
bash ./rosetta.sh
# 运行结束，获得分数后，由于Rosetta产生的格式不方便绘图，这里提供了一个Python脚本可以帮助我们把得分提取到csv文件中，需要提供三个参数：i. 结果的保存路径  ii. refile文件的路径  iii. csv文件的保存路径及名称
python ./rosetta_score2csv.py
# 通过这行以上命令我们就可以获得一个csv文件，里面保存了 Rosetta 的得分，如果我们需要进行绘图比对不同序列他们dG之间分数的差异，我们可以执行以下命令，并且在python内部文件修改好相关的参数：
python ./bin.py
