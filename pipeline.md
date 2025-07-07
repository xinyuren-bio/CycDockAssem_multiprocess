大家可以跟着这个流程进行操作：
### 1.下载pdb文件
```bash
cd ./Tools
python dowmload_pdb_from_rcsb.py --input_file $input_files output_dir $out_dir
```
### 2.将pdb文件切割转换为fragments
```bash
python pdb2fragments.py
```
### 3.执行对接


