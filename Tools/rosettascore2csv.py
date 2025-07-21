import os
import glob
import pandas as pd
import re

def extract_sc_data(file_path):
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()

        score_header_line = None
        score_data_line = None
        found_first_score_prefix = False 
        for i, line in enumerate(lines):
            line_stripped = line.strip()

            if line_stripped.startswith("SCORE:"):
                if score_header_line is None:
                    score_header_line = line_stripped
        
                elif score_data_line is None:
                    if any(char.isdigit() for char in line_stripped.replace("SCORE:", "")):
                        score_data_line = line_stripped
                       
                
            if score_header_line is not None and score_data_line is not None:

                break

        if not score_header_line or not score_data_line:
            print("警告: 文件 '{}' 中未找到完整的 SCORE 头部或数据行。".format(file_path))
            return None
        raw_headers = [h for h in re.findall(r'\S+', score_header_line.replace("SCORE:", ""))]
        
        raw_values = [v for v in re.findall(r'\S+', score_data_line.replace("SCORE:", ""))]
        
        if len(raw_headers) != len(raw_values):
            print("警告: 文件 '{}' 中 SCORE 头部和数据列数不匹配。头部 {} 列，数据 {} 列。跳过。".format(file_path, len(raw_headers), len(raw_values)))
            return None

        result_dict = {}
        for header, value in zip(raw_headers, raw_values):
            try:
                result_dict[header] = float(value)
            except ValueError:
                result_dict[header] = value 

        original_description_value = result_dict.get('description', '')

        if 'description' in result_dict and isinstance(result_dict['description'], str):
            current_description_value = result_dict['description']
            
            match = re.match(r'(.*?)_design_(\d+)(.*)', current_description_value)
            
            if match:
                pdb_basename_part = match.group(1) 
                resfile_index_part = match.group(2) 
                
                new_formatted_description = "{}_{}".format(pdb_basename_part, resfile_index_part)
                result_dict['description'] = new_formatted_description 
          
        result_dict['original_pdb_basename'] = original_description_value 
        
        return result_dict

    except Exception as e:
        print("错误: 处理文件 '{}' 时发生异常：{}".format(file_path, e))
        return None

def extract_sequence_from_resfile(resfile_path):
    sequence = []
    found_start = False
    
    if not os.path.isfile(resfile_path):
        return ""

    try:
        with open(resfile_path, 'r') as f:
            for line in f:
                line_stripped = line.strip()
                
                if not found_start:
                    if line_stripped == "start":
                        found_start = True
                    continue 

                parts = line_stripped.split()
                if len(parts) >= 1:
                    sequence.append(parts[-1])
        return "".join(sequence)
    except Exception as e:
        print("  警告: 从序列文件 '{}' 提取序列失败。错误：{}".format(resfile_path, e))
        return "" 

def process_all_sc_files(base_folder, sequence_resfile_folder, output_csv_path):
    all_extracted_data = []
    
    if not os.path.isdir(base_folder):
        print("错误: 基础文件夹 '{}' 不存在。请检查路径。".format(base_folder))
        return
    if not os.path.isdir(sequence_resfile_folder):
        print("错误: 序列文件文件夹 '{}' 不存在。请检查路径。".format(sequence_resfile_folder))
        return

    print("--- 开始处理 .sc 文件 ---")
    print("遍历文件夹: {}".format(base_folder))
    print("序列文件来源: {}".format(sequence_resfile_folder))

    subfolders = [d for d in os.listdir(base_folder) if os.path.isdir(os.path.join(base_folder, d))]

    if not subfolders:
        print("警告: 基础文件夹 '{}' 下没有找到任何子文件夹。".format(base_folder))
        return

    for subfolder_name in subfolders:
        subfolder_path = os.path.join(base_folder, subfolder_name)
        print("正在处理子文件夹: {}".format(subfolder_path))

        sc_files = glob.glob(os.path.join(subfolder_path, '*.sc'))

        if not sc_files:
            print("  警告: 子文件夹 '{}' 中没有找到任何 .sc 文件。".format(subfolder_name))
            continue

        for sc_file_path in sc_files:
            print("  提取文件: {}".format(os.path.basename(sc_file_path)))
            data = extract_sc_data(sc_file_path)
            
            if data:
                if 'description' in data: 
                    modified_description = data['description']
                    seq_resfile_path = os.path.join(sequence_resfile_folder, "{}.resfile".format(modified_description))
                    
                    extracted_sequence = extract_sequence_from_resfile(seq_resfile_path)
                    data['sequence'] = extracted_sequence 
                else:
                    data['sequence'] = "" 
                all_extracted_data.append(data)
    
    if not all_extracted_data:
        print("没有从任何 .sc 文件中提取到有效数据。未生成 CSV。")
        return
    df = pd.DataFrame(all_extracted_data)

    final_cols = []
    if 'description' in df.columns:
        final_cols.append('description')
    if 'sequence' in df.columns:
        final_cols.append('sequence')
    
    if 'total_score' in df.columns:
        final_cols.append('total_score')
    other_cols = [col for col in df.columns if col not in ['description', 'sequence', 'total_score']]
    final_cols.extend(sorted(other_cols))
    
    df = df[final_cols]

    try:
        df.to_csv(output_csv_path, index=False)
        print("\n--- 处理完成 ---")
        print("所有提取的数据已保存到: {}".format(output_csv_path))
    except Exception as e:
        print("错误: 保存 CSV 文件 '{}' 时发生异常：{}".format(output_csv_path, e))

if __name__ == "__main__":
    BASE_FOLDER_TO_SCAN = "/mnt/data/renxinyu/rosetta/pratice/TNfa/rosetta3rd_ouputs_protein" # 记录了结果的文件夹路径
    SEQUENCE_RESFILE_FOLDER = "/mnt/data/renxinyu/rosetta/pratice/TNfa/to3rosetta_resfiles_protein"  # resfile文件夹路径
    OUTPUT_CSV_REPORT = "3rd_proteinMPNN.csv"   # 生成后csv文件的保存路径

    process_all_sc_files(BASE_FOLDER_TO_SCAN, SEQUENCE_RESFILE_FOLDER, OUTPUT_CSV_REPORT)