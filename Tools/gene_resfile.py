import os
import re
import glob

def process_proteinmpnn_output_to_resfile(
    input_file_path,
    output_folder
):
    """
    从.fa文件中读取序列，
    并为每个序列生成一个.resfile文件
    输出文件名将基于输入文件的名称，后接序号和 .resfile 后缀。
    Args:
        input_file_path (str): 包含 ProteinMPNN 生成序列的输入文本文件路径。
        output_folder (str): 存放生成的新文件的目标文件夹路径。
    """
    if not os.path.exists(input_file_path):
        print(f"错误：输入文件 '{input_file_path}' 不存在。请检查路径。")
        return
    base_input_filename = os.path.splitext(os.path.basename(input_file_path))[0]
    sequence_pattern = re.compile(r'^[A-Z]+$')
    processed_sequence_count = 0
    print(f"正在读取输入文件：'{input_file_path}'...")
    
    try:
        with open(input_file_path, 'r', encoding='utf-8') as f_in:
            lines = f_in.readlines()
        for i in range(2, len(lines), 2): 
            if i + 1 < len(lines): 
                desc_line = lines[i].strip()
                sequence_line = lines[i+1].strip()

                if sequence_pattern.match(sequence_line): 
                    current_sequence = sequence_line
                    
                    processed_sequence_count += 1
                    output_filename = f"{base_input_filename}_{processed_sequence_count}.resfile"
                    output_file_path = os.path.join(output_folder, output_filename)

                    print(f"\n  正在为序列 '{current_sequence}' (来自描述 '{desc_line}') 生成文件：'{output_filename}'")

                    with open(output_file_path, 'w', encoding='utf-8') as f_out:
                        f_out.write("NATAA\n")
                        f_out.write("start\n")
                        for idx, amino_acid in enumerate(current_sequence):
                            f_out.write(f"{idx + 1} X PIKAA {amino_acid}\n")
                    
                    print(f"  文件 '{output_filename}' 生成成功。")
                else:
                    print(f"\n  警告: 跳过不符合序列模式的行 (行 {i+2}): '{sequence_line}'，或没有足够的行来形成描述-序列对。")
            else:
                print(f"\n  警告: 输入文件末尾存在不完整的描述-序列对，跳过行 {i+1} 及后续内容。")

    except Exception as e:
        print(f"处理输入文件 '{input_file_path}' 时发生严重错误: {e}")

    print(f"完成对文件 '{input_file_path}' 的处理。共成功生成 {processed_sequence_count} 个输出文件。")

if __name__ == "__main__":
    # 包含多个.fa文件的seqs文件夹。
    INPUT_MPNN_FILES_FOLDER = "/data1/home/renxinyu/mpnn/highmpnn/outputs/TNFa/to3rosetta/ProteinMPNN/seqs"
    # 用于存放生成的.resfile文件夹。
    OUTPUT_ROSETTA_RESFILES_FOLDER = "/data1/home/renxinyu/rosetta/outputs/to3rosetta/prroteinMPNN_resfiles"
    os.makedirs(OUTPUT_ROSETTA_RESFILES_FOLDER, exist_ok=True)
    # 2. 遍历输入文件夹中的所有文件
    input_files_list = glob.glob(os.path.join(INPUT_MPNN_FILES_FOLDER, '*.fa'))
    if not input_files_list:
        print(f"警告：在文件夹 '{INPUT_MPNN_FILES_FOLDER}' 中没有找到任何 .fa 文件。请检查路径和文件扩展名。")
    else:
        print(f"在 '{INPUT_MPNN_FILES_FOLDER}' 中找到 {len(input_files_list)} 个文件，准备开始批量处理...")
        for file_path in input_files_list:
            process_proteinmpnn_output_to_resfile(
                input_file_path=file_path,
                output_folder=OUTPUT_ROSETTA_RESFILES_FOLDER
            )
        print("\n所有文件批量处理完成。")