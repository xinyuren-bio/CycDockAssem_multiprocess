import os
import requests
import argparse
import re

def download_pdb(pdb_id_4char, output_dir, original_input_id=""):
    pdb_id_lower = pdb_id_4char.lower()
    pdb_id_upper = pdb_id_4char.upper()

    download_url = f"https://files.rcsb.org/download/{pdb_id_lower}.pdb"
    output_filename = f"{pdb_id_upper}.pdb"
    output_filepath = os.path.join(output_dir, output_filename)

    log_ref = f"(源自输入: {original_input_id})" if original_input_id else f"(PDB ID: {pdb_id_upper})"

    if os.path.exists(output_filepath):
        print(f"文件 {output_filepath} {log_ref} 已存在。跳过下载。")
        return

    try:
        response = requests.get(download_url, timeout=30)
        response.raise_for_status()

        with open(output_filepath, 'wb') as f:
            f.write(response.content)

    except Exception as e:
        print(f"处理PDB ID '{pdb_id_upper}' {log_ref} 时发生未知错误：{e}")

def main():
    parser = argparse.ArgumentParser(description="根据PDB ID列表从RCSB下载PDB文件。")
    parser.add_argument("input_file", help="包含PDB ID列表的输入文件路径。每行格式如 '12asA,4,5,...' 或仅 '1xyz'。")
    parser.add_argument("output_dir", help="用于存储下载的PDB文件的目录路径。")

    args = parser.parse_args()

    input_file_path = args.input_file
    output_directory_path = args.output_dir

    try:
        os.makedirs(output_directory_path, exist_ok=True)
        print(f"PDB文件将保存到: {os.path.abspath(output_directory_path)}")
    except OSError as e:
        print(f"错误：无法创建输出目录 '{output_directory_path}': {e}")
        return

    if not os.path.isfile(input_file_path):
        print(f"错误：找不到输入文件 '{input_file_path}'")
        return

    downloaded_pdb_files = set()

    print(f"\n正在从文件 '{input_file_path}' 读取PDB ID...")
    try:
        with open(input_file_path, 'r') as f:
            for line_number, line in enumerate(f, 1):
                line_content = line.strip()
                if not line_content or line_content.startswith('#'):
                    continue

                pdb_input_str = line_content.split(',')[0].strip()
                
                if pdb_input_str.isdigit():
                    continue
                
                pdb_id_4char = None
                chain_id = None
                is_valid_pdb_line = False

                if len(pdb_input_str) >= 4:
                    potential_pdb_id = pdb_input_str[:4]
                    if re.match(r"^[a-zA-Z0-9]{4}$", potential_pdb_id):
                        pdb_id_4char = potential_pdb_id
                        is_valid_pdb_line = True
                        if len(pdb_input_str) >= 5:
                            if re.match(r"^[a-zA-Z0-9]$", pdb_input_str[4]):
                                chain_id = pdb_input_str[4]
                            
                    else:
                        continue
                
                if not is_valid_pdb_line:
                    continue

                if pdb_id_4char:
                    pdb_id_to_download = pdb_id_4char.upper()
                    
                    if pdb_id_to_download not in downloaded_pdb_files:
                        download_pdb(pdb_id_to_download, output_directory_path, pdb_input_str)
                        downloaded_pdb_files.add(pdb_id_to_download)
                    else:
                        print(f"PDB文件 {pdb_id_to_download}.pdb (源自包含 '{pdb_input_str}' 的行) 之前已尝试下载/已存在。跳过。")

    except IOError as e:
        print(f"错误：无法读取输入文件 '{input_file_path}': {e}")
    except Exception as e:
        print(f"处理输入文件时发生未知错误: {e}")

    print("\n所有处理完成。")

if __name__ == "__main__":
    main()