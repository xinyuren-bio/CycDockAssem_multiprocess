import os
import subprocess
import argparse
import re

def run_assemble_cyc_command(assemble_executable, input_file_path, output_file_path, param1, param2, dockmodel_dir):
    command_list = [
        assemble_executable,
        input_file_path,
        output_file_path,
        str(param1),
        str(param2),
        dockmodel_dir
    ]

    try:
        result = subprocess.run(
            command_list,
            check=True,
            capture_output=True,
            text=True,
            encoding='utf-8'
        )
        print(f"✅ 成功处理文件: {os.path.basename(input_file_path)}. 输出到: {output_file_path}")
        if result.stdout:
            print(f"   stdout: {result.stdout.strip()}")
        if result.stderr:
            print(f"   stderr: {result.stderr.strip()}")
        return True
    except FileNotFoundError:
        print(f"🔴 错误: 找不到可执行文件 '{assemble_executable}'。请检查路径。")
        return False
    except subprocess.CalledProcessError as e:
        print(f"🔴 错误: 命令执行失败，处理文件 '{os.path.basename(input_file_path)}'。")
        print(f"   返回码: {e.returncode}")
        print(f"   stdout: {e.stdout.strip()}")
        print(f"   stderr: {e.stderr.strip()}")
        return False
    except Exception as e:
        print(f"🔴 错误: 处理文件 '{os.path.basename(input_file_path)}' 时发生未知错误: {e}")
        return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="批量执行 AssembleCyc 命令。")
    parser.add_argument("--assemble_exec", type=str, default="../utility/AssembleCyc",
                        help="AssembleCyc 可执行命令的完整路径 (例如: ../utility/AssembleCyc)。")
    parser.add_argument("--fraglinking_dir", type=str, required=True,
                        help="包含 fraglinking/ALK1fraglink_X 文件的输入目录路径。")
    parser.add_argument("--assembled_cyc_dir", type=str, required=True,
                        help="输出 AssembledCyc/ALK1Cyc_X.pdb 文件的目录路径。")
    parser.add_argument("--input_prefix", type=str, default="ALK1fraglink_",
                        help="输入文件名的前缀 (例如: ALK1fraglink_).")
    parser.add_argument("--input_suffix", type=str, default="",
                        help="输入文件名的后缀 (例如: '.txt' 或 '.pdb').")
    parser.add_argument("--output_prefix", type=str, default="Cyc_", # <--- 修改默认值为更通用的 "Cyc_"
                        help="输出文件名的前缀 (例如: ALK1Cyc_ 或 TestCyc_).") # <--- 更新帮助信息
    parser.add_argument("--param1", type=int, default=8,
                        help="第一个固定参数 (例如: 8)。")
    parser.add_argument("--param2", type=str, default="X",
                        help="第二个固定参数 (例如: X)。")
    parser.add_argument("--dockmodel_dir", type=str, required=True,
                        help="包含 SDOCK_Testbox1/2_*pdb 文件的 dockmodel 目录路径。")

    args = parser.parse_args()

    assemble_executable = os.path.expanduser(args.assemble_exec)
    fraglinking_dir = os.path.expanduser(args.fraglinking_dir)
    assembled_cyc_dir = os.path.expanduser(args.assembled_cyc_dir)
    input_prefix = args.input_prefix
    input_suffix = args.input_suffix
    output_prefix = args.output_prefix # <--- 现在这里会接收命令行传入的值
    param1 = args.param1
    param2 = args.param2
    dockmodel_dir = os.path.expanduser(args.dockmodel_dir)

    os.makedirs(assembled_cyc_dir, exist_ok=True)

    print(f"--- 开始批量执行 AssembleCyc 命令 ---")
    print(f"可执行文件: {assemble_executable}")
    print(f"输入目录: {fraglinking_dir}")
    print(f"输出目录: {assembled_cyc_dir}")
    print(f"输入文件前缀: {input_prefix}")
    print(f"输入文件后缀: {input_suffix}")
    print(f"输出文件前缀: {output_prefix}") # <--- 打印新的输出前缀
    print(f"固定参数: {param1}, {param2}")
    print(f"Dockmodel 目录: {dockmodel_dir}")
    print("-" * 30)

    file_pattern = re.compile(rf"^{re.escape(input_prefix)}(\d+){re.escape(input_suffix)}$")
    
    found_files_info = []
    try:
        for filename in os.listdir(fraglinking_dir):
            match = file_pattern.match(filename)
            if match and os.path.isfile(os.path.join(fraglinking_dir, filename)):
                file_number = int(match.group(1))
                found_files_info.append((file_number, filename))
        
        found_files_info.sort()
    except FileNotFoundError:
        print(f"🔴 错误: 输入目录 '{fraglinking_dir}' 未找到。")
    except Exception as e:
        print(f"🔴 错误: 读取输入目录 '{fraglinking_dir}' 时发生错误: {e}")
 
    if not found_files_info:
        print(f"🔴 警告: 在 '{fraglinking_dir}' 中未找到匹配 '{input_prefix}X{input_suffix}' 模式的文件。请检查输入目录和前缀/后缀。")

    total_files_to_process = len(found_files_info)
    print(f"将处理 {total_files_to_process} 个文件。")
    print("-" * 30)

    for i, (file_number, input_filename_base) in enumerate(found_files_info, 1):
        full_input_path = os.path.join(fraglinking_dir, input_filename_base)
        # 使用传入的 output_prefix 来构建输出文件名
        output_filename_base = f"{output_prefix}{file_number}.pdb" # <--- 使用动态的 output_prefix
        full_output_path = os.path.join(assembled_cyc_dir, output_filename_base)

        print(f"处理文件 {i}/{total_files_to_process}: {full_input_path}")
        run_assemble_cyc_command(assemble_executable, full_input_path, full_output_path, param1, param2, dockmodel_dir)
        print("-" * 20)

    print("--- 所有 AssembleCyc 任务执行完毕 ---")