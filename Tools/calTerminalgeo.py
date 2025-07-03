import os
import csv
import subprocess
import time
import sys
from multiprocessing import Pool
import psutil

# ==============================================================================
# ===                           1. 配置区                                  ===
# ==============================================================================
# 1. C文件可执行路径
C_PROCESSOR_PATH = "/data1/home/renxinyu/sdock/fragtermgeo/calTerminalgeo"
# 2. 输入PDB文件库的路径
INPUT_PDB_DIR = "/data1/home/renxinyu/data/sdock/3_fragments"
# 3. 新的输出目录，用于存储处理结果
OUTPUT_DIR = "/data1/home/renxinyu/data/sdock/geo_csv"
# 4. [续跑源1] 用于检查进度的旧格式文件所在目录
LEGACY_PROGRESS_DIR = "/data1/home/renxinyu/data/sdock/geo_3"
# 5. 进程数目
MAX_PROCESSES = 12
# 6. 结果文件名
MASTER_D_FILE = "d_all_3.csv"
MASTER_D_UPPER_FILE = "D_all_3.csv" # << [续跑源2] 程序也会检查此文件
# 7. 日志文件路径
LOG_FILE = "processing_log.txt"


def setup_logging(log_file):
    class Logger:
        def __init__(self, filename):
            self.terminal = sys.stdout
            self.log = open(filename, "a", encoding='utf-8')
        def write(self, message):
            self.terminal.write(message)
            self.log.write(message)
        def flush(self):
            self.terminal.flush()
            self.log.flush()
    sys.stdout = Logger(log_file)
    sys.stderr = sys.stdout

def load_all_completed_pdbs(legacy_dir_path, csv_file_path):
    completed_pdbs = set()
    print(f"[INFO] [1/2] 正在扫描旧格式目录 '{legacy_dir_path}' ...")
    if not os.path.isdir(legacy_dir_path):
        print(f"[WARN] 找不到指定的旧进度目录 '{legacy_dir_path}'。")
    else:
        for root, _, files in os.walk(legacy_dir_path):
            for filename in files:
                if filename.startswith('D'):
                    file_path = os.path.join(root, filename)
                    try:
                        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                            for line in f:
                                parts = line.strip().split()
                                if parts:
                                    completed_pdbs.add(parts[0])
                    except Exception as e:
                        print(f"\n[WARN] 读取旧进度文件 {file_path} 时出错: {e}")

    print(f"[INFO] 从旧格式目录中发现 {len(completed_pdbs)} 个已完成的PDB。")

    print(f"[INFO] [2/2] 正在扫描主CSV文件 '{csv_file_path}' ...")
    if not os.path.exists(csv_file_path):
        print(f"[INFO] 未找到主CSV文件，跳过检查。")
    else:
        try:
            with open(csv_file_path, 'r', encoding='utf-8', newline='') as f:
                reader = csv.reader(f)
                header = next(reader, None)
                if header and header[0] == 'pdb_filename':
                    for row in reader:
                        if row:
                            completed_pdbs.add(row[0])
                else:
                    print(f"[WARN] 主CSV文件 {csv_file_path} 格式不正确或为空。")
        except Exception as e:
            print(f"[ERROR] 读取主CSV文件 {csv_file_path} 失败: {e}")

    print(f"[INFO] 综合所有来源，总计发现 {len(completed_pdbs)} 个已完成的PDB。")
    return completed_pdbs


# --- 核心工作函数 (由子进程调用) ---
def process_single_pdb_file(file_info):
    relative_path, input_path = file_info
    filename_only = os.path.basename(input_path)
    len_num = 3

    command = [C_PROCESSOR_PATH, input_path, str(len_num)]

    try:
        process = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore'
        )
        cmd_output = process.stdout
    except subprocess.CalledProcessError as e:
        error_msg = f"C脚本执行失败 (退出码 {e.returncode}) on {filename_only}: {e.stderr.strip()}"
        return filename_only, False, error_msg
    except Exception as e:
        error_msg = f"未知错误 on {filename_only}: {e}"
        return filename_only, False, error_msg

    lines = cmd_output.strip().split('\n')
    if len(lines) < 3:
        error_msg = f"C脚本输出行数不足 ({len(lines)} lines) for {filename_only}."
        return filename_only, False, error_msg

    try:
        dist = float(lines[0].split()[-1])
        parts2 = lines[1].split()
        val_d_1, val_d_2, val_d_3 = float(parts2[2]), float(parts2[3]), float(parts2[4])
        parts3 = lines[2].split()
        val_D_1, val_D_2, val_D_3 = float(parts3[2]), float(parts3[3]), float(parts3[4])
    except (ValueError, IndexError) as e:
        error_msg = f"解析C脚本输出失败 for {filename_only}: {e}. 输出:\n{cmd_output.strip()}"
        return filename_only, False, error_msg

    results = [
        {'type': 'd', 'data': (filename_only, dist, val_d_1, val_d_2, val_d_3)},
        {'type': 'D', 'data': (filename_only, dist, val_D_1, val_D_2, val_D_3)}
    ]

    return filename_only, True, results


def main():
    setup_logging(LOG_FILE)
    try:
        original_stdout = sys.stdout.terminal
    except AttributeError:
        original_stdout = sys.stdout

    print(f"--- 脚本启动 ({time.strftime('%Y-%m-%d %H:%M:%S')}) ---")
    print(f"C 可执行文件: {os.path.abspath(C_PROCESSOR_PATH)}")
    print(f"输入目录: {os.path.abspath(INPUT_PDB_DIR)}")
    print(f"新输出目录: {os.path.abspath(OUTPUT_DIR)}")
    print(f"检查进度的旧目录: {os.path.abspath(LEGACY_PROGRESS_DIR)}")
    print(f"并行进程数: {MAX_PROCESSES}")
    print("-" * 60)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    master_d_path = os.path.join(OUTPUT_DIR, MASTER_D_FILE)
    master_D_path = os.path.join(OUTPUT_DIR, MASTER_D_UPPER_FILE)

    # 1. [续跑] MODIFIED: 调用综合函数，从两个来源加载进度
    completed_pdb_names = load_all_completed_pdbs(LEGACY_PROGRESS_DIR, master_D_path)

    # 2. 扫描输入目录，过滤掉已完成的
    files_to_process = []
    print("正在扫描输入PDB目录以收集待处理文件...", flush=True)
    for root, _, files in os.walk(INPUT_PDB_DIR):
        for file_name in files:
            if not (file_name.startswith(("tp_", "Tp_", "Pp_"))) or not file_name.endswith(".pdb"):
                continue
            
            if file_name not in completed_pdb_names:
                input_full_path = os.path.join(root, file_name)
                relative_path = os.path.relpath(input_full_path, INPUT_PDB_DIR)
                files_to_process.append((relative_path, input_full_path))

    total_to_process_count = len(files_to_process)
    if total_to_process_count == 0:
        print("\n所有符合条件的文件均已处理完成，无需任何操作。", flush=True)
        return

    print(f"扫描完成。发现需处理的新文件数: {total_to_process_count}", flush=True)
    
    processed_count = 0
    global_start_time = time.time()
    csv_header = ['pdb_filename', 'distance', 'value1', 'value2', 'value3']
    
    try:
        with open(master_d_path, 'a', encoding='utf-8', newline='') as f_d, \
             open(master_D_path, 'a', encoding='utf-8', newline='') as f_D:
            
            writer_d = csv.writer(f_d)
            writer_D = csv.writer(f_D)
            
            if f_d.tell() == 0: writer_d.writerow(csv_header)
            if f_D.tell() == 0: writer_D.writerow(csv_header)
            
            with Pool(processes=MAX_PROCESSES) as pool:
                for pdb_name, success, result_data in pool.imap_unordered(process_single_pdb_file, files_to_process):
                    processed_count += 1
                    
                    if success:
                        for result in result_data:
                            if result['type'] == 'd':
                                writer_d.writerow(result['data'])
                            elif result['type'] == 'D':
                                writer_D.writerow(result['data'])
                        f_d.flush()
                        f_D.flush()
                    else:
                        print(f"\n[FAIL] {result_data}", flush=True)

                    elapsed = time.time() - global_start_time
                    progress_percent = (processed_count / total_to_process_count) * 100
                    eta_str = "N/A"
                    if progress_percent > 1:
                        eta_seconds = (elapsed / progress_percent) * (100 - progress_percent)
                        eta_str = time.strftime('%H:%M:%S', time.gmtime(eta_seconds))
                    
                    progress_line = f"\r进度: {processed_count}/{total_to_process_count} ({progress_percent:.2f}%) | 耗时: {time.strftime('%H:%M:%S', time.gmtime(elapsed))} | 预计剩余: {eta_str}"
                    original_stdout.write(progress_line)
                    original_stdout.flush()

    except KeyboardInterrupt:
        print("\n\n[INFO] 检测到中断 (Ctrl+C)。程序将安全退出。", flush=True)
    except Exception as e:
        print(f"\n\n[FATAL] 发生严重错误: {e}", flush=True)
    finally:
        original_stdout.write("\n")
        original_stdout.flush()
        
        print(f"\n--- 脚本结束 ({time.strftime('%Y-%m-%d %H:%M:%S')}) ---")
        total_time = time.time() - global_start_time
        print(f"本次运行处理文件数: {processed_count}", flush=True)
        print(f"总计耗时: {time.strftime('%H:%M:%S', time.gmtime(total_time))}", flush=True)

if __name__ == "__main__":
    main()