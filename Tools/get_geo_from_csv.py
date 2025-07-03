import os
import csv
import sys
import time
from collections import defaultdict

# ==============================================================================
# ===                           1. 配置区                                  ===
# ==============================================================================

# 保存geo数据的CSV文件目录
INPUT_CSV_DIR = "/data1/home/renxinyu/data/geo_csv" 
# 输入的CSV文件名和对应的前缀
CSV_FILES_TO_PROCESS = {
    "d": "d_all_3.csv",
    "D": "D_all_3.csv"
}
# 输出GeoDict格式的目录
GEODICT_OUTPUT_DIR = "/data1/home/renxinyu/data/geo_3"
# 日志文件名
LOG_FILE = "conversion_log.txt"

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

def format_line_content(pdb_filename, distance, val1, val2, val3):
    formatted_pdb_name = f"{os.path.basename(pdb_filename):<17.17}"

    if distance < 10.0:
        spaces_and_val1_field = f"     {distance:5.3f}"
    else:
        spaces_and_val1_field = f"    {distance:6.3f}"
    
    formatted_values = f"{val1:9.3f}{val2:9.3f}{val3:9.3f}"

    return f"{formatted_pdb_name}{spaces_and_val1_field}{formatted_values}\n"

def main():
    setup_logging(LOG_FILE)
    
    print(f"--- 开始运行 ({time.strftime('%Y-%m-%d %H:%M:%S')}) ---")
    print(f"输入: {os.path.abspath(INPUT_CSV_DIR)}")
    print(f"输出: {os.path.abspath(GEODICT_OUTPUT_DIR)}")
    print("-" * 60)

    os.makedirs(GEODICT_OUTPUT_DIR, exist_ok=True)

    global_start_time = time.time()
    total_lines_processed = 0
    
    open_file_handles = {}

    try:
        for file_prefix, csv_filename in CSV_FILES_TO_PROCESS.items():
            input_csv_path = os.path.join(INPUT_CSV_DIR, csv_filename)

            if not os.path.exists(input_csv_path):
                print(f"[WARN] Input file not found, skipping: {input_csv_path}")
                continue

            print(f"\nProcessing file: '{input_csv_path}'...")
            
            with open(input_csv_path, 'r', encoding='utf-8', newline='') as f_in:
                reader = csv.reader(f_in)
                
                # Validate and skip the header
                header = next(reader, None)
                if not header or header[0] != 'pdb_filename':
                    print(f"[ERROR] Invalid or empty CSV file: {input_csv_path}. Skipping.")
                    continue

                # Process each data row in the CSV
                for i, row in enumerate(reader):
                    try:
                        # Unpack and convert data types from the CSV row
                        pdb_name, dist_str, v1_str, v2_str, v3_str = row
                        distance = float(dist_str)
                        val1 = float(v1_str)
                        val2 = float(v2_str)
                        val3 = float(v3_str)

                        output_filename = f"{file_prefix}{distance:.3f}"
                        output_filepath = os.path.join(GEODICT_OUTPUT_DIR, output_filename)

                        content_line = format_line_content(pdb_name, distance, val1, val2, val3)

                        if output_filepath not in open_file_handles:
                            open_file_handles[output_filepath] = open(output_filepath, 'a', encoding='utf-8')
                        
                        handle = open_file_handles[output_filepath]
                        handle.write(content_line)

                        total_lines_processed += 1
                        if (i + 1) % 25000 == 0:
                            print(f"\r  ...processed {i + 1} lines.", end="", flush=True)

                    except (ValueError, IndexError) as e:
                        print(f"\n[WARN] Skipping malformed row {i+2} in {csv_filename}: {row}. Reason: {e}")
                        continue
            
            print(f"\r  ...processed all lines in {csv_filename}. Done.     ")

    except Exception as e:
        print(f"\n[FATAL] A critical error occurred: {e}")
    finally:
        print("\nClosing all file handles...")
        closed_count = 0
        for handle in open_file_handles.values():
            handle.close()
            closed_count += 1
        
        print(f"Successfully closed {closed_count} files.")
        
        total_time = time.time() - global_start_time
        print(f"\n--- Script finished ({time.strftime('%Y-%m-%d %H:%M:%S')}) ---")
        print(f"Total lines processed and written: {total_lines_processed}")
        print(f"Total execution time: {time.strftime('%H:%M:%S', time.gmtime(total_time))}")

if __name__ == "__main__":
    main()