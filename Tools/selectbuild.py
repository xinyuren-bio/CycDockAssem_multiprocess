import os
import subprocess
import argparse
import re

def run_select_build_complex_command(select_build_executable, input_file_path, param1_score, param2_int, param3_score, param4_int, complex_model_dir, fixed_pdb_file, prefix):
    """
    执行 SelectBuildcomplex 命令的函数。
    """
    command_list = [
        select_build_executable,
        input_file_path,
        str(param1_score),
        str(param2_int),
        str(param3_score),
        str(param4_int),
        complex_model_dir, # 注意这里是目录
        fixed_pdb_file,
        prefix
    ]

    try:
        # 执行命令并捕获输出
        result = subprocess.run(
            command_list,
            check=True, # 如果返回非零退出码，则抛出 CalledProcessError
            capture_output=True,
            text=True,
            encoding='utf-8'
        )
        print(f"✅ Successfully processed file: {os.path.basename(input_file_path)}. Output to: {complex_model_dir}")
        if result.stdout:
            print(f"   stdout: {result.stdout.strip()}")
        if result.stderr:
            print(f"   stderr: {result.stderr.strip()}")
        return True
    except FileNotFoundError:
        print(f"🔴 Error: Executable '{select_build_executable}' not found. Please check the path.")
        return False
    except subprocess.CalledProcessError as e:
        print(f"🔴 Error: Command failed for file '{os.path.basename(input_file_path)}'.")
        print(f"   Return code: {e.returncode}")
        print(f"   stdout: {e.stdout.strip()}")
        print(f"   stderr: {e.stderr.strip()}")
        return False
    except Exception as e:
        print(f"🔴 Error: An unknown error occurred while processing file '{os.path.basename(input_file_path)}': {e}")
        return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Batch execution of SelectBuildcomplex command.")
    parser.add_argument("--select_build_exec", type=str, default="../utility/SelectBuildcomplex",
                        help="Full path to the SelectBuildcomplex executable (e.g., ../utility/SelectBuildcomplex).")
    parser.add_argument("--assembled_cyc_dir", type=str, required=True,
                        help="Input directory containing files like ALK1Cyc_X.pdb.")
    parser.add_argument("--complex_model_dir", type=str, required=True,
                        help="Output directory for complexmodel files (e.g., ../ALK1/complexmodel/).")
    parser.add_argument("--fixed_pdb_file", type=str, required=True,
                        help="Path to the fixed PDB file (e.g., ../ALK1/box1.pdb).")
    parser.add_argument("--param1_score", type=float, default=-55.0,
                        help="First score parameter (e.g., -55.0).")
    parser.add_argument("--param2_int", type=int, default=2,
                        help="First integer parameter (e.g., 2).")
    parser.add_argument("--param3_score", type=float, default=-55.0,
                        help="Second score parameter (e.g., -55.0).")
    parser.add_argument("--param4_int", type=int, default=3,
                        help="Second integer parameter (e.g., 3).")
    parser.add_argument("--input_prefix", type=str, default="ALK1Cyc_",
                        help="Prefix for input filenames (e.g., ALK1Cyc).")
    parser.add_argument("--input_suffix", type=str, default=".pdb",
                        help="Suffix for input filenames (e.g., .pdb).")

    args = parser.parse_args()

    select_build_executable = os.path.expanduser(args.select_build_exec)
    assembled_cyc_dir = os.path.expanduser(args.assembled_cyc_dir)
    complex_model_dir = os.path.expanduser(args.complex_model_dir)
    fixed_pdb_file = os.path.expanduser(args.fixed_pdb_file)
    param1_score = args.param1_score
    param2_int = args.param2_int
    param3_score = args.param3_score
    param4_int = args.param4_int
    input_prefix = args.input_prefix
    input_suffix = args.input_suffix

    # Ensure output directory exists
    os.makedirs(complex_model_dir, exist_ok=True)

    print(f"--- Starting batch execution of SelectBuildcomplex command ---")
    print(f"Executable: {select_build_executable}")
    print(f"Input directory: {assembled_cyc_dir}")
    print(f"Output directory: {complex_model_dir}")
    print(f"Fixed PDB file: {fixed_pdb_file}")
    print(f"Parameters: {param1_score}, {param2_int}, {param3_score}, {param4_int}")
    print("-" * 30)

    # Regex to match files like "ALK1Cyc_X.pdb"
    file_pattern = re.compile(rf"^{re.escape(input_prefix)}(\d+){re.escape(input_suffix)}$")
    
    found_files_info = []
    try:
        for filename in os.listdir(assembled_cyc_dir):
            match = file_pattern.match(filename)
            if match and os.path.isfile(os.path.join(assembled_cyc_dir, filename)):
                file_number = int(match.group(1))
                found_files_info.append((file_number, filename))
        
        found_files_info.sort() # Sort by file number to process in order
    except FileNotFoundError:
        print(f"🔴 Error: Input directory '{assembled_cyc_dir}' not found.")
    except Exception as e:
        print(f"🔴 Error: An error occurred while reading input directory '{assembled_cyc_dir}': {e}")

    if not found_files_info:
        print(f"🔴 Warning: No files matching '{input_prefix}X{input_suffix}' pattern found in '{assembled_cyc_dir}'. Please check input directory and prefix/suffix.")

    total_files_to_process = len(found_files_info)
    print(f"Will process {total_files_to_process} files.")
    print("-" * 30)

    for i, (file_number, input_filename_base) in enumerate(found_files_info, 1):
        full_input_path = os.path.join(assembled_cyc_dir, input_filename_base)

        print(f"Processing file {i}/{total_files_to_process}: {full_input_path}")
        run_select_build_complex_command(
            select_build_executable,
            full_input_path,
            param1_score,
            param2_int,
            param3_score,
            param4_int,
            complex_model_dir,
            fixed_pdb_file,
            input_prefix
        )
        print("-" * 20)

    print("--- All SelectBuildcomplex tasks completed ---")
