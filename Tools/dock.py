import multiprocessing
import subprocess
import os
import argparse
from tqdm import tqdm

def run_sdock_command(args_tuple):
    """
    Executes a single SDOCK command and captures its output.
    """
    command_args, task_id = args_tuple
    log_prefix = f"Task {task_id}"

    try:
        process = subprocess.Popen(
            command_args,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding='utf-8'
        )
        stdout, stderr = process.communicate()

        if process.returncode == 0:
            return (task_id, True, "Execution successful")
        else:
            error_details = []
            if stdout and stdout.strip():
                error_details.append(f"--- SDOCK Standard Output ---\n{stdout.strip()}")
            if stderr and stderr.strip():
                error_details.append(f"--- SDOCK Error Output ---\n{stderr.strip()}")
            if not error_details:
                full_error_output = "Program crashed with no standard or error output."
            else:
                full_error_output = "\n\n".join(error_details)
            error_message = (f"Command failed with return code {process.returncode}.\n\n"
                             f"{full_error_output}")
            tqdm.write(f"\n🔴 {log_prefix}: {error_message}\n")
            return (task_id, False, error_message)

    except FileNotFoundError:
        error_message = f"Error: Command '{command_args[0]}' not found. Check path."
        tqdm.write(f"\n🔴 {log_prefix}: {error_message}\n")
        return (task_id, False, error_message)
    except Exception as e:
        error_message = f"An unknown script error occurred: {e}"
        tqdm.write(f"\n🔴 {log_prefix}: {error_message}\n")
        return (task_id, False, error_message)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parallel SDOCK task runner.")

    # Core parameters
    parser.add_argument("--num_tasks", type=int, required=True,
                        help="Total number of SDOCK tasks/files to process. Also used as number of parallel processes.")
    parser.add_argument("--sdock_exec", type=str, required=False, default="/data1/home/renxinyu/sdock/SDOCK2.0-restrict/sdock",
                        help="Path to the SDOCK executable (e.g., ../SDOCK2.0-restrict/sdock).")
    parser.add_argument("--fixed_pdb", type=str, required=True,
                        help="Path to the fixed PDB file (e.g., ALK1box1.pdb).")
    parser.add_argument("--fragstruct_dir", type=str, required=True,
                        help="Path to the directory containing fragstruct_X files.")
    parser.add_argument("--fixed_wat_pdb", type=str, required=True,
                        help="Path to the fixed water PDB file (e.g., ALK1box1wat.pdb).")
    parser.add_argument("--fragwatmap_dir", type=str, required=True,
                        help="Path to the directory containing fragwatmap_X files.")
    parser.add_argument("--output_record_dir", type=str, required=True,
                        help="Path to the directory where output record files (e.g., box1recordf_X) will be saved.")
    parser.add_argument("--output_prefix", type=str, required=True,
                        help="Prefix for the output record files (e.g., box1recordf).")
    parser.add_argument("--common_params", type=str, default="-c 0.20 -x 12 -B 1 -p 0 -r ../SDOCK2.0-restrict/so3layer_648.qua -n 10 -d 1.3",
                        help="Common SDOCK parameters as a single quoted string.")
    
    args = parser.parse_args()

    # Assign parsed arguments to variables
    TOTAL_TASKS_TO_RUN = args.num_tasks
    PROCESSES_TO_USE = args.num_tasks # num_tasks and processes are the same
    SDOCK_EXECUTABLE = args.sdock_exec
    FIXED_PDB_FILE = args.fixed_pdb
    FRAGSTRUCT_DIR = args.fragstruct_dir
    FIXED_WAT_PDB_FILE = args.fixed_wat_pdb
    FRAGWATMAP_DIR = args.fragwatmap_dir
    OUTPUT_RECORD_DIR = args.output_record_dir
    OUTPUT_FILE_PREFIX = args.output_prefix
    COMMON_PARAMS_LIST = args.common_params.split() # Split common params string into a list

    # Ensure output directory exists
    os.makedirs(OUTPUT_RECORD_DIR, exist_ok=True)

    print(f"⚙️ Generating {TOTAL_TASKS_TO_RUN} tasks...")
    all_tasks_params = []
    for task_index in range(1, TOTAL_TASKS_TO_RUN + 1):
        # Construct full paths for fragstruct, fragwatmap, and output files
        fragstruct_file = os.path.join(FRAGSTRUCT_DIR, f'fragstruct_{task_index}')
        fragwatmap_file = os.path.join(FRAGWATMAP_DIR, f'fragwmap_{task_index}')
        output_file = os.path.join(OUTPUT_RECORD_DIR, f'{OUTPUT_FILE_PREFIX}_{task_index}')

        # Construct the SDOCK command for the current task
        current_command_args = [
            SDOCK_EXECUTABLE,
            FIXED_PDB_FILE,
            fragstruct_file,
            FIXED_WAT_PDB_FILE,
            fragwatmap_file,
        ]
        print(current_command_args)
        current_command_args.extend(COMMON_PARAMS_LIST)
        current_command_args.extend(["-o", output_file])
        
        all_tasks_params.append((current_command_args, task_index))

    total_tasks_to_process = len(all_tasks_params)
    print(f"✅ Successfully generated {total_tasks_to_process} SDOCK tasks.")

    if total_tasks_to_process == 0:
        print("🤷 No tasks to execute, exiting program.")
        sys.exit(0)

    print(f"🚀 Preparing to launch {total_tasks_to_process} tasks, using {PROCESSES_TO_USE} parallel processes.")
    all_results = []
    # Using a batch size of 1 for tqdm progress bar to update per task completion
    # as pool.imap_unordered yields results as they complete.
    BATCH_SIZE_FOR_PROGESS = 1 
    with multiprocessing.Pool(processes=PROCESSES_TO_USE) as pool:
        with tqdm(total=total_tasks_to_process, desc="SDOCK Overall Progress") as pbar:
            # imap_unordered is good for showing progress as tasks complete
            for result in pool.imap_unordered(run_sdock_command, all_tasks_params):
                all_results.append(result)
                pbar.update(1)

    print("\n🏁 All SDOCK tasks have been executed.")

    # Summarize results
    all_results.sort(key=lambda r: r[0])
    successful_tasks = [r for r in all_results if r[1]]
    failed_tasks = [r for r in all_results if not r[1]]

    print(f"\n--- Execution Summary ---")
    print(f"🟢 Success: {len(successful_tasks)}")
    print(f"🔴 Failed: {len(failed_tasks)}")

    if failed_tasks:
        print("\n--- Failed Task Details ---")
        for task_id, _, message in failed_tasks:
            summary = message.splitlines()[0] if message else "No detailed information"
            print(f"  - Task {task_id}: {summary}")