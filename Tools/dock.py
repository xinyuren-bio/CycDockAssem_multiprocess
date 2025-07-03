import multiprocessing
import subprocess
import time

def run_sdock_command(command_args, task_id):
    full_command = " ".join(command_args) # 将参数列表组合成一个字符串命令
    print(f"任务 {task_id}: 正在执行命令: {full_command}")
    try:
        result = subprocess.run(command_args, capture_output=True, text=True, check=True)
        if result.stderr:
            print(f"任务 {task_id} 标准错误:\n{result.stderr}")
    except Exception as e:
        print(f"任务 {task_id}: 发生未知错误: {e}")
    print(f"任务 {task_id}: 命令执行完毕。")

if __name__ == "__main__":
    # 定义可执行文件路径
    sdock_executable ='../../../CycAssem-main/SDOCK2.0-restrict/sdock'
    # 选择需要对接的蛋白
    fixed_pdb_file = '../2_gene_boxs2proteinandW/box1.pdb'
    # 选择对接代码的watermap文件
    fixed_wat_pdb_file = '../2_gene_boxs2proteinandW/box1wat.pdb'
    # 定义对接参数
    common_params = [
        "-c", "0.20",
        "-x", "12", 
        "-B", "1",
        "-p", "0",
        "-r", '../../../CycAssem-main/SDOCK2.0-restrict/so3layer_648.qua',
        "-n", "10", 
        "-d", "1.3"
    ]

    cpu_configs = [
        # CPU 1
        {"fragstruct": "fragstruct00", "fragwatmap": "fragwatmap00", "output_file": "box1recordf00"},
        # CPU 2
        {"fragstruct": "fragstruct01", "fragwatmap": "fragwatmap01", "output_file": "box1recordf01"},
        # CPU 3
        {"fragstruct": "fragstruct02", "fragwatmap": "fragwatmap02", "output_file": "box1recordf02"},
        # CPU 4
        {"fragstruct": "fragstruct03", "fragwatmap": "fragwatmap03", "output_file": "box1recordf03"},
        # CPU 5
        {"fragstruct": "fragstruct04", "fragwatmap": "fragwatmap04", "output_file": "box1recordf04"},
        # CPU 6
        {"fragstruct": "fragstruct05", "fragwatmap": "fragwatmap05", "output_file": "box1recordf05"},
    ]

    processes = []
    for i, config in enumerate(cpu_configs):
        current_command_args = [
            sdock_executable,
            fixed_pdb_file,
            config["fragstruct"],   
            fixed_wat_pdb_file,     
            config["fragwatmap"]    
        ]
        
        current_command_args.extend(common_params)

        current_command_args.extend(["-o", config["output_file"]])
        process = multiprocessing.Process(target=run_sdock_command, args=(current_command_args, i + 1))
        processes.append(process)
        process.start()
    for process in processes:
        process.join()

    print("\n所有SDOCK任务都已完成。")