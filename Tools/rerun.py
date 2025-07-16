import os
import argparse
import logging
import sys
import re # re 模块在这里实际未使用，可以移除，但保留无害

def setup_logging(log_dir):
    """
    配置日志记录器，将日志同时输出到文件和控制台。
    """
    log_file = os.path.join(log_dir, 'resume_prep_by_index.log')
    # 确保日志目录存在
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    
    # 移除可能存在的旧的日志处理器，避免重复输出
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
        
    logging.basicConfig(
        level=logging.INFO, # 默认日志级别为 INFO，可以通过修改此行查看更详细的 DEBUG 信息
        format='%(asctime)s [%(levelname)s] - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=[
            logging.FileHandler(log_file, mode='w', encoding='utf-8'), # 将日志写入文件 (覆盖模式 'w')
            logging.StreamHandler(sys.stdout) # 将日志输出到控制台
        ]
    )

def get_completed_task_ids_from_dockresult(dockresult_dir, original_recordf_path, num_tasks):
    """
    扫描指定的 dockresult 目录，并根据原始 box1recordf 文件中的预期输出路径，
    来判断哪些任务已经完成。
    
    参数:
        dockresult_dir (str): 存放最终对接结果文件的目录。
        original_recordf_path (str): 包含所有原始任务输出路径的 master 文件。
        num_tasks (int): 原始 master 文件中的总任务数量。
        
    返回:
        set: 包含所有已完成任务的全局 ID 集合。
    """
    completed_ids = set()
    logging.info(f"🔍 扫描目录: {dockresult_dir} 和原始文件: {original_recordf_path} 查找已完成的任务...")

    # 构建预期文件名到任务ID的映射，用于快速查找
    expected_output_filenames_to_id = {}
    try:
        with open(original_recordf_path, 'r', encoding='utf-8') as f:
            for task_id in range(1, num_tasks + 1):
                line = f.readline().strip()
                if not line:
                    logging.warning(f"⚠️ 原始 {original_recordf_path} 文件在行 {task_id} 处提前结束，总任务数可能不符。")
                    break
                # 从路径中提取文件名 (例如: /path/to/recordf_output_123 -> recordf_output_123)
                filename = os.path.basename(line)
                expected_output_filenames_to_id[filename] = task_id
    except FileNotFoundError:
        logging.error(f"❌ 错误：原始 recordf 文件未找到: {original_recordf_path}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"❌ 读取原始 recordf 文件时发生错误: {e}")
        sys.exit(1)

    # 遍历 dockresult 目录下的所有文件，检查哪些是已完成任务
    for filename in os.listdir(dockresult_dir):
        # 只有当文件在我们的预期文件名列表中时才处理
        if filename in expected_output_filenames_to_id:
            full_path = os.path.join(dockresult_dir, filename)
            # 确认文件确实存在且其大小大于 0 (非空文件)
            if os.path.exists(full_path) and os.path.getsize(full_path) > 0:
                task_id = expected_output_filenames_to_id[filename]
                completed_ids.add(task_id)
            # else:
            #    logging.debug(f"DEBUG: 文件 {full_path} 存在但为空或不存在，不计入完成任务。")
        # else:
        #    logging.debug(f"DEBUG: 文件 {filename} 在 {dockresult_dir} 中存在，但不在 '{os.path.basename(original_recordf_path)}' 的预期输出列表中。")

    logging.info(f"✅ 在 {dockresult_dir} 找到 {len(completed_ids)} 个已完成任务。")
    return completed_ids

def prepare_for_resume_by_index():
    """
    根据已完成的对接结果 (dockresult)，生成未完成任务的新的总输入文件，用于续跑。
    """
    parser = argparse.ArgumentParser(
        description="根据已完成的对接结果 (dockresult)，生成未完成任务的新的总输入文件。"
    )
    parser.add_argument("-d", "--docking_dir", type=str, 
                        default='/data1/home/renxinyu/sdock/Outputs/fcrn/dock',
                        help="对接文件所在的根目录。")
    parser.add_argument("--output_suffix", type=str, default="_resume_all",
                        help="新生成的未完成总文件名的后缀 (默认: _resume_all)。")
    
    args = parser.parse_args()

    # --- 配置所有相关文件和目录的路径 ---
    docking_base_path = args.docking_dir
    original_recordf_path = os.path.join(docking_base_path, 'box1recordf_1000W') # 原始的总输出记录文件
    original_fragstruct_path = os.path.join(docking_base_path, 'fragstruct_1000W') # 原始的总 fragstruct 文件
    original_fragwatmap_path = os.path.join(docking_base_path, 'fragwatmap_1000W') # 原始的总 fragwatmap 文件
    
    dockresult_dir = os.path.join(docking_base_path, 'dockresult') # 最终对接结果的存放目录

    # 新生成续跑文件的目录和名称
    output_dir = os.path.join(docking_base_path, 'resume_prepared_files1')
    os.makedirs(output_dir, exist_ok=True) # 确保输出目录存在
    
    new_recordf_path = os.path.join(output_dir, f'box1recordf{args.output_suffix}')
    new_fragstruct_path = os.path.join(output_dir, f'fragstruct{args.output_suffix}')
    new_fragwatmap_path = os.path.join(output_dir, f'fragwatmap{args.output_suffix}')

    # --- 设置日志 ---
    # 日志文件将保存在 'record' 子目录中
    setup_logging(os.path.join(docking_base_path, 'record'))

    logging.info("--- 🚀 开始准备续跑文件 (基于 dockresult 产物) ---")
    logging.info(f"📂 对接根目录: {docking_base_path}")
    logging.info(f"📂 最终结果目录: {dockresult_dir}")
    logging.info(f"🆕 新总文件输出目录: {output_dir}")
    logging.info(f"原始主文件: {original_recordf_path}")

    # 自动读取原始 recordf 文件的总行数作为 num_tasks
    num_tasks = 0
    try:
        with open(original_recordf_path, 'r', encoding='utf-8') as f:
            num_tasks = sum(1 for line in f)
        logging.info(f"📚 从原始 '{original_recordf_path}' 文件中自动检测到总任务数: {num_tasks}")
    except FileNotFoundError:
        logging.error(f"❌ 错误：原始 '{original_recordf_path}' 文件未找到。无法自动确定总任务数。请确保文件存在。", exc_info=True)
        sys.exit(1)
    except Exception as e:
        logging.error(f"❌ 读取原始 '{original_recordf_path}' 文件时发生错误，无法自动确定总任务数: {e}", exc_info=True)
        sys.exit(1)

    if num_tasks == 0:
        logging.warning(f"⚠️ 原始文件 '{original_recordf_path}' 为空或不包含任何任务。无需处理。")
        sys.exit(0)

    # 1. 获取所有已完成的任务ID (通过检查 dockresult 目录)
    completed_global_indices = get_completed_task_ids_from_dockresult(
        dockresult_dir, original_recordf_path, num_tasks
    )
    
    # 2. 读取原始总文件，筛选未完成任务的行
    uncompleted_recordf_lines = []
    uncompleted_fragstruct_lines = []
    uncompleted_fragwatmap_lines = []
    
    try:
        with open(original_recordf_path, 'r', encoding='utf-8') as f_recordf, \
             open(original_fragstruct_path, 'r', encoding='utf-8') as f_fragstruct, \
             open(original_fragwatmap_path, 'r', encoding='utf-8') as f_fragwatmap:
            
            # 遍历从 1 到 num_tasks 的所有全局任务ID
            for global_idx in range(1, num_tasks + 1):
                recordf_line = f_recordf.readline().strip()
                fragstruct_line = f_fragstruct.readline().strip()
                fragwatmap_line = f_fragwatmap.readline().strip()

                # 检查所有文件是否同步结束
                if not all([recordf_line, fragstruct_line, fragwatmap_line]):
                    logging.warning(f"⚠️ 原始输入文件在行 {global_idx} 处提前结束。可能文件行数不一致或提前到达文件末尾。")
                    break # 如果任何一个文件提前结束，则停止读取

                # 如果当前全局ID不在已完成任务的集合中，则将其行加入未完成列表
                if global_idx not in completed_global_indices:
                    uncompleted_recordf_lines.append(recordf_line)
                    uncompleted_fragstruct_lines.append(fragstruct_line)
                    uncompleted_fragwatmap_lines.append(fragwatmap_line)
                
    except FileNotFoundError as e:
        logging.error(f"❌ 错误：找不到原始总输入文件之一: {e}. 请检查路径和文件名。", exc_info=True)
        sys.exit(1)
    except Exception as e:
        logging.error(f"❌ 读取原始总文件时发生未知错误: {e}", exc_info=True)
        sys.exit(1)

    num_uncompleted = len(uncompleted_recordf_lines)
    logging.info(f"📊 识别到 {num_uncompleted} 个未完成任务需要重新运行。")

    if num_uncompleted == 0:
        logging.info("🎉 所有任务均已完成，无需生成新文件。")
        # 如果所有任务都已完成，尝试清理生成的空目录
        if os.path.exists(output_dir) and not os.listdir(output_dir): # 只有当目录存在且为空时才删除
            try:
                os.rmdir(output_dir)
                logging.info(f"🗑️ 已清理空目录: {output_dir}")
            except OSError as e:
                logging.warning(f"⚠️ 无法删除空目录 {output_dir}: {e}. 可能非空或权限问题。")
        sys.exit(0)

    # 3. 写入新的总未完成任务文件
    try:
        with open(new_recordf_path, 'w', encoding='utf-8') as f_recordf, \
             open(new_fragstruct_path, 'w', encoding='utf-8') as f_fragstruct, \
             open(new_fragwatmap_path, 'w', encoding='utf-8') as f_fragwatmap:
            
            for i in range(num_uncompleted):
                f_recordf.write(uncompleted_recordf_lines[i] + '\n')
                f_fragstruct.write(uncompleted_fragstruct_lines[i] + '\n')
                f_fragwatmap.write(uncompleted_fragwatmap_lines[i] + '\n')
        
        logging.info(f"✅ 成功生成新的未完成任务总文件：")
        logging.info(f"   - {new_recordf_path}")
        logging.info(f"   - {new_fragstruct_path}")
        logging.info(f"   - {new_fragwatmap_path}")
        logging.info("--- ✨ 续跑文件准备完成 ---")

    except Exception as e:
        logging.error(f"❌ 写入新文件时发生错误: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    prepare_for_resume_by_index()