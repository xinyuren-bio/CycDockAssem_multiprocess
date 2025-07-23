#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import math
import numpy as np
from tqdm import tqdm
import itertools
import concurrent.futures
import numpy as np
from scipy.spatial import cKDTree
import os

# --- Step 1 & 2: 导入库和定义常量 (Imports and Constants) ---
MAXCADIS = 13.5 * 13.5
MINCADIS = 4.7 * 4.7
AMMINCADIS = 3.4 * 3.4
TERMINALDISCAMAX3 = 7.5 * 7.5
TERMINALDISCAMAX4 = 10.5 * 10.5
MAXFRAGNUM = 300000

# --- Step 3: 创建 Protein 类 (Protein Class) ---
class Protein:
    def __init__(self):
        self.an = 0
        self.xyz = np.empty((0, 3), dtype=np.float32)
        self.name = ""
        self.score = 0.0
        # Initialize with -1 to indicate unset or invalid indices
        self.N, self.NC, self.C, self.CN, self.CAN, self.CAC = -1, -1, -1, -1, -1, -1

# --- Step 4: 转换几何计算函数 (Geometric Calculation Functions) ---
def getdx(x1, x2): return x1 - x2
def get_len(x): return np.linalg.norm(x)
def distance(x1, x2): return np.linalg.norm(x1 - x2)
def distance2(x1, x2): return np.sum((x1 - x2)**2)
def iprod(a, b): return np.dot(a, b)
def oprod(a, b): return np.cross(a, b)

def cos_angle(a, b):
    ip = np.dot(a, b)
    norm_product = np.linalg.norm(a) * np.linalg.norm(b)
    if norm_product == 0: return 0.0
    cos_val = ip / norm_product
    return np.clip(cos_val, -1.0, 1.0)

def cal_angle(x1, x2, x3):
    v12, v32 = getdx(x1, x2), getdx(x3, x2)
    # Ensure vectors are not zero length before calculating angle
    if np.linalg.norm(v12) == 0 or np.linalg.norm(v32) == 0:
        return 0.0 # Or raise an error, depending on desired behavior
    return math.degrees(math.acos(cos_angle(v12, v32)))

def cal_dih(x1, x2, x3, x4):
    x12, x32, x34 = getdx(x1, x2), getdx(x3, x2), getdx(x3, x4)
    m, n = oprod(x12, x32), oprod(x32, x34)
    # Ensure normal vectors are not zero length
    if np.linalg.norm(m) == 0 or np.linalg.norm(n) == 0:
        return 0.0 # Or handle as an error
    phi_rad = math.acos(cos_angle(m, n))
    sign = -1.0 if iprod(x12, n) < 0.0 else 1.0
    return math.degrees(sign * phi_rad)

# --- Step 5: 转换文件读取函数 (File Reading Function) ---
def getxyz(line):
    try:
        return [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    except (ValueError, IndexError):
        return None

def read_structure(pn):
    try:
        with open(pn, "r") as pf:
            p = Protein()
            p.name = os.path.basename(pn).strip() # Only store filename, not full path
            xyz_coords, atom_count, res_num = [], 0, 1
            for line in pf:
                if line.startswith("ATOM"):
                    coords = getxyz(line)
                    if coords is None: continue
                    xyz_coords.append(coords)
                    atom_type = line[12:16].strip()
                    atom_idx = atom_count # Store the current atom index

                    # Assign indices for critical atoms
                    if atom_type == "N":
                        if res_num == 1: p.N = atom_idx # First N
                        p.CN = atom_idx # C-terminus N (last N encountered for the peptide)
                    elif atom_type == "CA":
                        if res_num == 1: p.CAN = atom_idx # First CA
                        p.CAC = atom_idx # C-terminus CA (last CA encountered)
                    elif atom_type == "C":
                        if res_num == 1:
                            p.NC = atom_idx # N-terminus C (first C for the peptide)
                            res_num += 1 # Move to next residue after finding first C
                        p.C = atom_idx # C-terminus C (last C encountered)
                    atom_count += 1
                elif line.startswith("REMARK    score"):
                    try: p.score = float(line[15:])
                    except (ValueError, IndexError): p.score = 0.0
            p.an, p.xyz = atom_count, np.array(xyz_coords, dtype=np.float32)

            # Basic validation for critical indices after reading all atoms
            # If a PDB is very short or malformed, these might not be set.
            # We already initialized them to -1 in __init__.
            # It's better to ensure they point to valid atom positions if used.
            # For simplicity, we just rely on the -1 check in calculations.

            return p
    except FileNotFoundError:
        print(f"ERROR: Cannot open structure file {pn}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An error occurred while reading {pn}: {e}", file=sys.stderr)
        return None

def checkcollision(p1, p2):
    CLASH_DISTANCE_SQ = 2.5 * 2.5
    # Skip collision check if either protein has no atoms or critical atoms
    if p1.an == 0 or p2.an == 0:
        return True # Assuming no collision if no atoms (or handle as an error if empty files are truly invalid)

    tree1 = cKDTree(p1.xyz)
    tree2 = cKDTree(p2.xyz)

    # Find all pairs within CLASH_DISTANCE_SQ
    colliding_pairs_indices = tree1.query_ball_tree(tree2, math.sqrt(CLASH_DISTANCE_SQ))

    for i, j_list in enumerate(colliding_pairs_indices):
        if not j_list: # No clashes for this atom from p1
            continue
        for j in j_list:
            # Re-check for precise distance (query_ball_tree gives a sphere, not exact clash)
            if distance2(p1.xyz[i], p2.xyz[j]) < CLASH_DISTANCE_SQ:
                # Exclude specific atom types from collision check if they are terminal and connecting
                # Ensure the indices are valid before checking
                is_p1_c_or_cac = (i == p1.C and p1.C != -1) or (i == p1.CAC and p1.CAC != -1)
                is_p2_n_or_can = (j == p2.N and p2.N != -1) or (j == p2.CAN and p2.CAN != -1)

                is_p2_c_or_cac = (j == p2.C and p2.C != -1) or (j == p2.CAC and p2.CAC != -1)
                is_p1_n_or_can = (i == p1.N and p1.N != -1) or (i == p1.CAN and p1.CAN != -1)

                if (is_p1_c_or_cac and is_p2_n_or_can) or \
                   (is_p2_c_or_cac and is_p1_n_or_can):
                    continue
                return False # Collision detected
    return True # No significant collision found

def calculate_and_format_output(p1, p2, maxtotscore, maxlen):
    """
    Function to process a single pair of fragments.
    Returns the output line if criteria are met, otherwise None.
    This function contains the core logic for a single (p1, p2) pair.
    """
    # Score check
    if p1.score + p2.score > maxtotscore:
        return None

    # Collision check
    if not checkcollision(p1, p2):
        return None

    # Ensure critical atom indices for CA-CA distance are valid
    if p2.CAN == -1 or p1.CAC == -1: CAdis = float('inf')
    else: CAdis = distance2(p2.xyz[p2.CAN], p1.xyz[p1.CAC])

    if p1.CAN == -1 or p2.CAC == -1: CAdis2 = float('inf')
    else: CAdis2 = distance2(p1.xyz[p1.CAN], p2.xyz[p2.CAC])

    if not (AMMINCADIS <= CAdis <= MAXCADIS and AMMINCADIS <= CAdis2 <= MAXCADIS):
        return None

    isAM1, isAM2 = False, False
    dis, w, ang1, ang2 = 0, 0, 0, 0
    if CAdis < MINCADIS:
        # Check if indices are valid before access for AM1 calculations
        if p2.N != -1 and p1.C != -1:
            dis = distance2(p2.xyz[p2.N], p1.xyz[p1.C])
            if 1.0*1.0 < dis < 1.8*1.8:
                if all(idx != -1 for idx in [p2.CAN, p2.N, p1.C, p1.CAC]):
                    w = abs(cal_dih(p2.xyz[p2.CAN], p2.xyz[p2.N], p1.xyz[p1.C], p1.xyz[p1.CAC]))
                    ang1 = cal_angle(p2.xyz[p2.CAN], p2.xyz[p2.N], p1.xyz[p1.C])
                    ang2 = cal_angle(p2.xyz[p2.N], p1.xyz[p1.C], p1.xyz[p1.CAC])
                    if w >= 135 and 90 < ang1 < 150 and 90 < ang2 < 150:
                        isAM1 = True
                    else: return None
                else: return None # Invalid indices for dihedral/angle calculation
            else: return None
        else: return None # Invalid indices for distance calculation

    dis2, w2, ang12, ang22 = 0, 0, 0, 0
    if CAdis2 < MINCADIS:
        if p1.N != -1 and p2.C != -1:
            dis2 = distance2(p1.xyz[p1.N], p2.xyz[p2.C])
            if 1.0*1.0 < dis2 < 1.8*1.8:
                if all(idx != -1 for idx in [p1.CAN, p1.N, p2.C, p2.CAC]):
                    w2 = abs(cal_dih(p1.xyz[p1.CAN], p1.xyz[p1.N], p2.xyz[p2.C], p2.xyz[p2.CAC]))
                    ang12 = cal_angle(p1.xyz[p1.CAN], p1.xyz[p1.N], p2.xyz[p2.C])
                    ang22 = cal_angle(p1.xyz[p1.N], p2.xyz[p2.C], p2.xyz[p2.CAC])
                    if w2 >= 135 and 90 < ang12 < 150 and 90 < ang22 < 150:
                        isAM2 = True
                    else: return None
                else: return None # Invalid indices for dihedral/angle calculation
            else: return None
        else: return None # Invalid indices for distance calculation

    resn = 0
    if CAdis >= MINCADIS:
        if CAdis < TERMINALDISCAMAX3: resn = 1
        elif CAdis < TERMINALDISCAMAX4: resn = 2
        else: resn = 3
    if CAdis2 >= MINCADIS:
        if CAdis2 < TERMINALDISCAMAX3: resn += 1
        elif CAdis2 < TERMINALDISCAMAX4: resn += 2
        else: resn += 3
    if not (1 <= resn <= (maxlen - 4)): return None

    # --- Formatted Output ---
    output_line = f"{p1.name:>40s} {p2.name:>40s} {p1.score:8.3f}{p2.score:8.3f}  "
    if CAdis >= MINCADIS:
        # Ensure all required indices are valid before accessing for CAdis block
        if all(idx != -1 for idx in [p2.N, p2.CAN, p1.CAC, p1.C, p2.NC, p1.CN]):
            vals = [math.sqrt(CAdis),
                    cal_angle(p2.xyz[p2.N], p2.xyz[p2.CAN], p1.xyz[p1.CAC]),
                    cal_angle(p2.xyz[p2.CAN], p1.xyz[p1.CAC], p1.xyz[p1.C]),
                    cal_dih(p2.xyz[p2.N], p2.xyz[p2.CAN], p1.xyz[p1.CAC], p1.xyz[p1.C]),
                    cal_angle(p2.xyz[p2.NC], p2.xyz[p2.CAN], p1.xyz[p1.CAC]),
                    cal_angle(p2.xyz[p2.CAN], p1.xyz[p1.CAC], p1.xyz[p1.CN]),
                    cal_dih(p2.xyz[p2.NC], p2.xyz[p2.CAN], p1.xyz[p1.CAC], p1.xyz[p1.CN])]
            output_line += " 1: " + "".join(f"{v:8.3f}" for v in vals) + " "
        else: return None # Invalid indices for this block
    elif isAM1:
        vals = [math.sqrt(CAdis), math.sqrt(dis), w, ang1, ang2]
        output_line += " 0: " + "".join(f"{v:8.3f}" for v in vals) + "                 "

    if CAdis2 >= MINCADIS:
        # Ensure all required indices are valid before accessing for CAdis2 block
        if all(idx != -1 for idx in [p1.N, p1.CAN, p2.CAC, p2.C, p1.NC, p2.CN]):
            vals2 = [math.sqrt(CAdis2),
                    cal_angle(p1.xyz[p1.N], p1.xyz[p1.CAN], p2.xyz[p2.CAC]),
                    cal_angle(p1.xyz[p1.CAN], p2.xyz[p2.CAC], p2.xyz[p2.C]),
                    cal_dih(p1.xyz[p1.N], p1.xyz[p1.CAN], p2.xyz[p2.CAC], p2.xyz[p2.C]),
                    cal_angle(p1.xyz[p1.NC], p1.xyz[p1.CAN], p2.xyz[p2.CAC]),
                    cal_angle(p1.xyz[p1.CAN], p2.xyz[p2.CAC], p2.xyz[p2.CN]),
                    cal_dih(p1.xyz[p1.NC], p1.xyz[p1.CAN], p2.xyz[p2.CAC], p2.xyz[p2.CN])]
            output_line += "1: " + "".join(f"{v:8.3f}" for v in vals2) + " "
        elif isAM2: # This condition might be problematic if CAdis2 is >= MINCADIS AND isAM2 is True
            # This 'elif' implies it's either CAdis2 >= MINCADIS OR isAM2.
            # However, isAM2 condition check is within CAdis2 < MINCADIS block.
            # Keeping the original logic structure but be aware.
            vals2 = [math.sqrt(CAdis2), math.sqrt(dis2), w2, ang12, ang22]
            output_line += "0: " + "".join(f"{v:8.3f}" for v in vals2) + "                 "
    elif isAM2: # This handles the case where CAdis2 < MINCADIS and isAM2 is True
            vals2 = [math.sqrt(CAdis2), math.sqrt(dis2), w2, ang12, ang22]
            output_line += "0: " + "".join(f"{v:8.3f}" for v in vals2) + "                 "


    return output_line

def process_single_p1(p1_arg, p2_list_arg, maxtotscore_arg, maxlen_arg):
    """
    Worker function executed by each process.
    It takes one p1 fragment and iterates through p2_list with early exit.
    Returns a list of valid output lines found by this worker.
    """
    local_results = []
    required_p2_score_threshold = maxtotscore_arg - p1_arg.score

    for p2 in p2_list_arg:
        # Early exit condition: if p2's score alone is already too high,
        # then any subsequent p2 (since p2_list is sorted) will also be too high.
        if p1_arg.score + p2.score > maxtotscore_arg: # Double check the combined score
            break

        # Process the pair using the existing logic
        result_line = calculate_and_format_output(p1_arg, p2, maxtotscore_arg, maxlen_arg)
        if result_line:
            local_results.append(result_line)
    return local_results

def _process_p1_task_wrapper(args):
    """
    Helper function to unpack arguments for process_single_p1.
    Defined at top-level for pickling compatibility with multiprocessing.
    """
    return process_single_p1(*args)

# --- Step 6: 转换主逻辑 (Main Logic) ---
def main():
    """Main execution block, equivalent to C's main function."""
    if len(sys.argv) != 6:
        print(f"Usage: {sys.argv[0]} <fragment_list1> <fragment_list2> <max_total_score> <max_len> <output_file>", file=sys.stderr)
        sys.exit(1)

    frag_list_file1 = sys.argv[1]
    frag_list_file2 = sys.argv[2]
    maxtotscore = float(sys.argv[3])
    maxlen = int(sys.argv[4])
    output_file_path = sys.argv[5]

    # 读取第一个片段列表，并用 tqdm 显示进度
    p1_list = []
    try:
        with open(frag_list_file1, 'r') as f1:
            frag_paths1 = [line.strip() for line in f1 if line.strip()]
        for pdb_file in tqdm(frag_paths1, desc="Reading fragments 1"):
            p = read_structure(pdb_file)
            if p: # Only add if read_structure succeeded
                p1_list.append(p)
            if len(p1_list) >= MAXFRAGNUM:
                print(f"Warning: Max fragment limit reached for list 1: {MAXFRAGNUM}", file=sys.stderr)
                break
    except FileNotFoundError:
        print(f"ERROR: Cannot open file {frag_list_file1}", file=sys.stderr)
        sys.exit(1)

    # 读取第二个片段列表，并用 tqdm 显示进度
    p2_list = []
    try:
        with open(frag_list_file2, 'r') as f2:
            frag_paths2 = [line.strip() for line in f2 if line.strip()]
        for pdb_file in tqdm(frag_paths2, desc="Reading fragments 2"):
            p = read_structure(pdb_file)
            if p: # Only add if read_structure succeeded
                p2_list.append(p)
            if len(p2_list) >= MAXFRAGNUM:
                print(f"Warning: Max fragment limit reached for list 2: {MAXFRAGNUM}", file=sys.stderr)
                break
    except FileNotFoundError:
        print(f"ERROR: Cannot open file {frag_list_file2}", file=sys.stderr)
        sys.exit(1)

    # We are now submitting each p1 to a separate process as a "task".
    # Each process will then iterate through p2_list internally.
    total_p1_fragments = len(p1_list)

    try:
        with open(output_file_path, 'w') as outfile:
            # You can set a specific number of cores or use os.cpu_count()
            num_cores = 12 # Example: Set to 12 as per your output, or os.cpu_count() for max
            print(f"Using {num_cores} CPU cores for parallel processing.")

            with concurrent.futures.ProcessPoolExecutor(max_workers=num_cores) as executor:
                # Prepare arguments for each p1 fragment
                # We need to pass p2_list, maxtotscore, and maxlen to each worker
                # as they are needed in process_single_p1
                tasks_for_p1s = [(p1, p2_list, maxtotscore, maxlen) for p1 in p1_list]

                # Map the _process_p1_task_wrapper function to all tasks (each task is one p1 fragment)
                # tqdm wraps the future results, showing progress for each p1 processed
                results_iterator = tqdm(
                    executor.map(_process_p1_task_wrapper, tasks_for_p1s),
                    total=total_p1_fragments,
                    desc="Analyzing p1 fragments (parallel)"
                )

                # Each 'result' here is a list of output lines returned by one process_single_p1 call
                for results_from_one_p1 in results_iterator:
                    for line in results_from_one_p1:
                        outfile.write(line + "\n")

    except IOError as e:
        print(f"ERROR: Could not open output file '{output_file_path}' for writing: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred during parallel processing: {e}", file=sys.stderr)
        sys.exit(1)


    print(f"\nProcessing complete. Results have been written to '{output_file_path}'.")

if __name__ == "__main__":
    # This block is essential for multiprocessing on Windows
    # and good practice on other OS for robustness.
    import multiprocessing
    multiprocessing.freeze_support() # Recommended for pyinstaller/executables
    main()