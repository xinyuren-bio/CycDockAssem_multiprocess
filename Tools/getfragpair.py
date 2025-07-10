#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script is a Python conversion of the original C code.
It performs geometric analysis on protein fragments to find potential linkers.
The results are written to a specified output file.
"""

import sys
import math
import numpy as np

# --- (前面部分的代码保持不变) ---

# --- Step 1 & 2: 导入库和定义常量 (Imports and Constants) ---
MAXCADIS = 13.5 * 13.5
MINCADIS = 4.7 * 4.7
AMMINCADIS = 3.4 * 3.4
TERMINALDISCAMAX3 = 7.5 * 7.5
TERMINALDISCAMAX4 = 10.5 * 10.5
MAXFRAGNUM = 2000

# --- Step 3: 创建 Protein 类 (Protein Class) ---
class Protein:
    def __init__(self):
        self.an = 0
        self.xyz = np.empty((0, 3), dtype=np.float32)
        self.name = ""
        self.score = 0.0
        self.N, self.NC, self.C, self.CN, self.CAN, self.CAC = 0, 0, 0, 0, 0, 0

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
    return math.degrees(math.acos(cos_angle(v12, v32)))

def cal_dih(x1, x2, x3, x4):
    x12, x32, x34 = getdx(x1, x2), getdx(x3, x2), getdx(x3, x4)
    m, n = oprod(x12, x32), oprod(x32, x34)
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
            p.name = pn.strip()
            xyz_coords, atom_count, res_num = [], 0, 1
            for line in pf:
                if line.startswith("ATOM"):
                    coords = getxyz(line)
                    if coords is None: continue
                    xyz_coords.append(coords)
                    atom_type = line[12:16].strip()
                    if atom_type == "N":
                        if res_num == 1: p.N = atom_count
                        p.CN = atom_count
                    elif atom_type == "CA":
                        if res_num == 1: p.CAN = atom_count
                        p.CAC = atom_count
                    elif atom_type == "C":
                        if res_num == 1:
                            p.NC = atom_count
                            res_num += 1
                        p.C = atom_count
                    atom_count += 1
                elif line.startswith("REMARK    score"):
                    try: p.score = float(line[15:])
                    except (ValueError, IndexError): p.score = 0.0
            p.an, p.xyz = atom_count, np.array(xyz_coords, dtype=np.float32)
            return p
    except FileNotFoundError:
        print(f"ERROR: Cannot open structure file {pn}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while reading {pn}: {e}", file=sys.stderr)
        sys.exit(1)

def checkcollision(p1, p2):
    CLASH_DISTANCE_SQ = 2.5 * 2.5
    for ia in range(p1.an):
        for ja in range(p2.an):
            if distance2(p1.xyz[ia], p2.xyz[ja]) < CLASH_DISTANCE_SQ:
                if (ia in {p1.C, p1.CAC} and ja in {p2.N, p2.CAN}) or \
                   (ja in {p2.C, p2.CAC} and ia in {p1.N, p1.CAN}):
                    continue
                return False
    return True

# --- Step 6: 转换主逻辑 (Main Logic) ---
def main():
    """Main execution block, equivalent to C's main function."""
    # MODIFICATION: Check for 6 arguments now
    if len(sys.argv) != 6:
        print(f"Usage: {sys.argv[0]} <fragment_list1> <fragment_list2> <max_total_score> <max_len> <output_file>", file=sys.stderr)
        sys.exit(1)

    frag_list_file1 = sys.argv[1]
    frag_list_file2 = sys.argv[2]
    maxtotscore = float(sys.argv[3])
    maxlen = int(sys.argv[4])
    # MODIFICATION: Get the output file path from arguments
    output_file_path = sys.argv[5]

    # Read the first list of protein fragments
    p1_list = []
    try:
        with open(frag_list_file1, 'r') as f1:
            for line in f1:
                if pdb_file := line.strip():
                    p1_list.append(read_structure(pdb_file))
                    if len(p1_list) >= MAXFRAGNUM:
                        print(f"Warning: Max fragment limit reached for list 1: {MAXFRAGNUM}", file=sys.stderr)
                        break
    except FileNotFoundError:
        print(f"ERROR: Cannot open file {frag_list_file1}", file=sys.stderr)
        sys.exit(1)

    # Read the second list of protein fragments
    p2_list = []
    try:
        with open(frag_list_file2, 'r') as f2:
            for line in f2:
                if pdb_file := line.strip():
                    p2_list.append(read_structure(pdb_file))
                    if len(p2_list) >= MAXFRAGNUM:
                        print(f"Warning: Max fragment limit reached for list 2: {MAXFRAGNUM}", file=sys.stderr)
                        break
    except FileNotFoundError:
        print(f"ERROR: Cannot open file {frag_list_file2}", file=sys.stderr)
        sys.exit(1)
        
    # MODIFICATION: Open the output file and run the main loop inside the 'with' block
    try:
        with open(output_file_path, 'w') as outfile:
            # --- Main processing loop ---
            for p1 in p1_list:
                for p2 in p2_list:
                    # (The entire inner logic of the loop remains the same)
                    if p1.score + p2.score > maxtotscore:
                        continue
                    if not checkcollision(p1, p2):
                        continue
                    
                    CAdis = distance2(p2.xyz[p2.CAN], p1.xyz[p1.CAC])
                    CAdis2 = distance2(p1.xyz[p1.CAN], p2.xyz[p2.CAC])

                    if not (AMMINCADIS <= CAdis <= MAXCADIS and AMMINCADIS <= CAdis2 <= MAXCADIS):
                        continue
                    
                    isAM1, isAM2 = False, False
                    if CAdis < MINCADIS:
                        dis = distance2(p2.xyz[p2.N], p1.xyz[p1.C])
                        if 1.0*1.0 < dis < 1.8*1.8:
                            w = abs(cal_dih(p2.xyz[p2.CAN], p2.xyz[p2.N], p1.xyz[p1.C], p1.xyz[p1.CAC]))
                            ang1 = cal_angle(p2.xyz[p2.CAN], p2.xyz[p2.N], p1.xyz[p1.C])
                            ang2 = cal_angle(p2.xyz[p2.N], p1.xyz[p1.C], p1.xyz[p1.CAC])
                            if w >= 135 and 90 < ang1 < 150 and 90 < ang2 < 150:
                                isAM1 = True
                            else: continue
                        else: continue
                    
                    if CAdis2 < MINCADIS:
                        dis2 = distance2(p1.xyz[p1.N], p2.xyz[p2.C])
                        if 1.0*1.0 < dis2 < 1.8*1.8:
                            w2 = abs(cal_dih(p1.xyz[p1.CAN], p1.xyz[p1.N], p2.xyz[p2.C], p2.xyz[p2.CAC]))
                            ang12 = cal_angle(p1.xyz[p1.CAN], p1.xyz[p1.N], p2.xyz[p2.C])
                            ang22 = cal_angle(p1.xyz[p1.N], p2.xyz[p2.C], p2.xyz[p2.CAC])
                            if w2 >= 135 and 90 < ang12 < 150 and 90 < ang22 < 150:
                                isAM2 = True
                            else: continue
                        else: continue
                    
                    resn = 0
                    if CAdis >= MINCADIS:
                        if CAdis < TERMINALDISCAMAX3: resn = 1
                        elif CAdis < TERMINALDISCAMAX4: resn = 2
                        else: resn = 3
                    if CAdis2 >= MINCADIS:
                        if CAdis2 < TERMINALDISCAMAX3: resn += 1
                        elif CAdis2 < TERMINALDISCAMAX4: resn += 2
                        else: resn += 3
                    if not (1 <= resn <= (maxlen - 4)): continue
                    
                    # --- Formatted Output ---
                    output_line = f"{p1.name:>40s} {p2.name:>40s} {p1.score:8.3f}{p2.score:8.3f}  "
                    if CAdis >= MINCADIS:
                        vals = [math.sqrt(CAdis), 
                                cal_angle(p2.xyz[p2.N], p2.xyz[p2.CAN], p1.xyz[p1.CAC]),
                                cal_angle(p2.xyz[p2.CAN], p1.xyz[p1.CAC], p1.xyz[p1.C]),
                                cal_dih(p2.xyz[p2.N], p2.xyz[p2.CAN], p1.xyz[p1.CAC], p1.xyz[p1.C]),
                                cal_angle(p2.xyz[p2.NC], p2.xyz[p2.CAN], p1.xyz[p1.CAC]),
                                cal_angle(p2.xyz[p2.CAN], p1.xyz[p1.CAC], p1.xyz[p1.CN]),
                                cal_dih(p2.xyz[p2.NC], p2.xyz[p2.CAN], p1.xyz[p1.CAC], p1.xyz[p1.CN])]
                        output_line += " 1: " + "".join(f"{v:8.3f}" for v in vals) + " "
                    elif isAM1:
                        vals = [math.sqrt(CAdis), math.sqrt(dis), w, ang1, ang2]
                        output_line += " 0: " + "".join(f"{v:8.3f}" for v in vals) + "                 "
                    
                    if CAdis2 >= MINCADIS:
                        vals2 = [math.sqrt(CAdis2),
                                cal_angle(p1.xyz[p1.N], p1.xyz[p1.CAN], p2.xyz[p2.CAC]),
                                cal_angle(p1.xyz[p1.CAN], p2.xyz[p2.CAC], p2.xyz[p2.C]),
                                cal_dih(p1.xyz[p1.N], p1.xyz[p1.CAN], p2.xyz[p2.CAC], p2.xyz[p2.C]),
                                cal_angle(p1.xyz[p1.NC], p1.xyz[p1.CAN], p2.xyz[p2.CAC]),
                                cal_angle(p1.xyz[p1.CAN], p2.xyz[p2.CAC], p2.xyz[p2.CN]),
                                cal_dih(p1.xyz[p1.NC], p1.xyz[p1.CAN], p2.xyz[p2.CAC], p2.xyz[p2.CN])]
                        output_line += "1: " + "".join(f"{v:8.3f}" for v in vals2) + " "
                    elif isAM2:
                        vals2 = [math.sqrt(CAdis2), math.sqrt(dis2), w2, ang12, ang22]
                        output_line += "0: " + "".join(f"{v:8.3f}" for v in vals2) + "                 "
                    
                    # MODIFICATION: Print the final line to the specified file
                    print(output_line, file=outfile)

    except IOError as e:
        # If the output file cannot be opened, print an error to the screen
        print(f"ERROR: Could not open output file '{output_file_path}' for writing: {e}", file=sys.stderr)
        sys.exit(1)

    # Let the user know the process is complete.
    print(f"Processing complete. Results have been written to '{output_file_path}'.")

if __name__ == "__main__":
    main()