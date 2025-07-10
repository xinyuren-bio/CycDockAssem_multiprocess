import numpy as np
import argparse
import sys

# --- 数学计算函数 (无需改动) ---

def distance(p1: np.ndarray, p2: np.ndarray) -> float:
    """计算两点之间的欧几里得距离"""
    return np.linalg.norm(p1 - p2)

def calculate_angle(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    """
    计算由三个点p1-p2-p3构成的角度，p2为顶点。
    结果以度为单位。
    """
    v1 = p1 - p2
    v2 = p3 - p2
    
    # 向量长度为0的特殊情况
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    if norm_v1 == 0 or norm_v2 == 0:
        return 0.0

    cos_angle = np.dot(v1, v2) / (norm_v1 * norm_v2)
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    
    angle_rad = np.arccos(cos_angle)
    return np.rad2deg(angle_rad)

def calculate_dihedral(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray) -> float:
    """
    计算由四个点p1-p2-p3-p4构成的二面角。
    结果以度为单位。
    """
    v12 = p1 - p2
    v32 = p3 - p2
    v34 = p3 - p4

    m = np.cross(v12, v32)
    n = np.cross(v32, v34)

    # 向量长度为0的特殊情况
    norm_m = np.linalg.norm(m)
    norm_n = np.linalg.norm(n)
    if norm_m == 0 or norm_n == 0:
        return 0.0

    cos_phi = np.dot(m, n) / (norm_m * norm_n)
    cos_phi = np.clip(cos_phi, -1.0, 1.0)
    
    phi_rad = np.arccos(cos_phi)

    sign = np.sign(np.dot(v12, n))
    if sign == 0:
        sign = 1.0

    signed_phi_rad = sign * phi_rad
    return np.rad2deg(signed_phi_rad)

# --- 文件读取函数 (已修正) ---

def read_peptide_terminals_from_pdb(pdb_file: str, chain_length: int) -> dict:
    """
    从PDB文件中读取特定残基的原子坐标。
    此版本使用相对残基计数，不依赖PDB中的绝对编号。
    """
    atoms = {
        'N1': None, 'CA1': None, 'C1': None,
        'CA2': None,
        'Nn': None, 'CAn': None, 'Cn': None
    }

    relative_res_num = 0
    last_res_id = None # 用于跟踪上一个残基的标识符

    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if not line.startswith("ATOM"):
                    continue

                # PDB中残基的唯一标识符通常是链ID+残基编号+插入码 (这里简化为残基编号+名称部分)
                # C代码比较的是 line[17:26]，我们就用这个
                current_res_id = line[17:26]
                if current_res_id != last_res_id:
                    relative_res_num += 1
                    last_res_id = current_res_id

                atom_name = line[12:16].strip()
                
                try:
                    coords = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except ValueError:
                    continue

                # 使用相对残基编号进行判断
                if relative_res_num == 1:
                    if atom_name == 'N': atoms['N1'] = coords
                    elif atom_name == 'CA': atoms['CA1'] = coords
                    elif atom_name == 'C': atoms['C1'] = coords
                elif relative_res_num == 2 and atom_name == 'CA':
                    atoms['CA2'] = coords
                
                if relative_res_num == chain_length:
                    if atom_name == 'N': atoms['Nn'] = coords
                    elif atom_name == 'CA': atoms['CAn'] = coords
                    elif atom_name == 'C': atoms['Cn'] = coords

    except FileNotFoundError:
        print(f"错误: 无法打开蛋白质结构文件 {pdb_file}")
        sys.exit(1)
        
    return atoms

# --- 主函数 (无需改动) ---

def main():
    """主执行函数"""
    parser = argparse.ArgumentParser(
        description="计算多肽链末端之间的几何参数 (距离、角度、二面角)。"
    )
    parser.add_argument("pdb_file", type=str, help="输入的PDB文件名")
    parser.add_argument("length", type=int, help="蛋白质/多肽的残基总数")

    args = parser.parse_args()
    terminals = read_peptide_terminals_from_pdb(args.pdb_file, args.length)

    required_atoms = ['N1', 'C1', 'CA1', 'CAn', 'Cn', 'Nn']
    missing_atoms = [key for key in required_atoms if terminals[key] is None]
    
    if missing_atoms:
        print(f"错误: 未能在PDB文件中找到必需的原子: {', '.join(missing_atoms)}。请检查文件和肽链长度({args.length})是否正确。")
        sys.exit(1)

    N1, C1, CA1 = terminals['N1'], terminals['C1'], terminals['CA1']
    CAn, Cn, Nn = terminals['CAn'], terminals['Cn'], terminals['Nn']

    ca_dis = distance(CA1, CAn)
    n1_ca_can = calculate_angle(N1, CA1, CAn)
    ca1_can_cn = calculate_angle(CA1, CAn, Cn)
    dih1 = calculate_dihedral(N1, CA1, CAn, Cn)
    nn_can_ca1 = calculate_angle(Nn, CAn, CA1)
    can_ca1_c1 = calculate_angle(CAn, CA1, C1)
    dih2 = calculate_dihedral(Nn, CAn, CA1, C1)
    
    filename = args.pdb_file
    print(f"{filename} CA-CA {ca_dis:9.3f}")
    print(f"{filename} N1-Cn {n1_ca_can:9.3f}{ca1_can_cn:9.3f}{dih1:9.3f}")
    print(f"{filename} C1-Nn {nn_can_ca1:9.3f}{can_ca1_c1:9.3f}{dih2:9.3f}")


if __name__ == "__main__":
    main()