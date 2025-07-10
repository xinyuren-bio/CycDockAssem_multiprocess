import numpy as np
from scipy.spatial.transform import Rotation as R
import argparse

def rotate_protein_sequential(input_pdb_path, ref_atom1_coords, ref_atom2_coords, ref_atom3_coords, output_pdb_path):
    """
    使用与提供的C代码完全相同的串行旋转算法和浮点数精度来旋转蛋白质。
    """
    protein_atoms = []
    other_lines = []

    with open(input_pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    atom_info = {
                        'line': line,
                        'x': float(line[30:38]),
                        'y': float(line[38:46]),
                        'z': float(line[46:54]),
                        'original_line_start': line[:30],
                        'original_line_end': line[54:],
                    }
                    protein_atoms.append(atom_info)
                except ValueError:
                    other_lines.append(line)
            else:
                other_lines.append(line)

    if not protein_atoms:
        return

    # --- 关键改动：强制使用 np.float32 来匹配C代码的 `float` 类型 ---
    FLOAT_TYPE = np.float32
    
    all_atom_coords = np.array([[atom['x'], atom['y'], atom['z']] for atom in protein_atoms], dtype=FLOAT_TYPE)
    ref_coords = np.array([ref_atom1_coords, ref_atom2_coords, ref_atom3_coords], dtype=FLOAT_TYPE)

    center = (ref_coords[0] + ref_coords[1]) / 2.0
    
    coords_translated = all_atom_coords - center
    ref_coords_translated = ref_coords - center

    # 第一次旋转
    target_x_axis = np.array([1.0, 0.0, 0.0], dtype=FLOAT_TYPE)
    current_x_axis = ref_coords_translated[1] - ref_coords_translated[0]
    current_x_axis /= np.linalg.norm(current_x_axis)

    rot_axis_1 = np.cross(target_x_axis, current_x_axis)
    
    coords_rotated_1 = coords_translated
    ref_coords_rotated_1 = ref_coords_translated
    
    if np.linalg.norm(rot_axis_1) > 1e-6: # 使用1e-6以匹配float的精度
        rot_axis_1 /= np.linalg.norm(rot_axis_1)
        angle_1 = np.arccos(np.dot(target_x_axis, current_x_axis))
        rotation_1 = R.from_rotvec(rot_axis_1 * -angle_1)
        coords_rotated_1 = rotation_1.apply(coords_translated)
        ref_coords_rotated_1 = rotation_1.apply(ref_coords_translated)

    # 第二次旋转
    target_z_axis = np.array([0.0, 0.0, 1.0], dtype=FLOAT_TYPE)
    vec_to_ref2 = ref_coords_rotated_1[1]
    vec_to_ref3 = ref_coords_rotated_1[2]
    
    current_z_axis = np.cross(vec_to_ref2 / np.linalg.norm(vec_to_ref2), 
                              vec_to_ref3 / np.linalg.norm(vec_to_ref3))
    current_z_axis /= np.linalg.norm(current_z_axis)
    
    rot_axis_2 = np.cross(target_z_axis, current_z_axis)
    
    coords_final = coords_rotated_1
    ref_coords_final = ref_coords_rotated_1
    
    if np.linalg.norm(rot_axis_2) > 1e-6:
        rot_axis_2 /= np.linalg.norm(rot_axis_2)
        angle_2 = np.arccos(np.dot(target_z_axis, current_z_axis))
        rotation_2 = R.from_rotvec(rot_axis_2 * -angle_2)
        coords_final = rotation_2.apply(coords_rotated_1)
        ref_coords_final = rotation_2.apply(ref_coords_rotated_1)

    # 输出结果
    with open(output_pdb_path, 'w') as f:
        for line in other_lines:
            f.write(line)
        for i, atom_info in enumerate(protein_atoms):
            x_rot, y_rot, z_rot = coords_final[i]
            formatted_x = f"{x_rot:8.3f}"
            formatted_y = f"{y_rot:8.3f}"
            formatted_z = f"{z_rot:8.3f}"
            new_line = (
                atom_info['original_line_start'] +
                formatted_x + formatted_y + formatted_z +
                atom_info['original_line_end']
            )
            f.write(new_line.rstrip() + '\n')
            
    print("-" * 30)
    print(f"蛋白质已按C代码串行逻辑和32位浮点数精度旋转并保存到: {output_pdb_path}")
    print("旋转后参考原子的最终坐标:")
    print(f"box1       {ref_coords_final[0][0]:8.3f}{ref_coords_final[0][1]:8.3f}{ref_coords_final[0][2]:8.3f}")
    print(f"box2       {ref_coords_final[1][0]:8.3f}{ref_coords_final[1][1]:8.3f}{ref_coords_final[1][2]:8.3f}")
    print(f"3atm       {ref_coords_final[2][0]:8.3f}{ref_coords_final[2][1]:8.3f}{ref_coords_final[2][2]:8.3f}")
    
    box2_minus_box1 = ref_coords_final[1] - ref_coords_final[0]
    print(f"box2-box1  {box2_minus_box1[0]:8.3f}{box2_minus_box1[1]:8.3f}{box2_minus_box1[2]:8.3f}")

# main函数部分保持不变
if __name__ == "__main__":
    # ... (这部分不需要修改)
    parser = argparse.ArgumentParser(description='根据三个参考原子，使用与C代码兼容的串行算法和浮点数精度旋转蛋白质。')
    parser.add_argument('--input_pdb', type=str, required=True, help='输入蛋白质 PDB 文件的路径。')
    parser.add_argument('--output_pdb', type=str, required=True, help='输出旋转后 PDB 文件的路径。')
    parser.add_argument('--ref1', type=str, required=True, help='第一个参考原子的 X,Y,Z 坐标 (例如: "1.0,2.0,3.0")。')
    parser.add_argument('--ref2', type=str, required=True, help='第二个参考原子的 X,Y,Z 坐标 (例如: "4.0,5.0,6.0")。')
    parser.add_argument('--ref3', type=str, required=True, help='第三个参考原子的 X,Y,Z 坐标 (例如: "7.0,8.0,9.0")。')

    args = parser.parse_args()
    def parse_coords(coord_str):
        try: return [float(c) for c in coord_str.split(',')]
        except ValueError: raise argparse.ArgumentTypeError(f"无效的坐标格式: {coord_str}。")
    ref_atom1_coords = parse_coords(args.ref1)
    ref_atom2_coords = parse_coords(args.ref2)
    ref_atom3_coords = parse_coords(args.ref3)
    rotate_protein_sequential(args.input_pdb, ref_atom1_coords, ref_atom2_coords, ref_atom3_coords, args.output_pdb)