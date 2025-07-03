import numpy as np
from scipy.spatial.distance import cdist

def _calculate_steric_score(dist_matrix, config):
    score = 0.0
    accept_steric = config['ACCEPT_STERIC']
    if dist_matrix.size == 0: return 0.0
    score += np.sum(dist_matrix < 2.4) * (accept_steric * 2.0)
    mask = (dist_matrix >= 2.4) & (dist_matrix < 3.0)
    score += np.sum(1.0 - (dist_matrix[mask] - 2.4) / 0.6)
    return score

def check_steric(patch, p1, p2, sup_type, config):
    total_steric_score = 0.0
    accept_steric = config['ACCEPT_STERIC']
    accept_tersteric = 1.5
    
    # 使用 sup_type=0 的核心原子
    patch_core_coords = patch.coords[patch.get_core_indices(sup_type=0)]
    if patch_core_coords.size == 0:
        return accept_steric + 1.0, [1, 1]

    # --- 1. 主体碰撞检查 ---
    p2_clash_coords = p2.coords[p2.get_p2_clash_indices()]
    p1_clash_coords = p1.coords[p1.get_p1_clash_indices()]

    if p2_clash_coords.size > 0:
        dist_matrix_p2 = cdist(patch_core_coords, p2_clash_coords)
        total_steric_score += _calculate_steric_score(dist_matrix_p2, config)
        if total_steric_score > accept_steric: return total_steric_score, [1, 1]
    
    if p1_clash_coords.size > 0:
        dist_matrix_p1 = cdist(patch_core_coords, p1_clash_coords)
        total_steric_score += _calculate_steric_score(dist_matrix_p1, config)
        if total_steric_score > accept_steric: return total_steric_score, [1, 1]
    
    ter_flags = [1, 1]
    # 检查 patch vs p1的N端侧链
    p1_n_term_sc_coords = p1.coords[p1.get_n_term_sidechain_indices()]
    if p1_n_term_sc_coords.size > 0:
        dist_ter_p1 = cdist(patch_core_coords, p1_n_term_sc_coords)
        if np.any(dist_ter_p1 < 2.2):
            ter_flags[0] = 0
        else:
            mask = dist_ter_p1 < 2.8
            tersteric_score = np.sum(1.0 - (dist_ter_p1[mask] - 2.2) / 0.6)
            if tersteric_score > accept_tersteric:
                ter_flags[0] = 0

    # 检查 patch vs p2的C端侧链
    p2_c_term_sc_coords = p2.coords[p2.get_c_term_sidechain_indices()]
    if p2_c_term_sc_coords.size > 0:
        dist_ter_p2 = cdist(patch_core_coords, p2_c_term_sc_coords)
        if np.any(dist_ter_p2 < 2.2):
            ter_flags[1] = 0
        else:
            mask = dist_ter_p2 < 2.8
            tersteric_score = np.sum(1.0 - (dist_ter_p2[mask] - 2.2) / 0.6)
            if tersteric_score > accept_tersteric:
                ter_flags[1] = 0
    
    return total_steric_score, ter_flags


def check_steric2_np(patch, p1, p2, config):
    total_steric_score = 0.0
    accept_steric = config['ACCEPT_STERIC']
    accept_tersteric = 1.5

    # 1. 主体碰撞检查 (使用 sup_type=1 的核心原子)
    patch_core_coords = patch.coords[patch.get_core_indices(sup_type=1)]
    if patch_core_coords.size == 0: return accept_steric + 1.0, [0, 0]

    # 1a. patch vs p2 (使用 type2 的索引方法)
    p2_clash_coords = p2.coords[p2.get_p2_clash_indices_type2()]
    if p2_clash_coords.size > 0:
        dist_matrix_p2 = cdist(patch_core_coords, p2_clash_coords)
        total_steric_score += _calculate_steric_score(dist_matrix_p2, config)
        if total_steric_score > accept_steric: return total_steric_score, [0, 0]

    # 1b. patch vs p1 (使用 type2 的索引方法)
    p1_clash_coords = p1.coords[p1.get_p1_clash_indices_type2()]
    if p1_clash_coords.size > 0:
        dist_matrix_p1 = cdist(patch_core_coords, p1_clash_coords)
        total_steric_score += _calculate_steric_score(dist_matrix_p1, config)
        if total_steric_score > accept_steric: return total_steric_score, [0, 0]

    # 2. 末端精细碰撞检查
    ter_flags = [0, 0]

    # 2a. 检查 patch 和 p1 的连接端
    try:
        patch_term_N_idx = patch.residue_info[patch.residues[-1]]['atoms'].get('N')
        p1_term_N_idx = p1.residue_info[p1.residues[0]]['atoms'].get('N')
        d = np.linalg.norm(patch.coords[patch_term_N_idx] - p1.coords[p1_term_N_idx])
        
        if d < 0.8:
            ter_flags[0] = 1 # 距离足够近，进行下一步软碰撞检查
            p1_n_term_sc_coords = p1.coords[p1.get_n_term_sidechain_indices()]
            if p1_n_term_sc_coords.size > 0:
                dist_ter_p1 = cdist(patch_core_coords, p1_n_term_sc_coords)
                if np.any(dist_ter_p1 < 2.2): ter_flags[0] = 0
                else:
                    mask = dist_ter_p1 < 2.8
                    if np.sum(1.0 - (dist_ter_p1[mask] - 2.2) / 0.6) > accept_tersteric:
                        ter_flags[0] = 0
    except (KeyError, IndexError, TypeError):
        ter_flags[0] = 0

    # 2b. 检查 patch 和 p2 的连接端
    try:
        patch_term_C_idx = patch.residue_info[patch.residues[0]]['atoms'].get('C')
        p2_term_C_idx_in_p1 = p1.residue_info[p1.residues[p2.residue_count - 1]]['atoms'].get('C')
        d = np.linalg.norm(patch.coords[patch_term_C_idx] - p1.coords[p2_term_C_idx_in_p1])

        if d < 0.8:
            ter_flags[1] = 1 # 距离足够近，进行下一步软碰撞检查
            p2_c_term_sc_coords = p2.coords[p2.get_c_term_sidechain_indices()]
            if p2_c_term_sc_coords.size > 0:
                dist_ter_p2 = cdist(patch_core_coords, p2_c_term_sc_coords)
                if np.any(dist_ter_p2 < 2.2): ter_flags[1] = 0
                else:
                    mask = dist_ter_p2 < 2.8
                    if np.sum(1.0 - (dist_ter_p2[mask] - 2.2) / 0.6) > accept_tersteric:
                        ter_flags[1] = 0
    except (KeyError, IndexError, TypeError):
        ter_flags[1] = 0

    return total_steric_score, ter_flags


def dock_to_target(patch, target_coords, config):
    if patch.coords is None or target_coords is None or target_coords.size == 0:
        return 0.0
    patch_core_coords = patch.coords[patch.get_core_indices(sup_type=0)] # 假设用sup_type=0的核心
    if patch_core_coords.size == 0:
        return 0.0
    dist_matrix = cdist(patch_core_coords, target_coords)
    score = np.sum(1.0 - (dist_matrix[dist_matrix < 3.0] - 2.4) / 0.6)
    return score if score > 0 else 0.0

def _calculate_dock_score(fragment_coords, target_coords):
    if fragment_coords.size == 0 or target_coords.size == 0:
        return 0.0
    
    dist_matrix = cdist(fragment_coords, target_coords)
    
    # 吸引力项
    attr = np.sum((dist_matrix > 3.2) & (dist_matrix < 4.2))
    
    # 位阻项
    steric_mask = dist_matrix < 3.0
    steric = np.sum(1.0 - (dist_matrix[steric_mask] - 2.4) / 0.6)
    
    # C代码中的加权
    score = steric * 0.4 - attr * 0.1
    return score

def pepterdock(p1, p2, target_coords):
    scores = [0.0] * 4

    # score[0]: p1 N端侧链 vs target
    p1_n_term_sc = p1.coords[p1.get_n_term_sidechain_indices()]
    scores[0] = _calculate_dock_score(p1_n_term_sc, target_coords)

    # score[1]: p1 C端侧链 vs target
    p1_c_term_sc = p1.coords[p1.get_c_term_sidechain_indices()]
    scores[1] = _calculate_dock_score(p1_c_term_sc, target_coords)
    
    # score[2]: p2 N端侧链 vs target
    p2_n_term_sc = p2.coords[p2.get_n_term_sidechain_indices()]
    scores[2] = _calculate_dock_score(p2_n_term_sc, target_coords)
    
    # score[3]: p2 C端侧链 vs target
    p2_c_term_sc = p2.coords[p2.get_c_term_sidechain_indices()]
    scores[3] = _calculate_dock_score(p2_c_term_sc, target_coords)
    
    return scores