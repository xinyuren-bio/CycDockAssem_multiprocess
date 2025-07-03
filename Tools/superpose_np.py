import numpy as np

def fit_terminal(coords_ref, coords_mob):
    if coords_ref.shape != coords_mob.shape or coords_ref.shape[1] != 3:
        return None, None, 999.9 
    num_atoms = coords_ref.shape[0]
    
    # 1. 对两组坐标都进行中心化
    centroid_ref = np.mean(coords_ref, axis=0)
    centroid_mob = np.mean(coords_mob, axis=0)
    coords_ref_c = coords_ref - centroid_ref
    coords_mob_c = coords_mob - centroid_mob

    # 2. 计算协方差矩阵 H
    H = np.dot(coords_mob_c.T, coords_ref_c)

    # 3. SVD分解
    U, S, Vt = np.linalg.svd(H)

    # 4. 计算旋转矩阵 R
    R = np.dot(Vt.T, U.T)
    if np.linalg.det(R) < 0:
        Vt_corrected = Vt.T.copy()
        Vt_corrected[:, -1] *= -1
        R = np.dot(Vt_corrected, U.T)

    # 5. 计算最终用于对齐的平移向量
    t = centroid_ref - np.dot(R, centroid_mob)
    
    # 6. 计算RMSD
    coords_mob_aligned = np.dot(coords_mob, R.T) + t
    diff = coords_ref - coords_mob_aligned
    rmsd = np.sqrt(np.sum(diff * diff) / num_atoms)
    
    return R, t, rmsd