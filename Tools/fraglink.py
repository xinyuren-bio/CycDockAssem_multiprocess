#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
蛋白质片段链接程序 - Python版本 (预加载几何库 + 日志 + 并行)
"""

import os
import sys
import math
import numpy as np
import multiprocessing
import logging
import time
import datetime
import copy
from typing import List, Tuple, Dict
from tqdm import tqdm


ACCEPTSTERIC = 4.0; ACCEPTTERSTERIC = 1.5; ACCEPTDOCKSTERIC = 5.0; ACCEPTRMSD = 0.8; ACCEPTDOCKING = 2.0; DISDIFF = 0.5; ANG1DIFF = 15; ANG2DIFF = 15; DIHDIFF = 35; BESTPATCHNUM = 500; BESTNUM = 30; BESTCLUSTERNUM = 3; CULSTERRMSD = 0.6; TODEG = lambda A: A * 57.295779513

class Residue:
    def __init__(self): self.start = 0; self.CA = 0; self.C = 0; self.N = 0; self.end = 0; self.nam = ""
class Fragment:
    def __init__(self): self.an = 0; self.xyz = np.zeros((100, 3), dtype=np.float64); self.nam = [""] * 100; self.resn = 0; self.res = [Residue() for _ in range(5)]; self.N = 0; self.NC = 0; self.C = 0; self.CN = 0; self.CAN = 0; self.CAC = 0
class Geometry:
    @staticmethod
    def distance(x1: np.ndarray, x2: np.ndarray) -> float: return np.linalg.norm(x1 - x2)
    @staticmethod
    def distance2(x1: np.ndarray, x2: np.ndarray) -> float: return np.sum((x1 - x2) ** 2)
    @staticmethod
    def cal_angle(x1: np.ndarray, x2: np.ndarray, x3: np.ndarray) -> float:
        v1 = x1 - x2; v2 = x3 - x2; norm_v1 = np.linalg.norm(v1); norm_v2 = np.linalg.norm(v2)
        if norm_v1 == 0 or norm_v2 == 0: return 0.0
        cos_angle = np.clip(np.dot(v1, v2) / (norm_v1 * norm_v2), -1.0, 1.0)
        return TODEG(np.arccos(cos_angle))
    @staticmethod
    def cal_dih(x1: np.ndarray, x2: np.ndarray, x3: np.ndarray, x4: np.ndarray) -> float:
        v1 = x2 - x1; v2 = x3 - x2; v3 = x4 - x3; n1 = np.cross(v1, v2); n2 = np.cross(v2, v3); norm_n1 = np.linalg.norm(n1); norm_n2 = np.linalg.norm(n2)
        if norm_n1 == 0 or norm_n2 == 0: return 0.0
        cos_dih = np.clip(np.dot(n1, n2) / (norm_n1 * norm_n2), -1.0, 1.0)
        dih_rad = np.arccos(cos_dih)
        if np.dot(n1, v3) < 0: dih_rad = -dih_rad
        return TODEG(dih_rad)
class Superpose:
    @staticmethod
    def get_center(x: np.ndarray) -> np.ndarray: return np.mean(x, axis=0)
    @staticmethod
    def do_rot(x: np.ndarray, R: np.ndarray) -> np.ndarray: return x @ R.T
    @staticmethod
    def fit_terminal(target_coords: np.ndarray, mobile_coords: np.ndarray) -> Tuple[float, np.ndarray, np.ndarray]:
        center_target = Superpose.get_center(target_coords); center_mobile = Superpose.get_center(mobile_coords); target_centered = target_coords - center_target; mobile_centered = mobile_coords - center_mobile
        H = mobile_centered.T @ target_centered; U, S, Vt = np.linalg.svd(H); R = Vt.T @ U.T
        if np.linalg.det(R) < 0: Vt[-1, :] *= -1; R = Vt.T @ U.T
        mobile_rotated = Superpose.do_rot(mobile_centered, R); num_atoms = target_centered.shape[0]
        if num_atoms == 0: return 0.0, R, center_mobile
        rmsd = np.sqrt(np.sum((target_centered - mobile_rotated) ** 2) / num_atoms)
        return rmsd, R, center_mobile
class PDBReader:
    @staticmethod
    def getxyz(line: str) -> np.ndarray: return np.array([float(line[30+8*i : 38+8*i]) for i in range(3)], dtype=np.float64)
    @staticmethod
    def _read_fragment_internal(filename: str, backbone_only: bool) -> Fragment:
        p = Fragment(); n = 0; rn = 0; allowed_atoms = {"N", "CA", "C", "O"} if backbone_only else None
        try:
            with open(filename, 'r') as f: lines = f.readlines()
            for line in lines:
                if not line.startswith("ATOM"): continue
                if backbone_only and line[13:16].strip() not in allowed_atoms: continue
                if rn == 0 or line[17:26] != p.res[rn - 1].nam:
                    if rn >= len(p.res): break
                    p.res[rn].nam = line[17:26]; p.res[rn].start = n
                    if rn != 0: p.res[rn-1].end = n - 1
                    rn += 1
                if n >= len(p.xyz): break
                p.xyz[n] = PDBReader.getxyz(line); p.nam[n] = line[13:16].strip(); atom_name = line[13:16]; res_idx = rn - 1
                if atom_name == "N  ":
                    if rn == 1: p.N = n
                    p.res[res_idx].N = n; p.CN = n
                elif atom_name == "CA ":
                    if rn == 1: p.CAN = n
                    p.res[res_idx].CA = n; p.CAC = n
                elif atom_name == "C  ":
                    if rn == 1: p.NC = n
                    p.res[res_idx].C = n; p.C = n
                n += 1
            p.an = n; p.resn = rn
            if rn > 0: p.res[rn-1].end = n - 1
        except FileNotFoundError: return None
        return p
    @staticmethod
    def readfrag(filename: str) -> Fragment: return PDBReader._read_fragment_internal(filename, backbone_only=False)
    @staticmethod
    def readfrag_backbone_only(filename: str) -> Fragment: return PDBReader._read_fragment_internal(filename, backbone_only=True)
    @staticmethod
    def readtarget(filename: str) -> np.ndarray:
        target = [];
        try:
            with open(filename, 'r') as f:
                for line in f:
                    if line.startswith("ATOM"): target.append(PDBReader.getxyz(line))
        except FileNotFoundError: return None
        return np.array(target, dtype=np.float64)


class FragmentLinker:
    def __init__(self): self.pdb_reader = PDBReader(); self.patch_cache = {}
    def getgapgeo(self, p1: Fragment, p2: Fragment) -> Tuple:
        terxyz = np.zeros((2, 2, 4, 3)); terxyz[0,0,0,:]=p1.xyz[p1.N]; terxyz[0,0,1,:]=p1.xyz[p1.CAN]; terxyz[0,0,2,:]=p2.xyz[p2.CAC]; terxyz[0,0,3,:]=p2.xyz[p2.C]; terxyz[1,0,0,:]=p2.xyz[p2.N]; terxyz[1,0,1,:]=p2.xyz[p2.CAN]; terxyz[1,0,2,:]=p1.xyz[p1.CAC]; terxyz[1,0,3,:]=p1.xyz[p1.C]; terxyz[0,1,0,:]=p2.xyz[p2.CN]; terxyz[0,1,1,:]=p2.xyz[p2.CAC]; terxyz[0,1,2,:]=p1.xyz[p1.CAN]; terxyz[0,1,3,:]=p1.xyz[p1.NC]; terxyz[1,1,0,:]=p1.xyz[p1.CN]; terxyz[1,1,1,:]=p1.xyz[p1.CAC]; terxyz[1,1,2,:]=p2.xyz[p2.CAN]; terxyz[1,1,3,:]=p2.xyz[p2.NC]
        tmdis  = np.array([Geometry.distance(p1.xyz[p1.CAN], p2.xyz[p2.CAC]), Geometry.distance(p2.xyz[p2.CAN], p1.xyz[p1.CAC])]); ncaca  = np.array([Geometry.cal_angle(p1.xyz[p1.N], p1.xyz[p1.CAN], p2.xyz[p2.CAC]), Geometry.cal_angle(p2.xyz[p2.N], p2.xyz[p2.CAN], p1.xyz[p1.CAC])]); cacac  = np.array([Geometry.cal_angle(p1.xyz[p1.CAN], p2.xyz[p2.CAC], p2.xyz[p2.C]), Geometry.cal_angle(p2.xyz[p2.CAN], p1.xyz[p1.CAC], p1.xyz[p1.C])]); dih    = np.array([Geometry.cal_dih(p1.xyz[p1.N], p1.xyz[p1.CAN], p2.xyz[p2.CAC], p2.xyz[p2.C]), Geometry.cal_dih(p2.xyz[p2.N], p2.xyz[p2.CAN], p1.xyz[p1.CAC], p1.xyz[p1.C])]); ncaca2 = np.array([Geometry.cal_angle(p2.xyz[p2.CN], p2.xyz[p2.CAC], p1.xyz[p1.CAN]), Geometry.cal_angle(p1.xyz[p1.CN], p1.xyz[p1.CAC], p2.xyz[p2.CAN])]); cacac2 = np.array([Geometry.cal_angle(p2.xyz[p2.CAC], p1.xyz[p1.CAN], p1.xyz[p1.NC]), Geometry.cal_angle(p1.xyz[p1.CAC], p2.xyz[p2.CAN], p2.xyz[p2.NC])]); dih2   = np.array([Geometry.cal_dih(p2.xyz[p2.CN], p2.xyz[p2.CAC], p1.xyz[p1.CAN], p1.xyz[p1.NC]), Geometry.cal_dih(p1.xyz[p1.CN], p1.xyz[p1.CAC], p2.xyz[p2.CAN], p2.xyz[p2.NC])])
        return tmdis, ncaca, cacac, dih, ncaca2, cacac2, dih2, terxyz
    def centerterminal(self, terxyz0: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        tcenter = np.mean(terxyz0, axis=2); terxyz = terxyz0 - tcenter[:, :, np.newaxis, :]; return tcenter, terxyz
    def _vectorized_steric_check(self, coords1, coords2, d_max_sq, d_min_sq, penalty_d_min, penalty_d_offset):
        if coords1.shape[0] == 0 or coords2.shape[0] == 0: return 0.0
        dist_sq_matrix = np.sum((coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :]) ** 2, axis=2)
        steric = np.sum(dist_sq_matrix < d_min_sq) * ACCEPTSTERIC * 2
        if steric > ACCEPTSTERIC: return steric
        moderate_mask = (dist_sq_matrix >= d_min_sq) & (dist_sq_matrix < d_max_sq)
        if np.any(moderate_mask): steric += np.sum(1.0 - (np.sqrt(dist_sq_matrix[moderate_mask]) - penalty_d_min) / penalty_d_offset)
        return steric
    def checksteric(self, p: Fragment, p1: Fragment, p2: Fragment) -> Tuple[float, List[int]]:
        steric = 0.0; ter = [1, 1]; coords_p_internal = p.xyz[p.res[1].start : p.res[p.resn - 1].start]
        steric += self._vectorized_steric_check(coords_p_internal, p2.xyz[:p2.res[p2.resn - 1].start], 9.0, 5.76, 2.4, 0.6)
        if steric > ACCEPTSTERIC: return steric, ter
        steric += self._vectorized_steric_check(coords_p_internal, p1.xyz[p1.res[1].start : p1.an], 9.0, 5.76, 2.4, 0.6)
        if steric > ACCEPTSTERIC: return steric, ter
        tersteric = 0.0
        for ia in range(p.res[1].start, p.res[p.resn-1].start):
            for ja in range(4, p1.res[1].start):
                d2 = Geometry.distance2(p.xyz[ia], p1.xyz[ja]);
                if d2 < 4.84: ter[0] = 0; break
                if d2 < 7.84: tersteric += 1.0 - (math.sqrt(d2) - 2.2) / 0.6
                if tersteric > ACCEPTTERSTERIC: ter[0] = 0; break
            if ter[0] == 0: break
        tersteric = 0.0
        for ia in range(p.res[1].start, p.res[p.resn-1].start):
            for ja in range(p2.res[p2.resn-1].start + 4, p2.an):
                d2 = Geometry.distance2(p.xyz[ia], p2.xyz[ja]);
                if d2 < 4.84: ter[1] = 0; break
                if d2 < 7.84: tersteric += 1.0 - (math.sqrt(d2) - 2.2) / 0.6
                if tersteric > ACCEPTTERSTERIC: ter[1] = 0; break
            if ter[1] == 0: break
        return steric, ter
    def checksteric2(self, p: Fragment, p1: Fragment, p2: Fragment) -> Tuple[float, List[int]]:
        steric, ter = 0.0, [1, 1]; coords_p_terminal = p.xyz[p.NC : p.CN + 1]
        if p2.resn > 1:
            coords_p2 = np.delete(p2.xyz[:p2.res[p2.resn-1].start], p2.res[p2.resn-2].C, axis=0)
            steric += self._vectorized_steric_check(coords_p_terminal, coords_p2, 9.0, 5.76, 2.4, 0.6)
            if steric > ACCEPTSTERIC: return steric, ter
        if p1.resn > 1:
            coords_p1 = p1.xyz[p1.res[1].start + 1 : p1.an]
            steric += self._vectorized_steric_check(coords_p_terminal, coords_p1, 9.0, 5.76, 2.4, 0.6)
            if steric > ACCEPTSTERIC: return steric, ter
        if math.sqrt(Geometry.distance2(p.xyz[p.res[p.resn - 1].N], p1.xyz[p1.res[0].N])) >= 0.8: ter[0] = 0
        if math.sqrt(Geometry.distance2(p.xyz[p.res[0].C], p2.xyz[p2.res[p2.resn - 1].C])) >= 0.8: ter[1] = 0
        return steric, ter
    def pepterdock(self, target: np.ndarray, p1: Fragment, p2: Fragment) -> List[float]: return [0.0] * 4
    def dock2target(self, target: np.ndarray, p: Fragment) -> float:
        coords_p = p.xyz[p.res[1].start : p.res[p.resn-1].start]
        steric = self._vectorized_steric_check(coords_p, target, 9.0, 5.76, 2.4, 0.6)
        return min(steric, ACCEPTDOCKSTERIC)

    def readfraglib(self, geo_data: Dict, libname: str, typ: int, tmdis: float, ncaca: float, cacac: float, dih: float) -> List[Dict]:
        """修改后的版本，从内存中的 geo_data 字典读取数据，并修正了路径拼接的bug"""
        patches = []
        d_range = np.arange(tmdis - DISDIFF, tmdis + DISDIFF + 0.0005, 0.001)
        
        for d in d_range:
            # ========================= 关键修正 =========================
            # 先生成纯文件名，然后使用 os.path.join 来确保路径正确
            geo_filename = f"{'D' if typ == 0 else 'd'}{d:.3f}"
            disgeofn = os.path.join(libname, geo_filename)
            # ==========================================================

            lines = geo_data.get(disgeofn) # 从内存字典中获取文件内容
            if not lines: continue
            
            for line in lines:
                try:
                    a1, a2, h = float(line[27:36]), float(line[36:45]), float(line[45:54])
                    if abs(a1 - ncaca) > ANG1DIFF or abs(a2 - cacac) > ANG2DIFF: continue
                    
                    h_diff = h - dih
                    if h_diff < -180: h += 360
                    elif h_diff > 180: h -= 360
                    
                    if abs(h - dih) > DIHDIFF: continue
                    
                    score = ((d - tmdis) / DISDIFF) ** 2 + ((a1 - ncaca) / ANG1DIFF) ** 2 + \
                            ((a2 - cacac) / ANG2DIFF) ** 2 + ((h - dih) / DIHDIFF) ** 2
                    patches.append({'name': line[:17].strip(), 'score': score})
                except (ValueError, IndexError): continue
                
        patches.sort(key=lambda x: x['score'])
        return patches[:BESTPATCHNUM]

    def readpatchxyz(self, dirname: str, pnam: str, typ: int) -> Tuple[Fragment, np.ndarray]:
        if pnam in self.patch_cache: return self.patch_cache[pnam]
        filename = os.path.join(dirname, f"pre_{pnam}")
        p = self.pdb_reader.readfrag_backbone_only(filename)
        if p is None: return None, None
        terxyz = np.zeros((4, 3))
        if typ == 0: terxyz[0]=p.xyz[p.CN]; terxyz[1]=p.xyz[p.CAC]; terxyz[2]=p.xyz[p.CAN]; terxyz[3]=p.xyz[p.NC]
        else: terxyz[0]=p.xyz[p.N]; terxyz[1]=p.xyz[p.CAN]; terxyz[2]=p.xyz[p.CAC]; terxyz[3]=p.xyz[p.C]
        self.patch_cache[pnam] = (p, terxyz)
        return p, terxyz
    def buildcomplex(self, lpcenter: np.ndarray, trans: np.ndarray, rotR: np.ndarray, p: Fragment):
        mobile_coords = p.xyz[:p.an]; mobile_centered = mobile_coords - trans
        mobile_rotated = Superpose.do_rot(mobile_centered, rotR); p.xyz[:p.an] = mobile_rotated + lpcenter
    def getCAxyz(self, p: Fragment) -> np.ndarray:
        num_res = min(p.resn, 5); return np.array([p.xyz[p.res[ir].CA] for ir in range(num_res)])
    def CAcompare(self, resn: int, linkerCA: np.ndarray, N: int, bestCAn: List[int], bestCA: List[np.ndarray]) -> int:
        for ic in range(N):
            if resn != bestCAn[ic]: continue
            min_res = min(resn, len(bestCA[ic]))
            if min_res == 0: continue
            rmsd = np.sqrt(np.mean(np.sum((bestCA[ic][:min_res] - linkerCA[:min_res]) ** 2, axis=1)))
            if rmsd < CULSTERRMSD: return ic + 1
        return N + 1
    def cp2clusterCA(self, resn: int, linkerCA: np.ndarray, bestCA_slice: np.ndarray):
        num_to_copy = min(resn, len(bestCA_slice)); bestCA_slice[:num_to_copy] = linkerCA[:num_to_copy]
    def printpepstructure(self, id: int, frag1name: str, frag2name: str, fragdir: str, pepname: str, rmsd: float, steric: float, docking: float, tersteric: float, score: float, pepfname: str, p: Fragment, patchi: int, typ: int, ter: List[int], fragscore: List[float]):
        try:
            with open(pepfname, 'a') as f:
                f.write(f"TITLE   {patchi+1:2d} {id:2d} {typ+1:2d}\n"); f.write(f"REMARK  dock1 {ter[0]:2d} {fragscore[0]:8.3f} {frag1name}\n"); f.write(f"REMARK  dock2 {ter[1]:2d} {fragscore[1]:8.3f} {frag2name}\n"); f.write(f"REMARK  link {fragdir}/{pepname}\n"); f.write(f"REMARK  score {score:8.3f} {rmsd:8.3f}{steric:8.3f}{docking:8.3f}{tersteric:8.3f}\n")
                atno = 1
                for ir in range(p.resn):
                    loop_end = min(p.res[ir].start + 4, p.an)
                    for i in range(p.res[ir].start, loop_end):
                        f.write(f"ATOM  {atno:>5d}  {p.nam[i]:<4s}{p.res[ir].nam[:3]:>3s} A{p.res[ir].nam[6:]:>4s}   {p.xyz[i,0]:8.3f}{p.xyz[i,1]:8.3f}{p.xyz[i,2]:8.3f}\n"); atno += 1
                    if ir != 0 and ir != p.resn - 1:
                        side_chain_start = p.res[ir].start + 4; side_chain_end = min(p.res[ir].end + 1, p.an)
                        for i in range(side_chain_start, side_chain_end):        
                            f.write(f"ATOM  {atno:>5d}  {p.nam[i]:<4s}{p.res[ir].nam[:3]:>3s} A{p.res[ir].nam[6:]:>4s}   {p.xyz[i,0]:8.3f}{p.xyz[i,1]:8.3f}{p.xyz[i,2]:8.3f}\n"); atno += 1
                f.write("END\n")
        except (IOError, IndexError): pass

# =======================================================================
#  并行化执行与新功能
# =======================================================================

# 1. 新增：预加载所有几何库文件的函数
def preload_geo_data(libs: List[str]) -> Dict:
    logging.info("开始预加载几何库文件到内存...")
    geo_data = {}
    all_files_to_load = []
    for lib_dir in libs:
        if not os.path.isdir(lib_dir):
            logging.warning(f"几何库目录不存在: {lib_dir}")
            continue
        for filename in os.listdir(lib_dir):
            if filename.startswith('D') or filename.startswith('d'):
                all_files_to_load.append(os.path.join(lib_dir, filename))

    # 使用 tqdm 显示预加载进度
    for filepath in tqdm(all_files_to_load, desc="Preloading geo files"):
        try:
            with open(filepath, 'r') as f:
                geo_data[filepath] = f.readlines()
        except Exception as e:
            logging.warning(f"无法读取文件 {filepath}: {e}")

    total_size_mb = sum(sum(len(line) for line in lines) for lines in geo_data.values()) / (1024 * 1024)
    logging.info(f"预加载完成！共加载 {len(geo_data)} 个文件, 总大小: {total_size_mb:.2f} MB。")
    return geo_data

# 全局变量，用于在工作进程中初始化
WORKER_LINKER, WORKER_TARGET, WORKER_LIBS, WORKER_PATCH_DIR, WORKER_GEO_DATA = None, None, None, None, None

def init_worker(linker, target, libs, patch_dir, geo_data):
    """初始化每个工作进程的全局变量"""
    global WORKER_LINKER, WORKER_TARGET, WORKER_LIBS, WORKER_PATCH_DIR, WORKER_GEO_DATA
    WORKER_LINKER, WORKER_TARGET, WORKER_LIBS, WORKER_PATCH_DIR, WORKER_GEO_DATA = linker, target, libs, patch_dir, geo_data

def process_line(task_info):
    """处理单个输入行的函数，这是每个工作进程执行的核心任务"""
    line_num, line, output_prefix = task_info
    linker = copy.deepcopy(WORKER_LINKER)
    target, libs, patch_dir, geo_data = WORKER_TARGET, WORKER_LIBS, WORKER_PATCH_DIR, WORKER_GEO_DATA

    try:
        parts = line.strip().split()
        if len(parts) < 4: return (line_num, "Skipped_NotEnoughParts", 0)
        
        frg1nm, frg2nm, score1, score2 = parts[0], parts[1], parts[2], parts[3]
        fragscore = [float(score1), float(score2)]
        outname = f"{output_prefix}_{line_num}"
        if os.path.exists(outname): os.remove(outname)
        
        p1 = linker.pdb_reader.readfrag(frg1nm)
        p2 = linker.pdb_reader.readfrag(frg2nm)
        if p1 is None or p2 is None: return (line_num, "Skipped_FragReadError", 0)

        tmdis, ncaca, cacac, dih, ncaca2, cacac2, dih2, terxyz0 = linker.getgapgeo(p1, p2)
        tcenter, terxyz = linker.centerterminal(terxyz0)
        pepterdockscore = linker.pepterdock(target, p1, p2)
        
        all_valid_candidates = [[] for _ in range(2)]
        
        for ii in range(2):
            geos = [(tmdis[ii], ncaca[ii], cacac[ii], dih[ii], 0), (tmdis[ii], ncaca2[ii], cacac2[ii], dih2[ii], 1)]
            for current_tmdis, current_ncaca, current_cacac, current_dih, jj in geos:
                lib_configs = [(libs[0], 4.5, 7.5), (libs[1], 5.5, 10.5), (libs[2], 6.5, 13.5)]
                for il, (lib_dir, min_d, max_d) in enumerate(lib_configs):
                    if not (min_d < current_tmdis < max_d): continue
                    
                    patches = linker.readfraglib(geo_data, lib_dir, jj, current_tmdis, current_ncaca, current_cacac, current_dih)
                    
                    for patch_info in patches:
                        patch_name = patch_info['name']
                        p, pterxyz = linker.readpatchxyz(patch_dir, patch_name, jj)
                        if p is None: continue
                        rmsd, rotR, trans = Superpose.fit_terminal(terxyz[ii, jj], pterxyz)
                        if rmsd >= ACCEPTRMSD: continue
                        linker.buildcomplex(tcenter[ii, jj], trans, rotR, p)
                        p_for_check1, p_for_check2 = (p1, p2) if ii == 0 else (p2, p1)
                        if jj == 0: steric, tersteric = linker.checksteric(p, p_for_check1, p_for_check2)
                        else: steric, tersteric = linker.checksteric2(p, p_for_check1, p_for_check2)
                        if steric >= ACCEPTSTERIC: continue
                        docking = linker.dock2target(target, p)
                        if docking >= ACCEPTDOCKING: continue
                        terminalscore = 0.0
                        if ii == 0:
                            if tersteric[0] == 1: terminalscore += pepterdockscore[0]
                            if tersteric[1] == 1: terminalscore += pepterdockscore[3]
                        else:
                            if tersteric[0] == 1: terminalscore += pepterdockscore[2]
                            if tersteric[1] == 1: terminalscore += pepterdockscore[1]
                        score = rmsd * 4 + steric + docking + il * 2 + terminalscore
                        all_valid_candidates[ii].append({'score': score, 'p_object': p, 'rmsd': rmsd, 'steric': steric, 'docking': docking, 'terminalscore': terminalscore, 'tersteric': tersteric, 'name': patch_name, 'ii': ii, 'jj': jj})

        cluster_count = 0
        top_candidates = [[] for _ in range(2)]
        for ii in range(2):
            all_valid_candidates[ii].sort(key=lambda x: x['score'])
            top_candidates[ii] = all_valid_candidates[ii][:BESTNUM]

        bestCAn = [[0] * BESTNUM for _ in range(2)]; bestCA = [[np.zeros((5, 3)) for _ in range(BESTNUM)] for _ in range(2)]
        for ii in range(2):
            if not top_candidates[ii]: continue
            clusterN = 0
            for candidate_info in top_candidates[ii]:
                if clusterN >= BESTCLUSTERNUM: break
                p = candidate_info['p_object']; linkerCA = linker.getCAxyz(p)
                clusterid = linker.CAcompare(p.resn, linkerCA, clusterN, bestCAn[ii], bestCA[ii])
                if clusterid > clusterN:
                    linker.printpepstructure(clusterN + 1, frg1nm, frg2nm, patch_dir, candidate_info['name'], candidate_info['rmsd'], candidate_info['steric'], candidate_info['docking'], candidate_info['terminalscore'], candidate_info['score'], outname, p, candidate_info['ii'], candidate_info['jj'], candidate_info['tersteric'], fragscore)
                    bestCAn[ii][clusterN] = p.resn; linker.cp2clusterCA(p.resn, linkerCA, bestCA[ii][clusterN]); clusterN += 1
            cluster_count += clusterN

        if cluster_count > 0:
            return (line_num, "Success_FileGenerated", cluster_count)
        else:
            return (line_num, "Success_NoLinkersFound", 0)
    except Exception:
        return (line_num, "FAILED_UnhandledException", 0)

def main():
    if len(sys.argv) < 8 or len(sys.argv) > 9:
        print("用法: python fraglink_parallel_preload.py <fraglist> <target> <output_prefix> <lib3> <lib4> <lib5> <patchdir> [num_workers]")
        sys.exit(1)
    
    fraglist_file, target_file, output_prefix, lib3_dir, lib4_dir, lib5_dir, patch_dir = sys.argv[1:8]
    libs = [lib3_dir, lib4_dir, lib5_dir]
    
    cpu_count = os.cpu_count()
    if len(sys.argv) == 9:
        try: num_workers = int(sys.argv[8])
        except ValueError: print(f"错误：无效的进程数 '{sys.argv[8]}'。"); sys.exit(1)
    else: num_workers = cpu_count - 1 if cpu_count > 1 else 1

    log_filename = "/data1/home/renxinyu/sdock/systerm/fcrn/dock/fraglink.log"
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', handlers=[logging.FileHandler(log_filename, mode='w'), logging.StreamHandler(sys.stdout)])

    logging.info(f"脚本启动，将使用 {num_workers} 个进程。")
    logging.info(f"日志文件将保存在: {log_filename}")

    # 2. 调用预加载函数
    geo_data = preload_geo_data(libs)

    logging.info("主进程正在准备共享数据...")
    linker_template = FragmentLinker()
    target = linker_template.pdb_reader.readtarget(target_file)
    if target is None: logging.error(f"错误：无法读取目标文件 {target_file}"); sys.exit(1)
    
    try:
        with open(fraglist_file, 'r') as f: lines = f.readlines()
        tasks = [(i + 1, line, output_prefix) for i, line in enumerate(lines)]
        total_tasks = len(tasks)
        logging.info(f"任务准备完成，共 {total_tasks} 行待处理。")
    except FileNotFoundError:
        logging.error(f"错误：无法打开片段列表文件 {fraglist_file}"); sys.exit(1)

    start_time = time.time(); processed_count = 0
    log_interval = 500 
    logging.info("开始并行处理...")

    # 3. 将预加载的 geo_data 传递给工作进程
    with multiprocessing.Pool(processes=num_workers, initializer=init_worker, initargs=(linker_template, target, libs, patch_dir, geo_data)) as pool:
        results_iterator = pool.imap_unordered(process_line, tasks)

        for i, result in enumerate(results_iterator):
            processed_count = i + 1; line_num, status, cluster_count = result

            if "Success_FileGenerated" in status: logging.info(f"第 {line_num} 行处理成功，生成了 {cluster_count} 个结构。")
            elif "FAILED" in status: logging.warning(f"第 {line_num} 行处理失败，原因: {status}")
            
            if processed_count % log_interval == 0 or processed_count == total_tasks:
                elapsed_time = time.time() - start_time
                speed = processed_count / elapsed_time if elapsed_time > 0 else 0
                eta_seconds = (total_tasks - processed_count) / speed if speed > 0 else 0
                eta_formatted = str(datetime.timedelta(seconds=int(eta_seconds)))
                logging.info(f"进度: {processed_count}/{total_tasks} ({processed_count/total_tasks:.2%}) | 速度: {speed:.2f} 行/秒 | 剩余时间预估: {eta_formatted}")

    end_time = time.time()
    total_duration = str(datetime.timedelta(seconds=int(end_time - start_time)))
    logging.info(f"所有任务处理完成！总耗时: {total_duration}")

if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()