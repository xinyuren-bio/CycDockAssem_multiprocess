import os
import time
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import heapq
import glob
import traceback
import argparse
import logging
from data_structures import Fragment, FullAtomFragment
from geometry_np import distance, angle, dihedral
from superpose_np import fit_terminal
from scoring import check_steric, check_steric2_np, dock_to_target, pepterdock

NUM_ANGLE1 = 0
NUM_ANGLE2 = 0
NUM_DIH = 0
# ==============================================================================
# ===                           1. 配置区 (现在从命令行参数或默认值获取)      ===
# ==============================================================================
# CONFIG 字典将由 argparse 填充或使用默认值
CONFIG = {}

def calculate_gap_geometry(p1, p2):
    required_keys = ['N_N', 'CA_N', 'C_N', 'N_C', 'CA_C', 'C_C']
    for key in required_keys:
        if p1.key_atom_indices.get(key) is None or p2.key_atom_indices.get(key) is None: return None
    geo_params = {}
    p1_coords = {key: p1.coords[idx] for key, idx in p1.key_atom_indices.items()}
    p2_coords = {key: p2.coords[idx] for key, idx in p2.key_atom_indices.items()}
    geo_params['tmdis_0'] = distance(p1_coords['CA_N'], p2_coords['CA_C'])
    geo_params['ncaca_0'] = angle(p1_coords['N_N'], p1_coords['CA_N'], p2_coords['CA_C'])
    geo_params['cacac_0'] = angle(p1_coords['CA_N'], p2_coords['CA_C'], p2_coords['C_C'])
    geo_params['dih_0'] = dihedral(p1_coords['N_N'], p1_coords['CA_N'], p2_coords['CA_C'], p2_coords['C_C'])
    geo_params['ncaca2_0'] = angle(p2_coords['N_C'], p2_coords['CA_C'], p1_coords['CA_N'])
    geo_params['cacac2_0'] = angle(p2_coords['CA_C'], p1_coords['CA_N'], p1_coords['C_N'])
    geo_params['dih2_0'] = dihedral(p2_coords['N_C'], p2_coords['CA_C'], p1_coords['CA_N'], p1_coords['C_N'])
    geo_params['tmdis_1'] = distance(p2_coords['CA_N'], p1_coords['CA_C'])
    geo_params['ncaca_1'] = angle(p2_coords['N_N'], p2_coords['CA_N'], p1_coords['CA_C'])
    geo_params['cacac_1'] = angle(p2_coords['CA_N'], p1_coords['CA_C'], p1_coords['C_C'])
    geo_params['dih_1'] = dihedral(p2_coords['N_N'], p2_coords['CA_N'], p1_coords['CA_C'], p1_coords['C_C'])
    geo_params['ncaca2_1'] = angle(p1_coords['N_C'], p1_coords['CA_C'], p2_coords['CA_N'])
    geo_params['cacac2_1'] = angle(p1_coords['CA_C'], p2_coords['CA_N'], p2_coords['C_N'])
    geo_params['dih2_1'] = dihedral(p1_coords['N_C'], p1_coords['CA_C'], p2_coords['CA_N'], p2_coords['C_N'])
    
    # [PY_DEBUG] Print gap geometry results
    logging.info(f"[PY_DEBUG] GAP_GEOMETRY_0: tmdis={geo_params['tmdis_0']:.4f} ncaca={geo_params['ncaca_0']:.4f} cacac={geo_params['cacac_0']:.4f} dih={geo_params['dih_0']:.4f} ncaca2={geo_params['ncaca2_0']:.4f} cacac2={geo_params['cacac2_0']:.4f} dih2={geo_params['dih2_0']:.4f}")
    logging.info(f"[PY_DEBUG] GAP_GEOMETRY_1: tmdis={geo_params['tmdis_1']:.4f} ncaca={geo_params['ncaca_1']:.4f} cacac={geo_params['cacac_1']:.4f} dih={geo_params['dih_1']:.4f} ncaca2={geo_params['ncaca2_1']:.4f} cacac2={geo_params['cacac2_1']:.4f} dih2={geo_params['dih2_1']:.4f}")
    
    return geo_params

def find_linkers(geo_lib_dir, typ, target_tmdis, target_ncaca, target_cacac, target_dih, config):
    candidate_scores = []
    geo_lib_params = config['GEO_LIB_PARAMS']
    num_pass_dih = 0
    accepted_count = 0  # [PY_DEBUG] Counter
    
    # [PY_DEBUG] Print search parameters
    logging.info(f"[PY_DEBUG] find_linkers: lib={geo_lib_dir} type={typ} tmdis={target_tmdis:.4f} ncaca={target_ncaca:.4f} cacac={target_cacac:.4f} dih={target_dih:.4f}")
    
    for dish in range(int(geo_lib_params['DISDIFF'] * 2000 + 1)):
        d = target_tmdis + dish * 0.001 - geo_lib_params['DISDIFF']
        filename_prefix = 'D' if typ == 0 else 'd'
        geo_filepath = os.path.join(geo_lib_dir, f"{filename_prefix}{d:.3f}")
        if not os.path.exists(geo_filepath): continue
        try:
            with open(geo_filepath, 'r') as gf:
                for line in gf:
                    line = line.strip()
                    if len(line) < 54: continue
                    
                    parts1 = line[27:36].split()
                    if not parts1: continue
                    a1 = float(parts1[0])
                    
                    if abs(a1 - target_ncaca) > geo_lib_params['ANG1DIFF']:
                        continue

                    parts2 = line[36:45].split()
                    if not parts2: continue
                    a2 = float(parts2[0])
                    
                    if abs(a2 - target_cacac) > geo_lib_params['ANG2DIFF']:
                        continue
                    
                    parts3 = line[45:54].split()
                    if not parts3: continue
                    h = float(parts3[0])
                    
                    h_diff = h - target_dih
                    if h_diff < -180: h_diff += 360
                    if h_diff > 180: h_diff -= 360
                    if abs(h_diff) > geo_lib_params['DIHDIFF']:
                        continue
                    num_pass_dih += 1
                    patch_name = line[:17].strip()
                    score = (((a1 - target_ncaca) / geo_lib_params['ANG1DIFF'])**2 +
                             ((a2 - target_cacac) / geo_lib_params['ANG2DIFF'])**2 +
                             ((h_diff) / geo_lib_params['DIHDIFF'])**2)
                    candidate_scores.append((score, patch_name))
                    accepted_count +=1 # [PY_DEBUG] Increment counter
        
        except (ValueError, IndexError):
            continue
    
    logging.info(f"[PY_DEBUG] find_linkers_result: found {accepted_count} patches.")
    best_candidates = heapq.nsmallest(config['BESTPATCHNUM'], candidate_scores)
    return [name for score, name in best_candidates]

def perform_linker_search(geo_params, config):
    all_candidates_with_context = []
    geo_libs = {3: config['TRIPEP_GEO_DIR'], 4: config['TETRAPEP_GEO_DIR'], 5: config['PENTPEP_GEO_DIR']}
    for direction in [0, 1]:
        tmdis = geo_params[f'tmdis_{direction}']
        search_tasks = []
        if config['TERMINAL_DIS_CA_MIN'] < tmdis < config['TERMINAL_DIS_CA_MAX3']: search_tasks.append({'len_type': 3, 'lib_path': geo_libs[3]})
        if config['TERMINAL_DIS_CA_MIN'] < tmdis < config['TERMINAL_DIS_CA_MAX4']: search_tasks.append({'len_type': 4, 'lib_path': geo_libs[4]})
        if config['TERMINAL_DIS_CA_MIN'] < tmdis < config['TERMINAL_DIS_CA_MAX5']: search_tasks.append({'len_type': 5, 'lib_path': geo_libs[5]})
        for task in search_tasks:
            for sup_type in [0, 1]:
                params_key_suffix = f"_{direction}" if sup_type == 0 else f"2_{direction}"
                params = (
                    geo_params[f'tmdis_{direction}'], geo_params[f'ncaca{params_key_suffix}'],
                    geo_params[f'cacac{params_key_suffix}'], geo_params[f'dih{params_key_suffix}']
                )
                found_patches = find_linkers(task['lib_path'], sup_type, *params, config)
                for patch_name in found_patches:
                    all_candidates_with_context.append({
                        "name": patch_name, "direction": direction,
                        "sup_type": sup_type, "len_type": task['len_type']
                    })
    return all_candidates_with_context

def get_terminal_coords_for_superpose(p1, p2, patch, direction, sup_type):
    try:
        if sup_type == 0:
            mob_coords = np.vstack([
                patch.coords[patch.key_atom_indices['N_C']], patch.coords[patch.key_atom_indices['CA_C']],
                patch.coords[patch.key_atom_indices['CA_N']], patch.coords[patch.key_atom_indices['C_N']]
            ])
        else:
            mob_coords = np.vstack([
                patch.coords[patch.key_atom_indices['N_N']], patch.coords[patch.key_atom_indices['CA_N']],
                patch.coords[patch.key_atom_indices['CA_C']], patch.coords[patch.key_atom_indices['C_C']]
            ])
        if direction == 0:
            if sup_type == 0:
                ref_coords = np.vstack([
                    p1.coords[p1.key_atom_indices['N_N']], p1.coords[p1.key_atom_indices['CA_N']],
                    p2.coords[p2.key_atom_indices['CA_C']], p2.coords[p2.key_atom_indices['C_C']]
                ])
            else:
                ref_coords = np.vstack([
                    p2.coords[p2.key_atom_indices['N_C']], p2.coords[p2.key_atom_indices['CA_C']],
                    p1.coords[p1.key_atom_indices['CA_N']], p1.coords[p1.key_atom_indices['C_N']]
                ])
        else:
            if sup_type == 0:
                ref_coords = np.vstack([
                    p2.coords[p2.key_atom_indices['N_N']], p2.coords[p2.key_atom_indices['CA_N']],
                    p1.coords[p1.key_atom_indices['CA_C']], p1.coords[p1.key_atom_indices['C_C']]
                ])
            else:
                ref_coords = np.vstack([
                    p1.coords[p1.key_atom_indices['N_C']], p1.coords[p1.key_atom_indices['CA_C']],
                    p2.coords[p2.key_atom_indices['CA_N']], p2.coords[p2.key_atom_indices['C_N']]
                ])
        return ref_coords, mob_coords
    except (KeyError, TypeError, ValueError, IndexError):
        return None, None

def _calculate_rmsd(coords1, coords2):
    diff = coords1 - coords2
    return np.sqrt(np.mean(np.sum(diff * diff, axis=1)))

def cluster_results(all_best_candidates, config):
    if not all_best_candidates:
        return []

    sorted_candidates = sorted(all_best_candidates, key=lambda x: x['score'])
    
    final_representatives = []
    cluster_reps_ca_coords = []

    for candidate in sorted_candidates:
        current_patch = candidate['patch']
        current_ca_coords = current_patch.get_ca_coords()
        
        if current_ca_coords.size == 0:
            continue

        is_new_cluster = True
        for rep_ca_coords in cluster_reps_ca_coords:
            if current_ca_coords.shape[0] != rep_ca_coords.shape[0]:
                continue
            
            rmsd = _calculate_rmsd(current_ca_coords, rep_ca_coords)
            
            if rmsd < config['CULSTERRMSD']:
                is_new_cluster = False
                break
        
        if is_new_cluster:
            final_representatives.append(candidate)
            cluster_reps_ca_coords.append(current_ca_coords)
        
        if len(final_representatives) >= config['BESTCLUSTERNUM']:
            break
            
    return final_representatives

def process_pair(args):
    pair_index, frag1_path, frag2_path, frag_scores, target_coords_full, config = args
    
    # [PY_DEBUG] Start processing a new pair
    logging.info(f"\n[PY_DEBUG] START_PAIR: index={pair_index} frag1={frag1_path} frag2={frag2_path}")

    p1_geom = Fragment(frag1_path)
    p2_geom = Fragment(frag2_path)
    if not p1_geom.is_valid() or not p2_geom.is_valid(): return []
    
    p1_full = FullAtomFragment(frag1_path)
    p2_full = FullAtomFragment(frag2_path)
    if p1_full.coords is None or p2_full.coords is None: return []

    pepterdockscore = pepterdock(p1_full, p2_full, target_coords_full)
    # [PY_DEBUG] Print pepterdock scores
    logging.info(f"[PY_DEBUG] PEPTERDOCK: score0={pepterdockscore[0]:.4f} score1={pepterdockscore[1]:.4f} score2={pepterdockscore[2]:.4f} score3={pepterdockscore[3]:.4f}")


    geo_params = calculate_gap_geometry(p1_geom, p2_geom)
    if geo_params is None: return []
    
    linker_candidates = perform_linker_search(geo_params, config)
    if not linker_candidates: return []

    top_candidates = {0: [], 1: []}
    candidate_id_counter = 0

    for i, candidate_info in enumerate(linker_candidates):
        # [PY_DEBUG] Start evaluating a candidate
        logging.info(f"[PY_DEBUG] EVAL_CANDIDATE: direction={candidate_info['direction']} sup_type={candidate_info['sup_type']} len_type={candidate_info['len_type']} patch_name={candidate_info['name']}")
        
        patch = Fragment(os.path.join(config['FRAGLIB_DIR'], str('pre_' + candidate_info['name'])))
        print('coords', patch.coords)
        break
            # --- 在这里添加第一条日志 ---
        if not patch.is_valid():
            logging.info(f"[PY_DEBUG] REJECTED (invalid patch): {candidate_info['name']}")
            continue
        
        coords_ref_orig, coords_mob_orig = get_terminal_coords_for_superpose(p1_geom, p2_geom, patch, candidate_info['direction'], candidate_info['sup_type'])

        # --- 在这里添加第二条日志 ---
        if coords_ref_orig is None or coords_mob_orig is None:
            logging.info(f"[PY_DEBUG] REJECTED (missing atoms for superpose): {candidate_info['name']}")
            continue
        
        # coords_ref_orig, coords_mob_orig = get_terminal_coords_for_superpose(p1_geom, p2_geom, patch, candidate_info['direction'], candidate_info['sup_type'])
        # if coords_ref_orig is None or coords_mob_orig is None: continue
        
        ref_centroid = np.mean(coords_ref_orig, axis=0)
        mob_centroid = np.mean(coords_mob_orig, axis=0)
        R, t, rmsd = fit_terminal(coords_ref_orig - ref_centroid, coords_mob_orig - mob_centroid)
        logging.info(f"[PY_DEBUG] EVAL_RMSD: {rmsd:.4f}")
        
        if rmsd < config['ACCEPT_RMSD']:
            patch.coords = np.dot(patch.coords - mob_centroid, R.T) + ref_centroid
            
            direction = candidate_info['direction']
            sup_type = candidate_info['sup_type']
            steric_score = config['ACCEPT_STERIC'] + 1.0
            ter_flags = [0, 0]

            if direction == 0 and sup_type == 0:
                steric_score, ter_flags = check_steric(patch, p1_geom, p2_geom, sup_type, config)
            elif direction == 1 and sup_type == 0:
                steric_score, ter_flags = check_steric(patch, p2_geom, p1_geom, sup_type, config)
            elif direction == 0 and sup_type == 1:
                steric_score, ter_flags = check_steric2_np(patch, p1_geom, p2_geom, config)
            elif direction == 1 and sup_type == 1:
                steric_score, ter_flags = check_steric2_np(patch, p2_geom, p1_geom, config)
            
            logging.info(f"[PY_DEBUG] EVAL_STERIC: steric_score={steric_score:.4f} ter_flags={ter_flags}")

            if steric_score < config['ACCEPT_STERIC']:
                docking_score = dock_to_target(patch, target_coords_full, config)
                logging.info(f"[PY_DEBUG] EVAL_DOCKING: {docking_score:.4f}")

                if docking_score < config['ACCEPT_DOCKING']:
                    terminalscore = 0.0
                    if direction == 0:
                        if ter_flags[0] == 1: terminalscore += pepterdockscore[0]
                        if ter_flags[1] == 1: terminalscore += pepterdockscore[3]
                    else: # direction == 1
                        if ter_flags[0] == 1: terminalscore += pepterdockscore[2]
                        if ter_flags[1] == 1: terminalscore += pepterdockscore[1]
                    logging.info(f"[PY_DEBUG] EVAL_TERMSCORE: {terminalscore:.4f}")

                    linker_len_penalty = (candidate_info['len_type'] - 3) * 2.0
                    final_score = (rmsd * 4.0) + steric_score + docking_score + linker_len_penalty + terminalscore
                    logging.info(f"[PY_DEBUG] EVAL_FINALSCORE: final_score={final_score:.4f} (rmsd*4={(rmsd*4.0):.4f}, steric={steric_score:.4f}, docking={docking_score:.4f}, len_penalty={linker_len_penalty:.4f}, term_score={terminalscore:.4f})")

                    candidate_data = {
                        'name': candidate_info['name'], 'score': final_score, 'rmsd': rmsd,
                        'steric': steric_score, 'docking': docking_score, 'terminalscore': terminalscore,
                        'patch': patch, 'direction': direction, 'sup_type': sup_type, 'ter_flags': ter_flags
                    }

                    item_to_push = (-final_score, candidate_id_counter, candidate_data)
                    candidate_id_counter += 1

                    if len(top_candidates[direction]) < config['BESTNUM']:
                        heapq.heappush(top_candidates[direction], item_to_push)
                    else:
                        heapq.heappushpop(top_candidates[direction], item_to_push)
                        
    all_final_reps = []
    for direction in [0, 1]:
        direction_candidates = [item[2] for item in heapq.nlargest(config['BESTNUM'], top_candidates[direction])]
        if not direction_candidates:
            continue

        reps_for_direction = cluster_results(direction_candidates, config)
        
        for i, rep in enumerate(reps_for_direction):
            rep['cluster_id'] = i + 1
        
        all_final_reps.extend(reps_for_direction)

    if not all_final_reps: return []
    
    all_final_reps.sort(key=lambda x: (x['direction'], x['cluster_id']))
    
    return all_final_reps

def write_fraglink_output(final_candidates, config, job_args):
    if not final_candidates:
        return

    pair_index, frag1_path, frag2_path, frag_scores = job_args[0:4]
    
    os.makedirs(config['OUTPUT_DIR'], exist_ok=True)

    output_filename = os.path.join(
        config['OUTPUT_DIR'], 
        f"{config['OUTPUT_PREFIX']}_{pair_index}"
    )
    
    try:
        with open(output_filename, 'w') as f:
            for candidate in final_candidates:
                # [PY_DEBUG] Log before writing final output
                logging.info(f"[PY_DEBUG] WRITE_OUTPUT: pair_idx={pair_index} direction={candidate['direction'] + 1} cluster_id={candidate['cluster_id']} sup_type={candidate['sup_type'] + 1} score={candidate['score']:.4f} rmsd={candidate['rmsd']:.4f} steric={candidate['steric']:.4f} docking={candidate['docking']:.4f} terminalscore={candidate['terminalscore']:.4f} patch_name={candidate['name']}")
                
                cluster_id = candidate['cluster_id']
                patch = candidate['patch']
                
                f.write(f"TITLE   {candidate['direction'] + 1}  {cluster_id}  {candidate['sup_type'] + 1}\n")

                if candidate['direction'] == 0:
                    dock1_flag, dock1_score, dock1_path = candidate['ter_flags'][0], frag_scores[0], frag1_path
                    dock2_flag, dock2_score, dock2_path = candidate['ter_flags'][1], frag_scores[1], frag2_path
                else:
                    dock1_flag, dock1_score, dock1_path = candidate['ter_flags'][0], frag_scores[1], frag2_path
                    dock2_flag, dock2_score, dock2_path = candidate['ter_flags'][1], frag_scores[0], frag1_path
                
                f.write(f"REMARK  dock1 {dock1_flag:2d} {dock1_score:8.3f} {dock1_path}\n")
                f.write(f"REMARK  dock2 {dock2_flag:2d} {dock2_score:8.3f} {dock2_path}\n")
                f.write(f"REMARK  link {config['FRAGLIB_DIR']}/{candidate['name']}\n")
                f.write(f"REMARK  score   {candidate['score']:8.3f} {candidate['rmsd']:8.3f} {candidate['steric']:8.3f} {candidate['docking']:8.3f} {candidate['terminalscore']:8.3f}\n")
                
                atom_serial = 1
                for res_id in patch.residues:
                    res_info = patch.residue_info[res_id]
                    res_name = res_info['res_name']
                    chain_id = res_id[0] if res_id[0] else 'A'
                    res_seq = res_id[1]
                    for atom_idx in range(res_info['start_index'], res_info['end_index'] + 1):
                        atom_name = patch.atom_names[atom_idx]
                        x, y, z = patch.coords[atom_idx]
                        atom_line = (f"ATOM  {atom_serial:5d}  {atom_name:<4s}{res_name:>3s} {chain_id}{res_seq:4d}    "
                                     f"{x:8.3f}{y:8.3f}{z:8.3f}\n")
                        f.write(atom_line)
                        atom_serial += 1
                f.write("END\n")
    except IOError as e:
        logging.error(f"错误: 无法写入输出文件 '{output_filename}': {e}")
        with open(config['COMBINED_OUTPUT_LOG'], 'a') as f_err:
            f_err.write(f"IOError: 无法写入 {output_filename}: {e}\n")
            traceback.print_exc(file=f_err)
    except Exception as e:
        logging.error(f"错误: 写入文件 '{output_filename}' 时发生未知错误: {e}")
        with open(config['COMBINED_OUTPUT_LOG'], 'a') as f_err:
            f_err.write(f"未知错误: 写入 {output_filename}: {e}\n")
            traceback.print_exc(file=f_err)


def format_time(seconds):
    if seconds is None:
        return "N/A"
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return f"{int(h):02d}:{int(m):02d}:{int(s):02d}"


def main():
    import logging
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
        
    parser = argparse.ArgumentParser(description="Fraglink Python version for finding linkers.")
    parser.add_argument("--pair_list_file", type=str, required=True,
                        help="Path to the file containing pairs of fragment paths and scores.")
    parser.add_argument("--target_pdb_file", type=str, required=True,
                        help="Path to the target PDB file (e.g., box1.pdb).")
    parser.add_argument("--tripep_geo_dir", type=str, required=True,
                        help="Path to the tripeptide geometry directory.")
    parser.add_argument("--tetrapep_geo_dir", type=str, required=True,
                        help="Path to the tetrapeptide geometry directory.")
    parser.add_argument("--pentpep_geo_dir", type=str, required=True,
                        help="Path to the pentapeptide geometry directory.")
    parser.add_argument("--fraglib_dir", type=str, required=True,
                        help="Path to the fragment library directory.")
    parser.add_argument("--output_prefix", type=str, default="fraglink_pair",
                        help="Prefix for output linker files (e.g., 'fraglink_pair').")
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Directory to save output linker files.")
    parser.add_argument("--max_workers", type=int, default=os.cpu_count(),
                        help="Maximum number of parallel processes (default: CPU count).")
    parser.add_argument("--debug_mode", action='store_true',
                        help="Enable debug mode for more verbose output.")
    parser.add_argument("--combined_output_log", type=str, default="run_error_log.txt",
                        help="File to log errors and exceptions.")
    parser.add_argument("--progress_log_file", type=str, default="run_progress.log",
                        help="File to log progress messages.")

    args = parser.parse_args()

    global CONFIG
    CONFIG['DEBUG_MODE'] = args.debug_mode
    CONFIG['MAX_WORKERS'] = args.max_workers
    CONFIG['PAIR_LIST_FILE'] = args.pair_list_file
    CONFIG['TARGET_PDB_FILE'] = args.target_pdb_file
    CONFIG['TRIPEP_GEO_DIR'] = args.tripep_geo_dir
    CONFIG['TETRAPEP_GEO_DIR'] = args.tetrapep_geo_dir
    CONFIG['PENTPEP_GEO_DIR'] = args.pentpep_geo_dir
    CONFIG['FRAGLIB_DIR'] = args.fraglib_dir
    CONFIG['OUTPUT_PREFIX'] = args.output_prefix
    CONFIG['OUTPUT_DIR'] = args.output_dir
    CONFIG['COMBINED_OUTPUT_LOG'] = args.combined_output_log
    CONFIG['PROGRESS_LOG_FILE'] = args.progress_log_file
    
    CONFIG['ACCEPT_RMSD'] = 0.2
    CONFIG['ACCEPT_STERIC'] = 4
    CONFIG['ACCEPT_DOCKING'] = 2
    CONFIG['BESTPATCHNUM'] = 500
    CONFIG['BESTNUM'] = 30
    CONFIG['BESTCLUSTERNUM'] = 4
    CONFIG['CULSTERRMSD'] = 0.6
    CONFIG['GEO_LIB_PARAMS'] = {'DISDIFF': 0.5, 'ANG1DIFF': 15, 'ANG2DIFF': 15, 'DIHDIFF': 35.0} # Changed to match C defines
    CONFIG['TERMINAL_DIS_CA_MIN'] = 4.5
    CONFIG['TERMINAL_DIS_CA_MAX3'] = 7.5
    CONFIG['TERMINAL_DIS_CA_MAX4'] = 10.5
    CONFIG['TERMINAL_DIS_CA_MAX5'] = 13.5

    # In python version, I changed the logging level from INFO to DEBUG to make it more verbose
    log_level = logging.DEBUG if CONFIG['DEBUG_MODE'] else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=[
            logging.FileHandler(CONFIG['PROGRESS_LOG_FILE'], mode='w'),
            logging.StreamHandler()
        ]
    )


    logging.info("--- Fraglink Python版开始执行 ---")
    start_time = time.time()
    
    os.makedirs(CONFIG['OUTPUT_DIR'], exist_ok=True)
    logging.info(f"输出文件将保存在: {CONFIG['OUTPUT_DIR']}")

    logging.info(f"正在加载目标结构: {CONFIG['TARGET_PDB_FILE']}...")
    target_structure_full = FullAtomFragment(CONFIG['TARGET_PDB_FILE'])
    if target_structure_full.coords is None:
        logging.error("错误: 无法加载目标PDB文件，程序终止。")
        return
    target_coords_full = target_structure_full.coords
    logging.info("目标结构加载成功。")

    try:
        with open(CONFIG['PAIR_LIST_FILE'], 'r') as f:
            all_pairs = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        logging.error(f"错误: 找不到配对列表文件: {CONFIG['PAIR_LIST_FILE']}")
        return

    jobs = []
    logging.info(f"从 '{os.path.basename(CONFIG['PAIR_LIST_FILE'])}' 中加载了 {len(all_pairs)} 行，正在解析任务...")
    
    for i, line in enumerate(all_pairs):
        parts = line.split()
        if len(parts) < 4: 
            logging.warning(f"警告: 跳过无效行 (少于4部分): {line.strip()}")
            continue
        try:
            frag1_path = parts[0]; frag2_path = parts[1]
            frag_scores = [float(parts[2]), float(parts[3])]
            
            if not os.path.exists(frag1_path):
                logging.warning(f"警告: 片段文件 '{frag1_path}' 不存在，跳过配对任务 {i+1}。")
                continue
            if not os.path.exists(frag2_path):
                logging.warning(f"警告: 片段文件 '{frag2_path}' 不存在，跳过配对任务 {i+1}。")
                continue

            jobs.append((i + 1, frag1_path, frag2_path, frag_scores, target_coords_full, CONFIG))
        except ValueError as ve:
            logging.warning(f"警告: 解析行中的分数时出错，跳过行 {i+1}: {line.strip()} - {ve}")
            continue
        except Exception as e:
            logging.error(f"错误: 解析配对列表行时发生未知错误，跳过行 {i+1}: {line.strip()} - {e}")
            continue
    
    total_jobs = len(jobs)
    if total_jobs == 0:
        logging.error("错误: 未能从输入文件中解析出任何有效的配对任务。请检查PAIR_LIST_FILE内容。")
        return

    logging.info(f"成功解析了 {total_jobs} 个配对任务。")
    logging.info(f"使用最多 {CONFIG['MAX_WORKERS'] or os.cpu_count()} 个并行进程。")
    logging.info("-" * 40)

    job_count = 0
    try:
        with open(CONFIG['COMBINED_OUTPUT_LOG'], 'w') as f_err, ProcessPoolExecutor(max_workers=CONFIG['MAX_WORKERS']) as executor:
            future_to_job = {executor.submit(process_pair, job): job for job in jobs}
            
            for future in as_completed(future_to_job):
                job_args = future_to_job[future]
                pair_index = job_args[0]
                job_count += 1
                try:
                    final_candidates_from_pair = future.result()
                    if final_candidates_from_pair:
                        logging.info(f"进度: {job_count}/{total_jobs} | Pair {pair_index} 成功，找到 {len(final_candidates_from_pair)} 个结果，正在写入文件...")
                        write_fraglink_output(final_candidates_from_pair, CONFIG, job_args)
                    else:
                        logging.info(f"进度: {job_count}/{total_jobs} | Pair {pair_index} 完成，未找到符合要求的结果。")

                except Exception:
                    logging.error(f"!!!!!! 错误: 处理配对任务 {pair_index} 时发生严重错误 !!!!!!")
                    traceback.print_exc(file=f_err)
                
                elapsed_time = time.time() - start_time
                avg_time_per_job = elapsed_time / job_count
                jobs_remaining = total_jobs - job_count
                eta_seconds = avg_time_per_job * jobs_remaining
                
                eta_str = format_time(eta_seconds)
                elapsed_str = format_time(elapsed_time)
                
                logging.info(f"--- [进度: {job_count}/{total_jobs} | 已耗时: {elapsed_str} | 预计剩余: {eta_str}] ---")

    except Exception as e:
        logging.error(f"\n主进程发生严重错误: {e}")
        with open(CONFIG['COMBINED_OUTPUT_LOG'], 'a') as f_err:
            f_err.write(f"主进程错误: {e}\n")
            traceback.print_exc(file=f_err)
        return

    total_duration_str = format_time(time.time() - start_time)
    logging.info("-" * 40)
    logging.info(f"--- 所有任务执行完毕！总耗时: {total_duration_str} ---")

if __name__ == "__main__":
    main()
