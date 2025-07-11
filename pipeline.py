import os
import subprocess
import argparse
import yaml
import sys
import re
import glob
import numpy as np

def _run_command(command_parts, cwd=None, step_name="Unknown Step"):
    """
    Helper function to execute a shell command and handle its output.
    """
    full_command = ' '.join(command_parts)
    print(f"\n--- Running Command for {step_name} ---")
    print(f"Command: {full_command}")
    print(f"Working Directory: {os.path.abspath(cwd) if cwd else os.getcwd()}")

    try:
        process = subprocess.run(
            command_parts,
            cwd=cwd,
            check=True,
            capture_output=True,
            text=True,
            encoding='utf-8'
        )
        print(f"✅ {step_name} Succeeded.")
        if process.stdout:
            print(f"   Stdout:\n{process.stdout.strip()}")
        if process.stderr:
            print(f"   Stderr:\n{process.stderr.strip()}")
        return True, process.stdout.strip()
    except FileNotFoundError:
        print(f"🔴 {step_name} Failed: Command not found. Check executable path.")
        return False, f"Command not found: {command_parts[0]}"
    except subprocess.CalledProcessError as e:
        print(f"🔴 {step_name} Failed with error code {e.returncode}.")
        print(f"   Command: {e.cmd}")
        print(f"   Stdout:\n{e.stdout.strip()}")
        print(f"   Stderr:\n{e.stderr.strip()}")
        return False, f"Command failed: {e.cmd}\nStdout: {e.stdout.strip()}\nStderr: {e.stderr.strip()}"
    except Exception as e:
        print(f"🔴 {step_name} Failed with unexpected error: {e}")
        return False, f"Unexpected error: {e}"

def _create_dir(path, step_name):
    """Helper to create a directory."""
    try:
        os.makedirs(path, exist_ok=True)
        print(f"📁 {step_name}: Directory created/exists: {path}")
        return True
    except Exception as e:
        print(f"🔴 {step_name}: Failed to create directory {path}: {e}")
        return False

def _resolve_path(base_dir, path_template, replacements):
    """
    Resolves a path template by replacing placeholders and joining with base_dir.
    `replacements` is a dictionary like {'target_protein_name': 'ALK1', 'sdock_base_dir': 'SDOCK2.0-restrict/'}
    """
    if isinstance(path_template, str):
        resolved_path = path_template
        for key, value in replacements.items():
            if isinstance(value, str):
                resolved_path = resolved_path.replace(f"{{{{{key}}}}}", value)
        
        if os.path.isabs(resolved_path):
            return resolved_path
        else:
            return os.path.join(base_dir, resolved_path)
    return path_template


def run_step1_download_pdb(config, base_dir, replacements):
    step_name = "Step 1: Download PDB"
    if not config.get('enabled', True): return True, None
    
    script_path = _resolve_path(base_dir, config['script_path'], replacements)
    input_file = _resolve_path(base_dir, config['input_file'], replacements)
    output_dir = _resolve_path(base_dir, config['output_dir'], replacements)
    
    if not _create_dir(output_dir, step_name): return False, None

    command = [
        sys.executable, script_path,
        "--input_file", input_file,
        "--output_dir", output_dir
    ]
    cwd = os.path.dirname(script_path) 
    success, _ = _run_command(command, cwd=cwd, step_name=step_name)
    return success, None

def run_step2_pdb2fragments(config, base_dir, replacements):
    step_name = "Step 2: PDB to Fragments"
    if not config.get('enabled', True): return True, None

    script_path = _resolve_path(base_dir, config['script_path'], replacements)
    
    cwd = os.path.dirname(script_path)
    command = [sys.executable, script_path]
    
    success, _ = _run_command(command, cwd=cwd, step_name=step_name)
    return success, None

def run_step3_rotate_protein(config, base_dir, replacements):
    step_name = "Step 3: Rotate Protein"
    if not config.get('enabled', True): return True, None

    script_path = _resolve_path(base_dir, config['script_path'], replacements)
    input_pdb = _resolve_path(base_dir, config['input_pdb'], replacements)
    output_pdb = _resolve_path(base_dir, config['output_pdb'], replacements)

    command = [
        sys.executable, script_path,
        "--input_pdb", input_pdb,
        "--output", output_pdb,
        "--ref1", config['ref1'],
        "--ref2", config['ref2'],
        "--ref3", config['ref3']
    ]
    cwd = os.path.dirname(script_path)
    success, stdout_output = _run_command(command, cwd=cwd, step_name=step_name)
    
    box1_coords = None
    box2_coords = None

    if success and stdout_output:
        for line in stdout_output.splitlines():
            if line.strip().startswith("box1"):
                parts = line.strip().split()
                if len(parts) >= 4:
                    box1_coords = f"{parts[1]},{parts[2]},{parts[3]}"
            elif line.strip().startswith("box2"):
                parts = line.strip().split()
                if len(parts) >= 4:
                    box2_coords = f"{parts[1]},{parts[2]},{parts[3]}"
        
        if box1_coords and box2_coords:
            print(f"Extracted Box1 coords: {box1_coords}")
            print(f"Extracted Box2 coords: {box2_coords}")
            return True, {'box1_m_param': box1_coords, 'box2_m_param': box2_coords}
        else:
            print(f"🔴 {step_name} Failed: Could not extract Box1 or Box2 coordinates from output. Check rotate.py output format.")
            return False, None
    
    return success, None

def run_step4_generate_boxes_watmaps(config, base_dir, replacements, dynamic_params=None):
    step_name = "Step 4: Generate Boxes and Watmaps"
    if not config.get('enabled', True): return True, None

    sdock_restrict_dir = _resolve_path(base_dir, config['sdock_restrict_dir'], replacements)
    preprocess_exec = os.path.join(sdock_restrict_dir, "preprocess")
    watmap_exec = os.path.join(sdock_restrict_dir, "watmap")
    rotated_pdb = _resolve_path(base_dir, config['rotated_pdb'], replacements)
    
    box1_m_param = dynamic_params.get('box1_m_param') if dynamic_params and 'box1_m_param' in dynamic_params else config['box1_m_param']
    box2_m_param = dynamic_params.get('box2_m_param') if dynamic_params and 'box2_m_param' in dynamic_params else config['box2_m_param']
    print(f"Using Box1 M param: {box1_m_param}")
    print(f"Using Box2 M param: {box2_m_param}")
    if not box1_m_param or not box2_m_param:
        print(f"🔴 {step_name} Failed: Box M parameters are missing. Ensure Step 3 ran successfully or they are defined in config.")
        return False, None

    # Preprocess Box 1
    command_box1_preprocess = [
        preprocess_exec, rotated_pdb,
        "-o", _resolve_path(base_dir, config['box1_output_pdb'], replacements),
        "-a", os.path.join(sdock_restrict_dir, "ATM"),
        "-m", box1_m_param
    ]
    if not _run_command(command_box1_preprocess, step_name=f"{step_name} (Box1 Preprocess)")[0]: return False, None

    # Preprocess Box 2
    command_box2_preprocess = [
        preprocess_exec, rotated_pdb,
        "-o", _resolve_path(base_dir, config['box2_output_pdb'], replacements),
        "-a", os.path.join(sdock_restrict_dir, "ATM"),
        "-m", box2_m_param
    ]
    if not _run_command(command_box2_preprocess, step_name=f"{step_name} (Box2 Preprocess)")[0]: return False, None

    # Watmap Box 1
    command_box1_watmap = [
        watmap_exec,
        _resolve_path(base_dir, config['box1_output_pdb'], replacements),
        _resolve_path(base_dir, config['box1_watmap_output'], replacements)
    ]
    if not _run_command(command_box1_watmap, step_name=f"{step_name} (Box1 Watmap)")[0]: return False, None

    # Watmap Box 2
    command_box2_watmap = [
        watmap_exec,
        _resolve_path(base_dir, config['box2_output_pdb'], replacements),
        _resolve_path(base_dir, config['box2_watmap_output'], replacements)
    ]
    if not _run_command(command_box2_watmap, step_name=f"{step_name} (Box2 Watmap)")[0]: return False, None
    
    return True, None

def run_step5_process_fragments(config, base_dir, replacements):
    step_name = "Step 5: Process Fragments (Pre-processed)"
    if not config.get('enabled', True): return True, None

    preprocessed_fragments_path = replacements.get('preprocessed_fragments_path')
    water_fragments_path = replacements.get('water_fragments_path')

    if not preprocessed_fragments_path:
        print(f"🔴 {step_name}: Pre-processed fragments path not found in global configuration.")
        return False, None
    if not water_fragments_path:
        print(f"🔴 {step_name}: Water fragments path not found in global configuration.")
        return False, None

    if not os.path.isdir(preprocessed_fragments_path):
        print(f"🔴 {step_name}: Pre-processed fragments directory not found: {preprocessed_fragments_path}")
        return False, None
    if not os.path.isdir(water_fragments_path):
        print(f"🔴 {step_name}: Water-processed fragments directory not found: {water_fragments_path}")
        return False, None
    
    print(f"✅ {step_name}: Found pre-processed fragments at {preprocessed_fragments_path}")
    print(f"✅ {step_name}: Found water-processed fragments at {water_fragments_path}")
    
    return True, {
        'preprocessed_fragments_path': preprocessed_fragments_path,
        'water_fragments_path': water_fragments_path
    }

def run_step6_generate_frag_files(config, base_dir, replacements, dynamic_params=None):
    step_name = "Step 6: Generate Frag Files"
    if not config.get('enabled', True): return True, None

    tools_dir = _resolve_path(base_dir, config['tools_dir'], replacements)

    fraglib_input_dir_for_fragstruct = None
    if dynamic_params and 'preprocessed_fragments_path' in dynamic_params:
        fraglib_input_dir_for_fragstruct = dynamic_params['preprocessed_fragments_path']
    else:
        print(f"🔴 {step_name} Failed: Pre-processed fragments path not provided by Step 5 output.")
        return False, None

    # fragstruct.py
    fragstruct_script = os.path.join(tools_dir, "fragstruct.py")
    fragstruct_output_txt = _resolve_path(base_dir, config['fragstruct_output_txt'], replacements)
    command_fragstruct = [
        sys.executable, fragstruct_script,
        "--input_dir", fraglib_input_dir_for_fragstruct,
        "--output_txt", fragstruct_output_txt,
        "--tail", config['fragstruct_tail']
    ]
    if not _run_command(command_fragstruct, cwd=tools_dir, step_name=f"{step_name} (fragstruct.py)")[0]: return False, None

    # struct2wmap.py
    struct2wmap_script = os.path.join(tools_dir, "struct2wmap.py")
    fragwmap_output_file = _resolve_path(base_dir, config['fragwmap_output_file'], replacements)
    command_struct2wmap = [
        sys.executable, struct2wmap_script,
        "--input_file", fragstruct_output_txt,
        "--output_file", fragwmap_output_file,
        "--old_string", config['struct2wmap_old_string'],
        "--new_string", config['struct2wmap_new_string']
    ]
    if not _run_command(command_struct2wmap, cwd=tools_dir, step_name=f"{step_name} (struct2wmap.py)")[0]: return False, None

    # mkdir dockresult
    dockresult_dir = _resolve_path(base_dir, config['dockresult_dir'], replacements)
    if not _create_dir(dockresult_dir, step_name): return False, None

    # box1recordf.py
    recordf_script = os.path.join(tools_dir, "recordf.py")
    box1recordf_output_file = _resolve_path(base_dir, config['box1recordf_output_file'], replacements)
    box1recordf_prefix_string = _resolve_path(base_dir, config['box1recordf_prefix_string'], replacements)
    command_box1recordf = [
        sys.executable, recordf_script,
        "--input_file", fragstruct_output_txt,
        "--output_file", box1recordf_output_file,
        "--prefix_string", box1recordf_prefix_string
    ]
    if not _run_command(command_box1recordf, cwd=tools_dir, step_name=f"{step_name} (box1recordf.py)")[0]: return False, None

    # box2recordf.py
    box2recordf_output_file = _resolve_path(base_dir, config['box2recordf_output_file'], replacements)
    box2recordf_prefix_string = _resolve_path(base_dir, config['box2recordf_prefix_string'], replacements)
    command_box2recordf = [
        sys.executable, recordf_script,
        "--input_file", fragstruct_output_txt,
        "--output_file", box2recordf_output_file,
        "--prefix_string", box2recordf_prefix_string
    ]
    if not _run_command(command_box2recordf, cwd=tools_dir, step_name=f"{step_name} (box2recordf.py)")[0]: return False, None
    
    return True, None

def run_step7_perform_docking(config, base_dir, replacements):
    step_name = "Step 7: Perform Docking"
    if not config.get('enabled', True): return True, None

    tools_dir = _resolve_path(base_dir, config['tools_dir'], replacements)
    docking_base_dir = _resolve_path(base_dir, config['docking_base_dir'], replacements)
    split_n_parts = config['split_n_parts']
    dock_num_tasks = config['dock_num_tasks']
    sdock_executable = _resolve_path(base_dir, config['sdock_executable'], replacements)
    
    resolved_sdock_common_params_str = config['sdock_common_params']
    for key, value in replacements.items():
        if isinstance(value, str):
            resolved_sdock_common_params_str = resolved_sdock_common_params_str.replace(f"{{{{{key}}}}}", value)
    sdock_common_params = resolved_sdock_common_params_str.split()

    # Create directories for split files
    fragstruct_split_dir = os.path.join(docking_base_dir, 'fragstruct_')
    fragwmap_split_dir = os.path.join(docking_base_dir, 'fragwmap_')
    recordf_split_dir = os.path.join(docking_base_dir, 'recordf_')

    if not _create_dir(fragstruct_split_dir, step_name): return False, None
    if not _create_dir(fragwmap_split_dir, step_name): return False, None
    if not _create_dir(recordf_split_dir, step_name): return False, None

    split_script = os.path.join(tools_dir, "splitfragstruct.py")
    
    # Split fragstruct
    fragstruct_input = _resolve_path(base_dir, "{{target_protein_name}}/fragstruct", replacements)
    command_split_fragstruct = [
        sys.executable, split_script,
        "--input_file", fragstruct_input,
        "--output_prefix", os.path.join(fragstruct_split_dir, 'fragstruct'),
        "--n_parts", str(split_n_parts)
    ]
    if not _run_command(command_split_fragstruct, cwd=tools_dir, step_name=f"{step_name} (Split fragstruct)")[0]: return False, None

    # Split fragwmap
    fragwmap_input = _resolve_path(base_dir, "{{target_protein_name}}/fragwmap", replacements)
    command_split_fragwmap = [
        sys.executable, split_script,
        "--input_file", fragwmap_input,
        "--output_prefix", os.path.join(fragwmap_split_dir, 'fragwmap'),
        "--n_parts", str(split_n_parts)
    ]
    if not _run_command(command_split_fragwmap, cwd=tools_dir, step_name=f"{step_name} (Split fragwmap)")[0]: return False, None

    # Split box1recordf
    box1recordf_input = _resolve_path(base_dir, "{{target_protein_name}}/box1recordf", replacements)
    command_split_box1recordf = [
        sys.executable, split_script,
        "--input_file", box1recordf_input,
        "--output_prefix", os.path.join(recordf_split_dir, 'box1recordf'),
        "--n_parts", str(split_n_parts)
    ]
    if not _run_command(command_split_box1recordf, cwd=tools_dir, step_name=f"{step_name} (Split box1recordf)")[0]: return False, None

    # Split box2recordf
    box2recordf_input = _resolve_path(base_dir, "{{target_protein_name}}/box2recordf", replacements)
    command_split_box2recordf = [
        sys.executable, split_script,
        "--input_file", box2recordf_input,
        "--output_prefix", os.path.join(recordf_split_dir, 'box2recordf'),
        "--n_parts", str(split_n_parts)
    ]
    if not _run_command(command_split_box2recordf, cwd=tools_dir, step_name=f"{step_name} (Split box2recordf)")[0]: return False, None

    dock_script = os.path.join(tools_dir, "dock.py")

    # Dock box1
    command_dock_box1 = [
        sys.executable, dock_script,
        "--num_tasks", str(dock_num_tasks),
        "--fixed_pdb", _resolve_path(base_dir, config['box1_fixed_pdb'], replacements),
        "--fragstruct_dir", fragstruct_split_dir,
        "--fixed_wat_pdb", _resolve_path(base_dir, config['box1_fixed_wat_pdb'], replacements),
        "--fragwatmap_dir", fragwmap_split_dir,
        "--output_record_dir", recordf_split_dir,
        "--output_prefix", config['box1_output_prefix'],
        "--sdock_exec", sdock_executable,
        "--common_params", ' '.join(sdock_common_params)
    ]
    if not _run_command(command_dock_box1, cwd=tools_dir, step_name=f"{step_name} (Dock Box1)")[0]: return False, None

    # Dock box2
    command_dock_box2 = [
        sys.executable, dock_script,
        "--num_tasks", str(dock_num_tasks),
        "--fixed_pdb", _resolve_path(base_dir, config['box2_fixed_pdb'], replacements),
        "--fragstruct_dir", fragstruct_split_dir,
        "--fixed_wat_pdb", _resolve_path(base_dir, config['box2_fixed_wat_pdb'], replacements),
        "--fragwatmap_dir", fragwmap_split_dir,
        "--output_record_dir", recordf_split_dir,
        "--output_prefix", config['box2_output_prefix'],
        "--sdock_exec", sdock_executable,
        "--common_params", ' '.join(sdock_common_params)
    ]
    if not _run_command(command_dock_box2, cwd=tools_dir, step_name=f"{step_name} (Dock Box2)")[0]: return False, None

    return True, None

def run_step8_select_top_n_per_box(config, base_dir, replacements):
    step_name = "Step 8: Select Top N per Box"
    if not config.get('enabled', True): return True, None

    script_path = _resolve_path(base_dir, config['script_path'], replacements)
    top_n = config['top_n']

    # Box1
    box1_recordf = _resolve_path(base_dir, config['box1_recordf'], replacements)
    box1_output_path = _resolve_path(base_dir, config['box1_output_path'], replacements)
    command_box1 = [
        sys.executable, script_path,
        "--recordf", box1_recordf,
        "--output_path", box1_output_path,
        "--top_n", str(top_n)
    ]
    if not _run_command(command_box1, cwd=os.path.dirname(script_path), step_name=f"{step_name} (Box1)")[0]: return False, None

    # Box2
    box2_recordf = _resolve_path(base_dir, config['box2_recordf'], replacements)
    box2_output_path = _resolve_path(base_dir, config['box2_output_path'], replacements)
    command_box2 = [
        sys.executable, script_path,
        "--recordf", box2_recordf,
        "--output_path", box2_output_path,
        "--top_n", str(top_n)
    ]
    if not _run_command(command_box2, cwd=os.path.dirname(script_path), step_name=f"{step_name} (Box2)")[0]: return False, None
    
    return True, None

def run_step9_generate_complex_models(config, base_dir, replacements, dynamic_params=None):
    step_name = "Step 9: Generate Complex Models"
    if not config.get('enabled', True): return True, None

    utility_dir = _resolve_path(base_dir, config['utility_dir'], replacements)
    sdock_restrict_dir = _resolve_path(base_dir, config['sdock_restrict_dir'], replacements)
    genbuild_exec = os.path.join(utility_dir, "genBuildcommand")
    build_exec = os.path.join(sdock_restrict_dir, "build")
    dockmodel_dir = _resolve_path(base_dir, config['dockmodel_dir'], replacements)
    so3layer_qua = _resolve_path(base_dir, config['so3layer_qua'], replacements)

    if not _create_dir(dockmodel_dir, step_name): return False, None

    box1_m_param_step9 = "0,0,0"

    calculated_box2_m_param_step9 = config['box2_m_param']

    if dynamic_params and 'box1_m_param' in dynamic_params and 'box2_m_param' in dynamic_params:
        try:
            box1_coords_from_step3 = np.array([float(x) for x in dynamic_params['box1_m_param'].split(',')])
            box2_coords_from_step3 = np.array([float(x) for x in dynamic_params['box2_m_param'].split(',')])

            diff_coords = box1_coords_from_step3 - box2_coords_from_step3
            calculated_box2_m_param_step9 = f"{diff_coords[0]:.3f},{diff_coords[1]:.3f},{diff_coords[2]:.3f}"
            print(f"Computed Box2 M param for Step 9 (Box1_rot - Box2_rot): {calculated_box2_m_param_step9}")
        except Exception as e:
            print(f"🔴 {step_name}: Error calculating dynamic box2_m_param from Step 3 output: {e}")
            print(f"    Using default from config for Step 9: {config['box2_m_param']}")
            calculated_box2_m_param_step9 = config['box2_m_param']

    fragpair_output_file = _resolve_path(base_dir, config['fragpair_output_file'], replacements)
    # getfragpair_c_exec_path 变量已不再需要，因为 Tools/getfragpair.py 内部已处理其路径

    # genBuildcommand Box1
    box1_bestdocking = _resolve_path(base_dir, config['box1_bestdocking'], replacements)
    box1_gen_script = _resolve_path(base_dir, config['box1_gen_script'], replacements)
    box1_fixed_pdb = _resolve_path(base_dir, config['box1_fixed_pdb'], replacements)
    command_genbuild_box1 = [
        genbuild_exec, box1_bestdocking, box1_gen_script,
        build_exec, box1_fixed_pdb, dockmodel_dir,
        so3layer_qua, box1_m_param_step9
    ]
    if not _run_command(command_genbuild_box1, step_name=f"{step_name} (genBuildcommand Box1)")[0]: return False, None

    # genBuildcommand Box2
    box2_bestdocking = _resolve_path(base_dir, config['box2_bestdocking'], replacements)
    box2_gen_script = _resolve_path(base_dir, config['box2_gen_script'], replacements)
    box2_fixed_pdb = _resolve_path(base_dir, config['box2_fixed_pdb'], replacements)
    command_genbuild_box2 = [
        genbuild_exec, box2_bestdocking, box2_gen_script,
        build_exec, box2_fixed_pdb, dockmodel_dir,
        so3layer_qua, calculated_box2_m_param_step9
    ]
    if not _run_command(command_genbuild_box2, step_name=f"{step_name} (genBuildcommand Box2)")[0]: return False, None

    # Execute generated shell scripts
    if not _run_command(["bash", box1_gen_script], step_name=f"{step_name} (Execute genbox1frag.sh)")[0]: return False, None
    if not _run_command(["bash", box2_gen_script], step_name=f"{step_name} (Execute genbox2frag.sh)")[0]: return False, None

    # ls commands
    box1dockfrag_output = _resolve_path(base_dir, "{{target_protein_name}}/box1dockfrag", replacements)
    box2dockfrag_output = _resolve_path(base_dir, "{{target_protein_name}}/box2dockfrag", replacements)
    
    try:
        box1_files = glob.glob(os.path.join(dockmodel_dir, f"SDOCK_{replacements['target_protein_name']}box1*"))
        with open(box1dockfrag_output, 'w') as f:
            for file_path in sorted(box1_files):
                f.write(file_path + '\n')
        print(f"✅ {step_name}: ls output for box1 saved to {box1dockfrag_output}")

        box2_files = glob.glob(os.path.join(dockmodel_dir, f"SDOCK_{replacements['target_protein_name']}box2*"))
        with open(box2dockfrag_output, 'w') as f:
            for file_path in sorted(box2_files):
                f.write(file_path + '\n')
        print(f"✅ {step_name}: ls output for box2 saved to {box2dockfrag_output}")

        if not box1_files:
            print(f"⚠️ {step_name}: Warning: No files found matching 'SDOCK_{replacements['target_protein_name']}box1*' in '{dockmodel_dir}' for box1dockfrag.")
        if not box2_files:
            print(f"⚠️ {step_name}: Warning: No files found matching 'SDOCK_{replacements['target_protein_name']}box2*' in '{dockmodel_dir}' for box2dockfrag.")

    except Exception as e:
        print(f"🔴 {step_name}: Error during file listing for dockfrag: {e}")
        return False, None

    # --- MODIFIED: Call Python wrapper for getfragpair ---
    getfragpair_python_script = _resolve_path(base_dir, config['getfragpair_exec'], replacements) # Points to Tools/getfragpair.py

    command_getfragpair_python = [
        sys.executable, getfragpair_python_script,
        box1dockfrag_output, # Positional arg 1
        box2dockfrag_output, # Positional arg 2
        config['fragpair_score_arg'], # Positional arg 3
        config['fragpair_length_arg'], # Positional arg 4
        fragpair_output_file # Positional arg 5 (output file)
        # --c_exec_path 参数不再由 pipeline.py 传递，由 Tools/getfragpair.py 内部处理
    ]
    
    # Execute the Python wrapper script
    if not _run_command(command_getfragpair_python, cwd=os.path.dirname(getfragpair_python_script), step_name=f"{step_name} (getfragpair Python wrapper)")[0]: return False, None

    return True, None

def run_step10_find_linker(config, base_dir, replacements):
    step_name = "Step 10: Find Linker"
    if not config.get('enabled', True): return True, None

    tools_dir = _resolve_path(base_dir, config['tools_dir'], replacements)
    fraglinking_dir = _resolve_path(base_dir, config['fraglinking_dir'], replacements)
    fraglink_script = os.path.join(tools_dir, config['fraglink_script'])

    if not _create_dir(fraglinking_dir, step_name): return False, None

    # fraglink.py 的命令行参数
    pair_list_file = _resolve_path(base_dir, config['pair_list_file'], replacements)
    target_pdb_file = _resolve_path(base_dir, config['target_pdb_file'], replacements)
    tripep_geo_dir = _resolve_path(base_dir, config['tripep_geo_dir'], replacements)
    tetrapep_geo_dir = _resolve_path(base_dir, config['tetrapep_geo_dir'], replacements)
    pentpep_geo_dir = _resolve_path(base_dir, config['pentpep_geo_dir'], replacements)
    fraglib_dir = _resolve_path(base_dir, config['fraglib_dir'], replacements)
    output_prefix = config['output_prefix'].replace("{{target_protein_name}}", replacements['target_protein_name'])
    combined_output_log = _resolve_path(base_dir, config['combined_output_log'], replacements)
    progress_log_file = _resolve_path(base_dir, config['progress_log_file'], replacements)

    command_fraglink = [
        sys.executable, fraglink_script,
        "--pair_list_file", pair_list_file,
        "--target_pdb_file", target_pdb_file,
        "--tripep_geo_dir", tripep_geo_dir,
        "--tetrapep_geo_dir", tetrapep_geo_dir,
        "--pentpep_geo_dir", pentpep_geo_dir,
        "--fraglib_dir", fraglib_dir,
        "--output_prefix", output_prefix,
        "--output_dir", fraglinking_dir,
        "--max_workers", str(config.get('max_workers', os.cpu_count())),
        "--combined_output_log", combined_output_log,
        "--progress_log_file", progress_log_file
    ]
    if config.get('debug_mode', False):
        command_fraglink.append("--debug_mode")
    
    if not _run_command(command_fraglink, cwd=tools_dir, step_name=f"{step_name} (fraglink.py)")[0]: return False, None
    
    return True, None

def run_step11_assemble_cyc(config, base_dir, replacements):
    step_name = "Step 11: Assemble Cyc"
    if not config.get('enabled', True): return True, None

    tools_dir = _resolve_path(base_dir, config['tools_dir'], replacements)
    utility_dir = _resolve_path(base_dir, config['utility_dir'], replacements)
    assembled_cyc_dir = _resolve_path(base_dir, config['assembled_cyc_dir'], replacements)

    if not _create_dir(assembled_cyc_dir, step_name): return False, None

    # assemblecyc.py
    assemblecyc_script = os.path.join(tools_dir, config['assemblecyc_script'])
    assemblecyc_fraglinking_dir = _resolve_path(base_dir, config['assemblecyc_fraglinking_dir'], replacements)
    assemblecyc_output_dir = _resolve_path(base_dir, config['assemblecyc_output_dir'], replacements)
    assemblecyc_input_prefix = config['assemblecyc_input_prefix'].replace("{{target_protein_name}}", replacements['target_protein_name'])
    assemblecyc_input_suffix = config['assemblecyc_input_suffix']
    
    assemblecyc_output_prefix_resolved = config['assemblecyc_output_prefix'].replace("{{target_protein_name}}", replacements['target_protein_name'])

    dockmodel_dir_for_assemblecyc = _resolve_path(base_dir, config['dockmodel_dir'], replacements)

    command_assemblecyc = [
        sys.executable, assemblecyc_script,
        "--fraglinking_dir", assemblecyc_fraglinking_dir,
        "--assembled_cyc_dir", assemblecyc_output_dir,
        "--input_prefix", assemblecyc_input_prefix,
        "--input_suffix", assemblecyc_input_suffix,
        "--output_prefix", assemblecyc_output_prefix_resolved,
        "--dockmodel_dir", dockmodel_dir_for_assemblecyc
    ]
    if not _run_command(command_assemblecyc, cwd=tools_dir, step_name=f"{step_name} (assemblecyc.py)")[0]: return False, None

    # SelectBuildcomplex (C program)
    select_build_exec = _resolve_path(base_dir, config['select_build_exec'], replacements)
    select_build_input_dir = _resolve_path(base_dir, config['select_build_input_dir'], replacements)
    select_build_complex_model_dir = _resolve_path(base_dir, config['select_build_complex_model_dir'], replacements)
    select_build_fixed_pdb = _resolve_path(base_dir, config['select_build_fixed_pdb'], replacements)
    
    # --- 新增行：确保 SelectBuildcomplex 的输出目录存在 ---
    if not _create_dir(select_build_complex_model_dir, step_name=f"{step_name} (SelectBuildcomplex output dir)"): return False, None
    # --- 结束新增行 ---

    select_build_output_prefix = config['select_build_output_prefix'].replace("{{target_protein_name}}", replacements['target_protein_name'])

    input_prefix = config['select_build_input_prefix'].replace("{{target_protein_name}}", replacements['target_protein_name'])
    input_suffix = config['select_build_input_suffix']
    file_pattern = re.compile(rf"^{re.escape(input_prefix)}(\d+){re.escape(input_suffix)}$")
    
    found_files_info = []
    try:
        for filename in os.listdir(select_build_input_dir):
            match = file_pattern.match(filename)
            if match and os.path.isfile(os.path.join(select_build_input_dir, filename)):
                file_number = int(match.group(1))
                found_files_info.append((file_number, filename))
        found_files_info.sort()
    except Exception as e:
        print(f"🔴 {step_name}: Error reading input directory for SelectBuildcomplex: {e}")
        return False, None

    if not found_files_info:
        print(f"🔴 {step_name}: No files found for SelectBuildcomplex matching '{input_prefix}X{input_suffix}' in '{select_build_input_dir}'.")
        return False, None

    print(f"--- {step_name}: Running SelectBuildcomplex for {len(found_files_info)} files ---")
    for i, (file_number, input_filename_base) in enumerate(found_files_info, 1):
        full_input_path = os.path.join(select_build_input_dir, input_filename_base)
        
        command_select_build = [
            select_build_exec,
            full_input_path,
            str(config['select_build_param1_score']),
            str(config['select_build_param2_int']),
            str(config['select_build_param3_score']),
            str(config['select_build_param4_int']),
            select_build_complex_model_dir,
            select_build_fixed_pdb,
            select_build_output_prefix
        ]
        if not _run_command(command_select_build, step_name=f"{step_name} (SelectBuildcomplex {i}/{len(found_files_info)})")[0]: return False, None
        
    return True, None
def main():
    parser = argparse.ArgumentParser(description="Automated pipeline for SDOCK fragment docking and assembly.")
    parser.add_argument("config_file", type=str, help="Path to the YAML configuration file.")
    args = parser.parse_args()

    try:
        with open(args.config_file, 'r') as f:
            config = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"Error: Configuration file '{args.config_file}' not found.")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"Error: Invalid YAML format in '{args.config_file}': {e}")
        sys.exit(1)

    project_base_dir = os.path.expanduser(config.get('project_base_dir', os.getcwd()))
    if not os.path.isdir(project_base_dir):
        print(f"Error: project_base_dir '{project_base_dir}' specified in config does not exist.")
        sys.exit(1)
    
    target_protein_name = config.get('target_protein_name')
    if not target_protein_name:
        print("Error: 'target_protein_name' must be specified in the config.yml file.")
        sys.exit(1)

    sdock_base_dir = config.get('sdock_base_dir')
    if not sdock_base_dir:
        print("Error: 'sdock_base_dir' must be specified in the config.yml file.")
        sys.exit(1)

    preprocessed_fragments_path_raw = config.get('preprocessed_fragments_path')
    water_fragments_path_raw = config.get('water_fragments_path')

    if not preprocessed_fragments_path_raw:
        print("Error: 'preprocessed_fragments_path' must be specified in the config.yml file.")
        sys.exit(1)
    if not water_fragments_path_raw:
        print("Error: 'water_fragments_path' must be specified in the config.yml file.")
        sys.exit(1)

    print(f"--- Starting SDOCK Pipeline ---")
    print(f"Using project base directory: {project_base_dir}")
    print(f"Target Protein Name: {target_protein_name}")
    print(f"SDOCK Base Directory: {sdock_base_dir}")
    print(f"Pre-processed Fragments Path (from config): {preprocessed_fragments_path_raw}")
    print(f"Water Fragments Path (from config): {water_fragments_path_raw}")
    print(f"Loaded configuration from: {args.config_file}")
    print("-" * 40)

    pipeline_dynamic_params = {}
    
    global_replacements = {
        'target_protein_name': target_protein_name,
        'sdock_base_dir': sdock_base_dir,
        'preprocessed_fragments_path': _resolve_path(project_base_dir, preprocessed_fragments_path_raw, {'target_protein_name': target_protein_name}),
        'water_fragments_path': _resolve_path(project_base_dir, water_fragments_path_raw, {'target_protein_name': target_protein_name})
    }

    steps = [
        ("step1_download_pdb", run_step1_download_pdb),
        ("step2_pdb2fragments", run_step2_pdb2fragments),
        ("step3_rotate_protein", run_step3_rotate_protein),
        ("step4_generate_boxes_watmaps", run_step4_generate_boxes_watmaps),
        ("step5_process_fragments", run_step5_process_fragments),
        ("step6_generate_frag_files", run_step6_generate_frag_files),
        ("step7_perform_docking", run_step7_perform_docking),
        ("step8_select_top_n_per_box", run_step8_select_top_n_per_box),
        ("step9_generate_complex_models", run_step9_generate_complex_models),
        ("step10_find_linker", run_step10_find_linker),
        ("step11_assemble_cyc", run_step11_assemble_cyc),
    ]

    for step_key, step_func in steps:
        step_config = config.get(step_key)
        if step_config is None:
            print(f"Skipping {step_key}: Configuration not found in YAML.")
            continue
        if not step_config.get('enabled', True):
            print(f"Skipping {step_key}: Disabled in configuration.")
            continue

        print(f"\n===== Executing {step_key} =====")
        
        current_dynamic_params = None # 默认为 None

        # 根据步骤的依赖关系，决定是否传递 dynamic_params
        if step_key == "step4_generate_boxes_watmaps":
            current_dynamic_params = pipeline_dynamic_params.get('step3_output')
            success, result_data = step_func(step_config, project_base_dir, global_replacements, dynamic_params=current_dynamic_params)
        elif step_key == "step5_process_fragments":
            # Step 5 不需要 dynamic_params 作为输入，但它的输出需要被捕获
            success, result_data = step_func(step_config, project_base_dir, global_replacements) # 保持没有 dynamic_params 参数
        elif step_key == "step6_generate_frag_files":
            current_dynamic_params = pipeline_dynamic_params.get('step5_process_fragments_output') # 从这里获取 step5 的输出
            success, result_data = step_func(step_config, project_base_dir, global_replacements, dynamic_params=current_dynamic_params)
        elif step_key == "step9_generate_complex_models":
            current_dynamic_params = pipeline_dynamic_params.get('step3_output')
            success, result_data = step_func(step_config, project_base_dir, global_replacements, dynamic_params=current_dynamic_params)
        else:
            # 对于不需要特定 dynamic_params 的其他步骤
            success, result_data = step_func(step_config, project_base_dir, global_replacements)
        
        if not success:
            print(f"\n🔴 Pipeline FAILED at {step_key}. Exiting.")
            sys.exit(1)
        
        if result_data: # 只要 result_data 不为空，就存储它
            pipeline_dynamic_params[f'{step_key}_output'] = result_data

        print(f"===== {step_key} Completed Successfully =====")

    print("\n--- SDOCK Pipeline Completed Successfully! ---")

if __name__ == "__main__":
    main()
