import os
import subprocess
import argparse
import yaml
import sys
import re
import glob

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
            check=True,  # Raise CalledProcessError for non-zero exit codes
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

def run_step1_download_pdb(config, base_dir):
    step_name = "Step 1: Download PDB"
    if not config.get('enabled', True): return True, None
    
    script_path = os.path.join(base_dir, config['script_path'])
    input_file = os.path.join(base_dir, config['input_file'])
    output_dir = os.path.join(base_dir, config['output_dir'])
    
    if not _create_dir(output_dir, step_name): return False, None

    command = [
        sys.executable, script_path,
        "--input_file", input_file,
        "--output_dir", output_dir
    ]
    cwd = os.path.dirname(script_path) 
    success, _ = _run_command(command, cwd=cwd, step_name=step_name)
    return success, None

def run_step2_pdb2fragments(config, base_dir):
    step_name = "Step 2: PDB to Fragments"
    if not config.get('enabled', True): return True, None

    script_path = os.path.join(base_dir, config['script_path'])
    
    cwd = os.path.dirname(script_path)
    command = [sys.executable, script_path]
    
    success, _ = _run_command(command, cwd=cwd, step_name=step_name)
    return success, None

def run_step3_rotate_protein(config, base_dir):
    step_name = "Step 3: Rotate Protein"
    if not config.get('enabled', True): return True, None

    script_path = os.path.join(base_dir, config['script_path'])
    input_pdb = os.path.join(base_dir, config['input_pdb'])
    output_pdb = os.path.join(base_dir, config['output_pdb'])

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
            if line.strip().startswith("Box1 (参考原子1):"):
                parts = line.strip().split()
                if len(parts) >= 5:
                    box1_coords = f"{parts[3]},{parts[4]},{parts[5]}"
            elif line.strip().startswith("Box2 (参考原子2):"):
                parts = line.strip().split()
                if len(parts) >= 5:
                    box2_coords = f"{parts[3]},{parts[4]},{parts[5]}"
        
        if box1_coords and box2_coords:
            print(f"Extracted Box1 coords: {box1_coords}")
            print(f"Extracted Box2 coords: {box2_coords}")
            return True, {'box1_m_param': box1_coords, 'box2_m_param': box2_coords}
        else:
            print(f"🔴 {step_name} Failed: Could not extract Box1 or Box2 coordinates from output.")
            return False, None
    
    return success, None

def run_step4_generate_boxes_watmaps(config, base_dir, dynamic_params=None):
    step_name = "Step 4: Generate Boxes and Watmaps"
    if not config.get('enabled', True): return True, None

    sdock_restrict_dir = os.path.join(base_dir, config['sdock_restrict_dir'])
    preprocess_exec = os.path.join(sdock_restrict_dir, "preprocess")
    watmap_exec = os.path.join(sdock_restrict_dir, "watmap")
    rotated_pdb = os.path.join(base_dir, config['rotated_pdb'])
    
    box1_m_param = dynamic_params.get('box1_m_param') if dynamic_params else config['box1_m_param']
    box2_m_param = dynamic_params.get('box2_m_param') if dynamic_params else config['box2_m_param']

    if not box1_m_param or not box2_m_param:
        print(f"🔴 {step_name} Failed: Box M parameters are missing. Ensure Step 3 ran successfully or they are defined in config.")
        return False, None

    # Preprocess Box 1
    command_box1_preprocess = [
        preprocess_exec, rotated_pdb,
        "-o", os.path.join(base_dir, config['box1_output_pdb']),
        "-a", os.path.join(sdock_restrict_dir, "ATM"),
        "-m", box1_m_param
    ]
    if not _run_command(command_box1_preprocess, step_name=f"{step_name} (Box1 Preprocess)")[0]: return False, None

    # Preprocess Box 2
    command_box2_preprocess = [
        preprocess_exec, rotated_pdb,
        "-o", os.path.join(base_dir, config['box2_output_pdb']),
        "-a", os.path.join(sdock_restrict_dir, "ATM"),
        "-m", box2_m_param
    ]
    if not _run_command(command_box2_preprocess, step_name=f"{step_name} (Box2 Preprocess)")[0]: return False, None

    # Watmap Box 1
    command_box1_watmap = [
        watmap_exec,
        os.path.join(base_dir, config['box1_output_pdb']),
        os.path.join(base_dir, config['box1_watmap_output'])
    ]
    if not _run_command(command_box1_watmap, step_name=f"{step_name} (Box1 Watmap)")[0]: return False, None

    # Watmap Box 2
    command_box2_watmap = [
        watmap_exec,
        os.path.join(base_dir, config['box2_output_pdb']),
        os.path.join(base_dir, config['box2_watmap_output'])
    ]
    if not _run_command(command_box2_watmap, step_name=f"{step_name} (Box2 Watmap)")[0]: return False, None
    
    return True, None

def run_step5_process_fragments(config, base_dir):
    step_name = "Step 5: Process Fragments"
    if not config.get('enabled', True): return True, None

    fraglib_dir = os.path.join(base_dir, config['fraglib_dir'])
    
    if not os.path.isdir(fraglib_dir):
        print(f"🔴 {step_name}: Fragment library directory not found: {fraglib_dir}")
        return False, None

    # fragpreprocess.sh
    fragpreprocess_script = os.path.join(fraglib_dir, "fragpreprocess.sh")
    if not os.path.exists(fragpreprocess_script):
        print(f"🔴 {step_name}: fragpreprocess.sh not found at {fragpreprocess_script}")
        return False, None
    if not _run_command(["bash", fragpreprocess_script], cwd=fraglib_dir, step_name=f"{step_name} (Frag Preprocess)")[0]: return False, None

    # fragwater.sh
    fragwater_script = os.path.join(fraglib_dir, "fragwater.sh")
    if not os.path.exists(fragwater_script):
        print(f"🔴 {step_name}: fragwater.sh not found at {fragwater_script}")
        return False, None
    if not _run_command(["bash", fragwater_script], cwd=fraglib_dir, step_name=f"{step_name} (Frag Water)")[0]: return False, None
    
    return True, None

def run_step6_generate_frag_files(config, base_dir):
    step_name = "Step 6: Generate Frag Files"
    if not config.get('enabled', True): return True, None

    tools_dir = os.path.join(base_dir, config['tools_dir'])

    # fragstruct.py
    fragstruct_script = os.path.join(tools_dir, "fragstruct.py")
    fraglib_input_dir = os.path.join(base_dir, config['fraglib_input_dir'])
    fragstruct_output_txt = os.path.join(base_dir, config['fragstruct_output_txt'])
    command_fragstruct = [
        sys.executable, fragstruct_script,
        "--input_dir", fraglib_input_dir,
        "--output_txt", fragstruct_output_txt,
        "--tail", config['fragstruct_tail']
    ]
    if not _run_command(command_fragstruct, cwd=tools_dir, step_name=f"{step_name} (fragstruct.py)")[0]: return False, None

    # struct2wmap.py
    struct2wmap_script = os.path.join(tools_dir, "struct2wmap.py")
    fragwmap_output_file = os.path.join(base_dir, config['fragwmap_output_file'])
    command_struct2wmap = [
        sys.executable, struct2wmap_script,
        "--input_file", fragstruct_output_txt,
        "--output_file", fragwmap_output_file,
        "--old_string", config['struct2wmap_old_string'],
        "--new_string", config['struct2wmap_new_string']
    ]
    if not _run_command(command_struct2wmap, cwd=tools_dir, step_name=f"{step_name} (struct2wmap.py)")[0]: return False, None

    # mkdir dockresult
    dockresult_dir = os.path.join(base_dir, config['dockresult_dir'])
    if not _create_dir(dockresult_dir, step_name): return False, None

    # box1recordf.py
    recordf_script = os.path.join(tools_dir, "recordf.py")
    box1recordf_output_file = os.path.join(base_dir, config['box1recordf_output_file'])
    # Ensure prefix_string is correctly interpreted by recordf.py (absolute vs relative)
    box1recordf_prefix_string = os.path.join(base_dir, config['box1recordf_prefix_string'])
    command_box1recordf = [
        sys.executable, recordf_script,
        "--input_file", fragstruct_output_txt,
        "--output_file", box1recordf_output_file,
        "--prefix_string", box1recordf_prefix_string
    ]
    if not _run_command(command_box1recordf, cwd=tools_dir, step_name=f"{step_name} (box1recordf.py)")[0]: return False, None

    # box2recordf.py
    box2recordf_output_file = os.path.join(base_dir, config['box2recordf_output_file'])
    # Ensure prefix_string is correctly interpreted by recordf.py (absolute vs relative)
    box2recordf_prefix_string = os.path.join(base_dir, config['box2recordf_prefix_string'])
    command_box2recordf = [
        sys.executable, recordf_script,
        "--input_file", fragstruct_output_txt,
        "--output_file", box2recordf_output_file,
        "--prefix_string", box2recordf_prefix_string
    ]
    if not _run_command(command_box2recordf, cwd=tools_dir, step_name=f"{step_name} (box2recordf.py)")[0]: return False, None
    
    return True, None

def run_step7_perform_docking(config, base_dir):
    step_name = "Step 7: Perform Docking"
    if not config.get('enabled', True): return True, None

    tools_dir = os.path.join(base_dir, config['tools_dir'])
    docking_base_dir = os.path.join(base_dir, config['docking_base_dir'])
    split_n_parts = config['split_n_parts']
    dock_num_tasks = config['dock_num_tasks']
    sdock_executable = os.path.join(base_dir, config['sdock_executable'])
    sdock_common_params = config['sdock_common_params'].split()

    # Create directories for split files
    fragstruct_split_dir = os.path.join(docking_base_dir, 'fragstruct_')
    fragwmap_split_dir = os.path.join(docking_base_dir, 'fragwmap_')
    recordf_split_dir = os.path.join(docking_base_dir, 'recordf_')

    if not _create_dir(fragstruct_split_dir, step_name): return False, None
    if not _create_dir(fragwmap_split_dir, step_name): return False, None
    if not _create_dir(recordf_split_dir, step_name): return False, None

    split_script = os.path.join(tools_dir, "splitfragstruct.py")
    
    # Split fragstruct
    fragstruct_input = os.path.join(docking_base_dir, 'fragstruct')
    command_split_fragstruct = [
        sys.executable, split_script,
        "--input_file", fragstruct_input,
        "--output_prefix", os.path.join(fragstruct_split_dir, 'fragstruct'),
        "--n_parts", str(split_n_parts)
    ]
    if not _run_command(command_split_fragstruct, cwd=tools_dir, step_name=f"{step_name} (Split fragstruct)")[0]: return False, None

    # Split fragwmap
    fragwmap_input = os.path.join(docking_base_dir, 'fragwmap')
    command_split_fragwmap = [
        sys.executable, split_script,
        "--input_file", fragwmap_input,
        "--output_prefix", os.path.join(fragwmap_split_dir, 'fragwmap'),
        "--n_parts", str(split_n_parts)
    ]
    if not _run_command(command_split_fragwmap, cwd=tools_dir, step_name=f"{step_name} (Split fragwmap)")[0]: return False, None

    # Split box1recordf
    box1recordf_input = os.path.join(docking_base_dir, 'box1recordf')
    command_split_box1recordf = [
        sys.executable, split_script,
        "--input_file", box1recordf_input,
        "--output_prefix", os.path.join(recordf_split_dir, 'box1recordf'),
        "--n_parts", str(split_n_parts)
    ]
    if not _run_command(command_split_box1recordf, cwd=tools_dir, step_name=f"{step_name} (Split box1recordf)")[0]: return False, None

    # Split box2recordf
    box2recordf_input = os.path.join(docking_base_dir, 'box2recordf')
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
        "--fixed_pdb", os.path.join(base_dir, config['box1_fixed_pdb']),
        "--fragstruct_dir", fragstruct_split_dir,
        "--fixed_wat_pdb", os.path.join(base_dir, config['box1_fixed_wat_pdb']),
        "--fragwatmap_dir", fragwmap_split_dir,
        "--output_record_dir", recordf_split_dir, # This should be the 'record' subdir inside docking_base_path
        "--output_prefix", config['box1_output_prefix'],
        "--sdock_exec", sdock_executable,
        "--common_params", ' '.join(sdock_common_params) # Pass as single string
    ]
    if not _run_command(command_dock_box1, cwd=tools_dir, step_name=f"{step_name} (Dock Box1)")[0]: return False, None

    # Dock box2
    command_dock_box2 = [
        sys.executable, dock_script,
        "--num_tasks", str(dock_num_tasks),
        "--fixed_pdb", os.path.join(base_dir, config['box2_fixed_pdb']),
        "--fragstruct_dir", fragstruct_split_dir,
        "--fixed_wat_pdb", os.path.join(base_dir, config['box2_fixed_wat_pdb']),
        "--fragwatmap_dir", fragwmap_split_dir,
        "--output_record_dir", recordf_split_dir, # This should be the 'record' subdir inside docking_base_path
        "--output_prefix", config['box2_output_prefix'],
        "--sdock_exec", sdock_executable,
        "--common_params", ' '.join(sdock_common_params) # Pass as single string
    ]
    if not _run_command(command_dock_box2, cwd=tools_dir, step_name=f"{step_name} (Dock Box2)")[0]: return False, None

    return True, None

def run_step8_select_top_n_per_box(config, base_dir):
    step_name = "Step 8: Select Top N per Box"
    if not config.get('enabled', True): return True, None

    script_path = os.path.join(base_dir, config['script_path'])
    top_n = config['top_n']

    # Box1
    box1_recordf = os.path.join(base_dir, config['box1_recordf'])
    box1_output_path = os.path.join(base_dir, config['box1_output_path'])
    command_box1 = [
        sys.executable, script_path,
        "--file_list_path", box1_recordf, # getbestdocking.py expects file_list_path
        "--output_file_path", box1_output_path,
        "--top_n", str(top_n)
    ]
    if not _run_command(command_box1, cwd=os.path.dirname(script_path), step_name=f"{step_name} (Box1)")[0]: return False, None

    # Box2
    box2_recordf = os.path.join(base_dir, config['box2_recordf'])
    box2_output_path = os.path.join(base_dir, config['box2_output_path'])
    command_box2 = [
        sys.executable, script_path,
        "--file_list_path", box2_recordf, # getbestdocking.py expects file_list_path
        "--output_file_path", box2_output_path,
        "--top_n", str(top_n)
    ]
    if not _run_command(command_box2, cwd=os.path.dirname(script_path), step_name=f"{step_name} (Box2)")[0]: return False, None
    
    return True, None

def run_step9_generate_complex_models(config, base_dir):
    step_name = "Step 9: Generate Complex Models"
    if not config.get('enabled', True): return True, None

    utility_dir = os.path.join(base_dir, config['utility_dir'])
    sdock_restrict_dir = os.path.join(base_dir, config['sdock_restrict_dir'])
    genbuild_exec = os.path.join(utility_dir, "genBuildcommand")
    build_exec = os.path.join(sdock_restrict_dir, "build")
    dockmodel_dir = os.path.join(base_dir, config['dockmodel_dir'])
    so3layer_qua = os.path.join(base_dir, config['so3layer_qua'])

    if not _create_dir(dockmodel_dir, step_name): return False, None

    # genBuildcommand Box1
    box1_bestdocking = os.path.join(base_dir, config['box1_bestdocking'])
    box1_gen_script = os.path.join(base_dir, config['box1_gen_script'])
    box1_fixed_pdb = os.path.join(base_dir, config['box1_fixed_pdb'])
    command_genbuild_box1 = [
        genbuild_exec, box1_bestdocking, box1_gen_script,
        build_exec, box1_fixed_pdb, dockmodel_dir,
        so3layer_qua, config['box1_m_param']
    ]
    if not _run_command(command_genbuild_box1, step_name=f"{step_name} (genBuildcommand Box1)")[0]: return False, None

    # genBuildcommand Box2
    box2_bestdocking = os.path.join(base_dir, config['box2_bestdocking'])
    box2_gen_script = os.path.join(base_dir, config['box2_gen_script'])
    box2_fixed_pdb = os.path.join(base_dir, config['box2_fixed_pdb'])
    command_genbuild_box2 = [
        genbuild_exec, box2_bestdocking, box2_gen_script,
        build_exec, box2_fixed_pdb, dockmodel_dir,
        so3layer_qua, config['box2_m_param']
    ]
    if not _run_command(command_genbuild_box2, step_name=f"{step_name} (genBuildcommand Box2)")[0]: return False, None

    # Execute generated shell scripts
    if not _run_command(["bash", box1_gen_script], step_name=f"{step_name} (Execute genbox1frag.sh)")[0]: return False, None
    if not _run_command(["bash", box2_gen_script], step_name=f"{step_name} (Execute genbox2frag.sh)")[0]: return False, None

    # ls commands
    box1dockfrag_output = os.path.join(base_dir, 'ALK1/box1dockfrag') # Hardcoded as per user's example
    box2dockfrag_output = os.path.join(base_dir, 'ALK1/box2dockfrag') # Hardcoded as per user's example
    
    try:
        with open(box1dockfrag_output, 'w') as f:
            subprocess.run(["ls", os.path.join(dockmodel_dir, "SDOCK_ALK1box1*")], stdout=f, check=True, text=True)
        print(f"✅ {step_name}: ls output for box1 saved to {box1dockfrag_output}")

        with open(box2dockfrag_output, 'w') as f:
            subprocess.run(["ls", os.path.join(dockmodel_dir, "SDOCK_ALK1box2*")], stdout=f, check=True, text=True)
        print(f"✅ {step_name}: ls output for box2 saved to {box2dockfrag_output}")
    except subprocess.CalledProcessError as e:
        print(f"🔴 {step_name}: Failed to list files for dockfrag: {e}")
        return False, None
    except Exception as e:
        print(f"🔴 {step_name}: Error during ls command: {e}")
        return False, None

    # getfragpair.py (assuming it's a Python script accepting CLI args)
    getfragpair_script = os.path.join(base_dir, config['getfragpair_exec']) # Renamed from _exec to _script for clarity as it's Python
    fragpair_output_file = os.path.join(base_dir, config['fragpair_output_file'])

    command_getfragpair = [
        sys.executable, getfragpair_script,
        "--getfragpair_exec", getfragpair_script, # This seems redundant if getfragpair.py is the script itself.
                                                  # If getfragpair_exec in config points to the C program, then this needs re-evaluation.
                                                  # Assuming getfragpair_exec points to the Python wrapper script.
        "--parent_dir", os.path.dirname(box1dockfrag_output), # Assuming parent_dir is where box1dockfrag and box2dockfrag are
        "--prefix1", "box1dockfrag", # These prefixes don't match the original getfragpair.py
        "--prefix2", "box2dockfrag", # These prefixes don't match the original getfragpair.py
        "--num_files_per_group", "1", # Assuming these are single files, not groups
        "--score_arg", config['fragpair_score_arg'],
        "--length_arg", config['fragpair_length_arg'],
        "--output_file_name", os.path.splitext(os.path.basename(fragpair_output_file))[0] # Just the base name
    ]
    # Re-evaluating getfragpair.py call based on its last CLI arguments
    # The getfragpair.py script I generated expects:
    # --getfragpair_exec (path to C program)
    # --parent_dir (parent dir of box1_ and box2_ files)
    # --prefix1, --prefix2 (prefixes for numbered files)
    # --num_files_per_group (number of files)
    # --score_arg, --length_arg
    # --output_file_name

    # The user's example: ../utility/getfragpair ../ALK1/box1dockfrag ../ALK1/box2dockfrag -56 8 ../ALK1/ALK1dockfragpair
    # This looks like a direct call to the C program, not the Python wrapper.
    # Let's assume getfragpair_exec in config points to the C executable.

    getfragpair_c_exec = os.path.join(base_dir, config['getfragpair_exec']) # This should be the C executable
    command_getfragpair_c = [
        getfragpair_c_exec,
        box1dockfrag_output, # This is the file list, not a directory
        box2dockfrag_output, # This is the file list, not a directory
        config['fragpair_score_arg'],
        config['fragpair_length_arg'],
        fragpair_output_file
    ]
    if not _run_command(command_getfragpair_c, step_name=f"{step_name} (getfragpair C program)")[0]: return False, None

    return True, None

def run_step10_find_linker(config, base_dir):
    step_name = "Step 10: Find Linker"
    if not config.get('enabled', True): return True, None

    tools_dir = os.path.join(base_dir, config['tools_dir'])
    fraglinking_dir = os.path.join(base_dir, config['fraglinking_dir'])
    fraglink_script = os.path.join(tools_dir, config['fraglink_script'])

    if not _create_dir(fraglinking_dir, step_name): return False, None

    command_fraglink = [sys.executable, fraglink_script]
    if not _run_command(command_fraglink, cwd=tools_dir, step_name=f"{step_name} (fraglink.py)")[0]: return False, None
    
    return True, None

def run_step11_assemble_cyc(config, base_dir):
    step_name = "Step 11: Assemble Cyc"
    if not config.get('enabled', True): return True, None

    tools_dir = os.path.join(base_dir, config['tools_dir'])
    utility_dir = os.path.join(base_dir, config['utility_dir']) # Used for SelectBuildcomplex
    assembled_cyc_dir = os.path.join(base_dir, config['assembled_cyc_dir'])

    if not _create_dir(assembled_cyc_dir, step_name): return False, None

    # assemblecyc.py
    assemblecyc_script = os.path.join(tools_dir, config['assemblecyc_script'])
    assemblecyc_fraglinking_dir = os.path.join(base_dir, config['assemblecyc_fraglinking_dir'])
    assemblecyc_output_dir = os.path.join(base_dir, config['assemblecyc_output_dir'])
    command_assemblecyc = [
        sys.executable, assemblecyc_script,
        "--fraglinking_dir", assemblecyc_fraglinking_dir,
        "--assembled_cyc_dir", assemblecyc_output_dir
    ]
    if not _run_command(command_assemblecyc, cwd=tools_dir, step_name=f"{step_name} (assemblecyc.py)")[0]: return False, None

    # SelectBuildcomplex (C program)
    select_build_exec = os.path.join(base_dir, config['select_build_exec'])
    select_build_input_dir = os.path.join(base_dir, config['select_build_input_dir'])
    select_build_complex_model_dir = os.path.join(base_dir, config['select_build_complex_model_dir'])
    select_build_fixed_pdb = os.path.join(base_dir, config['select_build_fixed_pdb'])

    # Get all files matching the pattern in select_build_input_dir
    input_prefix = config['select_build_input_prefix']
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
            select_build_fixed_pdb
        ]
        if not _run_command(command_select_build, step_name=f"{step_name} (SelectBuildcomplex {i}/{len(found_files_info)})")[0]: return False, None
        
    return True, None

# --- Main Pipeline Execution ---
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
    
    print(f"--- Starting SDOCK Pipeline ---")
    print(f"Using project base directory: {project_base_dir}")
    print(f"Loaded configuration from: {args.config_file}")
    print("-" * 40)

    # Dictionary to store dynamic parameters passed between steps
    pipeline_dynamic_params = {}

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
        
        # Pass dynamic parameters if available for relevant steps
        if step_key == "step4_generate_boxes_watmaps":
            success, result_data = step_func(step_config, project_base_dir, dynamic_params=pipeline_dynamic_params.get('step3_output'))
        else:
            success, result_data = step_func(step_config, project_base_dir)
        
        if not success:
            print(f"\n🔴 Pipeline FAILED at {step_key}. Exiting.")
            sys.exit(1)
        
        # Store results for subsequent steps if needed
        if result_data:
            pipeline_dynamic_params[f'{step_key}_output'] = result_data

        print(f"===== {step_key} Completed Successfully =====")

    print("\n--- SDOCK Pipeline Completed Successfully! ---")

if __name__ == "__main__":
    main()