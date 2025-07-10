#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os

# Define constants
MAX_MODELS = 25

def select_mod(filename, dG, hbn, tot, best_n):
    """
    Reads a multi-model PDB file, selects the best models based on filtering criteria,
    and returns their indices.
    """
    models_data = []
    
    try:
        with open(filename, 'r') as f:
            current_score = 0.0
            current_hb = 0
            
            for line in f:
                if line.startswith("TITLE"):
                    # Start of a new model, create a placeholder
                    models_data.append({'score': 0.0, 'hb': 0})
                elif line.startswith("REMARK"):
                    parts = line[8:].strip().split()
                    if not parts:
                        continue
                    
                    # Ensure we have a model to assign data to
                    if not models_data:
                        continue
                        
                    if parts[0] == "score":
                        models_data[-1]['score'] = float(parts[1])
                    elif " ".join(parts[:3]) == "backbone HB No":
                        models_data[-1]['hb'] = int(parts[3])

    except FileNotFoundError:
        print(f"ERROR: Cannot open protein structure file {filename}", file=sys.stderr)
        sys.exit(1)
    
    # --- Selection Logic ---
    
    # First, find all candidates that pass the initial filters
    candidates = []
    for i, data in enumerate(models_data):
        score = data['score']
        hb = data['hb']
        
        if score <= dG and hb >= hbn:
            metric = score - (hb * 1.5)
            if metric <= tot:
                # Add the metric and original index to our list of candidates
                candidates.append({'metric': metric, 'index': i})
    
    # If there are fewer candidates than best_n, return them all
    if len(candidates) <= best_n:
        return [c['index'] for c in candidates]
        
    # Otherwise, sort the candidates by their metric (lower is better)
    # and take the top `best_n`
    candidates.sort(key=lambda x: x['metric'])
    
    best_indices = [c['index'] for c in candidates[:best_n]]
    
    return best_indices

def write_target(output_handle, target_path):
    """
    Reads a target PDB file and writes its ATOM records to an open file handle.
    """
    try:
        with open(target_path, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    # Truncate line at column 57 to mimic C's line[56]='\n'
                    # and write to the output file
                    output_handle.write(line[:56] + '\n')
    except FileNotFoundError:
        print(f"ERROR: Cannot open target structure file {target_path}", file=sys.stderr)
        sys.exit(1)

def get_output_name(input_path, output_dir, model_index):
    """
    Constructs a new filename for the output complex PDB file.
    """
    # Get the basename of the input file (e.g., 'Cyc_4.pdb')
    base = os.path.basename(input_path)
    # Get the name without the extension (e.g., 'Cyc_4')
    root, _ = os.path.splitext(base)
    
    # Create the new filename
    new_filename = f"{root}_{model_index + 1}.pdb"
    
    # Join with the output directory to get the full path
    return os.path.join(output_dir, new_filename)

def build_complex(filename, selected_indices, output_dir, target_path):
    """
    Reads the multi-model PDB again and writes out the selected models combined
    with the target structure.
    """
    if not selected_indices:
        print("No models were selected to build.")
        return

    try:
        with open(filename, 'r') as f:
            model_index = -1
            is_in_selected_model = False
            output_file_handle = None
            
            for line in f:
                if line.startswith("TITLE"):
                    model_index += 1
                    # Check if this new model is one of the ones we selected
                    if model_index in selected_indices:
                        is_in_selected_model = True
                        complex_filename = get_output_name(filename, output_dir, model_index)
                        
                        try:
                            # Open the new file for writing
                            output_file_handle = open(complex_filename, 'w')
                            # Write the target PDB first
                            write_target(output_file_handle, target_path)
                            output_file_handle.write("TER\n")
                            # Write the TITLE line of the current model
                            output_file_handle.write(line)
                        except IOError:
                            print(f"ERROR: Cannot create output file {complex_filename}", file=sys.stderr)
                            is_in_selected_model = False # Stop writing
                    else:
                        is_in_selected_model = False
                
                elif is_in_selected_model and output_file_handle:
                    # If we are in a selected model block, write the line
                    output_file_handle.write(line)
                    # If this is the end of the model, close the file
                    if line.startswith("END"):
                        output_file_handle.close()
                        output_file_handle = None # Reset handle
                        is_in_selected_model = False

    except FileNotFoundError:
        print(f"ERROR: Cannot open protein structure file {filename}", file=sys.stderr)
        sys.exit(1)


def main():
    """
    Main function to drive the script from the command line.
    """
    if len(sys.argv) != 8:
        print("Usage: python <script.py> <input_pdb> <dG_cutoff> <min_hbn> <tot_cutoff> <best_n> <output_dir> <target_pdb>", file=sys.stderr)
        sys.exit(1)
        
    # Parse command-line arguments
    input_pdb = sys.argv[1]
    dG_cutoff = float(sys.argv[2])
    min_hbn = int(sys.argv[3])
    tot_cutoff = float(sys.argv[4])
    best_n = int(sys.argv[5])
    output_dir = sys.argv[6]
    target_pdb = sys.argv[7]
    
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    # 1. Select the best models and get their indices
    selected_indices = select_mod(input_pdb, dG_cutoff, min_hbn, tot_cutoff, best_n)
    
    # 2. Build the complex PDB files for the selected models
    build_complex(input_pdb, selected_indices, output_dir, target_pdb)
    
    print(f"====> {len(selected_indices)} complex model(s) built for {input_pdb}")


if __name__ == "__main__":
    main()