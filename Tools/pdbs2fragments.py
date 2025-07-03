import os
import argparse
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
from tqdm import tqdm
import multiprocessing

warnings.simplefilter('ignore', PDBConstructionWarning)

def extract_fragment_from_pdb(pdb_filepath, structure_id, chain_id_to_extract,
                               start_res_num, end_res_num, output_pdb_filepath,
                               quiet=False):
    parser = PDBParser(QUIET=True)
    try:
        original_structure = parser.get_structure(structure_id, pdb_filepath)
    except FileNotFoundError:
        if not quiet: print(f"Error: PDB file not found at {pdb_filepath} (Worker)")
        return False # Indicate failure
    except Exception as e:
        if not quiet: print(f"Error parsing PDB file {pdb_filepath}: {e} (Worker)")
        return False

    sb = StructureBuilder()
    sb.init_structure(f"frag_{structure_id}_{chain_id_to_extract}_{start_res_num}_{end_res_num}")
    sb.init_model(0)
    sb.init_chain(chain_id_to_extract)

    residues_extracted_count = 0
    try:
        original_chain = original_structure[0][chain_id_to_extract]
        for res_num in range(start_res_num, end_res_num + 1):
            residue_id_tuple_standard = (' ', res_num, ' ')
            residue_to_add = None
            if residue_id_tuple_standard in original_chain:
                residue_to_add = original_chain[residue_id_tuple_standard]
            
            if residue_to_add:
                try:
                    sb.structure[0][chain_id_to_extract].add(residue_to_add.copy())
                    residues_extracted_count += 1
                except Exception as e_add:
                    if not quiet: print(f"Error adding residue {res_num} to fragment from {pdb_filepath}: {e_add} (Worker)")

        if residues_extracted_count == 0 :
            if not quiet and (end_res_num >= start_res_num) :
                 print(f"No residues extracted for {start_res_num}-{end_res_num} from {pdb_filepath}:{chain_id_to_extract} (Worker).")
            return False
        elif residues_extracted_count < (end_res_num - start_res_num + 1):
             if not quiet: print(f"Warning: Incomplete fragment {start_res_num}-{end_res_num} from {pdb_filepath}. Got {residues_extracted_count} (Worker).")
             # Decide if incomplete fragments should be saved or return False

    except KeyError:
        # Attempt to get available chain IDs if chain is not found
        available_chains_str = "N/A"
        try:
            available_chains_str = str([c.id for c in original_structure[0]])
        except Exception:
            pass # Ignore error if structure itself was problematic
        if not quiet: print(f"Chain '{chain_id_to_extract}' not found in {pdb_filepath} (model 0). Avail: {available_chains_str} (Worker)")
        return False
    except Exception as e:
        if not quiet: print(f"An error occurred accessing structure for {pdb_filepath}: {e} (Worker)")
        return False

    fragment_structure = sb.get_structure()
    io = PDBIO()
    io.set_structure(fragment_structure)
    try:
        io.save(output_pdb_filepath)
        return True
    except Exception as e:
        if not quiet: print(f"Error saving fragment PDB file {output_pdb_filepath}: {e} (Worker)")
        return False


def find_pdb_file(pdb_dir, pdb_id_4char_lc):
    potential_filenames = [
        f"{pdb_id_4char_lc}.pdb", f"{pdb_id_4char_lc}.pdb.gz",
        f"{pdb_id_4char_lc}.ent", f"{pdb_id_4char_lc}.ent.gz",
    ]
    for fname in potential_filenames:
        fpath = os.path.join(pdb_dir, fname)
        if os.path.exists(fpath): return fpath
    
    pdb_id_4char_uc = pdb_id_4char_lc.upper()
    potential_filenames_uc = [
        f"{pdb_id_4char_uc}.pdb", f"{pdb_id_4char_uc}.pdb.gz",
        f"{pdb_id_4char_uc}.ent", f"{pdb_id_4char_uc}.ent.gz",
    ]
    for fname in potential_filenames_uc:
        fpath = os.path.join(pdb_dir, fname)
        if os.path.exists(fpath): return fpath
    return None

def main():
    parser = argparse.ArgumentParser(description="Extract peptide fragments from PDB files using multiprocessing.")
    parser.add_argument("--input_file", required=True, help="Path to the index file.")
    parser.add_argument("--pdb_dir", required=True, help="Directory containing PDB files.")
    parser.add_argument("--output_dir", required=True, help="Directory to save extracted fragments.")
    parser.add_argument("--length", "-n", type=int, required=True, help="Length of peptide fragments.")
    parser.add_argument("--processes", "-p", type=int, default=os.cpu_count(), help="Number of worker processes to use (default: all CPU cores).")
    parser.add_argument("--quiet_errors", action="store_true", help="Suppress individual fragment extraction error/warning messages in workers.")

    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    fragment_length = args.length
    if fragment_length <= 0:
        print("Error: Fragment length (-n) must be a positive integer.")
        return

    # --- Prepare list of arguments for all tasks ---
    tasks_for_starmap = []
    
    print("Scanning input file to prepare tasks...")
    try:
        with open(args.input_file, 'r') as f_index_scan:
            for line_num, line in enumerate(f_index_scan):
                line = line.strip()
                if not line or line.startswith("#"): continue
                parts = line.split(',')
                if not parts or not parts[0]: continue
                
                pdb_chain_identifier = parts[0].strip()
                if len(pdb_chain_identifier) < 4: continue

                pdb_id_4char = pdb_chain_identifier[:4]
                chain_id_from_identifier = pdb_chain_identifier[4:]
                
                if not chain_id_from_identifier and len(pdb_chain_identifier) == 4:
                    if not args.quiet_errors: print(f"Warning: Task for '{pdb_chain_identifier}' (line {line_num+1}) has no explicit chain. Skipping.")
                    continue
                if not chain_id_from_identifier and len(pdb_chain_identifier) > 4: # Should not happen if format is PDBIDChain
                     if not args.quiet_errors: print(f"Warning: PDB identifier '{pdb_chain_identifier}' (line {line_num+1}) seems to imply a chain but it's empty. Skipping.")
                     continue


                chain_id_to_use = chain_id_from_identifier.upper()
                structure_id_for_biopython = pdb_id_4char.upper()
                pdb_id_for_file_search = pdb_id_4char.lower()

                pdb_filepath = find_pdb_file(args.pdb_dir, pdb_id_for_file_search)
                if not pdb_filepath:
                    if not args.quiet_errors: print(f"PDB file for ID '{pdb_id_4char}' not found. Skipping fragments for entry: {pdb_chain_identifier} (line {line_num+1}).")
                    continue # Skip all fragments for this PDB if its file is not found
                
                starting_residues_str = parts[1:]
                for res_str in starting_residues_str:
                    res_str = res_str.strip()
                    if not res_str: continue
                    try:
                        start_res_num = int(res_str)
                        end_res_num = start_res_num + fragment_length - 1
                        
                        output_filename = f"{pdb_chain_identifier}_{start_res_num}-{end_res_num}_len{fragment_length}.pdb"
                        output_pdb_filepath = os.path.join(args.output_dir, output_filename)
                        
                        tasks_for_starmap.append((
                            pdb_filepath,
                            structure_id_for_biopython,
                            chain_id_to_use,
                            start_res_num,
                            end_res_num,
                            output_pdb_filepath,
                            args.quiet_errors
                        ))
                    except ValueError:
                        if not args.quiet_errors: print(f"Warning: Invalid start_res_num '{res_str}' in {args.input_file} line {line_num+1}")
    except Exception as e:
        print(f"Error reading input file for preparing tasks: {e}")
        return

    if not tasks_for_starmap:
        print("No valid tasks to process after scanning the input file.")
        return
    
    print(f"Prepared {len(tasks_for_starmap)} total fragments to process using {args.processes} processes.")

    pool = multiprocessing.Pool(processes=args.processes)
    results = []
    try:
        results = list(tqdm(pool.starmap(extract_fragment_from_pdb, tasks_for_starmap), 
                            total=len(tasks_for_starmap), 
                            desc="Extracting Fragments", 
                            unit="fragment"))
    except Exception as e:
        print(f"An error occurred during multiprocessing: {e}")
    finally:
        pool.close() # Important to close the pool
        pool.join()  # Wait for all worker processes to finish

    successful_extractions = sum(1 for r in results if r is True) # Count True results
    failed_extractions = len(tasks_for_starmap) - successful_extractions
    
    print(f"\nFragment extraction process completed.")
    print(f"  Successfully extracted: {successful_extractions} fragments.")
    print(f"  Failed/Skipped extractions: {failed_extractions} fragments.")
    print(f"  Output directory: '{args.output_dir}'.")

if __name__ == "__main__":
    multiprocessing.freeze_support() 
    main()