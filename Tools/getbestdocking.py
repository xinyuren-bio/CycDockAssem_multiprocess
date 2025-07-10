import heapq
import os
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

def process_docking_file(file_path, filter_rank_enabled):
    results = []
    current_peptide = file_path
    pending_data = None

    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith("   Protein B: mobile"):
                    current_peptide = line.split()[-1]
                    print(current_peptide)
                    continue

                line_stripped = line.strip()
                if not line_stripped or line_stripped.startswith('%') or line_stripped.startswith('[source'):
                    continue
                
                try:
                    if ':' in line_stripped and pending_data:
                        score_str = line_stripped.split(':')[0].strip()
                        score = float(score_str)
                        model_no, model_rank = pending_data
                        
                        if not (filter_rank_enabled and model_rank > 1):
                            results.append((score, current_peptide, file_path, model_no))
                        
                        pending_data = None
                        continue

                    parts = line_stripped.split()
                    model_no = int(parts[0])
                    model_rank = int(parts[1])

                    pending_data = (model_no, model_rank)

                    if ':' in line_stripped:
                        before_colon = line_stripped.split(':')[0]
                        score_str = before_colon.split()[-1]
                        score = float(score_str)
                        
                        if not (filter_rank_enabled and model_rank > 1):
                            results.append((score, current_peptide, file_path, model_no))
                        
                        pending_data = None

                except (ValueError, IndexError):
                    pending_data = None
                    continue

    except FileNotFoundError:
        print(f"Warning: File {file_path} not found, skipping.")
    except Exception as e:
        print(f"Error: Unknown error processing file {file_path}: {e}")
        
    return results

def main():
    parser = argparse.ArgumentParser(description="Processes SDOCK output files to find top N results.")
    parser.add_argument("--recordf", type=str, required=True,
                        help="Path to the file containing a list of docking result file paths.")
    parser.add_argument("--output_path", type=str, required=True,
                        help="Path for the final output file where best results will be saved.")
    parser.add_argument("--top_n", type=int, default=50000,
                        help="Number of best results (Top N) to filter out (default: 50000).")
    parser.add_argument("--max_workers", type=int, default=None,
                        help="Number of threads to use for parallel processing (default: None, which means use as many as CPU cores).")
    parser.add_argument("--filter_rank", action='store_true',
                        help="Enable filtering to only include results with rank = 1.")

    args = parser.parse_args()

    recordf = args.recordf
    output_path = args.output_path
    top_n = args.top_n
    max_workers = args.max_workers
    filter_rank = args.filter_rank

    try:
        with open(recordf, 'r') as f:
            docking_files = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"Error: File list '{recordf}' not found. Please check configuration.")
        return

    print(f"Configuration loaded:")
    print(f"  - Reading file list from '{recordf}'.")
    print(f"  - Found {len(docking_files)} files to process.")
    print(f"  - Filtering for top {top_n} best results.")
    print(f"  - Results will be saved to '{output_path}'.")
    if not filter_rank:
        print("  - Note: Results with rank > 1 will NOT be filtered.")
    else:
        print("  - Only results with rank = 1 will be kept.")

    all_results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_file = {executor.submit(process_docking_file, path, filter_rank): path for path in docking_files}
        
        for i, future in enumerate(as_completed(future_to_file), 1):
            try:
                print(f"\rProcessing progress: {i}/{len(docking_files)}", end="")
                file_results = future.result()
                if file_results:
                    all_results.extend(file_results)
            except Exception as e:
                file_path_error = future_to_file[future]
                print(f"\nError: Task {file_path_error} raised an exception: {e}")
    
    print("\nAll files processed. Starting sorting and filtering...")

    if not all_results:
        print("Warning: No valid results found. Please check file content and code logic.")
        return
        
    best_results = heapq.nsmallest(top_n, all_results)

    try:
        with open(output_path, 'w') as f_out:
            for score, peptide_path, source_file_path, model_no in best_results:
                f_out.write(f"{peptide_path:<80} {source_file_path:<80} {model_no:>3} {score:8.3f}\n")
    except IOError as e:
        print(f"Error: Could not write to output file '{output_path}': {e}")

    print(f"Done! Top {len(best_results)} results saved to '{output_path}'")

if __name__ == "__main__":
    main()
