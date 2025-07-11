import heapq
import os
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

def process_docking_file(file_path, filter_rank_enabled):
    # Use a dictionary to store the best score for each unique combination of (peptide, file_path, model_no)
    # Key: (peptide, model_no) -> Value: (score, file_path)
    # We include file_path in the value so we can reconstruct the full tuple later.
    best_results_for_file = {}
    current_peptide = file_path # Initialize with file_path as a fallback

    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                # Update current_peptide from "Protein B: mobile" line
                if line.startswith("   Protein B: mobile"):
                    current_peptide = line.split()[-1]
                    continue # Move to the next line

                line_stripped = line.strip()
                if not line_stripped or line_stripped.startswith('%') or line_stripped.startswith('[source'):
                    continue # Skip empty lines, comments, and source lines

                # Attempt to parse lines containing model_no, model_rank, and score
                try:
                    # Check for lines where model_no, model_rank, and score are present
                    # Example: "1  1 -35.817 : some other info"
                    parts = line_stripped.split()
                    if len(parts) >= 3: # Ensure there are enough parts for model, rank, and score
                        model_no = int(parts[0])
                        model_rank = int(parts[1])
                        
                        # Find the score part. It's usually the third part if no ':'
                        # Or the part before ':' if ':' exists.
                        score_str_candidate = parts[2]
                        if ':' in line_stripped:
                            # If colon exists, take the part before it and then get the last number
                            before_colon = line_stripped.split(':')[0]
                            score_str = before_colon.strip().split()[-1]
                        else:
                            score_str = score_str_candidate # Assume it's the third part

                        score = float(score_str)
                        
                        # Apply rank filter
                        if not (filter_rank_enabled and model_rank > 1):
                            # Define a unique key for the dictionary
                            unique_key = (current_peptide, file_path, model_no)
                            
                            # If this combination is new, or if the current score is better than recorded
                            if unique_key not in best_results_for_file or score < best_results_for_file[unique_key][0]:
                                best_results_for_file[unique_key] = (score, model_no)
                                
                except (ValueError, IndexError):
                    # This block handles lines that don't fit the expected parsing pattern for scores
                    # and model numbers, such as pure "model_no model_rank" lines without a score,
                    # or other malformed lines.
                    continue # Skip to the next line if parsing fails

    except FileNotFoundError:
        print(f"Warning: File {file_path} not found, skipping.")
    except Exception as e:
        print(f"Error: Unknown error processing file {file_path}: {e}")
        
    # Convert the dictionary back to a list of tuples in the desired format
    # (score, peptide_path, source_file_path, model_no)
    final_results = []
    for (peptide_path, source_file_path, model_no), (score, _) in best_results_for_file.items():
        final_results.append((score, peptide_path, source_file_path, model_no))
        
    return final_results


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
                # The alignment in the C code was %s %s %2d %8.3f
                # For Python, we'll try to match the column widths.
                # peptide_path: Left-aligned, 80 characters
                # source_file_path: Left-aligned, 80 characters
                # model_no: Right-aligned, 3 characters
                # score: Right-aligned, 8.3f (8 characters total, 3 decimal places)
                f_out.write(f"{peptide_path:<80} {source_file_path:<80} {model_no:>3} {score:8.3f}\n")
    except IOError as e:
        print(f"Error: Could not write to output file '{output_path}': {e}")

    print(f"Done! Top {len(best_results)} results saved to '{output_path}'")

if __name__ == "__main__":
    main()