import math
import argparse
import os

def split_file(input_file, output_prefix, n_parts):
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    total_lines = len(lines)
    lines_per_part = math.ceil(total_lines / n_parts)
    
    for i in range(n_parts):
        start_idx = i * lines_per_part
        end_idx = min((i + 1) * lines_per_part, total_lines)
        
        output_file = f"{output_prefix}_{i+1}"
        with open(output_file, 'w') as f:
            f.writelines(lines[start_idx:end_idx])
        
        print(f"已创建文件 {output_file}，包含 {end_idx - start_idx} 行")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Splits a file into multiple parts.')
    parser.add_argument('--input_file', type=str, required=True, help='Path to the input file.')
    parser.add_argument('--output_prefix', type=str, required=True, help='Prefix for the output files (e.g., "output_part").')
    parser.add_argument('--n_parts', type=int, required=True, help='Number of parts to split the file into.')

    args = parser.parse_args()

    split_file(args.input_file, args.output_prefix, args.n_parts)
    