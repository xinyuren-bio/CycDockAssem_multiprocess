import re
import argparse

def extract_and_prepend_to_file(input_filename, output_filename, prefix_string):
    pattern = re.compile(r'pre_[tTP]p_(.{10})\.pdb')

    try:
        with open(input_filename, 'r') as infile:
            lines = infile.readlines()

        with open(output_filename, 'w') as outfile:
            for line in lines:
                match = pattern.search(line)
                if match:
                    extracted_string = match.group(1)
                    final_line = prefix_string + extracted_string + '\n'
                    outfile.write(final_line)
                else:
                    print(f"警告：行 '{line.strip()}' 未匹配到预期模式，已跳过。")
        print(f"提取并添加前缀完成！结果已保存到 '{output_filename}'。")

    except FileNotFoundError:
        print(f"错误：文件 '{input_filename}' 未找到。")
    except Exception as e:
        print(f"发生错误: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extracts a specific pattern from each line, prepends a string, and writes to an output file.')
    parser.add_argument('--input_file', type=str, required=True, help='Path to the input file.')
    parser.add_argument('--output_file', type=str, required=True, help='Path to the output file.')
    parser.add_argument('--prefix_string', type=str, required=True, help='The string to prepend to the extracted content.')

    args = parser.parse_args()

    extract_and_prepend_to_file(args.input_file, args.output_file, args.prefix_string)
