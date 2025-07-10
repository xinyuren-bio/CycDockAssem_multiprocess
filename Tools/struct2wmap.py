import argparse

def replace_text_in_file(input_filename, output_filename, old_string, new_string):
    try:
        with open(input_filename, 'r') as infile:
            lines = infile.readlines()

        with open(output_filename, 'w') as outfile:
            for line in lines:
                modified_line = line.replace(old_string, new_string)
                outfile.write(modified_line)
        print(f"替换完成！修改后的内容已写入 '{output_filename}'。")

    except FileNotFoundError:
        print(f"错误：文件 '{input_filename}' 未找到。")
    except Exception as e:
        print(f"发生错误: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Reads a file, replaces a string in each line, and writes to a new output file.')
    parser.add_argument('--input_file', type=str, required=True, help='Path to the input file.')
    parser.add_argument('--output_file', type=str, required=True, help='Path to the output file.')
    parser.add_argument('--old_string', type=str, required=True, help='The string to be replaced.')
    parser.add_argument('--new_string', type=str, required=True, help='The new string to replace with.')

    args = parser.parse_args()

    replace_text_in_file(args.input_file, args.output_file, args.old_string, args.new_string)
