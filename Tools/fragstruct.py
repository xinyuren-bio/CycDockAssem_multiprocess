import os
import argparse

def list_files_to_txt(directory_path, output_file_path, file_extension=None):
    try:
        with open(output_file_path, 'w') as outfile:
            for filename in os.listdir(directory_path):
                if file_extension and not filename.endswith(file_extension):
                    continue
                full_path = os.path.join(directory_path, filename)
                if os.path.isfile(full_path):
                    outfile.write(full_path + '\n')
        print(f"File paths successfully saved to: {output_file_path}")
    except FileNotFoundError:
        print(f"Error: Directory not found at {directory_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Lists files in a directory and saves their paths to a text file.')
    parser.add_argument('--input_dir', type=str, required=True, help='The path to the directory to scan.')
    parser.add_argument('--output_txt', type=str, required=True, help='The path to the output text file.')
    parser.add_argument('--tail', type=str, default=None, help='Only include files with this extension (e.g., ".pdb"). If None, all files are included.')

    args = parser.parse_args()

    list_files_to_txt(args.input_dir, args.output_txt, args.tail)
