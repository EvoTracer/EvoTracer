import os

def rename_faa_to_fasta(directory):
    # 检查目录是否存在
    if not os.path.exists(directory):
        print(f"Error: Directory '{directory}' does not exist.")
        return

    # 获取目录下的所有文件
    files = os.listdir(directory)
    count = 0
    
    print(f"Scanning directory: {directory}")

    for filename in files:
        if filename.endswith(".faa"):
            old_path = os.path.join(directory, filename)
            # 生成新文件名：去掉扩展名 .faa 加上 .fasta
            # 或者直接替换结尾，注意防止文件名本身包含.faa
            new_filename = filename[:-4] + ".fasta"
            new_path = os.path.join(directory, new_filename)
            
            try:
                os.rename(old_path, new_path)
                print(f"Renamed: {filename} -> {new_filename}")
                count += 1
            except OSError as e:
                print(f"Error renaming {filename}: {e}")

    print(f"Done. Renamed {count} files.")

if __name__ == "__main__":
    target_dir = "/home/wangyiwei/BactPG/BactPG2.0_20240204/test"
    rename_faa_to_fasta(target_dir)