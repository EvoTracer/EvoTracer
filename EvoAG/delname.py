import os

target_dir = "/home/vandark/mycode/BACT_AGmuti/input/gene"  # 修改为你的目标文件夹路径
pattern = "Salmonella_enterica_subsp._"

for filename in os.listdir(target_dir):
    if pattern in filename:
        new_name = filename.replace(pattern, "")
        old_path = os.path.join(target_dir, filename)
        new_path = os.path.join(target_dir, new_name)
        os.rename(old_path, new_path)
        print(f"Renamed: {filename} -> {new_name}")