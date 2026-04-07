import os
import sys
import re

def find_missing_files(tree_files, gene_folder):
    """
    Parses tree files to extract filenames, checks for their existence in a given folder,
    and reports any missing files.

    Args:
        tree_files (list): A list of paths to the tree files.
        gene_folder (str): The path to the folder containing the .fasta files.
    """
    all_names = set()

    # 1. 从所有树文件中提取菌株名称
    for tree_file in tree_files:
        try:
            with open(tree_file, 'r') as f:
                content = f.read()
                # 更新正则表达式，精确查找不包含括号、逗号或空白的字符串
                # 这符合Newick格式中叶子节点的定义
                names = re.findall(r"[^,():\s]+", content)
                
                # 过滤掉可能由分支长度等引入的纯数字或过短的字符串
                filtered_names = {name for name in names if not name.isdigit() and len(name) > 1}
                all_names.update(filtered_names)
        except FileNotFoundError:
            print(f"警告: 树文件未找到: {tree_file}")
        except Exception as e:
            print(f"读取文件 {tree_file} 时出错: {e}")

    if not all_names:
        print("未能从树文件中提取任何有效的菌株名称。")
        return

    print(f"从树文件中总共提取了 {len(all_names)} 个唯一的菌株名称。正在开始检查...")

    # 2. 检查每个文件是否存在
    missing_files = []
    for name in sorted(list(all_names)):
        # 构建预期的 .fasta 文件路径
        expected_file_path = os.path.join(gene_folder, f"{name}.fasta")
        
        # 检查文件是否存在
        if not os.path.exists(expected_file_path):
            missing_files.append(name)

    # 3. 报告结果
    if not missing_files:
        print("\n恭喜！所有文件都已在目标文件夹中找到。")
    else:
        print(f"\n--- 发现 {len(missing_files)} 个缺失的文件 ---")
        for missing_file in missing_files:
            print(missing_file)
        print("------------------------------------")

if __name__ == "__main__":
    # 检查是否提供了命令行参数
  

    # FASTA 文件所在的文件夹
    gene_directory = "/home/vandark/mycode/BACT_AGmuti/input/gene"
    
    # 从命令行获取树文件列表
    input_tree_files = ['/home/vandark/mycode/BACT_AGmuti/input/tree_file/enterica83+2.txt']

    find_missing_files(input_tree_files, gene_directory)