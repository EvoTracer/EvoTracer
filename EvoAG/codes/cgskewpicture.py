import os
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

# 计算GC skew
def gc_skew(seq, window_size=1000):
    skews = []
    for i in range(0, len(seq) - window_size + 1, window_size):
        window = seq[i:i + window_size]
        g = window.count('G')
        c = window.count('C')
        if g + c != 0:
            skew = (g - c) / (g + c)
        else:
            skew = 0
        skews.append(skew)
    return skews

# 绘制并保存GC skew图，带有基因长度标尺
def plot_gc_skew(skews, file_name, output_folder, genome_length):
    plt.figure(figsize=(8, 8))
    theta = np.linspace(0, 2 * np.pi, len(skews))
    radii = skews

    # 设置颜色：正值为绿色，负值为紫色
    colors = ['green' if val >= 0 else 'purple' for val in radii]

    # 绘制GC skew图
    ax = plt.subplot(111, polar=True)
    bars = ax.bar(theta, np.abs(radii), color=colors, width=0.1, bottom=1)
    ax.set_title(file_name, va='bottom', fontsize=20)
    
    # 添加基因长度标尺
    num_ticks = 5  # 标尺刻度数量
    tick_labels = [f"{i * genome_length / num_ticks / 1e6:.1f} Mbp" for i in range(num_ticks + 1)]
    ax.set_xticks(np.linspace(0, 2 * np.pi, num_ticks + 1))
    ax.set_xticklabels(tick_labels, fontsize=10)
    
    # 保存图像为PNG格式
    output_path = os.path.join(output_folder, f"{file_name}_gc_skew.png")
    plt.savefig(output_path, format="png")
    plt.close()  # 关闭图表以节省内存

# 主函数：遍历指定文件夹中的FASTA文件并生成图
def generate_gc_skew_plots(folder_path, output_folder, window_size=1000):
    # 创建输出文件夹（如果不存在）
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".fasta") or file_name.endswith(".fa"):
            file_path = os.path.join(folder_path, file_name)
            with open(file_path, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    skews = gc_skew(str(record.seq), window_size=window_size)
                    genome_length = len(record.seq)  # 获取基因长度
                    plot_gc_skew(skews, file_name, output_folder, genome_length)

# 使用指定文件夹路径和窗口大小
folder_path = "/home/vandark/mycode/output9.24"
output_folder = "/home/vandark/mycode/cgskew"
generate_gc_skew_plots(folder_path, output_folder)
