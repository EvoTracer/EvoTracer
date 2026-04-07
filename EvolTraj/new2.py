import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from Bio import SeqIO
import subprocess

# ================= 1. 全局配置区域 =================

# 根节点文件定义
ROOT_NODE = "aEC"
ROOT_FASTA = f"{ROOT_NODE}.fasta"
ROOT_GFF = f"{ROOT_NODE}.gff"

# 演化树结构 (Parent -> Child)
EVOLUTION_TREE = [
    ("aEC", "aEC1"), ("aEC1", "aB2G"), ("aEC1", "aAB1EFD"),
    ("aB2G", "aG"), ("aB2G", "aB2"), ("aAB1EFD", "aD"),
    ("aAB1EFD", "aAB1EF"), ("aAB1EF", "aAB1E"), ("aAB1E", "aE"),
    ("aAB1E", "aAB1"), ("aAB1", "aB1"), ("aAB1", "aA")
]

# 输出目录
OUTPUT_DIR = "Final_Large_SV_Analysis"

# [新增] 最小事件长度阈值
MIN_EVENT_LENGTH = 1000 

# ================= 2. 基础工具函数 =================

def ensure_dir(d):
    """创建输出目录"""
    if not os.path.exists(d): os.makedirs(d)

def parse_gff_all_features(gff_file):
    """解析 GFF 文件 (Genes & tRNAs)"""
    genes_db = []
    trnas_db = []
    
    if not os.path.exists(gff_file): return [], []
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            
            feature_type = parts[2]
            try:
                start, end = int(parts[3]), int(parts[4])
            except ValueError:
                continue
            
            info = parts[8]
            name = "unknown"
            if "Name=" in info: name = info.split("Name=")[1].split(";")[0]
            elif "ID=" in info: name = info.split("ID=")[1].split(";")[0]
            elif "gene=" in info: name = info.split("gene=")[1].split(";")[0]
            elif "product=" in info: name = info.split("product=")[1].split(";")[0]
            
            record = {'start': start, 'end': end, 'name': name}
            
            if feature_type in ['gene', 'CDS']:
                genes_db.append(record)
            if feature_type == 'tRNA':
                trnas_db.append(record)
                
    return genes_db, trnas_db

def find_overlapping_genes(event_start, event_end, genes_db):
    """查找落在事件区间内的所有基因"""
    hit_genes = set()
    for g in genes_db:
        if max(event_start, g['start']) <= min(event_end, g['end']):
            hit_genes.add(g['name'])
    
    if not hit_genes: return "Intergenic"
    return "; ".join(sorted(list(hit_genes)))

def run_mauve_3way(root, parent, child, out):
    """运行 3-way Mauve 比对"""
    if not all(os.path.exists(f) for f in [root, parent, child]):
        print(f"Error: Missing input FASTA files for {root}, {parent}, or {child}")
        return False
    
    cmd = ["progressiveMauve", "--output=" + out, root, parent, child]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError:
        print(f"Error running Mauve for {out}")
        return False

# ================= 3. 核心逻辑：事件解析 =================

def parse_xmfa_events(xmfa_file):
    """解析 XMFA 文件"""
    events = []
    with open(xmfa_file, 'r') as f: lines = f.readlines()
    block_lines = []
    for line in lines:
        line = line.strip()
        if line == '=':
            if block_lines: _process_block(block_lines, events)
            block_lines = []
        else:
            block_lines.append(line)
    if block_lines: _process_block(block_lines, events)
    return events 

def _process_block(lines, events):
    """处理单个 LCB Block"""
    seq_data, seq_starts = {}, {}
    curr_idx = -1
    
    for line in lines:
        if line.startswith('>'):
            parts = line.split()
            idx = int(parts[1].split(':')[0])
            try:
                s = int(parts[1].split(':')[1].split('-')[0])
                if s == 0: s = 1
            except: s = 1
            seq_starts[idx] = s
            seq_data[idx] = []
            curr_idx = idx
        elif curr_idx != -1:
            seq_data[curr_idx].append(line)

    if 1 not in seq_data or 2 not in seq_data or 3 not in seq_data: return

    s_root = "".join(seq_data[1])
    s_parent = "".join(seq_data[2])
    s_child = "".join(seq_data[3])
    max_len = max(len(s_root), len(s_parent), len(s_child))
    
    s_root = s_root.ljust(max_len, '-')
    s_parent = s_parent.ljust(max_len, '-')
    s_child = s_child.ljust(max_len, '-')

    curr_phys = seq_starts[1]
    last_anchor = curr_phys
    
    in_evt, evt_s, evt_t = False, -1, ""
    evt_seq_buffer = []

    for i in range(max_len):
        base_ec = s_root[i]
        base_a = s_parent[i]
        base_b = s_child[i]

        # 1. 坐标与锚点
        if base_ec != '-':
            map_pos = curr_phys
            last_anchor = curr_phys
            curr_phys += 1
        else:
            map_pos = last_anchor
        
        # 2. 事件判定 (B vs A)
        is_diff = False
        curr_t = ""
        curr_base = ""
        
        # Deletion: A有, B无
        if base_a != '-' and base_b == '-': 
            is_diff, curr_t = True, "Deletion"
            curr_base = base_a 
            
        # Insertion: A无, B有
        elif base_a == '-' and base_b != '-': 
            is_diff, curr_t = True, "Insertion"
            curr_base = base_b 
            
        # 3. 记录
        if is_diff:
            if not in_evt: 
                in_evt, evt_s, evt_t = True, map_pos, curr_t
                evt_seq_buffer = [curr_base]
            elif curr_t != evt_t:
                events.append({
                    'start': evt_s, 'end': map_pos, 'type': evt_t,
                    'sequence': "".join(evt_seq_buffer)
                })
                evt_s, evt_t = map_pos, curr_t
                evt_seq_buffer = [curr_base]
            else:
                evt_seq_buffer.append(curr_base)
        else:
            if in_evt:
                events.append({
                    'start': evt_s, 'end': map_pos, 'type': evt_t,
                    'sequence': "".join(evt_seq_buffer)
                })
                in_evt = False
                evt_seq_buffer = []
                
    if in_evt: 
        events.append({
            'start': evt_s, 'end': last_anchor, 'type': evt_t,
            'sequence': "".join(evt_seq_buffer)
        })

# ================= 4. 主流程 =================

def main():
    ensure_dir(OUTPUT_DIR)
    
    print(f"Loading Root Genome: {ROOT_FASTA}")
    if not os.path.exists(ROOT_FASTA):
        print("Error: Root fasta file not found.")
        return
    root_len = len(SeqIO.read(ROOT_FASTA, "fasta").seq)
    
    print(f"Parsing GFF: {ROOT_GFF}")
    genes_db, trna_list = parse_gff_all_features(ROOT_GFF)
    
    all_events_filtered = [] 
    csv_records = []           
    
    print(f"Analysis Configuration: Focusing on events > {MIN_EVENT_LENGTH} bp")

    # 1. 遍历演化树
    for parent, child in EVOLUTION_TREE:
        print(f"\nProcessing Branch: {parent} -> {child}")
        xmfa_path = os.path.join(OUTPUT_DIR, f"{parent}_{child}.xmfa")
        
        if not os.path.exists(xmfa_path):
            if not run_mauve_3way(ROOT_FASTA, f"{parent}.fasta", f"{child}.fasta", xmfa_path): 
                continue
        
        raw_events = parse_xmfa_events(xmfa_path)
        
        count_pass = 0
        for e in raw_events:
            s, e_end = e['start'], e['end']
            seq = e['sequence']
            actual_len = len(seq)
            
            # === [核心修改] 长度过滤 ===
            if actual_len <= MIN_EVENT_LENGTH:
                continue
            
            count_pass += 1
            all_events_filtered.append(e) # 用于绘图
            
            # CSV 记录
            affected = find_overlapping_genes(s, e_end, genes_db)
            csv_records.append({
                'Parent': parent,
                'Child': child,
                'Type': e['type'],
                'Start_EC': s,
                'End_EC': e_end,
                'Ref_Length_EC': max(0, e_end - s),
                'Actual_Seq_Length': actual_len, # 实际变异长度
                'Sequence': seq,
                'Affected_Genes': affected
            })
        print(f"  Found {len(raw_events)} events -> {count_pass} passed filter (> {MIN_EVENT_LENGTH}bp).")

    # 2. 导出 CSV (仅包含大片段)
    print("\nSaving Large SV Data to CSV...")
    df_events = pd.DataFrame(csv_records)
    csv_path = os.path.join(OUTPUT_DIR, "Large_SV_Details.csv")
    df_events.to_csv(csv_path, index=False)
    print(f"  Saved to: {csv_path}")

    # 3. 计算单点覆盖度 (基于过滤后的事件)
    print("Calculating Exact Coverage for Large SVs...")
    x_axis = np.arange(1, root_len + 2) 
    genome_coverage = np.zeros(root_len + 1, dtype=np.int32)
    
    for evt in all_events_filtered:
        s, e = evt['start'], evt['end']
        s = max(0, min(s, root_len))
        e = max(0, min(e, root_len))
        
        if s < e:
            genome_coverage[s:e] += 1
        elif s == e:
            genome_coverage[s] += 1
            
    y_axis = genome_coverage[1:]
    
    if len(y_axis) < len(x_axis):
        y_axis = np.append(y_axis, np.zeros(len(x_axis)-len(y_axis)))
    else:
        y_axis = y_axis[:len(x_axis)]

    # ================= 5. Matplotlib 静态绘图 =================
    print("Generating High-Res Plot for Large SVs...")
    
    plt.figure(figsize=(20, 6), dpi=300) 
    
    # A. 绘制主数据 (蓝色阶梯)
    plt.fill_between(x_axis, y_axis, step="post", 
                     color='#377eb8', alpha=0.6, linewidth=0)
    plt.step(x_axis, y_axis, where='post', color='#377eb8', linewidth=0.5)

    # B. 绘制 tRNA (墨绿色, 加粗)
    if trna_list:
        max_y = np.max(y_axis) if len(y_axis) > 0 else 5
        trna_starts = [t['start'] for t in trna_list]
        trna_heights = [max_y * 0.15] * len(trna_list)
        plt.bar(trna_starts, trna_heights, width=4000, 
                color='#006400', align='center', alpha=1.0, zorder=4)

    # C. 阈值线
    plt.axhline(y=1, color='red', linestyle='--', linewidth=1.5, zorder=3)

    # D. 装饰
    plt.title(f"Evolutionary Hotspots of Large SVs (>{MIN_EVENT_LENGTH/1000}kb, Anchor-based)", fontsize=16, weight='bold')
    plt.xlabel("Ancestral Coordinate (bp)", fontsize=14)
    plt.ylabel("Exact Event Frequency", fontsize=14)
    
    ax = plt.gca()
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    plt.xlim(0, root_len)
    plt.ylim(bottom=0)
    
    # 图例
    blue_patch = mpatches.Patch(color='#377eb8', alpha=0.6, label=f'Large SV Count (>{MIN_EVENT_LENGTH}bp)')
    green_patch = mpatches.Patch(color='#006400', label='tRNA Sites')
    red_line = Line2D([0], [0], color='red', linestyle='--', linewidth=1.5, label='Threshold = 2')
    plt.legend(handles=[blue_patch, green_patch, red_line], loc='upper right', fontsize=12)
    
    plt.tight_layout()
    
    png_path = os.path.join(OUTPUT_DIR, "Large_SV_Hotspots.png")
    pdf_path = os.path.join(OUTPUT_DIR, "Large_SV_Hotspots.pdf")
    
    plt.savefig(png_path, dpi=300)
    plt.savefig(pdf_path)
    print(f"Done. Plots saved to:\n  {png_path}")

if __name__ == "__main__":
    main()