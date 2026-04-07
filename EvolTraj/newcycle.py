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

ROOT_NODE = "aEC"
ROOT_FASTA = f"{ROOT_NODE}.fasta"
ROOT_GFF = f"{ROOT_NODE}.gff"

EVOLUTION_TREE = [
    ("aEC", "aEC1"), ("aEC1", "aB2G"), ("aEC1", "aAB1EFD"),
    ("aB2G", "aG"), ("aB2G", "aB2"), ("aAB1EFD", "aD"),
    ("aAB1EFD", "aAB1EF"), ("aAB1EF", "aAB1E"), ("aAB1E", "aE"),
    ("aAB1E", "aAB1"), ("aAB1", "aB1"), ("aAB1", "aA")
]

OUTPUT_DIR = "Final_Large_SV_Analysis"
MIN_EVENT_LENGTH = 1000 

# ================= 2. 基础工具函数 =================

def ensure_dir(d):
    if not os.path.exists(d): os.makedirs(d)

def parse_gff_all_features(gff_file):
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
    hit_genes = set()
    for g in genes_db:
        if max(event_start, g['start']) <= min(event_end, g['end']):
            hit_genes.add(g['name'])
    if not hit_genes: return "Intergenic"
    return "; ".join(sorted(list(hit_genes)))

def run_mauve_3way(root, parent, child, out):
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

        if base_ec != '-':
            map_pos = curr_phys
            last_anchor = curr_phys
            curr_phys += 1
        else:
            map_pos = last_anchor
        
        is_diff = False
        curr_t, curr_base = "", ""
        
        if base_a != '-' and base_b == '-': 
            is_diff, curr_t = True, "Deletion"
            curr_base = base_a 
        elif base_a == '-' and base_b != '-': 
            is_diff, curr_t = True, "Insertion"
            curr_base = base_b 
            
        if is_diff:
            if not in_evt: 
                in_evt, evt_s, evt_t = True, map_pos, curr_t
                evt_seq_buffer = [curr_base]
            elif curr_t != evt_t:
                events.append({'start': evt_s, 'end': map_pos, 'type': evt_t, 'sequence': "".join(evt_seq_buffer)})
                evt_s, evt_t = map_pos, curr_t
                evt_seq_buffer = [curr_base]
            else:
                evt_seq_buffer.append(curr_base)
        else:
            if in_evt:
                events.append({'start': evt_s, 'end': map_pos, 'type': evt_t, 'sequence': "".join(evt_seq_buffer)})
                in_evt = False
                evt_seq_buffer = []
                
    if in_evt: 
        events.append({'start': evt_s, 'end': last_anchor, 'type': evt_t, 'sequence': "".join(evt_seq_buffer)})

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
            
            if actual_len <= MIN_EVENT_LENGTH:
                continue
            
            count_pass += 1
            all_events_filtered.append(e)
            
            affected = find_overlapping_genes(s, e_end, genes_db)
            csv_records.append({
                'Parent': parent, 'Child': child, 'Type': e['type'],
                'Start_EC': s, 'End_EC': e_end,
                'Ref_Length_EC': max(0, e_end - s),
                'Actual_Seq_Length': actual_len,
                'Sequence': seq, 'Affected_Genes': affected
            })
        print(f"  Found {len(raw_events)} events -> {count_pass} passed filter (> {MIN_EVENT_LENGTH}bp).")

    print("\nSaving Large SV Data to CSV...")
    df_events = pd.DataFrame(csv_records)
    csv_path = os.path.join(OUTPUT_DIR, "Large_SV_Details.csv")
    df_events.to_csv(csv_path, index=False)
    
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

    # ================= 5. Matplotlib SCI-Quality Circular Plot =================
    print("Generating High-Res Circular Plot for Large SVs...")
    
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
    plt.rcParams['svg.fonttype'] = 'none' 
    
    fig, ax = plt.subplots(figsize=(12, 12), subplot_kw={'projection': 'polar'}, dpi=300)
    
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    
    # 极速数据降维算法
    changes = np.where(np.diff(y_axis) != 0)[0]
    xs_reduced = [0]
    ys_reduced = [y_axis[0]]
    for idx in changes:
        xs_reduced.extend([idx, idx])
        ys_reduced.extend([y_axis[idx], y_axis[idx+1]])
    xs_reduced.extend([len(y_axis)-1])
    ys_reduced.extend([y_axis[-1]])
    xs_reduced = np.array(xs_reduced)
    ys_reduced = np.array(ys_reduced)
    
    theta = 2 * np.pi * (xs_reduced / root_len)
    max_freq = max(ys_reduced) if len(ys_reduced) > 0 and max(ys_reduced) > 0 else 1
    
    # 径向放大系数，拉长刻度间距
    y_scale = max(3.0, 15.0 / max_freq) 
    inner_hole_radius = 8.0 
    R_outer = inner_hole_radius + max_freq * y_scale 
    
    ax.set_ylim(0, R_outer + R_outer * 0.1)

    # 核心填充区
    r_data = R_outer - ys_reduced * y_scale
    ax.fill_between(theta, R_outer, r_data, color='#00A087', alpha=0.85, zorder=3, linewidth=0)
    ax.plot(theta, r_data, color='#007a67', linewidth=0.5, zorder=4)
    
    # 标尺对齐 0 刻度
    axis_angle = 0.0  
    ax.plot([axis_angle, axis_angle], [R_outer - max_freq * y_scale, R_outer], color='#2c3e50', linewidth=2, zorder=10)
    ax.text(axis_angle, R_outer - max_freq * y_scale - 1.5, "Frequency", ha='center', va='top', 
            fontsize=12, fontweight='bold', color='#2c3e50', zorder=10)

    for y_val in range(1, int(max_freq) + 1):
        r_val = R_outer - y_val * y_scale
        circle_theta = np.linspace(0, 2*np.pi, 1000)
        ax.plot(circle_theta, np.full_like(circle_theta, r_val), color='#cccccc', linestyle='-', linewidth=0.6, zorder=1)
        ax.plot(circle_theta, np.full_like(circle_theta, r_val), color='white', linestyle='--', linewidth=0.8, alpha=0.8, zorder=3.5)
        
        tick_width = 0.02
        ax.plot([axis_angle - tick_width, axis_angle + tick_width], [r_val, r_val], color='#2c3e50', linewidth=2, zorder=10)
        ax.text(axis_angle + 0.04, r_val, str(y_val), ha='left', va='center', 
                fontsize=11, fontweight='bold', color='#2c3e50', 
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1.5), zorder=11)

    # 外圈基因组骨架
    circle_theta = np.linspace(0, 2*np.pi, 1000)
    ax.plot(circle_theta, np.full_like(circle_theta, R_outer), color='#2c3e50', linewidth=2.5, zorder=5)
    
    # 阈值线
    ax.plot(circle_theta, np.full_like(circle_theta, R_outer - 1 * y_scale), color='#E64B35', linestyle='-', linewidth=1.5, zorder=6)
    
    ax.set_frame_on(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    
    # === [核心修改区：按 0.5 Mb 步长生成刻度] ===
    unit, divisor = "Mb", 1_000_000
    step_bp = 500_000  # 0.5 Mb = 500,000 bp
    
    tick_positions = np.arange(0, root_len, step_bp)
    tick_angles = 2 * np.pi * (tick_positions / root_len)
    
    tick_length = R_outer * 0.025
    text_offset = R_outer * 0.06

    for angle, pos in zip(tick_angles, tick_positions):
        ax.plot([angle, angle], [R_outer, R_outer + tick_length], color='#2c3e50', linewidth=1.5, zorder=5)
        rot = np.degrees(angle)
        if 90 < rot < 270: 
            rot += 180
            
        # 强制格式化为 1 位小数，比如 0.0 Mb, 0.5 Mb, 1.0 Mb
        label_text = f"{pos/divisor:.1f} {unit}"
        ax.text(angle, R_outer + text_offset, label_text, ha='center', va='center', 
                rotation=rot, rotation_mode='anchor', fontsize=12, color='#333333', fontweight='medium')
    # ==========================================

    # 自定义刻度 (保持高亮的 3.70 和 3.95)
    custom_points = [3.70, 3.95]
    for val in custom_points:
        bp_pos = val * divisor
        if bp_pos <= root_len: 
            custom_angle = 2 * np.pi * (bp_pos / root_len)
            ax.plot([custom_angle, custom_angle], [R_outer, R_outer + tick_length * 1.5], color='#E64B35', linewidth=2.5, zorder=7)
            
            rot = np.degrees(custom_angle)
            if 90 < rot < 270: 
                rot += 180
                
            label_text = f"{val:.2f}"
            ax.text(custom_angle, R_outer + text_offset * 1.3, label_text, ha='center', va='center', 
                    rotation=rot, rotation_mode='anchor', fontsize=12, color='#E64B35', fontweight='bold', zorder=8)

    # 中心标题
    ax.text(0, 0, f"Evolutionary Hotspots\n> {MIN_EVENT_LENGTH/1000} kb SVs", 
            ha='center', va='center', fontsize=16, fontweight='bold', color='#2c3e50')
    
    # 图例配置
    blue_patch = mpatches.Patch(color='#00A087', alpha=0.85, label=f'SV Frequency')
    red_line = Line2D([0], [0], color='#E64B35', linestyle='-', linewidth=1.5, label='Threshold = 1')
    ax.legend(handles=[blue_patch, red_line], loc='lower right', bbox_to_anchor=(1.15, -0.05), 
              fontsize=12, frameon=False)
    
    # 保存输出
    plt.tight_layout()
    png_path = os.path.join(OUTPUT_DIR, "Large_SV_Hotspots_Circos.png")
    pdf_path = os.path.join(OUTPUT_DIR, "Large_SV_Hotspots_Circos.pdf")
    svg_path = os.path.join(OUTPUT_DIR, "Large_SV_Hotspots_Circos.svg")
    
    plt.savefig(png_path, dpi=600, bbox_inches='tight') 
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.savefig(svg_path, format='svg', bbox_inches='tight', transparent=True)
    
    print(f"Done. High-Res Circular Plots saved to:\n  {png_path}\n  {pdf_path}\n  {svg_path}")

if __name__ == "__main__":
    main()