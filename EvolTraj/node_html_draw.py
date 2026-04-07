import pandas as pd
import plotly.graph_objects as go
import os
import numpy as np
from Bio import SeqIO

# ==========================================
# --- 1. 配置区域 ---
# ==========================================

CSV_FILE = '/home/wangyiwei/nodeacientmauve/paintwithtrna/Large_SV_Details.csv'
OUTPUT_HTML = 'Figure_5_Separated.html'
OUTPUT_SVG = 'Figure_5_Separated_Editable.svg'
GBK_DIR = '/home/wangyiwei/nodeacientmauve/allnodeGBK'

# --- 视觉参数 ---
NODE_FILL = '#66c2a5'
NODE_BORDER = '#FFFFFF'
NODE_TEXT_COLOR = '#1A1A1A'
NODE_SIZE = 55
EDGE_COLOR = '#B0B0B0'
EDGE_WIDTH = 2.5
CANVAS_SIZE = 900

# 坐标定义
NODE_COORDS = {
    "aEC":      (5,   55), "aEC1":     (15,  55),
    "aB2G":     (30,  80), "aG":       (50,  92), "aB2":      (55,  75),
    "aAB1EFD":  (25,  35), "aD":       (35,  15), "aAB1EF":   (40,  35),
    "aAB1E":    (55,  35), "aE":       (60,  55), "aAB1":     (70,  30),
    "aB1":      (85,  40), "aA":       (75,  55)
}

# ==========================================
# --- 2. 核心函数 (不变) ---
# ==========================================

def get_bezier_curve(p0, p1, num_points=100):
    x0, y0 = p0; x1, y1 = p1
    dist = abs(x1 - x0) / 1.5
    if dist < 10: dist = 10
    if x1 < x0:
        ctrl_1 = (x0 + 10, y0); ctrl_2 = (x1 - 10, y1)
    else:
        ctrl_1 = (x0 + dist, y0); ctrl_2 = (x1 - dist, y1)
    t = np.linspace(0, 1, num_points)
    x = (1-t)**3 * x0 + 3*(1-t)**2 * t * ctrl_1[0] + 3*(1-t)*t**2 * ctrl_2[0] + t**3 * x1
    y = (1-t)**3 * y0 + 3*(1-t)**2 * t * ctrl_1[1] + 3*(1-t)*t**2 * ctrl_2[1] + t**3 * y1
    return x, y

def format_bp(bp):
    if bp >= 1_000_000: return f"{bp/1_000_000:.2f}Mb"
    elif bp >= 1_000:   return f"{bp/1_000:.1f}kb"
    else:               return f"{bp}bp"

gbk_cache = {}
def load_gbk(node_name):
    if node_name in gbk_cache: return gbk_cache[node_name]
    for d in [GBK_DIR, '.']:
        for ext in ['.gbk', '.gb', '.gbf']:
            path = os.path.join(d, f"{node_name}{ext}")
            if os.path.exists(path):
                try:
                    r = next(SeqIO.parse(path, "genbank"))
                    gbk_cache[node_name] = r
                    return r
                except: pass
    return None

def get_real_gene_count(row):
    p, c, type_, seq = row['Parent'], row['Child'], row['Type'], row['Sequence']
    if type_ == 'Deletion':
        rec = load_gbk(p)
        if not rec: return 0
        s, e = row['Start_EC'], row['End_EC']
        count = 0
        for f in rec.features:
            if f.type == 'CDS':
                if s <= (f.location.start + f.location.end)/2 <= e: count += 1
        return count
    elif type_ == 'Insertion':
        rec = load_gbk(c)
        if not rec or not isinstance(seq, str) or len(seq) < 10: return 0
        full_seq = str(rec.seq).upper()
        idx = full_seq.find(seq.upper())
        if idx != -1:
            s, e = idx, idx + len(seq)
            count = 0
            for f in rec.features:
                if f.type == 'CDS':
                    if s <= (f.location.start + f.location.end)/2 <= e: count += 1
            return count
    return 0

# ==========================================
# --- 3. 这里的逻辑完全变了: 每一个元素都独立添加 ---
# ==========================================

def main():
    if not os.path.exists(CSV_FILE):
        print(f"Error: 找不到 {CSV_FILE}")
        return

    # --- 数据处理 ---
    df = pd.read_csv(CSV_FILE)
    df['Real_Gene_Count'] = df.apply(get_real_gene_count, axis=1)

    edge_stats = {}
    for _, row in df.iterrows():
        key = (row['Parent'], row['Child'])
        if key not in edge_stats:
            edge_stats[key] = {'gain_g':0, 'gain_l':0, 'loss_g':0, 'loss_l':0}
        if row['Type'] == 'Insertion':
            edge_stats[key]['gain_g'] += row['Real_Gene_Count']
            edge_stats[key]['gain_l'] += row['Actual_Seq_Length']
        elif row['Type'] == 'Deletion':
            edge_stats[key]['loss_g'] += row['Real_Gene_Count']
            edge_stats[key]['loss_l'] += row['Actual_Seq_Length']

    fig = go.Figure()

    # ==========================
    # 1. 独立绘制每一条连线 (Lines)
    # ==========================
    for (p, c), stat in edge_stats.items():
        if p in NODE_COORDS and c in NODE_COORDS:
            x0, y0 = NODE_COORDS[p]
            x1, y1 = NODE_COORDS[c]
            xc, yc = get_bezier_curve((x0, y0), (x1, y1))
            
            # 【重点】每次循环都调用 add_trace，而不是把它们合并
            # 这样在 SVG 里，这一条线就是一个独立的 <path> 标签
            fig.add_trace(go.Scatter(
                x=xc, y=yc,
                mode='lines',
                line=dict(color=EDGE_COLOR, width=EDGE_WIDTH),
                name=f"Line_{p}_to_{c}",  # 给个名字，AI图层面板里能看到
                hoverinfo='none',
                showlegend=False,
                opacity=1.0
            ))
            
            # 标注文字
            lines = []
            if stat['gain_l'] > 0:
                lines.append(f"<span style='color:#E64B35'>+{stat['gain_g']} ({format_bp(stat['gain_l'])})</span>")
            if stat['loss_l'] > 0:
                lines.append(f"<span style='color:#3498DB'>-{stat['loss_g']} ({format_bp(stat['loss_l'])})</span>")
            
            if lines:
                label_html = "<br>".join(lines)
                mid_idx = int(len(xc) * 0.5)
                # 简单避让参数
                dy = yc[mid_idx+1] - yc[mid_idx-1]
                ax, ay = 10, 25
                if c == 'aA': ax, ay = -30, -30
                elif c == 'aB1': ax, ay = 30, 20
                elif dy < 0: ax, ay = 10, -25

                # 文字本身就是独立对象，不需要特殊处理
                fig.add_annotation(
                    x=xc[mid_idx], y=yc[mid_idx],
                    text=label_html,
                    showarrow=False,
                    ax=ax, ay=ay,
                    bgcolor="rgba(255,255,255,0.7)",
                    font=dict(size=10, family="Arial")
                )

    # ==========================
    # 2. 独立绘制每一个节点 (Nodes)
    # ==========================
    # 将光晕和圆圈彻底分开绘制，方便你删除光晕或移动圆圈
    
    for node, (nx, ny) in NODE_COORDS.items():
        # 2.1 独立的光晕 (如果不想要光晕，在 AI 里可以直接选中删掉)
        fig.add_trace(go.Scatter(
            x=[nx], y=[ny],
            mode='markers',
            name=f"Halo_{node}", # 独立的图层名
            marker=dict(size=NODE_SIZE + 10, color='rgba(102, 194, 165, 0.4)'),
            showlegend=False,
            hoverinfo='none'
        ))

        # 2.2 独立的实心圆
        fig.add_trace(go.Scatter(
            x=[nx], y=[ny],
            mode='markers', # 纯圆圈
            name=f"Circle_{node}",
            marker=dict(
                size=NODE_SIZE, 
                color=NODE_FILL,
                line=dict(color=NODE_BORDER, width=2)
            ),
            showlegend=False,
            hoverinfo='none'
        ))
        
        # 2.3 独立的文字
        # 我们用 annotation 代替 Scatter text，这样在 SVG 里文字和圆圈完全分离，不会被编组
        fig.add_annotation(
            x=nx, y=ny,
            text=node,
            showarrow=False,
            font=dict(color=NODE_TEXT_COLOR, size=11, family="Arial", weight="bold"),
            name=f"Text_{node}"
        )

    # ==========================
    # 3. 布局设置 (洁癖模式)
    # ==========================
    fig.update_layout(
        width=CANVAS_SIZE,
        height=CANVAS_SIZE,
        plot_bgcolor='rgba(0,0,0,0)', # 透明背景
        paper_bgcolor='rgba(0,0,0,0)',
        margin=dict(l=0, r=0, t=0, b=0), # 零边距
        xaxis=dict(visible=False, showgrid=False, zeroline=False),
        yaxis=dict(visible=False, showgrid=False, zeroline=False),
        showlegend=False
    )

    print(f">>> 导出 SVG: {OUTPUT_SVG}")
    try:
        fig.write_image(OUTPUT_SVG, width=CANVAS_SIZE, height=CANVAS_SIZE)
        print(">>> 完成！现在SVG中每条线和圆都是独立的路径。")
    except Exception as e:
        print(f"Export Error: {e}")
        fig.write_html(OUTPUT_HTML)

if __name__ == "__main__":
    main()