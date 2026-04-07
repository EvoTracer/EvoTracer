import pandas as pd
import os
import sys

# ================= 配置区域 =================

# 默认尝试读取的文件路径 (优先读取合并后的，如果没有则读取未合并的)
DEFAULT_FILES = [
    "Final_Large_SV_Analysis_Merged/Merged_Large_SV_Details.csv",
    "Final_Large_SV_Analysis/Large_SV_Details.csv",
    "Large_SV_Details.csv" # 同时也查找当前目录
]

# ================= 核心功能函数 =================

def load_data():
    """尝试加载 CSV 数据文件"""
    for f in DEFAULT_FILES:
        if os.path.exists(f):
            print(f"成功加载数据文件: {f}")
            return pd.read_csv(f)
    
    print("错误: 未找到结果 CSV 文件。")
    print(f"请确保以下任一文件存在: {DEFAULT_FILES}")
    return None

def query_region(df, query_start, query_end):
    """
    查询指定区间内的事件
    判定标准: 事件区间与查询区间存在重叠 (Overlap)
    """
    # 筛选逻辑: max(Event_Start, Query_Start) <= min(Event_End, Query_End)
    # 简化为: (Event_Start <= Query_End) AND (Event_End >= Query_Start)
    mask = (df['Start_EC'] <= query_end) & (df['End_EC'] >= query_start)
    result = df[mask].copy()
    
    return result

def print_stats(result_df, start, end):
    """打印统计结果"""
    print("\n" + "="*50)
    print(f"  查询区间: aEC {start} bp - {end} bp")
    print("="*50)
    
    count = len(result_df)
    print(f"\n[1] 总事件数: {count}")
    
    if count == 0:
        print("\n在该区间内未检测到任何 >1KB 的结构变异事件。")
        return

    # 统计节点分布
    print("\n[2] 节点(分支)事件统计:")
    # 按照 Parent -> Child 分组计数
    branch_stats = result_df.groupby(['Parent', 'Child']).size().reset_index(name='Events_Count')
    branch_stats = branch_stats.sort_values(by='Events_Count', ascending=False)
    
    print("-" * 40)
    print(f"{'Parent':<10} -> {'Child':<10} | {'Count':<5}")
    print("-" * 40)
    for _, row in branch_stats.iterrows():
        print(f"{row['Parent']:<10} -> {row['Child']:<10} | {row['Events_Count']:<5}")
    
    # 详细列表
    print("\n[3] 详细事件列表:")
    print("-" * 100)
    # 选取关键列进行展示
    cols = ['Parent', 'Child', 'Type', 'Start_EC', 'End_EC', 'Actual_Seq_Length', 'Affected_Genes']
    
    # 格式化输出
    # 如果没有 Affected_Genes 列，做个兼容处理
    if 'Affected_Genes' not in result_df.columns:
        result_df['Affected_Genes'] = 'N/A'

    header = f"{'Branch':<15} | {'Type':<12} | {'Range (EC)':<18} | {'Length':<8} | {'Genes'}"
    print(header)
    print("-" * 100)
    
    for _, row in result_df.iterrows():
        branch = f"{row['Parent']}->{row['Child']}"
        rng = f"{row['Start_EC']}-{row['End_EC']}"
        # 截断过长的基因列表
        genes = str(row['Affected_Genes'])
        if len(genes) > 40: genes = genes[:37] + "..."
        
        print(f"{branch:<15} | {row['Type']:<12} | {rng:<18} | {row['Actual_Seq_Length']:<8} | {genes}")
    print("-" * 100)

# ================= 主程序 =================

def main():
    print(">>> 结构变异(SV) 区域查询工具 <<<")
    
    # 1. 加载数据
    df = load_data()
    if df is None:
        return

    while True:
        print("\n请输入查询坐标 (输入 q 退出):")
        try:
            inp_start = input("  起始位置 (Start): ").strip()
            if inp_start.lower() == 'q': break
            if not inp_start: continue
            
            inp_end = input("  结束位置 (End):   ").strip()
            if inp_end.lower() == 'q': break
            if not inp_end: continue

            start = int(inp_start)
            end = int(inp_end)
            
            if start > end:
                print("错误: 起始位置不能大于结束位置。")
                continue
                
            # 2. 执行查询
            results = query_region(df, start, end)
            
            # 3. 展示结果
            print_stats(results, start, end)
            
        except ValueError:
            print("错误: 请输入有效的整数坐标。")
        except KeyboardInterrupt:
            print("\n退出。")
            break

if __name__ == "__main__":
    main()