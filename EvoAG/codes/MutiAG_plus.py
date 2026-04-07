import re
import os
import sys
import shutil
import random
import string
import subprocess
from multiprocessing import Pool, Manager
threads = int(sys.argv[1])
#树读取函数

def read_tree(treefile_path):
    pattern1 = re.compile(r'\(([^()\s,]+),([^()\s,]+)\),\(([^()\s,]+),([^()\s,]+)\)')  # 二叉树规则
    pattern2 = re.compile(r'\(([^()\s,]+),([^()\s,]+)\),([^()\s,]+)')  # 三叉树规则1
    pattern3 = re.compile(r'([^()\s,]+),\(([^()\s,]+),([^()\s,]+)\)')  # 三叉树规则2

    file = open(treefile_path, 'r', encoding='utf-8')
    lines = file.readlines()
    lines_str = ''.join(lines)
    matches1 = pattern1.findall(lines_str)
    matches2 = [('@p1__L_@' + m[0], m[1], m[2]) for m in pattern2.findall(lines_str)]
    matches3 = [('@p3__L_@' + m[0], m[1], m[2]) for m in pattern3.findall(lines_str)]

    combined_matches = matches1 + matches2 + matches3
    return combined_matches , lines_str
# 分割函数
def cut(combined_matches):
    grouped_matches = [combined_matches[i:i+4] for i in range(0, len(combined_matches), 4)]
    return grouped_matches
# BactAG1功能集成
def run_bactag1(input1, input2, input3, output):
    log_file = os.path.join(log_dir, f"{output}.log")
    # 为每个进程创建独立的工作目录
    process_dir = os.path.join('temp', output)
    os.makedirs(process_dir, exist_ok=True)
    
    # 复制必要的文件到进程目录
    fastafiles = [f"{input1}.fasta", f"{input2}.fasta", f"{input3}.fasta"]
    for fasta in fastafiles:
        src = os.path.join('input/gene', fasta)
        dst = os.path.join(process_dir, fasta)
        if os.path.exists(src):
            shutil.copy2(src, dst)
        else:
            # 如果在input/gene中找不到，尝试在output文件夹中寻找
            src_output = os.path.join('output', fasta)
            if os.path.exists(src_output):
                shutil.copy2(src_output, dst)
    
    # 复制perl脚本到进程目录
    perl_scripts = ['seqExtFromMultiGenome.pl', 'buchong.py', 'tiaozhen.py']
    for script in perl_scripts:
        for possible_dir in ['codes/Tool/', 'temp']:
            src = os.path.join(possible_dir, script)
            if os.path.exists(src):
                dst = os.path.join(process_dir, script)
                shutil.copy2(src, dst)
                break
    
    with open(log_file, 'w') as fh:
        sys.stdout = fh
        sys.stderr = fh
        try:
            os.chdir(process_dir)
            print(f"Successfully changed directory to {process_dir}")
        except Exception as e:
            print(f"Cannot change directory to {process_dir}: {e}")
            sys.exit(1)

        def run_command(command):
            print(f"Running command: {command}")
            
            result = subprocess.run(command, shell=True, text=True, capture_output=True)
            print(result.stdout)
            if result.stderr:
                print(result.stderr)

        print("(1) AG backbone inference and refinement ...")
        print("Alignment ...")
        run_command(f"progressiveMauve --output={input1}.vs.{input2}.xmfa --output-guide-tree={input1}.vs.{input2}.guide_tree --backbone-output={input1}.vs.{input2}.backbone {input1}.fasta {input2}.fasta")
        
        print("Formatting alignment backbone file ...")
        run_command(f"progBackbonePrep {input1}.vs.{input2}.backbone >{input1}.vs.{input2}.backbone.txt")
        
        print("Reordering the homologous blocks of the backbone file ...")
        run_command(f"homBlkReorder {input1}.vs.{input2}.backbone.txt >{input1}.vs.{input2}.homblk.txt")
        
        print("Parsing the orthologous blocks ...")
        run_command(f"orthoParsing {input1}.vs.{input2}.homblk.txt >{input1}.vs.{input2}.orthBlk.txt")
        
        print("Combining the short discontinuous orthologous blocks (< 1000 bp) ...")
        run_command(f"orthoCombine {input1}.vs.{input2}.orthBlk.txt 1000 >{input1}.vs.{input2}.combined_orthBlk.txt")
        
        print("Refining the AG backbone blocks with the genomes of multiple strains ...")
        run_command(f"orthJoin1 {input1}.vs.{input2}.combined_orthBlk.txt {input1}.vs.{input2}.combined_orthBlk.txt {input1} {input1} >{input1}.vs.{input2}.Joint.txt")
        
        print("Retrieving the AG backbone sequence ...")
        run_command(f"homBB {input1}.vs.{input2}.Joint.txt {input1} >Houtenae_BG0_HOM{output}.txt")
        run_command(f"orthBB Houtenae_BG0_HOM{output}.txt >Houtenae_BG0_ORTH{output}.txt")
        run_command(f"perl seqExtFromMultiGenome.pl {input1} {input2} Houtenae_BG0_ORTH{output}.txt >Houtenae_BG0_ORTH{output}.fasta")
        print("The first step finished ...\n\n")

        print("(2) Recursively patching AG with the genomes of representative Houtenae strains and the orthologous relationship between them and the neighbor ancient genome ...")
        print("1-st round Alignment ...")
        run_command(f"progressiveMauve --output={input1}.vs.Houtenae_BG0_ORTH{output}.xmfa --output-guide-tree={input1}.vs.Houtenae_BG0_ORTH{output}.guide_tree --backbone-output={input1}.vs.Houtenae_BG0_ORTH{output}.backbone {input1}.fasta Houtenae_BG0_ORTH{output}.fasta")
        run_command(f"progressiveMauve --output={input1}.vs.{input3}.xmfa --output-guide-tree={input1}.vs.{input3}.guide_tree --backbone-output={input1}.vs.{input3}.backbone {input1}.fasta {input3}.fasta")
        
        print("Formatting alignment backbone file ...")
        run_command(f"progBackbonePrep {input1}.vs.Houtenae_BG0_ORTH{output}.backbone >{input1}.vs.Houtenae_BG0_ORTH{output}.backbone.txt")
        run_command(f"progBackbonePrep {input1}.vs.{input3}.backbone >{input1}.vs.{input3}.backbone.txt")
        
        print("Reordering the homologous blocks of the backbone file ...")
        run_command(f"homBlkReorder {input1}.vs.Houtenae_BG0_ORTH{output}.backbone.txt >{input1}.vs.Houtenae_BG0_ORTH{output}.homblk.txt")
        run_command(f"homBlkReorder {input1}.vs.{input3}.backbone.txt >{input1}.vs.{input3}.homblk.txt")
        
        print("Parsing the orthologous blocks ...")
        run_command(f"orthoParsing {input1}.vs.Houtenae_BG0_ORTH{output}.homblk.txt >{input1}.vs.Houtenae_BG0_ORTH{output}.orthBlk.txt")
        run_command(f"orthoParsing {input1}.vs.{input3}.homblk.txt >{input1}.vs.{input3}.orthBlk.txt")
        
        print("Patching ...")
        run_command(f"patching {input1}.vs.Houtenae_BG0_ORTH{output}.orthBlk.txt {input1}.vs.{input3}.orthBlk.txt Houtenae_BG0_ORTH{output}.txt {input1} {input3} {input3}_added >Houtenae_BG1_ORTH1{output}.txt")
        run_command(f"python buchong.py {input1}.vs.Houtenae_BG0_ORTH{output}.orthBlk.txt {input1}.vs.{input3}.orthBlk.txt Houtenae_BG1_ORTH1{output}.txt Houtenae_BG1_ORTH2{output}.txt")
        run_command(f"python tiaozhen.py Houtenae_BG1_ORTH2{output}.txt Houtenae_BG1_ORTH{output}.txt {input1} {input3} {input2}")
        
        print("Retrieving the AG sequence ...\n\n")
        run_command(f"perl seqExtFromMultiGenome.pl {input1} {input2} Houtenae_BG1_ORTH{output}.txt >Houtenae_BG1_ORTH{output}.fasta")

        print("2-nd round Alignment ...")
        run_command(f"progressiveMauve --output={input2}.vs.Houtenae_BG1_ORTH{output}.xmfa --output-guide-tree={input2}.vs.Houtenae_BG1_ORTH{output}.guide_tree --backbone-output={input2}.vs.Houtenae_BG1_ORTH{output}.backbone {input2}.fasta Houtenae_BG1_ORTH{output}.fasta")
        run_command(f"progressiveMauve --output={input2}.vs.{input3}.xmfa --output-guide-tree={input2}.vs.{input3}.guide_tree --backbone-output={input2}.vs.{input3}.backbone {input2}.fasta {input3}.fasta")
        
        print("Formatting alignment backbone file ...")
        run_command(f"progBackbonePrep {input2}.vs.Houtenae_BG1_ORTH{output}.backbone >{input2}.vs.Houtenae_BG1_ORTH{output}.backbone.txt")
        run_command(f"progBackbonePrep {input2}.vs.{input3}.backbone >{input2}.vs.{input3}.backbone.txt")
        
        print("Reordering the homologous blocks of the backbone file ...")
        run_command(f"homBlkReorder {input2}.vs.Houtenae_BG1_ORTH{output}.backbone.txt >{input2}.vs.Houtenae_BG1_ORTH{output}.homblk.txt")
        run_command(f"homBlkReorder {input2}.vs.{input3}.backbone.txt >{input2}.vs.{input3}.homblk.txt")
        
        print("Parsing the orthologous blocks ...")
        run_command(f"orthoParsing {input2}.vs.Houtenae_BG1_ORTH{output}.homblk.txt >{input2}.vs.Houtenae_BG1_ORTH{output}.orthBlk.txt")
        run_command(f"orthoParsing {input2}.vs.{input3}.homblk.txt >{input2}.vs.{input3}.orthBlk.txt")
        
        print("Patching ...")
        run_command(f"patching {input2}.vs.Houtenae_BG1_ORTH{output}.orthBlk.txt {input2}.vs.{input3}.orthBlk.txt Houtenae_BG1_ORTH{output}.txt {input2} {input3} {input3}_added >Houtenae_BG2_ORTH1{output}.txt")
        run_command(f"python buchong.py {input2}.vs.Houtenae_BG1_ORTH{output}.orthBlk.txt {input2}.vs.{input3}.orthBlk.txt Houtenae_BG2_ORTH1{output}.txt Houtenae_BG2_ORTH2{output}.txt")
        run_command(f"python tiaozhen.py Houtenae_BG2_ORTH2{output}.txt Houtenae_BG2_ORTH{output}.txt {input2} {input3} {input1}")
        
        print("Retrieving the AG sequence ...\n\n")
        run_command(f"perl seqExtFromMultiGenome.pl {input2} {input1} Houtenae_BG2_ORTH{output}.txt >Houtenae_BG2_ORTH{output}.fasta")
        source_file = f"Houtenae_BG2_ORTH{output}.fasta"
        destination_file = f"{output}.fasta"
        shutil.copy(source_file, destination_file)
        new_location = '../../output/'
        os.makedirs(new_location, exist_ok=True)
        shutil.copy(destination_file, os.path.join(new_location, destination_file))
        print(f"File copied to {new_location}")

    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    os.chdir(original_dir)  
# BactAG2功能集成
def run_bactag2(input1, input2, input3, input4, output):
    # 创建日志文件
    log_file = os.path.join(log_dir, f"{output}.log")
    # 为每个进程创建独立的工作目录
    process_dir = os.path.join('temp', output)
    os.makedirs(process_dir, exist_ok=True)
    
    # 复制必要的文件到进程目录
    fastafiles = [f"{input1}.fasta", f"{input2}.fasta", f"{input3}.fasta", f"{input4}.fasta"]
    for fasta in fastafiles:
        src = os.path.join('input/gene', fasta)
        dst = os.path.join(process_dir, fasta)
        if os.path.exists(src):
            shutil.copy2(src, dst)
        else:
            # 如果在input/gene中找不到，尝试在output文件夹中寻找
            src_output = os.path.join('output', fasta)
            if os.path.exists(src_output):
                shutil.copy2(src_output, dst)
    
    # 复制perl脚本到进程目录
    perl_scripts = ['seqExtFromMultiGenome.pl', 'buchong.py', 'tiaozhen.py']
    for script in perl_scripts:
        for possible_dir in ['codes/Tool/', 'temp']:
            src = os.path.join(possible_dir, script)
            if os.path.exists(src):
                dst = os.path.join(process_dir, script)
                shutil.copy2(src, dst)
                break

    with open(log_file, 'w') as fh:
        # 重定向标准输出和标准错误
        sys.stdout = fh
        sys.stderr = fh
    
        # 切换到指定目录
        try:
            os.chdir(process_dir)
            print(f"Successfully changed directory to {process_dir}")
        except Exception as e:
            print(f"Cannot change directory to {process_dir}: {e}")
            sys.exit(1)

        #执行系统命令并打印输出
        def run_command(command):
            print(f"Running command: {command}")
            # 替换工具路径为绝对路径
            
            
            result = subprocess.run(command, shell=True, text=True, capture_output=True)
            print(result.stdout)
            if result.stderr:
                print(result.stderr)

        #系统命令
        print("(1) AG backbone inference and refinement ...")
        print("Alignment ...")
        run_command(f"progressiveMauve --output={input1}.vs.{input2}.xmfa --output-guide-tree={input1}.vs.{input2}.guide_tree --backbone-output={input1}.vs.{input2}.backbone {input1}.fasta {input2}.fasta")
        
        print("Formatting alignment backbone file ...")
        run_command(f"progBackbonePrep {input1}.vs.{input2}.backbone >{input1}.vs.{input2}.backbone.txt")
        
        print("Reordering the homologous blocks of the backbone file ...")
        run_command(f"homBlkReorder {input1}.vs.{input2}.backbone.txt >{input1}.vs.{input2}.homblk.txt")
        
        print("Parsing the orthologous blocks ...")
        run_command(f"orthoParsing {input1}.vs.{input2}.homblk.txt >{input1}.vs.{input2}.orthBlk.txt")
        
        print("Combining the short discontinuous orthologous blocks (< 1000 bp) ...")
        run_command(f"orthoCombine {input1}.vs.{input2}.orthBlk.txt 1000 >{input1}.vs.{input2}.combined_orthBlk.txt")
        
        print("Refining the AG backbone blocks with the genomes of multiple strains ...")
        run_command(f"orthJoin1 {input1}.vs.{input2}.combined_orthBlk.txt {input1}.vs.{input2}.combined_orthBlk.txt {input1} {input1} >{input1}.vs.{input2}.Joint.txt")
        
        print("Retrieving the AG backbone sequence ...")
        run_command(f"homBB {input1}.vs.{input2}.Joint.txt {input1} >Houtenae_BG0_HOM{output}.txt")
        run_command(f"orthBB Houtenae_BG0_HOM{output}.txt >Houtenae_BG0_ORTH{output}.txt")
        run_command(f"perl seqExtFromMultiGenome.pl {input1} {input2} Houtenae_BG0_ORTH{output}.txt >Houtenae_BG0_ORTH{output}.fasta")
        print("The first step finished ...\n\n")

        print("(2) Recursively patching AG with the genomes of representative Houtenae strains and the orthologous relationship between them and the neighbor ancient genome ...")
        print("1-st round Alignment ...")
        run_command(f"progressiveMauve --output={input1}.vs.Houtenae_BG0_ORTH{output}.xmfa --output-guide-tree={input1}.vs.Houtenae_BG0_ORTH{output}.guide_tree --backbone-output={input1}.vs.Houtenae_BG0_ORTH{output}.backbone {input1}.fasta Houtenae_BG0_ORTH{output}.fasta")
        run_command(f"progressiveMauve --output={input1}.vs.{input3}.xmfa --output-guide-tree={input1}.vs.{input3}.guide_tree --backbone-output={input1}.vs.{input3}.backbone {input1}.fasta {input3}.fasta")
        
        print("Formatting alignment backbone file ...")
        run_command(f"progBackbonePrep {input1}.vs.Houtenae_BG0_ORTH{output}.backbone >{input1}.vs.Houtenae_BG0_ORTH{output}.backbone.txt")
        run_command(f"progBackbonePrep {input1}.vs.{input3}.backbone >{input1}.vs.{input3}.backbone.txt")
        
        print("Reordering the homologous blocks of the backbone file ...")
        run_command(f"homBlkReorder {input1}.vs.Houtenae_BG0_ORTH{output}.backbone.txt >{input1}.vs.Houtenae_BG0_ORTH{output}.homblk.txt")
        run_command(f"homBlkReorder {input1}.vs.{input3}.backbone.txt >{input1}.vs.{input3}.homblk.txt")
        
        print("Parsing the orthologous blocks ...")
        run_command(f"orthoParsing {input1}.vs.Houtenae_BG0_ORTH{output}.homblk.txt >{input1}.vs.Houtenae_BG0_ORTH{output}.orthBlk.txt")
        run_command(f"orthoParsing {input1}.vs.{input3}.homblk.txt >{input1}.vs.{input3}.orthBlk.txt")
        
        print("Patching ...")
        run_command(f"patching {input1}.vs.Houtenae_BG0_ORTH{output}.orthBlk.txt {input1}.vs.{input3}.orthBlk.txt Houtenae_BG0_ORTH{output}.txt {input1} {input3} {input3}_added >Houtenae_BG1_ORTH1{output}.txt")
        run_command(f"python buchong.py {input1}.vs.Houtenae_BG0_ORTH{output}.orthBlk.txt {input1}.vs.{input3}.orthBlk.txt Houtenae_BG1_ORTH1{output}.txt Houtenae_BG1_ORTH2{output}.txt")
        run_command(f"python tiaozhen.py Houtenae_BG1_ORTH2{output}.txt Houtenae_BG1_ORTH{output}.txt {input1} {input3} {input2}")
        
        print("Retrieving the AG sequence ...\n\n")
        run_command(f"perl seqExtFromMultiGenome.pl {input1} {input2} Houtenae_BG1_ORTH{output}.txt >Houtenae_BG1_ORTH{output}.fasta")

        print("2-nd round Alignment ...")
        run_command(f"progressiveMauve --output={input2}.vs.Houtenae_BG1_ORTH{output}.xmfa --output-guide-tree={input2}.vs.Houtenae_BG1_ORTH{output}.guide_tree --backbone-output={input2}.vs.Houtenae_BG1_ORTH{output}.backbone {input2}.fasta Houtenae_BG1_ORTH{output}.fasta")
        run_command(f"progressiveMauve --output={input2}.vs.{input3}.xmfa --output-guide-tree={input2}.vs.{input3}.guide_tree --backbone-output={input2}.vs.{input3}.backbone {input2}.fasta {input3}.fasta")
        
        print("Formatting alignment backbone file ...")
        run_command(f"progBackbonePrep {input2}.vs.Houtenae_BG1_ORTH{output}.backbone >{input2}.vs.Houtenae_BG1_ORTH{output}.backbone.txt")
        run_command(f"progBackbonePrep {input2}.vs.{input3}.backbone >{input2}.vs.{input3}.backbone.txt")
        
        print("Reordering the homologous blocks of the backbone file ...")
        run_command(f"homBlkReorder {input2}.vs.Houtenae_BG1_ORTH{output}.backbone.txt >{input2}.vs.Houtenae_BG1_ORTH{output}.homblk.txt")
        run_command(f"homBlkReorder {input2}.vs.{input3}.backbone.txt >{input2}.vs.{input3}.homblk.txt")
        
        print("Parsing the orthologous blocks ...")
        run_command(f"orthoParsing {input2}.vs.Houtenae_BG1_ORTH{output}.homblk.txt >{input2}.vs.Houtenae_BG1_ORTH{output}.orthBlk.txt")
        run_command(f"orthoParsing {input2}.vs.{input3}.homblk.txt >{input2}.vs.{input3}.orthBlk.txt")
        
        print("Patching ...")
        run_command(f"patching {input2}.vs.Houtenae_BG1_ORTH{output}.orthBlk.txt {input2}.vs.{input3}.orthBlk.txt Houtenae_BG1_ORTH{output}.txt {input2} {input3} {input3}_added >Houtenae_BG2_ORTH1{output}.txt")
        run_command(f"python buchong.py {input2}.vs.Houtenae_BG1_ORTH{output}.orthBlk.txt {input2}.vs.{input3}.orthBlk.txt Houtenae_BG2_ORTH1{output}.txt Houtenae_BG2_ORTH2{output}.txt")
        run_command(f"python tiaozhen.py Houtenae_BG2_ORTH2{output}.txt Houtenae_BG2_ORTH{output}.txt {input2} {input3} {input1}")
        
        print("Retrieving the AG sequence ...\n\n")
        run_command(f"perl seqExtFromMultiGenome.pl {input2} {input1} Houtenae_BG2_ORTH{output}.txt >Houtenae_BG2_ORTH{output}.fasta")
        
        print("1-st round Alignment ...")
        run_command(f"progressiveMauve --output={input1}.vs.Houtenae_BG2_ORTH{output}.xmfa --output-guide-tree={input1}.vs.Houtenae_BG2_ORTH{output}.guide_tree --backbone-output={input1}.vs.Houtenae_BG2_ORTH{output}.backbone {input1}.fasta Houtenae_BG2_ORTH{output}.fasta")
        run_command(f"progressiveMauve --output={input1}.vs.{input4}.xmfa --output-guide-tree={input1}.vs.{input4}.guide_tree --backbone-output={input1}.vs.{input4}.backbone {input1}.fasta {input4}.fasta")
        
        print("Formatting alignment backbone file ...")
        run_command(f"progBackbonePrep {input1}.vs.Houtenae_BG2_ORTH{output}.backbone >{input1}.vs.Houtenae_BG2_ORTH{output}.backbone.txt")
        run_command(f"progBackbonePrep {input1}.vs.{input4}.backbone >{input1}.vs.{input4}.backbone.txt")
        
        print("Reordering the homologous blocks of the backbone file ...")
        run_command(f"homBlkReorder {input1}.vs.Houtenae_BG2_ORTH{output}.backbone.txt >{input1}.vs.Houtenae_BG2_ORTH{output}.homblk.txt")
        run_command(f"homBlkReorder {input1}.vs.{input4}.backbone.txt >{input1}.vs.{input4}.homblk.txt")
        
        print("Parsing the orthologous blocks ...")
        run_command(f"orthoParsing {input1}.vs.Houtenae_BG2_ORTH{output}.homblk.txt >{input1}.vs.Houtenae_BG2_ORTH{output}.orthBlk.txt")
        run_command(f"orthoParsing {input1}.vs.{input4}.homblk.txt >{input1}.vs.{input4}.orthBlk.txt")
        
        print("Patching ...")
        run_command(f"patching1 {input1}.vs.Houtenae_BG2_ORTH{output}.orthBlk.txt {input1}.vs.{input4}.orthBlk.txt Houtenae_BG2_ORTH{output}.txt {input1} {input4} {input4}_added >Houtenae_BG3_ORTH1{output}.txt")
        run_command(f"python buchong.py {input1}.vs.Houtenae_BG2_ORTH{output}.orthBlk.txt {input1}.vs.{input4}.orthBlk.txt Houtenae_BG3_ORTH1{output}.txt Houtenae_BG3_ORTH2{output}.txt")
        run_command(f"python tiaozhen.py Houtenae_BG3_ORTH2{output}.txt Houtenae_BG3_ORTH{output}.txt {input1} {input4} {input2}")
        
        print("Retrieving the AG sequence ...\n\n")
        run_command(f"perl seqExtFromMultiGenome.pl {input1} {input2} Houtenae_BG3_ORTH{output}.txt >Houtenae_BG3_ORTH{output}.fasta")
        
        print("2-nd round Alignment ...")
        run_command(f"progressiveMauve --output={input2}.vs.Houtenae_BG3_ORTH{output}.xmfa --output-guide-tree={input2}.vs.Houtenae_BG3_ORTH{output}.guide_tree --backbone-output={input2}.vs.Houtenae_BG3_ORTH{output}.backbone {input2}.fasta Houtenae_BG3_ORTH{output}.fasta")
        run_command(f"progressiveMauve --output={input2}.vs.{input4}.xmfa --output-guide-tree={input2}.vs.{input4}.guide_tree --backbone-output={input2}.vs.{input4}.backbone {input2}.fasta {input4}.fasta")
        
        print("Formatting alignment backbone file ...")
        run_command(f"progBackbonePrep {input2}.vs.Houtenae_BG3_ORTH{output}.backbone >{input2}.vs.Houtenae_BG3_ORTH{output}.backbone.txt")
        run_command(f"progBackbonePrep {input2}.vs.{input4}.backbone >{input2}.vs.{input4}.backbone.txt")
        
        print("Reordering the homologous blocks of the backbone file ...")
        run_command(f"homBlkReorder {input2}.vs.Houtenae_BG3_ORTH{output}.backbone.txt >{input2}.vs.Houtenae_BG3_ORTH{output}.homblk.txt")
        run_command(f"homBlkReorder {input2}.vs.{input4}.backbone.txt >{input2}.vs.{input4}.homblk.txt")
        
        print("Parsing the orthologous blocks ...")
        run_command(f"orthoParsing {input2}.vs.Houtenae_BG3_ORTH{output}.homblk.txt >{input2}.vs.Houtenae_BG3_ORTH{output}.orthBlk.txt")
        run_command(f"orthoParsing {input2}.vs.{input4}.homblk.txt >{input2}.vs.{input4}.orthBlk.txt")
        
        print("Patching ...")
        run_command(f"patching1 {input2}.vs.Houtenae_BG3_ORTH{output}.orthBlk.txt {input2}.vs.{input4}.orthBlk.txt Houtenae_BG3_ORTH{output}.txt {input2} {input4} {input4}_added >Houtenae_BG4_ORTH1{output}.txt")
        run_command(f"python buchong.py {input2}.vs.Houtenae_BG3_ORTH{output}.orthBlk.txt {input2}.vs.{input4}.orthBlk.txt Houtenae_BG4_ORTH1{output}.txt Houtenae_BG4_ORTH2{output}.txt")
        run_command(f"python tiaozhen.py Houtenae_BG4_ORTH2{output}.txt Houtenae_BG4_ORTH{output}.txt {input2} {input4} {input1}")
        
        print("Retrieving the AG sequence ...\n\n")
        run_command(f"perl seqExtFromMultiGenome.pl {input2} {input1} Houtenae_BG4_ORTH{output}.txt >Houtenae_BG4_ORTH{output}.fasta")
        
        source_file = f"Houtenae_BG4_ORTH{output}.fasta"
        destination_file = f"{output}.fasta"
        shutil.copy(source_file, destination_file)
        new_location = '../../output/'
        os.makedirs(new_location, exist_ok=True)
        shutil.copy(destination_file, os.path.join(new_location, destination_file))
    # 恢复标准输出和标准错误
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    os.chdir(original_dir) 
# 单进程函数
def process_group(group, ID_string, output):
    count = len(group)
    if count == 3:  # 三叉树
        if group[0].startswith('@p1__L_@'):
            group = (group[0].replace('@p1__L_@', '', 1), group[1], group[2])
            run_bactag1(group[0], group[1], group[2], ID_string)
            result = group[0] + "+" + group[1] + " OutsidE " + group[2] + " = " + ID_string
            strre = "("+ group[0] +","+ group[1] +")"
        elif group[0].startswith('@p3__L_@'):
            group = (group[0].replace('@p3__L_@', '', 1), group[1], group[2])
            run_bactag1(group[1], group[2], group[0], ID_string)
            result = group[1] + "+" + group[2] + " OutsidE " + group[0] + " = " + ID_string
            strre = "("+ group[1] +","+ group[2] +")"
    elif count == 4:  # 四叉树
        run_bactag2(group[0], group[1], group[2],group[3], ID_string)
        result = group[0] + "+" + group[1] + " OutsidE " + group[2] + "+" + group[3] + " = " + ID_string
        strre = "("+ group[0] +","+ group[1] +")" 

    output.put((result, strre))  
# 工作函数
def worker_wrapper(args):
    group, ID_string, output = args
    process_group(group, ID_string, output)
# 主函数，多进程
def main(num_threads, treefile_path, IDfilepath):
    while True:
        combined_matches , lines_str= read_tree(treefile_path)
        if not combined_matches:
            print("Fin")
            sys.exit()  # 终止
        print(combined_matches)
        grouped_matches = cut(combined_matches)
        
        with Manager() as manager:
            output = manager.Queue()
            generated_ids = set()
            
            # 读取ID文件内容
            with open(IDfilepath, 'r') as IDfileread:
                content = IDfileread.read()
            
            # 展开grouped_matches并计算总组数
            flattened_groups = [group for groups in grouped_matches for group in groups]
            total_groups = len(flattened_groups)
            
            # 生成total_groups个不同的随机字符串
            while len(generated_ids) < total_groups:
                ID_string = ''.join(random.choice(string.ascii_letters + string.digits) for _ in range(7))
                if ID_string not in generated_ids and ID_string not in content:
                    generated_ids.add(ID_string)
            
            # 将生成的ID字符串转换为列表
            generated_ids = list(generated_ids)
            
            # 配对每个组和一个唯一的随机字符串
            group_seed_pairs = [
                (group, ID_string, output)
                for (group, ID_string) in zip(flattened_groups, generated_ids)
            ]

            # 使用进程池进行并行处理
            with Pool(processes=num_threads) as pool:
                pool.map(worker_wrapper, group_seed_pairs)
            
            # 从队列中读取结果并写入文件
            with open(IDfilepath, 'a') as bactIDfile:
                while not output.empty():
                    result, strre = output.get()
                    print(f"Result: {result}, Strre: {strre}") 
                    bactIDfile.write(result + '\n')
                    # 替换treefile中的匹配项为ID_string
                    lines_str = lines_str.replace(strre, result.split(' = ')[-1])
            # 将修改后的内容写回treefile
            with open(treefile_path, 'w', encoding='utf-8') as file:
                file.write(lines_str)
# 主函数运行
if __name__ == "__main__":
    num_threads = threads # 线程大小后面设置为外部输入
    directory = "input/tree_file"
    files = os.listdir(directory)
    if len(files) == 1:
        treefile_path = os.path.join(directory, files[0])
    else:
        print("There are some problems with the treefile. Please check it")

    fastafilefolder_path = "input/gene/"
    IDfilepath = 'codes/bactID.txt'
    Toolpath = 'codes/Tool/'
    folder_path = "temp/"
    original_dir = os.getcwd() 
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    else:
        for filename1 in os.listdir(folder_path):
            file_path = os.path.join(folder_path, filename1)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(f"Wrong when clean {file_path} : {e}")

    # 不再需要复制工具文件到temp目录，因为使用绝对路径
    # 只复制FASTA文件到temp目录（如果需要的话）
    # 注意：由于每个进程会创建自己的子目录并复制文件，这里可能不需要复制
    
    log_dir = 'log/'
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    main(num_threads, treefile_path, IDfilepath)


    