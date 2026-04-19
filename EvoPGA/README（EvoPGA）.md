# EvoPGA (Evoerial Pan-Genome Annotation)

[English](#english) | [中文](#chinese)

---

<h2 id="english">English</h2>

**EvoPGA** is a versatile software package designed for extracting gene annotations from GenBank files and mapping them with Pan-Genome (PG) UIDs (e.g., Salmonella 26-genome PG) using mutual best hit alignments for PGAG or RAST annotation formats.

Recently, the core parsers and annotators have been rewritten entirely in **pure Go** (`Evopga` cli tool), eliminating previous Python/Go scattered scripts and hard-coded path dependencies.

### Directory Structure
* `bin/`: Contains the compiled executable binary `Evopga`.
* `Evopga/`: The pure Go source code for the main CLI tool.
* `codes/`: Legacy scripts (`gbParse.py`, `PGA.go`, etc.) kept for reference purposes.
* `EvoCG1.0/`: Associated Evoerial core-genome alignment tool (CG).
* `test/` & `input/`: Test datasets and reference files (e.g., `26_PG.txt`, `nA_AG.gbk`).

### Build Instructions
The core tool `Evopga` requires no external dependencies. To rebuild the binary from source:
```bash
cd Evopga
go build -o ../bin/Evopga .
```

### Usage
The single `Evopga` binary inside `bin/` provides two subcommands:

1. **parse**: Extract feature sequence elements (CDS, rRNA, tRNA, ncRNA, misc_feature) from a given GenBank file.
   ```bash
   ./bin/Evopga parse ./test/nA_AG.gbk > output_PGAG.tab.txt
   ```

2. **annotate**: Read the mutual best hit `.mutbest.filt.txt` files and annotate the specific gene table with PG UIDs.
   ```bash
   ./bin/Evopga annotate -pg ./test/26_PG.txt -tab ./output_PGAG.tab.txt -mode PGAG -strain nA_AG -mutbestDir ./EvoCG1.0/result/out_mutbest_filt > nA_AG_annotated.tab.txt
   ```

To run the whole pipeline automatically (parsing > CG > annotating), use the built-in pipeline:
```bash
./bin/Evopga pipeline -gbk ./test/nA_AG.gbk -pg ./test/26_PG.txt -seq ./test/seq -out ./output
```

---

<h2 id="chinese">中文 (Chinese)</h2>

**EvoPGA** 是一个自动化细菌泛基因组注释工具包，旨在解析 GenBank 文件，并根据相互最佳序列比对（Mutual Best Hit）结果为基因注释表（PGAG 或 RAST 格式）自动添加相应的泛基因组簇 UID 注释（如 沙门氏菌 26 基因组 PG 聚类）。

近期，其核心解析与标注代码已被完全重写为一个**无前置依赖纯 Go 语言开发的命令行工具** (`Evopga`)，从而解决了过去老代码对 Python 的依赖及严重硬编码路径导致无法跨设备移植的问题。

### 目录结构
* `bin/`: 存放本软件包编译得到的核心可执行文件 `Evopga`。
* `Evopga/`: 现代 Go 语言环境重写构建的核心源码目录。
* `codes/`: 早期的 Python 与 Go 原型脚本（如 `gbParse.py`、`PGA.go`），仅作留档参考。
* `EvoCG1.0/`: 连用的细菌核心基因组提取工具。
* `test/` 与 `input/`: 测试数据集及配置文件与参比 PG 集（如 `26_PG.txt`, `nA_AG.gbk`）。

### 编译说明
新版的 `Evopga` 由纯 Go 语言编写，没有其他第三方包依赖，非常轻量化。您可以在任意系统重新编译：
```bash
cd Evopga
go build -o ../bin/Evopga .
```

### 使用指南
打包好的二进制程序 `Evopga` 位于 `bin/` 下，有两个子命令：

1. **parse** (解析): 从 GenBank 文件提取 CDS, rRNA, tRNA, ncRNA, misc_feature 等元件以供后续使用。
   ```bash
   ./bin/Evopga parse ./test/nA_AG.gbk > output_PGAG.tab.txt
   ```

2. **annotate** (注释): 利用相互最佳序列比对，赋予查询基因表对应的 PG UID。
   ```bash
   ./bin/Evopga annotate -pg ./test/26_PG.txt -tab ./output_PGAG.tab.txt -mode PGAG -strain nA_AG -mutbestDir ./EvoCG1.0/result/out_mutbest_filt > nA_AG_annotated.tab.txt
   ```

如需自动执行全流程管线，只需在根目录执行下面这条命令（替换为您自己的输入路径即可自动调用内建的工作流及CG核心组工具）：
```bash
./bin/Evopga pipeline -gbk ./test/nA_AG.gbk -pg ./test/26_PG.txt -seq ./test/seq -out ./output
```