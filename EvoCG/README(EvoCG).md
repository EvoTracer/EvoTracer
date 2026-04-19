# BactCG

BactCG 是一个用于基因组和核心基因（Core Gene）分析的快速处理管线工具包。本项目将原先分散的 Python 脚本、Perl 脚本以及 BactCG 2.0 的各组件模块统一重构、编译并打包成了一个独立的高性能 Go 命令行程序。由于其使用了通用的相对路径体系，它非常适合分发装箱，在不同设备间可以零修改无缝迁移。

## 环境依赖 (Dependencies)

在运行完整流水线之前，请确保系统中已正确安装以下外部工具：

- **NCBI BLAST+** (`blastn` / `blastp` 需要在系统环境变量 `$PATH` 中)
- **CD-HIT** (默认使用 `$PATH` 中的 `cd-hit`，若无则尝试检测当前相对路径：`Tool/cd-hit-v4.6.7-2017-0501/cd-hit`)
- **Clustalw2** (默认使用 `$PATH` 中的 `clustalw2`，若无则尝试检测当前相对路径：`Tool/clustalw2`)

## 编译与安装

如果你需要重新编译本项目，只需进入项目目录并执行：

```bash
go mod tidy
go build -o bactcg .
```

这将在当前目录生成一个名为 `bactcg` 的可执行文件。

## 使用说明 (Usage)

基本命令格式：

```bash
./bactcg [command] [flags]
```

### 1. 运行核心自动化流水线 (替代 1.cd-hit_and_BactCG.sh)

`run` 命令是整个工具的核心启动项，它会自动检查环境变量和相对路径依赖，并交互式询问用户是否开启体积累积过滤(QC)模块。另外，你能借由参数覆盖原先硬编码的一系列分析指标，如 CD-HIT 阈值或参考菌株：

```bash
# 最简执行（采用默认参数运行）
./bactcg run 

# 携带参数运行示例：
./bactcg run \
  --ref "BMU_04865" \
  --cd-c 0.7 \
  --cd-s 0.7 \
  --cg1 0.8 \
  --cg2 0.9 \
  -i "input/seq/prot_file/" \
  -o "output/CG_results/"
```

**可用的流水线参数：**
* `-r, --ref`：BactCG 指定的参考菌株名字
* `--cd-c`：cd-hit `-c` 阈值参数（默认 `"0.7"`）
* `--cd-s`：cd-hit `-s` 阈值参数（默认 `"0.7"`）
* `--cg1`：BactCG 第一个参数（默认 `"0.8"`）
* `--cg2`：BactCG 第二个参数（默认 `"0.9"`）
* `-i, --input`：输入的蛋白质全量 `fasta` 文件夹相对路径
* `-o, --output`：流水线统一的结果输出目录相对路径

### 2. 序列与数据清洗命令

替代了早期的 Python 清洗和质控脚本。

- **清理格式空行和泛头 (对应 `del.py` 与 `remove_empty_lines.py`)**
  ```bash
  ./bactcg clean -i ./dataset/1
  ```

- **体积质控与过滤 (对应 `0-1.filterate_file_M1.py`)**
  计算指定目录中文件的平均大小，剔除并删除低于平均体积 90% 的蛋白文件，将剔除 ID 保存在 `remove_id.txt`。
  ```bash
  ./bactcg filter -i input/seq/prot_file/
  ```

- **多序列比对与格式转换 (对应 `clustal.pl` 与 `gcg2meg.pl`)**
  对指定文件夹中的 `.fa` 文件进行 Clustalw2 比对生成 `.gcg` 格式件，并将其转换为对应的 `.meg` (MEGA) 文件。
  ```bash
  ./bactcg clustal -d "output/CG_results/2.result"
  ```

### 3. 原生 BactCG 2.0 内置算法模块

所有原版 BactCG 2.0 的功能也已被完美内部整合。您可以以子命令的形式直接调用它们，且不用再理会任何路径冲突问题。

```bash
./bactcg cg            # 运行 BactCG2.0 主算法逻辑
./bactcg bestpicker    # 运行 BestPicker 筛选逻辑
./bactcg batchprot     # 批量处理 ProtAcc 到 locTag
./bactcg combsingle    # 运行 combSingleGeneOrtho
./bactcg lenext        # 运行长度扩展逻辑 (lenExt)
./bactcg mutbest       # 运行 mutbest 处理
./bactcg protacc       # 处理单文件 ProtAcc 到 locTag
./bactcg simcov        # 运行相似度同源过滤 (simcovfilter)
```

## 帮助信息

若想查询任意命令的详细参数或帮助信息，请追加 `-h` 或 `--help`：

```bash
./bactcg --help
./bactcg filter --help
```
