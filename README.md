# EvoTracer 计算引擎部署与安装指南

EvoTracer 引擎支持在本地工作站或高性能计算集群（HPC）上通过命令行（CLI）运行，提供从序列比对、AOG 重建到动态进化图谱（Evo1DGR / EvolTraj）生成的端到端分析。

##  方案一：Docker 容器化部署（推荐）

使用 Docker 部署可以跳过所有依赖包的配置，是最为推荐、且保证环境一致性的方案（尤其适合无 root 权限或存在多版本冲突的服务器）。

### 1. 拉取官方镜像

在已安装 Docker 的终端中执行以下命令，获取最新的 EvoTracer 引擎镜像：


```bash
docker pull evotracer/evotracer:latest
```
------

### 2. 验证安装

检查镜像是否能够正常输出帮助文档：

```bash
docker run --rm evotracer/evotracer:latest evotracer --help
```
------

### 3. 投递任务示例

运行分析时，需要将宿主机（本地）的数据目录挂载（-v）到容器内部。例如，将本地的
/path/to/my_data 挂载到容器的 /data 目录：

```bash
docker run --rm -v /path/to/my_data:/data evotracer/evotracer:latest \
  evotracer pipeline \
  --input /data/input_genomes \
  --outdir /data/results \
  --module cg，ag，traj，pg，pga，1dgr \
  --threads 16
```
------

##  方案二：源码编译与安装部署（适用二次开发）

 

如果需要在集群的特定节点上运行，或者需要对 EvoTracer 的算法模块进行深度定制，可以通过源码配合 Conda 虚拟环境进行安装。

 

### 1. 系统环境准备

确保系统已安装 Git 以及 Miniconda3 或 Anaconda3以及mauve，Mega11。

------


### 2. 克隆源代码仓库

```bash
git clone https://github.com/EvoTracer /EvoTracer.git
cd EvoTracer
```
------



### 3. 创建并激活 Conda 隔离环境

EvoTracer 的核心环境依赖树已写入 environment.yml（包含必要的比对工具、和 Python/C++ 依赖）。

```bash
# 创建环境（自动安装依赖）
conda env create -f environment.yml
# 激活 EvoTracer 工作环境
conda activate evotracer_env
```
------


### 4. 完整分析测试运行

使用自带的 E. coli Demo 测试集验证安装

```bash
evotracer pipeline \
    --input ./example/Ecoli_demo_genomes.zip \
    --annotation ./example/Ecoli_annotation.gbk \
    --module cg，ag，traj，pg，pga，1dgr \
    --outdir ./test_output \
    --threads 8
```
------
## 进阶使用：独立分析与可视化工具

> **注意**：在运行以下脚本前，请确保您已激活 Conda 虚拟环境（`conda activate evotracer_env`）

### 工具一：EvolTraj工具
EvolTraj 专门用于描绘特定现存谱系或古老节点沿着进化轨迹的序列获得（Insertion）与丢失（Deletion）情况。它通过高精度的坐标锚定技术，帮助研究者直观地追踪基因组在演化尺度上的动态波动。


* **使用步骤：**
1. **修改脚本配置**：打开 `EvolTraj.py`，根据您的实际数据修改顶部「全局配置区域」的变量：
   - `ROOT_NODE`：根节点名称（需提供对应的 `.fasta` 文件）。
   - `EVOLUTION_TREE`：定义演化树结构的列表。
   - `MIN_EVENT_LENGTH`：结构变异的长度阈值（默认 1000）。
2. **准备文件**：确保工作目录下拥有相应的基因组 `.fasta` 文件与注释 `.gbk` 文件。
3. **运行脚本**：
```bash
   python EvolTraj.py
```
------
### 工具二：EvoCG工具
**EvoCG** 是用于快速处理基因组和核心基因（Core Gene）的分析管线。它将序列清洗（去除空行和非法字符）、cd-hit 聚类以及核心基因提取流程整合为高性能的 Go 命令行程序，适合处理大规模蛋白质序列集。

* **使用步骤**：
1. **数据准备**：将所有待分析的蛋白质全量 `fasta` 文件放入指定文件夹。
2. **配置参数**：确定参考菌株名称（`--ref`）以及 cd-hit 的相似度阈值（默认 `0.7`）。
3. **运行分析**：
```bash
./bactcg run --ref "Ref_Name" -i "./input/seq/" -o "./output/CG_results/"
```
* **输出结果**：生成核心基因集、比对结果序列以及经过格式清洗的标准化序列文件。

---

### 工具三：EvoPG工具
**EvoPG** 提供了优化的单文件版泛基因组分析管线。它无缝集成了序列计数、映射格式化、BLAST 过滤以及代表性序列生成等辅助流程，完全消除了对外部多个复杂子程序的依赖。

* **使用步骤**：
1. **数据准备**：确保输入目录包含所有的基因组蛋白质数据集（`.fasta` 或 `.faa`）。
2. **设置阈值**：确定分析所需的相似度阈值（例如 `0.7` 代表 70% 一致性）。
3. **运行分析**：
```bash
./BactPG <序列目录> <相似度阈值> <线程数>
```
* **输出结果**：输出包含 Pan-Genome 聚类矩阵、各基因簇的代表序列以及详细的序列映射统计表。

---

### 工具四：EvoAG工具
**EvoAG** 是一套用于细菌祖先基因组推断与重建的集成软件包。它通过高度并行的计算架构，显著提升了从系统发育树推演内部节点（Ancestral Nodes）序列的性能。

* **使用步骤**：
1. **准备输入**：系统发育树目录、各节点基因序列目录、菌株 ID 映射表 `bactID.txt`。
2. **运行分析**：
```bash
./BactAG -tree ./input/tree_file -gene ./input/gene -id ./codes/bactID.txt
```
* **输出结果**：生成重建后的古老节点基因组序列文件，及其推导出的进化事件（获得/丢失）原始记录。

---

### 工具五：EvoPGA工具
**EvoPGA** 能够解析 GenBank 文件，并根据相互最佳比对（Mutual Best Hit）结果，为每个基因自动添加对应的泛基因组簇 UID 注释。

* **使用步骤**：
1. **数据准备**：准备待注释的 `.gbk` 文件以及对应的泛基因组聚类结果（`.txt`）。
2. **运行自动化流水线**：
```bash
./bin/bactpga pipeline -gbk ./input/nA_AG.gbk -pg ./input/26_PG.txt -seq ./input/seq -out ./output
```
* **输出结果**：生成带有完整 PG UID 注释信息的基因表及更新后的 GenBank 注释文件。

---

### 工具六：Evo1DGR工具
**Evo1DGR** 专门用于生成细菌基因组的一维演化图谱。通过解析由 EvoAG 获得的古基因组生成关系，它能够精准定位特定的古基因组片段，为后续的可视化提供坐标支撑。

* **使用步骤**：
1. **环境配置**：在 `OneDGR` 目录下运行 `pip install -e .` 以激活全局命令 `onedgr`。
2. **运行命令**：
```bash
onedgr <Ancient_Node_ID> --bactid ./bactID.txt -f ./ag_output -o MyResults -w 8
```
* **输出结果**：生成包含指定节点一维映射坐标的 `.txt` 文件（这是 **EvoFragAnn** 的核心输入）。

---
### 工具七：EvoFragAnn工具

EvoFragAnn 利用 EvoPGA 的分析结果，对 Evo1DGR 和 EvolTraj 产生的序列片段进行自动化的功能注释。它能够将复杂的序列变异映射回基因组背景，生成带有交互功能的极坐标 Mosaic HTML 图谱，支持点击查看具体片段包含的基因详细信息。
**使用步骤：**

1. **数据准备**：在您的工作目录（即脚本中配置的 `BASE_DIR`）下准备两个文件夹：
   - `1DGR_en/`：存放 1DGR 流程生成的 `.txt` 映射结果文件。
   - `GBK_en/`：存放对应菌株的标准 `.gbk` 注释文件。
   
2. **修改脚本配置**：使用文本编辑器打开 `EvoFragAnn.py`，根据实际情况修改顶部「用户配置区域」的变量：
   - `BASE_DIR`：修改为您存放上述数据的基础工作目录绝对路径。
   - `NAME_MAP` 和 `RELATIONSHIPS`：根据您的实际序列命名习惯或原始 ID 映射表，更新节点名称和演化分支关系。
   
3. **运行分析**：在终端中执行以下命令运行脚本：
```bash
python EvoFragAnn.py
```
------

## 附录：硬件配置要求参考

   最低配置：4 CPU 核心, 16 GB 内存（适用于小批量菌株）。

   推荐配置：16+ CPU 核心, 64+ GB 内存（适用于大批量菌株）。

   存储空间：视输入基因组队列文件大小和数量而定，建议预留至少 50 GB 以上的空闲读写空间以应对临时比对文件的生成。