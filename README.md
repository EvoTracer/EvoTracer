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

## 附录：硬件配置要求参考

   最低配置：4 CPU 核心, 16 GB 内存（适用于小批量菌株）。

   推荐配置：16+ CPU 核心, 64+ GB 内存（适用于大批量菌株）。

   存储空间：视输入基因组队列文件大小和数量而定，建议预留至少 50 GB 以上的空闲读写空间以应对临时比对文件的生成。
