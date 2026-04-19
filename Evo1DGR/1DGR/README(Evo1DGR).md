# OneDGR  运行指南

本文档提供了 OneDGR 程序的安装说明、命令行示例及参数说明。

## 环境要求与安装

- **环境要求**：Python >= 3.6（目前仅依赖标准库，无需额外安装第三方库）。
- **安装方法**：本工具已配置 `setup.py`，支持直接安装为本地包。在当前目录（即包含 `OneDGR` 文件夹的目录）下，运行以下命令进行安装：

```bash
pip install -e ./OneDGR
```

安装成功后，系统会注册全局命令 `onedgr`，可以直接进行调用。

## 运行命令

**方式一：使用全局快捷命令（推荐，需先完成安装）**

```bash
onedgr obYMMbN --Evoid /home/hcd/1DGR/EvoID.txt -f /home/hcd/AG/output -o MyResults -w 8
```

**方式二：通过 Python 脚本直接运行**

在使用此方式前，请确保您当前的工作目录位于包含 `OneDGR/run_onedgr.py` 的父目录下：

```bash
python OneDGR/run_onedgr.py obYMMbN --Evoid /home/hcd/1DGR/EvoID.txt -f /home/hcd/AG/output -o MyResults -w 8
```

## 参数说明

- `obYMMbN`:输入感兴趣的古代基因组ID
- `--Evoid`: EvoAG获得的古基因组生成关系（例如：`/home/hcd/1DGR/EvoID.txt`）。
- `-f`: 指定输入文件夹路径（例如：`/home/hcd/AG/output`）。
- `-o`: 指定输出结果保存的目录名称或前缀（此处为 `MyResults`）。
- `-w`: 指定运行的线程数或工作进程数（此处为 `8`，表示使用 8 个线程以加速计算）。

## 注意事项

- 运行此命令前，请确保您当前的工作目录位于包含 `OneDGR/run_onedgr.py` 的父目录下。
- 请检查 `--Evoid` 和 `-f` 提供的绝对路径在您的文件系统中是否真实存在且具有读取权限。
- 确保已配置好运行所需的 Python 环境及相关依赖。

---

# OneDGR User Guide

This document provides installation instructions, command-line examples, and parameter descriptions for the OneDGR program.

## Environmental Requirements & Installation

- **Requirements**: Python >= 3.6 (Currently depends only on the standard library, no third-party packages required).
- **Installation**: This tool is configured with a `setup.py` and can be installed as a local package. In the current directory (the directory containing the `OneDGR` folder), run the following command:

```bash
pip install -e ./OneDGR
```

Once successfully installed, the system will register a global command `onedgr`, which can be called directly.

## Run Commands

**Method 1: Using the global command (Recommended, requires installation first)**

```bash
onedgr obYMMbN --Evoid /home/hcd/1DGR/EvoID.txt -f /home/hcd/AG/output -o MyResults -w 8
```

**Method 2: Running the Python script directly**

Before using this method, ensure your current working directory is the parent directory containing `OneDGR/run_onedgr.py`:

```bash
python OneDGR/run_onedgr.py obYMMbN --Evoid /home/hcd/1DGR/EvoID.txt -f /home/hcd/AG/output -o MyResults -w 8
```

## Parameters 

- `obYMMbN`: Input the ancient genome ID of interest.
- `--Evoid`: The ancient genome generation relationship obtained by EvoAG (e.g., `/home/hcd/1DGR/EvoID.txt`).
- `-f`: Specify the input folder path (e.g., `/home/hcd/AG/output`).
- `-o`: Specify the output directory name or prefix (here it is `MyResults`).
- `-w`: Specify the number of threads or worker processes to run (here it is `8`, meaning 8 threads are used to accelerate computation).

## Notes

- Before running this command natively using Python, ensure your current working directory is the parent directory of `OneDGR/run_onedgr.py`.
- Please check that the absolute paths provided for `--Evoid` and `-f` actually exist in your file system and have read permissions.
- Ensure the required Python environment is properly configured.