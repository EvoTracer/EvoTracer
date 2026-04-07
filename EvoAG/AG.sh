#!/bin/bash
#此为BactAG的总控脚本，可进行从输入到输出的全流程操作
#设置线程数
threads=24

# 设置要检查的目录路径
IDfile="codes/bactID.txt" 
AGoutput="BactAG-output$(date '+%Y-%m-%d %H:%M:%S')"
output_dir="output"
log_dir="log"
# 创建日志目录
if [ -d "$log_dir" ]; then
    echo "Directory exists: $log_dir. Deleting it..."
    rm -rf "$log_dir"
    echo "Directory deleted: $log_dir"
fi

# 创建新的 output 文件夹
echo "Creating a new directory: $log_dir"
mkdir "$log_dir"


# 检查文件是否存在
if [ -f "$IDfile" ]; then
    echo "File exists: $IDfile. Deleting it..."
    rm "$IDfile"
    echo "File deleted: $IDfile"
else
    echo "File does not exist: $IDfile"
fi

# 创建一个新的文件并写入默认内容
echo "Creating a new file: $IDfile"
echo "# This is BactID  $(date '+%Y-%m-%d %H:%M:%S')" > "$IDfile"

# 检查 output 文件夹是否存在
if [ -d "$output_dir" ]; then
    echo "Directory exists: $output_dir. Deleting it..."
    rm -rf "$output_dir"
    echo "Directory deleted: $output_dir"
fi

# 创建新的 output 文件夹
echo "Creating a new directory: $output_dir"
mkdir "$output_dir"

#运行主程序

python codes/MutiAG_plus.py $threads

# 运行完毕，输出结果

mkdir "$AGoutput"
cp -r "$output_dir" "$AGoutput"
cp "$IDfile" "$AGoutput"
cp -r input/tree_file "$AGoutput"
cp -r "$log_dir" "$AGoutput"

rm -rf "$output_dir"
rm -rf "$IDfile"
rm -rf "$log_dir"
echo "All done! All Output files are in $AGoutput"