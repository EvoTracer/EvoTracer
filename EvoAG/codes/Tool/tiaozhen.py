import sys

input_file = sys.argv[1]
output_file = sys.argv[2]
param1 = sys.argv[3]
param2 = sys.argv[4]
param3 = sys.argv[5]

def process_file(input_file, output_file, param1, param2, param3):
    lines = []

    # 读取输入文件内容
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # 检查是否存在包含特定表头的行，如果存在则移除
    lines = [line for line in lines if line.strip() != "Source_Strain\tSource_st\tSource_end\tAnnotation\tHom_strain"]

    if len(lines) < 1:
        print("输入文件内容不符合预期格式")
        return

    # 处理每一行数据
    processed_lines = []
    for line in lines:
        parts = line.strip().split('\t')

        if len(parts) < 5:
            continue

        # 如果第一列是数字，并且既不等于 param1 也不等于 param3，则进行替换
        try:
            if parts[0].isdigit() and parts[0] != param1 and parts[0] != param3:
                # 构建新的行内容
                new_line = f"{param1}\t{parts[0]}\t{parts[1]}\t{param2}\t{param2}_added"
                processed_lines.append(new_line)
            else:
                processed_lines.append('\t'.join(parts))
        except ValueError:
            processed_lines.append('\t'.join(parts))

    # 写入到输出文件
    with open(output_file, 'w') as f:
        f.write("Source_Strain\tSource_st\tSource_end\tAnnotation\tHom_strain\n")  # 重新写入表头
        for line in processed_lines:
            f.write(line + '\n')

if __name__ == "__main__":
    # 检查命令行参数是否符合预期
    if len(sys.argv) != 6:
        print("Usage: python script.py <input_file> <output_file> <param1> <param2> <param3>")
        sys.exit(1)

    # 调用函数处理文件
    process_file(input_file, output_file, param1, param2, param3)
