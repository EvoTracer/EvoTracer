import sys

def insert_into_file(file_A, file_B, file_C, file_D):
    lines_to_insert_before = []
    lines_to_insert_after = []
    
    # 读取file_B文件的第一行和最后一行数据
    with open(file_B, 'r') as f_b:
        first_line_b = f_b.readline().strip().split()
        first_column_b = int(first_line_b[0])

        f_b.seek(0)  # 重置文件指针到文件开头
        last_lines_b = f_b.readlines()[-1].strip().split()
        last_column_b = int(last_lines_b[1])

    # 处理file_A文件，根据逻辑处理行内容
    with open(file_A, 'r') as f_a:
        for line in f_a:
            columns = line.strip().split()
            first_column_a = int(columns[0])
            second_column_a = int(columns[1])
            
            if second_column_a > first_column_b and first_column_a < first_column_b:
                columns[1] = str(first_column_b)  # 修改第二列元素为B文件第一行的第一列元素
            
            modified_line = ' '.join(columns)
            
            if first_column_a < first_column_b:
                lines_to_insert_before.append(modified_line)
            elif first_column_a > last_column_b:
                lines_to_insert_after.append(modified_line)

    # 读取file_C文件的内容
    with open(file_C, 'r') as f_c:
        lines_c = f_c.readlines()
        
    # 将处理后的行和file_C文件的内容按顺序组合
    new_lines_d = lines_to_insert_before + lines_c + lines_to_insert_after

    # 将结果写入file_D文件
# 将结果写入文件D
    with open(file_D, 'w') as f_d:
        for line in new_lines_d:
            f_d.write('\t'.join(line.split()) + '\n')


if __name__ == "__main__":
    # 检查命令行参数是否符合预期
    if len(sys.argv) != 5:
        print("Usage: python script.py <file_A> <file_B> <file_C> <file_D>")
        sys.exit(1)

    # 获取命令行参数
    file_B = sys.argv[1]
    file_A = sys.argv[2]
    file_C = sys.argv[3]
    file_D = sys.argv[4]

    # 调用函数处理文件
    insert_into_file(file_A, file_B, file_C, file_D)