#Modl1虽然可以批量循环，但是会损失一定的数据，具体参与运行的菌株请查看temp/AMR_list/

#对输入文件做一个过滤，规则如下：
#1.input的list、prot_file、nuc_file三者ID取交集，复制prot_file、nuc_file到temp问文件夹
#2.对temp文件夹的prot_file：计算prot_file文件大小均值，把所有小于均值90%大小的文件放弃，相应nuc_file中也删除，记录相应ID
#3.利用此时的temp文件夹的prot_file的ID，更新Input的list，并将该总list拆成每个子list存在temp/AMR_list文件夹中
#同时，对一个耐药分类进行计数统计：小于5，则提示不再运行程序（因为5倍交叉验证要求至少5个样本），反之则利用新的AMR_list.txt、prot_file、nuc_file运行程序
#输出计数统计结果


import shutil
import function_common as func
import os

#文件夹检查 or 新建
#output部分
func.check_dir(r"output/CG_results/")


#输入路径
input_prot_file = r"input/seq/prot_file/"


#2.获取prot文件夹内所有文件的文件大小平均值，小于平均值90%的文件舍去
total_size = 0  # 文件夹内所有文件的总大小
file_count = 0  # 文件夹内文件的数量

# 遍历文件夹中的所有文件并计算它们的大小
for dirpath, dirnames, filenames in os.walk(input_prot_file):
    for filename in filenames:
        file_path = os.path.join(dirpath, filename)
        try:
            # 获取文件大小并将其添加到总大小中
            size = os.path.getsize(file_path)
            total_size = total_size + size
            file_count = file_count + 1
        except OSError:
            # 如果无法获取文件大小，则跳过该文件
            pass
# 计算平均文件大小
if file_count > 0:
    average_size = total_size / file_count
else:
    average_size = 0
cutoff_size = average_size * 0.9 #分割线文件大小

#记录所有文件大小小于cutoff_size的文件id
small_file_id_list = []
for dirpath_new, dirnames_new, filenames_new in os.walk(input_prot_file):
    for filename_new in filenames_new:
        file_path_new = os.path.join(dirpath_new, filename_new)
        size_new = os.path.getsize(file_path_new)
        if size_new < cutoff_size:
            small_file_id = filename_new.replace(".fasta","")
            small_file_id_list.append(small_file_id) #记录错误的id
            os.remove(file_path_new) #删除错误的文件(蛋白文件)


#输出小于cutoff_size的文件id列表，输出为txt
with open(r"remove_id.txt", "w") as f:
    for i in small_file_id_list:
        f.write(i + "\n")

            

