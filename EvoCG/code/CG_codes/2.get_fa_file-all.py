
#对每一行CG家族蛋白，建一个文件，取第一个菌株的蛋白名作为文件名，后缀为.fa，
#在这个文件里储存113个菌株该家族的各自蛋白序列，注释行直接改为相应菌株名。

import os
from itertools import islice
from operator import add

#输入
input_dir = r"output/CG_results/1.cd-hit_output/1.cd-hit_fatsa/" #输入的 经过去冗余的 faa蛋白质文件
output_file = r"output/CG_results/all-strain-together/2.all_prot-sequence.fa" #输出的合并后的总蛋白质文件
output_dir = r"output/CG_results/all-strain-together/2.result/" #输出的 存储fa文件集合的文件夹
CG_file = r"output/CG_ALL.txt" #输入的 CG.tab.txt文件
rename_CG_file = r"output/CG_results/all-strain-together/2.rename_CG.tab.txt" #输出的 重命名之后的CG.tab.txt文件

os.mkdir(r'output/CG_results/all-strain-together')

#1.重命名CG.tab.txt文件
line_output_list = []  #存储重命名后的信息
with open(CG_file) as input_file:
    #处理第一行（即索引行）
    first_line = next(input_file) #取出第一行（即索引行）
    first_line_list = first_line.replace("\n","").split("\t") #消除换行符 \t分隔 
    index_list = []  
    for pp1 in first_line_list:
        index_list.append(pp1 + "---")
    #处理第二行以后的内容（即索引行）
    for pp2 in input_file:
        prot_list = pp2.replace("\n","").split("\t") #消除换行符 \t分隔 
        output_line = list(map(add,index_list,prot_list)) #CG文件的索引和菌名结合 (菌株名-原始蛋白ID)
        line_output_list.append(output_line)
#输出重命名之后的CG.tab.txt文件
str_new = "\t"
with open(rename_CG_file,"a+",encoding="utf-8") as t:
    #输出
    t.write(first_line) #表头
    for qq in line_output_list:
        line_out = str_new.join(qq)
        t.write(line_out + "\n")


#2.蛋白ID重命名 并合并为一个总文件

seq_dict = {} #存储序列字典表
for root, dirs, files in os.walk(input_dir):#文件夹内文件批量循环   
    for name in files:
       
       index_name = name.split(".fasta")[0] #取出文件名
       file_path = os.path.join(root, name) #获取文件路径

       with open(file_path) as p:
           for line in p:
               if ">" in line:
                   seq_name = line.replace("\n","").split(">")[1]
                   seq_name = ">" + index_name + "---" + seq_name  #构造蛋白序列名字(即：>菌株名|原始蛋白ID)
                   seq_dict[seq_name] = ""
               elif line == "\n":
                   pass
               else:
                   seq_dict[seq_name] = seq_dict[seq_name] + line.replace("\n","") #成功读取fasta 文件
#输出all_prot-sequence.fa
with open(output_file,"a+",encoding="utf-8") as all_seq:
    for key_1,value_1 in seq_dict.items():
        str_all_seq = key_1 + "\n" + value_1 + "\n"
        all_seq.write(str_all_seq)



#3.开始扫描 重命名之后的CG.tab.txt文件，生成fa文件集合

#先生成一个文件夹result 存储fa文件集合
try:
    os.makedirs(output_dir)
except FileExistsError:
    pass

#开始读取 重命名之后的CG.tab.txt文件
input_file_new_CG = open(rename_CG_file) 
for line_new_CG in islice(input_file_new_CG, 1, None):  #跳过第一行表头
   #获取fa文件名
   line_new_CG = line_new_CG.replace("\n","") #消除换行符
   line_new_CG_list = line_new_CG.split("\t") #按照分隔符划分
   output_file_name = line_new_CG_list[0] #第一个名字作为fa文件名
    #扫描重命名之后的CG.tab.txt文件每一行
   for element in line_new_CG_list: 

      key_2 = ">" + element
      value_2 = seq_dict[key_2]
      #输出该蛋白质文件
      with open(output_dir + output_file_name + ".fa","a+",encoding="utf-8") as t:
         t.write(key_2 + "\n")
         t.write(value_2 + "\n")
