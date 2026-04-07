#主要目的是将相同位点相同氨基酸删除后，依次拼接序列，获得SNP文件
#但是从 先合并序列（out.meg）再完成同位点删除 计算量太大
#因此改变策略 先处理每个mega文件的同位点删除（生成new_mega），再依次合并 减少计算量
#在这里，先获取全部蛋白的全mega文件


import os
from itertools import islice

input_dir_old_mega = r"output/CG_results/all-strain-together/3.mega/" #输入的 多重序列比对后的mega文件夹
output_dir_new_mega = r"output/CG_results/all-strain-together/4.SNP_mega/" #输出的 存储同位点删除后的mega文件夹
output_new_mega_file = r"output/CG_results/all-strain-together/4.SNP_mage.meg" #输出的 拼接后的mega文件
output_log_file = r"output/CG_results/all-strain-together/4.log_order_len.txt" #输出的 日志信息：拼接的顺序和拼接的长度


#现在进行第一步 先处理每个mega文件的同位点删除（生成new_mega）
for root, dirs, files in os.walk(input_dir_old_mega):#文件夹内文件批量循环   
    for name in files:

        #以每个文件为单位 开始处理
        file_path = os.path.join(root, name)

        # 构造文件字典
        prot_dict = {}
        with open(file_path) as p1:
            for line_1 in islice(p1,4,None):
                if line_1.startswith("#"):
                    faa_name_1 = line_1.replace("\n","")
                    prot_dict[faa_name_1] = ""
                else:
                    line_1 = line_1.replace("\n","")
                    prot_dict[faa_name_1] = prot_dict[faa_name_1] + line_1

        #开始匹配相同位置的氨基酸 并 记录删除的位点信息
        line_one = list(prot_dict.values())[0] #文件第一条序列
        seq_len = len(line_one) #获取序列初始长度

        del_num_list = [] #存储需要删除的位点坐标信息
        del_num = 0
        for i in line_one:  #依次读取序列的每个碱基 去与其它序列匹配
            for key1,value1 in prot_dict.items():
                aa_false = value1[del_num]  #每个序列 该位置的碱基
                if i == aa_false:
                    ideol = True
                else:
                    ideol = False
                    break
            if ideol == True:  #说明该位点序列完全相同 记录删除信息
                del_num_list.append(del_num)
            del_num = del_num + 1

        #开始删除del_num_list 中记录的位置
        #删除会导致移码，所以选择将不需要删除的序列 保留需要的序列
        prot_dict_new = {}
        for key2,value2 in prot_dict.items():
            str_seq = ""
            for i in range( len(value2) ):
                if i not in del_num_list:
                    str_seq = str_seq + value2[i]
            prot_dict_new[key2] = str_seq

        #创建输出文件夹  存储同位点删除后的mega文件
        try:
            os.makedirs(output_dir_new_mega)
        except FileExistsError:
            pass
        
        #保存结果
        with open(output_dir_new_mega + "new-" + name,"a+",encoding="utf-8") as output_line:
            
            #表头
            output_line.write("#mega\n")
            output_line.write("!Title CG;\n")
            output_line.write("!Format DataType=Protein indel=-;\n")
            output_line.write("\n")

            for key_3,value_3 in prot_dict_new.items():
                output_line.write(key_3 + "\n")
                output_line.write(value_3 + "\n") 
                output_line.write("\n")



#现在进行第二步 对所有new_mega文件依次合并（只要属于同一菌株就合并）同时，新mega文件里面序列ID改为菌株名
#另外还要输出拼接顺序和new_mega文件序列长度

prot_sequence_dict = {}  #存储最终的拼接后的mega信息
prot_log_dict = {} #存储日志信息：拼接的顺序和拼接的长度

num = 0
for root, dirs, files in os.walk(output_dir_new_mega):#在同位点删除后的new_mega文件夹内文件批量循环   
    for name in files:

        num = num + 1

        #每个文件读取
        file_path = os.path.join(root, name)

        #第一个文件
        if num == 1:
            with open(file_path) as p1:
                for line_1 in islice(p1,4,None):
                    if line_1.startswith("#"):
                        faa_name_1 = line_1.replace("\n","").split("---")[0]
                        prot_sequence_dict[faa_name_1] = ""
                    else:
                        line_1 = line_1.replace("\n","")
                        prot_sequence_dict[faa_name_1] = prot_sequence_dict[faa_name_1] + line_1
                        seq_1 = prot_sequence_dict[faa_name_1]
                        log_num = len(seq_1)
        #第二个文件开始进行拼接 
        else:
            with open(file_path) as p2:
                for line_2 in islice(p2,4,None):
                    if line_2.startswith("#"):
                        faa_name_2 = line_2.replace("\n","").split("---")[0]
                        prot_sequnece = prot_sequence_dict[faa_name_2]
                        len_0 = len(prot_sequnece)
                    else:
                        add_seqence = line_2.replace("\n","") 
                        prot_sequnece = prot_sequnece + add_seqence
                        prot_sequence_dict[faa_name_2] = prot_sequnece
                        len_1 = len(prot_sequnece)
        
            #获取此次拼接进来的文件的序列长度
            log_num = len_1 - len_0

        #记录拼接进来的文件名和序列长度   
        prot_log_dict[name] = log_num


#输出 拼接后的mega文件
with open(output_new_mega_file,"a+",encoding="utf-8") as output_line:
    
    #表头
    output_line.write("#mega\n")
    output_line.write("!Title CG;\n")
    output_line.write("!Format DataType=Protein indel=-;\n")
    output_line.write("\n")

    for key_2,value_2 in prot_sequence_dict.items():
        output_line.write(key_2 + "\n")
        output_line.write(value_2 + "\n") 
        output_line.write("\n")


#输出日志信息：拼接的顺序和拼接的长度
with open(output_log_file,"a+",encoding="utf-8") as output_line_3:
    for key_3,value_3 in prot_log_dict.items():
        str_3 = str(key_3) + "\t" + str(value_3) + "\n"
        output_line_3.write(str_3)