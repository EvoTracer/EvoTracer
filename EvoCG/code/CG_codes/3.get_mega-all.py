#对所有蛋白质文件（.fa）使用clustalw2进行对比，并格式化为mega文件输出
#多线程版本

import subprocess
from joblib import Parallel,delayed
import function_common as func
import os


#文件夹检查 or 新建
func.check_dir(r"output/CG_results/all-strain-together/3.gcc/")
func.check_dir(r"output/CG_results/all-strain-together/3.mega/")


input_fa_dir = r"output/CG_results/all-strain-together/2.result/"
output_gcc_dir = r"output/CG_results/all-strain-together/3.gcc/"
output_mage_fir = r"output/CG_results/all-strain-together/3.mega/"


#获取strain_list子集菌株名字
# AMR_strain_dict = func.get_AMR_strain_dict(r"temp/AMR_list/")
# strain_list = list(AMR_strain_dict.keys())


# for strain_name in strain_list:
#     path1 = input_fa_dir + strain_name + "/" #定位到菌株的fa文件夹
path1 = input_fa_dir
# #造输出路径
# func.check_dir(output_gcc_dir + strain_name + "/")
# func.check_dir(output_mage_fir + strain_name + "/")

#获取输入文件夹 input_fa_dir 的fa文件集合
fa_file_list = []
fa_files = os.listdir(path1) 
for i in fa_files:
    if i.endswith(".fa"):
        fa_file_list.append(i)

#定义函数：通过clustalw2获取mega文件
def get_mega(fa_file_name):

    name = fa_file_name.replace(".fa","") #获取去除.fa后缀的文件名

    input_fa_file = path1 + fa_file_name
    output_gcc_file = output_gcc_dir +  "/" + name + ".gcc"
    output_mega_file = output_mage_fir + "/"+ name + ".meg"

    subprocess.getstatusoutput(f"Tool/clustalw2 -infile={input_fa_file} -type=protein -output=gcg -outfile={output_gcc_file} -pwmatrix=GONNET -pwgapopen=10 -pwgapext=0.1 -gapopen=10 -gapext=0.2 -gapdist=4 -align")
    subprocess.getstatusoutput(f"perl Tool/gcg2meg.pl {output_gcc_file} {output_mega_file}")


# #套用多线程 (但是是乱序输出)
Parallel(n_jobs=-2)([delayed(get_mega)(fa_file_name) for fa_file_name in fa_file_list] ) 





