#各个菌株内部：CD-HIT（蛋白） 去冗余+家族聚类 
#之后对去冗余之后的蛋白质跑BactCG
#此程序只对cd hit的阈值>40%生效


input_dir="input/seq/prot_file/" #输入的蛋白质文件夹
output_dir="output/CG_results/" #输出文件夹
cd_hit_path='Tool/cd-hit-v4.6.7-2017-0501/cd-hit' #cd-hit脚本路径

# BactCG_test="Tool/BactCG2.0/test" #BactCG 输入文件路径
BactCG_tool="Tool/BactCG2.0/CG"  #BactCG 脚本路径
BactCG_fasta_file_name='REF' #BactCG 指定的参考菌株名字 （变）

#BactCG_fasta_file_name="BMU_04865" #BactCG 指定的参考菌株名字 （变）

cd_hit_c=0.7 #cd-hit -c的参数（变）
cd_hit_s=0.7 #cd-hit -s的参数（变）
BactCG_1=0.8 #BactCG第一个参数（变）
BactCG_2=0.9 #BactCG第二个参数（变）


for name in $(ls $input_dir)
do 	

#进行cd-hit处理
#创建输出文件夹
output_dir_finial=$output_dir"1.cd-hit_output"
mkdir -p $output_dir_finial;
#执行cd-hit脚本
$cd_hit_path -i $input_dir$name -o $output_dir_finial"/"$(basename -s .fasta $name)"_cd-hit.fasta" -c $cd_hit_c -M 0 -T 0 -s $cd_hit_s -d 1000000

#将输出的cd-hit fasta文件copy放到一个文件夹里面
#创建输出文件夹
output_dir_cd_hit=$output_dir_finial"/1.cd-hit_fatsa"
mkdir -p $output_dir_cd_hit;

#复制文件 移动一份到BactCG_test
old_faa_file_name=$output_dir_finial"/"$(basename -s .fasta $name)"_cd-hit.fasta"
new_faa_file_name=$output_dir_cd_hit"/"$(basename -s .fasta $name)".fasta"
cp $old_faa_file_name $new_faa_file_name
done;


#跑BactCG
Tool/BactCG2.0/CG output/CG_results/1.cd-hit_output/1.cd-hit_fatsa $BactCG_fasta_file_name $BactCG_1 $BactCG_2 ;
#取出CG.tab.txt文件
cp Tool/BactCG2.0/result/cg_result/CG.tab.txt  output/CG_ALL.txt

#删除BactCG相应文件
#删除result文件夹
rm -rf Tool/BactCG2.0/result


