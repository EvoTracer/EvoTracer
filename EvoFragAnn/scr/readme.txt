1-bg_LT2_geneid4go.r 根据背景基因组的go文件，构建gene和go_term的对应关系。
go文件下载方式：https://geneontology.org/docs/download-go-annotations/
1-bg_ecoli_k12_geneid4go.r 因为ecoli k12 没有找到go文件，但是有在线数据库
library(org.EcK12.eg.db) # 大肠杆菌 K-12 专用包
所以根据在线包上的数据下载构建
2-ecoli_go_batch_file.R 批量进行go enrich，输入是文件夹，遍历文件夹下所有文件逐个处理
2-go_enrich.R 对单个输入基因文件进行处理

因为功能富集本质是超几何检验，输入数据要跟背景进行对应桥接。
发现输入数据和背景的基因名其实是一样的：
以“基因名”为桥梁
LT2 (GAF文件): 包含 GO term <--> 基因名 (如 bcfE) 的对应关系。
PGAP (GFF文件): 包含 你的ID (如 pgaptmp_000002) <--> 基因名 (如 arcA) 的对应关系。
连接: 通过 基因名 将两者串联，构建出 GO term <--> 你的ID 的背景文件。
    所以逻辑是1-构建背景：包含 GO term <--> 基因名 (如 bcfE) 的对应关系。
    输入文件就是一列gene（基因名 (如 bcfE) ）
    2-进行桥联，统计检验，功能富集

ps:之前做的操纵子组文章，都是输入菌株跟背景是同一个菌株，所以可以直接通过locustag进行对应
    