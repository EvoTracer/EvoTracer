library(tidyverse)
library(clusterProfiler)

# 1. 读取 GAF 文件
# 注意：quote = "" 防止文件中有奇奇怪怪的引号导致读取错误
gaf_data <- read.table("/home/sunleting/storage/202602214-ag_function_enrich/input/69.S_typhimurium_ATCC_700720.goa.txt", 
                       sep = "\t", comment.char = "!", quote = "", fill = TRUE, stringsAsFactors = FALSE)

# 2. 提取 TERM2GENE 关系表
# 【重要修正】：这里改用 V3 (Symbol)，而不是 V10 (Description)
lt2_term2gene <- gaf_data %>%
  select(GO_ID = V5, Gene_Symbol = V3) %>% 
  filter(Gene_Symbol != "" & GO_ID != "") %>% 
  distinct() # 去重

# 检查一下长什么样，应该是：
# GO:0007155  bcfE
head(lt2_term2gene)

write.csv(lt2_term2gene, 
            file = "/home/sunleting/storage/202602214-ag_function_enrich/input/LT2_background_TERM2GENE.csv", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = TRUE)