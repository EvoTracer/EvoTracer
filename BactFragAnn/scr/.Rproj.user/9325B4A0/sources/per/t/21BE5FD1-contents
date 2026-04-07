library(tidyverse)
library(clusterProfiler)
BiocManager::install("org.EcK12.eg.db")
library(org.EcK12.eg.db) # 大肠杆菌 K-12 专用包
library(KEGGREST)

# =================================================
# 1. 构建 GO 背景 (利用 org.EcK12.eg.db)
# =================================================
# 提取 GO -> SYMBOL 的对应关系
ec_go_bg <- AnnotationDbi::select(org.EcK12.eg.db, 
                                  keys = keys(org.EcK12.eg.db, keytype="SYMBOL"), 
                                  columns = c("GO", "SYMBOL"), 
                                  keytype = "SYMBOL") %>%
  filter(ONTOLOGY == "BP") %>% 
  dplyr::select(GO_ID = GO, Gene_Symbol = SYMBOL) %>% # 使用 dplyr::select 重命名列
  filter(!is.na(GO_ID) & !is.na(Gene_Symbol)) %>%
  distinct()


write.csv(ec_go_bg, "/home/sunleting/storage/202602214-ag_function_enrich/input/Ecoli_GO_Background.csv", row.names=F)


# =================================================
# 2. 构建 KEGG 背景 (利用 KEGGREST 抓取 eco)
# 获取 eco (E. coli K-12 MG1655) 的通路信息
# eco_link <- keggLink("pathway", "eco")
# eco_genes <- keggList("eco")

# eco_kegg_bg <- data.frame(GeneID = names(eco_link), PathwayID = eco_link) %>%
#   left_join(data.frame(GeneID = names(eco_genes), Raw_Desc = eco_genes), by="GeneID") %>%
#   mutate(Gene_Symbol = str_extract(Raw_Desc, "^[^;]+")) %>% # 提取基因名
#   mutate(PathwayID = str_remove(PathwayID, "path:")) %>%
#   select(PathwayID, Gene_Symbol) %>%
#   distinct()
# write.csv(eco_kegg_bg, "/home/sunleting/storage/202602214-ag_function_enrich/input/Ecoli_KEGG_Background.csv", row.names=F)