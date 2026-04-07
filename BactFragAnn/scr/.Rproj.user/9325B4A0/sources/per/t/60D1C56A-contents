getwd()
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(stringr)

#读取手动准备好的背景基因集
gene_GO <- read.csv('/Users/vas/工作数据存放/202602214-ag_function_enrich/input/Ecoli_GO_Background.csv',
                    header = TRUE,
                    stringsAsFactors = FALSE)
colnames(gene_GO)
colnames(gene_GO) <- c("term","gene")
colnames(gene_GO)
#gene_GO <- gene_GO[!duplicated(gene_GO$gene), ]

#为直接注释补充为间接注释
#term2gene <- buildGOmap(gene_GO)
term2gene <- gene_GO[, c("term", "gene")]  # 确保顺序：term在前，gene在后
go2term_df <- clusterProfiler::go2term(term2gene$term)#获取GO术语名称 #将GoId转换为GoTerm
#将GoId转换为GoOnt
go2ont <- go2ont(term2gene$term)

input_dir <- "/Users/vas/工作数据存放/202602214-ag_function_enrich/input/ecoli_region"
# file_list <- list.files(input_dir, pattern = "\\.xlsx$", full.names = TRUE)
file_list <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

# 主输出目录
main_output_dir <- "/Users/vas/工作数据存放/202602214-ag_function_enrich/output/ecoli_go/region"
dir.create(main_output_dir, showWarnings = FALSE, recursive = TRUE)

#读取基因列表文件中的基因名称
for (file in file_list) { 
  file_base<- tools::file_path_sans_ext(basename(file))
   accessory_genes <- read.csv(file)$Gene
  #accessory_genes <- read.csv(file)$REFSEQ
  
  #GO 富集分析
  accessory_go_rich <- enricher(gene = accessory_genes,  #待富集的基因列表
                                TERM2GENE = term2gene,  #背景基因集
                                TERM2NAME = go2term_df, 
                                pAdjustMethod = 'BH',  #指定 p 值校正方法
                                pvalueCutoff = 1,  #指定 p 值阈值（可指定 1 以输出全部）
                                qvalueCutoff = 1)  #指定 q 值阈值（可指定 1 以输出全部）
  
  head(accessory_go_rich)
  # 检查目标基因与背景集的匹配情况
  matched_genes <- accessory_genes[accessory_genes %in% gene_GO$gene]
  # cat("目标基因总数:", length(accessory_genes), "\n")
  # cat("在背景集中匹配的基因数:", length(matched_genes), "\n")
  # cat("匹配比例:", round(length(matched_genes)/length(accessory_genes)*100, 1), "%\n")
  # 
  
  # 创建匹配结果数据框
  matching_stats <- data.frame(
    Metric = c("Input file", "Total target genes", "Matched genes in background", "Matching percentage"),
    Value = c(
      basename(file),
      length(accessory_genes),
      length(matched_genes),
      paste0(round(length(matched_genes)/length(accessory_genes)*100, 1), "%")
    )
  )
  
  # 创建匹配基因列表
  matched_gene_list <- data.frame(Matched_Genes = matched_genes)
  
  # 为当前文件创建专属输出目录
  file_output_dir <- file.path(main_output_dir, file_base)
  dir.create(file_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 输出匹配结果
  stats_file <- file.path(file_output_dir, paste0(file_base, "_matching_stats.txt"))
  write.table(matching_stats, stats_file, sep = "\t", row.names = FALSE, quote = FALSE)
  genes_file <- file.path(file_output_dir, paste0(file_base, "_matched_genes.txt"))
  write.table(matched_gene_list, genes_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("\n处理文件:", basename(file), "\n")
  cat("目标基因总数:", length(accessory_genes), "\n")
  cat("在背景集中匹配的基因数:", length(matched_genes), "\n")
  cat("匹配比例:", round(length(matched_genes)/length(accessory_genes)*100, 1), "%\n")
  cat("匹配统计已保存至:", stats_file, "\n")
  cat("匹配基因列表已保存至:", genes_file, "\n")
  
  
  #再把 GO Ontology 信息添加在上述 GO 富集结果中
  accessory_GO <- accessory_go_rich@result
  accessory_GO <- merge(accessory_GO, go2ont, by.x = 'ID', by.y = "go_id", all = FALSE)
  
  
  #输出
  out_file <- file.path(file_output_dir, paste0(file_base, "_go_enrich.txt"))
  write.table(accessory_GO, out_file, 
              sep = '\t', row.names = FALSE, quote = FALSE)
  
  # --- 新增：专门过滤出 P < 0.05 的结果用于绘图 ---
  accessory_go_rich_plot <- accessory_go_rich
  # 提取 pvalue < 0.05 的行 (如果你想用校正后的P值，请将 pvalue 改为 p.adjust)
  accessory_go_rich_plot@result <- accessory_go_rich_plot@result[accessory_go_rich_plot@result$pvalue < 0.05, ]
  
  # 检查过滤后是否还有数据，如果没有则跳过当前文件的绘图，防止代码崩溃
  if (nrow(accessory_go_rich_plot@result) == 0) {
    cat("  -> 注意：该文件没有 pvalue < 0.05 的富集结果，已跳过绘图步骤。\n")
    next
  }
  ####以下为绘图
  
  # --- 新增：专门过滤出 P < 0.05 的结果用于绘图，并按显著性排序 ---
  accessory_go_rich_plot <- accessory_go_rich
  
  # 1. 提取出结果数据框
  res_df <- accessory_go_rich_plot@result
  
  # 2. 过滤 pvalue < 0.05 的行 (如果想更严格，可把 pvalue 改为 p.adjust)
  res_df <- res_df[res_df$pvalue < 0.05, ]
  
  # 3. 严格按 pvalue 从小到大排序 (显著性从高到低)
  res_df <- res_df[order(res_df$pvalue), ]
  
  # 4. 把排序和过滤后的数据重新赋给画图对象
  accessory_go_rich_plot@result <- res_df
  
  # 检查过滤后是否还有数据，如果没有则跳过当前文件的绘图，防止代码崩溃
  if (nrow(accessory_go_rich_plot@result) == 0) {
    cat("  -> 注意：该文件没有 pvalue < 0.05 的富集结果，已跳过绘图步骤。\n")
    next
  }
  # -----------------------------------------------
  
  # 1. Barplot (柱状图)
  # clusterProfiler 会按传入对象的顺序提取 top N (这里的 showCategory = 13)
  bar_plot_all <- barplot(accessory_go_rich_plot,
                          showCategory = 13,
                          color = "p.adjust",
                          title = "GO enrichment barplot") +
    scale_color_gradient(high="#D53142",low="#345097")
  
  ggsave(filename = file.path(file_output_dir, paste0(file_base, "_go_barplot.pdf")), 
         plot = bar_plot_all,
         width = 8, height = 6)
  
  # 2. Cnetplot (基因-概念网络图)
  # 默认展示前 5 个最显著的 term（可以自己加 showCategory = N 修改）
  cnet_plot <- cnetplot(accessory_go_rich_plot) 
  ggsave(filename=file.path(file_output_dir,paste0(file_base, "_go_cnetplot.pdf")),
         plot=cnet_plot,
         width = 11, 
         height = 9)
  
  # 3. Emapplot (富集功能网络图)
  # 必须先计算过滤+排序后各个Term之间的相似度
  accessory_go_rich_sim <- pairwise_termsim(accessory_go_rich_plot) 
  emap_plot <- emapplot(accessory_go_rich_sim) 
  
  ggsave(filename=file.path(file_output_dir,paste0(file_base, "_go_emapplot.pdf")),
         plot=emap_plot,
         width = 11, 
         height = 9)
  
  # 4. Dotplot (气泡图)
  # 默认通常按 GeneRatio 排序，通过强制 orderBy = "pvalue" 或因为前面已经排好序，确保最显著的在上面
  dot_plot_all <- dotplot(accessory_go_rich_plot, 
                          showCategory = 20, 
                          orderBy = "pvalue", # 强制按 pvalue 排序气泡图
                          title = "GO enrichment") +
    theme_minimal()
  
  ggsave(filename = file.path(file_output_dir, paste0(file_base, "_go_dotplot.pdf")),
         plot=dot_plot_all,
         width = 9, 
         height = 7)
}