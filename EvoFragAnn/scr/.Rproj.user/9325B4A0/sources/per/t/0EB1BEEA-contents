# 加载必要的R包
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(stringr)
library(dplyr)

# ==========================================
# 设置您的 Windows 基础工作路径 (注意这里用的是正斜杠 /)
base_dir <- "/Users/vas/工作数据存放/202602214-ag_function_enrich"
# ==========================================

# 1. 读取手动准备好的背景基因集
# 拼接完整的输入文件路径
bg_file <- file.path(base_dir, "input", "LT2_background_TERM2GENE.csv")

gene_GO <- read.csv(bg_file,
                    header = TRUE,
                    stringsAsFactors = FALSE)

# 规范化列名
colnames(gene_GO) <- c("term","gene")

# 2. 准备GO注释和映射文件
term2gene <- gene_GO[, c("term", "gene")]  
go2term_df <- clusterProfiler::go2term(term2gene$term) 
go2ont <- go2ont(term2gene$term)

# 3. 设置输入和输出目录
input_dir <- file.path(base_dir, "input", "sal_1")
file_list <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

# 创建主输出目录
main_output_dir <- file.path(base_dir, "output", "sal_go", "np")

dir.create(main_output_dir, showWarnings = FALSE, recursive = TRUE)

# 4. 循环处理每个基因列表文件
for (file in file_list) { 
  file_base <- tools::file_path_sans_ext(basename(file))
  accessory_genes <- read.csv(file)$gene
  
  # 5. 执行GO富集分析
  accessory_go_rich <- enricher(gene = accessory_genes,
                                TERM2GENE = term2gene,
                                TERM2NAME = go2term_df,
                                pAdjustMethod = 'BH',
                                pvalueCutoff = 1,
                                qvalueCutoff = 1)
  
  # 6. 检查并统计匹配情况
  matched_genes <- accessory_genes[accessory_genes %in% gene_GO$gene]
  
  matching_stats <- data.frame(
    Metric = c("Input file", "Total target genes", "Matched genes in background", "Matching percentage"),
    Value = c(
      basename(file),
      length(accessory_genes),
      length(matched_genes),
      paste0(round(length(matched_genes)/length(accessory_genes)*100, 1), "%")
    )
  )
  
  matched_gene_list <- data.frame(Matched_Genes = matched_genes)
  
  # 为当前文件创建专属输出目录
  file_output_dir <- file.path(main_output_dir, file_base)
  dir.create(file_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 7. 输出统计信息和全量表格
  stats_file <- file.path(file_output_dir, paste0(file_base, "_matching_stats.txt"))
  write.table(matching_stats, stats_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  genes_file <- file.path(file_output_dir, paste0(file_base, "_matched_genes.txt"))
  write.table(matched_gene_list, genes_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("\n处理文件:", basename(file), "\n")
  cat("匹配比例:", round(length(matched_genes)/length(accessory_genes)*100, 1), "%\n")
  
  accessory_GO <- accessory_go_rich@result
  accessory_GO <- merge(accessory_GO, go2ont, by.x = 'ID', by.y = "go_id", all = FALSE)
  
  out_file <- file.path(file_output_dir, paste0(file_base, "_go_enrich.txt"))
  write.table(accessory_GO, out_file, sep = '\t', row.names = FALSE, quote = FALSE)
  
  # 8. 绘图部分 (只保留 p.adjust < 0.05)
  accessory_go_rich_filtered <- dplyr::filter(accessory_go_rich, p.adjust < 0.05)
  
  if (nrow(accessory_go_rich_filtered@result) > 0) {
    bar_plot_all <- barplot(accessory_go_rich_filtered,
                            showCategory = 13,
                            color = "p.adjust",
                            title = "GO enrichment barplot (p.adjust < 0.05)") +
      scale_color_gradient(high = "#D53142", low = "#345097")
    
    ggsave(filename = file.path(file_output_dir, paste0(file_base, "_go_barplot.pdf")), 
           plot = bar_plot_all,
           width = 8, height = 6)
    
    cat("绘制并保存柱状图成功。\n")
  } else {
    cat(">>> 提示: 文件", basename(file), "中没有 p.adjust < 0.05 的显著富集条目，已跳过绘图。\n")
  }
}