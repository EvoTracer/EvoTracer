getwd()
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(stringr)

#读取手动准备好的背景基因集
gene_GO <- read.csv('input/LT2_background_TERM2GENE.csv',
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
#将GoId转换为GoTerm
#go2term <- go2term(term2gene$GO)
# 将GoId转换为GoOnt
go2ont <- go2ont(term2gene$term)

#读取基因列表文件中的基因名称
accessory_genes <- read.delim('/home/sunleting/storage/202602214-ag_function_enrich/input/sal_1/Sal_NP_List.csv', 
                              stringsAsFactors = FALSE)$Gene

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
cat("目标基因总数:", length(accessory_genes), "\n")
cat("在背景集中匹配的基因数:", length(matched_genes), "\n")
cat("匹配比例:", round(length(matched_genes)/length(accessory_genes)*100, 1), "%\n")

#再把 GO Ontology 信息添加在上述 GO 富集结果中
accessory_GO <- accessory_go_rich@result
accessory_GO <- merge(accessory_GO, go2ont, by.x = 'ID', by.y = "go_id", all = FALSE)

#输出
write.table(accessory_GO, "/home/sunleting/storage/202602214-ag_function_enrich/output/go/sal_go_enrichment_results.txt", 
            sep = '\t', row.names = FALSE, quote = FALSE)


#例如:
#clusterProfiler 包里的一些默认作图方法，例如
# dotplot(accessory_go_rich,
#         showCategory=10,
#         color="p.adjust",
#         title = "accessory"
#         ) +
#   scale_color_gradient(high="#D53142",low="#345097") +
#   facet_grid(~Ontology)


bar_plot_all<-barplot(accessory_go_rich,
                      showCategory = 10,
                      color = "p.adjust",
                      title = "GO enrichment barplot") +
  scale_color_gradient(high="#D53142",low="#345097")
# 
output_dir<-"/home/sunleting/storage/202602214-ag_function_enrich/output/go"
ggsave(filename = file.path(output_dir, "sal_go_barplot_allterms.pdf"), 
       plot = bar_plot_all,  # 替换bar_plot为你的图形对象
       width = 8, 
       height = 6)

# cnet_plot<-cnetplot(accessory_go_rich) #网络图展示富集功能和基因的包含关系
# # emapplot(accessory_go_rich) #网络图展示各富集功能之间共有基因关系
# ggsave(filename=file.path(output_dir,"sal_go_cnetplot.pdf"),
#        plot=cnet_plot,
#        width = 20, 
#        height =9)

dot_plot_all<-dotplot(accessory_go_rich, showCategory = 20, title = "GO enrichment (all terms)") +
  theme_minimal()
ggsave(filename=file.path(output_dir,"sal_go_dotplot.pdf"),
       plot=dot_plot_all,
       width = 9, 
       height = 7)

dot_plot_20<-dotplot(accessory_go_rich, showCategory = sum(accessory_go_rich@result$p.adjust < 0.2)) +
  ggtitle("Significant GO Terms (p.adj < 0.2)")
# ggsave(filename=file.path(output_dir,"1-4pan_go_dotplot_20.pdf"), 
#        dotplot(accessory_go_rich, showCategory = sum(accessory_go_rich@result$p.adjust < 0.2)) +
#          ggtitle("Significant GO Terms (p.adj < 0.2)"),
#        width = 9, height = 7)


# 默认画前 20 个 term（按 p.adjust 排序）
barplot(accessory_go_rich, 
        showCategory = 15, 
        title = "GO Enrichment - Barplot",
        font.size = 12)
ggsave(filename=file.path(output_dir,"sal_go_barplot_top20.pdf"), 
       plot = last_plot(), 
       width = 8, height = 6)

