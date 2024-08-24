rm(list = ls())
load(file = 'step1output.Rdata')

library(rtracklayer)
library(dplyr)
library(data.table)
library(fgsea)
library(ggplot2)

# 读取GTF文件
gtf_file <- "data/gencode.v46.chr_patch_hapl_scaff.annotation.gtf"
gtf_data <- import(gtf_file)

# 将GTF数据转换为数据框
gtf_df <- as.data.frame(gtf_data)

# 提取基因名和Ensemble ID
gene_info <- gtf_df %>%
  filter(type == "gene") %>%
  dplyr::select(gene_id = gene_id, gene_name = gene_name)

# 确保基因ID没有版本号
gene_info$gene_id <- sub("\\..*", "", gene_info$gene_id)

# 读取重要特征数据
important_features <- read.csv("data/feature_importance.csv")

# 为重要特征数据框添加symbol列
important_features <- important_features %>%
  left_join(gene_info, by = c("feature" = "gene_id"))

# 定义函数读取GMT文件
read_gmt <- function(file) {
  gmt_lines <- readLines(file)
  gmt_list <- list()
  
  for (line in gmt_lines) {
    parts <- unlist(strsplit(line, "\t")) # 按制表符分割每一行
    gene_set_name <- parts[1] # 基因集名称
    genes <- parts[3:length(parts)] # 跳过第二个空列，从第三个元素开始提取基因ID列表
    
    # 去除空字符串
    genes <- genes[genes != ""]
    
    gmt_list[[gene_set_name]] <- genes # 将基因集名称作为键，基因ID列表作为值
  }
  
  return(gmt_list)
}

# 读取GMT文件
gmt_file <- "data/combined_ensembl.gmt"
gmt_data <- read_gmt(gmt_file)

names(gmt_data)[names(gmt_data) == "PLK1 signaling events Homo sapiens e5e87977-6194-11e5-8ac5-06603eb7f303"] <- "PID PLK1 PATHWAY"

# 创建排名列表
res <- DEGs
ranks <- res$log2FoldChange
names(ranks) <- rownames(res)

# 按log2FoldChange排序
ranks <- sort(ranks, decreasing = TRUE)
set.seed(17)
fgsea_res <- fgsea(pathways = gmt_data, 
                  stats    = ranks,
                  eps      = 0.0,
                  minSize  = 15,
                  maxSize  = 500)

# 修改Source列的筛选条件
fgsea_res$source <- ifelse(grepl("R-HSA-|Homo sapiens", fgsea_res$pathway), "Reactome",
                              ifelse(grepl("WP", fgsea_res$pathway), "WikiPathways",
                                     ifelse(grepl("_|outerRadialGlia", fgsea_res$pathway), "Literature",
                                            ifelse(grepl("UP|DN|\\.V[1-9]|CAHOY|PID PLK1 PATHWAY", fgsea_res$pathway) & !grepl("DNA", fgsea_res$pathway), "MsigDB",
                                                   "KEGG"))))
fgsea_res <- fgsea_res[fgsea_res$padj < 0.01, ]

# 提取important_features的前20行，并选择feature列（Ensembl基因ID）
top20_features <- important_features %>%
  dplyr::slice(1:20) %>%
  dplyr::select(feature)

# 创建两个空列用于存储结果
fgsea_res$top20_genes_in_pathway <- NA
fgsea_res$num_top20_genes_in_pathway <- 0

# 对fgsea_res中的每一行进行操作
for (i in 1:nrow(fgsea_res)) {
  # 获取当前通路的基因集
  pathway_genes <- gmt_data[[fgsea_res$pathway[i]]]
  
  # 查找前20个基因中哪些在当前通路中
  matching_genes <- top20_features$feature[top20_features$feature %in% pathway_genes]
  
  # 存储匹配基因和数量
  fgsea_res$top20_genes_in_pathway[i] <- paste(matching_genes, collapse = ", ")
  fgsea_res$num_top20_genes_in_pathway[i] <- length(matching_genes)
}

# 补充gene_name信息到top20_genes_in_pathway列
fgsea_res$top20_genes_in_pathway <- sapply(fgsea_res$top20_genes_in_pathway, function(genes) {
  if (genes != "") {
    gene_ids <- unlist(strsplit(genes, ", "))
    gene_names <- gene_info$gene_name[match(gene_ids, gene_info$gene_id)]
    return(paste(gene_names, collapse = ", "))
  } else {
    return(NA)
  }
})

# 拆分数据框
split_data <- split(fgsea_res, fgsea_res$source)

# 输出每个拆分后的数据框
kegg_df <- split_data[["KEGG"]]
msigdb_df <- split_data[["MsigDB"]]
reactome_df <- split_data[["Reactome"]]
wikipathways_df <- split_data[["WikiPathways"]]
literature_df <- split_data[["Literature"]]

# 自定义 KEGG
library(stringr)
kegg_df$pathway <- str_to_title(kegg_df$pathway)

# 自定义 MsigDB
msigdb_df <- split_data[["MsigDB"]]
msigdb_df <- msigdb_df %>%
  filter(pathway %in% c("CAHOY NEURONAL", "CAHOY OLIGODENDROCUTIC", "PID PLK1 PATHWAY"))
# msigdb_df <- msigdb_df %>%
#   mutate(pathway = case_when(
#     pathway == "CAHOY NEURONAL" ~ "Genes up-regulated in neurons",
#     pathway == "CAHOY OLIGODENDROCUTIC" ~ "Genes up-regulated in oligodendrocytes",
#     pathway == "PID PLK1 PATHWAY" ~ "PLK1 signaling events",
#     TRUE ~ pathway
#   ))

# 自定义 Reactome
reactome_df <- split_data[["Reactome"]]
reactome_df$pathway <- sub("(.*)\\s[^\\s]*$", "\\1", reactome_df$pathway) %>% 
  sub("^(.)", "\\U\\1", ., perl = TRUE)

small_words <- c("a", "an", "and", "as", "at", "but", "by", "for", "in", "nor", "of", "on", "or", "so", "the", "to", "up", "with")

title_case_with_exceptions <- function(text, exceptions) {
  text <- str_to_title(text)
  
  for (word in exceptions) {
    text <- gsub(paste0("\\b", str_to_title(word), "\\b"), word, text, ignore.case = TRUE)
  }
  
  return(text)
}

reactome_df$pathway <- sapply(reactome_df$pathway, title_case_with_exceptions, exceptions = small_words)

reactome_df <- reactome_df %>%
  mutate(pathway = case_when(
    pathway == "Post Nmda Receptor Activation Events" ~ "Post NMDA Receptor Activation Events",
    TRUE ~ pathway
  ))


# 自定义 WikiPathways
wikipathways_df$pathway <- sub("(.*\\s)([^\\s]+)$", "\\1(\\2)", wikipathways_df$pathway) %>%
  sub("^(.)", "\\U\\1", ., perl = TRUE)

wikipathways_df <- wikipathways_df %>%
  mutate(pathway = case_when(
    pathway == "Neovascularisation processes (WP4331)" ~ "Neovascularisation Processes (WP4331)",
    TRUE ~ pathway
  ))


# 重新合并所有拆分后的数据框
new_fgsea_res <- rbind(kegg_df, reactome_df, wikipathways_df, msigdb_df)

# 创建气泡图
source_order <- c("KEGG", "Reactome", "MsigDB", "WikiPathways")
new_fgsea_res$source <- factor(new_fgsea_res$source, levels = source_order)
bubble_plot <- ggplot(new_fgsea_res, aes(x = NES, y = reorder(pathway, NES), size = size, color = padj)) +
  geom_point(alpha = 0.8, stroke = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey25", size = 0.5) + 
  scale_size_continuous(name = "No. genes", range = c(3, 8)) +
  scale_color_gradient(name = "p-adjusted", low = "#177F97", high = "#B7F1B2") +
  facet_grid(source ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),  # Remove left y-axis text
    axis.title.y = element_blank(), # Remove left y-axis title
    axis.text.y.right = element_text(size = 12, face = "bold"),
    axis.title.y.right = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    strip.text.y = element_text(size = 8, colour = "white", face = "bold"), 
    strip.background = element_rect(colour = "#1B1B1B", size = 1.5, fill = "#1B1B1B"),
    panel.border = element_rect(colour = "#1B1B1B", size = 1.5, fill = NA),
    panel.grid.major = element_line(size = 0.5, linetype = "dashed", colour = scales::alpha("grey70", 0.5)),
    panel.grid.minor = element_line(size = 0.25, linetype = "dashed", colour = scales::alpha("grey80", 0.3)),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway"
  ) +
  scale_y_discrete(position = "right") + # Move y-axis to the right
  guides(
    color = guide_colorbar(
      frame.linewidth = 0.5 # Set the border width
    )
  )

print(bubble_plot)

######## 火山图 ########
# 加载必要的库
library(ggplot2)
library(ggrepel)
library(dplyr)

# 设置阈值
logFC_t = 0.585
padj_t = 0.05

# 绘制基础火山图
p <- ggplot(data = DEGs, 
            aes(x = log2FoldChange, 
                y = -log10(padj))) +
  geom_point(alpha=0.6, size=4, 
             aes(color=Change), stroke=0.9) +
  scale_color_manual(values=c("#3C79C0", "grey","#CC1C2F")) +
  geom_vline(xintercept=c(-logFC_t, logFC_t), lty=4, col="black", linewidth=0.8) +
  geom_hline(yintercept = -log10(padj_t), lty=4, col="black", linewidth=0.8) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_blank(),  # 移除图例标题
        legend.text = element_text(size = 14, face = "bold")
  )

# 将行名转换为列
DEGs_id <- DEGs %>%
  rownames_to_column(var = "ensemble_id")

# 将DEGs数据框与gene_info合并，添加symbol列
DEGs_id <- DEGs_id %>%
  left_join(gene_info, by = c("ensemble_id" = "gene_id"))

# 筛选重要基因
for_label <- DEGs_id %>% 
  filter(gene_name %in% important_features$gene_name[1:20])


# 添加标签并保存火山图
volcano_plot <- p +
  geom_point(size = 4, shape = 1, data = for_label, stroke=0.9) +
  ggrepel::geom_label_repel(
    aes(label = gene_name),
    data = for_label,
    fontface = "bold",
    color = "black",
    size = 4,  # 设置label的字体大小，例如3
    box.padding = unit(0.35, "lines"),
    label.padding = unit(0.3, "lines"),
    label.size = 0.5,  # 设置标签框的粗细
    fill = alpha("white", 0.8)  # 设置标签框的填充颜色和透明度
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    panel.border = element_rect(colour = "#1B1B1B", size = 1.5, fill = NA),
    panel.grid.major = element_line(size = 0.5, linetype = "dashed", colour = scales::alpha("grey70", 0.5)),
    panel.grid.minor = element_line(size = 0.25, linetype = "dashed", colour = scales::alpha("grey80", 0.3))
  )

# 显示火山图
print(volcano_plot)

ggsave("volcano.png", plot = volcano_plot, width = 9, height = 6, dpi = 600)

library(patchwork)

# 假设 volcano_plot 和 bubble_plot 是你要合并的两个图
combined_plot <- volcano_plot + bubble_plot + 
  plot_layout(ncol = 2, widths = c(2.5, 1))  # 设置比例，2:1 表示火山图宽度是气泡图的两倍

# 显示合并后的图
print(combined_plot)

# 保存图像
ggsave("combined_plot.png", plot = combined_plot, width = 15, height = 6, dpi = 600)

ggsave("bubble_plot.png", plot = bubble_plot, width = 9.5, height = 7.5, dpi = 600)

# 保留指定的列
df <- new_fgsea_res[, c("pathway", "padj", "NES", "size", "source")]
write.csv(df, file = "df.csv", row.names = FALSE)


# Assuming the list column is named 'list_col'
literature_df$top20_genes_in_pathway <- sapply(literature_df$top20_genes_in_pathway, toString)
literature_df$leadingEdge <- sapply(literature_df$leadingEdge, toString)
filtered_df <- literature_df[literature_df$num_top20_genes_in_pathway > 2, ]

# Write the filtered data frame to a CSV file
write.csv(filtered_df, file = "literature_df_filtered.csv", row.names = FALSE)
