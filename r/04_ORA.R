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
                   minSize  = 5,
                   maxSize  = 500)

fgsea_res <- fgsea_res[fgsea_res$padj < 0.01, ]


# 选择前 20 个基因
top_genes <- important_features$feature[1:20]

# 定义函数进行超几何检验
ora_test <- function(gene_set, genes_of_interest, background_genes) {
  # gene_set: 一个基因集中的基因列表
  # genes_of_interest: 感兴趣的基因（如前20个基因）
  # background_genes: 背景基因集合（所有基因）
  
  overlap <- length(intersect(gene_set, genes_of_interest))
  set_size <- length(gene_set)
  total_genes <- length(background_genes)
  target_size <- length(genes_of_interest)
  
  p_value <- phyper(overlap - 1, set_size, total_genes - set_size, target_size, lower.tail = FALSE)
  return(p_value)
}

# 准备背景基因集合
background_genes <- rownames(DEGs)

# 对每个基因集进行ORA分析
ora_results <- data.frame(GeneSet = character(), PValue = numeric(), stringsAsFactors = FALSE)

for (gene_set_name in names(gmt_data)) {
  gene_set <- gmt_data[[gene_set_name]]
  p_value <- ora_test(gene_set, top_genes, background_genes)
  
  ora_results <- rbind(ora_results, data.frame(GeneSet = gene_set_name, PValue = p_value))
}

# 使用Benjamini-Hochberg方法校正p值
ora_results$AdjustedPValue <- p.adjust(ora_results$PValue, method = "BH")

# 添加其他参数列
ora_results$GeneSetSize <- sapply(gmt_data[ora_results$GeneSet], length)
ora_results$Overlap <- sapply(ora_results$GeneSet, function(gs) length(intersect(gmt_data[[gs]], top_genes)))
ora_results$EnrichmentFactor <- ora_results$Overlap / (ora_results$GeneSetSize * length(top_genes) / length(background_genes))

# 添加具体的overlap基因列，使用gene_name
ora_results$OverlapGenes <- sapply(ora_results$GeneSet, function(gs) {
  # 找到重叠的gene_id
  overlap_ids <- intersect(gmt_data[[gs]], top_genes)
  
  # 使用gene_info数据框将gene_id转换为gene_name
  overlap_names <- gene_info %>%
    filter(gene_id %in% overlap_ids) %>%
    pull(gene_name)
  
  # 将重叠基因名称用逗号连接成字符串
  paste(overlap_names, collapse = ", ")
})

# 根据校正后的p值排序
ora_results <- ora_results %>% arrange(AdjustedPValue)

# 查看显著结果
significant_results <- ora_results %>% filter(AdjustedPValue < 0.01)
print(significant_results)

# 修改Source列的筛选条件
significant_results$source <- ifelse(grepl("R-HSA-|Homo sapiens", significant_results$GeneSet), "Reactome",
                           ifelse(grepl("WP", significant_results$GeneSet), "WikiPathways",
                                  ifelse(grepl("_|outerRadialGlia", significant_results$GeneSet), "Literature",
                                         ifelse(grepl("UP|DN|\\.V[1-9]|CAHOY|PID PLK1 PATHWAY", significant_results$GeneSet) & !grepl("DNA", significant_results$GeneSet), "MsigDB",
                                                "KEGG"))))

significant_results <- significant_results %>%
  mutate(GeneSet = case_when(
    GeneSet == "Glial Cell Differentiation WP2276" ~ "Glial Cell Differentiation (WP2276)",
    TRUE ~ GeneSet
  ))



