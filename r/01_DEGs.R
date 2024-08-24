rm(list = ls())
load(file = 'step0output.Rdata')

library(DESeq2)
library(dplyr)

######## 函数：差异表达分析 ########
DGE_analysis <- function(colData, countData, padj_t = 0.05, logFC_t = 0.585) {
  
  ## 过滤低表达数据
  keep_genes <- sapply(unique(colData$Stage), function(stage) {
    stage_samples <- colData$ID[colData$Stage == stage]
    stage_counts <- countData[, stage_samples]
    rowSums(stage_counts > 0) >= floor(0.50 * length(stage_samples))
  })
  countData <- countData[rowSums(keep_genes) > 0, ]
  print(paste("After Filtering Genes:", 
              paste(dim(countData)[1], "Genes", dim(countData)[2], "Patients")))
  
  ## DESeq2分析
  colData <- data.frame(row.names = colnames(countData), 
                        condition = factor(colData$Stage, levels = c("Primary", "Recurrent")), 
                        subject = factor(gsub("_P|_R", "", colData$ID)))
  suppressMessages({
    dds <- DESeqDataSetFromMatrix(countData = countData, 
                                  colData = colData, 
                                  design = ~ subject + condition)
    dds <- DESeq(dds)
  })
  
  ## 提取结果并按padj排序
  res <- results(dds, contrast = c("condition", "Recurrent", "Primary"))
  res <- as.data.frame(res)
  res <- arrange(res, padj)
  res <- na.omit(res)
  
  ## 根据padj和log2FoldChange筛选显著上下调基因
  down <- (res$padj < padj_t) & (res$log2FoldChange < -logFC_t)
  up <- (res$padj < padj_t) & (res$log2FoldChange > logFC_t)
  res$Change <- ifelse(down, "Down", ifelse(up, "Up", "Stable"))
  
  ## 输出显著上下调和稳定基因的数量
  change_table <- table(res$Change)
  print(paste("Number of DEGs: Upregulated =", change_table["Up"], 
              "Downregulated =", change_table["Down"], 
              "Stable =", change_table["Stable"])) 
  
  return(res)
}

## 执行差异表达分析
colData <- rbind(colData_train, colData_test)
countData <- cbind(countData_train, countData_test)

DEGs <- DGE_analysis(colData_test, countData_test)


## 保存数据
save(DEGs, file = "step1output.Rdata")
