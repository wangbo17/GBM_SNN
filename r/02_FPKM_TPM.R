rm(list = ls())
load(file = 'step0output.Rdata')

library(edgeR)

######## 函数：FPKM 标准化 ########
calculate_log_fpkm <- function(colData, countData, gene_lengths) {
  dge <- DGEList(counts = countData, group = factor(colData$Stage, levels = c("Primary", "Recurrent")))
  cpm_values <- cpm(dge, normalized.lib.sizes = FALSE)
  gene_lengths <- gene_lengths[match(rownames(countData), gene_lengths$ID), ]
  gene_lengths_kb <- gene_lengths$Length / 1000
  fpkm_values <- sweep(cpm_values, 1, gene_lengths_kb, FUN = "/")
  log_fpkm_values <- log2(fpkm_values + 1)
  return(t(log_fpkm_values))
}

log_fpkm_test <- calculate_log_fpkm(colData_test, countData_test, gene_lengths)
log_fpkm_train <- calculate_log_fpkm(colData_train, countData_train, gene_lengths)

## 保存数据
write.csv(log_fpkm_test, 'data/fpkm_test.csv')
write.csv(log_fpkm_train, 'data/fpkm_train.csv')


######## 函数：TPM 标准化 ########
calculate_log_tpm <- function(colData, countData, gene_lengths) {
  gene_lengths <- gene_lengths[match(rownames(countData), gene_lengths$ID), ]
  gene_lengths_kb <- gene_lengths$Length / 1000
  rpk_values <- sweep(countData, 1, gene_lengths_kb, FUN = "/")
  scaling_factors <- colSums(rpk_values)
  tpm_values <- sweep(rpk_values, 2, scaling_factors, FUN = "/") * 1e6
  log_tpm_values <- log2(tpm_values + 1)
  return(t(log_tpm_values))
}

log_tpm_test <- calculate_log_tpm(colData_test, countData_test, gene_lengths)
log_tpm_train <- calculate_log_tpm(colData_train, countData_train, gene_lengths)

## 保存数据
write.csv(log_tpm_test, 'data/tpm_test.csv')
write.csv(log_tpm_train, 'data/tpm_train.csv')

##未转置
calculate_log_tpm <- function(colData, countData, gene_lengths) {
  gene_lengths <- gene_lengths[match(rownames(countData), gene_lengths$ID), ]
  gene_lengths_kb <- gene_lengths$Length / 1000
  rpk_values <- sweep(countData, 1, gene_lengths_kb, FUN = "/")
  scaling_factors <- colSums(rpk_values)
  tpm_values <- sweep(rpk_values, 2, scaling_factors, FUN = "/") * 1e6
  log_tpm_values <- log2(tpm_values + 1)
  return(log_tpm_values)
}

calculate_log_fpkm <- function(colData, countData, gene_lengths) {
  dge <- DGEList(counts = countData, group = factor(colData$Stage, levels = c("Primary", "Recurrent")))
  cpm_values <- cpm(dge, normalized.lib.sizes = FALSE)
  gene_lengths <- gene_lengths[match(rownames(countData), gene_lengths$ID), ]
  gene_lengths_kb <- gene_lengths$Length / 1000
  fpkm_values <- sweep(cpm_values, 1, gene_lengths_kb, FUN = "/")
  log_fpkm_values <- log2(fpkm_values + 1)
  return(log_fpkm_values)
}

colData <- rbind(colData_test, colData_train)
countData <- cbind(countData_test, countData_train)
log_tpm_countData <- calculate_log_tpm(colData, countData, gene_lengths)
log_fpkm_countData <- calculate_log_fpkm(colData, countData, gene_lengths)

## 保存数据
save(colData, countData, log_tpm_countData, log_fpkm_countData,
     file = "step2output.Rdata")
