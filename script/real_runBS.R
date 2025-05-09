#!/usr/bin/env Rscript

read.lister <- function(file) {
  dat <- read.table(
    file,
    skip = 1,
    row.names = NULL,
    col.names = c("chr", "pos", "strand", "context", "M", "Cov"),
    colClasses = c("character", "integer", "character", "character",
                   "integer", "integer"))
  
  dat <- dat[dat$context == "CG", ]
  dat$context <- NULL
  dat$chr <- paste("chr", dat$chr, sep = "")
  tmp <- dat[dat$strand == "+",]
  BS.forward <- BSseq(
    pos = tmp$pos,
    chr = tmp$chr,
    M = as.matrix(tmp$M, ncol = 1),
    Cov = as.matrix(tmp$Cov, ncol = 1),
    sampleNames = "forward")
  tmp <- dat[dat$strand == "-",]
  BS.reverse <- BSseq(
    pos = tmp$pos - 1L,
    chr = tmp$chr,
    M = as.matrix(tmp$M, ncol = 1),
    Cov = as.matrix(tmp$Cov, ncol = 1),
    sampleNames = "reverse")
  BS <- combine(BS.forward, BS.reverse)
  BS <- collapseBSseq(BS, group = c("a", "a"), type = "integer")
  BS
}

# ========== Step 0: 读取命令行参数 ==========
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Please provide sample parameters, such as: Rscript bsmooth_dmrs.R /mnt/data/samples_20")
}
#sample_tag <- args[1]
input_dir <- args[1]
g1 <- args[2]
g2 <- args[3]
threads <- args[4]

input_dir <- "/home/user/liyang/project/methDmr/real-data/prostate_cancer/GSE158927/split/result/formatted_BSmooth"
g1 <- "lethal"
g2 <- "normal"
threads <- 8

# ========== Step 1: 加载包 ==========
suppressMessages(library(bsseq))

# ========== Step 2: 构建路径 ==========
sample_tag <- basename(input_dir)
output_dir <- file.path(input_dir, "BSmooth_result")
#print(output_dir)  
meta_file <- file.path(output_dir, "in_bs.list")

if (!file.exists(meta_file)) {
  stop("Unable to find sample list file：", meta_file)
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ========== Step 3: 读取样本信息 ==========
meta <- read.csv(meta_file, header = FALSE, stringsAsFactors = FALSE, sep = '\t')
colnames(meta) <- c("Sample", "Group", "Directory")

# ========== Step 4: 构建 BSseq List ==========
bs_list <- list()
for (i in seq_len(nrow(meta))) {
  bs <- read.lister(meta$Directory[i])  # 你自己定义的 read.lister 函数应提前存在或source()
  name <- paste(meta$Group[i], meta$Sample[i], sep = "_")
  sampleNames(bs) <- name
  bs_list[[name]] <- bs
}

# ========== Step 5: 合并样本 ==========
bs_all <- Reduce(combine, bs_list)
group_vector <- c(rep(g1, sum(meta$Group == g1)),
                  rep(g2, sum(meta$Group == g2)))
pData(bs_all)$Rep <- group_vector

# ========== Step 6: BSmooth ==========
bs_smooth <- BSmooth(bs_all)

# ========== Step 7: t统计量 ==========
group1 <- sampleNames(bs_smooth)[pData(bs_smooth)$Rep == g1]
group2 <- sampleNames(bs_smooth)[pData(bs_smooth)$Rep == g2]

bs_tstat <- BSmooth.tstat(
  BSseq = bs_smooth, 
  group1 = group1, 
  group2 = group2,
  estimate.var = "group2",
  local.correct = FALSE,
  mc.cores = threads,
  verbose = TRUE
)

# ========== Step 8: 检测 DMR ==========
dmrs <- dmrFinder(bs_tstat, cutoff = c(-4.6, 4.6), stat = "tstat")
dmrs <- dmrFinder(bs_tstat, qcutoff = c(0.05, 0.95), stat = "tstat" )
dmrs_filter <- subset(dmrs, n >= 5 & abs(meanDiff) >= 0.1)

# ========== Step 9: 保存结果 ==========
write.table(dmrs, file.path(output_dir, "dmrs.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(dmrs_filter, file.path(output_dir, "dmrs_filter.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("R detection completed, results have been saved to: ", output_dir, "\n")


