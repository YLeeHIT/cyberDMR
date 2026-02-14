library(DSS)
library(bsseq)

# ===============================
# 1. Specify the input directory
# ===============================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript run_dss_auto.R <input_dir>")
  }

input_dir <- args[1]

# ===============================
# 2. Automatically search for sample files
# ===============================

control_files <- list.files(
  input_dir,
  pattern = "^noh_sorted_control_.*\\.tsv$",
  full.names = TRUE
  )

treat_files <- list.files(
  input_dir,
  pattern = "^noh_sorted_treatment_.*\\.tsv$",
  full.names = TRUE
  )

if (length(control_files) == 0 || length(treat_files) == 0) {
  stop("No control or treatment files found.")
  }  

cat("Control samples:\n")
print(control_files)

cat("Treatment samples:\n")
print(treat_files)

all_files <- c(control_files, treat_files)
sample_names <- basename(all_files)

# ===============================
# 3. Read data
# ===============================

read_one <- function(f) {
  df <- read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
    col.names = c("chr", "pos", "N", "X")
    )
  df
  }

bs_list <- lapply(all_files, read_one)

BSobj <- makeBSseqData(bs_list, sampleNames = sample_names)

# ===============================
# 4. DML inspection
# ===============================

group <- c(
  rep("Control", length(control_files)),
  rep("Treat", length(treat_files))
  )

dmlTest <- DMLtest(
  BSobj,
  group1 = which(group == "Control"),
  group2 = which(group == "Treat"),
  smoothing = TRUE,
  smoothing.span = 500
  )

# ===============================
# 5. Call DMR (default parameter)
# ===============================

dmrs <- callDMR(dmlTest)

# ===============================
# 6. output result
# ===============================

out_file <- file.path(input_dir, "DSS_DMR_result.tsv")

write.table(
  dmrs,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
  )

cat("DMR calling finished.\n")
cat("Output:", out_file, "\n")

  
