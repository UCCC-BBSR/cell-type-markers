# make a reference matrix from a Seurat object
# This script reads a Seurat object from an RDS file, generates a reference matrix using the `make_ref_matrix` function from the `clustifyr` package,
# and then saves the resulting reference matrix to a specified output RDS file.
# The script uses the `optparse` package to handle command-line arguments for input and output file paths, as well as optional parameters for label column, assay, and slot.

# example usage:
# Rscript make-refmatrix-from-seurat.R -i seurat.rds -o refmatrix.rds -a RNA -l cell_type -s data

args <- commandArgs(trailingOnly = TRUE)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Input Seurat RDS file", metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", help = "Output reference matrix RDS file", metavar = "FILE"),
  make_option(c("-l", "--label_col"), type = "character", default = "cell_type", help = "Metadata column with labels [default %default]"),
  make_option(c("-a", "--assay"), type = "character", default = "RNA", help = "Assay to use [default %default]"),
  make_option(c("-s", "--slot"), type = "character", default = "data", help = "Slot to use [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

library("clustifyr")
library("Seurat")


so <- readRDS(opt$input)
ref_mat <- make_ref_matrix(
  so,
  label_col = opt$label_col,
  assay = opt$assay,
  slot = opt$slot
)
saveRDS(ref_mat, file = opt$output)