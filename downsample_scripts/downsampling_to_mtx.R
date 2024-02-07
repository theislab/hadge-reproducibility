library(Seurat)
library(DropletUtils)
library(argparse)
library(dplyr)


parser <- ArgumentParser("Parameters for 10x matrix")
parser$add_argument("--hto_h5", help = "Path to HTO H5 ")
parser$add_argument("--hto_mtx_path", help = "Path to save HTO mtx. ")

args <- parser$parse_args()

hto_h5 <- Read10X_h5(args$hto_h5)

write10xCounts(args$hto_mtx_path,hto_h5, gene.type = "Gene Expression",version="3",)

print("----------------------------------")
print(args$hto_mtx_path)
print("10x counts folder created")