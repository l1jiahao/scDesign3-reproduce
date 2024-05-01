# library(scDesign3)
# library(SingleCellExperiment)
# library(ggplot2)
# theme_set(theme_bw())

# ### Read csv data
# # Load gene expression data
# counts <- read.csv("./csv_format/muris_gene_expression_countT.csv", row.names = 1)
# # csv 读取的时候会把 colnames 中的 '-' 读取成 '.', 现在再次转换一下
# colnames(counts) <- gsub("\\.", "-", colnames(counts))
# colnames(counts) <- gsub("^X([0-9-]+)$", "\\1", colnames(counts))

# # Load metadata
# metadata <- read.csv("./csv_format/muris_obs.csv", row.names = 1)

# counts <- counts[, rownames(metadata)]
# # Ensure column names in counts match row names in metadata
# if(!all(colnames(counts) == rownames(metadata))) stop("Misalignment between counts and metadata")

# # Check for NA values
# if(any(is.na(counts))) stop("NA values in counts")
# if(any(is.na(metadata))) stop("NA values in metadata")

# # Confirm matching of columns and rows
# # Create SingleCellExperiment object
# sce <- SingleCellExperiment(
#   assays = list(counts=counts),
#   colData = metadata)

# # Convert celltype to factor ensuring no NA values exist
# # 检查 celltype 列是否含有 NA 值
# if(any(is.na(metadata$celltype))) {stop("NA values in celltype column of metadata")}
# print(dim(sce))
# print(summary(colData(sce)))
# # print(table(metadata$celltype))
# colData(sce)$celltype <- as.factor(colData(sce)$celltype)
# # logcounts(sce) <- counts(sce)

library(scDesign3)
library(SingleCellExperiment)
library(data.table)  # 加载 data.table 包
library(Matrix)      # 加载 Matrix 包处理稀疏矩阵
library(ggplot2)
theme_set(theme_bw())

# 使用 fread 读取 CSV 数据，并立即进行必要的列选择和类型转换
counts <- fread("./csv_format/muris_gene_expression_countT.csv", data.table = FALSE)
rownames(counts) <- counts[[1]]  # 假设第一列为行名
counts <- counts[,-1]            # 移除第一列

# 转换列名中的 '.' 为 '-'
colnames(counts) <- gsub("\\.", "-", colnames(counts))
colnames(counts) <- gsub("^X([0-9-]+)$", "\\1", colnames(counts))

metadata <- fread("./csv_format/muris_obs.csv", data.table = FALSE)
rownames(metadata) <- metadata[[1]]
metadata <- metadata[,-1]

# 使用 Matrix 包将 counts 转换为稀疏矩阵，减少内存使用
sparse_counts <- as(counts, "sparseMatrix")

counts <- NULL  # 不再需要 counts 数据框，释放其内存
gc()            # 调用垃圾回收

# 创建 SingleCellExperiment 对象
sce <- SingleCellExperiment(
  assays = list(counts = sparse_counts),
  colData = metadata
)


set.seed(123)
example_simu <- scdesign3(
    sce = sce,
    assay_use = "counts",
    celltype = "celltype",
    pseudotime = NULL,
    spatial = NULL,
    other_covariates = NULL,
    family_use = "nb",
    n_cores = 1,
    usebam = FALSE,
    mu_formula = "celltype",
    sigma_formula = "celltype",
    corr_formula = "celltype",
    copula = "vine",
    DT = TRUE,
    pseudo_obs = FALSE,
    return_model = FALSE,
    nonzerovar = TRUE,
    nonnegative = TRUE
)

saveRDS(example_simu, file = "./output/muris_gen.rds")
