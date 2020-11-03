options(repr.plot.width = 16, repr.plot.height = 10)
.libPaths()

library(dplyr)
library(Seurat)

# Provide directory of 10X data
gse_id <- "GSE132080"
base_dir <- file.path("data", "perturbseq")
data_dir <- file.path(base_dir, gse_id)

processed_output_file <- file.path(base_dir, paste0(gse_id, "_processed_matrix.tsv.gz"))
processed_output_file

# Initialize the Seurat object
perturbseq_data <- Seurat::Read10X(data.dir = data_dir)

perturbseq <- Seurat::CreateSeuratObject(
    counts = perturbseq_data,
    project = "crispri",
    min.cells = 3,
    min.features = 200
)

perturbseq

# Identify proportion of cells mapping to the mitochondrial genome
perturbseq[["percent.mt"]] <- Seurat::PercentageFeatureSet(perturbseq, pattern = "^MT-")

Seurat::VlnPlot(perturbseq, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Visualize relationships between counts
plot1 <- Seurat::FeatureScatter(perturbseq, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- Seurat::FeatureScatter(perturbseq, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter cells with too few and too much RNA, and cells that have outlier percentage of mito
perturbseq <- subset(
    perturbseq,
    subset = nFeature_RNA > 2000 & nFeature_RNA < 7000 & percent.mt < 15
)

# Normalize
perturbseq <- Seurat::NormalizeData(perturbseq, normalization.method = "LogNormalize", scale.factor = 10000)

# Find the top 2,000 most highly variable genes
perturbseq <- Seurat::FindVariableFeatures(perturbseq, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(perturbseq), 10)

top10

# Plot feature expression average and variance
plot1 <- VariableFeaturePlot(perturbseq)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Z score scale data and extract
all.genes <- rownames(perturbseq)
perturbseq <- ScaleData(perturbseq, features = all.genes)

scaled_perturbseq_data <- GetAssayData(object = perturbseq, slot = "scale.data")

dim(scaled_perturbseq_data)

# Subset to selected genes and output to file
tibble::rownames_to_column(
    as.data.frame(
        scaled_perturbseq_data[VariableFeatures(perturbseq), ]
    ),
    var = "gene"
) %>%
    readr::write_tsv(processed_output_file)
