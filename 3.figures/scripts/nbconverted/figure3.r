suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

source(file.path("utils", "themes.R"))

# Load perturbseq results
perturbseq_results_dir <- file.path("../1.calculate-metrics/perturb-seq/results")

gse_id <- "GSE132080"
results_file <- file.path(perturbseq_results_dir, paste0(gse_id, "_grit.tsv"))

output_dir <- "figures"

# Load bulk data
bulk_grit_cols <- readr::cols(
    perturbation = readr::col_character(),
    group = readr::col_character(),
    grit = readr::col_double(),
    id = readr::col_character(),
    sequence = readr::col_character(),
    gene = readr::col_character(),
    gamma_day5 = readr::col_double(),
    gamma_day10 = readr::col_double(),
    relative_activity_day5 = readr::col_double(),
    relative_activity_day10 = readr::col_double()
)

bulk_grit_df <- readr::read_tsv(results_file, col_types = bulk_grit_cols)

# Order for plotting
bulk_grit_df$gene <- factor(bulk_grit_df$gene, levels = unique(bulk_grit_df$gene))

print(dim(bulk_grit_df))
head(bulk_grit_df)

# Load single cell grit results
sc_results_file <- file.path(perturbseq_results_dir, paste0(gse_id, "_single_cell_grit.tsv.gz"))

sc_grit_cols <- readr::cols(
    perturbation = readr::col_character(),
    group = readr::col_character(),
    grit = readr::col_double(),
    grit_gene = readr::col_character(),
    grit_guide = readr::col_character()
)

sc_df <- readr::read_tsv(sc_results_file, col_types = sc_grit_cols)

# Load UMAP embeddings
sc_embeddings_file <- file.path(perturbseq_results_dir, paste0(gse_id, "_single_cell_umap_embeddings.tsv.gz"))

sc_embeddings_cols <- readr::cols(
    Metadata_cell_identity = readr::col_character(),
    Metadata_cell_barcode = readr::col_character(),
    Metadata_guide_identity = readr::col_character(),
    Metadata_read_count = readr::col_double(),
    Metadata_UMI_count = readr::col_double(),
    Metadata_coverage = readr::col_double(),
    Metadata_gemgroup = readr::col_double(),
    Metadata_good_coverage = readr::col_logical(),
    Metadata_number_of_cells = readr::col_double(),
    Metadata_gene_identity = readr::col_character(),
    Metadata_barcode = readr::col_character(),
    Metadata_sequence = readr::col_character(),
    umap_0 = readr::col_double(),
    umap_1 = readr::col_double(),
    grit_gene = readr::col_character()
)

sc_embeddings_df <- readr::read_tsv(sc_embeddings_file, col_types = sc_embeddings_cols)

# Merge single cell data
sc_df <- sc_embeddings_df %>%
    dplyr::right_join(
        sc_df,
        by = c("Metadata_cell_identity" = "perturbation", "grit_gene" = "grit_gene")
    ) %>%
    dplyr::full_join(
        bulk_grit_df,
        by = c("Metadata_guide_identity" = "perturbation", "Metadata_gene_identity" = "gene"),
        suffix = c("", "_bulk_activity"),
        keep = TRUE
)

sc_df$gene = factor(sc_df$gene, levels = unique(bulk_grit_df$gene))

sc_df$Metadata_gene_identity = factor(
    sc_df$Metadata_gene_identity, levels = c("neg", paste(unique(bulk_grit_df$gene)))
)

print(dim(sc_df))
head(sc_df)

table(sc_df$gene)

table(sc_df$Metadata_gene_identity)

# Load single cell gene expression data
perturbseq_data_dir <- file.path("../0.download-data/data/perturbseq/")

gene_exp_file <- file.path(perturbseq_data_dir, paste0(gse_id, "_final_analytical.tsv.gz"))

sc_gene_cols <- readr::cols(
    .default = readr::col_double(),
    Metadata_cell_identity = readr::col_character(),
    Metadata_cell_barcode = readr::col_character(),
    Metadata_guide_identity = readr::col_character(),
    Metadata_good_coverage = readr::col_logical(),
    Metadata_gene_identity = readr::col_character(),
    Metadata_barcode = readr::col_character(),
    Metadata_sequence = readr::col_character()
)

sc_gene_exp_df <- readr::read_tsv(gene_exp_file, col_types = sc_gene_cols)

print(dim(sc_gene_exp_df))
head(sc_gene_exp_df)

bulk_gene_gg = (
    ggplot(bulk_grit_df, aes(y = relative_activity_day5, x = grit))
    + geom_point(size = 1.2, alpha = 0.8)
    + xlab("Grit (bulk)")
    + ylab("Relative Activity (Day 5)")
    + facet_wrap("~gene", scales = "free_x")
    + custom_grit_theme
    + theme(axis.text = element_text(size = 6))
    + labs(tag = "a")
)

bulk_gene_gg

top_genes <- c("HSPA5", "GATA1", "RPL9", "GINS1", "HSPA9")

top_gene_embeddings <- list()
gene_ggs <- list()

for (gene in top_genes) {
    gene_embedding_df <- sc_df %>%
        dplyr::filter(grit_gene == !!gene) %>%
        dplyr::mutate(grit_facet_label = paste("RA:", round(relative_activity_day5, 3)))

    gene_neg_ctrl <- "Neg. Ctrl"
    gene_embedding_df[gene_embedding_df$Metadata_gene_identity == "neg", "grit_facet_label"] <- gene_neg_ctrl

    factor_order <- c(gene_neg_ctrl, paste("RA:", unique(round(sort(gene_embedding_df$relative_activity_day5), 3))))

    gene_embedding_df$grit_facet_label <- factor(
        gene_embedding_df$grit_facet_label, levels = factor_order
    )
    
    top_gene_embeddings[[gene]] <- gene_embedding_df
    
    gene_gg = (
        ggplot(gene_embedding_df, aes(x = umap_0, y = umap_1))
        + geom_point(aes(fill = grit), size = 1.8, pch = 21, stroke = 0, alpha = 0.5)
        + facet_grid("grit_facet_label~grit_gene", scales = "free", as.table = TRUE)
        + custom_grit_theme
        + xlab("")
        + ylab("")
        + theme(
            plot.margin = margin(t = 0.25, r = 0, b = 0.25, l = -0.7, "cm"),
            axis.text = element_text(size = 6),
        )
    )
    
    if (gene == "HSPA5") {
        gene_gg <- gene_gg + ylab("UMAP 1") + labs(tag = "b")
    }
    
    gene_ggs[[gene]] <- gene_gg
}

top_gene_embeddings <- dplyr::bind_rows(top_gene_embeddings)

print(dim(top_gene_embeddings))
head(top_gene_embeddings)

grit_range <- range(top_gene_embeddings$grit, na.rm = TRUE)

sc_grit_gg <- wrap_plots(gene_ggs, ncol = length(gene_ggs)) +
    plot_layout(guides = "collect") & 
    scale_colour_continuous(limits = grit_range) &
    xlab("UMAP 0") &
    scale_fill_gradientn(
        name = "scGrit",
        colors = scgrit_palette(10),
        values = scales::rescale(c(max(grit_range) + 1, 3, 0)),
        breaks = c(0, 1, 3, 5, 7),
        limits = grit_range
    )

sc_grit_gg

patchwork_plot <- (bulk_gene_gg | sc_grit_gg) + plot_layout(widths = c(0.6, 1))

for (extension in c(".png", ".pdf")) {
    output_file <- file.path(output_dir, paste0("figure3", extension))
    ggsave(output_file, patchwork_plot, width = 12, height = 6.2)
}

patchwork_plot
