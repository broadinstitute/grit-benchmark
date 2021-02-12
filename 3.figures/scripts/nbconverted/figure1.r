suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(magick))

source(file.path("utils", "themes.R"))

# Set file names
output_dir <- "figures"

cell_health_dir <- file.path("../1.calculate-metrics/cell-health/results")
ceres_dir <- file.path("../2.compare-metrics/cell-health/results")

grit_file <- file.path(cell_health_dir, "cell_health_grit.tsv")
non_replicate_cor_file <- file.path(cell_health_dir, "cell_health_nonreplicate_95thpercentile.tsv")
reprod_file <- file.path(cell_health_dir, "cell_health_replicate_reproducibility.tsv")
ceres_file <- file.path(ceres_dir, "cell_health_grit_ceres.tsv")

# Model predictions
commit <- "7d63d4a43014c757fd0c77c0fd1c19540f17cc3d"
count_file <- paste0("https://raw.githubusercontent.com/broadinstitute/cell-health/", commit, "/1.generate-profiles/tables/cell_count_summary.tsv")

sc_plate <- "SQ00014613"
sc_grit_file <- file.path(cell_health_dir, paste0("cellhealth_single_cell_grit_", sc_plate, "_chr2.tsv.gz"))
sc_umap_file <- file.path(cell_health_dir, paste0("cellhealth_single_cell_umap_embeddings_", sc_plate, "_chr2.tsv.gz"))

# Load data
grit_cols <- readr::cols(
    perturbation = readr::col_character(),
    group = readr::col_character(),
    grit = readr::col_double(),
    cell_line = readr::col_character(),
    barcode_control = readr::col_character(),
    cor_method = readr::col_character()
)

non_rep_cols <- readr::cols(
    cell_line = readr::col_character(),
    similarity_metric = readr::col_double()
)

rep_cols <- readr::cols(
    cell_line = readr::col_character(),
    group = readr::col_character(),
    perturbation = readr::col_character(),
    median_replicate_correlation = readr::col_double(),
    median_control_correlation = readr::col_double()
)

count_cols <- readr::cols(
    gene_name = readr::col_character(),
    pert_name = readr::col_character(),
    cell_line = readr::col_character(),
    cell_count = readr::col_double()
)

ceres_cols <- readr::cols(
    perturbation = readr::col_character(),
    group = readr::col_character(),
    grit = readr::col_double(),
    cell_line = readr::col_character(),
    barcode_control = readr::col_character(),
    cor_method = readr::col_character(),
    query = readr::col_character(),
    `_id` = readr::col_double(),
    DepMap_ID = readr::col_character(),
    stripped_cell_line_name = readr::col_character(),
    ceres_score = readr::col_double(),
    grit_mean = readr::col_double()
)

grit_df <- readr::read_tsv(grit_file, col_types = grit_cols)
nonrep_df <- readr::read_tsv(non_replicate_cor_file, col_types = non_rep_cols)
rep_df <- readr::read_tsv(reprod_file, col_types = rep_cols)
count_df <- readr::read_tsv(count_file, col_types = count_cols)
ceres_df <- readr::read_tsv(ceres_file, col_types = ceres_cols)

# Load single cell data
sc_grit_cols <- readr::cols(
    perturbation = readr::col_character(),
    group = readr::col_character(),
    grit = readr::col_double(),
    gene = readr::col_character(),
    guide = readr::col_character(),
    grit_gene = readr::col_character(),
    grit_guide = readr::col_character()
)

sc_umap_cols <- readr::cols(
    Metadata_cell_identity = readr::col_character(),
    Metadata_Plate = readr::col_character(),
    Metadata_Well = readr::col_character(),
    Metadata_gene_name = readr::col_character(),
    Metadata_pert_name = readr::col_character(),
    Metadata_broad_sample = readr::col_character(),
    Metadata_cell_line = readr::col_character(),
    Metadata_TableNumber = readr::col_character(),
    Metadata_ImageNumber = readr::col_double(),
    Metadata_ObjectNumber_cytoplasm = readr::col_double(),
    Metadata_Cytoplasm_Parent_Cells = readr::col_double(),
    Metadata_Cytoplasm_Parent_Nuclei = readr::col_double(),
    Metadata_ObjectNumber_cells = readr::col_double(),
    Metadata_ObjectNumber = readr::col_double(),
    umap_0 = readr::col_double(),
    umap_1 = readr::col_double(),
    grit_gene = readr::col_character()
)

sc_grit_df <- readr::read_tsv(sc_grit_file, col_types = sc_grit_cols)
sc_umap_df <- readr::read_tsv(sc_umap_file, col_types = sc_umap_cols)

url <- "https://raw.githubusercontent.com/broadinstitute/grit-benchmark/2a4b0edc3d96b67b59d5071a4c3abb201725ac37/media/Figure1A_Simple_Interpretation.png"
panel_a_gg <- magick::image_ggplot(
    magick::image_scale(
        magick::image_read(url),
        "3500x3500"
    ), interpolate = TRUE
) + labs(tag = "a") + theme(plot.margin = unit(c(-10, 0, -10, 0), "cm"))

panel_a_gg

grit_focus_df <- grit_df %>%
    dplyr::filter(barcode_control == "cutting_control", cor_method == "pearson") %>%
    dplyr::left_join(rep_df, by = c("perturbation", "group", "cell_line")) %>%
    tidyr::drop_na()

head(grit_focus_df)

mean_grit_df <- grit_focus_df %>%
    dplyr::group_by(cell_line) %>%
    dplyr::mutate(mean_grit = round(mean(grit), 3)) %>%
    dplyr::select(cell_line, mean_grit) %>%
    unique() %>%
    dplyr::ungroup() %>%
    dplyr::arrange(cell_line)

colnames(mean_grit_df) <- c("Cell line", "grit (mean)")

mean_grit_df

table_grob_theme <- gridExtra::ttheme_default(
    core = list(
        fg_params = list(cex = 0.6),
        padding = grid::unit.c(unit(1, "mm"), unit(1, "mm"))
    ),
    colhead = list(
        fg_params = list(cex = 0.7),
        padding = grid::unit.c(unit(3, "mm"), unit(2, "mm"))
        )
)

grit_mean_grob <- tableGrob(mean_grit_df, theme = table_grob_theme, rows = NULL)

mean_grit_gg <- (
    ggplot(grit_focus_df, aes(x = grit, fill = cell_line))
    + geom_density(alpha = 0.5)
    + scale_fill_manual(name = "Cell line", values = cell_line_colors)
    + annotation_custom(
        grit_mean_grob,
        xmin = 2,
        xmax = 3,
        ymin = 0.55,
        ymax = 0.75
    )
    + xlab("Grit")
    + ylab("Density")
    + custom_grit_theme
    + labs(tag = "b")
)

mean_grit_gg

replicate_reprod_panel <- (
    ggplot(grit_focus_df, aes(x = median_control_correlation, y = median_replicate_correlation))
    + geom_point(
        aes(size = grit, fill = grit, shape = cell_line),
        alpha = 0.7
    )
    + geom_point(
        data = grit_focus_df %>% dplyr::filter(grit < 1),
        aes(size = grit, shape = cell_line),
        fill="grey",
        alpha = 0.8
    )
    + scale_fill_gradient(name = "Grit", high = "yellow", low = "red")
    + scale_shape_manual(name = "Cell line", values = c(21, 23, 25))
    + scale_size_continuous(guide = FALSE, range = c(0.2, 2.5))
    + geom_vline(xintercept = 0, linetype = "dashed", color = "red")
    + geom_hline(yintercept = 0, linetype = "dashed", color = "red")
    + geom_hline(yintercept = mean(nonrep_df$similarity_metric), linetype = "dashed", color = "blue")
    + ylab("Median replicate correlation")
    + xlab("Median control correlation")
    + custom_grit_theme
    + labs(tag = "c")
)

replicate_reprod_panel

ceres_figure_df <- ceres_df %>%
    dplyr::filter(barcode_control == "cutting_control", cor_method == "pearson") %>%
    dplyr::left_join(
        count_df,
        by = c("perturbation" = "pert_name", "group" = "gene_name", "cell_line" = "cell_line")
    )

head(ceres_figure_df)

ceres_panel_gg = (
    ggplot(
        ceres_figure_df,
        aes(x = grit, y = ceres_score)
    )
    + geom_point(aes(fill = log10(cell_count), size = log10(cell_count), shape = cell_line))
    + geom_vline(xintercept = 1, linetype = "dashed", color = "blue")
    + annotate("rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf, alpha = .2)
    + scale_shape_manual(name = "Cell line", values = c(21, 23, 25))
    + scale_fill_gradient(name = "log10\ncell count", low = "white", high = "red")
    + scale_size_continuous(guide = FALSE, range = c(0.2, 2.5))
    + xlab("Grit")
    + ylab("CERES")
    + custom_grit_theme
    + labs(tag = "d")
)

ceres_panel_gg

gene <- "ITGAV"

focus_plot_df <- sc_umap_df %>%
    dplyr::filter(grit_gene == !!gene) %>%
    dplyr::left_join(
        sc_grit_df %>%
            dplyr::filter(grit_gene == !!gene),
        by = c("Metadata_cell_identity" = "perturbation", "grit_gene" = "grit_gene")
    )

head(focus_plot_df)

control_perts <- focus_plot_df %>%
    dplyr::filter(Metadata_gene_name == "Chr2") %>%
    dplyr::select(Metadata_pert_name) %>%
    unique() %>%
    dplyr::pull(Metadata_pert_name)

gene_perts <- focus_plot_df %>%
    dplyr::filter(Metadata_gene_name == !!gene) %>%
    dplyr::select(Metadata_pert_name) %>%
    unique() %>%
    dplyr::pull(Metadata_pert_name)

pert_order <- c(gene_perts, sort(control_perts))

focus_plot_df$Metadata_pert_name <- factor(
    focus_plot_df$Metadata_pert_name, levels = pert_order
)

focus_plot_df <- focus_plot_df %>%
    dplyr::filter(Metadata_pert_name %in% pert_order)

gene_umap_gg = (
    ggplot(focus_plot_df, aes(x = umap_0, y = umap_1))
    + geom_point(aes(fill = grit), pch = 21, size = 0.6, stroke = 0.1, alpha = 0.7)
    + facet_wrap("~Metadata_pert_name", nrow = 2)
    + scale_fill_gradientn(
        name = "scGrit",
        colors = scgrit_palette(10),
        values = scales::rescale(c(2, 1.4, 1.1)),
        limits = c(-0.7, 2.1)
    )
    + xlab("UMAP 0")
    + ylab("UMAP 1")
    + custom_grit_theme
    + labs(tag = "e")
)

gene_umap_gg

patchwork_plot <- (
    (
        panel_a_gg / ceres_panel_gg
    ) + plot_layout(heights = c(0.5, 1)) | (
        (
            mean_grit_gg + replicate_reprod_panel
        ) / gene_umap_gg
    )
) + plot_layout(widths = c(1, 2))


for (extension in c(".png", ".pdf")) {
    output_file <- file.path(output_dir, paste0("figure1", extension))
    ggsave(output_file, patchwork_plot, width = 12, height = 6.5)
}

patchwork_plot
