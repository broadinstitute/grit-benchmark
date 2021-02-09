# For consistent themes and colors in the grit benchmark paper

library(ggplot2)
library(RColorBrewer)

# Cell health cell line colors
cell_line_colors <- c(
  "A549" = "#861613",
  "ES2" = "#1CADA8",
  "HCC44" = "#2A364D"
)

# Standard ggplot theme
custom_grit_theme <- (
    theme_bw()
    + theme(
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        strip.background = element_rect(color = "black", fill = "#fdfff4"),
        strip.text = element_text(size = 7)
    )
)

scgrit_palette <- colorRampPalette(rev(brewer.pal(9, "PuBuGn")))