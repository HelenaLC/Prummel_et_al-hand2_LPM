# load required packages
suppressPackageStartupMessages({
    library(circlize)
    library(ComplexHeatmap)
    library(ggplot2)
    library(matrixStats)
    library(scales)
    library(viridis)
})

# load utility functions
source(file.path("code", "utils.R"))

# load RNA-seq data
z <- read.csv(file.path("data", "Rehrauer18.csv"), 
    header = TRUE, row.names = 1, check.names = FALSE)
z <- as.matrix(z)

# construct row & column annotations
col_anno <- columnAnnotation(
    foo = anno_block(
        gp = gpar(fill = hue_pal()(3)),
        labels = c(
            "crocid tum" = "Crocidotlie-exposed\ntumor (mesothelioma)",
            "crocid" = "Crocidolite-exposed\nmesothelium",
            "sham" = "Sham-exposed mesothelium"),
        labels_gp = gpar(col = "white", fontsize = 6)))

df <- data.frame(
    min = log10(rowMins(z)+1),
    max = log10(rowMaxs(z)+1))
qs <- quantile(unlist(df), c(0,.5,1))
anno_pal <- colorRamp2(qs, c("grey95", "gold", "red3"))
row_anno <- rowAnnotation(df = df,
    simple_anno_size = unit(4.5, "mm"),
    col = list(min = anno_pal, max = anno_pal),
    gp = gpar(col = "white"),
    show_legend = c(TRUE, FALSE),
    annotation_name_gp = gpar(fontsize = 8),
    annotation_legend_param = list(
        min = list(
            title = "min / max\ncount (log10)",
            title_gp = gpar(fontsize = 8),
            labels_gp = gpar(fontsize = 6))))

mat <- t(scale(t(log10(z+1))))
min <- floor(min(mat)*2)/2
max <- ceiling(max(mat)*2)/2
hm_pal <- c("navy", "royalblue", "grey95", "yellowgreen", "darkgreen")
hm_pal <- colorRamp2(c(min, min/2, 0, max/2, max), hm_pal)

# plot heatmap & save to .pdf
hm <- Heatmap(
    matrix = mat,
    col = hm_pal,
    row_names_side = "left",
    top_annotation = col_anno,
    right_annotation = row_anno,
    rect_gp = gpar(col = "white"),
    row_names_gp = gpar(fontsize = 8),
    heatmap_legend_param = list(
        title = "scaled normalized\ncount (log10)",
        title_gp = gpar(fontsize = 8),
        labels_gp = gpar(fontsize = 6)),
    column_split = colnames(z),
    column_title = NULL,
    show_column_names = FALSE)

hm <- grid.grabExpr(draw(hm))
fn <- file.path("figures", "Rehrauer18_heatmap.pdf")
ggsave(fn, hm, width = 15, height = 10.5, units = "cm")
