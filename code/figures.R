# load packages
suppressMessages({
    library(circlize)
    library(ComplexHeatmap)
    library(dplyr)
    library(ggplot2)
    library(magrittr)
    library(muscat)
    library(RColorBrewer)
    library(Seurat)
    library(SingleCellExperiment)
})

# load data & metadata
sce <- readRDS(file.path("output", "LPM-sce.rds"))
cluster_anno <- readxl::read_xlsx(file.path("data", "cluster_anno.xlsx"))

cluster_cols <- c(
    "unknown" = "grey",
    "endoderm_1" = brewer.pal(9, "Blues")[8],
    "endoderm_2" = brewer.pal(9, "Blues")[6],
    "endoderm_3" = brewer.pal(9, "Blues")[4],
    "hand2-high_1" = brewer.pal(9, "Reds")[8],
    "hand2-high_2" = brewer.pal(9, "Reds")[6],
    "hand2-high_3" = brewer.pal(9, "Reds")[4],
    "pharyngeal/head mesoderm_1" = brewer.pal(9, "Greens")[8],
    "pharyngeal/head mesoderm_2" = brewer.pal(9, "Greens")[6],
    "pharyngeal/head mesoderm_3" = brewer.pal(9, "Greens")[4],
    "anterior vasculature" = "gold",
    "heart field_1" = brewer.pal(9, "Purples")[8],
    "heart field_2" = brewer.pal(9, "Purples")[6],
    "heart field_3" = "plum",
    "posterior hemangioblasts and kidney_1" = "mediumaquamarine",
    "posterior hemangioblasts and kidney_2" = "aquamarine")

all(cluster_anno$label %in% names(cluster_cols))

anno <- setNames(cluster_anno$label, cluster_anno$cluster_number)
anno <- anno[as.character(sce$cluster_id)]
sce$cluster_id <- droplevels(factor(anno, levels = names(cluster_cols)))

dr <- do.call(cbind, reducedDims(sce))
colnames(dr) <- gsub("*\\.", "", colnames(dr))

df <- data.frame(t(logcounts(sce)), colData(sce), dr, check.names = FALSE)
df <- df[sample(nrow(df), nrow(df)), ]
df[, rownames(sce)] <- t(muscat:::.scale(t(df[, rownames(sce)])))

# UMAP colored by cluster ID ---------------------------------------------------
.plot_umap <- function(df, col, pal) {
    if (missing(pal)) pal <- cluster_cols
    p <- ggplot(df, aes_string(col = sprintf("`%s`", col), 
        x = "UMAP_1", y = "UMAP_2")) + geom_point(size = 0.2) + 
        guides(color = guide_legend(override.aes = list(size = 2))) +
        theme_void(base_size = 10) + theme(aspect.ratio = 1, 
            legend.key.size = unit(2, "mm"))
    if (is.factor(df[[col]])) {
        p + scale_color_manual(values = pal)
    } else {
        #qs <- quantile((es <- df[, col]), c(.01, .99))
        #es[es < qs[1]] <- qs[1]
        #es[es > qs[2]] <- qs[2]
        p + scale_color_gradientn(NULL,
            colors = c("black", "darkcyan", "white", "orange", "brown"))
            #theme(legend.position = "none")
    }
}

p <- .plot_umap(df, "cluster_id")
ggsave(file.path("figures", "UMAP-cluster_id.pdf"), p,
    width = 12, height = 6, units = "cm",
    dpi = 300, useDingbats = FALSE)

# cluster-markers dotplot ------------------------------------------------------
mgs <- readxl::read_xlsx(file.path("data", "marker_genes.xlsx"))
ss <- strsplit(mgs$symbol, ", ")
mgs <- setNames(ss, mgs$cluster)

cs <- split(colnames(sce), sce$cluster_id)
fs <- unique(unlist(mgs))
fs <- fs[!is.na(fs)]
gs <- sapply(fs, grep, rownames(sce))

rmv <- sapply(gs, length) == 0
fs <- fs[!rmv]
gs <- unlist(gs[!rmv])

ps <- sapply(cs, function(i) 
    rowMeans(counts(sce)[gs, i, drop = FALSE] > 0)) %>% 
    set_rownames(fs)
ms <- sapply(cs, function(i) 
    rowMeans(logcounts(sce)[gs, i, drop = FALSE])) %>% 
    set_rownames(fs)
ms0 <- muscat:::.scale(ms)

d <- dist(t(ms))
h <- hclust(d)
col_o <- h$order
row_o <- unlist(mgs[col_o])
row_o <- make.names(row_o)
row_o <- row_o[!duplicated(row_o)]
row_o <- row_o[row_o %in% rownames(ms)]
row_o <- match(row_o, fs, nomatch = 0)

gg_df <- cbind(
    reshape2::melt(ms0), 
    p = reshape2::melt(ps)$value) %>% 
    mutate_if(is.factor, as.character)

p <- ggplot(gg_df, aes(x = Var1, y = Var2, col = value, size = p)) +
    geom_point() +
    scale_size_continuous(range = c(0, 3)) +
    scale_color_gradientn("scaled mean\nexpression", breaks = seq(0, 1, 0.2), 
        colors = c("black", "darkcyan", "white", "orange", "brown")) +
    guides(size = guide_legend(order = 1, "expression\nfrequency")) +
    scale_x_discrete(limits = rownames(ms)[row_o], position = "top", expand = c(0, 0.5)) +
    scale_y_discrete(limits = colnames(ms)[rev(col_o)], expand = c(0, 0.5)) +
    coord_equal(clip = "off") + theme_minimal(base_size = 8) + theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 0, vjust = 1),
        legend.key.size = unit(3, "mm"))

ggsave(file.path("figures", "marker_dotplot.pdf"), p,
    width = 16, height = 6, units = "cm",
    dpi = 300, useDingbats = FALSE)

# UMAPs colored by expression of canonical markers -----------------------------
gs <- sapply(unlist(mgs), grep, colnames(df), value = TRUE)
gs <- gs[sapply(gs, function(g) length(g) > 0)]
for (g in unique(gs)) {
    p <- .plot_umap(df, g)
    ggsave(file.path("figures", "UMAPs", sprintf("umap-%s.pdf", g)), p,
        width = 5, height = 5, units = "cm",
        dpi = 300, useDingbats = FALSE)
}

# violins ----------------------------------------------------------------------
fs <- unlist(mgs)
fs <- fs[!is.na(fs)]
fs <- unique(c("mCherry", "hand2", fs))
gs <- sapply(fs, grep, colnames(df), value = TRUE)

for (g in gs) {
    max <- ceiling(max(df[, g])*0.5)/0.5
    p <- ggplot(df, aes_string(x = "cluster_id", y = sprintf("`%s`", g), 
        col = "cluster_id", fill = "cluster_id")) +
        geom_violin(size = 0.4, alpha = 0.4) + 
        geom_boxplot(size = 0.2, width = 0.2, outlier.size = 0, fill = "white") +
        scale_color_manual(values = cluster_cols) +
        scale_fill_manual(values = cluster_cols) +
        scale_x_discrete(NULL, expand = c(0, 1)) +
        scale_y_continuous("expression", limits = c(0, max), expand = c(0, 0.25)) +
        theme_bw() + theme_linedraw(base_size = 8) + ggtitle(g) +
        theme(aspect.ratio = 1/2, legend.position = "none",
            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(size = 0.2, color = "lightgrey"))
    
    ggsave(file.path("figures", "violins", sprintf("%s.pdf", g)), p,
        width = 9, height = 5, units = "cm",
        dpi = 300, useDingbats = FALSE)
}

