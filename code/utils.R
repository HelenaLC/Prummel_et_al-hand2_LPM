# canonical marker genes -------------------------------------------------------
mgs_canonical <- readLines(file.path("data", "marker_genes.txt"))
mgs_canonical <- strsplit(mgs_canonical, ", ")[[1]]

# cluster color palette --------------------------------------------------------
library(RColorBrewer)
cluster_pal <- c(
    "endoderm_1" = brewer.pal(9, "Blues")[8],
    "endoderm_2" = brewer.pal(9, "Blues")[6],
    "endoderm_3" = brewer.pal(9, "Blues")[4],
    "hand2-high_1" = brewer.pal(9, "Reds")[9],
    "hand2-high_2" = brewer.pal(9, "Reds")[7],
    "hand2-high_3" = brewer.pal(9, "Reds")[5],
    "hand2-high_4" = brewer.pal(9, "Reds")[3],
    "head mesoderm_1" = brewer.pal(9, "Greens")[8],
    "head mesoderm_2" = brewer.pal(9, "Greens")[5],
    "anterior vasculature" = "gold",
    "cardiopharyngeal_1" = brewer.pal(9, "Purples")[8],
    "cardiopharyngeal_2" = brewer.pal(9, "Purples")[5],
    "heart field" = "plum",
    "hemangioblasts and kidney_1" = "lightseagreen",
    "hemangioblasts and kidney_2" = "aquamarine")

# 0-to-1 scaling using (1%)-/(99%)-quantiles as boundaries ---------------------
.scale <- function(x) {
    if (!is(x, "matrix")) x <- as.matrix(x)
    qs <- rowQuantiles(x, probs = c(0.01, 0.99))
    qs <- matrix(qs, ncol = 2)
    x <- (x - qs[, 1]) / (qs[, 2] - qs[, 1])
    x[x < 0] <- 0
    x[x > 1] <- 1
    return(x)
}

# dotplot with rows = clusters & columns = genes
# where points are colored by mean scaled expression & 
# sized by expression frequency, e.g., mean(count > 0)
#   - so: a SeuratObject
#   - gs: features to plot
#   - os: list of row and column orders
.dotplot <- function(so, gs) {
    p <- DotPlot(so, features = gs)
    m <- acast(p$data, features.plot ~ id, value.var = "avg.exp.scaled")
    fun <- function(m) order.dendrogram(as.dendrogram(hclust(dist(m))))
    ro <- colnames(m)[fun(t(m))]
    co <- rownames(m)[fun(m)]
    p + scale_y_discrete(limits = ro, expand = c(0, 0.5)) +
        scale_x_discrete(limits = co, expand = c(0, 0.5),
            labels = so[["RNA"]]@meta.features[co, "symbol"]) +
        scale_color_gradientn(
            breaks = seq(floor(min(m)), ceiling(max(m))),
            limits = c(floor(min(m)*4)/4, ceiling(max(m)*4)/4), 
            colors = c("royalblue", "grey95", "yellowgreen", "darkgreen")) +
        scale_size_continuous(labels = formatC(seq(0,1,0.25), 2, format = "f")) +
        guides(size = guide_legend("expression\nfrequency"),
            color = guide_colorbar("scaled mean\nexpression", order = 1)) +
        coord_equal() + theme_void() + theme(
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            axis.text.y = element_text(hjust = 1))
}
        