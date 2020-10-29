# summarize zUMIs counts into single count matrix

# list samples
fns <- list.files(file.path("data", "raw"), full.names = TRUE)
ids <- gsub("\\..*", "", basename(fns))
names(fns) <- ids <- sapply(strsplit(ids, "_"), .subset, 1)

# get UMI counts
y <- lapply(ids, function(id) {
    x <- readRDS(fns[[id]])
    x <- x$exons$umicounts
    x <- as.matrix(x)
    x <- as.data.frame(x)
    colnames(x) <- paste(id, colnames(x), sep = "_")
    x$GENEID <- rownames(x)
    return(x)
})
y <- Reduce(function(...) dplyr::full_join(..., by = "GENEID"), y)

# clean up
y[is.na(y)] <- 0
rownames(y) <- y$GENEID
y$GENEID <- NULL
y <- as.matrix(y)

# write to .rds
saveRDS(y, file.path("data", "UMI_counts_merged.rds"))