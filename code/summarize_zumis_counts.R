# Summarize zUMIs counts into a single count matrix

# list samples
sample_ids <- list.dirs(file.path("..", "data"))
sample_ids <- c("SN1_AHY3WGBGX3_S5", "SN2_AHY3WGBGX3_S6",
                "SN3_AHY3WGBGX3_S7", "SN4_AH2TW3BGX5_S1")
names(sample_ids) <- sapply(strsplit(sample_ids, "_"), .subset, 1)

## Get UMI counts
counts <- Reduce(function(...) dplyr::full_join(..., by = "GENEID"), 
                 lapply(sample_ids, function(w) {
                   x <- readRDS(paste0("zUMIs/", w, "/zUMIs_output/expression/zUMIs_output/expression/", w, ".dgecounts.rds"))
                   x <- as.data.frame(as.matrix(x$exon$umicounts))
                   colnames(x) <- paste0(strsplit(w, "_")[[1]][1], "_", colnames(x))
                   x$GENEID = rownames(x)
                   x
                 })
)

## Clean up
counts[is.na(counts)] <- 0
rownames(counts) <- counts$GENEID
counts$GENEID <- NULL
counts <- as.matrix(counts)
dim(counts)

## Save 
saveRDS(counts, file="zUMIs/umicounts_merged.rds")

