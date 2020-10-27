args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(cdna)
print(ncrna)
print(custom_transcripts)
print(outfile)

suppressPackageStartupMessages(library(Biostrings))

cdna <- readDNAStringSet(cdna)
ncrna <- readDNAStringSet(ncrna)
custom_transcripts <- readDNAStringSet(custom_transcripts)
fasta <- c(cdna, ncrna, custom_transcripts)

tx2gene <- data.frame(t(sapply(names(fasta), function(nm) {
  tmp <- strsplit(nm, " ")[[1]]
  tx <- tmp[1]
  gene <- gsub("gene:", "", tmp[grep("^gene:", tmp)])
  c(tx = tx, gene = gene)
})), stringsAsFactors = FALSE)
rownames(tx2gene) <- NULL
head(tx2gene)

write.table(tx2gene, file = gsub("rds$", "txt", outfile), 
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
saveRDS(tx2gene, file = outfile)

sessionInfo()
date()