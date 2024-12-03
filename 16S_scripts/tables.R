# load sequence table
sequencetable.nochim <- readRDS("sequencetable.nochim.rds")

# create ASV fasta file
# store sequences and asv headers to use in manipulating tables
asv_seqs <- colnames(sequencetable.nochim)
asv_headers <- vector(dim(sequencetable.nochim)[2], mode = "character")

# format headers for fasta file
for (i in 1:dim(sequencetable.nochim)[2]) {
        asv_headers[i] <- paste(">ASV", i, sep = "_")
}

# create and write fasta table with sequences for each ASV
asv_fasta <- c(rbind(asv_headers, asv_seqs))
saveRDS(asv_fasta, "ASV_fasta.rds")
write(asv_fasta, "ASVs.fa")

# create counts table by transposing sequence table and changing row names
asv_tab <- t(sequencetable.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASV_counts.tsv", sep = "\t", quote = F, col.names = NA)

# create taxonomy table
# remove > from headers
asv_headers_tax <- sub('.', '', asv_headers)

taxonomy <- readRDS("taxonomy.rds")

# change row names from sequences to ASV numbers
rownames(taxonomy) <- asv_headers_tax

write.table(taxonomy, "ASV_taxonomy.tsv", sep = "\t", quote = F, col.names = NA)

module load R
Rscript --no-save ~/scripts/tables.R
