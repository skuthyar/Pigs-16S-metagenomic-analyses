# decontaminate seqs: run through the 4.0.2 version of R instead of default 
library(dplyr) 
packageVersion("dvplyr")
library(decontam)
packageVersion("decontam")

# load tables
counts <- read.table(file = 'ASV_counts.tsv', sep = '\t', header = TRUE, row.names = 1)
tax_table <- read.table(file = 'ASV_taxonomy.tsv', sep = '\t', header = TRUE, row.names = 1)
asv_fasta <- readRDS("ASV_fasta.rds")

# counts colnames should be a list of the samples, tax_table colnames should be kingdom phylum class etc...
print("Counts columns")
colnames(counts)
print("Tax columns")
colnames(tax_table)

# create a list of T/F values where TRUE = negatives, FALSE = samples or positives (from sample.names)
# rep is the number of samples that are TRUE or FALSE in order 
vector_for_decontam <- c(rep(FALSE,107), rep(TRUE,4), rep(FALSE,1), rep(TRUE,2), rep(FALSE,54), rep(TRUE,2), rep(FALSE,1), rep(TRUE,1), rep(FALSE,9), rep(TRUE,1), rep(FALSE,1), rep(TRUE,1), rep(FALSE,40), rep(TRUE,2), rep(FALSE,49), rep(TRUE,6), rep(FALSE,24), rep(TRUE,16), rep(FALSE,51))

# determine contaminants and store the ASV numbers identified as contam_asvs
contam_df <- isContaminant(t(counts), neg = vector_for_decontam)
# table output: how many of the ASVs were contaminants
print("Contaminants")
table(contam_df$contaminant)
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])

# display the taxonomy of each of the ASVs that were identified
tax_table[row.names(tax_table) %in% contam_asvs, ]

saveRDS(contam_asvs, "contam_asvs.rds")

# remove chloroplast, mitochondria, no phylum ASVs

tax_table <- read.table(file = 'ASV_taxonomy.tsv', sep = '\t', header = TRUE, row.names = 1)
contam_asvs <- readRDS("contam_asvs.rds")
asv_fasta <- readRDS("ASV_fasta.rds")
counts <- read.table(file = 'ASV_counts.tsv', sep = '\t', header = TRUE, row.names = 1)

chloro <- row.names(filter(tax_table, tax_table$Phylum == "Cyanobacteria" & tax_table$Order == "Chloroplast"))
mito <- row.names(filter(tax_table, tax_table$Class == "Alphaproteobacteria" & tax_table$Family == "Mitochondria"))
nobac <- row.names(filter(tax_table, tax_table$Kingdom != "Bacteria" | is.na(tax_table$Phylum)))
to.remove = c(contam_asvs, chloro, mito, nobac)
saveRDS(to.remove, "asvs_to_remove.rds")
print("Number of ASVs to remove")
length(to.remove)

# remove contaminant ASVs from the ASV fasta file
contam_indices <- which(asv_fasta %in% paste0(">", to.remove))
dont_want <- sort(c(contam_indices, contam_indices + 1))
asv_fasta_decontam <- asv_fasta[- dont_want]
write(asv_fasta_decontam, "ASVs_decontam.fa")
saveRDS(asv_fasta_decontam, "ASVs_decontam.rds")

# remove contaminant ASVs from counts table
counts_decontam <- counts[!row.names(counts) %in% to.remove, ]
write.table(counts_decontam, "ASV_counts_decontam.tsv", sep = "\t", quote=F, col.names = NA)

# remove contaminant ASVs from taxonomy table
tax_table_decontam <- tax_table[!row.names(tax_table) %in% to.remove, ]
write.table(tax_table_decontam, "ASV_taxonomy_decontam.tsv", sep = "\t", quote=F, col.names = NA)

export MODULEPATH=/projects/builder-group/jpg/modulefiles/applications:$MODULEPATH
module load R/4.0.2
Rscript --no-save ~/scripts/decontam.R
