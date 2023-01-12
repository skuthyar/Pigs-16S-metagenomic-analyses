Code associated with "Pig microbial and immune traits respond to host evolutionary history and ecology across multiple scales."

Description: These files will allow you to replicate the bioinformatic and statistical analyses we used to assess the impact of domestication on the pig gut microbiome and immune state. 

The following files should be used in order to process raw 16S rRNA gene sequencing data: 
1. Cutadapt to remove primers
2. Dada2 to learn errors, dereplicate, infer ASVs, merge amplicons, make sequence tables, and assign taxonomy
3. Decontam to remove contaminants, mitochondria, and chloroplasts
4. Tables to create final counts and taxonomy tables
5. Tree to create a phylogenetic tree

The output will be an ASV decontaminated counts table, an ASV taxonomy table, and a phylogenetic tree. The R markdown file (16S_analyses.Rmd) will allow you to replicate any statistical analyses and visualizations for 16S rRNA gene sequencing. The excel file, fecal_metadata_Apr2022.xlsx, contains all the metadata needed for these analyses.

###############################################################################################

The following files should be used in order to process raw shotgun metagenomic data: 
1. Quality filtering and trimming using Fastp
2. Map sequences against reference genome
3. Assemble to contigs using megahit
4. Gene prediction with prodigal
5. Dereplication with cd-hit
6. Align to UniProt and taxonomic classification
7. Annotation with EggNog, CARD, and VFDB
8. Salmon gene quantification

Outputs will be annotation files and quantification files. The .RMD file will allow for replication of statistical analyses and visualizations. The excel file, "dom_cat_metagenomic_metadata.xlsx," contains all the metadata associated with the samples.

The SNP.R file along with the plink.vcf and SNP_metadata.xlsx files can be used to replicate all SNP analyses. 
