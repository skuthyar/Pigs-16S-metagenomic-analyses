# Snakefile

# Define paths and resources
FASTQ_DIR = "/projects/ps-aspenlab/data/PATH_TO_FASTQ_FILES"
OUTPUT_DIR = "/projects/ps-aspenlab/LAST_NAME/PROJECT"
SAMPLES_FILE = "/projects/ps-aspenlab/samples.txt"
TRIMMED_DIR = f"{OUTPUT_DIR}/trimmed"
FILTERED_DIR = f"{OUTPUT_DIR}/filtered"
SAVED_OUTPUT_DIR = "/projects/ps-aspenlab/PATH_TO_SAVED_OUTPUT_FOLDER"
SILVA_DB = "~/databases/silva_nr99_v138.1_wSpecies_train_set.fa.gz"

# Rule to prepare environment
rule prepare_conda_env:
    output:
        "conda_env_ready.txt"
    shell:
        """
        source ~/.bashrc
        conda activate cutadaptenv
        echo "conda environment activated" > {output}
        """

# Rule to unzip fastq files
rule unzip_fastq:
    input:
        fastq_files=expand(f"{FASTQ_DIR}/*.fastq.gz",)
    output:
        unzipped_fastqs=expand(f"{FASTQ_DIR}/{wildcards.sample}_R1_001.fastq", sample=SAMPLES_FILE)
    shell:
        """
        gunzip {input.fastq_files}
        """

# Rule to create sample list file
rule create_samples_list:
    input:
        fastq_files=expand(f"{FASTQ_DIR}/*_R1_001.fastq")
    output:
        samples_txt=SAMPLES_FILE
    shell:
        """
        ls {input.fastq_files} | cut -f 1 -d "_" > {output}
        """

# Rule to create output directories
rule create_output_dirs:
    output:
        trimmed_dir=TRIMMED_DIR,
        filtered_dir=FILTERED_DIR
    shell:
        """
        mkdir -p {output.trimmed_dir} {output.filtered_dir}
        """

# Rule for cutadapt trimming
rule trim_reads:
    input:
        r1=f"{FASTQ_DIR}/{wildcards.sample}_L001_R1_001.fastq",
        r2=f"{FASTQ_DIR}/{wildcards.sample}_L001_R2_001.fastq"
    output:
        r1_trimmed=f"{TRIMMED_DIR}/{wildcards.sample}_trimmed_R1.fastq",
        r2_trimmed=f"{TRIMMED_DIR}/{wildcards.sample}_trimmed_R2.fastq"
    params:
        forward_adapter="GTGCCAGCMGCCGCGGTAA",
        reverse_adapter="GGACTACHVGGGTWTCTAAT"
    shell:
        """
        cutadapt -g {params.forward_adapter} -G {params.reverse_adapter} \
        -o {output.r1_trimmed} -p {output.r2_trimmed} \
        {input.r1} {input.r2} || exit 1
        """

# Rule to check if trimming was successful
rule check_trimming:
    input:
        trimmed_r1=f"{TRIMMED_DIR}/{wildcards.sample}_trimmed_R1.fastq",
        trimmed_r2=f"{TRIMMED_DIR}/{wildcards.sample}_trimmed_R2.fastq"
    output:
        "trim_check_{wildcards.sample}.txt"
    shell:
        """
        if [[ -s {input.trimmed_r1} && -s {input.trimmed_r2} ]]; then
            echo "Trimming successful for {wildcards.sample}" > {output}
        else
            echo "Trimming failed for {wildcards.sample}" > {output}
            exit 1
        fi
        """

# Rule to gzip sequence files back
rule gzip_fastq:
    input:
        fastqs=expand(f"{FASTQ_DIR}/*.fastq",)
    output:
        gzipped_fastqs=expand(f"{FASTQ_DIR}/*.fastq.gz",)
    shell:
        """
        gzip {input.fastqs} || exit 1
        """

# R script quality control plotting
rule plot_quality:
    input:
        trimmed_r1=expand(f"{TRIMMED_DIR}/{wildcards.sample}_trimmed_R1.fastq", sample=SAMPLES_FILE),
        trimmed_r2=expand(f"{TRIMMED_DIR}/{wildcards.sample}_trimmed_R2.fastq", sample=SAMPLES_FILE)
    output:
        "quality_profiles.pdf"
    script:
        "scripts/plot_quality.R"  # Assuming R script exists to generate quality plots

# Rule to filter and trim using dada2
rule filter_trim_dada2:
    input:
        trimmed_r1=expand(f"{TRIMMED_DIR}/{wildcards.sample}_trimmed_R1.fastq", sample=SAMPLES_FILE),
        trimmed_r2=expand(f"{TRIMMED_DIR}/{wildcards.sample}_trimmed_R2.fastq", sample=SAMPLES_FILE)
    output:
        filtered_r1=f"{FILTERED_DIR}/{wildcards.sample}_filtered_R1.fastq",
        filtered_r2=f"{FILTERED_DIR}/{wildcards.sample}_filtered_R2.fastq"
    shell:
        """
        Rscript -e '
        library(dada2);
        filtFs <- c("{input.trimmed_r1}");
        filtRs <- c("{input.trimmed_r2}");
        filterAndTrim(filtFs, "{output.filtered_r1}", filtRs, "{output.filtered_r2}", maxEE=c(2,2), rm.phix=TRUE, minLen=175, truncLen=c(230,230));
        ' || exit 1
        """

# Rule for DADA2 error learning
rule learn_errors:
    input:
        filtered_r1=f"{FILTERED_DIR}/{wildcards.sample}_filtered_R1.fastq",
        filtered_r2=f"{FILTERED_DIR}/{wildcards.sample}_filtered_R2.fastq"
    output:
        err_r1=f"{OUTPUT_DIR}/{wildcards.sample}_err_r1.rds",
        err_r2=f"{OUTPUT_DIR}/{wildcards.sample}_err_r2.rds"
    shell:
        """
        Rscript -e '
        library(dada2);
        errF <- learnErrors("{input.filtered_r1}", multithread = TRUE);
        errR <- learnErrors("{input.filtered_r2}", multithread = TRUE);
        saveRDS(errF, "{output.err_r1}");
        saveRDS(errR, "{output.err_r2}");
        ' || exit 1
        """

# Rule to run DADA2 pipeline (dereplication, inference, merging, etc.)
rule dada2_pipeline:
    input:
        err_r1=f"{OUTPUT_DIR}/{wildcards.sample}_err_r1.rds",
        err_r2=f"{OUTPUT_DIR}/{wildcards.sample}_err_r2.rds",
        filtered_r1=f"{FILTERED_DIR}/{wildcards.sample}_filtered_R1.fastq",
        filtered_r2=f"{FILTERED_DIR}/{wildcards.sample}_filtered_R2.fastq"
    output:
        merged_amplicons=f"{OUTPUT_DIR}/{wildcards.sample}_merged_amplicons.rds",
        sequence_table=f"{OUTPUT_DIR}/{wildcards.sample}_sequence_table.rds"
    shell:
        """
        Rscript -e '
        library(dada2);
        errF <- readRDS("{input.err_r1}");
        errR <- readRDS("{input.err_r2}");
        filtFs <- c("{input.filtered_r1}");
        filtRs <- c("{input.filtered_r2}");
        derepF <- derepFastq(filtFs, verbose = TRUE);
        derepR <- derepFastq(filtRs, verbose = TRUE);
        dadaF <- dada(derepF, err = errF, pool = "pseudo", multithread = TRUE);
        dadaR <- dada(derepR, err = errR, pool = "pseudo", multithread = TRUE);
        merged_amplicons <- mergePairs(dadaF, derepF, dadaR, derepR, verbose = TRUE);
        saveRDS(merged_amplicons, "{output.merged_amplicons}");
        sequencetable <- makeSequenceTable(merged_amplicons);
        saveRDS(sequencetable, "{output.sequence_table}");
        ' || exit 1
        """

# Add more rules for the rest of the pipeline as necessary

# Final cleanup and taxonomy assignment
rule final_steps:
    input:
        sequence_table=f"{OUTPUT_DIR}/{wildcards.sample}_sequence_table.rds"
    output:
        final_output="final_output.txt"
    shell:
        """
        Rscript -e '
        library(dada2);
        sequence_table <- readRDS("{input.sequence_table}");
        taxonomy <- assignTaxonomy(sequence_table, "{SILVA_DB}", multithread = TRUE);
        saveRDS(taxonomy, "taxonomy.rds");
        write.table(taxonomy, "ASV_taxonomy.tsv", sep = "\t", quote = F, col.names = NA);
        ' || exit 1
        """

# Makefile based on your input
rule all:
    input:
        "final_output.txt"

