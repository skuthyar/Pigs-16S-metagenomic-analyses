# Activate conda and cutadapt, if there are problems the culprit is probably one of these two lines
source ~/.bashrc
conda activate cutadaptenv

# Change directory to easily access fastq sequence files
cd /projects/ps-aspenlab/data/<<PATH_TO_FASTQ_FILES>>

# Unzip sequence files
gunzip *.fastq.gz

# Creates text file with list of all the sample names
# This line will extract the first field using delimiter "_" of all files ending in _R1_001.fastq, if files are in different format may need adjusting
ls *_R1_001.fastq | cut -f 1 -d "_" > <<PATH_TO_SAVED_OUTPUT_FOLDER>>/samples.txt

# Change to lab directory and create output file in personal folder
cd /projects/ps-aspenlab
mkdir <<LAST_NAME>>/<<PROJECT>>/trimmed
mkdir <<LAST_NAME>>/<<PROJECT>>/filtered

# Works for samples with the file names in format: SAMPLEID_SAMPLENUMBER_L001_R(1/2)_001.fastq, may need to adjust
for SAMPLEID in $(cat <<PATH_TO_SAMPLES_FILE>>);

do
        echo "On sample: $SAMPLEID"
        cutadapt -g GTGCCAGCMGCCGCGGTAA -G GGACTACHVGGGTWTCTAAT -o <<PATH_TO_TRIMMED>>/trimmed/${SAMPLEID}_trimmed_R1.fastq -p <<PATH_TO_TRIMMED>>/trimmed/${SAMPLEID}_trimmed_R2.fastq <<PATH_TO_FASTQ_FILES>>/${SAMPLEID}_*_L001_R1_001.fastq <<PATH_TO_FASTQ_FILES>>/${SAMPLEID}_*_L001_R2_001.fastq
done

# Re-zip sequence files
gzip /projects/ps-aspenlab/<<PATH_TO_FASTQ_FILES>>/*.fastq

# final output is a set of trimmed fastq files in a "trimmed" folder
