# assemble contigs using megahit 
conda activate megahitenv

for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID"
        megahit -1 ~/processed_seqs/${SAMPLEID}_final_R1.fastq.gz -2 ~/processed_seqs/${SAMPLEID}_final_R2.fastq.gz --num-cpu-threads 8 --min-contig-len 500 -o ~/megahit/${SAMPLEID}_megahit
done     

# build index for contigs 
module load bowtie2

for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID"
        bowtie2-build ~/megahit/${SAMPLEID}_megahit/final.contigs.fa megahit_index
done

# align sequencing reads to assembled contigs
module load bowtie2

for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID"
        bowtie2 -p 15 -t -x megahit_index -1 ~/processed_seqs/${SAMPLEID}_final_R1.fastq.gz -2 ~/processed_seqs/${SAMPLEID}_final_R2.fastq.gz -S ${SAMPLEID}_aligned.sam
done

# extract unmapped sequence  
module load samtools
module load bedtools

for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID"
        samtools view -b -f 12 -F 256 ${SAMPLEID}_aligned.sam > ${SAMPLEID}_unmapped.bam
        samtools sort -n ${SAMPLEID}_unmapped.bam -o ${SAMPLEID}_unmapped_sorted.bam
        bedtools bamtofastq -i ${SAMPLEID}_unmapped_sorted.bam -fq ${SAMPLEID}_unmapped_r1.fastq -fq2 ${SAMPLEID}_unmapped_r2.fastq
done

# coassemble unmapped reads
conda activate megahitenv

for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID"
        megahit -1 ${SAMPLEID}_unmapped_r1.fastq -2 ${SAMPLEID}_unmapped_r2.fastq --min-count 2 --k-min 27 --k-max 87 --k-step 10 --num-cpu-threads 20 --min-contig-len 500 -o ~/megahit/${SAMPLEID}_unmapped_megahit
done

# merge all contigs 
for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID"
        cat ${SAMPLEID}_megahit/final.contigs.fa ${SAMPLEID}_unmapped_megahit/final.contigs.fa > ${SAMPLEID}_all_final_contigs.fasta
done
