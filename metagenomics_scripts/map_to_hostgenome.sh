# map to pig reference genome  
module load bwa
module load samtools 

for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID" 
        bwa mem -t 8 -T 30 ~/pigref ~/clean/${SAMPLEID}_R1_001.clean.fastq.gz ~/clean/${SAMPLEID}_R2_001.clean.fastq.gz -R '@RG\tID:${SAMPLEID}_join\tPL:illumina\tSM:${SAMPLEID}_join' > ~/${SAMPLEID}.sam
done 

# convert sam to bam and sort
module load samtools 

for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID"
        samtools view -h -b -S ${SAMPLEID}.sam > ${SAMPLEID}.bam
        samtools sort ${SAMPLEID}.bam -o ${SAMPLEID}_sorted.bam
        samtools index ${SAMPLEID}_sorted.bam 
done

# extract unmapped sequence 
module load samtools

for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID"
        samtools view -b -f 12 -F 256 ${SAMPLEID}.bam > ${SAMPLEID}_filter.bam
done

# sort bam file
module load samtools

for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID"
        samtools sort -n ${SAMPLEID}_filter.bam -o ${SAMPLEID}_filter_sorted.bam 
done

# bam to fastq 
module load bedtools

for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID"
        bedtools bamtofastq -i ${SAMPLEID}_filter_sorted.bam -fq ${SAMPLEID}_r1.fastq -fq2 ${SAMPLEID}_r2.fastq
        mv ${SAMPLEID)_r1.fastq ~/processed_seqs
        mv ${SAMPLEID)_r2.fastq ~/processed_seqs
done

