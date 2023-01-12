for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID" 
        fastp -i ~/raw/${SAMPLEID}_S*_L001_R1_001.fastq.gz -I ~/raw/${SAMPLEID}_S*_L001_R2_001.fastq.gz --cut_by_quality3 -W 4 -M 20 -n 5 -c -l 50 -w 3 -o ~/clean/${SAMPLEID}_R1_001.clean.fastq.gz -O ~/clean/${SAMPLEID}_R2_001.clean.fastq.gz
done
