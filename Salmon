#build Salmon index 
source activate salmon

salmon index -t ~/allsamples.drep_nucl.fna -i ~/allsamples.salmon_nucl_index

#Salmon quantification
source activate salmon

for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID"
        salmon quant -i ~/allsamples.salmon_nucl_index -l A -1 ~/processed_seqs/${SAMPLEID}_final_R1.fastq.gz -2 ~/processed_seqs/${SAMPLEID}_final_R2.fastq.gz -p 8 --validateMappings -o ~/${SAMPLEID}_quant_nucl
done
