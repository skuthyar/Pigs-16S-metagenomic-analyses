conda activate prodigalenv

for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID"
        prodigal -i ~/megahit/${SAMPLEID}_all_final_contigs.fasta -a ~/prodigal/${SAMPLEID}_all_final_protein.faa -d ~/prodigal/${SAMPLEID}_all_final_nucl.fna -o ~/prodigal/${SAMPLEID}_all_final.gff -f gff -p meta
done
