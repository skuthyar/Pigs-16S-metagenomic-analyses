# dereplication at protein level
for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID"
        awk '/>/{sub(">","&"FILENAME"--");sub(/\.faa/,x)}1' ${SAMPLEID}_all_final_protein.faa > ${SAMPLEID}_prot_appended.faa
done    

source activate cd-hit

for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID"
        cd-hit -i ~/prodigal/${SAMPLEID}_prot_appended.faa -o ~/${SAMPLEID}_all_final_derep_protein.faa -c 1.00 -n 5 -M 80000 -d 0 -T 16
        cat ~/${SAMPLEID}_all_final_derep_protein.faa|grep "^>"|awk -F ' ' '{print $1}'|awk -F '>' '{print $2}' >~/PIG_geneID.list
        cat ~/${SAMPLEID}_all_final_derep_protein.faa > ~/all.prot_appended.faa
done

# dereplication at nucleotide level
for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID"
        awk '/>/{sub(">","&"FILENAME"--");sub(/\.faa/,x)}1' ${SAMPLEID}_all_final_nucl.fna > ${SAMPLEID}_nucl_appended.fna
done

source activate cd-hit

for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID" 
        cd-hit-est -i ~/prodigal/${SAMPLEID}_nucl_appended.fna -o ~/${SAMPLEID}_derep_nucl.fna -c 0.95 -n 10 -d 0 -M 16000 -T 8
done        

cat ~/*_derep_nucl.fna > ~/allsamples.drep_nucl.fna
