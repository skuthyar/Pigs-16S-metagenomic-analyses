# CARD for proteins
source activate rgi

for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID"
        rgi main -i ~/all.prot_appended_clean.faa -o ~/${SAMPLEID}_protein_faa.card -t protein -a DIAMOND -n 6 --local --include_loose --clean --debug
done

# CARD for nucleotide file
source activate rgi

for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID"
        rgi main -i ~/allsamples.drep_nucl.fna -o ~/${SAMPLEID}_nucl.card -a DIAMOND -n 6 --local --include_loose --clean --debug
done
