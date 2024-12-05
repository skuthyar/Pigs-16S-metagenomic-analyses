# VFDB for protein 
module load blast

for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID"
        blastp -query ~/all.prot_appended_clean.faa -db ~/VFDB/VFDB_setB_pro_edited.fas -out ~/${SAMPLEID}.vfdb.tab -evalue 1e-5 -outfmt 6 -num_threads 6 -num_alignments 5 
done

# VFDB for nucl
module load blast

for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID"
      blastn -query ~/allsamples.drep_nucl.fna -db ~/VFDB/VFDB_setB_nt.fas -out ~/${SAMPLEID}.nucl.vfdb.tab -evalue 1e-5 -outfmt 6 -num_threads 6 -num_alignments 5 
done
