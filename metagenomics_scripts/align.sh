# Align protein sequence of gene catalog to the Uniprot TrEMBL database
conda activate diamondenv

for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID"
        diamond blastp -q ~/derep_protein/${SAMPLEID}_all_final_derep_protein.faa -d uniprot_trembl.dmnd -p 8 --sensitive -o ${SAMPLEID}_all_final_derep_protein.faa.default.diamond2uniprot_trembl    
done

# Taxonomic classification based on the LCA algorithms 
conda activate basta_py3

for SAMPLEID in $(cat ~/samples.txt);

do
        echo "On sample: $SAMPLEID"
        basta sequence -l 25 -i 80 -e 0.00001 -m 3 -b 1 -p 60 ~/${SAMPLEID}_all_final_derep_protein.faa.default.diamond2uniprot_trembl ~/${SAMPLEID}_all_final_derep_protein.faa.diamond2uniprot_trembl.dmnd.lca.out ~/taxonomy/prot
done
