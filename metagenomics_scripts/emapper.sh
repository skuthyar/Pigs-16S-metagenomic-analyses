# protein level
emapper.py -i ~/all.prot_appended.faa --dmnd_db data/eggnog_proteins.dmnd --data_dir ~/eggNOG/eggnog-mapper-2.1.7/data/ -o ~/all.protein.eggout --cpu 6 --matrix BLOSUM62 --seed_ortholog_evalue 1e-5 --dbtype seqdb -m diamond

# nucleotide level 
emapper.py -i ~/allsamples.drep_nucl.fna --itype CDS --translate --dmnd_db data/eggnog_proteins.dmnd --data_dir ~/eggNOG/eggnog-mapper-2.1.7/data/ -o ~/all.nucl.eggout --cpu 6 --matrix BLOSUM62 --seed_ortholog_evalue 1e-5 --dbtype seqdb -m diamond 
