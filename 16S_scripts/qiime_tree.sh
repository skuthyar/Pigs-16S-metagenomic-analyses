source ~/.bashrc
conda activate qiime2-2021.2

#import fasta file and convert to qiime2 object (.qza)
#mb_filtered_seqs.fna is a fasta file with all ASVs and ids
qiime tools import \
  --input-path ASVs_decontam.fa \
  --output-path mb_filtered_seqs.qza \
  --type 'FeatureData[Sequence]'

#run multiple alignment
qiime alignment mafft \
  --i-sequences mb_filtered_seqs.qza \
  --o-alignment aligned-mb-seqs.qza

#mask alignment to reduce ambiguity
qiime alignment mask \
  --i-alignment aligned-mb-seqs.qza \
  --o-masked-alignment masked-aligned-mb-seqs.qza

qiime phylogeny fasttree \
  --i-alignment masked-aligned-mb-seqs.qza \
  --o-tree fasttree-tree.qza

qiime phylogeny midpoint-root \
  --i-tree fasttree-tree.qza \
  --o-rooted-tree fasttree-tree-rooted.qza

#export rooted tree
qiime tools export \
  --input-path fasttree-tree-rooted.qza \
  --output-path mb-rooted-fasttree_filtered.tre
  
