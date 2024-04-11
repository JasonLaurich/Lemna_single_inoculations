#! /bin/bash

#some useful commands
#conda init
#conda info --envs

#start qiime
conda activate qiime2-2023.2

#import data (just Churchill and Wellspring field microbiomes)
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_microbiome_fastqmanifest.csv \
  --input-format PairedEndFastqManifestPhred33 \
  --output-path /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_imported_demux.qza

#summarize imported data
qiime demux summarize \
  --i-data /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_imported_demux.qza \
  --o-visualization /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_imported_demux_view.qzv

#trim reads
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_imported_demux.qza \
  --p-cores 40 \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --o-trimmed-sequences /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_imported_trim.qza 

##merge paired-end reads
qiime vsearch merge-pairs \
  --i-demultiplexed-seqs /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_imported_trim.qza \
  --o-merged-sequences /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_imported_join.qza

##quality filter
qiime quality-filter q-score \
  --i-demux /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_imported_join.qza \
  --o-filtered-sequences /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_imported_filtered.qza \
  --o-filter-stats /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_imported_filter-stats.qza

##demultiplex summmarize to see what truncation should be
qiime demux summarize \
--i-data /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_imported_filtered.qza \
--o-visualization /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_imported_filter.qzv

#run deblur
qiime deblur denoise-16S \
--i-demultiplexed-seqs /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_imported_filtered.qza \
--p-trim-length 400 \
--p-sample-stats \
--o-table /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_deblurtable.qza \
--o-representative-sequences /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_deblurrep-seqs.qza \
--o-stats /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_deblur-stats.qza

#summarize and visualize each output file 
qiime feature-table summarize \
  --i-table /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_deblurtable.qza \
  --o-visualization /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_deblurtable.qzv

qiime feature-table tabulate-seqs \
  --i-data /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_deblurrep-seqs.qza \
  --o-visualization /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_deblurrep-seqs.qzv

qiime deblur visualize-stats \
  --i-deblur-stats /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_deblur-stats.qza \
  --o-visualization /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_deblur-stats.qzv

#assign taxonomy
#Downloaded 3-Feb-2024:
#2022.10.backbone.full-length.fna.qza
#2022.10.taxonomy.asv.nwk.qza
#2022.10.phylogeny.asv.nwk.qza

qiime greengenes2 non-v4-16s \
--i-table /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_deblurtable.qza \
--i-sequences /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_deblurrep-seqs.qza \
--i-backbone /symbiont/megan.frederickson/Lemna_microbiome/2022.10.backbone.full-length.fna.qza \
--o-mapped-table /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_gg2.qza \
--o-representatives /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_gg2.fna.qza

qiime greengenes2 taxonomy-from-table \
      --i-reference-taxonomy /symbiont/megan.frederickson/Lemna_microbiome/2022.10.taxonomy.asv.nwk.qza \
      --i-table /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_gg2.qza \
      --o-classification /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_gg2.taxonomy.qza

#filter mitochondria and chloroplasts
qiime taxa filter-table \
  --i-table /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_gg2.qza \
  --i-taxonomy /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_gg2.taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well-no-mitochondria-no-chloroplast.qza

qiime tools export \
  --input-path /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well-no-mitochondria-no-chloroplast.qza \
  --output-path /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_exported-feature-table

biom convert \
  -i /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_exported-feature-table/feature-table.biom \
  -o /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_exported-feature-table/feature_table.tsv \
  --to-tsv

qiime taxa filter-seqs \
--i-sequences /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_gg2.fna.qza \
--i-taxonomy /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_gg2.taxonomy.qza \
--p-include p__ \
--p-exclude mitochondria,chloroplast \
--o-filtered-sequences /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well-seqs-with-phyla-no-mitochondria-chloroplast.qza

qiime tools export \
 --input-path /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well-seqs-with-phyla-no-mitochondria-chloroplast.qza \
 --output-path /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_taxonomy

#make barplot
qiime taxa barplot \
--i-table /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well-no-mitochondria-no-chloroplast.qza \
--i-taxonomy /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well_gg2.taxonomy.qza \
--o-visualization /symbiont/megan.frederickson/Lemna_microbiome/church-well/taxa-bar-plots.qzv

#Single strain Sanger sequences

#Trimmed, joined, and cleaned in Geneious

#NConcatenate into single fasta file
cat /symbiont/megan.frederickson/jason-single-strain/'raw data'/fasta/* >> /symbiont/megan.frederickson/jason-single-strain/'raw data'/fasta/combined.fna

#Import to qiime2

qiime tools import \
  --input-path /symbiont/megan.frederickson/jason-single-strain/'raw data'/fasta/combined.fna \
  --output-path /symbiont/megan.frederickson/jason-single-strain/ss_imported_sequences_edited.qza \
  --type 'FeatureData[Sequence]'

#qiime tools peek /symbiont/megan.frederickson/jason-single-strain/ss_imported_sequences_edited.qza

#Merge with 16S short read seuquences
qiime feature-table merge-seqs \
    --i-data /symbiont/megan.frederickson/Lemna_microbiome/church-well/church_well-seqs-with-phyla-no-mitochondria-chloroplast.qza /symbiont/megan.frederickson/jason-single-strain/ss_imported_sequences_edited.qza \
    --o-merged-data /symbiont/megan.frederickson/Lemna_microbiome/church-well/merged-data.qza

#Export
 qiime tools export \
  --input-path /symbiont/megan.frederickson/Lemna_microbiome/church-well/merged-data.qza \
  --output-path /symbiont/megan.frederickson/Lemna_microbiome/church-well/merged-data

#align sequences
qiime alignment mafft \
  --i-sequences /symbiont/megan.frederickson/Lemna_microbiome/church-well/merged-data.qza \
  --o-alignment /symbiont/megan.frederickson/Lemna_microbiome/church-well/aligned-sequences.qza

#qiime feature-table tabulate-seqs --help

qiime feature-table tabulate-seqs \
    --i-data /symbiont/megan.frederickson/Lemna_microbiome/church-well/aligned-sequences.qza \
    --o-visualization /symbiont/megan.frederickson/Lemna_microbiome/church-well/rep-seqs.qzv
  
#make tree
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences /symbiont/megan.frederickson/Lemna_microbiome/church-well/merged-data.qza \
  --o-alignment /symbiont/megan.frederickson/Lemna_microbiome/church-well/merged_alignment.qza \
  --o-masked-alignment /symbiont/megan.frederickson/Lemna_microbiome/church-well/merged_masked_alignment.qza \
  --o-tree /symbiont/megan.frederickson/Lemna_microbiome/church-well/merged_tree_unrooted.qza \
  --o-rooted-tree /symbiont/megan.frederickson/Lemna_microbiome/church-well/merged_tree_rooted.qza

 qiime tools export \
  --input-path /symbiont/megan.frederickson/Lemna_microbiome/church-well/merged_tree_rooted.qza \
  --output-path /symbiont/megan.frederickson/Lemna_microbiome/church-well/exported-merged-tree

