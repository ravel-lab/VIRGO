#!/bin/tcsh

#setenv DEREP_FASTA /local/scratch/bma/IMG/vdb/gene_index/total_4_13_14.derep.faa
setenv DEREP_FASTA /local/scratch/bma/vdb/AA_seq/nucleotide_derep_100.fas
#setenv DEREP_FASTA /local/scratch/bma/vdb/AA_seq/protein_match_derep100.faa 

#/local/scratch/bma/IMG/vdb/crabtree/bing-clustering/makeClusterFASTAFiles.pl $DEREP_FASTA jaccard-clusters.out final-clusters.cog cluster_fasta_files
/local/scratch/bma/vdb/AA_seq/vdb/makeClusterFASTAFiles.pl $DEREP_FASTA jaccard.0.55.out final-clusters.cog cluster_fasta_files_nt
