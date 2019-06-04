#!/bin/tcsh

set fasta=$1
setenv DEREP_FASTA $fasta 

mkdir cluster_fasta_files_jacs
makeClusterFASTAFiles.pl $DEREP_FASTA jaccard-clusters.out final-clusters.cog cluster_fasta_files_jacs 
