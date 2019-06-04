#!/bin/tcsh

setenv PERL_MOD_DIR /usr/local/projects/ergatis/package-latest/lib/perl5
set pi=80
set ls=0.6

#Run Jaccard clustering using modified clusterBsmlPairwiseAlignments script:
clusterBsmlPairwiseAlignments-mod.pl \
 --bsmlSearchList=wu-blastp.bsml.list \
 --log=test-jaccard.log \
 --linkscore=$ls \
 --percent_identity=$pi \
 --percent_coverage=60 \
 --p_value=1e-5 \
 --outfile=jaccard-clusters.out \
 --cluster_path=/usr/local/projects/ergatis/package-latest/bin/cluster

