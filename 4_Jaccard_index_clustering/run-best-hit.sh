#!/bin/tcsh
# Run bidirectional blast analysis on Jaccard clusters
setenv PERL_MOD_DIR /usr/local/projects/ergatis/package-latest/lib/perl5


# run on full set of BSML docs
./CogBsmlLoader-mod.pl \
 --bsmlSearchList=wu-blastp.bsml.list \
 --jaccardClusters=jaccard-clusters.out \
 --outfile=all-best-hits.btab \     
 --pvalcut=1e-5 \    #1e-5 \
 --coverageCutoff=70 \
 --log=run-cog-bsml-loader.log


# do reciprocal best hit analysis
setenv WDIR /local/scratch/bma/IMG/vdb/crabtree/bing-clustering
/usr/local/projects/ergatis/package-latest/bin/best_hit -i ./all-best-hits.btab -j 0.0 > ./final-clusters.cog
