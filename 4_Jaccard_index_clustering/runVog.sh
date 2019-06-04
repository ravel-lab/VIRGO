#!/bin/bash

## link the tab and bsml list here
fasta=$1 ##the fasta file that is AA that was used originally for all-way blastp
cp /home/bma/.stuffbox/script/VOG/0_jaccard/z* .
cp /home/bma/.stuffbox/script/VOG/0_jaccard/postVog.sh .
mkdir -p logs
#./z01_jac_clustering.sh
qsub -cwd -b y -l mem_free=45G -P jravel-lab -q all.q -V -e logs/ -o logs/ ./z01_jac_clustering.sh
#./z02_run_best_hit.sh
qsub -cwd -b y -l mem_free=4G -P jravel-lab -q all.q -V -e logs/ -o logs/ ./z02_run_best_hit.sh
#./z03_make_fasta_files.sh $fasta > cluster.log

./postVog.sh
