#!/bin/bash


for X in cluster_fasta_files_jacs/*.fa; do number=`grep "^>" $X | wc -l`; name=`basename $X .fa`; echo -e "$name\t$number"; done > family.count

for X in cluster_fasta_files_jacs/*.fa; do name=`basename $X .fa`; echo $name; grep "^>" $X; done > tmp

awk -F"\t" '{if ($1 !~ /^>/) clusterID=$1; if ($1 ~ /^>/) print clusterID"\t"$0;}' tmp | sed 's/>//g' > cluster_gene.tbl

cp family.count multi.count 

sed '/^singletons/d' multi.count > tmp0
echo -e "singletons\t1" >> tmp0
mv tmp0 multi.count
##change singleton count to 1 in multi.count

awk -F"\t" 'NR==FNR{a[$1]=$2;next} {if ($1 in a) print $0"\t"a[$1]; else print $0"\tNA"}' multi.count cluster_gene.tbl > tmp

awk -F"\t" 'NR==FNR{a[$1]=$2;next} {if ($2 in a) print $0"\t"a[$2]; else print $0"\tNA"}' /local/projects-t4/bma/CHARM/vdb/AA_seq/species_name/total.species.vers2 tmp > tmp1

awk -F"\t" 'NR==FNR{a[$7]=$2;next} {if ($2 in a) print $0"\t"a[$2]; else print $0"\tNA"}' /local/scratch/bma/VOG/5_nonRedundant/clst.MG.95.ID tmp1 > tmp2

awk -F"\t" 'NR==FNR{a[$1]=$2;next} {if ($2 in a) print $0"\t"a[$2]; else print $0"\tNA"}'  /local/projects-t4/bma/CHARM/vdb/gene_index/gene_annotation/total_gene.anno tmp2 > tmp3

awk -F"\t" 'NR==FNR{a[$2]=$5;next} {if ($2 in a) print $0"\t"a[$2]; else print $0"\tNA"}' /local/scratch/bma/VOG/2_annotation/eggnog/aa.annotation tmp3 > tmp4

awk -F"\t" 'NR==FNR{a[$2]=$6;next} {if ($2 in a) print $0"\t"a[$2]; else print $0"\tNA"}' /local/scratch/bma/VOG/2_annotation/eggnog/aa.annotation tmp4 > tmp5

echo -e "family\tgene\tsize\tspecies\tGC\tannotation\teggNogg\tCOGannot" > family.short.info.txt
cat tmp5 >> family.short.info.txt

awk -F"\t" 'NR==FNR{a[$1]=$0;next} {if ($2 in a) print $0"\t"a[$2]; else print $0"\tNA"}' /local/scratch/bma/VOG/2_annotation/all.merge.tbl.txt tmp5 > tmp6

#awk -F"\t" '{print $2}' tmp5 > tmp6
#for X in `cat tmp6`; do grep -w "^$X" /local/scratch/bma/VOG/2_annotation/all.merge.tbl.txt; done > tmp7

echo -e "family\tgene\tsize\tspecies\tGC\tannotation\teggNogg\tCOGannot\tPC\tCDD\tGO\tGene3D\tHamap\tInterpro\tMobiDBLite\tPIRSF\tPRINTS\tPfam\tProDomProSitePatterns\tProSiteProfiles\tSFLD\tSMART\tSUPERFAMILY\tTIGRFAM" > family.long.info.txt

cat tmp6 >> family.long.info.txt


rm tmp*


