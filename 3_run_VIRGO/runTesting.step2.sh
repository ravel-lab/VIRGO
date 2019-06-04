#!/bin/bash

display_usage() {
        echo "Argument:" 
        echo -e "\nUsage:\n-p working-path\n-d VIRGO-path\n" 
        }
        if [  $# -le 2 ]
        then
                display_usage
                exit 1
        fi
        if [[ ( $# == "--help") ||  $# == "-h" ]]
        then
                display_usage
                exit 0
        fi
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -p|--working-path)
    WPATH="$2"
    shift # past argument
    shift # past value
    ;;
    -d|--VIRGO-path)
    DIR="$2"
    shift # past argument
    shift # past value
    ;;
    --default)
    DEFAULT=YES
    shift # past argument
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" 

echo "WORKING PATH   = ${WPATH}"
echo "VIRGO path   = ${DIR}"

cd ${WPATH}/
dir=${DIR}/1_VIRGO
echo -e "sample\tgene count" > summary.geneRichness.txt
for X in *.out; do name=`basename $X .out`; count=`cat $X | wc -l`; echo -e "$name\t$count"; done >> summary.geneRichness.txt
for X in *.out; do name=`basename $X .out`; awk -F"\t" 'NR==FNR{a[$2]=$3;next} {if ($1 in a) print $0"\t"a[$1]}' $dir/1.taxon.tbl.txt $X | awk -F"\t" -v name="${name}" '{a[$4]++;b[$4]+=$2}END{ for (i in a) print name"\t"i"\t"a[i]"\t"b[i]}' | sort -t$'\t' -k3nr,3 > $name.tmp; done 
for X in *.tmp; do name=`basename $X .tmp`; awk -v OFS="\t" 'NR==FNR{a = a + $4;next} {perc = ($4/a)*100;print $0"\t"perc }' $X $X > $name.tmp2; rm $X; done
cat *.tmp2 > tmp
awk 'BEGIN {FS=OFS="\t"} NR>1 {total[$1$2]=$3;id[$1]++;name[$2]++} END {n=asorti(name,copy1);printf "PC"; for (i=1;i<=n;i++) {printf "%s%s",FS,copy1[i]}; print ""; m=asorti(id,copy2); for (j=1;j<=m;j++) {printf "%s",copy2[j]; for (k=1;k<=n;k++) {printf "%s%s",FS,total[copy2[j]copy1[k]]}; print ""}}' tmp | awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = 0 }; 1' > summary.Count.txt 
awk 'BEGIN {FS=OFS="\t"} NR>1 {total[$1$2]=$4;id[$1]++;name[$2]++} END {n=asorti(name,copy1);printf "PC"; for (i=1;i<=n;i++) {printf "%s%s",FS,copy1[i]}; print ""; m=asorti(id,copy2); for (j=1;j<=m;j++) {printf "%s",copy2[j]; for (k=1;k<=n;k++) {printf "%s%s",FS,total[copy2[j]copy1[k]]}; print ""}}' tmp | awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = 0 }; 1' > summary.Abundance.txt
awk 'BEGIN {FS=OFS="\t"} NR>1 {total[$1$2]=$5;id[$1]++;name[$2]++} END {n=asorti(name,copy1);printf "PC"; for (i=1;i<=n;i++) {printf "%s%s",FS,copy1[i]}; print ""; m=asorti(id,copy2); for (j=1;j<=m;j++) {printf "%s",copy2[j]; for (k=1;k<=n;k++) {printf "%s%s",FS,total[copy2[j]copy1[k]]}; print ""}}' tmp | awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = 0 }; 1' > summary.Percentage.txt
rm tmp *.tmp2
for X in *.out; do name=`basename $X .out`; echo -e "Sample\tAnnotation\tnumGene\tnumReads" > $name.annotation.txt; awk -F"\t" 'NR==FNR{a[$2]=$6;next} {if ($1 in a) print $0"\t"a[$1]}' $dir/3.eggnog.NOG.txt $X | awk -F"\t" -v name="${name}" '{a[$4]++;b[$4]+=$2}END{ for (i in a) print name"\t"i"\t"a[i]"\t"b[i]}' | sort -t$'\t' -k3nr,3 >> $name.annotation.txt; done
