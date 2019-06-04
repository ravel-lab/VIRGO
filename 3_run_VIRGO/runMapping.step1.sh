#!/bin/bash

display_usage() { 
	echo "Argument:" 
	echo -e "\nUsage:\n-r read.fastq\n-p prefix \n-d VIRGO-path\n" 
	} 
	if [  $# -le 3 ] 
	then 
		display_usage
		exit 1
	fi 
	if [[ ( $# == "--help") ||  $# == "-h" ]] 
	then 
		display_usage
		exit 0
	fi 

echo "Arguments are all good !!!"

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -r|--reads)
    READS="$2"
    shift # past argument
    shift # past value
    ;;
    -p|--prefix)
    PREFIX="$2"
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

echo "reads file  = ${READS}"
echo "sample prefix  = ${PREFIX}"
echo "VIRGO path   = ${DIR}"

dir=`echo $PWD`
output=$dir/temp_mapping
mkdir -p $output
reads=2
ref_db=$DIR/0_db
dir_annot=$DIR/1_VIRGO
ref=$ref_db/VIRGO
bowtie -p 16 -l 25 --fullref --chunkmbs 200 --best --strata -m 20 --mm $ref --suppress 2,4,5,6,7,8 $READS $output/$PREFIX.reads2ref
awk -F"\t" '{a[$2]++;}END{ for (i in a) print i"\t"a[i]}' $output/$PREFIX.reads2ref | sort -t$'\t' -k2nr,2 | awk -F"\t" -v reads=$reads '$2 >= reads' | awk -F"\t" 'NR==FNR{a[$1]=$2;next} ($1 in a) {print $0"\t"a[$1]}' $dir_annot/0.geneLength.txt - > $output/$PREFIX.out
rm $output/$PREFIX.reads2ref 
