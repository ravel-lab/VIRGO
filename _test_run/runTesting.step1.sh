#!/bin/bash

display_usage() { 
	echo "Argument:" 
	echo -e "\nUsage:\n-1 read1.fastq(.gz)\n-2 read2.fastq(.gz)\n-p prefix \n-d VIRGO-path\n" 
	} 
	if [  $# -le 4 ] 
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
    -1|--read1)
    READ1="$2"
    shift # past argument
    shift # past value
    ;;
    -2|--read2)
    READ2="$2"
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

echo "read file 1  = ${READ1}"
echo "read file 2  = ${READ2}"
echo "sample prefix    = ${PREFIX}"
echo "VIRGO path   = ${DIR}"

dir=`echo $PWD`
output=$dir/temp_testsample
mkdir -p $output
reads=2
seqtk sample -s=33 $READ1 1000000 > $output/$PREFIX.1.fq
seqtk sample -s=33 $READ2 1000000 > $output/$PREFIX.2.fq
db_dir=$DIR/0_db
dir_annot=$DIR/1_VIRGO

ref=$db_dir/VIRGO
echo "$ref"
bowtie -p 16 -l 25 --fullref --chunkmbs 200 --best --strata -m 20 --mm $ref --suppress 2,4,5,6,7,8 $output/$PREFIX.1.fq $output/$PREFIX.1.reads2ref
bowtie -p 16 -l 25 --fullref --chunkmbs 200 --best --strata -m 20 --mm $ref --suppress 2,4,5,6,7,8 $output/$PREFIX.2.fq $output/$PREFIX.2.reads2ref
cat $output/$PREFIX.1.reads2ref $output/$PREFIX.2.reads2ref | awk -F"\t" '{a[$2]++;}END{ for (i in a) print i"\t"a[i]}' - | sort -t$'\t' -k2nr,2 | awk -F"\t" -v reads=$reads '$2 >= reads' | awk -F"\t" 'NR==FNR{a[$1]=$2;next} ($1 in a) {print $0"\t"a[$1]}' $dir_annot/0.geneLength.txt - > $output/$PREFIX.out
rm $output/$PREFIX.1.reads2ref $output/$PREFIX.2.reads2ref $output/$PREFIX.1.fq $output/$PREFIX.2.fq 

