Tutorial of VIRGO & VOG


Author: Bing Ma 
email: igs.bma@gmail.com


### please contact the author if there is any question or suggestion on the optimal practical use of the database. Thank you and enjoy!


### prerequisite Package: 
### Please have the following programs added in your path!
# blast 2.2.3+
# bowtie version 1.1.2+
# GNU Awk 3.1+
# seqtk 1.0+


### database structure explanation
# VIRGO bowtie db
# VIRGO blastn db
# VOG script


### Figure: relational database of VIRGO and VOG


######################################################################
### MODULE 1: Quick and dirty test run!###############################
######################################################################


### This module is designed to test run RAW/sequencer-delivered metagenome or metatranscriptome samples, so you will know whether sample are what's expected


## step 0: set up your working directory
mkdir -p /your/working/directory
cd /your/working/directory


## step 1
## input the paired end reads to -1 and -2 for a sample, either in fastq format or a zip fastq.gz format 


Argument:

Usage:
-1 read1.fastq(.gz)
-2 read2.fastq(.gz)
-p prefix 

runTesting.step1.sh -1 /path/to/read1.fastq(.gz) -2 /path/to/read2.fastq(.gz) -p prefix samplePrefix

# for example: 
/path/to/runTesting.step1.sh -1 /local/sequence/illumina/data/2016/160822_K00134_0079_AHF33JBBXX/Unaligned/Project_XXXX/Sample_IL100073457/IL100073457_S169_L008_R1_001.fastq.gz -2 /local/sequence/illumina/data/2016/160822_K00134_0079_AHF33JBBXX/Unaligned/Project_XXXX/Sample_IL100073457/IL100073457_S169_L008_R2_001.fastq.gz -p sampleName


## step 2
runTesting.step2.sh /path/to/output/of/step1/
# for example:
/local/projects-t2/M8910/VOG_demo/2_run_VIRGO/runTesting.step2.sh /your/working/path/


## step 3 
# output includes:
# 1. summary.Abundance.txt: the number of reads per species
# 2. summary.Count.txt: the number of genes per species
# 3. summary.Percentage.txt: normalized abundance in percentage per species (sum up to 100)
# 4. summary.geneRichness.txt: number of genes per sample
# 5. annotation.txt: EggNOG annotation file per sample
# 6. Sample.out: the list of non-redundant genes with gene length


######################################################################
### MODULE 2: Mapping to VIRGO non-redundant gene database ###########
######################################################################


### this is module is to map reads from a metagenome or metatranscriptome sample to VIRGO and present the results
### !!! please note VIRGO does NOT support paired-end mapping, please merge the paired end reads alone or together with single-end reads into one fastq file.

## step 0. Prepare your working directory
mkdir -p /your/working/directory
cd /your/working/directory

## step 1
## input the read file in fastq format to -r option for a sample
## either quality assessed or raw fastq file would be fine

Argument:

Usage:
-r reads.fastq
-p prefix 


# for example: 
runMapping.step1.sh -r sample1.fastq -p samplePrefix


## step 2: run reads mapping
runMapping.step2.sh /path/to/output/of/step1/
# for example:
runMapping.step2.sh /your/working/path/temp_mapping/


## 
# output includes:
# 1. summary.Abundance.txt: the number of reads per species
# 2. summary.Count.txt: the number of genes per species
# 3. summary.Percentage.txt: normalized abundance in percentage per species (sum up to 100)
# 4. summary.geneRichness.txt: number of genes per sample
# 5. summary.NR.abundance.txt: the number of reads per non-redundant gene
# 6. gene.lst.txt: the list of non-redundant genes with gene length
# 7. EggNOG.annotation.txt: EggNOG annotation file per sample
# 8. EC.annotation.txt: list of non-redundant genes with EC number
# 9. GC.txt: list of non-redundant genes with gene count category (HGC: high gene count, LGC: low gene count)
# 10. geneProduct.txt: list of non-redundant genes with gene product annotation
# 11. Kegg.module.annotation.txt: list of non-redundant genes with KEGG module annotation, including module ID and annotation
# 12. Kegg.ortholog.annotation.txt: list of non-redundant genes with KEGG ortholog (KO) annotation, including ortholog ID and annotation
# 13. Kegg.pathway.annotation.txt: list of non-redundant genes with KEGG pathway annotation, including pathway ID and annotation
# 14. proteinFamily.annotation.txt: list of non-redundant genes with protein family annotation, from the database of CDD, GO, Gene3D, Hamap, Interpro, MobiDBLite, PIRSF, PRINTS, Pfam, ProDom, ProSitePatterns, ProSiteProfiles, SFLD, SMART, SUPERFAMILY, TIGRFAM
# 15. rxn.annotation.txt: list of non-redundant genes with KEGG reaction


######################################################################
### MODULE 3: BLAST to VIRGO database ################################
######################################################################

## VIRGO was also made blastable in both nucleotide and amino acid. Just run your normal blastn or blastp command
## # blast 2.2.3+ compatible 


######################################################################  
### MODULE 4. Use jaccard clustering to generate protein cluster
######################################################################

# step 0
Run all-vs-all BLASTP 

# step 1
#Run Jaccard clustering using modified clusterBsmlPairwiseAlignments script:

./run-jaccard.sh

Script/program:
  clusterBsmlPairwiseAlignments-mod.pl

Inputs: 
  o wu-blastp.bsml.list (list of BLASTP BSML files)
  o location of executable that generates clusters from polypeptide pairs
  o various parameters (linkscore, percent_identity, percent_coverage, p_value)

Outputs:
  o jaccard-clusters.out(.gz)

clusterBsmlPairwiseAlignments-mod.pl:
  o reads and filters all of the BLASTP data and generates a list of polypeptide pairs
  o writes the polypeptide pairs to a file (20808.pairs.gz)
  o invokes /usr/local/projects/ergatis/package-latest/bin/cluster to convert the pairs
    into clusters

Note that all of clustering/filtering parameters are applied during the pair generation
step, not the pair clustering step.

Example: 
clusterBsmlPairwiseAlignments-mod.pl \
 --bsmlSearchList=wu-blastp.bsml.list \
 --log=test-jaccard.log \
 --linkscore=0.6 \
 --percent_identity=80 \
 --percent_coverage=60 \
 --p_value=1e-5 \
 --outfile=jaccard-clusters.out \
 --cluster_path=/path/to/cluster/files/

# step 3. Run bidirectional best hit analysis on Jaccard clusters.

./run-best-hit.sh (includes steps 3a and 3b, shown below)

3a. Convert BLASTP BSML to a single btab file of best hits.

Script/program:
  o CogBsmlLoader-mod.pl

Inputs:
  o wu-blastp.bsml.list (list of BLASTP BSML files)
  o jaccard-clusters.out from Jaccard clustering step
  o parameters (minimum p-value and minimum % coverage)

Outputs:
  o summary btab file (/path/to/btab)
    (which includes the Jaccard cluster info. in the last column)

3b. Merge all clusters with reciprocal best hits.

Script/program:
  o best_hit.pl

Inputs:
  o all-best-hits.btab from step 3a
  o Minimum Jaccard coefficient (for additional cluster pruning if set to > 0)

Outputs:
  o a final cluster file (final-clusters.cog.gz)

4. Generate FASTA files.

Script/program:
  o makeClusterFASTAFiles.pl

Inputs:
  o The original polypeptide file (e.g., total.faa)
  o jaccard-clusters.out from step 2
  o final-clusters.cog from step 3b
  o The path to an output directory

Outputs:
  o A FASTA file for each cluster of size > 1
  o A singletons.fa file that contains all the clusters of size 1
  o A cluster size histogram (on stdout - see make-fasta-files-2.log.gz)



