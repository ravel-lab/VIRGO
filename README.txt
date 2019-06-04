DISCLOSURE: VIRGO/VOG is free for Academia as open source. It is not free for industry. 


Tutorial of VIRGO (human vaginal non-redundant gene catalog) & VOG (vaginal orthologous groups)


Author: Bing Ma 
email: igs.bma@gmail.com


### please contact the author if there is any question or suggestions. Thank you! Enjoy!


### prerequisite Package ###
### Please have the following programs added to your path!
# blast 2.2.3+
# bowtie version 1.1.2+
# GNU Awk 3.1+
# seqtk 1.0+


### database structure ###
## _test_run contains the test dataset and test results
## 0_db contains the searchable database and sequences
# VIRGO bowtie version 1 database
# VIRGO bowtie version 2 database
# VIRGO BLASTN database
# VIRGO BLASTP database
## 1_VIRGO contains the annotation for all genes
# 0.geneLength.txt contains nucleotide basepair length      
# 1.taxon.tbl.txt contains taxonomy and gene length
# 2.A.proteinFamily.txt contains annotation from 17 proteins databases such as TIGRFAM, Interpro, SUPERFAMILY, CDD, Pfam, one source for a gene per line
# 2.B.proteinFamily.txt contains annotation from 17 proteins databases such as TIGRFAM, Interpro, SUPERFAMILY, CDD, Pfam, tabulated for one gene per line
# 3.eggnog.NOG.txt contains eggNOG annotation
# 4.GC.txt contains high or low gene count information, indicating a gene that is >95% of the times seen in a high or low gene count community
# 5.EC.txt contains EC number (Enzyme Commission number)
# 6.product.txt contains gene product annotation 
# 7.rxn.txt KEGG reaction annotatio
# 8.A.kegg.ortholog.txt  KEGG ortholog annotation
# 8.B.kegg.module.txt KEGG module annotation
# 8.C.kegg.pathway.txt KEGG pathway annotation
## 2_fungal_phage  contains the annotation for 10,908 fungal and 15,965 phage genes
## 3_run_VIRGO contains scripts (used below) to run metagenomes or metatranscriptome samples
## 4_Jaccard_index_clustering contains script to perform jaccard clustering to generate protein families
## 5_VOG contains the VIRGO to VOG and annotations for each of the protein family
# 0.VOG.VIRGO.tbl.txt  contains the VIRGO ID to VOG ID
# 1.VOG.size.txt  # the size of the VOG protein family
# 2.VOG.GC.txt  # gene count 
# 3.VOG.funcat.txt # functional category based on eggNOG
# 4.VOG.taxonomy.txt  # consensus taxonomy 
# 5.VOG.alignmentScore.txt # the alignment score indicate quality of the protein family alignment 

######################################################################
### MODULE 0: Test run VIRGO #########################################
######################################################################

### It is to get you a quick idea of the output, how to run, etc. Feel free to mess with it!

In directory _test_run contains the test dataset, 

> cd _test_run

## specify the path that VIRGO was downloaded to, it should contain the database structure such as /path/to/VIRGO/0_db, and /path/to/VIRGO/1_VIRGO
> ./runTesting.step1.sh -1 sub1.fq -2 sub2.fq -p test -d /path/to/VIRGO/
> ./runTesting.step2.sh -p temp_testsample -d /path/to/VIRGO/
# Done test run for the quick run! You should see outputs of all files as below! 
> ls temp_testsample/
summary.Abundance.txt  summary.Percentage.txt    test.1.fq         test.2.fq            test.out
summary.Count.txt      summary.geneRichness.txt  test.1.reads2ref  test.annotation.txt

>./runMapping.step1.sh -r sub1.fq -p sub1 -d /path/to/VIRGO/
>./runMapping.step1.sh -r sub2.fq -p sub2 -d /path/to/VIRGO/
>./runMapping.step2.sh -p temp_mapping/ -d /path/to/VIRGO/
>ls temp_mapping/
sub1.EC.annotation.txt             sub1.reads2ref                     sub2.out
sub1.EggNOG.annotation.txt         sub1.rxn.annotation.txt            sub2.proteinFamily.annotation.txt
sub1.GC.txt                        sub2.EC.annotation.txt             sub2.rxn.annotation.txt
sub1.gene.lst.txt                  sub2.EggNOG.annotation.txt         summary.Abundance.txt
sub1.geneProduct.txt               sub2.GC.txt                        summary.Count.txt
sub1.kegg.module.annotation.txt    sub2.gene.lst.txt                  summary.NR.abundance.txt
sub1.kegg.ortholog.annotation.txt  sub2.geneProduct.txt               summary.Percentage.txt
sub1.kegg.pathway.annotation.txt   sub2.kegg.module.annotation.txt    summary.geneRichness.txt
sub1.out                           sub2.kegg.ortholog.annotation.txt
sub1.proteinFamily.annotation.txt  sub2.kegg.pathway.annotation.txt


######################################################################
### MODULE 1: To explore a single gene of interest ###################
######################################################################

### It is to quickly pull annotations for your favorite gene, just using a single VIRGO geneID!!

## Step 1. For a gene of interest, just perform a BLASTN or BLASTP or your favorite similarity search against our blastable/searchable database to find similarity hit

# For example, you can use BLAST. Below is a pseudo-code to run blast
> blastn -query nucleotide.fasta -db 0_db/VIRGO -outfmt 6 -out yourGene.nucleotide.txt 
> blastp -query protein.fasta -db 0_db/VIRGO -outfmt 6 -out yourGene.protein.txt 

## Step 2. After you found a VIRGO gene hit to your favorite gene (for example V1000057 is the top hit with high sequence similarity), by "grep" annotation files of VIRGO, you will get lots of information in a second! For example, for the hit gene V1000057 you just found, you can pull out its length, taxonomy, GO, TIGRFAM, Interpro, SUPERFAMILY, CDD, Pfam, eggNOG, EC number, gene product, KEGG ortholog, KEGG module, KEGG pathway etc.

>grep -w "V1000057" 1_VIRGO/*

0.geneLength.txt:V1000057	990
1.taxon.tbl.txt:Cluster_219398	V1000057	Lactobacillus_iners	990
2.A.proteinFamily.txt:V1000057	GO	; Phosphoribosyltransferase domain; Phosphoribosyltransferase-like; Ribose-phosphate diphosphokinase; Ribose-phosphate pyrophosphokinase, N-terminal domain
2.A.proteinFamily.txt:V1000057	TIGRFAM	;  ribose-phosphate diphosphokinase
2.A.proteinFamily.txt:V1000057	Interpro	; Phosphoribosyltransferase domain; Phosphoribosyltransferase-like; Ribose-phosphate diphosphokinase; Ribose-phosphate pyrophosphokinase, N-terminal domain
2.A.proteinFamily.txt:V1000057	SUPERFAMILY	; PRTase-like; PRTase-like
2.A.proteinFamily.txt:V1000057	CDD	; PRTases_typeI
2.A.proteinFamily.txt:V1000057	Pfam	; N-terminal domain of ribose phosphate pyrophosphokinase; Phosphoribosyl synthetase-associated domain
2.B.proteinFamily.txt:V1000057	; PRTases_typeI	; Phosphoribosyltransferase domain; Phosphoribosyltransferase-like; Ribose-phosphate diphosphokinase; Ribose-phosphate pyrophosphokinase, N-terminal domain			; Phosphoribosyltransferase domain; Phosphoribosyltransferase-like; Ribose-phosphate diphosphokinase; Ribose-phosphate pyrophosphokinase, N-terminal domain			; N-terminal domain of ribose phosphate pyrophosphokinase; Phosphoribosyl synthetase-associated domain						; PRTase-like; PRTase-like	;  ribose-phosphate diphosphokinase
3.eggnog.NOG.txt:Cluster_219398	V1000057	PRS	map00030,map00230,map01100,map01110,map01120,map01230	F	Phosphoribosyl pyrophosphate synthase	COG0462
5.EC.txt:Cluster_219398	V1000057	ec:2.7.6.1		 ribose-phosphate diphosphokinase; ribose-phosphate pyrophosphokinase; PRPP synthetase; phosphoribosylpyrophosphate synthetase; PPRibP synthetase; PP-ribose P synthetase; 5-phosphoribosyl-1-pyrophosphate synthetase; 5-phosphoribose pyrophosphorylase; 5-phosphoribosyl-alpha-1-pyrophosphate synthetase; phosphoribosyl-diphosphate synthetase; phosphoribosylpyrophosphate synthase; pyrophosphoribosylphosphate synthetase; ribophosphate pyrophosphokinase; ribose-5-phosphate pyrophosphokinase
6.product.txt:Cluster_219398	V1000057	Ribose-phosphate pyrophosphokinase
8.A.kegg.ortholog.txt:V1000057	K00948	ljf:FI9785_163	dvvi:GSVIVP00018977001  GSVIVT00018977001; assembled CDS; K00948 ribose-phosphate pyrophosphokinase [EC:2.7.6.1]
8.B.kegg.module.txt:V1000057	M00005	PRPP biosynthesis, ribose 5P => PRPP
8.C.kegg.pathway.txt:V1000057	map00230	Purine metabolism

## Step 3. You can further pull out the nucleotide and amino acid sequences

> sed -n '/V1000057/{n;p;}' 0_db/NT.fasta 
ATGAACAAAGAAGAATTACATGAATTATTACATCCAATGAAGTTAATTGGGTTGGGTGGCAACCAGACATTGGCTAAGCGTATTGCAAAAGCTTTGAATAAAGAGTTACTTGAGACAGCTGTTAAGCATTTTAGTGATGGTGAAATTCAAGTTAATATTACAGAAACTGTTAGAGGCTGTGATGTTTATGTTATTCAATCTATTCAAGATCCAGTCAATGAGAATTTCATGGAACTAATGATTGTTTTGGATGCTTTAAATAGAGCTTCAGCACACAGTGTGACAGTAGTTGTTCCTTATATGGCTTATTCTCGTCAAGATACTAAAAATAGATCTCGTGAGCCTATTACTGCTAAATTATTGGCTAATTTGCTACAATTAACTAGAATTGATGATGTAGTTGCTTTGGATTTACATGCATCTCAAATTCAAGGATTTTACAATATACCTGTAGATCATCTTCATGCTTTACCAATTCTTGCTCAGTATTTTTTAGATAATGGTATTGTAGATAGAGATTCTGATTCCGTAGTGGTAGTTTCTCCAGATCATTCAGGTGCTAAATTGGCAAGGAACTTTGGTTCATATTTCAATGCACCAATAGCTATTGTTGATCAACGAGGTGCTCGTTATGATACTGAAGTGCATGATATGATAGGTGATGTTAAGGATAAAACTGTAATTATTGTAGATGATTTGATTGATACTGGTTCTAGAATATCGTCATCTACCAAATCAGTATTGGCAGCTGGTGCGAAAAAAGTGTACGTGGCTGCAACTCATGCCTTGTTATCTCAAAATGCTATTGATGTTTTAAATCAATTAAATTTGGAACAAATAGTTGTAACTGATACTATAGAGCATAAAAAATATCCTGATAGGATAGTACGGTTATCTGTAGATAGATTGCTTGCCAAAGGTATAGATTATATCTACAATGATAAATCTATTCATCAAATTTTTGATGAACAAAATAGAATACACCATTAA

> sed -n '/V1000057/{n;p;}' 0_db/AA.fasta 
MNKEELHELLHPMKLIGLGGNQTLAKRIAKALNKELLETAVKHFSDGEIQVNITETVRGCDVYVIQSIQDPVNENFMELMIVLDALNRASAHSVTVVVPYMAYSRQDTKNRSREPITAKLLANLLQLTRIDDVVALDLHASQIQGFYNIPVDHLHALPILAQYFLDNGIVDRDSDSVVVVSPDHSGAKLARNFGSYFNAPIAIVDQRGARYDTEVHDMIGDVKDKTVIIVDDLIDTGSRISSSTKSVLAAGAKKVYVAATHALLSQNAIDVLNQLNLEQIVVTDTIEHKKYPDRIVRLSVDRLLAKGIDYIYNDKSIHQIFDEQNRIHH


######################################################################
### MODULE 2: Explore a list of genes of interest#####################
######################################################################

### It is to quickly pull out information for a list of genes of interest, using a list of VIRGO ID

## Step 1. Generate a list of genes of VIRGO ID. For example, a list file with VIRGO geneID

> cat lst
V1891288
V1891259
V1891167
V1891137
V1891104

## Step 2. Pull annotations for this list of genes, such as taxonomy, functional annotation, eggNOG, KEGG, etc.
> awk -F"\t" 'NR==FNR{a[$2]=$3;next} ($1 in a) {print $0"\t"a[$1]}' 1_VIRGO/1.taxon.tbl.txt lst
V1891288	Lactobacillus_crispatus
V1891259	Lactobacillus_crispatus
V1891167	Lactobacillus_crispatus
V1891137	Lactobacillus_crispatus
V1891104	Lactobacillus_crispatus


> awk -F"\t" 'NR==FNR{a[$1]=$0;next} ($1 in a) {print a[$1]}' 1_VIRGO/2.B.proteinFamily.txt lst
V1891288		; Amino acid/polyamine transporter I			; Amino acid/polyamine transporter I				; Amino acid permease		
V1891259		; Diphosphomevalonate decarboxylase; Diphosphomevalonate/phosphomevalonate decarboxylase; GHMP kinase N-terminal domain; GHMP kinase, C-terminal domain; Ribosomal protein S5 domain 2-type fold; Ribosomal protein S5 domain 2-type fold, subgroup			; Diphosphomevalonate decarboxylase; Diphosphomevalonate/phosphomevalonate decarboxylase; GHMP kinase N-terminal domain; GHMP kinase, C-terminal domain; Ribosomal protein S5 domain 2-type fold; Ribosomal protein S5 domain 2-type fold, subgroup		; (Full/Desc.) diphosphomevalonate decarboxylase		; GHMP kinases C terminal; GHMP kinases N terminal domain						; Ribosomal protein S5 domain 2-like; GHMP Kinase, C-terminal domain	;  diphosphomevalonate decarboxylase
V1891137		; AAA domain; Exopolysaccharide synthesis protein; P-loop containing nucleoside triphosphate hydrolase; Polymerase/histidinol phosphatase-like	; P-loop containing nucleotide triphosphate hydrolases		; AAA domain; Exopolysaccharide synthesis protein; P-loop containing nucleoside triphosphate hydrolase; Polymerase/histidinol phosphatase-like			; AAA domain						; P-loop containing nucleoside triphosphate hydrolases; PHP domain-like	;  capsular exopolysaccharide family
V1891104		; Bacteriophage/Gene transfer agent portal protein			; Bacteriophage/Gene transfer agent portal protein				; Phage portal protein		

> awk -F"\t" 'NR==FNR{a[$2]=$0;next} ($1 in a) {print a[$1]}'  1_VIRGO/3.eggnog.NOG.txt lst
Cluster_519635	V1891288	YBEC		E	amino acid	COG0531
Cluster_213474	V1891259	MVAD	map00900,map01100,map01110	I	diphosphomevalonate decarboxylase	COG3407
Cluster_94013	V1891137	YWQD		M	Capsular exopolysaccharide family	COG0489
Cluster_686323	V1891104			S	Phage Portal Protein	COG4695

# Step 3. Pull out gene count information. For example, 3 out of the 5 genes from the list were most likely seen in low gene count (LGC) vaginal communities. It makes sense right? L. crispatus are often seen dominating a vaginal community with low community complexity. The other 2 genes can be seen in LGC or HGC vaginal communities. 
> awk -F"\t" 'NR==FNR{a[$2]=$0;next} ($1 in a) {print a[$1]}'  ../1_VIRGO/4.GC.txt lst
Cluster_519635	V1891288	LGC
Cluster_54462	V1891167	LGC
Cluster_686323	V1891104	LGC


######################################################################
### MODULE 3: Quick and dirty test run!###############################
######################################################################

### This module is designed to test run RAW/sequencer-delivered metagenome or metatranscriptome samples, so you will know whether samples are what's expected for QC purpose. 

## Step 1. feed the paired end reads to -1 and -2 for a sample, either in fastq format or a zip fastq.gz format 

Argument:

Usage:
-1 read1.fastq(.gz)
-2 read2.fastq(.gz)
-p prefix 
-d path to VIRGO

./runTesting.step1.sh -1 /path/to/read1.fastq(.gz) -2 /path/to/read2.fastq(.gz) -p prefix -d /full/path/to/VIRGO/

# for example. The output should be a sampleName.out txt file
/path/to/runTesting.step1.sh -1 Read1.fastq.gz -2 Read2.fastq.gz -p sampleName -d /virgo/path/

# The output should be a sampleName.out txt file. The 1st column is the VIRGO gene ID, 2nd column is the # reads mapped onto the VIRGO database, 3rd column is the gene length. File was sorted by the # reads mapped. 
> head sampleName.out 
V1593031	2425	3663
V1607456	566	390
V1420626	348	1017
V1684181	327	1442
V1316648	324	480
V1541216	274	1386
V1064785	269	303
V1872750	246	1179
V1885045	227	1254
V1573770	226	291

## Step 2: repeat for all your sample sets. Below is a pseudo-code
> for sample in *.fq; do ./runTesting.step1.sh -1 $sample.1.fq -2 $sample.2.fq -p $sample -d /virgo/path/; done

## Step 3: summarize the stats of multiple samples
./runTesting.step2.sh -p /path/to/output/of/step1/ -d /path/to/VIRGO/

# for example:
./runTesting.step2.sh -p temp_testsample/ -d /path/to/VIRGO/
 
# output includes:
# 1. summary.Abundance.txt: the number of reads per species
# 2. summary.Count.txt: the number of genes per species
# 3. summary.Percentage.txt: normalized abundance in percentage per species (sum up to 100)
# 4. summary.geneRichness.txt: number of genes per sample
# 5. annotation.txt: EggNOG annotation file per sample
# 6. Sample.out: the list of non-redundant genes with gene length


######################################################################
### MODULE 4: Mapping to VIRGO non-redundant gene database ###########
######################################################################

###  FULL run module. This is module is to map reads from a metagenome or metatranscriptome sample to VIRGO and present the results
### !!! please note VIRGO does NOT yet support paired-end mapping, please merge the paired end reads alone or together with single-end reads into one fastq file.

## Step 1:  feed the fastq read file of a sample to -r. Please note either quality assessed or raw fastq file would be fine

Argument:

Usage:
-r reads.fastq
-p prefix 
-d path to VIRGO

# for example: 
./runMapping.step1.sh -r sample1.fastq -p samplePrefix -d /full/path/to/VIRGO/

# The output should be a sampleName.out txt file. The 1st column is the VIRGO gene ID, 2nd column is the # reads mapped onto the VIRGO database, 3rd column is the gene length. File was sorted by the # reads mapped. 
> head sampleName.out 
V1593031	2425	3663
V1607456	566	390
V1420626	348	1017
V1684181	327	1442
V1316648	324	480
V1541216	274	1386
V1064785	269	303
V1872750	246	1179
V1885045	227	1254
V1573770	226	291

## Step 2: repeat for all your sample sets. Below is a pseudo-code
> for sample in *.fq; do ./runMapping.step1.sh -r sample1.fastq -p $sample -d /virgo/path/; done

## step 3: summarize the stats of multiple samples
./runMapping.step2.sh -p /path/to/output/of/step1/ -d /path/to/VIRGO/

# for example:
./runMapping.step2.sh -p temp_mapping/ -d /path/to/VIRGO/

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
### MODULE 5: BLAST to VIRGO database ################################
######################################################################

## VIRGO was also made blastable in both nucleotide and amino acid. Just run your normal blastn or blastp command
## # blast 2.2.3+ compatible 


######################################################################
### MODULE 6: For advanced user ######################################
######################################################################

### A basic understanding of the bash script and simple mapping result 

# parameters you can tune:
1. "reads = " # default value is 2. 
# in step 1 scripts, the cutoff of the number of reads for calling a gene hit. 
2. "bowtie -l " # default value is 25
# in step 1 scripts, the -l is the seed length for read mapping. You can go lower to get more hits or higher to get more strict in mapping quality. 
3. "--VIRGO-path" # the full path to the VIRGO directory 
# hard code the path so that you don't need to specify every time you run it.
4. output directory and output names 
# yes the scripts hard-coded the output names, the functional annotation. Feel free to name it better. 
5. Annotation file
# With the original fasta file of the VIRGO database, feel free to generate your own annotation!

######################################################################  
### MODULE 7. Use jaccard clustering to generate protein cluster######
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

######################################################################  
### MODULE 8. Generate protein family from VIRGO results##############
######################################################################

## The idea is simple. Now you have a list of VIRGO for your sample, you simply pull the protein family information 

# for example, the first 2 genes are in the same protein family named "1", this protein family has 2 genes in it. 
> head 0.VOG.VIRGO.tbl.txt 
1	V1616449
1	V1199137
10	V1891647
10	V1633691
100	V1861392
100	V1870003
1000	V1074024
1000	V1015725
10000	V1651986
10000	V1663564

