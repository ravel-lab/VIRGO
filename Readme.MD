# Tutorial for VIRGO  & VOG

##### Author: Bing Ma
##### Email: igs.bma@gmail.com
&nbsp;
**VIRGO** is a non-redundant gene catalog for the microbial communities that inhabit the human vagina. It was built using a combination of metagenomic and urogential isolate genome sequences. VIRGO can be used to characterize the **taxonomic** and/or **functional** composition of vaginal  **metagenomic** and **metatranscriptomic** datasets. Available functional annotations include: KEGG, COG, eggNOG, CDD, GO, Pfam, TIGRfam, and EC number.
&nbsp;
**VOG** is a comprehensive collection of orthlogous genes for the vaginal microbiota as defined by a novel jaccard clustering method. This complementary database of amino acid sequences can be used to improve functional annotation, comparative genomics and evolution of vaginal orthologous protein families. 
&nbsp;

#### The complete VIRGO and VOG databases, along with the set of vaginal metagenomes used to build and test VIRGO are available for download at: http://virgo.igs.umaryland.edu/ 



*Please contact the author if you have any questions or suggestions. Thank you and Enjoy!*

---
**Required dependencies** 

> BLAST 2.2.3+ 
> Bowtie 1.1.2+
> GNU Awk 3.1+
> seqtk 1.0+

---
**Database structure**

- **_test_run** contains the test dataset and results
- **0_db** contains the searchable database files and sequences
    - VIRGO Bowtie 1 database
    - VIRGO Bowtie 2 database
    - VIRGO BLASTN database
    - VIRGO BLASTP database
- **1_VIRGO** contains the annotation files for the genes
    - ***0.geneLength*** contains gene length in bp
    - ***1.taxon.table*** contains taxonomic information
    - ***2.A.proteinFamily.txt*** contains annotations from 17 sources, one source per line
    - ***2.B.proteinFamily.txt*** contains annotations from 17 sources, one gene per line
    - ***3.eggnog.txt*** contains eggNOG annotations
    - ***4.GC.txt*** contains high or low gene count information indicating whether a gene has a preference (>95%) for high or low gene count communities
    - ***5.EC.txt*** contains Enzyme Commission number (EC number)
    - ***6.product.txt*** contains gene product annotations
    - ***7.rxn.txt*** contains KEGG reaction annotations
    - ***8.A.kegg.ortholog.txt*** contains KEGG ortholog annotations
    - ***8.B.kegg.module.txt*** contain the KEGG module annotations
    - ***8.C.kegg.pathway.txt*** contains the KEGG pathway annotations
- **2_fungal_phage** contains the annotations for 10,908 fungal and 15,965 bacteriophage genes
- **3_run_VIRGO** contains scripts used to analyze metagenomic and metatranscriptomic datasets 
- **4_Jaccard_index_clustering** contains the script used to perform the jaccard clustering to identify VOGs
- **5_VOG contains** the VOG sequences and linker file to VIRGO genes, contents includes:
    - ***0.VOG.VIRGO.tbl.txt*** link between VIRGO and VOG ids
    - ***1.VOG.size.txt*** size of the VOG protein family
    - ***2.VOG.GC.txt*** gene count annotations
    - ***3.VOG.funcat.txt*** eggNOG functional category
    - ***4.VOG.taxonomy.txt*** consensus taxonomy
    - ***5.VOG.alignmentScore.txt*** alignment score for the protein family

---

### Module 0: Test run VIRGO
This module outlines how test VIRGO to ensure it is working as intended.

**In directory _test_run contains the test dataset:** 

> cd _test_run

**Run step1 and step2 scripts making sure to specify the path that VIRGO was downloaded to (it should contain the database structure such as /path/to/VIRGO/0_db, and /path/to/VIRGO/1_VIRGO)**
> ./runTesting.step1.sh -1 sub1.fq -2 sub2.fq -p test -d /path/to/VIRGO/
> ./runTesting.step2.sh -p temp_testsample -d /path/to/VIRGO/

**You should see outputs of all files as below!**

> ls temp_testsample/

- summary.Abundance.txt
- summary.Percentage.txt
- summary.Count.txt      
- summary.geneRichness.txt 
- test.1.reads2ref  
- test.annotation.txt
- test.1.fq         
- test.2.fq            
- test.out

---

### MODULE 1: Mapping to VIRGO non-redundant gene database ####
This module demonstrates how to map reads from a metagenome or metatranscriptome sample to VIRGO and present the results.

*Please note VIRGO does NOT yet support paired-end mapping, please merge the paired end reads alone or together with single-end reads into one fastq file.*

***runMapping.Step1.sh***  
**Argument:**

>**-r** fastq file containing reads as single end
>**-p** sample prefix 
>**-d** path to VIRGO

**Example:**
>**./runMapping.step1.sh** **-r** sample1.fastq **-p** samplePrefix **-d** /full/path/to/VIRGO/

**The output should be a sampleName.out text file. The 1st column is the VIRGO gene ID, 2nd column is the number of reads mapped onto the VIRGO database, the 3rd column is the gene length. File is sorted by the number of reads mapped column.**

> **head** sampleName.out 
```
V1593031  2425  3663
V1607456  566 390
V1420626  348 1017
V1684181  327 1442
V1316648  324 480
V1541216  274 1386
V1064785  269 303
V1872750  246 1179
V1885045  227 1254
V1573770  226 291
```
**Repeat for all your sample sets. Below is an example in pseudo-code**

> for sample in *.fq; do **./runMapping.step1.sh** **-r** sample1.fastq **-p** $sample **-d** /virgo/path/; done

**Summarize the stats of multiple samples**

**./runMapping.step2.sh** **-p** /path/to/output/of/step1/ **-d** /path/to/VIRGO/

**Output includes:**

- summary.Abundance.txt: the number of reads per species
- summary.Count.txt: the number of genes per species
- summary.Percentage.txt: normalized abundance in percentage per species (sum up to 100)
- summary.geneRichness.txt: number of genes per sample
- summary.NR.abundance.txt: the number of reads per non-redundant gene
- gene.lst.txt: the list of non-redundant genes with gene length
- EggNOG.annotation.txt: EggNOG annotation file per sample
- EC.annotation.txt: list of non-redundant genes with EC number
- GC.txt: list of non-redundant genes with gene count category (HGC: high gene count, LGC: low gene count
- geneProduct.txt: list of non-redundant genes with gene product annotation
- Kegg.module.annotation.txt: list of non-redundant genes with KEGG module annotation, including module ID and annotation
- Kegg.ortholog.annotation.txt: list of non-redundant genes with KEGG ortholog (KO) annotation, including ortholog ID and annotation
- Kegg.pathway.annotation.txt: list of non-redundant genes with KEGG pathway annotation, including pathway ID and annotation
- proteinFamily.annotation.txt: list of non-redundant genes with protein family annotation, from the database of CDD, GO, Gene3D, Hamap, Interpro, MobiDBLite, PIRSF, PRINTS, Pfam, ProDom, ProSitePatterns, ProSiteProfiles, SFLD, SMART, SUPERFAMILY, TIGRFAM
- rxn.annotation.txt: list of non-redundant genes with KEGG reaction

---
### MODULE 2: For advanced user ###

*Requires  abasic understanding of the bash script and the simplified mapping result* 

**Parameters you can tune:**

- **"reads** = " default value is 2. : *in step 1 scripts, the cutoff of the number of reads for calling a gene hit.**
- **"bowtie -l"** default value is 25 : *in step 1 scripts, the -l is the seed length for read mapping. You can go lower to get more hits or higher to get more strict in mapping quality.*
- **"--VIRGO-path"** the full path to the VIRGO directory : *hard code the path so that you don't need to specify every time you run it.*
- **output directory and output names** : *Yes the scripts hard-coded the output names, the functional annotation. Feel free to name it better.* 
- **Annotation file** : *Feel free to generate your own annotations for the original fasta file of the VIRGO database*

---
### MODULE 3: To explore a single gene of interest ###

It is easy to quickly pull annotations for your favorite gene, using a single VIRGO geneID! After you found a VIRGO gene hit to your favorite gene, by using "grep" with the annotation files of VIRGO, you will get lots of information, For example, for the gene V1000057, you can pull out its length, taxonomy, GO, TIGRFAM, Interpro, SUPERFAMILY, CDD, Pfam, eggNOG, EC number, gene product, KEGG ortholog, KEGG module, KEGG pathway etc.

>**grep** **-w** "V1000057" 1_VIRGO/*

- **0.geneLength.txt** :V1000057  990
- **1.taxon.tbl.txt** :Cluster_219398 V1000057  Lactobacillus_iners 990
- **2.A.proteinFamily.txt** :V1000057 GO  ; Phosphoribosyltransferase domain; Phosphoribosyltransferase-like; Ribose-phosphate diphosphokinase; Ribose-phosphate pyrophosphokinase, N-terminal domain
- **2.A.proteinFamily.txt** :V1000057 TIGRFAM ;  ribose-phosphate diphosphokinase
- **2.A.proteinFamily.txt** :V1000057 Interpro  ; Phosphoribosyltransferase domain; Phosphoribosyltransferase-like; Ribose-phosphate diphosphokinase; Ribose-phosphate pyrophosphokinase, N-terminal domain
- **2.A.proteinFamily.txt** :V1000057 SUPERFAMILY ; PRTase-like; PRTase-like
- **2.A.proteinFamily.txt** :V1000057 CDD ; PRTases_typeI
- **2.A.proteinFamily.txt** :V1000057 Pfam  ; N-terminal domain of ribose phosphate pyrophosphokinase; Phosphoribosyl synthetase-associated domain
- **2.B.proteinFamily.txt** :V1000057 ; PRTases_typeI ; Phosphoribosyltransferase domain; Phosphoribosyltransferase-like; Ribose-phosphate diphosphokinase; Ribose-phosphate pyrophosphokinase, N-terminal domain     ; Phosphoribosyltransferase domain; Phosphoribosyltransferase-like; Ribose-phosphate diphosphokinase; Ribose-phosphate pyrophosphokinase, N-terminal domain     ; N-terminal domain of ribose phosphate pyrophosphokinase; Phosphoribosyl synthetase-associated domain            ; PRTase-like; PRTase-like  ;  ribose-phosphate diphosphokinase
- **3.eggnog.NOG.txt** :Cluster_219398  V1000057  PRS map00030,map00230,map01100,map01110,map01120,map01230 F Phosphoribosyl pyrophosphate synthase COG0462
- **5.EC.txt** :Cluster_219398  V1000057  ec:2.7.6.1     ribose-phosphate diphosphokinase; ribose-phosphate pyrophosphokinase; PRPP synthetase; phosphoribosylpyrophosphate synthetase; PPRibP synthetase; PP-ribose P synthetase; 5-phosphoribosyl-1-pyrophosphate synthetase; 5-phosphoribose pyrophosphorylase; 5-phosphoribosyl-alpha-1-pyrophosphate synthetase; phosphoribosyl-diphosphate synthetase; phosphoribosylpyrophosphate synthase; pyrophosphoribosylphosphate synthetase; ribophosphate pyrophosphokinase; ribose-5-phosphate pyrophosphokinase
- **6.product.txt** :Cluster_219398 V1000057  Ribose-phosphate pyrophosphokinase
- **8.A.kegg.ortholog.txt** :V1000057 K00948  ljf:FI9785_163  dvvi:GSVIVP00018977001  GSVIVT00018977001; assembled CDS; K00948 ribose-phosphate pyrophosphokinase [EC:2.7.6.1]
- **8.B.kegg.module.txt** :V1000057 M00005  PRPP biosynthesis, ribose 5P => PRPP
- **8.C.kegg.pathway.txt** :V1000057  map00230  Purine metabolism

---

### MODULE 5: Exploring a list of genes of interest###

It is possible to pull out information for a list of genes of interest, using a list of VIRGO ID

**Step 1. Generate a list of genes by their VIRGO ID. e.g**

> cat lst

```
V1891288
V1891259
V1891167
V1891137
V1891104
```

**Step 2. Pull annotations for this list of genes, such as taxonomy, functional annotation, eggNOG, KEGG, etc.**

> awk -F"\t" 'NR==FNR{a[$2]=$3;next} ($1 in a) {print $0"\t"a[$1]}' 1_VIRGO/1.taxon.tbl.txt lst
```
V1891288  Lactobacillus_crispatus
V1891259  Lactobacillus_crispatus
V1891167  Lactobacillus_crispatus
V1891137  Lactobacillus_crispatus
V1891104  Lactobacillus_crispatus
```

---

### MODULE 6: BLAST to VIRGO database ###

**VIRGO was also made blastable in both nucleotide and amino acid. Just run your normal blastn or blastp command**

blast 2.2.3+ compatible 

---
### MODULE 7. Using jaccard clustering to generate protein cluster###

**Step 1**
Run all-vs-all BLASTP on set of pro 

**Step 2**
Run Jaccard clustering using modified clusterBsmlPairwiseAlignments script

>./run-jaccard.sh

Script/programs:
  
 *clusterBsmlPairwiseAlignments-mod.pl*

**Inputs:** 
  - Wu-blastp.bsml.list (list of BLASTP BSML files)
  - Location of executable that generates clusters from polypeptide pairs
  - Various parameters (linkscore, percent_identity, percent_coverage, p_value)

**Outputs:**
  - Jaccard-clusters.out(.gz)

*clusterBsmlPairwiseAlignments-mod.pl:*

  - Reads and filters all of the BLASTP data and generates a list of polypeptide pairs
  - Writes the polypeptide pairs to a file (20808.pairs.gz)
  - Invokes /usr/local/projects/ergatis/package-latest/bin/cluster to convert the pairs
    into clusters

*Note that all of clustering/filtering parameters are applied during the pair generation
step, not the pair clustering step.*

Example: 

>clusterBsmlPairwiseAlignments-mod.pl \
 --bsmlSearchList=wu-blastp.bsml.list \
 --log=test-jaccard.log \
 --linkscore=0.6 \
 --percent_identity=80 \
 --percent_coverage=60 \
 --p_value=1e-5 \
 --outfile=jaccard-clusters.out \
 --cluster_path=/path/to/cluster/files/

**Step 3. Run bidirectional best hit analysis on Jaccard clusters.**

> ./run-best-hit.sh 

**3a. Convert BLASTP BSML to a single btab file of best hits.**

**Script/program:**
  *CogBsmlLoader-mod.pl*

**Inputs:**
  - Wu-blastp.bsml.list (list of BLASTP BSML files)
  - Jaccard-clusters.out from Jaccard clustering step
  - Parameters (minimum p-value and minimum % coverage)

**Outputs:**
  - Summary btab file (/path/to/btab)
    *(which includes the Jaccard cluster info. in the last column)*

**3b. Merge all clusters with reciprocal best hits.**

**Script/program:**
  *best_hit.pl*

**Inputs:**
  - All-best-hits.btab from step 3a
  - Minimum Jaccard coefficient (for additional cluster pruning if set to > 0)

**Outputs:**
  -A final cluster file (final-clusters.cog.gz)

**4. Generate FASTA files.**

**Script/program:**
  *makeClusterFASTAFiles.pl*

**Inputs:**
  - The original polypeptide file (e.g., total.faa)
  - Kaccard-clusters.out from step 2
  - Final-clusters.cog from step 3b
  - The path to an output directory

**Outputs:**
  - A FASTA file for each cluster of size > 1
  - A singletons.fa file that contains all the clusters of size 1
  - A cluster size histogram (on stdout - see make-fasta-files-2.log.gz)

