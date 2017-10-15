* **Allen Lab pipeline**
* Jessica Blanton
* Last updated: Jan 2017

# 16S amplicon analysis pipeline


### Outline of pipeline
   **Getting the files and metadata organized**
   
1. Obtain and verify data 
2. Create metadata mapping file
3. Make softlinks of raw data to working directory

   **Preprocessing your sequence files**

4. Trim Sequences with trimmomatic 
5. Merge Paired-Read sequence files with PEAR
6. Filter sequences by qscore and size, with renaming option
7. Remove primer sequences with the cutadapt program
8. Filter sequences for size, remove homopolymers, & change headers back again
9. Split chimeric sequences from no-chimeras with the vsearch program

   **Data Processing**

10. Pick otus (& much more)
11. Assign taxonomy to representative sequences with Qiime using the RDP classifier, .80 bootstrap cutoff
12. Build/modify final biom table 
13. Package and export essential files 

Downstream Processing: software recommendations


### Outline of Databases used
- **rdp.gold** 
	- Download from Broad: https://sourceforge.net/projects/microbiomeutil/files/microbiomeutil-r20110519.tgz/download)
- **SILVAv128**: 
	- Download current Qiime-formatted releases: https://www.arb-silva.de/download/archive/qiime/)
	
	- **NOTE: I have created custom databases which include novel sequence groups which are currently not in SILVA, or that have corrected taxonomic strings.
	- Otu picking: `/rep_set/rep_set_all/97/97_otus.fasta`
	- Taxonomic assignment: `/taxonomy_all/99/
majority_taxonomy_7_levels.txt`
	- alignment for tree building: `/rep_set/rep_set_all/99/99_otus_aligned.fasta`

### Programs used
- Qiime v1.9.1 (command line version)
- Cutadapt 1.11
- Vsearch v2.3.0_linux_x86_64


##### Occasional errors with poor documentation 
- Parallel scripts in Qiime often max out memory of any computer, and cause processes to hang
- Sometimes it is best to import a json formatted biom file into R through phyloseq, rather than the more update hd5f format

---
# Running the Pipeline:

## Getting the files and metadata organized
*an ounce of planning...*

#### 1. Obtain and verify data 
Know this often requires disentanglingl it from other collaborator's datasets, or partial runs at the sequencing Core

- Are the number of files and names correct?

#### 2. Create metadata mapping file
*This is strongly recommended even if not using a program that immediately requires it!*

Include at least the following essentials:

- Exact Name of datafile e.g. "JB009"
- Name of source specimen
- Descriptive name of sample set e.g. "SL12_K-HL"
- Experimental conditions
- Primer names for the target-specific PCR

#### 3. Make softlinks of raw data to working directory
This protects your raw-data archive

## Preprocessing your sequence files

#### 4. Trim Sequences with trimmomatic 
Optional step, only if data is particulary bad, or will be analyzed as single reads instead of paired

#### 5. Merge Paired-Read sequence files with PEAR
```
pear -f 021515_FGM3PE300/LP46*R1* -r 021515_FGM3PE300/LP46*R2* \
-o V34_pearmerge/LP46_merged \
-q 35 -t 50 -j 12 >> V34_pearmerge/output.txt
```

#### 6. Filter sequences by qscore and size, with renaming option
- **qiime** (recommended for ease of pipeline integration)

	\* or *
- cutadapt
- vsearch

**Preparing mapping files for** ```split\_libraries_fastq.py```:

The Qiime pipeline expects *only* multiplexed data. For demultiplexed data as is received from a MiSeq run, you essentially must work around this by treating each sample like an independent dataset at this stage.  This is somewhat onorous as indicated in this pipeline; code should be re-written using loops.

Make master mapfile 

(if renaming desired, indicate so here by entering filenames in the "run_prefix" column, and desired name in "SampleID" column)

```
nano V34_addition_mapfile.txt
		#SampleID	BarcodeSequence	LinkerPrimerSequence	ReversePrimer	InputFileName	run_prefix	EpiLum_portion	MidHindTotalGutPortion	DissectionMethod_FieldUCSBLabJBLabUCSB	Catch	run_number	SourceSample	Species	Description
		LP44		GTGYCAGCMGCCGCGGTAA	CCGYCAATTYMTTTRAGTTT	LP44.extendedFrags.fastq	LP44	lumenal	total	LabLP	SJ16	FGM3	SJ16A	Scomber.japonicus	SJ16A.L.D.Jstore
		LP45		GTGYCAGCMGCCGCGGTAA	CCGYCAATTYMTTTRAGTTT	LP45.extendedFrags.fastq	LP45	lumenal	total	LabLP	SJ16	FGM3	SJ16B	Scomber.japonicus	SJ16B.L.D.Jstore
		LP46		GTGYCAGCMGCCGCGGTAA	CCGYCAATTYMTTTRAGTTT	LP46.extendedFrags.fastq	LP46	lumenal	total	LabLP	SJ16	FGM3	SJ16D	Scomber.japonicus	SJ16D.L.D.Jstore
		LP43		GTGYCAGCMGCCGCGGTAA	CCGYCAATTYMTTTRAGTTT	LP43.extendedFrags.fastq	LP43	lumenal	total	LabLP	SJ16	FGM3	SJ16F	Scomber.japonicus	SJ16F.L.D.Jstore
```
Qiime-validate master mapfile, passing flags to ignore barcode (-b) and primer sequences (-p)

```
source /opt/miniconda3/bin/activate qiime1
validate_mapping_file.py -m *mapfile*.txt -o mapfile_out/ -p -b
```
Check that there is one a single expected error relating to barcodes:
"If no barcodes are present, and the added_demultiplex_field option isn't used, only a single SampleID can be present.   no location"

```
cat mapfile_out/*.log
# Errors and warnings are written as a tab separated columns, with the first column showing the error or warning, and the second column contains the location of the error or warning, written as row,column, where 0,0 is the top left header item (SampleID).  Problems not specific to a particular data cell will be listed as having 'no location'.
# Errors -----------------------------
# If no barcodes are present, and the added_demultiplex_field option isn't used, only a single SampleID can be present.   no location
# Warnings ---------------------------
```
Split validated mapfile to 1 per sample

```
mkdir per_sample_maps/
grep "#SampleID" *mapfile*.txt > per_sample_maps/LP43map.txt | grep "LP43" *mapfile*.txt >> per_sample_maps/LP43map.txt
grep "#SampleID" *mapfile*.txt > per_sample_maps/LP44map.txt | grep "LP44" *mapfile*.txt >> per_sample_maps/LP44map.txt
grep "#SampleID" *mapfile*.txt > per_sample_maps/LP45map.txt | grep "LP45" *mapfile*.txt >> per_sample_maps/LP45map.txt
grep "#SampleID" *mapfile*.txt > per_sample_maps/LP46map.txt | grep "LP46" *mapfile*.txt >> per_sample_maps/LP46map.txt

```
Filter sequences by qscore ≥ 20, and optionally rename sequences with the Qiime script ```split_libraries_fastq.py```

```
cd /home/jess/projects/amplicons/amplicons_053016/053016_V34/122716_DNAcomp_add
nohup split_libraries_fastq.py \
-i /home/jess/projects/raw_reads/miseq/V34_addition_pearmerge/LP43_merged.assembled.fastq,\
/home/jess/projects/raw_reads/miseq/V34_addition_pearmerge/LP44_merged.assembled.fastq,\
/home/jess/projects/raw_reads/miseq/V34_addition_pearmerge/LP45_merged.assembled.fastq,\
/home/jess/projects/raw_reads/miseq/V34_addition_pearmerge/LP46_merged.assembled.fastq \
--sample_id LP43,LP44,LP45,LP46 \
-o combined_seqs \
-m per_sample_maps/LP43map.txt,per_sample_maps/LP44map.txt,per_sample_maps/LP45map.txt,per_sample_maps/LP46map.txt \
-q 19 \
--barcode_type 'not-barcoded' &
```

#### 7. Remove primer sequences with the cutadapt program

```
mkdir PrimersRemoved/

cutadapt \
   -g 341F=^CCTACGGGNGGCWGCAG \
   combined_seqs/seqs.fna \
   | cutadapt \
   -a 785R=GGATTAGATACCCBDGTAGT$ \
   - > PrimersRemoved/primertrimmed_withempties.fasta;

# Remove the empty sequences that sometimes appear:
# Removing primers with cutadapt sometimes results in sequences with length=0, which will throw errors downstream.

nohup awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' \
PrimersRemoved/primertrimmed_withempties.fasta > PrimersRemoved/primertrimmed.fasta &
```
Check histograms before filtering seqs again
	
	less combined_seqs/histograms.txt

#### 8. Filter sequences for size, remove homopolymers, & change headers back again
using `split_libraries.py`

```
nohup split_libraries.py \
   -m *mapfile*.txt \
   -f PrimersRemoved/primertrimmed.fasta \
   -b 0  \
   -p \
   -l 200 \
   -L 600 \
   -o PrimersRemoved/ \
   -j "run_prefix" &
```
#### 9. Split chimeric sequences from no-chimeras with the vsearch program

* **NOTE: use "--fasta_width 0" to make sure sequence output is a single line (not wrapped)
* **NOTE: this step can take a long time!!!

```
vsearch --uchime_ref PrimersRemoved/seqs.fna \
--db /home/jess/scratch/ribosomal_databases/qiimechimeraDB_RDBgold9/rdp_gold.fa \
--chimeras sorted_chimera_out.fna --nonchimeras seqs_chimeras_filtered.fna \
--threads 10 --minh 0.2 --mindiv 1.5 \
--log vsearch_chimeras.log --fasta_width 0

```

## Data Processing
The following steps things constitute decisions for downstream analysis and in many cases choices will be are specific to your amplicon data!


#### 10. Pick otus (& much more)
1. **Qiime (Uclust default)**
1. Vsearch (USEARCH algorithm) 
1. Swarm clustering
1. Deblur high resolution method
1. DaDa2 high resolution method


**Qiime Open-reference clustering (used here mainly out of convenience)**
This script includes many tasks which could be performed separately:

1. Pick otus
2. Pick rep set
3. Filter singleton OTUs
4. Pynast_align/filter failures from rep set
5. Filter alignment gaps
5. Build tree from final alignment
6. Assign taxonomy (recommend doing this separately if using resource-intensive RDPclassifier algorithm)
7. Build final biom table)


Make parameter file pointing to appropriate SILVA databases for alignment (and taxonomic classifications if performing this step here). 

*Alternatively, clustering with the open reference workflow script can be done with skipping taxonomic assignment or alignment/tree building. These resource-consuming steps may be done separately if needed.*

- **NOTE: it is recommended that you avoid using Qiime's parallel scripts and multithreading, as they have a tendency to multiply the workload too much, causing most computers to hang

```
echo 'pick_otus:similarity 0.98
assign_taxonomy:id_to_taxonomy_fp $FP/Silva_123_Custom_taxamap_JBv9.txt
assign_taxonomy:reference_seqs_fp /$FP/99_Silva_123_Custom_repset_JBv8.txt
align_seqs:template_fp $FP/Silva_111_post/rep_set_aligned/97_Silva_111_rep_set.fasta
align_seqs:min_length 200

' >parameters_open_98pick_otus_notax.txt;
```
Run pick-otus workflow without taxonomy assignment

```
pick_open_reference_otus.py \
-i VsearchChecked/seqs_chimeras_filtered.fna \
-o picked_98ORotus_silva \
-m uclust \
--suppress_taxonomy_assignment \
-p parameters_open_98pick_otus_notax.txt
```
Check number of OTUs formed by counting the representative sets selected

	grep ">" rep_set.fna -c



#### 11. Assign taxonomy to representative sequences with Qiime using the RDP classifier, .80 bootstrap cutoff
This is remarkably memory intensive- the bigger the DB, the more pairwise comparisons per repseq!  Takes a long time, if using SILVAv123 or higher, I run this with more memory on the triton supercomputer.  Alternative is to run it with a less comprehensive DB.

```
assign_taxonomy.py \
-m rdp -c 0.8 \
--rdp_max_memory 1000000 \
-i rep_set.fna \
--id_to_taxonomy_fp $FP/qiimecustomDB/Silva_123_Custom_taxamap_JBv9.txt \
--reference_seqs_fp $FP/qiimecustomDB/99_Silva_123_Custom_repset_JBv8.txt
```

#### 12. Build/modify final biom table 
Final table should have no singletons, be filtered of non-16S sequences, and have a taxonomy assigned

```
biom add-metadata \
-i otu_table_mc2_no_pynast_failures.biom \
--observation-metadata-fp rdp_assigned_taxonomy/rep_set_tax_assignments.txt \
-o otu_table_mc2_no_pynast_failures_rdptax.biom \
--sc-separated taxonomy \
--observation-header OTUID,taxonomy
```
Summarize your biom table and verify that it makes sense
- correct number of samples
- expected readcounts

```
biom summarize-table -i otu_table_mc2_w_rdptax_nopynast_failures.biom -o table_summary.txt  
```

#### 13. Compress and export essential files 

The objective here is to bundle up the files that an outside researcher can perform all possible downstream analysis, short of testing additional OTU clustering methods. Ideally, you will also include your exact preprocessing commands from this pipeline- even if you're the one doing the analyses, keeping a copy of your exact commands with your files will make future troubleshooting possible.

Create a tarball of all necessary files

```
tar -cvf 16S_011717.tar otu_table_mc2_no_pynast_failures.biom rep_set.tre V34_addition_mapfile.txt rep_set.fna 16S_011717_cmds.txt
gzip 16S_011717.tar

```
File key:

|Export file|example filename|
|----|------------|
|Count+taxonomy table (biom format) |`otu_table_mc2_no_pynast_failures.biom`|
|Treefile |`rep_set.tre`|
|Master mapping file|`V34_addition_mapfile.txt `|
|Representative sequences for OTUs|`rep_set.fna`|
|Your preprocessing commands|`16S_011717_cmds.txt`|

### Fin.

--- 
# Downstream Processing: software recommendations 

While Qiime supports an increasing number analysis scripts, the packages noted below have been found to provide more customizable and agile methods for analysis. 
Notably, the *Phyloseq* R package is excellent for data management, as all data from a study is stored in, and accessed from a single S4 object.

Moving into R also allows for easier use the numerous statistics and graphing packages available.  Furthermore, other useful programs which are not currently integrated into qiime are made easier to use using R to manage and transform data. A similar case can be made for up-and-coming python programs, so the recommendations here are solid but certainly the only game in town.

|Program|Runs as|Features
|-------|-------|-------|
Phyloseq v1.19.1 | R package| data management, extensive diversity analyses, wrappered plotting functions|
metagenomeSeq v1.16.0 |R package| cumulative-sum scaling (CSS) method to normalize data
DESeq2 1.16.1| R package| Differential abundance method based on the negative binomial distribution
SPARCC| python2.x script| python module for computing correlations in compositional data

