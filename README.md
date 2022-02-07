# 'DNA damage independent inhibition of NF-κB transcription by anthracyclines' - Chora *et al.*, 2022
## RNA-Seq Analysis

This analysis is part of this published pre-print: https://www.biorxiv.org/content/10.1101/2020.04.27.065003v2

This repository contains transcriptomic data from murine Bone Marrow Derived Macrophages (BMDM's) stimulated with LPS and treated with PBS, Epirubicin and Aclarubicin

--------------------------------------------------------------------------------------------------------------------------------------------------------------------
All samples were processed the following way:

#### Quality Control
This quality control is done using the program **Fastqc** (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/); considering this program produces a report per file, we can also use a program to merge all the reports into one - **MultiQC** (https://multiqc.info/).

```
mkdir fastqc_reports
mkdir mQC_reports

fastqc raw/*.txt.gz -o fastqc_reports -t 5
```
#### Alignment
For the alignment steps, we're going to use **STAR Alignment** (https://github.com/alexdobin/STAR)

###### 1) Create the Reference Genome Indexes

First of all, create a new folder to store the genome indexes and also create a shortcut for it:
```
mkdir gen_index
export gen_index=./gen_index
```

Then, you can assess it by just using:

```
cd $gen_index
```

To align our samples, we need the reference genome indexes, that are just the corresponding reference genome (.fasta) and its corresponding annotation (.gtf). We used the genome and correspondig annotation for the version 89.

```
wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

wget ftp://ftp.ensembl.org/pub/release-89/gtf/mus_musculus/Mus_musculus.GRCm38.89.gtf.gz
gunzip Mus_musculus.GRCm38.97.gtf.gz
```
 Now that we have the files, we proceed to use STAR with the option of  `genomeGenerate`

```
STAR --runThreadN 10 --runMode genomeGenerate --genomeDir $gen_index --genomeFastaFiles $gen_index/Mus_musculus.GRCm38.dna.primary_assembly.fa --sjdbGTFfile $gen_index/Mus_musculus.GRCm38.97.gtf 

```
Options explained:

- `--runThreadN 10` Nr. of Cores used
- `--runMode genomeGenerate` Argument for the program to know what is going to run
- `--genomeDir $gen_index/` Path to the genome indexes
- `--genomeFastaFiles $gen_index/Mus_musculus.GRCm38.dna.primary_assembly.fa` Path to the fasta file of the reference genome
- `--sjdbGTFfile $gen_index/Mus_musculus.GRCm38.95.gtf` Path to the gtf file of the reference genome

Afterwards, you can proceed to the alignment step. We are going to use the program qualimap (http://qualimap.bioinfo.cipf.es/) & the log files of STAR to assess the quality of this step.

The alignment step will perform a cycle but, before, we need to create the several required folders to store the necessary files:
```
mkdir aligned
mkdir qualimap_aligned

```
Now, we go to the directory where the merged files are and proceed with the cycle:
```
cd raw/
export gen_index=../../gen_index #this step helps in case of using a different window
```
```
for f in *.txt.gz;
do STAR --genomeDir $gen_index --readFilesIn $f --readFilesCommand zcat --sjdbGTFfile $gen_index/Mus_musculus.GRCm38.97.gtf --quantMode GeneCounts --runThreadN 24 --outFileNamePrefix ../aligned/"$f"_ --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx;
qualimap rnaseq -bam ../aligned/"$f"_*.bam -gtf $gen_index/Mus_musculus.GRCm38.97.gtf -outdir ../qualimap_aligned/$f --java-mem-size=40G;
rm -rf ../aligned/$f*.bam; done
```

This cycle will perform the alignment and perform the corresponding report for the BAM file.

After checking the qualimap report and the log.final.out file for each sample if everything went as expected during the alignment, we proceeded to the R Analysis. The steps corresponding to this analysis are stored in the 'R Analysis' folder.
