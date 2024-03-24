## Training Dataset
### Reference genome
```
wget -O my_sequence.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta&retmode=text"

```
### V3 primer schemes: https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019
```
wget https://github.com/artic-network/artic-ncov2019/blob/master/primer_schemes/nCoV-2019/V3/nCoV-2019.primer.bed
```

## Install workflow toolkits
### Tools required: porechop, trimmomatic, ivar, medaka [tabix, pbgzip, bcftools, samtools, minimap2]
```
sudo apt update
sudo apt-get install python3-pip
sudo apt-get install build-essential
pip install --upgrade pip
pip install --upgrade pip setuptools
conda create -n SARS python=3.10
conda activate SARS
pip install medaka
conda install {required tools} -c bioconda

```
### Trim adapters from the reads
```
porechop -i file.fastq -o file_arm.fastq
```
### Trimmomatic is a widely used tool in bioinformatics for trimming and filtering raw sequencing reads
It calculates the average quality score in a sliding window of a specified size along each read and trims bases from the end of the read if the average quality score falls below a certain threshold
the SLIDINGWINDOW option is used with the following syntax: SLIDINGWINDOW:<windowSize>:<requiredQuality>
SE indicates that you're processing single-end reads, -phred33 specifies the quality encoding of the input FASTQ file
```
trimmomatic SE -threads 23 -phred33 file_arm.fastq SLIDINGWINDOW:50:10 MINLEN:100

```
### Map sequences to reference: 
-ax map-ont: Specifies the mapping mode suitable for long reads
```
minimap2 -ax map-ont ./reference.fasta file_arm.fastq -o file_aln.sam

```
### 
Use samtools to convert the sam to bam, then sort the reads
Use ivar to trim primers from your sequencing reads(primer bed file used depends on the primer vesrion used to sequence)
Then sort and index with samtools
```
samtools view -bS file_aln.sam > file_aln.bam
samtools sort file_aln.bam > file_sorted.bam
ivar trim -e -i file_sorted.bam -b ./nCoV-2019.primer.bed -p file_primertrim.bam 
samtools sort file_primertrim.bam -o file_primertrim_sorted
samtools index file_primertrim_sorted. bam
```
### Create consensus genome using ivar; usually samtools mpileup output is piped to ivar consensus
-aa: This option tells samtools mpileup to output all positions, including those with zero coverage. Without this option, samtools mpileup would only output positions covered by reads
-A: This option disables the BAQ (Base Alignment Quality) computation. BAQ is used to improve alignment quality by taking into account nearby base qualities. 
-d 0: This option sets the minimum read coverage depth for a position to be considered. In this case, setting -d 0 means that all positions, even those with zero coverage, will be included in the pileup output.
-Q 0: This option sets the minimum base quality for a base to be considered. By setting -Q 0, you include all bases regardless of their base quality scores.
```
 samtools mpileup -aa -A -d 0 -Q 0 file_primertrim_sorted.bam | ivar consensus -p file_consensus.fa

```
### Polish consensus genome using medaka
```
medaka consensus file_consensus.fa file_polished.fasta  -model r941_min_high_g360 -batch 200 -threads 4

