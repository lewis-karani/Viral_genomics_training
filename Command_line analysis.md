## Training Dataset
### Reference genome
```
wget -O my_sequence.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta&retmode=text"

```
### Reverse reads
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR209/005/ERR2093245/ERR2093245_2.fastq.gz
```

## Install workflow toolkits
### Tools required: porechop, trimmomatic, minimap2, samtools, ivar, medaka [tabix, pbgzip, bbtools, samtools, ]
```
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
```
### Extract the content of tar file
```
tar -vxzf sratoolkit.tar.gz
```
### Export to PATH
```
export PATH=$PATH:$PWD/sratoolkit.3.1.0-ubuntu64/bin
```
### Refresh the terminal
```
source ~/.bashrc
```
### Verify successful install
```
fastq-dump
```
