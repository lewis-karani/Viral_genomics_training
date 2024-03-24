## Training Dataset
### Forward reads
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR209/005/ERR2093245/ERR2093245_1.fastq.gz
```
### Reverse reads
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR209/005/ERR2093245/ERR2093245_2.fastq.gz
```

## Install Sra-toolkit
### Fetch the tar file from the canonical location at NCBI:
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
