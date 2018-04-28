Cluster commands – do all analyses in scratch (500gb)
Login to SSH
```Tsch
ssh -X rshukla@scclogin.camhres.ca
```
Then login to development node 1 or 2 (dev01 or dev02)
```Tsch
ssh -X dev01
```
-X stands for X server, mostly used for any any program having GUI interphase
Check for avilable modules
```Tsch
module avail
```
Load the private module RNASEQ/1.0.0 (has all the tools required for alignment and alignment QC and further anlysis)
```Tsch
module load RNASEQ/1.0.0
```
create new screen with session name
```Tsch
screen -S your_session_name
```
BASESPACE
Create “basespace” folder: 

```Tsch
mkdir basespace
```
Link basemount here: 
```Tsch
basemount basespace
```
Change directories:
```Tsch
cd basespace
```
Search for all .fastq.gz:
```Tsch
find ./ -name "*.fastq.gz" 
```
Find/Copy to new the target directory
```Tsch
find ./ -name "*.fastq.gz" -exec cp \{\} <Target directory path> \;
```
FASTQC
Make sure RNASEQ module is loaded: module load RNASEQ/1.0.0
Load Java 1.8: module load JAVA/1.8.0
Run FastQC: fastqc –O “<target directory>” *(or *.fastq.gz)

REFERENCE GENOME AND ANNOTATIONS
Get reference genome from: ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/
Download to “reference” folder: cd <reference dir>
wget -c ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
Get annotations from: ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/
wget –c ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.gtf.gz
Extract both reference genome and annotations file:
gunzip Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.86.gtf.gz
Make new folder: GRCh38 – subfolders: hisat2, fasta, genes
Move annotations(.gtf) to genes: mv “file” <destination>
Move reference genome to fasta
Run his: **30 mins
hisat2-build -p 16 “reference genome file” <hisat2 dir>
Get preliminary alignment statistics
Unpaired (one off)
hisat2 -x  $INDEX/XY --known-splicesite-infile $REFERENCE/splicesites.tsv -p 8 -U F1R10_S10_L001_R1_001.fastq.gz,F1R10_S10_L002_R1_001.fastq.gz | samtools view -bS - | samtools sort > $BAM_U/F1R10.bam
Paired (one off)
hisat2 -x  $INDEX/XY --known-splicesite-infile $REFERENCE/splicesites.tsv -p 8 -1 F1R10_S10_L001_R1_001.fastq.gz,F1R10_S10_L002_R1_001.fastq.gz -2 F1R10_S10_L001_R2_001.fastq.gz,F1R10_S10_L002_R2_001.fastq.gz | samtools view -bS - | samtools sort > $BAM_P/F1R10.bam
Parallel Processing (to run all subjects, do this)
Open “hisat_cmd_gen.bash”
Put proper file directories for sources/output etc
Change 2nd line commented SRC_DIR
Put in main directory
Load parallel processing module: module load GNU_PARALLEL
Create cmd.list: ./hisat_cmd_gen.bash > cmd.list
Run parallel processing: qbatch    -c4 -j4 --ppj 8  -w 24:00:00 cmd.list




