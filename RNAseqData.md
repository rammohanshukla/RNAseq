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
Load the private module RNASEQ/1.0.0 (has all the tools required for alignment, alignment QC and further anlysis)
```Tsch
module load RNASEQ/1.0.0
```
create new screen with session name
```Tsch
screen -S your_session_name
```
Download fastq.gz files from BASESPACE using Basemount
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
Run FASTQC on fastq.gz (needs Java 1.8) 
```Tsch
module load JAVA/1.8.0
```
Run FastQC and save the output to <target directory>
```Tsch
fastqc –O “<target directory>” *.fastq.gz
```
Read more on FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/





