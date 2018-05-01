Make sub folders for Alignments Index, reference genome, and annotations within GRCh38 folder
```Tsch
mkdir -p GRCh38/{hisat2,fasta,genes}
```
make shortcut for each folder
```Tsch
export REFERENCE=/Path for fasta folder/
export GTF=/Path for genes folder/ 
export INDEX=/Path for Alignment index/
export BAM_P=/bg01/homescinet/scinet/course/ss2017/15_rnaseq/15-RNASeq/Questions/BAMFiles_Paired
export BAM_P=/bg01/homescinet/scinet/course/ss2017/15_rnaseq/15-RNASeq/Questions/BAMFiles_UnPaired
```

BRING REFERENCE GENOME AND ANNOTATIONS

reference genome from: ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/

Get annotations from: ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/

Download to $REFERENCE
```Tsch
cd $REFERENCE
```
```Tsch
wget -c ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
```
Extract both file:
```Tsch
gunzip Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
```
```Tsch
cd $GTF
```
```Tsch
wget –c ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.gtf.gz
```
Extract file:
```Tsch
gunzip Homo_sapiens.GRCh38.86.gtf.gz
```
```Tsch
cd $REFERENCE  
```
```Tsch
hisat2-build -p 16 <reference genome file> $INDEX
```
Extract splice sites
```Tsch
cd $GTF
```
```Tsch
hisat2_extract_splice_sites.py $GTF <target_dir/splicesites.tsv>
```
Run hisat2

cd to directory with FastQ files
```Tsch
hisat2-build -p 16 $REFERENCE/genomefile <hisat2 dir>
```
Get preliminary alignment statistics

Unpaired (one off)
```Tsch
hisat2 -x  $INDEX/XY --known-splicesite-infile $REFERENCE/splicesites.tsv -p 8 -U F1R10_S10_L001_R1_001.fastq.gz,F1R10_S10_L002_R1_001.fastq.gz | samtools view -bS - | samtools sort > $BAM_U/F1R10.bam
```
Paired (one off)
```Tsch
hisat2 -x  $INDEX/XY --known-splicesite-infile $REFERENCE/splicesites.tsv -p 8 -1 F1R10_S10_L001_R1_001.fastq.gz,F1R10_S10_L002_R1_001.fastq.gz -2 F1R10_S10_L001_R2_001.fastq.gz,F1R10_S10_L002_R2_001.fastq.gz | samtools view -bS - | samtools sort > $BAM_P/F1R10.bam
```
Parallel Processing (to run all subjects, do this)

Open “hisat_cmd_gen.bash”

Put proper file directories for sources/output etc

Change 2nd line commented SRC_DIR

Put in main directory

Load parallel processing module: 
```Tsch
module load GNU_PARALLEL
```
Create cmd.list: 
```Tsch
./hisat_cmd_gen.bash > cmd.list
```
Run parallel processing: 
```Tsch
qbatch -c4 -j4 --ppj 8  -w 24:00:00 cmd.list
```
Post-alignment stats (% aligned)

Create directory for output in bam folder: mkdir flagstat_output

Run “flagstat_miseq”: 
```Tsch
./flagstat_miseq.bash > cmd.list
```
Edit .bash file to make sure directories are linked correctly

Run cmd.list: 
```Tsch
./cmd.list (may need to change permissions: chmod 744 cmd.list)
```
Output of each file will be in flagstat_output

To get a convenient list of mapping rates, run flagstat_extract script in flagstat_output directory: 
```Tsch
./flagstat_extract.bash 
```
If viewing alignment output in genomic alignment viewer (IGV) we need to index the files 

Index bam files: 
```Tsch
samtools index –b <filename>
```
