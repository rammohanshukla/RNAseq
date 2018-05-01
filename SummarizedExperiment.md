```R
library(DESeq2)
library(Rsamtools);
library("GenomicFeatures");
library("GenomicAlignments");
library("BiocParallel");
library("pheatmap");
library("RColorBrewer");
library("PoiClaClu")
```
```R
getwd()
setwd(“target dir”)
list.files()
```

o	Set “search rate” for bam files (how many reads at a time it will look for)

```R
files<-list.files(pattern=".bam$");
file.exists(files)
bamfiles <- BamFileList(files, yieldSize=200000)
```
o	Prepare GTF file for counting, change file name depending on species

```R
seqinfo(bamfiles[1])
(txdb <- makeTxDbFromGFF("Homo_sapiens.GRCh38.86.gtf", format="gtf", circ_seqs=character())) 
(ebg <- exonsBy(txdb, by="gene"))
```
 
 #this is the counting step – run overnight (MAKE SURE SCREEN WAS RUN BEFOREHAND)
```R
register(MulticoreParam(workers=10))
se <- summarizeOverlaps(features=ebg, reads=bamfiles,mode="Union",singleEnd=FALSE,ignore.strand=TRUE,fragments=TRUE)
```


Once run, these commands let you look at the dimensions of the summarized experiment (se – aka your dataset)

```R
se 
```
```R
Output: 
class: RangedSummarizedExperiment 
dim: 58051 90 
metadata(0): #There is no metadata file 
assays(1): counts #There is 1 assay file (actual data) 
rownames(58051): ENSG00000000003 ENSG00000000005 ... ENSG00000283698
  ENSG00000283699
rowRanges metadata column names(0):
colnames(90): F1R10.bam F1R11.bam ... F8R5.bam F8R6.bam
colData names(0): #There is no Ptable
```
```R
dim(se) 
assayNames(se)
head(assay(se), 3) 
colSums(assay(se)) # this provides the total number of reads mapped to exons for each sample
round(colSums(assay(se)/1e6,1),2) #this gives the same information but in millions
mean(round(colSums(assay(se)/1e6,1),2)) #average number of reads mapped to exons across all samples  
RR <- round(colSums(assay(se)/1e6,1),2) #puts mapped reads, in millions, into the object RR – this can be renamed to your liking
write.csv(RR,file=”RRDepth.csv”) # outputs the object RR to a .csv, this also provides the order of the samples – important for integrating with the pTable, make sure your pTable has the names and order in its rows as in the assay file. Save as pTable.csv and upload to .bam directory of cluster.
```
#These next functions are to do checks – make sure the names of columns match up with the sampleTable (pTable)
```R
rowRanges(se)
str(metadata(rowRanges(se))) 
colData(se) 
```
```R
pTable<-read.csv("pTable.csv") 
colnames(se)
rownames(pTable) #check that the results from these two are the same
pTable<-DataFrame(pTable) 
(colData(se) <- DataFrame(pTable)) 
colData(se) 
se$Cell_type <- relevel(se$Cell_type, "Pyr") #not needed now, but later for differential expression
round(colSums(assay(se))/1e6, 1) #check that the column sums are still correct
```

#To get metadata
Go to: https://www.ensembl.org/biomart
-Select Human Genes (GRCh38.p12), and filter by Gene stable IDs
-Select Attributes: Gene stable ID, Transcript stable ID, Gene description, Strand, Gene name, Transcript length (including UTRs and CDS), Gene % GC content, Transcription start site (TSS), Chromosome/scaffold name, Transcript count
-Select any others you are interested in (there is a plethora of data to capture, if you are interested)
-Download, unzip using gunzip, and convert to .csv if required (may download as .txt even if .csv is selected)	


#To merge metadata with se
#Since not all genes will have annotations, we will separate the gene list into with and without annotations, apply the relevant annotations, recombine the files, then merge with se
```R
featureData <- rownames(se)
sum(duplicated(featureData)) #if this is not 0, there are duplicates
write.csv(featureData,file=”featureData.csv”)
mData <- read.csv(“metadata.csv”)
mData = mdata[!duplicated(mData$ Gene.stable.ID),]
noAnnotations <- setdiff(featureData, mData$Gene.stable.ID)
withAnnotations <- intersect(featureData, mData$Gene.stable.ID)
mTable1 <- mData[which(withAnnotations %in% mData$Gene.stable.ID),]
mTable2 <- as.data.frame(matrix(, nrow=486, ncol=10))
colnames(mTable2) <- colnames(mTable1)
mTable2$Gene.stable.ID <- noAnnotations
colnames(mTable1)==colnames(mTable2) # check that these are the same
x <- rbind(mTable1,mTable2)
dim(x) #check that the number of genes (rows) and metadata (columns) is correct
sum(duplicated(x$Gene.stable.ID)) #if this is not 0, there are duplicates
y <- x[match(as.factor(featureData), x$Gene.stable.ID),]
sum(as.factor(featureData)== y$Gene.stable.ID)
rownames(y) <- y$Gene.stable.ID
mcols(se) <- cbind(mcols(se), y)
save(se,file=”se.rData”)
```
