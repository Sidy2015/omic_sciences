



# Introduction to FASTQ 
The first step in NGS sequencing is to read the FASTQ files and align them. We're not going to do that because the FASTQ files are too big -- but we will walk through the steps. There are many many software programs for alignment. If you want to run a truncated version to see what FASTQ files look like, you can download a single sample in a UNIX shell that we will see later here:


## Links for this experiment

Study information at the Sequence Read Archive:

http://www.ncbi.nlm.nih.gov/Traces/sra/?study=SRP033351

Himes et al paper at PubMed Central:

http://www.ncbi.nlm.nih.gov/pubmed/24926665 

Example sample table stored in our course repo on github:

https://github.com/genomicsclass/labs/blob/master/course5/airway_sample_table.csv

Details on creating such a sample table from SRA and GEO:

http://www.bioconductor.org/packages/release/data/experiment/vignettes/airway/inst/doc/airway.html

The European Nucleotide Archive (EMBL-EBI):

http://www.ebi.ac.uk/ena 

The Sequence Read Archive (NCBI):

http://www.ncbi.nlm.nih.gov/sra/ 

## Fastq file commands

Downloading from the ENA:

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508_2.fastq.gz

Alias for ls:

alias ll='ls -lGh'

Unzipping:

gunzip *.fastq.gz


Looking at the FASTQ files:

less SRR1039508_1.fastq
wc -l SRR1039508_1.fastq ### number of lines


Quality control with [fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

Navigate to your FASTQ directory:

fastqc --noextract SRR1039508_1.fastq SRR1039508_2.fastq --outdir=/Users/gurinina/Dropbox/Documents/Rprojects/2022_GENOMICS_COURSE/fastqc

```{r}

library(fastqcr)

qc.dir = "/Users/gurinina/Dropbox/Documents/Rprojects/2022_GENOMICS_COURSE/fastqc"

qc <- qc_aggregate(qc.dir)
library(dplyr)
  tmp = qc %>%
    group_by(sample, module, status) %>%    
    filter(status %in% c("WARN", "FAIL")) %>%
    arrange(sample)
  # Summary of qc
summary(qc)
qc_stats(qc)
qc_fails(qc)#: Displays samples or modules that failed.
qc_warns(qc)#: Displays samples or modules that warned.
qc_problems(qc)#: Union of qc_fails() and qc_warns(). Display which samples or modules that failed or warned.
  qc_fails(qc, "module")
  qc_warns(qc, "module")
  qc_problems(qc, "module")
qc_problems(qc, "module", compact = FALSE)
qc_fails(qc, "sample")

##### THIS ONE IS MOST IMPORTANT; BUT ALL OF THIS IS IN THE HTML FILE
qc_report(qc.path=qc.dir, 
          result.file="reportFile", preview = TRUE)

qc_report(qc.dir, result.file = "/Users/gurinina/Dropbox/Documents/Rprojects/2022_GENOMICS_COURSE/fastq",
          experiment = "RNA seq hg19")
qc.file <- qc_aggregate(qc.dir)
qc <- qc_read(qc.file)
names(qc)
qc_plot(qc, "Per sequence GC content")

qc_plot(qc, "Per base sequence quality")

qc_plot(qc, "Per sequence quality scores")

qc_plot(qc, "Sequence duplication levels")

library(ShortRead)

fastqFile = "/Users/gurinina/Dropbox/Documents/Rprojects/2022_GENOMICS_COURSE/fastq/SRR1039508_1.fastq"
# read fastq file
fq = readFastq(fastqFile)

reads = sread(fq)[1:1000]
reads[2]
reverseComplement(reads[2])
alphabetFrequency(reads[1:5] )
w = vcountPattern("GCTGGGCG",reads)
ww = which(w!=0)
strsplit(reads[ww],"GCTGGGCG")

# get quality scores per base as a matrix
qPerBase = as(quality(fq[1:1e+06]), "matrix")


# get number of bases per read that have quality score below 20
# we use this
qcount = rowSums( qPerBase <= 20) 

# Number of reads where all Phred scores >= 20
length(qcount[qcount == 0])
## class: ShortReadQ
## length: 10699 reads; width: 72 cycles

```
~ 1.3 % of the reads have all 63 bp >= 20, is that enough?
what is the cutoff?
do we filter reads here or wait for the exon match to do it?
what is a per tile quality score, how to read the QC reports?
The FASTQ file is 22e+06 reads long, you could write a loop and filter out the poor quality reads
We can write out the filtered fastq file with the ShortRead::writeFastq() function
```{r}
# write out fastq file with only reads where all 
# quality scores per base are above 20
writeFastq(fq[qcount == 0], 
           paste(fastqFile, "Qfiltered", sep="_"))
           
```