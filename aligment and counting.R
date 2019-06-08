###############################################################################################
# RNA-seq analysis in R
# Alignment and Counting
# source: https://bioinformatics-core-shared-training.github.io/RNAseq-R/align-and-count.nb.html
################################################################################################
library(Rsubread)

##We create the list of all the reads we want to align
fastq.files <- list.files(path = "./data", pattern = ".fastq.gz$", full.names = TRUE)
fastq.files

###In order to build an index you need to have the fasta file (.fa), which can be downloaded from the UCSC genome browser. 
#buildindex(basename="/home/anna/R_exercises/RNAseq/chr1_mm10",reference="/home/anna/R_exercises/RNAseq/chr1.fa")

##The index is already built!!

## Now we can align our 12 fastq.gz files using the align command.
##The default output file name is the filename with “.subread.BAM” added at the end.

dir('data')

test <- align(index = "data/chr1_mm10", readfile1= fastq.files)#We have 4 index files in data
#directory. "index" arguments takes the chr1_mm10 part they have in common 
test

##############
# CHALLENGE  #
##############

##########################################################################################################
# Try aligning the fastq files allowing multi-mapping reads (set unique = FALSE), 
# and allowing for up to 6 “best” locations to be reported (nBestLocations = 6). 
# Specify the output file names (bam.files.multi) by substituting “.fastq.gz” with “.multi.bam” 
# so we don’t overwrite our unique alignment bam files.
# Look at the proportion of reads mapped and see if we get any more reads mapping by specifying a 
# less stringent criteria.
##########################################################################################################

test1 <- align(index = "data/chr1_mm10", readfile1= fastq.files, unique = FALSE, 
              nBestLocations = 6)
test1

# get a summary of the proportion of reads that mapped to the reference genome using 
# the propmapped function.

bam.files <- list.files(path = "./data", pattern = ".BAM$", full.names = TRUE)
bam.files

props <- propmapped(files=bam.files)
props

####################
# QUALITY CONTROL  #
####################

# Extract quality scores
qs <- qualityScores(filename="data/SRR1552450.fastq.gz",nreads=100)
# Check dimension of qs
dim(qs)
# Check first few elements of qs with head
head(qs)
# the overall distribution of quality scores across the 100 reads,
boxplot(qs)

######################
#     CHALLENGE      #
######################

#Extract quality scores for SRR1552451.fastq.gz for 50 reads
qs50 <- qualityScores(filename="data/SRR1552451.fastq.gz",nreads=50)
#Plot a boxplot of the quality scores for SRR1552451.fastq.gz
boxplot(qs50)

############
# COUNTING #
############

#count the mapped reads across the mouse genone. (featureCounts contains 
#built-in annotation for mouse (mm9, mm10) and human (hg19) genome assemblies 
#(NCBI refseq annotation).)
fc <- featureCounts(bam.files, annot.inbuilt="mm10")

# See what slots are stored in fc
names(fc)

## Take a look at the featurecounts stats
fc$stat

## Take a look at the dimensions to see the number of genes
dim(fc$counts)
## Take a look at the first 6 lines
head(fc$counts)

################
#  CHALLENGE  ##
################

#Redo the counting over the exons, rather than the genes
exon <- featureCounts(bam.files, annot.inbuilt="mm10", useMetaFeatures = FALSE)
#Check the dimension of the counts slot to see how much larger it is.
dim(exon$counts)

# Using your “.multi.bam” files, redo the counting over genes, allowing for multimapping reads 
# (specify countMultiMappingReads = TRUE), calling the object fc.multi. Check the stats.


