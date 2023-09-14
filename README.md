# 2023_MacroevolTetraopes

Custom scripts were employed for the analysis of the Macroevolution of Milkweed Longhorn Beetles Project.

We did not develop new software for this study; hence, we offer the scripts for conducting analyses and provide references to other sources utilized.

For each section, we provide a folder with a Jupyter Notebook with a detailed explanation of the analysis, and the scripts used are provided. Additionally, whenever applicable, we provide the output files and reports from the different software employed:

## Index
- Quality control
- Genome assembly
- Ultraconserved Element (UCE) custom probe set design
- UCE matrix generation and partitioning
- Phylogenomics
- Divergence Time Estimation
- Ancestral range estimation

### Quality control

```
#### Fastp ####
### Script based on Matthew Van Dam's script

## create a clean dir for your clean reads to go
dir.create("clean") #do once
#make list of file names and paste together commands for fastp as above example 
f1 = list.files(getwd(), ".*_1.fq.gz", recursive=TRUE, full.names=TRUE)
f2 = list.files(getwd(), ".*_2.fq.gz", recursive=TRUE, full.names=TRUE)
f1=f1[-21] # this removes any files that may match, change number to remove right one if need be
f2=f2[-21] # this removes any files that may match, change number to remove right one if need be

# need to match the regular expression of your file names e.g. /home/ngutierrez/array1/Tetraopini/0.untrimmed/Nayeli_Pool1_CKDL220004503-1a-7UDI293-AK900_HGC7NDSX3_L4_1.fq.gz

c1 = gsub("ed/", "ed/clean/", f1 )
c2 = gsub("ed/", "ed/clean/", f2 )

# Fastp

c1u = gsub(".fq.gz", ".unpaired.fq.gz", c1)
c2u = gsub(".fq.gz", ".unpaired.fq.gz", c2)
htmlrep = gsub(".unpaired.fq.gz",".html",c1u)
cmd = paste("fastp --in1 ", f1, " --in2 ", f2, " --out1 ", c1," --out2 ", c2, " --unpaired1 ", c1u, " --unpaired2 ", c2u, " --thread 16 --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --html ", htmlrep, sep="") # --detect_adapter_for_pe

mclapply(cmd, system, mc.cores=getOption("mc.cores", 4)) # do only once !!!!!!!!!!!!, this runs 16 threads X 4 = 64 total threads

cat1 = c1u
cat2 = c2u
catted = sub("_1.unpaired.fq.gz", "_unpaired.fq.gz", c1u)
cmd = paste("cat", cat1, cat2, ">>", catted)
lapply(cmd, system) #run only once !!!!!!!!!!!

## not in R, remove all the unpaied r1 r2 files keeping only the catted ones
### you will want to check find . -type f -name '*.unpaired.fq.gz' before running the line below as this can delete many files be very careful, check to make sure you are in the correct directory
#find . -type f -name '*.unpaired.fq.gz' -exec rm {} \;

```
