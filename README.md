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

Supplementary
- In silico test of the UCE Coleoptera probes on the Tetraopes genomes
- 

## Quality control

```
# Fastp 
# Based on Matthew Van Dam's script

# Create a clean dir for your clean reads to go
dir.create("clean") #do once
# Make list of file names and paste together commands for fastp as above example 
f1 = list.files(getwd(), ".*_1.fq.gz", recursive=TRUE, full.names=TRUE)
f2 = list.files(getwd(), ".*_2.fq.gz", recursive=TRUE, full.names=TRUE)
f1=f1[-21] # this removes any files that may match, change number to remove right one if need be
f2=f2[-21] # this removes any files that may match, change number to remove right one if need be

# Need to match the regular expression of your file names e.g. /home/ngutierrez/array1/Tetraopini/0.untrimmed/Nayeli_Pool1_
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

# Not in R, remove all the unpaied r1 r2 files keeping only the catted ones
find . -type f -name '*.unpaired.fq.gz' -exec rm {} \;

```

### Genome assembly

```
# Spades
# Based on Matthew Van Dam's script

# Define c1, c2, and catted again in 0.trimmed directory
f1 = list.files(getwd(), ".*_1.fq.gz", recursive=TRUE, full.names=TRUE)
f2 = list.files(getwd(), ".*_2.fq.gz", recursive=TRUE, full.names=TRUE)
f1=f1[-21] # this removes any files that may match, change number to remove right one if need be
f2=f2[-21] # this removes any files that may match, change number to remove right one if need be

c1 = gsub("ed/", "ed/clean/", f1 )
c2 = gsub("ed/", "ed/clean/", f2 )

# Removes the funky names parts
outs = gsub("_CKDL2.*","_spades_vanilla", catted)

# Make Spades commands for denovo assembly then paste them into two different screens and or CPUs
cmd = paste("/array1/ngutierrez/software/anaconda3/envs/r_env/bin/spades.py -k 21,33,55,77,99,127 --pe1-1 ", c1, " --pe1-2 ", c2, " --pe1-s ", catted, " -o ", outs, " -m 700 -t 16")
cmd = paste("spades.py -k 21,33,55,77,99,127 --pe1-1 ", c1, " --pe1-2 ", c2, " --pe1-s ", catted, " -o ", outs, " -m 700 -t 16")
cmd = paste("spades.py -k 21,33,55,77,99,127 --pe1-1 ", c1, " --pe1-2 ", c2, " --pe1-s ", catted, " -o ", outs, " -m 700 -t 16")

# Paste the lines into a screen session e.g. screen -S some_name , to re-enter screen -r some_name
cmd = cmd[1:10] # for one cpu
cmd = cmd[11:20] # for the other cpu ... until you have all the assemblies done
spades.py -k 21,33,55,77,99,127 --pe1-1  /home3/mvandam/Xerces/clean/EL10_CKDL210001104-1a_HFLHVCCX2_L4_1.fq.gz  --pe1-2  /home3/mvandam/Xerces/clean/EL10_CKDL210001104-1a_HFLHVCCX2_L4_2.fq.gz  --pe1-s  /home3/mvandam/Xerces/clean/EL10_CKDL210001104-1a_HFLHVCCX2_L4_unpaired.fq.gz  -o  /home3/mvandam/Xerces/clean/EL10_spades_vanilla  -m 700 -t 36

spades.py -k 21,33,55,77,99,127 --pe1-1  /home/ngutierrez/array1/Tetraopini/0.untrimmed/clean/Nayeli_Pool1_1.fq.gz  --pe1-2  /home/ngutierrez/array1/Tetraopini/0.untrimmed/clean/Nayeli_Pool1_2.fq.gz  --pe1-s  /home/ngutierrez/array1/Tetraopini/0.untrimmed/clean/Nayeli_Pool1_unpaired.fq.gz  -o  /home/ngutierrez/array1/Tetraopini/0.untrimmed/clean/Nayeli_Pool1_spades_vanilla  -m 700 -t 36

```

### Assembly evaluation with QUAST

```
for f in *.fasta; do quast --eukaryote "$f" -o "$f".quast; done

```
## 

```
# Phyluce 1.7.1
# Count of matches that we recovered to UCE loci in the probe set, and extract all of the “good” loci to a monolithic FASTA 
phyluce_assembly_get_match_counts \
    --locus-db /home/ngutierrez/array1/Tetraopini/tutorial4_allspp/in-silico-lastz/probe.matches.sqlite \
    --samples 39 \
    --taxon-list-config in-silico-coleoptera-taxon-sets.conf \
    --taxon-group 'nohtesta' \
    --output taxon-sets/no-htesta/insilico-incomplete.conf \
    --log-path log \
    --incomplete-matrix

# Change to the taxon-sets/all directory
cd taxon-sets/no-htesta
# Make a log directory to hold our log files - this keeps things neat
mkdir log

# Extract the FASTA information for each locus into a monolithic FASTA file:
phyluce_assembly_get_fastas_from_match_counts \
    --contigs ../../coleoptera-genome-fasta \
    --locus-db ../../in-silico-lastz/probe.matches.sqlite \
    --match-count-output insilico-incomplete.conf \
    --output insilico-incomplete.fasta \
    --incomplete-matrix insilico-incomplete.incomplete \
    --log-path log

# Align the conserved locus data
phyluce_align_seqcap_align \
    --input insilico-incomplete.fasta \
    --output mafft \
    --taxa 39 \
    --incomplete-matrix \
    --cores 18 \
    --no-trim \
    --output-format fasta \
    --log-path log

# Trim the conserved locus alignments
phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
    --alignments mafft \
    --output mafft-gblocks \
    --b1 0.5 \
    --b4 8 \
    --cores 18 \
    --log log

# Remove the locus names from each alignment
phyluce_align_remove_locus_name_from_files \
    --alignments mafft-gblocks \
    --output mafft-gblocks-clean \
    --cores 18 \
    --log-path log

# Get stats across the aligned loci
phyluce_align_get_align_summary_data \
    --alignments mafft-gblocks-clean \
    --cores 28 \
    --log-path log

# Get uce mean mength stats (taken from tutorial II)
# Exploding the monolithic FASTA file
cd /home/ngutierrez/array1/Tetraopini/tutorial4_allspp/taxon-sets/no-htesta

# By taxon
phyluce_assembly_explode_get_fastas_file \
    --input insilico-incomplete.fasta \
    --output exploded-fastas \
    --by-taxon

# Then
for i in exploded-fastas/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv;
done

# Generate an incomplete matrix
# 50 
phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-gblocks-clean \
    --taxa 39 \
    --output mafft-gblocks-50p \
    --percent 0.50 \
    --cores 18 \
    --log log
# 75
phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-gblocks-clean \
    --taxa 39 \
    --output mafft-gblocks-70p \
    --percent 0.75 \
    --cores 18 \
    --log log
# 90 
phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-gblocks-clean \
    --taxa 39 \
    --output mafft-gblocks-90p \
    --percent 0.90 \
    --cores 18 \
    --log log

# Getting input data for pre-partitioning
# Need to get the alignment in nexus format
# 50 
phyluce_align_concatenate_alignments \
    --alignments mafft-gblocks-50p \
    --output mafft-gblocks-50p-raxml-nexus \
    --log-path log --nexus
# 75 
phyluce_align_concatenate_alignments \
    --alignments mafft-gblocks-75p \
    --output mafft-gblocks-75p-raxml-nexus \
    --log-path log --nexus
# 90 
phyluce_align_concatenate_alignments \
    --alignments mafft-gblocks-90p \
    --output mafft-gblocks-90p-raxml-nexus \
    --log-path log --nexus

# Edit the file 
cat mafft-gblocks-50p-raxml-nexus.nexus | sed 's/uce-/uce_/g' > mafft-gblocks-50p-raxml-nexus-edited.nexus
cat mafft-gblocks-65p-raxml-nexus.nexus | sed 's/uce-/uce_/g' > mafft-gblocks-75p-raxml-nexus-edited.nexus
cat mafft-gblocks-90p-raxml-nexus.nexus | sed 's/uce-/uce_/g' > mafft-gblocks-90p-raxml-nexus-edited.nexus

# Getting input data for raxml
# Setup the PHYLIP-formatted files for raxm
phyluce_align_concatenate_alignments \
    --alignments mafft-gblocks-50p \
    --output mafft-gblocks-50p-raxml \
    --log-path log --phylip
```
