#!/bin/bash
while getopts g:1:2:s:t:m: option
do
case "${option}"
in
g) GENOME_DB=${OPTARG};;
1) READONE=${OPTARG};;
2) READTWO=${OPTARG};;
s) SAMPLE=${OPTARG};;
t) THREADS=${OPTARG};;
m) MAPFILE=${OPTARG};;
esac
done

## USAGE ##
# generates a bam file from fastq files
# filters the bam file to contigs of interest
# can be replaced by any filtered bam file
# options:
# -g path to bwa genome database name
# -1 path to read 1 file
# -2 path to read 2 file
# -s sample name/output prefix
# -t threads
# -m mapfile with names of contigs of interest, as written in the fasta file used to make bwa genome database


# cutadapt uses specific adapters here, will be replaced with trimmomatic which predicts adapters automatically eventually
# nextseq-trim option is necessary for nextseq data
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --nextseq-trim=20 -o tmp.trimmed.seqprep.${SAMPLE}_R1.fastq -p tmp.trimmed.seqprep.${SAMPLE}_R2.fastq ${READONE} ${READTWO}

# map reads
bwa mem -t ${THREADS} ${GENOME_DB} tmp.trimmed.seqprep.${SAMPLE}_R1.fastq tmp.trimmed.seqprep.${SAMPLE}_R2.fastq -o tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam

# generate bam file from sam, this should probably just be passed as an option to bwa mem
samtools view -S -b tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam > ${SAMPLE}.mergedandpe.bwamem.bam

# sort and index for filtering
samtools sort ${SAMPLE}.mergedandpe.bwamem.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.bam
samtools index ${SAMPLE}.sorted.mergedandpe.bwamem.bam

# filter using names of contigs in mapfile
samtools view -b ${SAMPLE}.sorted.mergedandpe.bwamem.bam $(cat ${MAPFILE} | tr "\n" " ") > filtered.sorted.${SAMPLE}.bam

# remove tmp files and unfiltered files
# tmp command should probably be more specific, or use a tmp directory
rm tmp*
rm ${SAMPLE}.mergedandpe.bwamem.bam