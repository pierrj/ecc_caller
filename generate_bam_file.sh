#!/bin/bash
#MIT License
#
#Copyright (c) 2020 Pierre Michel Joubert
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
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
cutadapt -j ${THREADS} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --nextseq-trim=20 -o tmp.trimmed.seqprep.${SAMPLE}_R1.fastq -p tmp.trimmed.seqprep.${SAMPLE}_R2.fastq ${READONE} ${READTWO}

# map reads
bwa mem -q -a -t ${THREADS} ${GENOME_DB} tmp.trimmed.seqprep.${SAMPLE}_R1.fastq tmp.trimmed.seqprep.${SAMPLE}_R2.fastq -o tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam

# generate bam file from sam, this should probably just be passed as an option to bwa mem
samtools view -S -b tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam > ${SAMPLE}.mergedandpe.bwamem.bam

# sort and index for filtering
samtools sort -@ ${THREADS} ${SAMPLE}.mergedandpe.bwamem.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.bam
samtools index ${SAMPLE}.sorted.mergedandpe.bwamem.bam

# filter using names of contigs in mapfile
samtools view -b ${SAMPLE}.sorted.mergedandpe.bwamem.bam $(cat ${MAPFILE} | tr "\n" " ") > filtered.sorted.${SAMPLE}.bam

# filter to only primary alignments
samtools view -b -F 256 filtered.sorted.${SAMPLE}.bam > no_secondary.filtered.sorted.${SAMPLE}.bam

# filter to only uniquely mapped alignments
samtools view -b -q 1 no_secondary.filtered.sorted.${SAMPLE}.bam > uniq.filtered.sorted.${SAMPLE}.bam

# sort by read name for subsampling reads later
samtools sort -n -@ ${THREADS} filtered.sorted.${SAMPLE}.bam > multimapped.filtered.name_sorted.${SAMPLE}.bam


# remove tmp files and unfiltered files
# tmp command should probably be more specific, or use a tmp directory
rm tmp*
rm ${SAMPLE}.mergedandpe.bwamem.bam
rm filtered.sorted.${SAMPLE}.bam
rm ${SAMPLE}.sorted.mergedandpe.bwamem.bam