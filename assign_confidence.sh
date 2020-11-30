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
while getopts m:s:t:b:r: option
do
case "${option}"
in
m) MAPFILE=${OPTARG};;
s) SAMPLE=${OPTARG};;
t) THREADS=${OPTARG};;
b) FILTERED_BAMFILE=${OPTARG};;
r) CONFIRMED_SPLITREADS=${OPTARG};;
esac
done

## USAGE ##
# this script assigns confidence to previously confirmed eccDNA forming regions
# this script mostly serves as a wrapper/preparation for the two python scripts
# see converage_confirm_nodb.py for more details on specific confidence cutoffs
# options
# -m mapfile with names of contigs of interest, as written in the fasta file used to make bwa genome database
# -s sample name/output prefix
# -t threads
# -b bamfile of mapped reads, filtered to contigs of interest
# -r confirmed split read file, output of call_ecc_regions.sh

# get chrom/scaffold count from mapfile
# generate bam file with scaffolds renamed to chrom/scaffold number for compatability
chrom_count=$(wc -l ${MAPFILE} | awk '{print $1}')
for (( i = 1 ; i < ${chrom_count}+1; i++)); do echo $i ; done > tmp.chrom_count
paste tmp.chrom_count ${MAPFILE} > tmp.chrom_count_and_names
samtools view -H ${FILTERED_BAMFILE} | awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{ if ($2 ~ /^SN/ && substr($2, 4) in a) {print $1, "SN:"a[substr($2,4)], $3} else {print $0}}' tmp.chrom_count_and_names - |\
    samtools reheader - ${FILTERED_BAMFILE} > renamed.filtered.sorted.${SAMPLE}.bam

samtools index renamed.filtered.sorted.${SAMPLE}.bam

# regenerate 0 indexed split read bed file with scaffold numbers for compatability
# parallel.confirmed is input for cluster_eccs.py
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_count_and_names ${CONFIRMED_SPLITREADS} > parallel.plusone.confirmed
awk -v OFS='\t' '{print $1-1, $2, $3}' parallel.plusone.confirmed > parallel.confirmed

# uses hierarchical clustering to merge highly similar confirmed eccDNA forming regions together
# see cluster_eccs.py for more info
python ${ECC_CALLER_PYTHON_SCRIPTS}/cluster_eccs.py ${SAMPLE} ${chrom_count} 10

# split confirmed eccDNA forming regions into chunks to be used with GNU parallel
# then assign confidence to eccDNAs based off split read counts and read coverage
# see coverage_confirm_nodb.py for more details
shuf merged.confirmed > shuf.merged.confirmed
split --number=l/${THREADS} --numeric-suffixes=1 shuf.merged.confirmed merged.confirmed
parallel -j ${THREADS} --link python ${ECC_CALLER_PYTHON_SCRIPTS}/coverage_confirm_nodb.py ${SAMPLE} {} renamed.filtered.sorted.${SAMPLE}.bam ::: $(seq -w 1 ${THREADS})

# put parallel chunks back together
cat $(find . -maxdepth 1 -name "ecccaller_output.${SAMPLE}.details.tsv*" | xargs -r ls -1 | tr "\n" " ") > ecccaller_output.${SAMPLE}.details.tsv
cat $(find . -maxdepth 1 -name "ecccaller_output.${SAMPLE}.bed*" | xargs -r ls -1 | tr "\n" " ") > ecccaller_output.${SAMPLE}.bed

# rename output files to original chrom/scaffold names
paste ${MAPFILE} tmp.chrom_count > tmp.chrom_names_and_count
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_names_and_count ecccaller_output.${SAMPLE}.details.tsv > ecccaller_output.${SAMPLE}.renamed.details.tsv
awk -v OFS='\t' 'NR==FNR{a[$2]=$1;next}{$1=a[$1];}1' tmp.chrom_names_and_count ecccaller_output.${SAMPLE}.bed > ecccaller_output.${SAMPLE}.renamed.bed

# clean up tmp files
rm ecccaller_output.${SAMPLE}.details.tsv*
rm ecccaller_output.${SAMPLE}.bed*
rm parallel.confirmed
rm parallel.plusone.confirmed
rm renamed.filtered.sorted.${SAMPLE}.bam
rm renamed.filtered.sorted.${SAMPLE}.bam.bai
rm merged.confirmed*

## should probably sort outputs at the end