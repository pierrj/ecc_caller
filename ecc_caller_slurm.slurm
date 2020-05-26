#!/bin/bash
#SBATCH --job-name=doesnt_matter
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
cd $SLURM_SUBMIT_DIR
cd ${sample}
samtools view -b -q 1 sorted.mergedandpe.${sample}_bwamem.bam $(cat /global/scratch/users/pierrj/references/guy11_genome_baoetal2017.contignames | tr "\n" " ") > guy11.sorted.mergedandpe.${sample}_bwamem.bam
samtools view -f 81 -F 4 guy11.sorted.mergedandpe.${sample}_bwamem.bam > tmp.reverseread1.mergedandpe.${sample}_bwamem.sam
samtools view -f 145 -F 4 guy11.sorted.mergedandpe.${sample}_bwamem.bam > tmp.reverseread2.mergedandpe.${sample}_bwamem.sam
samtools view -f 65 -F 20 guy11.sorted.mergedandpe.${sample}_bwamem.bam > tmp.forwardread1.mergedandpe.${sample}_bwamem.sam
samtools view -f 129 -F 20 guy11.sorted.mergedandpe.${sample}_bwamem.bam > tmp.forwardread2.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.reverseread1.mergedandpe.${sample}_bwamem.sam tmp.reverseread1.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.reverseread2.mergedandpe.${sample}_bwamem.sam tmp.reverseread2.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.forwardread1.mergedandpe.${sample}_bwamem.sam tmp.forwardread1.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${sample}_bwamem.sam
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.forwardread2.mergedandpe.${sample}_bwamem.sam tmp.forwardread2.mergedandpe.${sample}_bwamem.sam > tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${sample}_bwamem.sam
cat  tmp.samechromosome.exactlytwice.reverseread1.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.reverseread2.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread1.mergedandpe.${sample}_bwamem.sam tmp.samechromosome.exactlytwice.forwardread2.mergedandpe.${sample}_bwamem.sam > samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.bam

samtools view guy11.sorted.mergedandpe.${sample}_bwamem.bam | awk '{ if (($2 == 83 || $2 == 147 ) && $9 > 0) print $0 ; else if (($2 == 99 || $2 == 163) && $9 <0) print $0}' | cat <(samtools view -H mergedandpe.${sample}_bwamem.bam) - | samtools view -b -h - > tmp.outwardfacing.mergedandpe.${sample}_bwamem.bam
bedtools bamtobed -i tmp.outwardfacing.mergedandpe.${sample}_bwamem.bam | sort -k 4,4 > tmp.outwardfacing.mergedandpe.${sample}_bwamem.bed
mv tmp.outwardfacing.mergedandpe.${sample}_bwamem.bed tmp.outwardfacing.mergedandpe.${sample}_bwamem.bed.old
awk 'BEGIN {OFS="\t"}; {print $1,$2,$3,substr($4, 1, length($4)-2),$5,$6}' tmp.outwardfacing.mergedandpe.${sample}_bwamem.bed.old > tmp.outwardfacing.mergedandpe.${sample}_bwamem.bed.old.trimmed
awk 'NR==FNR{a[$4]++; next} a[$4]==2' tmp.outwardfacing.mergedandpe.${sample}_bwamem.bed.old.trimmed tmp.outwardfacing.mergedandpe.${sample}_bwamem.bed.old.trimmed > outwardfacing.mergedandpe.${sample}_bwamem.bed

bedtools genomecov -d -ibam guy11.sorted.mergedandpe.${sample}_bwamem.bam | awk ' /M/ {print $0}' > genomecoverage.mergedandpe.${sample}_bwamem.bed

rm tmp.*

ipcluster start -n $SLURM_NTASKS --cluster-id="slurm-${SLURM_JOBID}" &
sleep 300
ipython /global/home/users/pierrj/git/python/ecc_caller_python.py samechromosome.exactlytwice.all.mergedandpe.${sample}_bwamem.bam outwardfacing.mergedandpe.${sample}_bwamem.bed genomecoverage.mergedandpe.${sample}_bwamem.bed ${sample} 56
ipcluster stop --cluster-id="slurm-${SLURM_JOBID}"