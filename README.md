<!-- ecc_caller -->
## ecc_caller

ecc_caller is a pipeline designed to take Illumina NovaSeq paired end sequencing reads generated from isolated circular DNA and call eccDNA forming regions.

The pipeline uses information from two structural read variants to call eccDNA forming regions: split reads and opposite facing read pairs.

It also assigns a confidence score based off the number of split reads corresponding to each region as well as its coverage.

## Requirements

Conda installation instructions coming soon!

ecc_caller requires the following software (more details and versions coming soon). Make sure these are added to PATH:
* cutadapt
* BWA MEM
* samtools
* bedtools
* GNU parallel
* Picard tools

ecc_caller contains Python scripts which were written in Python 3.8 and requires the following modules:
* sys
* collections
* itertools
* pandas
* numpy
* csv
* scipy
* statistics
* subprocess

<!-- USAGE EXAMPLES -->
## Usage

Once all software is installed with paths set and the git repository has been cloned make sure to set this path variable so the wrapper scripts know where to find the python scripts and your picard.jar
   ```sh
   export ECC_CALLER_PYTHON_SCRIPTS=/path/to/git/repo/python_scripts/
   export ECC_CALLER_PICARD=/path/to/picard/
   ```

ecc_caller uses a mapfile to work with any input genome. This mapfile should be a list of the first field of all fasta entries (scaffolds/chromosomes) of interest for the analysis. This can be all scaffolds in the original genome file or only scaffolds of interest (i.e. excluding mitochondria).

Example command to create mapfile:
   ```sh
   grep '>' guy11_genome_baoetal2017_with_70-15_mito.fasta | grep -v mito | awk '{print substr($1,2)}' > mapfile
   ```

Before mapping, remember to create a bwa index of your genome:
   ```sh
   bwa index -p genome_bwa genome.fasta 
   ```

Then remove adapters and map sequencing reads
   ```sh
   generate_bam_file.sh -g genome_bwa -1 R1.fastq -2 R2.fastq -s output_name -t n_threads -m mapfile
   ```
   
Next call eccDNA forming regions
   ```sh
   call_ecc_regions.sh -m mapfile -s output_name -t n_threads -b uniq.filtered.sorted.output_name.bam -q multimapped.filtered.name_sorted.output_name.bam
   ```
   
Finally, assign confidence to each eccDNA forming region
   ```sh
   assign_confidence.sh -m mapfile -s output_name -t n_threads -b no_secondary.filtered.sorted.output_name.bam -r output_name.confirmedsplitreads.bed
   ```
   
## Outputs

ecc_caller currently outputs many useful files for analysis and QC. File names need to be cleaned up a bit, update coming soon.

A few important files to note:
* filtered.sorted.output_name.bam - output from bwa mem, sorted, and filtered to scaffolds of interest
* lengthfiltered.merged.splitreads.output_name.renamed.bed - bed file containing split reads to be used for ecc calling
* outwardfacing.output_name.bed - bed file containing all outward facing read pairs to be used for ecc calling
* outputname.confirmedsplitreads.bed - output from call_ecc_regions.sh, split reads confirmed by the presence of an outward facing read pair
* ecccaller_splitreads.outputname.tsv - output from cluster_eccs.py, some clustering is preformed for highly similar eccDNA forming regions, this file shows which representative eccDNA forming region was called for each cluster of eccDNA forming regions
* ecccaller_output.output_name.renamed.details.tsv - tsv containing location of eccDNA forming regions as well as their confidence score, reason for confidence score and depth of coverage (see below)
* ecccaller_output.output_name.renamed.bed - bed file of all eccDNA forming regions, colored by confidence score (see below)

## Confidence Scores for eccDNA Forming Regions

coverage_confirm_nodb.py assigns confidence scores to each eccDNA region. These criteria can currently only be modified by changing the python script but I will eventually include this as a command line argument. These are as follows:
* lowq - low quality, greater than 5% of the region is uncovered by mapped reads and/or only one split read is associated with this region (colored red in bed file output)
* conf - medium confidence, 2 split reads associated with this region, coverage of the region is less than double that of the flanking regions of equal sizes (colored yellow in bed file output)
* hconf - high confidence, 2 split reads and coverage of the region is more than double that of the flanking regions of equal sizes OR 3 or more split reads associated with this region (colored green in bed file output)

## Acknowledgements

Initial framework for the pipeline was inspired by the methods described in the following paper but all code in this repository is original:

Møller, H.D., Mohiyuddin, M., Prada-Luengo, I. et al. Circular DNA elements of chromosomal origin are common in healthy human somatic tissue. Nat Commun 9, 1069 (2018). https://doi.org/10.1038/s41467-018-03369-8

Thank you to Ksenia Krasileva and all members of the KVK lab for their feedback on progress on this pipeline.

<!-- LICENSE -->
## License

ecc_caller is freely available under the MIT license

<!-- CONTACT -->
## Contact

Please feel free to contact me directly with any questions

Pierre Joubert - [@pmjoubert](https://twitter.com/pmjoubert) - pierrj@berkeley.edu

Project Link: [https://github.com/pierrj/ecc_caller](https://github.com/pierrj/ecc_caller)
