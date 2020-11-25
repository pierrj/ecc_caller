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

Once all software is installed with paths set and the git repository has been cloned make sure to set this path variable so the wrapper scripts know where to find the python scripts 
   ```sh
   export ECC_CALLER_PYTHON_SCRIPTS=/path/to/git/repo/python_scripts
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
   call_ecc_regions.sh -m mapfile -s output_name -t n_threads -b filtered.sorted.output_name.bam
   ```
   
Finally, assign confidence to each eccDNA forming region
   ```sh
   assign_confidence_nodb.sh -m mapfile -s output_name -t n_threads -b filtered.sorted.output_name.bam -r output_name.confirmedsplitreads.bed
   ```

<!-- LICENSE -->
## License

ecc_caller is freely available under the MIT license

<!-- CONTACT -->
## Contact

Please feel free to contact me directly with any questions

Pierre Joubert - [@pmjoubert](https://twitter.com/pmjoubert) - pierrj@berkeley.edu

Project Link: [https://github.com/pierrj/ecc_caller](https://github.com/pierrj/ecc_caller)
