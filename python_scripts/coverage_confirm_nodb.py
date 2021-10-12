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
import sys
import csv
import numpy as np
import statistics
import subprocess

## USAGE ##
# this python script uses read coverage information as well as split read count to assign confidence to eccDNA regions
# lowq eccDNA regions either have read coverage less than 95% or only one split read
# conf eccDNA regions have more than two split reads but are not more covered than neighboring regions
# hconf eccDNA regions have more than three split reads OR two split reads and are more covered than neighboring regions
# options:
# "output_name" - sample name/output prefix
# "output_number" - used to keep track of GNU parallel chunks
# "bam_file" - input bamfile of mapped reads

output_name = str(sys.argv[1])
output_number = str(sys.argv[2])
bam_file = str(sys.argv[3])

# read in output from splitread calls, input must be named merged.confirmed with chunk number for parallelization
with open('merged.confirmed'+output_number) as merged:
    merged_reader = csv.reader(merged, delimiter = '\t')
    flat_merged_list = [[int(row[0]), int(row[1]), int(row[2]), int(row[3])] for row in merged_reader]

# function that does the confidence assignments
def confidence_check(ecc):
    # get coverage of confirmed ecc region
    region_len = ecc[2] - ecc[1]
    beforestart = ecc[1] - region_len # get coordinates of regions before and after
    afterstart =  ecc[2] + 2 + region_len
    if beforestart > 0: # deal with what happens if the region before includes the end of the scaffold
        region = [ecc[0], beforestart, afterstart]
    else:
        region = [ecc[0], ecc[1], afterstart]
    samtools_get_region = ['samtools', 'depth', '-a', '-d 8000', '-r', str(region[0]+1)+':'+str(region[1])+'-'+str(region[2]), bam_file] # get region coverage
    sp = subprocess.Popen(samtools_get_region, shell= False, stdout=subprocess.PIPE)
    region_cov = [int(line.decode("utf-8").strip().split("\t")[2]) for line in sp.stdout] # translate samtools output to list
    if beforestart > 0: # deal with what happens if the region before includes the end of the scaffold
        ecc_region_cov = region_cov[region_len+1:((2 * region_len+1)+1)]
        before_region_cov = region_cov[0:region_len+1]
        after_region_cov = region_cov[(2 * region_len)+2:]
    else:
        ecc_region_cov = region_cov[0:region_len+1]
        after_region_cov = region_cov[region_len+1:((2 * region_len+1)+1)]
    mean_region = round(statistics.mean(ecc_region_cov), 2) # get mean coverage for target region, as well as regions before and after
    if beforestart > 0:
        mean_before = round(statistics.mean(before_region_cov), 2)
    else:
        mean_before = 'N/A'
    if len(after_region_cov) == region_len + 1:
        mean_after = round(statistics.mean(after_region_cov), 2)
    else:
        mean_after = 'N/A'
    # write coverage string containing the means of the regions before, within, and after the confirmed ecc
    coverage_string = str(mean_before)+';'+str(mean_region)+';'+str(mean_after)
    # if less than 95% of the ecc region is covered than the ecc is low confidence
    if ecc_region_cov.count(0) / len(ecc_region_cov) > 0.05:
        ecc.append('lowq')
        ecc.append('incomplete_coverage')
        ecc.append(coverage_string)
        return ecc
    if ecc[3] == 1:
        ecc.append('lowq')
        ecc.append('only_one_splitread')
        ecc.append(coverage_string)
        return ecc
    # if the ecc before or after regions fall beyond the borders of the chromosome than the ecc is medium confidence if it only has support from 2 split reads
    if beforestart <= 0 and ecc[3] < 3:
        ecc.append('conf')
        ecc.append('too_close_to_start')
        ecc.append(coverage_string)
        return ecc
    if len(after_region_cov) != region_len + 1 and ecc[3] <3:
        ecc.append('conf')
        ecc.append('too_close_to_end')
        ecc.append(coverage_string)
        return ecc
    if beforestart <= 0 and ecc[3] >= 3:
        ecc.append('hconf')
        ecc.append('splitreads')
        ecc.append(coverage_string)
        return ecc
    if len(after_region_cov) != region_len + 1 and ecc[3] >= 3:
        ecc.append('hconf')
        ecc.append('splitreads')
        ecc.append(coverage_string)
        return ecc
    # if the mean coverage of the eccDNA is twice the size of the before and after regions OR the eccDNA is supported by three or more splitreads then eccDNA is high confidence
    if mean_region >= (2*mean_before) and mean_region >= (2*mean_after) and ecc[3] < 3:
        ecc.append('hconf')
        ecc.append('coverage')
        ecc.append(coverage_string)
        return ecc
    if (mean_region < (2*mean_before) or mean_region < (2*mean_after)) and ecc[3] >= 3:
        ecc.append('hconf')
        ecc.append('splitreads')
        ecc.append(coverage_string)
        return ecc
    if mean_region >= (2*mean_before) and mean_region >= (2*mean_after) and ecc[3] >= 3:
        ecc.append('hconf')
        ecc.append('splitreads_and_coverage')
        ecc.append(coverage_string)
    else:
        ecc.append('conf')
        ecc.append('coverage_and_splitreads_too_low')
        ecc.append(coverage_string)
    return ecc

confidence_flat_merged_list = list(map(confidence_check, flat_merged_list))

# write tsv which includes each eccDNA location as well as their confidence and supporting number of split reads
with open('ecccaller_output.' + output_name + '.details.tsv'+output_number, 'w', newline = '') as bed:
    w = csv.writer(bed, delimiter = '\t')
    for i in range(len(confidence_flat_merged_list)):
        row = [confidence_flat_merged_list[i][0]+1, confidence_flat_merged_list[i][1], confidence_flat_merged_list[i][2], confidence_flat_merged_list[i][3], confidence_flat_merged_list[i][4],confidence_flat_merged_list[i][5], confidence_flat_merged_list[i][6]]
        w.writerow(row)

# write bed version of tsv with details missing and color coded by confidence level
with open('ecccaller_output.' + output_name + '.bed'+output_number, 'w', newline = '') as bed:
    w = csv.writer(bed, delimiter = '\t')
    for i in range(len(confidence_flat_merged_list)):
        if confidence_flat_merged_list[i][4] == 'lowq':
            color = '255,0,0'
        if confidence_flat_merged_list[i][4] == 'conf':
            color = '255,255,0'
        if confidence_flat_merged_list[i][4] == 'hconf':
            color = '0,255,0'
        row = [confidence_flat_merged_list[i][0]+1, confidence_flat_merged_list[i][1], confidence_flat_merged_list[i][2], i, 0, '+', confidence_flat_merged_list[i][1], confidence_flat_merged_list[i][2], color ]
        w.writerow(row)

# write bed version of all splitreads from conf and hconf eccdnas
with open('ecccaller_output.splitreads.' + output_name + '.bed'+output_number, 'w', newline = '') as bed:
    w = csv.writer(bed, delimiter = '\t')
    for i in range(len(confidence_flat_merged_list)):
        row = [confidence_flat_merged_list[i][0]+1, confidence_flat_merged_list[i][1], confidence_flat_merged_list[i][2]]
        if confidence_flat_merged_list[i][4] != "lowq":
            for k in range(confidence_flat_merged_list[i][3]):
                w.writerow(row)