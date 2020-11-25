import sys
from collections import Counter
from itertools import compress
import pandas as pd
import numpy as np
import csv

## USAGE ##
# this script uses information from splitreads (putative eccDNA forming regions) and opposite facing read pairs to confirm eccDNA forming regions
# options:
# "split_read_file" - splitread bed file (putative eccDNA forming regions) where scaffold names are  scaffold numbers
# "outwardfacing_read_file" - outward facing read pairs file where scaffold names are scaffold numbers
# "output_name" - sample name/output prefix
# "scaffold_number" - number of chroms/scaffolds in input bed file
# "output_number" - used to keep track of GNU parallel chunks

split_read_file = str(sys.argv[1])
outwardfacing_read_file = str(sys.argv[2])
output_name = str(sys.argv[3])
scaffold_number = int(sys.argv[4])
output_number = str(sys.argv[5])

# open putative ecc list and index to speed up confirming eccs
with open(split_read_file, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    eccloc_list = []
    for row in file_reader:
        ecc_loc = [int(row[0]) - 1, int(row[1]), int(row[2])]
        eccloc_list.append(ecc_loc)

# open opposite facing discordant read file as a dictionary of numpy arrays for each scaffold
outward_facing = pd.read_csv(outwardfacing_read_file, sep='\t', usecols = [0,1,2,3,4], names = ['chrom', 'start', 'end', 'end_inner', 'start_inner'])
numpy_dict = {}
for key in range(scaffold_number):
    numpy_dict[key] = outward_facing.loc[outward_facing['chrom'] == key + 1].to_numpy()

def confirmeccs(ecc): # using numpy arrays and only looking at opposite facing read pairs within the same scaffold makes this very fast
    mask1 = numpy_dict[ecc[0]][:, 1] >= ecc[1]
    mask2 = numpy_dict[ecc[0]][:, 2] <= ecc[2] 
    mask_total = np.logical_and(mask1, mask2)
    masked = numpy_dict[ecc[0]][mask_total, :] # filter to opposite read pairs within the putative eccDNA forming region
    masked[:,3] -= ecc[1]
    masked[:,4] = ecc[2]-masked[:,4]
    final =  masked[masked[:,3] + masked[:, 4] <= 500, :] # does proximity filtering based off an estimated insert size of 400 + 25%, this may need to be tweaked per library
    if final.size == 0:
        return False
    return True

# get true/false list if each ecc is confirmed, then compress only keeps where true is in the list
yesornoeccs = list(map(confirmeccs, eccloc_list))
confirmedeccs = list(compress(eccloc_list, yesornoeccs))

# write confirmed eccs to file, number is here to put parallelized chunks back together
with open('parallel.confirmed'+output_number, 'w', newline = '') as confirmed:
    w = csv.writer(confirmed, delimiter = '\t')
    w.writerows(confirmedeccs)