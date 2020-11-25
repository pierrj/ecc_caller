import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster
import csv
import sys

## USAGE ##
# this python script uses hierarchical clustering to merge highly similar eccDNA forming regions together
# usually this means that the start and end coordinates of the eccDNA forming regions are within ~10 bp of each other
# options
# "outputname" - sample name/output prefix
# "scaffold_number" - number of chroms/scaffolds in input bed file
# "max_d" - max depth for defining clusters of eccDNA forming regions, number can be decreased for more conservative clustering or increased for less conservative

outputname = str(sys.argv[1])
scaffold_number = int(sys.argv[2])
max_d = int(sys.argv[3])

# function calculates distance of each eccDNA from the mean eccDNA of the cluster
def get_distance_from_point(row):
    start_distance = abs(row['start'] - start_mean)
    end_distance = abs(row['end'] - end_mean)
    return start_distance + end_distance

# read in input bed file (currently must be named "parallel.confirmed" with scaffold names as 0 indexed numbers)
parallel_confirmed = pd.read_csv("parallel.confirmed", sep = '\t', names = ['chrom', 'start', 'end'])

# merge all entries that are indentical and count their number, these are the starting number of split reads per eccDNA forming region
parallel_confirmed = parallel_confirmed.groupby(parallel_confirmed.columns.tolist()).size().reset_index().\
    rename(columns={0:'splitreads'})

# sort entries by scaffold and make individual dataframes per scaffold to speed things up
parallel_confirmed_dict = {}
for key in range(scaffold_number):
    parallel_confirmed_dict[key] = parallel_confirmed.loc[parallel_confirmed['chrom'] == key]

representative_eccs= []
representative_variants = {}
for chrom in range(scaffold_number): # this is much faster than looking through all of the entries for all scaffolds at once
    scaffold_subset = parallel_confirmed_dict[chrom]
    scaffold_subset_startend =  scaffold_subset[["start", "end"]]
    if scaffold_subset_startend.empty: # some scaffolds don't have any entries
        continue
    else:
        linkage_for_scaffold = linkage(scaffold_subset_startend) # hierarchical clustering
    clusters = fcluster(linkage_for_scaffold, max_d, criterion='distance') # make clusters based off input distance cutoff
    for i in range(1, np.amax(clusters)+1): # loop through each cluster
        boolean_array = clusters==i
        cluster_subset = scaffold_subset[boolean_array]
        if len(cluster_subset) == 1: # many clusters will have just one entry
            clustered = 'no' # useful for detailed look at which entries were merged later
            rep_start = cluster_subset.iloc[0]['start']
            rep_end = cluster_subset.iloc[0]['end']
            split_read_count = cluster_subset.iloc[0]['splitreads']
            point_of_interest = [chrom, rep_start, rep_end, split_read_count ,clustered]
            representative_eccs.append(point_of_interest)
        else:
            clustered = 'yes' # useful for detailed look at which entries were merged later
            start_mean = cluster_subset.mean(axis=0)['start'] # get mean eccDNA, this is the center of the cluster
            end_mean = cluster_subset.mean(axis=0)['end']
            cluster_subset['distance_from_center'] = cluster_subset.apply(get_distance_from_point, axis=1) # get distance from center for all entries
            rep_start = int(cluster_subset[cluster_subset.distance_from_center == cluster_subset.distance_from_center.min()].iloc[0]['start'])
            rep_end = int(cluster_subset[cluster_subset.distance_from_center == cluster_subset.distance_from_center.min()].iloc[0]['end']) # get start and end coordinate of entry closest to center
            split_read_count = cluster_subset['splitreads'].sum() # sum split reads in the cluster
            point_of_interest = [chrom, rep_start, rep_end, split_read_count ,clustered]
            representative_eccs.append(point_of_interest)
            cluster_subset = cluster_subset[["chrom","start", "end", "splitreads"]]
            cluster_subset_listoflists = cluster_subset.values.tolist()
            cluster_subset_tupleoftuples = tuple(tuple(sub) for sub in cluster_subset_listoflists)
            representative_variants[tuple(point_of_interest)] = cluster_subset_tupleoftuples # store all entries in the cluster under the name of the chosen representative entry for detailed look later

# write output to "merged.confirmed", this is the name of the file that coverage_confirm_nodb.py uses
with open('merged.confirmed', 'w', newline = '') as merged:
    w = csv.writer(merged, delimiter = '\t')
    w.writerows(representative_eccs)

# write tsv where one column is the chosen variant for each cluster and the second column is all of the entries that were merged into that cluster
with open('ecccaller_splitreads.' + outputname + '.tsv', 'w', newline="") as variants_dict:
    w = csv.writer(variants_dict, delimiter = '\t')
    for key, value in representative_variants.items():
        w.writerow([key, value])