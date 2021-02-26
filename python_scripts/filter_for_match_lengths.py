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
import re
import math

## USAGE ##
# this python script filters candidate split reads so that only pairs of split reads that have matching match lengths are kept

input_file = str(sys.argv[1])
output_file = str(sys.argv[2])

with open(input_file, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    with open(output_file, 'w', newline = '') as confirmed:
        w = csv.writer(confirmed, delimiter = '\t')
        for row1 in file_reader:
            row2 = next(file_reader)
            matches_row1 = re.findall(r'(\d+)([A-Z]{1})', row1[5])
            matches_sums_row1 = {'M': 0, 'other': 0}
            for i in range(len(matches_row1)):
                if matches_row1[i][1] == 'M':
                    matches_sums_row1['M'] += int(matches_row1[i][0])
                else:
                    matches_sums_row1['other'] += int(matches_row1[i][0])
            matches_row2 = re.findall(r'(\d+)([A-Z]{1})', row2[5])
            matches_sums_row2 = {'M': 0, 'other': 0}
            for i in range(len(matches_sums_row2)):
                if matches_row2[i][1] == 'M':
                    matches_sums_row2['M'] += int(matches_row2[i][0])
                else:
                    matches_sums_row2['other'] += int(matches_row2[i][0])
            read_length_matches = matches_sums_row1['M'] + matches_sums_row2['M']
            read_length_nonmatches = matches_sums_row1['other'] + matches_sums_row2['other']
            if math.isclose(read_length_matches,read_length_nonmatches, abs_tol=10):
                w.writerow(row1)
                w.writerow(row2)