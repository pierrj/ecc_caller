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
import csv
import sys


## USAGE ##
# simple script that fixes chunks so that read names never occur in more than one chunk
# first argument is input chunk
# second argument is chunk number of input chunk

file_name = str(sys.argv[1])
output_number = str(sys.argv[2])

with open('split_line_fix.'+output_number, 'w', newline = '') as confirmed:
    w = csv.writer(confirmed, delimiter = '\t')
    if output_number == 1:
        pass
    else:
        with open(file_name, newline = '') as file:
            file_reader = csv.reader(file, delimiter = '\t')
            for row in file_reader:
                current_line = row
                first_read = current_line[3]
                while current_line[3] == first_read:
                    w.writerow(current_line)
                    current_line = next(file_reader)
                break