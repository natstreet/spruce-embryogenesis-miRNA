#!/hpc2n/eb/software/Python/3.11.5-GCCcore-13.2.0/bin/python
# modified from codes of Dr Ahsan Z Rizvi
# to read fastq, by Fai Theerarat Kochakarn
#########################
import numpy
import sys
import gzip
#########################

seq1=list()
with gzip.open(sys.argv[1], 'rt') if sys.argv[1].endswith('.gz') else open(sys.argv[1]) as file:
    for i, line in enumerate(file):
        if i % 4 == 1:  # line 2 of every 4 lines
           line2 = line.strip()  # remove newline
           N = len(line2)
           if N > 16 and N < 25: # here I put N to > 16 and < 25 because len() is now after removing newline
              seq1.append(line2) 

uniq,counts = numpy.unique(seq1,return_counts=True)

#####################################

with open(sys.argv[1]+'.dedup.fa', 'w') as f:
    for i, sequence in enumerate(uniq):
        string='>read'+str("{:08d}".format(i))+'_x'+str(counts[i])+'\n'
        f.write(string)
        f.write(sequence+'\n')
