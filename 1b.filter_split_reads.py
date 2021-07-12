#!/share/ClusterShare/thingamajigs/jamtor/local/lib/miniconda3/envs/snkenv/bin/python

import sys
import os
import re
import pysam

project_name = sys.argv[1]
sample_name = sys.argv[2]

home_dir = '/share/ScratchGeneral/jamtor/'
project_dir = home_dir + 'projects/' + project_name + '/'
result_dir = project_dir + 'results/'
in_dir = result_dir + 'BWA_and_picard/int_bams/' + sample_name + '/'
out_dir = result_dir + 'BWA_and_picard/bams/' + sample_name + '/'

# define input and output files:
primbam_name = in_dir + '/' + sample_name + '.markdups.primary.bam'
splitbam_name = in_dir + '/' + sample_name + '.markdups.split.supp.bam'
bam_name = in_dir + '/' + sample_name + '.markdups.bam'
out_name = out_dir + re.sub('markdups', 'consensus', os.path.basename(bam_name))

# define filenames using pysam:
primbam = pysam.Samfile(primbam_name, 'rb')
splitbam = pysam.Samfile(splitbam_name, 'rb')
bam = pysam.Samfile(bam_name, 'rb')
outbam = pysam.Samfile(out_name, 'wb', template = bam)
# fetch primary alignment qnames:
primary_qnames = []
for read in primbam.fetch(until_eof=True):
    primary_qnames.append(read.query_name)
print('Number of primary alignments = ', len(primary_qnames))
# fetch readnames of split reads without primary alignments:
throw_split = []
total_split = 0
for read in splitbam.fetch(until_eof=True):
    total_split += 1
    qname = read.query_name
    if qname not in primary_qnames:
        throw_split.append(qname)
print('Number of split alignments = ', total_split)
print('Number of split alignments without matching primaries = ', len(throw_split))
# define counters for reads processed and written:
n = 0
w = 0
t = 0
# filter out split reads without primary alignments:
for read in bam.fetch(until_eof=True):
    n += 1
    qname = read.query_name
    if qname not in throw_split:
        w += 1
        outbam.write(read)
    else:
        t += 1
print('Number alignments written to output = ', w)
print('Number alignments thrown = ', t)

primbam.close()
splitbam.close()
bam.close()
outbam.close()