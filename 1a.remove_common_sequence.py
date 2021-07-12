import pysam
import sys
import re
import os

project_name = 'hepatoblastoma'

home_dir = '/share/ScratchGeneral/jamtor/'
project_dir = home_dir + 'projects/' + sys.argv[1] + '/'
fq_dir = project_dir + 'raw_files/'
result_dir = project_dir + 'results/'
bam_dir = result_dir + 'BWA_and_picard/bams/'

# check whether bam file is provided:
if len(sys.argv) == 3:
	assert sys.argv[2].endswith('.bam')

# make output filehandle:
inbam_name = sys.argv[2]
#inbam_name = bam_dir + '409_004_D9YWF_GTAGAGGA-CTCTCTAT_L001.mapped_and_UMI.bam'
outbam_name = bam_dir + re.sub('.initial_mapped_and_UMI.bam$', '', os.path.basename(inbam_name)) + \
    '/' + re.sub('initial_mapped_and_UMI.bam$', 'uncollapsed.bam', os.path.basename(inbam_name))

print('Input bam = ', inbam_name)
print('Output bam = ', outbam_name)

# define filenames using pysam:
inbam = pysam.Samfile(inbam_name, 'rb')
outbam = pysam.Samfile(outbam_name, 'wb', template=inbam)

# define counters for reads processed and written:
n = 0
w = 0

print('Processing reads...')

for read in inbam.fetch(until_eof=True):
	n += 1
	# fetch ZB tag:
	umi1 = read.get_tag('ZB')
	assert umi1 is not None
	# remove last 11 characters:
	umi1_fix = re.sub('...........$', '', umi1)
	# replace original umi:
	read.set_tag('ZB', umi1_fix, value_type='Z')
	# fetch RX tag:
	umi2 = read.get_tag('RX')
	assert umi2 is not None
	# remove last 11 characters:
	umi2_fix = re.sub('[A-Z]-', '', 
		re.sub('...........$', '', umi1)
	)
	# replace original umi:
	read.set_tag('RX', umi2_fix, value_type='Z')

	# if length = 12, write to file:
	if len(umi1_fix) == 12 & len(umi2_fix) == 12:
		if '-' not in umi2_fix:
			w += 1
			outbam.write(read)
		else:
			print('Hyphen still in RX UMI sequence: ', umi2_fix, ' read=', read)
			break
	else:
		print(
			'UMI sequence wrong length: ZB=', umi1_fix, ' RX=', umi2_fix, ' read=', read
		)

print('All reads processed, checked and written')

inbam.close()
outbam.close()