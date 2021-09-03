#!/bin/bash

### from qiagen DNA pipeline github, https://github.com/qiaseq/qiaseq-dna ###

# make sure docker engine is running!

# fastqs must be named as 'sample_name_R1_001.fastq.gz' and 
# 'sample_name_R2_001.fastq.gz'

# define and create directories:
#sample_name=$1
#capture_id=$2
sample_name="324_044_D9YW9_GCTCATGA-CTCTCTAT_L001"
capture_id="CDHS-33412Z-324"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/hepatoblastoma"
cont_dir="$project_dir/containers"
ref_dir="$project_dir/refs"
raw_dir="$project_dir/raw_files/$sample_name"
working_dir="$project_dir/results/smcounter2/$sample_name"
data_dir="$project_dir/data"
code_dir="$project_dir/code"

mkdir -p $working_dir
mkdir -p $cont_dir
mkdir -p $code_dir

# copy fastqs to working dir:
\cp $raw_dir/*.fastq.gz $working_dir

# update params file with sample_name and copy to working_dir:
\cp $code_dir/qiaseq-dna/misc_workflow/run_consensus.params.txt \
  $working_dir/template_params.txt
cat $working_dir/template_params.txt | \
  sed "s/insert_sample_name/$sample_name/g" | \
  sed "s/insert_capture_id/$capture_id/g" \
  > $working_dir/run_consensus.params.txt

# copy primer and roi files to working_dir:
\cp $ref_dir/$capture_id.primers.txt $working_dir
\cp $ref_dir/$capture_id.roi.bed $working_dir

# run pipeline through container:
singularity exec \
  -B $code_dir:/srv/qgen/code/ \
  -B $data_dir:/srv/qgen/data/ \
  -B $working_dir:/srv/qgen/example/ \
  --pwd /srv/qgen/example/ \
  $cont_dir/qiaseq-dna_latest.sif python /srv/qgen/code/qiaseq-dna/misc_workflow/run_consensus.py \
  $sample_name \
  run_consensus.params.txt \
   > run.log 2>&1

# remove fastqs from working_dir::
rm $working_dir/*.fastq.gz