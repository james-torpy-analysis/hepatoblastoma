#!/bin/bash

### from qiagen DNA pipeline github, https://github.com/qiaseq/qiaseq-dna ###

# make sure docker engine is running!

# fastqs must be named as 'sample_name_R1_001.fastq.gz' and 
# 'sample_name_R2_001.fastq.gz'

# define and create directories:
sample_name=$1
capture_id=$2
genome=$3
#sample_name="324_032_D9HGF_TAGGCATG-CTCTCTAT_L001"
#capture_id="CDHS-33412Z-324"
#genome="hg19"

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
\cp $code_dir/qiaseq-dna/run_sm_counter_v2.params.txt \
  $working_dir/template_params.txt

if [ "$genome" == "hg38" ]; then
  cat $working_dir/template_params.txt | \
    sed "s/insert_sample_name/$sample_name/g" | \
    sed "s/insert_capture_id/$capture_id/g" | \
    sed "s/insert_genome_path/hg38\/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna/" | \
    sed "s/insert_snpEff_genome/GRCh38.92/" | \
    sed "s/insert_genome_version/hg38/" \
    > $working_dir/run_sm_counter_v2.params.txt
else
  cat $working_dir/template_params.txt | \
    sed "s/insert_sample_name/$sample_name/g" | \
    sed "s/insert_capture_id/$capture_id/g" | \
    sed "s/insert_genome_path/hg19\/ucsc.hg19.fa/" | \
    sed "s/insert_snpEff_genome/GRCh37.75/" | \
    sed "s/insert_genome_version/hg19/" \
    > $working_dir/run_sm_counter_v2.params.txt
fi

# copy primer and roi files to working_dir:
\cp $ref_dir/$capture_id.primers_$genome.txt $working_dir
\cp $ref_dir/$capture_id.roi_$genome.bed $working_dir

# run pipeline through container:
singularity exec \
  -B $code_dir:/srv/qgen/code/ \
  -B $data_dir:/srv/qgen/data/ \
  -B $working_dir:/srv/qgen/example/ \
  --pwd /srv/qgen/example/ \
  $cont_dir/qiaseq-dna_latest.sif python /srv/qgen/code/qiaseq-dna/run_qiaseq_dna.py \
  run_sm_counter_v2.params.txt \
  v2 \
  single \
  $sample_name > run.log 2>&1

# remove fastqs from working_dir::
rm $working_dir/*.fastq.gz