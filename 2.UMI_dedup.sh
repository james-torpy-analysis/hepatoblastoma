#!/bin/bash

source ~/.bashrc

sample_name=$1
#sample_name="324_057_DB674_CGTACTAG-CTCTCTAT_L001"

printf "\n\n"
echo "project name is hepatoblastoma"
printf "\n"

printf "\n\n"
echo "sample name is $sample_name"
printf "\n"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/hepatoblastoma"
result_dir="$project_dir/results"
in_dir="$result_dir/smcounter2/$sample_name"
bam_dir="$result_dir/picard/$sample_name"
stats_dir="$result_dir/read_stats/$sample_name"
script_dir="$project_dir/scripts"

picard_dir="$home_dir/local/bin"
fgbio_dir="$home_dir/local/lib/fgbio/target/scala-2.13"

mkdir -p $bam_dir
mkdir -p $stats_dir

cd $in_dir

# call snk conda env for samtools/python:
conda activate snkenv


##################################################################################################################################
### 1. Sort bam by queryname ###
##################################################################################################################################

if [ ! -f $bam_dir/$sample_name.sorted.by.query.bam ]; then

  samtools sort -n -o $bam_dir/$sample_name.align.sorted.by.query.bam \
    $in_dir/$sample_name.align.bam

else

  printf "\n"
  echo "bam_dir/$sample_name.sorted.by.query.bam already exists, skipping to dedup step..."
  printf "\n"
  
fi;


##################################################################################################################################
### 2. Collapse bam using Picard MarkDuplicates ###
##################################################################################################################################

if ! ( [ -f $bam_dir/$sample_name.dedup.sorted.by.coord.bam.bai ] && [ -f $bam_dir/$sample_name.dedup.sorted.by.coord.bam ] ); then

  printf "\n\n"
  echo "--------------------------------------------------"
  echo "collapsing UMIs with Picard MarkDuplicates..."
  echo "--------------------------------------------------"
  printf "\n"
  
  java -jar $picard_dir/picard.jar MarkDuplicates \
    I=$bam_dir/$sample_name.align.sorted.by.query.bam \
    O=$bam_dir/$sample_name.dedup.bam \
    M=$sample_name.dedup_metrics.txt \
    ASSUME_SORT_ORDER=queryname \
    BARCODE_TAG="mi" \
    REMOVE_DUPLICATES=True
  
  samtools sort -o $bam_dir/$sample_name.dedup.sorted.by.coord.bam $bam_dir/$sample_name.dedup.bam
  samtools index $bam_dir/$sample_name.dedup.sorted.by.coord.bam

  rm $bam_dir/$sample_name.align.sorted.by.query.bam
  rm $bam_dir/$sample_name.dedup.bam

  printf "\n"
  echo "output file is $bam_dir/$sample_name.dedup.bam..."
  printf "\n"

else

  printf "\n"
  echo "$bam_dir/$sample_name.dedup.sorted.by.coord.bam already exists, exiting..."
  printf "\n"
  
fi;


