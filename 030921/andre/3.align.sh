#!/bin/bash

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/hepatoblastoma"
result_dir="$project_dir/andre/trimmed"
genome_dir="$project_dir/genome"
in_dir="$result_dir/no_barcode"
out_dir="$project_dir/andre/alignment"

mkdir -p $out_dir
cd $out_dir

for inF in $in_dir/*R1*.fastq.gz; do

  echo $inF

  cSample=$(basename $inF | cut -d "_" -f 2)_tr_bc
  cFile=$(basename $inF)
  cDir=$(dirname $inF)
  cFile2=$(echo $cFile | sed 's/_R1/_R2/g')

  echo $cSample $cFile $cFile2

  bwa mem -t 16 -Ma $genome_dir/GRCh37.p13.genome.fa $cDir/$cFile $cDir/$cFile2 | 
    samtools view -O BAM -o ${cSample}.bam

  samtools fixmate -m ${cSample}.bam ${cSample}.bam_
  samtools sort ${cSample}.bam_ | samtools markdup -m s - ${cSample}.bam
  rm ${cSample}.bam_
  samtools index ${cSample}.bam

done
