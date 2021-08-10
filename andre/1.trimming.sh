#!/bin/bash

imagename="https://depot.galaxyproject.org/singularity/cutadapt:3.4--py39h38f01e4_1"

home_dir="/share/ScratchGeneral/jamtor"
raw_dir="$home_dir/projects/hepatoblastoma/raw_files"
result_dir="/share/ScratchGeneral/jamtor/projects/hepatoblastoma/andre/trimmed"

for inF in $raw_dir/*/*R1*fastq.gz; do

  cFile=$(basename $inF)
  cFile2=$(echo $cFile | sed 's/_R1/_R2/g' )
  cID=$(echo $cFile | sed "s/_R1.fastq.gz//")
  
  ls $raw_dir/$cID/$cFile $raw_dir/$cID/$cFile2
  
  f1out=$result_dir/$cID/$(basename $cFile | sed "s/.fastq.gz//")
  f2out=$result_dir/$cID/$(basename $cFile2 | sed "s/.fastq.gz//")

  mkdir -p $result_dir/$cID
  
  singularity run \
    -B $result_dir/$cID:$result_dir/$cID \
    -B $raw_dir/$cID:$raw_dir/$cID \
    ${imagename} cutadapt \
    -j 4 \
    -a CTGTCTCTTATACACATCT \
    -A CAAAACGCAATACTGTACATTCTGTCTCTTATACACATCT \
    -o $f1out.tr.fastq.gz \
    -p $f2out.tr.fastq.gz \
    -m 20 \
    $raw_dir/$cID/$cFile $raw_dir/$cID/$cFile2 > $f1out.log

done
