#!/bin/bash

home_dir="/share/ScratchGeneral/jamtor"
result_dir="/share/ScratchGeneral/jamtor/projects/hepatoblastoma/andre/trimmed"
out_dir="$result_dir/fasta"

mkdir -p $out_dir

for inF in $result_dir/no_barcode/*R1*.fastq.gz; do

  cSample=$(basename $inF | cut -d "_" -f 2)
  echo $cSample

  inF2=$(echo $inF | perl -pe 's/_R1/_R2/' )

  zcat $inF $inF2 | \
    perl -e 'while(<>){ if ($.%4==1){ s/^@/>/; print  } if ($.%4==2){ print }   }' \
    > $out_dir/$cSample.fa

done
