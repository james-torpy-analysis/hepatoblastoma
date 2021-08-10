#!/bin/bash

imagename="https://depot.galaxyproject.org/singularity/blat:35--2"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/hepatoblastoma"
out_dir="$project_dir/andre/blat"
fasta_dir="$project_dir/andre/trimmed/fasta"
script_dir="$project_dir/scripts/andre"
filt_dir="$project_dir/andre/alignment"

mkdir -p $out_dir
mkdir -p $filt_dir
cd $out_dir

for inF in $fasta_dir/*.fa; do

  cSample=$(basename $inF | cut -f 2 -d _ | cut -f 1 -d .)
  echo $cSample

  singularity run \
    -B $fasta_dir:$fasta_dir \
    -B $out_dir:$out_dir \
    ${imagename} \
    blat $out_dir/CTNNB1_loc.fa \
    $inF \
    $out_dir/$cSample.loc.psl \
    -t=dna \
    -q=dna \
    -noHead

  perl  $script_dir/5.filter_blat.pl $project_dir/andre $cSample # | head

done

cat *.tab | grep -P "^[S0-9]" > calls.summary
cat *.tab | grep -vP "^[S0-9]" > evid.summary