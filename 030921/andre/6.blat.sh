#!/bin/bash

imagename="https://depot.galaxyproject.org/singularity/blat:35--2"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/hepatoblastoma"
fasta_dir="$project_dir/andre/trimmed/fasta"
ref_dir="$project_dir/refs"
out_dir="$project_dir/andre/blat"
script_dir="$project_dir/scripts/andre"

mkdir -p $out_dir
cd $out_dir

for inF in $fasta_dir/*.fa; do

  cSample=$(basename $inF | cut -f 2 -d _ | cut -f 1 -d .)
  echo $cSample

#  singularity run \
#    -B $fasta_dir:$fasta_dir \
#    -B $out_dir:$out_dir \
#    ${imagename} \
#    blat $out_dir/CTNNB1_loc.fa \
#    $inF \
#    $out_dir/$cSample.loc.psl \
#    -t=dna \
#    -q=dna \
#    -noHead

  
  blat $ref_dir/CTNNB1_loc.fa \
    $inF \
    $out_dir/$cSample.loc.psl \
    -t=dna \
    -q=dna \
    -noHead

  perl  $script_dir/temp.pl $project_dir/andre $cSample # | head

done

cat *.tab | grep -P "^[S0-9]" > calls.summary
cat *.tab | grep -vP "^[S0-9]" > evid.summary