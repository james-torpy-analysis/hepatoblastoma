#!/bin/bash

#sample_name="409_021_DBV4V_CTCTCTAC-CTCTCTAT_L001"
#filename="/share/ScratchGeneral/jamtor/projects/ewing_ctDNA/results/svaba/bwa/markduplicates/$sample_name/$sample_name.svaba.unfiltered.sv.vcf"

sample_name=$1
filename=$2

in_dir="${filename%/*}"
grep \## $filename > $in_dir/$sample_name.svaba.unfiltered.sv.formatted.vcf
grep -v \## $filename > $in_dir/$sample_name.svaba.unfiltered.sv.temp.vcf

awk '{print $1, $2, $3, $4, $5, $6, $7, $8}' OFS='\t' $in_dir/$sample_name.svaba.unfiltered.sv.temp.vcf >> $in_dir/$sample_name.svaba.unfiltered.sv.formatted.vcf
