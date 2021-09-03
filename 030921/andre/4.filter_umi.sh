#!/bin/bash

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/hepatoblastoma"
in_dir="$project_dir/andre/alignment"
out_dir="$in_dir/only_UMIs"

mkdir -p $out_dir
cd $out_dir

for inB in $in_dir/*.bam; do

  echo $inB

  perl -e '$outf=$ARGV[0]; $outf=~ s/[.]lumpy//g;
  
  $cSample=$1 if $ARGV[0] =~ /([S0-9]+)_tr_bc/;
  
  open(OUT, " | samtools view -b > $cSample.UMIs.bam ") || die "file not found";
  
  open(IN1, "samtools view -h $ARGV[0] | ") || die "file not found";
  
  while(<IN1>){ chomp; @_=split("\t",$_);
  
  if (/^@/){print OUT "$_\n" ;next}
  
  ($id, $chr, $st)=($_[0], $_[2], $_[3]);
  
  if ($id =~ /_([ATCG]+)$/){
  
   if (exists($umi_2_chr_2_pos{$1}{$chr})){
  
    if( abs($umi_2_chr_2_pos{$1}{$chr}-$st)<500){
  
     $umi_2_chr_2_pos{$1}{$chr}=$st; 
  
     if (!exists($whitel{$id})){ $c++; next; } # print OUT "#$_\n";
  
    }
  
   }
  
   $umi_2_chr_2_pos{$1}{$chr}=$st;
  
   $whitel{$id}++;
  
   print OUT "$_\n" ;
  
  }
  
  }close(IN1);
  
  close(OUT); print STDERR "c:$c\n";
  
  `samtools index $cSample.UMIs.bam`;
  
  ' $inB | head

done
