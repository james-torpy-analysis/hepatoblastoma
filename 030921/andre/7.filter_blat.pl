



# 1.	matches - Number of bases that match that aren't repeats
# 2.	misMatches - Number of bases that don't match
# 3.	repMatches - Number of bases that match but are part of repeats
# 4.	nCount - Number of "N" bases
# 5.	qNumInsert - Number of inserts in query
# 6.	qBaseInsert - Number of bases inserted in query
# 7.	tNumInsert - Number of inserts in target
# 8.	tBaseInsert - Number of bases inserted in target
# 9.	strand - "+" or "-" for query strand. For translated alignments, second "+"or "-" is for target genomic strand.
# 10.	qName - Query sequence name
# 11.	qSize - Query sequence size.
# 12.	qStart - Alignment start position in query
# 13.	qEnd - Alignment end position in query
# 14.	tName - Target sequence name
# 15.	tSize - Target sequence size
# 16.	tStart - Alignment start position in target
# 17.	tEnd - Alignment end position in target
# 18.	blockCount - Number of blocks in the alignment (a block contains no gaps)
# 19.	blockSizes - Comma-separated list of sizes of each block. If the query is a protein and the target the genome, blockSizes are in amino acids. See below for more information on protein query PSLs.
# 20.	qStarts - Comma-separated list of starting positions of each block in query
# 21.	tStarts - Comma-separated list of starting positions of each block in target

$stem=shift(@ARGV);
$cSample=shift(@ARGV);

$chr=3;
$minBlockSize_master=12;

#######
$_path="$stem/blat/*$cSample.loc.psl";
@t=glob($_path);
die "no such file $_path" if scalar(@t)==0;
$inPsl=$t[0];
$outBed=$inPsl;
$outBed=~ s/psl$/f1.bed/;

$outTab=$inPsl;
$outTab=~ s/psl$/f1.tab/;

##### read fastq reads
$_path="$stem/trimmed/no_barcode/*/3*_$cSample\_*fastq.gz";
$count=0;
foreach $f (glob($_path)){
print "reading $f ...\n";
$count++;
#Perl_e_2input_files

open(IN1, "zcat $f | ") || die "file not found";
while($a1=<IN1> and $a2=<IN1> and $a3=<IN1> and $a4=<IN1>){
$id2fq_read{$count}{$id}="$a1$a2$a3$a4";
$id=$1 if $a1 =~ /^@([^ ]+) /;chomp($a2);
#print "$id # $a2 \n";
$id2read{$count}{$id}=$a2;
}close(IN1);

}
die "!2 fastq files with $_path $count " if $count!=2;




$_path="$stem/alignment/only_UMIs/$cSample.UMIs.bam";
@t=glob($_path);
die "no such file $_path" if scalar(@t)==0;
$inBamUMI=$t[0];


$_path="$stem/alignment/*$cSample\_tr_bc.bam";
@t=glob($_path);
die "no such file $_path" if scalar(@t)==0;
$inBam=$t[0];


%lng;
%all;







print STDERR "$inPsl\n";
open(BED, ">$outBed") || die " $outBed out file not found";
open(TAB, ">$outTab") || die " $outTab out file not found";
open(IN1, "<$inPsl") || die " $inPsl file not found";

print BED "#gffTags\n";

while(<IN1>){ chomp; @_=split("\t",$_);


next if $_[17]==1; # must be split read
# next if $_[8]==3; # gap on ref must be at least 8 bases

chop($_[18]); @bl=split(",",$_[18]);
$is_long=scalar(grep {$_*1 >= $minBlockSize_master } @bl) < scalar(@bl) ? 0:1;

#print $_[18]." ".scalar(grep {$_*1 > $minBlockSize_master } @bl).">".scalar(@bl)."\n";

##### select master split reads with at least 13 bases on either side


$_[13]="3"; $_[14]=198022430; $_[15]+=41265000; $_[16]+=41265000;
chop($_[20]); @st=split(",",$_[20]); @st = map { $_=$_+41265000-$_[15] } @st;

$_[20]=join(",",@st);#.","; 

$UMI="_";
if ($_[9] =~ /_([ATCG]+)$/){
	$UMI=$1;
	if (exists($umi_2_chr_2_pos{$UMI}{$_[13]})){
		if( abs($umi_2_chr_2_pos{$UMI}{$_[13]}-$_[15])<500  ){
			$stat{dupl_UMI}++;
			$umi_2_chr_2_pos{$UMI}{$_[13]}=$_[15];
			next;
		}
	}
	$umi_2_chr_2_pos{$UMI}{$_[13]}=$_[15];			
}
#print join("\t",@_)."\n";
$col=!$is_long ? "#b3b3b3":"#c34343";
#$st=($_[15]+$st[0]+$bl[0]-1); $en=($_[15]+$st[1]-2); $len=$en-$st+1; next if $len<=7;
$st=($_[15]+$st[0]+$bl[0]); $en=($_[15]+$st[1]-1); $len=$en-$st+1; next if $len<=7;

$del="3:$st-$en";

print BED join("\t",($_[13],$_[15],$_[16],"Name2=".$_[9].";q_gap_count=$_[4];q_gap_bases=$_[5];t_gap_count=$_[6];t_gap_bases=$_[7];block_sizes=$_[18];SVLEN=$len;DEL=3:$st-$en;", 60, $_[8], $_[15], $_[16], $col, $_[17], $_[18], $_[20]))."\n";

$sr={};
$$sr{st}=$st;
$$sr{en}=$en;
$$sr{DEL}="3:$st-$en";
$$sr{UMI}=$UMI;
$$sr{read_id}=$_[9];
$$sr{q_gap_count}=$_[4];
$$sr{q_gap_bases}=$_[5];
$$sr{t_gap_count}=$_[6];
$$sr{t_gap_bases}=$_[7];
$$sr{block_sizes}=$_[18];
$$sr{SVLEN}=$len;
$$sr{is_long}=$is_long;


if ($is_long){
	$lng{$del}{st}=$st;
	$lng{$del}{en}=$en;
	$lng{$del}{DEL}="3:$st-$en";
	$lng{$del}{SVLEN}=$len;
	push @{$lng{$del}{SR}}, $sr;
	$lng{$del}{UMIs}++ if $UMI ne "_";
	$lng{$del}{no_UMIs}++ if $UMI eq "_";	

}

push @{$all{$del}{SR}}, $sr;
$all{$del}{UMIs}++ if $UMI ne "_";
$all{$del}{no_UMIs}++ if $UMI eq "_";

}close(IN1);
close(BED);

$lng_c=0;

foreach $del (sort { scalar(@{$all{$b}{SR}}) <=> scalar(@{$all{$a}{SR}})    } keys %lng){
	
	$lng_c++;
	$st=$lng{$del}{st};
	$en=$lng{$del}{en};
	$nb_SR_UMI_long=exists($lng{$del}{UMIs}) ? $lng{$del}{UMIs}:0;
	$nb_SR_no_UMI_long=exists($lng{$del}{no_UMIs}) ? $lng{$del}{no_UMIs}:0;
	
	$nb_SR_UMI=exists($all{$del}{UMIs}) ? $all{$del}{UMIs}:0;
	$nb_SR_no_UMI=exists($all{$del}{no_UMIs}) ? $all{$del}{no_UMIs}:0;
		
	######### get depth, only use UMI to determine depth (else VAF unreal)
	if ($nb_SR_UMI>0){
	$vaf_base_on="UMI_reads";
	#print "nb_SR_UMI: $nb_SR_UMI, $inBamUMI\n";
	$cov_up=get_depth_UMI($st-1,$st-1,$inBamUMI);
	$cov_dn=get_depth_UMI($en+1,$en+1,$inBamUMI);
	$VAF_up=$nb_SR_UMI/$cov_up;
	$VAF_dn=$nb_SR_UMI/$cov_dn;
	
	}else{
	#print "nb_SR_no_UMI: $nb_SR_no_UMI, $inBam\n";
	$vaf_base_on="all_reads";
	$cov_up=get_depth($st-1,$st-1,$inBam);
	$cov_dn=get_depth($en+1,$en+1,$inBam);
	$VAF_up=$nb_SR_no_UMI/$cov_up;
	$VAF_dn=$nb_SR_no_UMI/$cov_dn;
	
	}
	$VAF_avg=($VAF_up+$VAF_dn)/2;
	
	$seq_up=get_seq($st-1-100,$st-1);
	$seq_dn=get_seq($en+1,$en+1+100);
	######### compute depth
	
	######### frameshift
	$cds2_ovl=overlap($st,$en,41265560,41265572);
	$cds3_ovl=overlap($st,$en,41266017,41266244);
	$cds4_ovl=overlap($st,$en,41266445,41266698);
	
	$sum_cds=$cds2_ovl+$cds3_ovl+$cds4_ovl;
	$fr_shift=($sum_cds)%3!=0 ? "yes":"no";
	$fr_shift="-" if ($sum_cds)==0;

	######### annotate master
	$annot="";
	
	## confidence
	if($nb_SR_UMI>=5){
		$confidence = "high";
	}elsif($nb_SR_UMI>1){
		$confidence = "medium";
	}else{
		$confidence = "low";	
	}
	$annot.= "$confidence confident, ";

	## secondary
	$annot.=($lng_c==1) ? 'primary, ':'secondary, ';
	
	$br_seq=substr($seq_up,-5).substr($seq_dn,0,5);
	
	%pos2base={};
	$coord_confirmed_in_read=0;
	foreach $sr (@{$all{$del}{SR}}){
		$r1=$id2read{1}{$$sr{read_id}};
		$r2=$id2read{2}{$$sr{read_id}};
		
		$seq_p1='';
		$seq_p2='';
		if (index($r1,$br_seq)>-1){
			$stcoord=index($r1,$br_seq);			
			$seq_p1=sprintf "%+150s", substr($r1,0,$stcoord+5);
			$seq_p2=substr($r1,$stcoord+5);
			#pos_stat($seq_p1,$seq_p2);
			$br_seq_res='r1';
			$coord_confirmed_in_read++;
		}elsif(index($r2,$br_seq)>-1){
			$stcoord=index($r2,$br_seq);
			$seq_p1=sprintf "%+150s", substr($r2,0,$stcoord+5);
			$seq_p2=substr($r2,$stcoord+5);
			#pos_stat($seq_p1,$seq_p2);
			$br_seq_res='r2';
			$coord_confirmed_in_read++;
		}else{
			$br_seq_res='-';
		}
		print TAB join("\t",("",$cSample, $lng_c, $$sr{is_long} ? "long":"short", $$sr{UMI}, $$sr{DEL}, $$sr{read_id}, $$sr{q_gap_count}, $$sr{q_gap_bases}, $$sr{t_gap_count}, $$sr{t_gap_bases}, $$sr{block_sizes}, $$sr{SVLEN}, $r1, $r2, $br_seq, $br_seq_res, $seq_p1, $seq_p2))."\n";
	}


	print TAB join("\t", $cSample, $lng_c, $lng{$del}{DEL}, $seq_up, $seq_dn, $coord_confirmed_in_read, $lng{$del}{SVLEN}, $nb_SR_UMI, $nb_SR_no_UMI, $nb_SR_UMI_long, $nb_SR_no_UMI_long, 
								$cds2_ovl, $cds3_ovl, $cds4_ovl, $fr_shift, $cov_up, $cov_dn, $VAF_up, $VAF_dn, $VAF_avg, $vaf_base_on, $annot
								 )."\n";
		
# 	if ($confidence ne "low"){
# 		
# 		$rsteam="$stem/blat/assembly/reads/$cSample\_$lng_c";
# 		open(FQ1, ">$rsteam.read1.fastq") || die " $rsteam out file not found";
# 		open(FQ2, ">$rsteam.read2.fastq") || die " $rsteam out file not found";
# 
# 		%read_to_extract={};
# 		foreach $sr (@{$all{$del}{SR}}){
# 			$read_to_extract{$$sr{read_id}}++;
# 			print FQ1 $id2fq_read{1}{$$sr{read_id}};
# 			print FQ2 $id2fq_read{2}{$$sr{read_id}};
# 		}
# 		extract_reads(\%read_to_extract,$lng_c,$inBam);
# 		
# 		close(FQ1);
# 		close(FQ2);
# 		
# 	}

	
}
close(TAB);


sub overlap{ my($R1,$R2,$Q1,$Q2)=@_; ($QuerB,$QuerE) = sort {$a <=> $b} ($Q1,$Q2); ($RefB,$RefE)   = sort {$a <=> $b} ($R1,$R2);   $ovlLen=(sort {$a <=> $b} ($QuerE,$RefE))[0]-(sort {$b <=> $a} ($QuerB,$RefB))[0]+1; $returnVal=($ovlLen>0)? $ovlLen:0;return $returnVal;}

# sub pos_stat{
# 	my($r1,$r2)=@_;
# 
# 	@array = split(//, reverse($r1) );
# 	my $i=0;
# 	foreach my $base (@array){
# 		$i++;
# 		last if $base eq "";
# 		$pos2base{1}{$i}{$base}++;
# 	}
# 
# 
# 	@array = split(//, $r2);
# 	my $i=0;
# 	foreach my $base (@array){
# 		$i++;
# 		$pos2base{2}{$i}{$base}++;
# 	}
# }
# 
# sub pos_stat_out{
# 	
# 	for my $i (0..150){
# 		last if !exists()
# 		pos2base{1}{$i}
# 	}
# }


sub extract_reads{
my($read_to_extract,$lng_c,$inBam)=@_;
$cmd="samtools view -h -r 3:41238915-41283944 $inBam | ";

open(BAMOUT, " | samtools view -b - > $stem/blat/assembly/$cSample\_$lng_c.bam ") || die " $outBed out file not found";

open(BAM, $cmd) || die "file not found";
 while(<BAM>){ 
 if (/^@/){print BAMOUT $_; next;}
 chomp; @_=split("\t",$_); 
 next if !exists($$read_to_extract{$_[0]});
 print BAMOUT "$_\n"; 
}close(BAM);
close(BAMOUT);

}

sub get_depth{
my($st,$en,$inBam)=@_;
$cmd="samtools depth -a -r $chr:$st-$en $inBam | ";
#print STDERR " $cmd \n";
open(BAM, $cmd) || die "file not found";
 while(<BAM>){ chomp; @_=split("\t",$_); 
 #print "$$v[0] $$v[$i] $_[1] $_[2]\n";
 return $_[2];
 
}close(BAM);

}

sub get_seq{
my($st,$en)=@_;
$cmd="samtools faidx /g/data/gd7/resources/gatk-resource-bundle/2.8/b37/human_g1k_v37_decoy.fasta $chr:$st-$en | ";
#print STDERR " $cmd \n";
$seq='';
open(FASTA, $cmd) || die "file not found";
 while(<FASTA>){ chomp; next if />/; $seq.=$_;
}close(FASTA);
return $seq;
}


sub get_depth_UMI{
my($st,$en,$inBam)=@_;
$cmd="samtools mpileup -Q 0 -r $chr:$st-$en $inBam | cut -f 1,2,4 | ";
#print STDERR " $cmd \n";
open(BAM, $cmd) || die "file not found";
 while(<BAM>){ chomp; @_=split("\t",$_); 
 #print "$$v[0] $$v[$i] $_[1] $_[2]\n";
 return $_[2];
 
}close(BAM);

}






































