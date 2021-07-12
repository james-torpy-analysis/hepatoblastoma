#!/bin/bash

source ~/.bashrc

project_name=$1
sample_name=$2
#project_name="hepatoblastoma"
#sample_name="324_001_DB674_TAAGGCGA-CTCTCTAT_L001"

printf "\n\n"
echo "project name is $project_name"
printf "\n"

printf "\n\n"
echo "sample name is $sample_name"
printf "\n"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project_name"
fq_dir="$project_dir/raw_files/$sample_name"
result_dir="$project_dir/results"
collapse_dir="$project_dir/results/collapsing"
int_dir="$result_dir/BWA_and_picard/int_bams/$sample_name"
fq_out_dir="$result_dir/BWA_and_picard/fastqs/$sample_name"
bam_dir="$result_dir/BWA_and_picard/bams/$sample_name"
genome_dir="$project_dir/genome"
stats_dir="$result_dir/read_stats/$sample_name"
script_dir="$project_dir/scripts"

picard_dir="$home_dir/local/bin"
fgbio_dir="$home_dir/local/lib/fgbio/target/scala-2.13"

mkdir -p $collapse_dir
mkdir -p $int_dir
mkdir -p $fq_out_dir
mkdir -p $bam_dir
mkdir -p $stats_dir

cd $collapse_dir

# call snk conda env for samtools/python:
conda activate snkenv

if [ ! -f $fq_out_dir/$sample_name.withUMI.fastq ]; then

  printf "\n\n"
  echo "--------------------------------------------------"
  echo "converting fastq files to one unmapped bam file..."
  echo "--------------------------------------------------"
  printf "\n"
  java -jar $picard_dir/picard.jar FastqToSam \
    FASTQ=$fq_dir/$sample_name\_R1.fastq.gz \
    FASTQ2=$fq_dir/$sample_name\_R2.fastq.gz \
    O=$int_dir/$sample_name.unmapped.bam \
    SM=sample
  
  printf "\n\n"
  echo "--------------------------------------------------"
  echo "extracting UMI and common sequence from R2..."
  echo "--------------------------------------------------"
  printf "\n"
  java -jar $fgbio_dir/fgbio-1.4.0-468a843-SNAPSHOT.jar ExtractUmisFromBam \
    --input=$int_dir/$sample_name.unmapped.bam \
    --output=$int_dir/$sample_name.unmapped.withUMI.bam \
    --read-structure=1M149T 23M152T \
    --molecular-index-tags=ZA ZB \
    --single-tag=RX
  
  printf "\n\n"
  echo "--------------------------------------------------"
  echo "converting back to single fastq file..."
  echo "--------------------------------------------------"
  printf "\n"
  java -jar $picard_dir/picard.jar SamToFastq I=$int_dir/$sample_name.unmapped.withUMI.bam \
    F=$fq_out_dir/$sample_name.withUMI.fastq \
    INTERLEAVE=true

else

  printf "\n"
  echo "$fq_out_dir/$sample_name.withUMI.fastq already exists, skipping to first alignment step..."
  printf "\n"
  
fi;


##################################################################################################################################
### 1. Uncollapsed bam ###
##################################################################################################################################

if [ ! -f $bam_dir/$sample_name.uncollapsed.bam ]; then

  printf "\n\n"
  echo "--------------------------------------------------"
  echo "aligning fastq..."
  echo "--------------------------------------------------"
  printf "\n"
  bwa mem -p -t 5 $genome_dir/GRCh37.p13.genome.fa \
    $fq_out_dir/$sample_name.withUMI.fastq > $int_dir/$sample_name.initial_mapped.bam
  
  # create reference dƒictionary if needed:
  #java -jar $picard_dir/picard.jar CreateSequenceDictionary -R $genome_dir/GRCh37.p13.genome.fa
  
  printf "\n\n"
  echo "--------------------------------------------------"
  echo "merging UMIs and common sequences to bam file..."
  echo "--------------------------------------------------"
  printf "\n"
  java -jar $picard_dir/picard.jar MergeBamAlignment \
    UNMAPPED=$int_dir/$sample_name.unmapped.withUMI.bam \
    ALIGNED=$int_dir/$sample_name.initial_mapped.bam \
    O=$int_dir/$sample_name.initial_mapped_and_UMI.bam \
    R=$genome_dir/GRCh37.p13.genome.fa \
    SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 \
    ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true
  
  printf "\n\n"
  echo "--------------------------------------------------"
  echo "removing fake read 1 UMI and common sequence..."
  echo "--------------------------------------------------"
  printf "\n"
  
  # remove read 1 UMI and common sequence using python script:
  python $script_dir/1a.remove_common_sequence.py \
    $project_name $int_dir/$sample_name.initial_mapped_and_UMI.bam 

else

  printf "\n"
  echo "$bam_dir/$sample_name.uncollapsed.bam already exists, skipping to UMI collapsing step..."
  printf "\n"
  
fi;

# isolate discordant reads:
samtools view -h -F 1294 $bam_dir/$sample_name.uncollapsed.bam | \
  samtools view -bh > $bam_dir/$sample_name.uncollapsed.discordant.bam
samtools index $bam_dir/$sample_name.uncollapsed.discordant.bam

# isolate split reads:
samtools view -h -f 2048 $bam_dir/$sample_name.uncollapsed.bam | \
  samtools view -bh > $bam_dir/$sample_name.uncollapsed.split.bam
samtools index $bam_dir/$sample_name.uncollapsed.split.bam


##################################################################################################################################
### 2. Collapse bam using Picard MarkDuplicates ###
##################################################################################################################################

if [ ! -f $int_dir/$sample_name.markdups.bam.bai ]; then

  printf "\n\n"
  echo "--------------------------------------------------"
  echo "collapsing UMIs with Picard MarkDuplicates..."
  echo "--------------------------------------------------"
  printf "\n"
  
  java -jar $picard_dir/picard.jar MarkDuplicates \
    I=$bam_dir/$sample_name.uncollapsed.bam \
    O=$int_dir/$sample_name.markdups.bam \
    M=$sample_name.markdups_metrics.txt \
    BARCODE_TAG="rx" \
    REMOVE_DUPLICATES=True
  
  samtools index $int_dir/$sample_name.markdups.bam

else

  printf "\n"
  echo "$int_dir/$sample_name.markdups.bam already exists, skipping to read stratification step..."
  printf "\n"
  
fi;


##################################################################################################################################
### 3. Separate primary, discordant and split reads ###
##################################################################################################################################

printf "\n\n"
echo "--------------------------------------------------"
echo "preparing files for split read duplicate removal..."
echo "--------------------------------------------------"
printf "\n"

# filter and index primary reads:
samtools view -h -F 2048 $int_dir/$sample_name.markdups.bam | \
  samtools view -bh > $int_dir/$sample_name.markdups.primary.bam
samtools index $int_dir/$sample_name.markdups.primary.bam

# filter and index discordant reads:
samtools view -h -F 1294 $int_dir/$sample_name.markdups.bam | \
  samtools view -bh > $int_dir/$sample_name.markdups.discordant.bam
samtools index $int_dir/$sample_name.markdups.discordant.bam

# filter split reads:
samtools view -hf 2048 $int_dir/$sample_name.markdups.bam | \
  samtools view -bh > $int_dir/$sample_name.markdups.split.supp.bam
  samtools index $int_dir/$sample_name.markdups.split.supp.bam


##################################################################################################################################
### 4. Remove split reads without matching primary alignments, and separate reads ###
##################################################################################################################################

printf "\n\n"
echo "---------------------------------------------------------------------------------------"
echo "removing split read duplicates and writing to $bam_dir/$sample_name.consensus.bam..."
echo "---------------------------------------------------------------------------------------"
printf "\n"

python $script_dir/1b.filter_split_reads.py $project_name $sample_name


printf "\n\n"
echo "--------------------------------------------------"
echo "indexing final bam files..."
echo "--------------------------------------------------"
printf "\n"

# index bam:
samtools index $bam_dir/$sample_name.consensus.bam

# filter and index concordant pairs:
samtools view -hf 0x2 $bam_dir/$sample_name.consensus.bam | \
  samtools view -bh > $bam_dir/$sample_name.consensus.concordant.pairs.bam
samtools index $bam_dir/$sample_name.consensus.concordant.pairs.bam

# filter and index discordant reads:
samtools view -h -F 1294 $bam_dir/$sample_name.consensus.bam | \
  samtools view -bh > $bam_dir/$sample_name.consensus.discordant.bam
samtools index $bam_dir/$sample_name.consensus.discordant.bam

# filter split reads:
samtools view -h -f 2048 $bam_dir/$sample_name.consensus.bam | \
  samtools view -bh > $bam_dir/$sample_name.consensus.split.bam
samtools index $bam_dir/$sample_name.consensus.split.bam




