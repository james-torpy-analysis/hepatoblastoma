home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/hepatoblastoma"
raw_dir="$project_dir/raw_files"

samplenames=( "324_070" )

for s in ${samplenames[@]}; do

  replicates=( $(ls $raw_dir | grep $s | sed "s/\///") )
  mkdir -p $raw_dir/$s\_combined

  for r in ${replicates[@]}; do

    gunzip $raw_dir/$r/$r*.fastq.gz

    wc $raw_dir/$r/$r*R1.fastq >> $raw_dir/$s\_merge.log
    wc $raw_dir/$r/$r*R2.fastq >> $raw_dir/$s\_merge.log

    cat $raw_dir/$r/$r*R1.fastq >> $raw_dir/$s\_combined/$s\_combined_R1.fastq
    cat $raw_dir/$r/$r*R2.fastq >> $raw_dir/$s\_combined/$s\_combined_R2.fastq

    wc $raw_dir/$s\_combined/$s\_combined_R1.fastq >> $raw_dir/$s\_merge.log
    wc $raw_dir/$s\_combined/$s\_combined_R2.fastq >> $raw_dir/$s\_merge.log

  done

  gzip $raw_dir/$s\_combined/*

done