#!/bin/bash
#Arguments:
# -t (int) - specify number of threads (default: 8)
# -H (host_genome.fa/.fna/.fasta - can be gzipped)
# -i (folder containing reads)


Threads=8

while getopts 'i:H:t:' opt ; do
	case $opt in
    	i) Input=$OPTARG ;;
    	H) Host=$OPTARG ;;
    	t) Threads=$OPTARG ;;
    esac
done

cd "$Input"
for R1 in *R1*.fastq.gz
do
  echo "$R1"
  R2=$(echo "$R1" | sed 's/_R1_/_R2_/')
  echo "$R2"
  sfile1=$(basename "$R1")
  sfile2=$(basename "$R2")
  samplename=${sfile1%%_S*}

  bwa mem -t "$Threads" "$Host" "$R1" "$R2" | samtools view -@ "$Threads" -b -f 4 -o "$Input"/"$samplename"_nonHost.bam -;
  samtools sort -@ "$Threads" -n -O BAM "$Input"/"$samplename"_nonHost.bam | samtools fastq -@ "$Threads" -1 "$Input"/"$samplename"_NH1.fastq -2 "$Input"/"$samplename"_NH2.fastq -s "$Input"/"$samplename"_NHS.fastq -;
  
done
