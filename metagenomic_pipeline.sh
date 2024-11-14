#!/bin/bash
#Arguments:
# -t (int) - specify number of threads (default: 8)
# -i (dir) - directory containing raw reads
# -d (dir) - kraken2 database directory
# -H (file) - host genome file

#dependencies: BWA, Samtools, BBTools (lilypad), SPAdes, viralverify (SPAdes add-on), pigz, Kraken2, KrakenTools, QUAST, HMMER (hmmsearch), BLAST
#required: Pfam-A (https://www.ebi.ac.uk/interpro/download/pfam/), some sort of blastdb either virus-nr or custom

Threads=16


while getopts 'i:t:d:H:' opt ; do
	case $opt in
    	i) Input=$OPTARG ;;
    	t) Threads=$OPTARG ;;
    	d) DB=$OPTARG ;;
    	H) Host=$OPTARG ;;
    esac
done

cd "$Input"
for R1 in *_R1*.fastq.gz
do
  echo "$R1"
  R2=$(echo "$R1" | sed 's/_R1_/_R2_/')
  echo "$R2"
  sfile1=$(basename "$R1")
  sfile2=$(basename "$R2")
  samplename=${sfile1%%_R1*}
  
  mkdir "$samplename"_results
  mkdir "$samplename"_results/kraken
  mkdir "$samplename"_results/SPAdes
  mkdir "$samplename"_results/verify
  mkdir "$samplename"_results/QUAST
  mkdir "$samplename"_results/lilypad
  mkdir "$samplename"_results/lilypad/QUAST
  mkdir "$samplename"_results/lilypad/verify
  
  viral1="$samplename"_results/"$samplename"_viral_1.fastq
  viral2="$samplename"_results/"$samplename"_viral_2.fastq
  NH1="$samplename"_results/"$samplename"_viral_NH1.fastq
  NH2="$samplename"_results/"$samplename"_viral_NH2.fastq
  
  #should be fine assuming kraken was installed via conda/apt
  kraken2 --db "$DB" --threads "$Threads" --confidence 0.1 --output "$samplename"_results/kraken/"$samplename".kraken --report "$samplename"_results/kraken/"$samplename".report --paired --gzip-compressed "$R1" "$R2"
  
  #should be fine assuming kronatools was installed via conda/apt
  ktImportTaxonomy -t 5 -m 3 -o "$samplename"_results/kraken/"$samplename".html "$samplename"_results/kraken/"$samplename".report
  
  #will require you to either add KrakenTools to $PATH, or edit and change to your install location
  extract_kraken_reads.py -k "$samplename"_results/kraken/"$samplename".kraken -r "$samplename"_results/kraken/"$samplename".report -t 2759 2157 2 2731619 2788787 10878 9606 --exclude --include-children --fastq-output -s "$R1" -s2 "$R2" -o "$viral1" -o2 "$viral2"
  
  rm "$samplename"_results/kraken/"$samplename".kraken
  
  #host removal
  bwa mem -t "$Threads" "$Host" "$viral1" "$viral2" | samtools view -@ "$Threads" -b -f 4 -o "$samplename"_results/"$samplename"_viral_nonHost.bam -;
  samtools sort -@ "$Threads" -n -O BAM "$samplename"_results/"$samplename"_viral_nonHost.bam | samtools fastq -@ "$Threads" -1 "$NH1" -2 "$NH2" -s "$samplename"_results/"$samplename"_viral_NHS.fastq
  
  #same as above, will work if SPAdes is in your $PATH - otherwise change
  spades.py -t "$Threads" --meta -k 21,33,55,77,99 -1 "$NH1" -2 "$NH2" -o "$samplename"_results/SPAdes
  
  
  SPAdes="$samplename"_results/SPAdes
  
  
  #zip up potentially huge reads
  pigz "$viral1"
  pigz "$viral2"
  
  #SPAdes Quality
  quast.py -o "$samplename"_results/QUAST "$SPAdes"/contigs.fasta
  
  cp "$SPAdes"/contigs.fasta "$samplename"_results/lilypad/contigs.fasta
  
  #lilypad prep
  bwa index "$samplename"_results/lilypad/contigs.fasta
  
  bwa mem -t 16 "$samplename"_results/lilypad/contigs.fasta "$R1" "$R2" | samtools fixmate -@ 16 -O bam - "$samplename"_results/lilypad/scaffolds.bam
  
  #lilypad - add own or path
  lilypad.sh in="$samplename"_results/lilypad/scaffolds.bam ref="$samplename"_results/lilypad/contigs.fasta out="$samplename"_results/lilypad/scaffolds.fasta
  
  #lilypad quality
  quast.py -o "$samplename"_results/lilypad/QUAST "$samplename"_results/lilypad/scaffolds.fasta

  #if you have Pfam-A and blastdb stored elsewhere please change this path - it is easier to keep them in the same folder. if you use a different blastdb please chane the "--db" flag below
  #viralverify on SPAdes + lilypad output
  cd /home/$USER/mnt/GAPDCII/blastdbs
  
  #make sure this also matches your install location
  /home/$USER/viralVerify/bin/viralverify -f "$Input"/"$SPAdes"/scaffolds.fasta -t "$Threads" --hmm Pfam-A.hmm.gz --db virus_all -o "$Input"/"$samplename"_results/verify
  
  /home/$USER/viralVerify/bin/viralverify -f "$Input"/"$samplename"_results/lilypad/scaffolds.fasta -t "$Threads" --hmm Pfam-A.hmm.gz --db virus_all -o "$Input"/"$samplename"_results/lilypad/verify
  
  cd "$Input"
  
  #cleanup huge BAM files
  rm "$samplename"_results/lilypad/scaffolds.bam
  rm "$samplename"_results/"$samplename"_viral_nonHost.bam
  
done

echo "Complete"
