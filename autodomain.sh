#!/bin/bash
# Dependencies:
# HMMsearch (HMMER)
# Easel - usually packaged with HMMER

# Required inputs:
# Fasta file containing one or more protein sequences.
# HMM file (can be .gz) of your desired target domains. These can be downloaded from InterPro. Files containing multiple domains can be used, but this is not the intended function.

VERSION=1.0

usage()
{ echo "
autodomain.sh Version $VERSION
Usage: bash $0 [options] -H <hmm database> -i <fasta file> -o <output directory>
	-h Displays this message
	-m Minimum domain length (aa), use to filter partial outputs [optional, but recommended]
	-t Number of threads/cores to use [optional] default = all
";
}

if [[ $1 == "-h" ]]; then
usage
exit 0
fi

Threads=$(grep -c ^processor /proc/cpuinfo)

while getopts 'i:H:t:o:m:h:' opt ; do
	case $opt in
		i) Input=$OPTARG ;;
		H) hmm=$OPTARG ;;
		t) Threads=$OPTARG ;;
		o) Outdir=$OPTARG ;;
		m) minlength=$OPTARG ;;
		\?) usage
		exit
		;;
	esac
done

echo -e "\nChecking options..."

#error handling stuff
if [[ $Input == "" ]] || [[ $hmm == "" ]] || [[ $Outdir == "" ]]; then
echo "ERROR: One or more arguments missing!"
usage
exit 1
fi

if [[ -d "$Outdir" ]]; then
echo "Output directory: $Outdir"
else
echo "ERROR: $Outdir is not a directory!"
usage
exit
fi

if [[ -f "$Input" ]] && [[ "$Input" =~ ".fa" ]]; then
echo "Input: $Input"
else
echo "ERROR: $Input is not a fasta file!"
usage
exit
fi

if [[ -f "$hmm" ]] && [[ "$hmm" =~ ".hmm" ]]; then
echo "HMM: $hmm"
else
echo "ERROR: $hmm is not a HMM file!"
usage
exit
fi

if [[ $minlength == "" ]]; then
echo "Filtering not active."
else
echo "Filtering at $minlength minimum length."
fi

echo -e "Proceeding... \n"

#actual code
filename=$(basename "$Input")
name=${filename%%.fa*}

esl-sfetch --index "$Input"

cd "$Outdir"

hmmsearch --noali --cpu "$Threads" -o "$name"_hmmout --domtblout "$name"_domtbl "$hmm" "$Input"
grep -v "^#" "$name"_domtbl | awk '{print $1"_"$4, $20, $21, $1}' | esl-sfetch -Cf "$Input" - > "$name"_domainhits.fasta
sed -i 's/ /_/g' "$name"_domainhits.fasta

#optional part for filtering
if [[ "$minlength" != "" ]]; then
echo "Filtering output for sequences >""$minlength"
seqtk seq -L "$minlength" "$name"_domainhits.fasta > "$name"_domainhits-filtered.fasta
echo "Filtered sequence hits can be found at""$name""_domainhits-filtered.fasta"
else
echo "Sequence hits can be found at ""$name""_domainhits.fasta"
fi
