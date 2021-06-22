#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage $0 <file with seq pairs to align>"
	echo "  will align all pairs from input file with water"
	exit 1
fi

seqspath=disordered_sequences_splitted_fasta
regspath=disordered_regions_splitted_fasta
outfn=$(echo $1|sed -e 's/lists/results_water/' -e 's/\.dat/\.txt/')
rm -f $outfn

echo -e "#disprotid1\disprotid2\id1\tid2\tlenght1\tlenght2\tidentity\tsimilarity\tgaps\tscore\tsequence1\tsequence2\treg1\treg2" > $outfn
#echo -e "#disprotid\tuniprotid\tregid\tlenght1\tlenght2\tidentity\tsimilarity\tgaps\tscore\tsequence1\tsequence2\tstart\tend\treg1\treg2" > $outfn # regseq

[ $(hostname) != 'quaoar' ] && module load emboss/6.6.0 # needs to be loaded on remote only
#module load emboss 2>/dev/null  # to avoid printing errors in local

while read l; do
	seq1=$(echo "$l"|cut -f1)
	seq2=$(echo "$l"|cut -f2)

	ab=$(water -outfile=/dev/stdout -asequence=$seqspath/$seq1.fasta -bsequence=$seqspath/$seq2.fasta  -gapopen=10 -gapextend=0.5 | parser/parse.py)
	ba=$(water -outfile=/dev/stdout -asequence=$seqspath/$seq2.fasta -bsequence=$seqspath/$seq1.fasta  -gapopen=10 -gapextend=0.5 | parser/parse.py)

#	ab=$(water -outfile=/dev/stdout -asequence=$seqspath/$seq1.fasta -bsequence=$regspath/$seq2.fasta  -gapopen=10 -gapextend=0.5 | parser/parse.py) # regseq
#	echo "$ab" >>$outfn
	abid=$(echo "$ab"|cut -f 5|bc -l)
        baid=$(echo "$ba"|cut -f 5|bc -l)
	if (( $(echo "$abid >= $baid" |bc -l) )); then
		echo "$ab" >>$outfn

	else
		echo "$ba" >>$outfn
	fi

done <<< "$(cat "$1")"
