#!/bin/bash

WORKDIR=/Volumes/Bennett_BACKUP/Research/curtobacterium/curto-population-genomics/roary_coregenome/coregenes90

# in order to calculate dN/dS ratios for each gene, we need to import into R
# BUT each alignment needs to be a multiple of 3 nt

rm -f $WORKDIR/checkgenes.txt
rm -f $WORKDIR/removed_genes.txt

mkdir -p $WORKDIR/goodalignments

for f in *.ffn.aln 
do

	gene=${f%.total.ffn.aln}
	numgen=1

	count-bp-in-seqs.sh $f > /dev/null 2>&1

	bpcount=$(grep -v ">" ${gene}.total.ffn-bp-count.txt | head -n28 | sort | uniq)

	#check to make sure alignments look good
	check=$(grep -v ">" ${gene}.total.ffn-bp-count.txt | head -n28 | sort | uniq | wc -l)

	if [ "$check" -eq "1" ]
		then
			:
		else 
			echo $gene >> $WORKDIR/removed_genes.txt
			# rm -f $WORKDIR/$gene.*
			echo "${gene} removed!"

	fi

	# check to make sure open reading frames are good to go!
	if (( $bpcount % 3 == 0 ))           # no need for brackets
		then
			mv $f $WORKDIR/goodalignments
			echo "$gene is good!"
		else
			echo "${gene} Invalid"
			echo $gene >> $WORKDIR/checkgenes.txt
	fi

	rm -f ${gene}.total.ffn-bp-count.txt

done
