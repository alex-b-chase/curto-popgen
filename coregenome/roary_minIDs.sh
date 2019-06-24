#!/bin/bash

BASE=/Volumes/Bennett_BACKUP/Research/curtobacterium/curto-population-genomics

cd $BASE/genomes

# # for ROARY, need .gff3 files, so:
# # convert each file to ffn and faa files for the coregenome analysis
# for f in *.gb 
# do
# 	bp_genbank2gff3.pl $f 
# 	gb2fasta.py $f ffn
# 	gb2fasta.py $f faa

# done

cp *.gff $BASE/roary_coregenome

cd $BASE/roary_coregenome

for minID in {70,75,80,85,90,95,97,99}
do

	cd $BASE

	rm -rf coregenome${minID}

	roary -f coregenome${minID} -p 4 -i ${minID} *.gff

	cd coregenome${minID}

	create_pan_genome_plots.R
	python $BASE/roary_plots.py accessory_binary_genes.fa.newick gene_presence_absence.csv 

done

rm -f *.gff
