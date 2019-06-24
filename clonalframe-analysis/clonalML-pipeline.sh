#!/bin/bash

WORKDIR=/Volumes/Bennett_BACKUP/Research/curtobacterium/curto-population-genomics/mauve_new
OUTDIR=$WORKDIR/clonalframe

mkdir -p $OUTDIR
cd $WORKDIR

mv newnames.txt $OUTDIR/

threshold=1500

# command strips out variable regions from the alignment to leave only core alignment blocks longer than ${threshold} nt bp
stripSubsetLCBs mauve_alignment.xmfa mauve_alignment.xmfa.bbcols core_alignment${threshold}.xmfa ${threshold}

xmfa2fasta.pl --file core_alignment${threshold}.xmfa > $OUTDIR/core_alignment${threshold}.fa

cd $OUTDIR

# mauve exports the full_alignment file as a iterative genome name
# go to the .log file and get the order
# also check that MMLR14_002 and MMLR14_014 are SUPER clsoely related for the .guide_tree file to double check
fasta-rename.py core_alignment${threshold}.fa newnames.txt core_alignment${threshold}.new.fa

# run clonalframe using the MAUVE core genome alignment and the
# PROKKA/ROARY coregenome tree generated using my own scripts and RAxML tree
# tree needs to be read in and exported as .nwk file for clonalframe
# make sure tip nodes and newnames.txt as the EXACT same

# # inference of recombination in bacterial genomes
# # raxml_coregenes90.nwk == output from ROARY ortholog calling
# # raxml_corealignment.nwk == output from MAUVE coregenome alignment (core_alignment${threshold}.new.fa)
# # both trees generated using:  $ raxml -s $ALIGNMENT -m GTRGAMMA -x 100 -# 100 -p 4321 -f a -T 64

# ClonalFrameML raxml_coregenes90.nwk core_alignment${threshold}.new.fa curtocore${threshold}
# ClonalFrameML raxml_corealignment.nwk core_alignment${threshold}.new.fa corealign${threshold}

# # one the output of clonalframe
# Rscript /Users/alexchase/software/ClonalFrameML/src/cfml_results.R curtocore${threshold}
# Rscript /Users/alexchase/software/ClonalFrameML/src/cfml_results.R corealign${threshold}



# now by each population
search-fasta.py -i core_alignment1500.new.fa -m population2.txt -o core_alignment1500.pop2.fa
raxml -s core_alignment1500.pop2.fa -m GTRGAMMA -n corealignment_pop2 -x 100 -# 100 -p 4321 -f a -T 2
ClonalFrameML raxml_corealignment_pop2.nwk core_alignment1500.pop2.fa corealign1500_pop2

search-fasta.py -i core_alignment1500.new.fa -m population4.txt -o core_alignment1500.pop4.fa
raxml -s core_alignment1500.pop4.fa -m GTRGAMMA -n corealignment_pop4 -x 100 -# 100 -p 4321 -f a -T 2
ClonalFrameML raxml_corealignment_pop4.nwk core_alignment1500.pop4.fa corealign1500_pop4

search-fasta.py -i core_alignment1500.new.fa -m population3.txt -o core_alignment1500.pop3.fa
raxml -s core_alignment1500.pop3.fa -m GTRGAMMA -n corealignment_pop3 -x 100 -# 100 -p 4321 -f a -T 2
ClonalFrameML raxml_corealignment_pop3.nwk core_alignment1500.pop3.fa corealign1500_pop3

