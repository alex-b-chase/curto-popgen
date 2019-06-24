#!/bin/bash

# look for ICEs and prophage DNA in genomes to see if they are near the regions of interest

GENOMEDIR=/Volumes/Bennett_BACKUP/Research/curtobacterium/curto-population-genomics/genomes
HGTDIR=/Users/alexchase/software/phiSpyNov11_v2.3
ICEDIR=$HGTDIR/ICEberg
PHASTDIR=$HGTDIR/phaster 

# cd $ICEDIR

# makeblastdb -in ICEberg_seq.fas  -dbtype 'nucl' -out ICEberg

# cd $GENOMEDIR

# mkdir -p $GENOMEDIR/hgtelements

# for f in *.fasta
# do

# 	genome=${f%.bgpipe.fasta}

# 	blastn -db $ICEDIR/ICEberg -query $f -outfmt "6 qseqid sseqid pident evalue qstart qend qlen slen" \
# 	-max_target_seqs 2 -evalue .00001 -num_threads 4 > ${genome}.temp.blastICE.txt

# 	awk -v v1="$genome" '{print v1, $0}' FS="\t" OFS="\t" ${genome}.temp.blastICE.txt > ${genome}.blastICE.txt
# 	rm ${genome}.temp.blastICE.txt

# done

# echo -e "genome\tcontig\tICEhit\tpid\tevalue\tstartpos\tendpos\tcontiglen\tICElength" > $GENOMEDIR/hgtelements/total.ICE.txt
# cat *.blastICE.txt | sed 's/_final//g' >> $GENOMEDIR/hgtelements/total.ICE.txt
# rm *.blastICE.txt


# for g in *.gb
# do

# 	genome1=${g%.bgpipe.gb}
# 	genome=$(echo $genome1 | sed 's/_final//g')

# 	python $HGTDIR/genbank_to_seed.py $g $GENOMEDIR/hgtelements/${genome}/
# 	python $HGTDIR/phiSpy.py -i $GENOMEDIR/hgtelements/${genome}/ -o $GENOMEDIR/hgtelements/${genome}phi/

# 	rm -r $GENOMEDIR/hgtelements/${genome}/

# 	awk -F"\t" '$10 == "1" { print $0 }' $GENOMEDIR/hgtelements/${genome}phi/prophage_tbl.txt |
# 	awk -v v1="$genome" '{print v1, $0}' FS="\t" OFS="\t" > $GENOMEDIR/hgtelements/${genome}.phispy.txt
# 	rm -r $GENOMEDIR/hgtelements/${genome}phi/

# done

# echo -e "geneID\tannotation\tcontig\tstartpos\tendpos\tgenepos\trank\tstatus\tpp\tfinalstatus" > total.phiSpy.txt
# cat *.phispy.txt >> total.phiSpy.txt 





