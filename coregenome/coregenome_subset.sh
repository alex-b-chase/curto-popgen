#!/bin/bash


usage(){
	echo "$(basename "$0") [-h] 

This program to take output from ROARY and subset core identified core genomes from outdirXX


In other words, you cannot have .fa files that represent both nucleotide genomes and translated amino acid genomes

correct use:
	coregenome_subset.sh -i minID 

where:
	--help or -h  			show this help text
	--identity or -i 		min BLASTp ID that ROARY generated (output files must be named 'outdir{minID}')
"
}


# get user defined options
if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit 1
fi


# get user input for other parameters
while test $# -gt 0; do
		case "$1" in
				-i)
						shift
						if test $# -gt 0; then
								export IDENTITY=$1
						else
								echo "no identity specified"
								exit 1
						fi
						shift
						;;
				--identity*)
						export IDENTITY=`echo $1 | sed -e 's/^[^=]*=//g'`
						shift
						;;
				*)
						break
						;;
		esac
done




### going to subset core genome genes from Curtobacterium elevation genomes + Frigoribacterium outgroup

BASE=/Volumes/Bennett_BACKUP/Research/curtobacterium/curto-population-genomics
BASEDIR=$BASE/roary_coregenome
GENOMEDIR=$BASE/genomes

### use the XX% BLASTp minID based on the results in: 
### $BASEDIR/roary_ortho_results.xlsx

# get number of genomes for all
numgen=$(ls $GENOMEDIR/*.gff | wc -l)

minID=$IDENTITY 

# rm -f $BASEDIR/coregenes${minID}.txt

# ################################################################################################
# # go through the core genes output from ROARY
# # create reference files for each core gene across all genomes
# # output the results into a gene-specific reference file
# ################################################################################################

# csv2txt.py < $BASEDIR/coregenome${minID}/gene_presence_absence.csv > $BASEDIR/coregenes${minID}.temp.txt

# head -n1 $BASEDIR/coregenes${minID}.temp.txt | cut -f1,3,4,5,15- > $BASEDIR/coregenes${minID}.txt

# # check that isolates and sequences with that gene are the same and export those genes only
# awk -F"\t" -v var="$numgen" '(NR>1) && ($4 == var ) ' $BASEDIR/coregenes${minID}.temp.txt | \
# awk -F"\t" -v var="$numgen" '(NR>1) && ($5 == var ) ' | \
# cut -f1,3,4,5,15- >> $BASEDIR/coregenes${minID}.txt

# cut -f1,2 $BASEDIR/coregenes${minID}.txt > $BASEDIR/coregeneIDs${minID}.txt

# rm -f $BASEDIR/coregenes${minID}.temp.txt


# # transpose the file
# cut -f5- $BASEDIR/coregenes${minID}.txt | tr -d '"' | awk '
# { 
#     for (i=1; i<=NF; i++)  {
#         a[NR,i] = $i
#     }
# }
# NF>p { p = NF }
# END {    
#     for(j=1; j<=p; j++) {
#         str=a[1,j]
#         for(i=2; i<=NR; i++){
#             str=str" "a[i,j];
#         }
#         print str
#     }
# }' | tr -d '\r' > $BASEDIR/coregenes${minID}.ids.txt

# rm -rf $BASEDIR/tmp
# mkdir -p $BASEDIR/tmp

# while read genome genes
# do

# 	echo $genes | tr ' ' '\n' > $BASEDIR/tmp/$genome.txt

# 	echo -ne "starting with ${genome}"

# 	while read geneID
# 	do

# 		gene=$(grep -w ${geneID} $BASEDIR/coregenes${minID}.txt | cut -f1)

# 		if [[ $geneID == *".p01" ]]
# 		then

# 			locusID=$(echo $geneID | cut -f1 -d'.')

# 		elif [[ $geneID == *".r01" ]]
# 		then

# 			locusID=$(echo $geneID | cut -f1 -d'.')

# 		else
# 			locusID=$(grep "ID=${geneID};" $BASEDIR/coregenome${minID}/fixed_input_files/${genome}.gff | \
# 			grep -o 'Parent=.*' | cut -f2 -d'=' | cut -f1 -d'.')

# 		fi

# 		echo -e "${genome}\t${locusID}" >> $BASEDIR/tmp/$gene.ids.txt

# 	done < $BASEDIR/tmp/$genome.txt
	

# 	echo -e "\t - done!"

# done < $BASEDIR/coregenes${minID}.ids.txt



echo ""
echo "Starting the organization by each core gene"
echo ""


################################################################################################
# go through the reference ids.txt files to find the genome file and sequence corresponding
# to the AA sequence of the core genes
# output the results into a core gene specific AA file with all seqs in there
################################################################################################

rm -rf $BASEDIR/coregenes${minID}
mkdir -p $BASEDIR/coregenes${minID}

cd $BASEDIR/tmp/

find . -type f -name "*.ids.txt" -exec mv {} $BASEDIR/coregenes${minID} \;

# rm -rf $BASEDIR/tmp/

cd $BASEDIR/coregenes${minID}

# # remove old files before beginning
# rm -f *.total.aln
# rm -f *.total.fa 
# rm -f *.fasta 
# rm -f *.faa 
# rm -f *.temp.faa
# rm -f *.temp.txt

# # rm -f $BASE/removed_genes.txt

# loop through the lines with the hit core genes and subset by the specific core gene
for f in *.ids.txt
do
	filename=${f%.ids.txt}

	echo -ne "Processing core gene ${filename}\t"

	# check to make sure these are not paralogs
	check=$(cut -f1 $f | sort | uniq -c | grep -v '^ *1 ')

	if [ -z "$check" ]
		then
			:
		else 
			echo "$filename is a paralog. Deleting..."
			echo $filename >> $BASEDIR/paralogs.txt 

			rm $f 
			break
	fi

	echo -ne "...subsetting\t"

	while read genome sequence 
	do
		genfile=${genome}.ffn
		# genfile=${genome}.faa
		sequence2=$(echo $sequence | cut -f1 -d' ')

		# some genomes may not exist for some reason or were not included 
		# so if genome was deleted, do nothing (:)
		if [ -f $GENOMEDIR/$genfile ]
			then 
				# awk -v "seq=$sequence2" -v RS='>' '$1 == seq {print RS $0}' $BASE/prokka/$genome/$genfile | \
				# sed "s/[>].*$/>$genome/" >> $BASEDIR/coregenes${minID}/$filename.ids.faa

				echo ${sequence2} > $BASEDIR/coregenes${minID}/$filename.$genome.temp.txt

				filterbyname.sh in=$GENOMEDIR/$genfile \
				out=$BASEDIR/coregenes${minID}/$filename.$genome.temp.ffn \
				include=t names=$BASEDIR/coregenes${minID}/$filename.$genome.temp.txt \
				substring=name ow=t > /dev/null 2>&1

				sed "s/[>].*$/>$genome/" $BASEDIR/coregenes${minID}/$filename.$genome.temp.ffn >> $BASEDIR/coregenes${minID}/$filename.ids.ffn
			else
				:
		fi

		rm $BASEDIR/coregenes${minID}/$filename.$genome.temp.ffn
		rm $BASEDIR/coregenes${minID}/$filename.$genome.temp.txt

	done < $f

	cat $BASEDIR/coregenes${minID}/$filename.ids.ffn | tr -d "[ -%,;\(\):=\.\\\[]\"\']" | sed "s/\*//g" > $BASEDIR/coregenes${minID}/$filename.total.fas 

	rm $BASEDIR/coregenes${minID}/$filename.ids.ffn

	echo -ne "...aligning\t"

	check2=$(cat $BASEDIR/coregenes${minID}/$filename.total.fas | wc -l)
	if [ $check2 -gt $numgen ]
		then
			# align each gene
			clustalo -i $BASEDIR/coregenes${minID}/$filename.total.fas -o $BASEDIR/coregenes${minID}/$filename.total.ffn.aln --force

		else
			echo $filename >> $BASE/removed_genes.txt
			rm -f $BASEDIR/coregenes${minID}/$filename.*
			echo -e "- removed!"
	fi

	echo -ne "...checking\t"

	if [ -f $BASEDIR/coregenes${minID}/$filename.total.ffn.aln ]
		then

		check3=$(grep '>' $BASEDIR/coregenes${minID}/$filename.total.ffn.aln | wc -l)

		if [[ $check3 == $numgen ]]
			then
				:
			else 
				echo $filename >> $BASE/removed_genes.txt
				rm -f $BASEDIR/coregenes${minID}/$filename.*
				echo -e "- removed!"
		fi
	fi

	# rm $BASEDIR/coregenes${minID}/$filename.ids.faa
	echo -e "- done!"

done

# rm *.ids.txt

################################################################################################
# build a consensus reference tree by combining all 30 marker genes into one file
# construct the tree
################################################################################################

# combine the fasta files together for tree construction keeping the proteins correctly ordered
# all individual alignment files must have the same header

# too many files to do at once so break into pieces
# for letter in {a..z} 
# do
# 	echo $letter
# 	catfasta2phyml.pl -f ${letter}*.total.aln > ${letter}.total.sub.aln

# 	if [[ -s ${letter}.total.sub.aln ]]
# 	  then
# 	  	:
# 	  else
# 	  	rm -f ${letter}.total.sub.aln
# 	fi
# done

# echo "done!"

catfasta2phyml.pl -f *.ffn.aln > concat.aligned.coregenes${minID}.fa

echo "starting tree building!"

# build the tree using the aligned file and output 100 bootstrap replications
# kept crashing, running on the other lab computer - WORKED so run separately!

rm -f RAxML*
rm -rf $BASEDIR/coregenes${minID}/master_tree/
mkdir $BASEDIR/coregenes${minID}/master_tree/

mv concat.aligned.coregenes${minID}.fa $BASEDIR/coregenes${minID}/master_tree/
cd $BASEDIR/coregenes${minID}/master_tree/


# try and run on HPC, no way will this work on desktop
raxmlHPC-PTHREADS -s concat.aligned.coregenes${minID}.fa -m PROTGAMMAWAG -n coregenes${minID} -x 100 -# 100 -p 4321 -f a -T 8



















