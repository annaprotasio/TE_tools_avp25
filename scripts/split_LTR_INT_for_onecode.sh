#!/bin/bash

if [ $# -ne 1 ]
then
    echo -e "\nusage: $0 <LTR.fa> \n"
    echo -e "DESCRIPTION: 	This script uses Blast to find the LTR subparts. Blast needs to be previously installed and found in the path\n\n!!! THIS SCRIPT COMES WITH ABSOLUTLEY NO GUARANTEES - PLEASE CHECK YOUR OUTPUT AGAINST YOUR INPUT TO MAKE SURE IT BEHAVES AS EXPECTED !!!\n"

    echo -e "DEPENDENCIES: 	Requires installation of BLAST+ \"https://anaconda.org/bioconda/blast\" \n"
    echo -e "INPUT:       	<LTR.fa>    		A FASTA file with the LTR elements sequences"
    echo -e "OUTPUT:	    writes a GFF to output\n"     

    exit
fi

fasta=$1

#samp_name=$(echo "$fastq" | cut -f 1 -d '.')


# The first step is to make the self-blast

# check if database has been made


if [ -f "$fasta.nin" ]; then
    echo "$fasta.nin exists."
else 
    echo "$fasta.nin does not exist. Building database"
    makeblastdb -in $fasta -dbtype nucl
fi

echo "running self-blast"

# get current time so subsequent blast results don't overwrite each other
current_time=$(date "+%Y.%m.%d-%H.%M.%S") 

# perform blast - `qlen` is added as last field because we use it for comparison later on. 
blastn -query $fasta -db $fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"  >  blast.$current_time

# parse blast 
echo "###gff-annotation" > LTR_annotation.gff
cat blast.$current_time | awk '{OFS="\t";
	if ( ($1==$2) &&  	 # query name == subject name
	($4<1000) && 		 # match length < 1000; finds the LTR subpart 
	($7<20) && 			 # qstart < 20; start of 5LTR subpart match
	($9>1000) && 		 # sstart > 1000; matches start of 3LTR
	($10 >= ($13-10))) 	 # send >= (query length - 10nt)
	## build GFF
	{print $1,"manual","LTR",1,$8,".","+",".","Superfamily \""$1 "\"; Family \""$1"\"; Type \"LTR\"; Subpart \"5LTR\";\n"$1,"manual","LTR",$8+1,$9-1,".","+",".","Superfamily \""$1 "\"; Family \""$1"\"; Type \"INT\"; Subpart \"INT\";\n"$1,"manual","LTR",$9,$10,".","+",".","Superfamily \""$1 "\"; Family \""$1"\"; Type \"LTR\"; Subpart \"3LTR\";"}}'>> LTR_annotation.gff

# make bed file

cat LTR_annotation.gff | awk '{OFS="\t"; print $1,$4-1,$5,$1"-"$NF,$6,$7}' |  tr -d '";' | sed '/^#/d' > LTR_annotation.bed

#rm blast.$current_time

# report summary to user. 

SEQS_IN=$(grep -c ">" $fasta)
SEQS_OUT=$(awk '{print $1}' LTR_annotation.bed | sort -u | wc -l)

if [ $SEQS_IN != $SEQS_OUT ]
then
    echo -e "\nWARNING: There are $SEQS_IN named elements in your input file BUT $SEQS_OUT named elements in the output. Please check manually"
    else
    echo -e "\nSUCCESS! There are $SEQS_IN named elements in your input file AND $SEQS_OUT named elements in the output."
fi    






