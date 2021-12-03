#!/bin/bash

if [ $# -ne 1 ]
then
    echo -e "\nusage: $0 <LTR.fa> \n"
    echo -e "DESCRIPTION: 	This script uses Blast to find the LTR subparts. Blast needs to be previously installed and found in the path\n"

    echo -e "DEPENDENCIES: 	Requires installation of fastq-tools \"https://anaconda.org/bioconda/blast\" \n"
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

current_time=$(date "+%Y.%m.%d-%H.%M.%S") 

blastn -query $fasta -db $fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"  >  blast.$current_time

# check that the length of the LTR element corresponds to the last coordinate of the 3' LTR

entries=$(awk '{print $1}' < blast.$current_time | wc -l)
unique=$(awk '{print $1}' < blast.$current_time | sort -u | wc -l)

if [ $entries != $unique ]
then
    echo -e "\nMore than 1 possible output for at least one entry or one entry missing from output, please check output manually"
fi    

# parse blast 
echo "###gff-annotation" > LTR_annotation.gff
cat blast.$current_time | awk '{OFS="\t";if ( ($1==$2) &&  ($4<1000) && ($7<10) && ($9>1000) && ($10==$13) ) {print $1,"manual","LTR",$7,$8,".","+",".","Superfamily \""$1 "\"; Family \""$1"\"; Type \"LTR\"; Subpart \"5LTR\";\n"$1,"manual","LTR",$8+1,$9-1,".","+",".","Superfamily \""$1 "\"; Family \""$1"\"; Type \"INT\"; Subpart \"INT\";\n"$1,"manual","LTR",$9,$10,".","+",".","Superfamily \""$1 "\"; Family \""$1"\"; Type \"LTR\"; Subpart \"3LTR\";"}}'>> LTR_annotation.gff

# make bed file

cat LTR_annotation.gff | awk '{OFS="\t"; print $1,$4-1,$5,$1"-"$NF,$6,$7}' |  tr -d '";' | sed '/^#/d' > LTR_annotation.bed

rm blast.$current_time





