#!/bin/bash


if [ $# -ne 1 ]
then
    echo -e "\nusage: $0 <file.gff> \n"
    echo -e "DESCRIPTION: This script takes a GFF file and makes an indexed file. This can be used in various downstream analyses, for example, display GFF annotation in Artemis. "
    echo -e "OUTPUT:      produces a <file.sorted.gff.gz"
    echo -e "             produces a <fasta.in.blast.flank.bed.fa> bedtools getfasta \n"
    echo -e "REQUIRES:    tabix installation"
    		
     
    exit
fi

gff_in=$1
name=$(echo "$gff_in" | cut -f 1 -d '.')

(grep ^"#" $gff_in; grep -v ^"#" $gff_in | sort -k1,1 -k4,4n) | bgzip > $name.sorted.gff.gz;
tabix -p gff $name.sorted.gff.gz;