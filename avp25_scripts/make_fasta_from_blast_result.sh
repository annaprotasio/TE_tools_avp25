#!/bin/bash


if [ $# -ne 4 ]
then
    echo -e "\nusage: $0 <genome.fa> <fasta.in> <min_length> <flank>\n"
    echo -e "DESCRIPTION: This script takes a fasta sequence <fasta.in>, blasts it to the <genome.fa>, recovers locations with alingment lenght > <min_length>, prints them as bed file, extends bed coordinates <flank> bases in each direction,  and makes fasta from that BED."
    echo -e "OUTPUT    :  produces a <fasta.in.bed> which is the blast results BED file file"
    echo -e "          :  produces a <fasta.in.blast.flank.bed> which is the extended BED"
    echo -e "          :  produces a <fasta.in.blast.flank.bed.fa> bedtools getfasta \n"
    echo -e "REQUIRES  :  file with chromosome lengths with same name as genome, as in <genome.fa.length>, and found in the same dir as genome (use \"samtools faixd\" and then \"cut -f 1,2\""
    echo -e "          :  <genome.fa> must be formatted for blast with \"makeblastdb\" "		
     
    exit
fi

genome=$1
fasta_in=$2
out=`basename $2`
min_length=$3
flank=$4

blastn -query $fasta_in -db $genome -outfmt "6 qseqid sseqid pident length mismatch qstart qend sstart send sstrand" -evalue 0.0000000000000000000000000001 > $out.blast.o

awk -v "ml=$min_length" '{OFS="\t"; if ($4 > ml ) { if ($10~/plus/) {print $2, $8, $9, $1, $3, "+"} else {print $2, $9, $8, $1, $3, "-"}}}' < $out.blast.o > $out.blast.bed

bedtools slop -s  -i $out.blast.bed -g $genome.length -b $flank > $out.blast.flank.bed

bedtools getfasta -fi $genome -fo $out.blast.bed.fa  -bed $out.blast.flank.bed -s 

