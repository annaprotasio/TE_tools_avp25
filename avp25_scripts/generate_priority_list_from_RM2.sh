#!/bin/bash

if [ $# -ne 3 ]
then
    echo -e "\nusage: $0 <RM2output.fa> <genome.fa> <pfam_db_dir> \n"
    echo -e "DESCRIPTION: 	This script runs a little pipeline that collects RT domains from peptide sequences, makes and MAFF alignment, selects informative blocks with GBLOCKS and calculates a tree with IQTREE.\n"

    echo -e "DEPENDENCIES: 	Requires installation of pfam_scan.pl, HMMER (hmmsearch), cd-hit, mafft, gblocks and iqtree\n"

    echo -e "INPUT:       	<RM2output.fa>		output from RM2 often with the ending `rm2.db-families.fa` "
    echo -e "             	<genome.fa>			genome used to predict the library"
    echo -e "             	<pfam_db_dir>		path to the Pfam database directory\n"

    echo -e "OUTPUT:	    A table\n"     

    exit
fi

rmout=$1
genome=$2
pfamdb=$3

# P2 reduce redundancy

cd-hit-est -i $rmout -o cdhit.fa -c 0.8

# P1 extract info from headers from cd-hit-output

perl rm2_fams2table.pl cdhit.fa # makes cdhit.fa.tab

sort cdhit.fa.tab | awk '{OFs="\t"; $NF=""; print $0}' > col1.txt

# P4 obtain sequence length

samtools faidx cdhit.fa

sort cdhit.fa.fai | awk '{print $2}' > col2.txt

sort cdhit.fa.fai | awk '{print $1}' > fam_names.txt

# P5 blast hits

makeblastdb -in $genome -dbtype nucl

blastn -query cdhit.fa -db chr_dummy.fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" | awk '{OFS="\t"; if ($3 >= 80 && (($4/$13) > 0.5 ) ) {print $0,$4/$13}  }' > genome.blast.o

cat fam_names.txt genome.blast.o | sed 's/\#/ /g' | awk '{print $1}' | sort | uniq -c | awk '{print $1-1}' > col3.txt

# P3 predict domains with Pfam

getorf -sequence cdhit.fa -outseq cdhit.orf -minsize 300

pfam_scan.pl -fasta cdhit.orf -dir $Pfamdb | awk '{if ($6~/^PF/) {print $1}}' | sed 's/\#/ /1' | awk '{print $1}' | sort > pf.domains.count

grep '>' cdhit.fa | sed 's/\#/ /1;s/^>//g' | awk '{print $1}' >> pf.domains.count

cat pf.domains.count | sort | uniq -c | awk '{print $1-1}' > col4.txt


# paste all outputs

paste -d "\t" col1.txt col2.txt col3.txt col4.txt > final_priority.table.tab






