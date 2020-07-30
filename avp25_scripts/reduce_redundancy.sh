#!/bin/bash

if [ $# -ne 2 ]
then
    echo -e "\nusage: $0 <input.fasta> <min_identity>\n"
    echo -e "DESCRIPTION: This script takes an alignment in fasta format and prints a consensus with given parameters. It's a wrapper."
    echo -e "OUTPUT:      produces a <consensus.fa>\n"
     
    exit
fi

fasta_in=$1
min_id=$2

source activate rm2

# Remove INT and LTR redundancies

cd-hit-est -i $fasta_in -o $fasta_in.cd-hit-est-c$min_id -c $min_id

# produce ORFs from the reduced file

getorf -reverse N -minsize 300 -find 1 -sequence $fasta_in.cd-hit-est-c$min_id -outseq $fasta_in.cd-hit-est-c$min_id.pep

