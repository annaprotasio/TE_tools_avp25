#!/bin/bash

if [ $# -ne 1 ]
then
    echo -e "\nusage: $0 <fastq.gz> <ngs_sample_size> <TE_lib.fa> \n"
    echo -e "DESCRIPTION: 	This script uses RepeatMasker to map NGS reads to a TE library.\n"
    echo -e "             	Only reads with at least 80% of their length mapped to the TE are kept.\n\n"

    echo -e "DEPENDENCIES: 	Requires installation of fastq-tools \"https://anaconda.org/bioconda/fastq-tools\" \n"
    echo -e "			 	Requires installation of fastx_toolkit \"https://anaconda.org/bioconda/fastx_toolkit\" \n"  #fastx_toolkit   

    echo -e "INPUT:       	<fastq.gz>    		NGS reads in FASTQ format\n"
    echo -e "             	<TE_lib.fa>   		TE library, ideally manually curated, in FASTA format\n"
    echo -e "             	<ngs_sample_size> 	integer, millions of reads to sample from fastq file\n\n"

    echo -e "OUTPUT:	    shoots a dotplot into X11\n"     

    exit
fi

fastq=$1
sample_size=$2
TE_lib=$3

readsM=$2*1000000
samp_name=$(echo "$fastq" | cut -f 1 -d '.')


# The first step is to sample the fastq file. We use the fastq-sample script from fastq-tools "https://anaconda.org/bioconda/fastq-tools"

echo "sampling $readsM reads from fastq"

fastq-sample -n $readsM -o $samp_name -s 1 $fastq

# change the name of the fastq-sample output so we know it is compressed

mv $samp_name.fastq $samp_name.fastq.gz

# transform fastq into fasta 

fastq_to_fasta -i $samp_name.fastq.gz -o $samp_name.fasta

# Run repeatmasker against the library

RepearMasker -a -nolow -no_is -lib $TE_lib $samp_name.fasta

# Parse RepeatMasker.out

cat $samp_name.fasta.out | awk '{if ($5~/ERR/) { if ($7-$6 > 80) {print $10}}}' | sort | uniq -c | awk '{OFS="\t"; print $2,$1}' > sampe_name.TE.readcount

