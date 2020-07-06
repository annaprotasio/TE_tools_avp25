These instructions describe how to produce annotation of TEs into the genome using a curated library, RepeatMasker and further processing using Onecodetofindthemall and some of my scripts (which are accessible here).  

Summary of steps:

## 1. Have a manually curated library of repeats. 

This should be in FASTA format and the LTR and INT subparts of LTR elements should be places separately. See "Annex 1" for an example. Typically, the manually curated library will have the full elements as FASTA, and this is correct to have. But, for the use of Onecodetofindthemall, it is necessary to run RepeatMasker with the LTR elements split into its subparts. Therefore, it is useful to have a BED file, describing the locations of each LTR subpart and then use `bedtools getfasta` to extract the sequences. Given the full manually curated library (without splitting into subparts) a BED file will look like this:

```
Boudicca#LTR.edit.cons.fa    0    320    Boudicca_LTR#LTR    0    +
Boudicca#LTR.edit.cons.fa    321    4997    Boudicca#LTR    0    +
Fugitive_cat_31_2_31.final.cons.fa    0    315    Fugitive_LTR#LTR    0    +
Fugitive_cat_31_2_31.final.cons.fa    316    4532    Fugitive#LTR    0    +
W2#SINE    0    710    W2#SINE    0    +
Merlin_Sm1m#DNA    0    1113    Merlin_Sm1m#DNA    0    +
RTESj-like_LINE_rnd-1-fam-73.fa    0    3950    RTESj-like#LINE    0    +
```
Note that the 1st column refers to the name of the chromosome in the full manually curated library file. Then:

```
bedtools getfasta -fi full_manually_curated_library.fa -fo june2020_cons_lib.fa -bed june2020_cons_lib.bed -name
```

The output of this, should look like Annex 1. 

## 2. Run RepeatMasker with the curated library

RepeatMasker finds the repeats that are described in the library in the target sequence, in this case, a genome. If you use RepeatMasker please cite the [original work](http://www.repeatmasker.org/faq.html#faq3). Better done in a HPC since it takes about 12 hours in a standard laptop for a ~400 Mb genome (8 cores, mac). The basic command is:

```
RepeatMasker -no_is -nolow -dir ./ -lib <path_to_library> <path_to_genome>
```

The command to run it in Sanger HPC is as follows.

```
module load tetools/1.1
bsub -q normal -R "select[type==X86_64 && mem > 8000] rusage[mem=8000]" -n16 -M8000 -o rpmk.o -e rpmk.e -J rpmk  "RepeatMasker -no_is -nolow -pa 16 -s -poly -excln -dir ./ -lib june2020_cons_lib.fa genome_reference.fa"
```

- no\_is		Skips bacterial insertion element check
- nolow		Does not mask low_complexity DNA or simple repeats
- s			Slow search; 0-5% more sensitive, 2-3 times slower than default
- poly		Reports simple repeats that may be polymorphic (in file.poly)
- excln		Calculates repeat densities (in .tbl) excluding runs of >=20 N/Xs the query
- lcambig	Probably better for nematode/tapeworm species. Outputs ambiguous DNA transposon fragments using a lower case name. All other repeats are listed in upper case. Ambiguous fragments match multiple repeat elements and can only be called based on flanking repeat information. (not used)




## 3. Onecodetofindthemall - Build dictionary step.

This set of perl scripts described in [Bailly-Bechet et al 2014](https://mobilednajournal.biomedcentral.com/articles/10.1186/1759-8753-5-13) is extremely useful to parse RepeatMasker outputs. Please refer to the article and site the original paper if you use any of the `Onecodetofindthemall` scripts. See the [README](https://github.com/annaprotasio/TE_tools_avp25/blob/master/Onecodetofindthemall/README) but basically, this script takes the output of RepeatMasker and matches the internal subpart and LTR subpart of a given LTR element, hence why they were separated in the initial curated library FASTA file. The <build> file generated only has information for the LTR elements

```
perl Onecodetofindthemall/build_dictionary.pl --rm genome_reference.fa.out  > genome_reference.fa.out.build
```

## 4. Onecodetofindthemall - Parse RepeatMasker Annotation.

This script is going to generate the annotation. Why is this better than just the RepeatMasker output? Main points: 
- for each feature found in the genome, it includes the percentage of the reference element that the given feature has. This is important when wanting to remove small fractions from the elements found in the manually curated 'reference' library. 
- builds 'complete' LTR elements, by joining LTR subparts with the corresponding internal subparts if they are within a given number of bases, the default being ~1kb. 
- incorporates a filtering for strict matches, this means that it will only report features that are > 80bases long and that share 80% identity with the reference element. Alongside downstream filtering for 80% of the length if the reference (reported in A), it is possible to apply the 80-80-80 rule. 

This script has dependencies. In particular, there is the length file which needs to be generated in order to run this step. The best why to do this is to use the `samtools faidx` tool on the file generated to run RepeatMasker, that is, the one that looks like Annex 1, and keep the first and second columns of the output.

```
samtools faidx june2020_cons_lib.fa
```
generates a `split_LTR_library.fa.fai`, from which we extract 1st and 2nd columns. 

```
cut -f 1,2 june2020_cons_lib.fa.fai > june2020_cons_lib.length
```

Other parameters can be added, full list of options can be found [here](https://github.com/annaprotasio/TE_tools_avp25/blob/master/Onecodetofindthemall/README) 

To run the parser:

```
perl Onecodetofindthemall/one_code_to_find_them_all.pl --rm genome_reference.fa.out --ltr genome_reference.fa.out.build  --length split_LTR_library.length
```

The output of the script is a series of CSV files, three per chromosome, with detailed annotation of all elements found. The summary files are *elem_sorted.csv and the summary lines in these files start with `###`. For example: 

A perfect, fully assembled Copia element:

```
###2442/41712/2442 0.089 0.000 0.000 chr2L 2294097 2299243 5145 C Copia_LTR LTR/Copia NA NA NA 845/846/847 3 1.000 
2442 0.0 0.0 0.0 chr2L 2294097 2294372 276 C Copia_LTR LTR/Copia (0) 276 1 845 1 
41712 0.1 0.0 0.0 chr2L 2294373 2298967 4593 C Copia_I LTR/Copia (0) 4593 1 846 1 
2442 0.0 0.0 0.0 chr2L 2298968 2299243 276 C Copia_LTR LTR/Copia (0) 276 1 847 1 
```

The top line is the summary/assembly of the three elements underneath. The last value in the ### line is the percentage_reference. 

## 5. Generate valid GFF with final annotation

In order to obtain an annotation file that can be used in standard counting algorithms such as `featureCounts` or `htseq-counts`, it is necessary to produce a GFF. This can be done with `https://github.com/annaprotasio/TE_tools_avp25/blob/master/onecode-2-GFF.pl`. This script takes all the *elem\_sorted.csv, combines them, re-defined solo_LTRs so they can be treated separatelly, filters for a user-defined minimum percentage reference (percentage - length - that the feature has with respect to the reference) and puts results in a user-defined file name. 

Example run:
```
perl onecode-2-GFF.pl 0.8 june2020_cons_lib.TE_annot.gff 
``` 

__________
Annex 1 - FASTA file of manually curated library for RepeatMasker annotation and posterior use with Onecodetofindthemall. Notice that for LTR elements, the sequences are split. "Boudicca\_LTR#LTR" represents that LTR portion of the LTR-element while "Boudicca#LTR" represents the INT or internal portion of the element.
```
>Boudicca_LTR#LTR
TGTAGCTGTAAATAATTCCCCAAAATCAATAGCTCATTAAGTGTTCGACATTGGATGCTGACCACGTTGTTTAAGTGGACTTGTTTAACTGAAGGCTTTCGGTCATGCATCCGTACCAGATCACGTTTCAACACGGCGTGGGACACAGCTCCTACGTTTCATTAGCCAGCTTATAGAAAAGGTCTCTCGATATATTGGTAACGAATCAAATATATCTTTCCTGCCTCTTCTCGTCTGACTTCTGATTTCTATCTGACTGGGAACGATCGAATATAGAAGGTGCTACTGGTAGCTAGTTGTAAAACTAGAATATCCGTTGCATCATC
>Boudicca#LTR
TTGGTGAGACCGACGTGATCTGGTACACCAAGTCAATATACTGTGAGTCTGTTCCATTTTGGAATATTATTATTCATTGTTATTCTTTATTTGAATTTTATATATTGTGATATACTTTTTTTCGCGAATTTTTTTTCCCGGATATATTTTAATTAGACATTTTTCGTTATCTATATTACTATGACGGAACACTCACCCAAGGTTTTTAAGTTAAAGTCTATTCCCCCTATTTCGATTCAGTTAACGCCATTTTGGCCCGACAATATCGAGTCCTGGTTCTGCTATGCAGAAGCCGATTTTTGCATGCACGGAATCACGGACTCGCGCACACGTTTCCTCGCGGTAGTTAAGGCGTTACCGCGCGAGTTTAACAGGTACGTAACACCTAGTATGTTTACTAGTGATGTTCCCGAACCTTACGAAACGCTTAAACTATCGATTTTAAAACGCGGAGACCTAACTGATCGACAACGGTTAGACCAACTCCTTAACAATATCGACTTGCAACATGGTTCTGCGACGGACATGTTGCTTAGAATGAGAGAAGTCATCGGCCAACGAACTTTCGATGATGGCCTATTTAGACAACTTTTCTTGTCTAAACT
>Penelope#PLE
AAGCTGTGAAAAATTCAAATACTTCACTTCTAAACATCAAAATCACCCACCTGAGCTACAAATCTTCTCCACCATCTCAAAATCATACGAGCGTAGTACGCTTCTCGTTGTTCGAAAATGGGAGAATTTCGCCAAAAGGCAAGCTTCCACTCGAGAACAACTATATTTTCTCCACGAGTGCATTCGCCAGCATGTTCTGCCAACAAGTGTAAGATACAGACCTCCAATCAACTGCCCATTAGGATGGAAAATAGCTGAGCAAAGTGGAAGAAAAATGGTTCATCTAATGATAACGGACGCACATTACAGGATCCGCAAATACCAACATGTCATAGAGGACCAAAAGAGGAAGCTACAACAGACACTCACACCAGAACACTTCGAAACCTTGACAACAGCCATCGAAACATTATGCAAACGACACAGAGAAAAAAGAAGAAACATCCTTAAAAAGAAAATAATCAAAATCTCAAAACCAACTATAGAAGAAAACACAAACAACAGCGTGGTCAACCTATCAAAACAACCACTAACAAACACAGAAGAACAGCTACTCCGAAAAGGATTAAATTACAATACAAGCGATGCCTCAAAACTAGAATTCTACGCAGCACTAGAGTCATCACTCAAAACTTCCGGAATAAGTGAGGAAACGCAACAGGAAATAAGACAATCTATAGTACCACTAACTCAACGGATGAAGGGCAACAACCAACTGACGTCACAAGAACAAACAGCTTTGAAGAAACTCAGAACCAGAAAAG
```
