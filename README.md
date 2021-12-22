# TE_tools_avp25
Transposable Elements tools collected by Anna 

## Summary of Tools / files

### TE_annotation.md
Is a short guideline of how to annotate repeats in the genome, starting with a manually curated library of repeats. The aim is to obtain a GTF annotation of repeats in the genome that can be used for downstream analyses such as TE expression, small RNA mapping, etc. 

### onecode-2-GFF.pl
Takes the output from Onecodetofindthemall (https://mobilednajournal.biomedcentral.com/articles/10.1186/1759-8753-5-13), which are a set of *.elem_sorted.csv files, and makes a GFF annotation file that can be used to quantify reads from mapping to genome. The output is highly customised for our needs working on the /Schistosoma mansoni/ genome transposable elements annotation. However, this script is also highly "customisable" and the fileds can be changed to suit each user. 

Example use:
`cd path-to-where-I-ran-onecode`

`perl onecode-2-GFF.pl 0.8 out.gff`

### RMout2GTF.pl
My own parser of RepeatMasker *.out file to GTF format. Removes simple repeats and any annotation of less than N base pairs (defined by the user). 

Example use:
`perl RMout2GTF.pl yourspecies_repeatmasker.fa.out 100` # will keep only those annotations of > 100 base paires

