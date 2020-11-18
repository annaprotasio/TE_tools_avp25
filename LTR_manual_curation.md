# Manual curation of LTR retrotransposons

LTR retroelements are formed by three parts: the 5' LTR subpart, the internal subpart and the 3' LTR subpart. In the vast majority of cases, both LTR subparts are very similar to each other, or perhaps identical. This characteristic in the anatomy of the LTR retroelements allow us to use bioinformatic approached to find their location, based on the similarity of the termini. In addition, the internal part, in most cases, will have protein coding potential, and will encode for conserved domains that are required to fulfil the life cycle of the element. More introductory information on LTR retrotransposons can be found in [Wikipedia](https://en.wikipedia.org/wiki/LTR_retrotransposon).

In my experience, I started manual curation of TEs (transposable elements) with LTR retroelements, because they are probably the easiest. There is the need for some coding skills, but the scripts you will required are already written (inser link). For an introduction on UNIX command line see this excellent [Unix/Linux tutorial for beginners](http://www.ee.surrey.ac.uk/Teaching/Unix/index.html) published by the University of Sussex and this other tutorial on the use of the three horses of the apocalypse: [AWK, SED and GREP](https://www-users.york.ac.uk/~mijp1/teaching/2nd_year_Comp_Lab/guides/grep_awk_sed.pdf) published by the University of York. 


## Software you will need.

In order to assist with the different steps of the manual curation, there are a number of softwares that come in handy. Here is the list:

| program  | run mode  |  example |  use | need |
|---|---|---|---|
|  dotplot  | GUI  | gepard   |  pairwise alignment | essential |
| blast  | commandline  | blast+  | homology search  | essential |
| bedtools | commandline | bedtools | extract sequences from genome and multi-task | essential |
| Alignment editor | GUI  | AliView  | view and edit an alignment  |  essential |
| various sequence manipulation tools | commandline | EMBOSS | multi-task | essential |
| genome viewer | GUI | Artemis | genome annotation | optional | 
| conserve protein domain finder | commandline/ web | Pfam_scan.pl/interproscan | finds protein domains in ORFs | optional |


## Steps

In a nutshell, the general steps for doing manual curation of LTR retroelements is as follows and can be replicated for each predicted LTR family. 

1. use consensus sequence to find locations in the genome where the consensus is found (recommended program, blast or RepeatMasker output)
2. based on the location of the hit, extend the coordinates as desired, up and downstream of the hit, to capture the LTR subparts (various `bedtools` programs)
3. align the sequences in an alignment viewer and edit the alignement (AliView)
4. produce a consensus sequence of the alignment (EMBOSS `cons` program)
5. check for the appearance of LTR subparts (Dotplot - gepard or other). 
6. repeat steps 2-5 using the new consensus and extending to the blast hit to include more flanking sequence until the LTRs are found. 
7. translate the new consensus and look for the long ORF. 
8. look for conserved protein coding domains in the ORF.
9. visualise the LTR retroelement and annotate its features (optional but recommended). Once all are done, it will be possible to collect all annotations into one GFF that can be used in many ways. 
10. reduce redundancy and find related LTR retroelements using `cd-hit-est`. Collect all sequences for LTR and INT subparts (using my scrip `convert_embls_to_gff_and_bed.pl` and `bedtools getfasta`) and run `cd-hit-est -i <in.fasta> -o cd-hit-est_output.o -c 0.80`, note to include the identity threshold of 0.80, as described in the 80-80-80 rule (add REF)





## In detail 

After getting the RepeatModeler2 output, there will be a list of all the families that were identified. The first step is to look at their names and their classification:

`grep ">" your_org-families.fa`

You can then extract the names of the families that are classified as LTR. Notice that there are two groups of consensus sequences and are tagged either as "Type=LTR" or "Type=INT" and they refer to the subpart that they represent. 

`grep "\#LTR" your_org-families.fa > LTR-families.ids`

I like to have to have a CSV (comma separated value) file with the names and the classification (INT/LTR) so I can prioritise and keep notes:

`grep "\#LTR" your_org-families.fa | grep Type | awk '{print $1","$3}' | sed 's/>//g' > ltr-families.csv`

This `ltr-families.csv` can be opened and edited in Excel or Google docs. 



## Alternative protocols

1. in order to reduce or prioritise families for manual curation, it is possible to predict conserved domain in bulk. For that, for each conserved family, predict the translation of in all 6 frames `transeq -sequence <in_fasta> -outseq <out_pep> -frame 6` and then predict any conserved domains (for example using `pfam_scan.pl` from the [EBI](https://www.ebi.ac.uk/Tools/pfa/pfamscan/), also available through `conda`). It is then possible to select only those sequences that have coding potential for an LTR retroelement conserved domain. 

