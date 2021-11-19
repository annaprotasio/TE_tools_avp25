# Distance Analysis of RT domains in LTR elements

- List of reference sequences and their accession numbers. We use capital letters for the sequence names to facilitate the recognition in the distance trees later on. 

```
grep ">" ref.seqs.fasta | cut -d " " -f 1,2 | sort | sed 's/>//g;s/$/;/g'
```

BEL AAB03640.1;
BLASTOPIA CAA81643.1;
CER-1 AAA50456.1;
CSRN1 AAK07486.1;
MAG X17219.1\_1;
PAO XP\_028041443.1;
SUSHI TWW62587.1;
TOM CAA80824.1;
TY3 Q7LHG5|YI31B\_YEAST;
ULYSES CAA39967.1;
YOYO Q17318|Q17318_CERCA;

- Note that the *S.mansoni* LTRs are in nucleotide format. We use `getorf` from EMBOSS package to obtain peptide sequences formed by more than 300 nucleotides. We also combine all peptide sequences from both lists

`getorf -sequence LTR_bypart_INT.fasta -outseq LTR_bypart_INT.pep.fa -minsize 300 -reverse n` #code

`cat ref.seqs.fasta LTR_bypart_INT.pep.fa > all.LTRs.seqs.fa` #code


- We are now ready to search for conserved domains using the Pfam script `pfam_scan.pl` and the database of profile-HMMs (see [https://pfam-docs.readthedocs.io/en/latest/ftp-site.html](https://pfam-docs.readthedocs.io/en/latest/ftp-site.html))

`pfam_scan.pl -fasta all.LTRs.seqs.fa -dir ~/Documents/bin/Pfam/Pfam_db/ -outfile all.LTRs.seqs.pfam`#code

- At this point, it is useful to transform the PFAM output into a 6-col BED file. We then extract only the RVT_1 domains by using the PF accession no. 

`pfam_2_bed.sh all.LTRs.seqs.pfam all.LTRs.seqs.pfam.bed` #code
`grep PF00078 all.LTRs.seqs.pfam.bed > RT_01.bed` #code


|TE\_name|start|end|Pfam\_domain|name|strand|
|---------|---|---|----------|----|--|
|CER-1|1071|1231|PF00078.28|RVT\_1|+|
|CSRN1|494|653|PF00078.28|RVT\_1|+|
|SUSHI|581|741|PF00078.28|RVT\_1|+|
|TOM|210|375|PF00078.28|RVT\_1|+|
|TY3|665|823|PF00078.28|RVT\_1|+|
|YOYO|238|403|PF00078.28|RVT\1|+|
|ULYSES|415|514|PF00078.28|RVT\_1|+|
|Aut4\_INT\_1|655|801|PF00078.28|RVT\_1|+|
|Boudicca\_INT\_3|231|390|PF00078.28|RVT\_1|+|
|Fugitive\_INT\_2|552|711|PF00078.28|RVT\_1|+|
|Batuque\_INT\_1|503|649|PF00078.28|RVT\_1|+|
|Oxum\_INT\_1|548|706|PF00078.28|RVT\_1|+|
|Oxossi\_INT\_2|525|684|PF00078.28|RVT\_1|+|
|Columbina\_INT\_1|554|704|PF00078.28|RVT\_1|+|
|Oxala\_INT\_1|503|661|PF00078.28|RVT\_1|+|
|Pierrot\_INT\_1|541|699|PF00078.28|RVT\_1|+|
|Oxaim\_INT\_1|505|663|PF00078.28|RVT\_1|+|
|Aut1\_INT\_2|541|700|PF00078.28|RVT\_1|+|
|Aut3\_INT\_1|557|701|PF00078.28|RVT\_1|+|
|Aut5\_INT\_1|535|695|PF00078.28|RVT\_1|+|
|Aut6\_INT\_1|544|702|PF00078.28|RVT\_1|+|
|Saci2\_INT\_2|523|683|PF00078.28|RVT\_1|+|
|Saci3\_INT\_2|305|464|PF00078.28|RVT\_1|+|
|Saci4\_INT\_2|230|389|PF00078.28|RVT\_1|+|
|Saci5\_INT\_1|464|622|PF00078.28|RVT_1|+|


- This seach did not recognise the RT domains in the Bel/Pao elements. To fix this, we can generate our own HMMs from a sequence alignment and used this custom HMM to find more distantly related RT domains. First, we index the fasta sequence using `samtools faidx` [REF] (otherwise the next step won't work). Second, we pull the parts of the sequences that match the RT domain, using the function `getfasta` from `bedtools` [REF]. Third, we align the sequences with clustal Omega [REF]. 

`samtools faidx all.LTRs.seqs.fa`  #code

`bedtools getfasta -fi all.LTRs.seqs.fa -bed RT_01.bed -fo RT_01.fa` #code

`clustalO -i RT_01.fa -o  RT_01.align.clu` #code

- Typically, the alignment is manually curated / altered to remove uninformative fields. But in this case, there is no need for it. We can then use this alignment `RT_01.align.clu` to generate an profile-HMM using `hmmbuild` [REF].

`hmmbuild RT_custHMM01.stk RT_01.align.clu` #code

- And we can use this new HMM (stored in a stockolm format) to search the peptide sequences again. 

`hmmsearch --domtblout RT_custHMM01.out.doms RT_custHMM01.stk all.LTRs.seqs.fa` #code

- We observe that some new hits appear, such as Saci1, which was not identified with the HMMs from Pfam. There are also some very high e-value matches that are not significant so we can filter them out at the same time that we make our BED format file. 

`grep "RT_01.align" RT_custHMM01.out.doms | awk '{OFS="\t"; if ($12<0.001) {print $1, $20, $21}}' > RT_custHMM01.out.doms.bed` #code

- We now repeat what we did previously: we use the BED file to extract the peptide sequences from the original FASTA, and we align them with clustalO. 

`bedtools getfasta -fi all.LTRs.seqs.fa -bed RT_custHMM01.out.doms.bed -fo RT_custHMM01.out.doms.fasta` #code

`clustalO -i RT_custHMM01.out.doms.fasta -o  RT_custHMM01.out.doms.clu --force` #code

- The alignment looks really good. However, we are still missing 4 Bel/Pao RT domains. Therefore, we iterate the generation of the profile HMM from this alignment, that now, will capture the sequences arising from including Saci1 into the model. 

`hmmbuild RT_custHMM02.stk RT_custHMM01.out.doms.clu` #code

- And we can use this new HMM (stored in a stockolm format) to search the peptide sequences again. 

`hmmsearch --domtblout RT_custHMM02.out.doms RT_custHMM02.stk all.LTRs.seqs.fa` #code

- We inspect the matches and find all Bel/Pao sequences. Note that this was achieved by including only one representative from this family in the `hmmbuild` step. Let's filter again for e-value when we create the BED file.

`grep "RT_custHMM01.out.doms" RT_custHMM02.out.doms | awk '{OFS="\t"; if ($12<0.001) {print $1, $20, $21}}' > RT_custHMM02.out.doms.bed` #code

- Extract the FASTA from the original peptide sequence file and align.

`bedtools getfasta -fi all.LTRs.seqs.fa -bed RT_custHMM02.out.doms.bed -fo RT_custHMM02.out.doms.fasta` #code

`clustalO -i RT_custHMM02.out.doms.fasta -o  RT_custHMM02.out.doms.clu --force` #code

- Now that we are happy with the alignment, we can use a more sensitive aligner and produce a phylip output that we can use to build a tree. 

`source activate rm2` #code
`prank -d=RT_custHMM02.out.doms.fasta +F -o=RT_custHMM02.out.doms.prk` #code

- We use [IQTREE](http://www.iqtree.org/) to build a consensus tree. 

`iqtree -s RT_custHMM02.out.doms.prk.best.phy -m rtREV -B 1000 -alrt 1000 -bnni --prefix 01_tree -redo`

Explained:
`-m rtREV` is the model for reverse transcriptase substitutions.
`-B 1000` specifies the number of ultra fast bootstrap. Note that this approach produces a bootstrap support value that has different interpretation to the traditional bootstrap value (see this [help](http://www.iqtree.org/doc/Frequently-Asked-Questions#how-do-i-interpret-ultrafast-bootstrap-ufboot-support-values) page).
`-alrt 1000` performs the SH-aLRT test, **one would typically start to rely on the clade if its SH-aLRT >= 80% and UFboot >= 95%.**
`-bnni` reduce the risk of overestimating branch supports with UFBoot due to severe model violations.


- We can see that the pHMM starts the RT domain in the 2nd block of conservation (Xiong & Eickbush 1990). It is possible to extend the alignment further up and downstream to include more sequence and look for the blocks of conervation found in RT (as per publication cites above). 

`awk '{OFS="\t"; print $1,$2-40,$3+30}' RT_custHMM02.out.doms.bed > RT_custHMM02.out.doms.mod.bed` 

`bedtools getfasta -fi all.LTRs.seqs.fa -bed RT_custHMM02.out.doms.mod.bed -fo RT_custHMM02.out.doms.mod.fa` 

`clustalO -i RT_custHMM02.out.doms.mod.fa -o RT_custHMM02.out.doms.mod.clu --force` 






