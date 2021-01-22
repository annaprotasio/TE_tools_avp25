# Distance Analysis of RT domains in PLE elements

- After extensive manual curation of the prospective PLE element, the list of remaining sequences are:

Acamarachi
Aracar
Arizaro
Cercyon
Chachacomani
Inca
Licancabur
Mururata
Palpana
Perere-0
Tocorpuri
Uturunku
Wallatani
Yaypuri


- List of reference sequences and their accession numbers. We use capital letters for the sequence names to facilitate the recognition in the distance trees later on. 

```
>Penelope_Dv AAA92124.2 ORF1 [Drosophila virilis] (from U49102)
>Poseidon?_Od AAT48673.1 pol protein [Oikopleura dioica] (from AY634216)
>FOLCA A0A226F598 (from Pfam db)
>Xena_Tr AAK58879.1 putative reverse transcriptase [Takifugu rubripes]
>XP_031758619.1 uncharacterized protein LOC116410997 [Xenopus tropicalis]
>XP_008123603.1 PREDICTED: uncharacterized protein LOC103282659 [Anolis carolinensis]
```

- Additional sequences to somehow root the tree later:
```
>TERT_RAT sp|Q673L6|TERT_RAT 
>TERT_YEAST sp|Q06163|TERT_YEAST 
```


- We are now ready to search for conserved domains using the Pfam script `pfam_scan.pl` and the database of profile-HMMs (see [https://pfam-docs.readthedocs.io/en/latest/ftp-site.html](https://pfam-docs.readthedocs.io/en/latest/ftp-site.html))

```
cat PLE.Sman.cdhit.tseq.pep reference.ple.pep telomerases.pep > sman.ref.telo.pep #code

pfam_scan.pl -fasta sman.ref.telo.pep -dir ~/Documents/bin/Pfam/Pfam_db/ -outfile sman.ref.telo.pep.pfam #code

```

- At this point, it is useful to transform the PFAM output into a 6-col BED file. We then extract only the RVT_1 domains by using the PF accession no. 

`pfam_2_bed.sh sman.ref.telo.pep.pfam > sman.ref.telo.pep.pfam.bed ` #code
`grep RVT_1 sman.ref.telo.pep.pfam.bed > sman.ref.telo.RVT_1.bed`

- This seach did not recognise the RT domains in the **XP\_031758619 and XP\_008123603 or  the yeast telomerase** elements. To fix this, we can generate our own HMMs from a sequence alignment and used this custom HMM to find more distantly related RT domains. First, we index the fasta sequence using `samtools faidx` [REF] (otherwise the next step won't work). Second, we pull the parts of the sequences that match the RT domain, using the function `getfasta` from `bedtools` [REF]. Third, we align the sequences with clustal Omega [REF]. 

```
samtools faidx sman.ref.telo.pep   #code
bedtools getfasta -fi sman.ref.telo.pep -bed sman.ref.telo.RVT_1.bed -fo sman.ref.telo.RVT_1.pep #code
mafft sman.ref.telo.RVT_1.pep > sman.ref.telo.RVT_1.pep.maf #code
```

It is clear from the alignment that we need to use more sequence from the flanks of the entries. Let's add 30 AA each way:

```
cat sman.ref.telo.pep.fai | awk '{OFS="\t"; print $1,$2}' > sman.ref.telo.pep.length #code
bedtools slop -i sman.ref.telo.RVT_1.bed -b 30 -g sman.ref.telo.pep.length > RVT1_f30.bed #code
bedtools getfasta -fi sman.ref.telo.pep -bed RVT1_f30.bed -fo RVT1_f30.pep #code
mafft  RVT1_f30.pep >  RVT1_f30.pep.maf #code
```

- The alignment is manually edited so it includes only the sequence boundaries described in the RT alignment of figure 2 from Arkipova 2006. This is used to generate own RT hmm.

`hmmbuild ownRTf30.stk RVT1_f30.pep.maf` #code

- And we can use this new HMM (stored in a stockolm format) to search the peptide sequences again plus the telomere sequences.


```
hmmsearch --domtblout sman.ref.telo.ownRTf30.out ownRTf30.stk sman.ref.telo.pep #code
hmmsearch_2_bed.sh sman.ref.telo.ownRTf30.out > sman.ref.telo.ownRTf30.out.bed
```

This approach finds all the RT domains in all the sequences. We can now build an alignment to analyse their phylogenetic distances, although two need to be removed :

Wallatani#PLE\_1	66	101	RVT1\_f30.pep	3.7e-11	+
Chachacomani#PLE\_2	174	207	RVT1\_f30.pep	7.2e-11	+


```
samtools faidx sman.ref.telo.pep
bedtools getfasta -fi sman.ref.telo.pep -bed sman.ref.telo.ownRTf30.out.bed -fo sman.ref.telo.ownRTf30.out.pep
mafft sman.ref.telo.ownRTf30.out.pep > sman.ref.telo.ownRTf30.out.pep.maf
```

We can also make an alignment with PRANK:
```
source activate rm2 #code
prank -d=sman.ref.telo.ownRTf30.out.pep +F -o=sman.ref.telo.ownRTf30.out.pep.prk #code
```

- We use [IQTREE](http://www.iqtree.org/) to build a consensus tree. 

`iqtree -s sman.ref.telo.ownRTf30.out.pep.prk.best.phy -m rtREV -B 1000 -alrt 1000 -bnni --prefix 01_tree`

Explained:
`-m rtREV` is the model for reverse transcriptase substitutions.
`-B 1000` specifies the number of ultra fast bootstrap. Note that this approach produces a bootstrap support value that has different interpretation to the traditional bootstrap value (see this [help](http://www.iqtree.org/doc/Frequently-Asked-Questions#how-do-i-interpret-ultrafast-bootstrap-ufboot-support-values) page).
`-alrt 1000` performs the SH-aLRT test, **one would typically start to rely on the clade if its SH-aLRT >= 80% and UFboot >= 95%.**
`-bnni` reduce the risk of overestimating branch supports with UFBoot due to severe model violations.

Tree looks pretty good. 

Let's edit the names of the entries so the tree looks even prettier!


