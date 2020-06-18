# TE_tools_avp25
 Transposable Elements tools

## Summary of Tools

### onecode-2-GFF.pl
Takes the output from Onecodetofindthemall (https://mobilednajournal.biomedcentral.com/articles/10.1186/1759-8753-5-13), which are a set of *.elem_sorted.csv files, and makes a GFF annotation file that can be used to quantify reads from mapping to genome. The output is highly customised for our needs working on the /Schistosoma mansoni/ genome transposable elements annotation. However, this script is also highly "customisable" and the fileds can be changed to suit each user. 

Example use:
`cd path-to-where-I-ran-onecode`
`perl onecode-2-GFF.pl 0.8 out.gff`