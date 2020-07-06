#!/usr/bin/perl -w
use strict;
use warnings;

my $arg_count = @ARGV;
if ($arg_count != 2) {
  &print_usage;
  die "Please specify 2 arguments!";
}

my $perc_ref =  shift;		#minimum perc of reference required for each TE, comes from args
my $outfile  =  shift;		#name of file where output is to be save, comes from args

open( OUT, ">$outfile");
print OUT "\#\#\#gff-annotation.\n";


for my $file (glob '*.elem_sorted.csv') {
    open my $FH, '<', $file or die "$file: $!";
    open( OUT, ">>$outfile");
    while (<$FH>) {
        next if /^\n/ ;
        next unless /^\#/ ;
        my @line = split /\s+/, $_ ;
        
        # filter for percentage identity
        next if ($line[16] < $perc_ref);
       
        # redefine the strand annotation as +/-
        my $strand = "";
    	if ( $line[8] eq 'C') {
		$strand = "-";
    	}
    	else {$strand = "+";
    	}    	
    	
    	# define if the feature is a solo LTR or not
    	my $element = "";
    	if ( 
    		($line[9] =~ /.+LTR/ ) && 		# is the element field, such as Boudicca_LTR 
    		($line[10] eq 'LTR') &&			# is the Order fiels, such as LTR, LINE, CR1, etc. 
    		($line[15] == 1)				# this field described the number of features used to make the element. If == 1, plus LTR, then it's a solo LTR
    		) { $element = $line[9] . "_solo"
    		} else {
    		$element = $line[9]
    		}
    	print OUT "$line[4]\tOneCode\trepeat_region\t$line[5]\t$line[6]\t$line[1]\t$strand\t.\trepeat_id \"$element\"; Order \"$line[10]\"; Perc_ref \"$line[16]\";\n";  	 
        
        
    }
}
close OUT or die $!;


sub print_usage {
  
  print "Usage: $0 <perc_ref> <file_GFF>\n\n";
  print "Parser script. Takes the output from \"Onecodetofindthemall\" (https://mobilednajournal.biomedcentral.com/articles/10.1186/1759-8753-5-13) and makes a GFF/GTF that can be used for quantification. This script discriminates between soloLTRs, providing a distict feature id.\n\n";
  print "IMPORTANT:    $0 must be called from the same directory where the output of Onecode is, as it will look for *.elem_sorted.csv files ... hence, there is no unpit declaration in the command callin\n\n";
  print "INPUT:        <perc_ref>            float - for each TE feature, Onecode reports the percentage of the reference. This value sets the minimim. For strict parameter, use 0.8\n";
  print "INPUT:        <file_GFF>            Name of GFF file where the script should save the annotation. The file is standard GTF/GFF with the 9th field having bespoke annotation\n\n";
  print "Authored by Anna V. Protasio. Report bugs and constructive feedback to avp25\@cam.ac.uk. GNU GENERAL PUBLIC LICENSE. \n\n";
  
  }

exit