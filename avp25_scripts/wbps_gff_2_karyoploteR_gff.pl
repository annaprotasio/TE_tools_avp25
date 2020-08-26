#!/usr/bin/perl
#created by Anna V. Protasio (University of Cambridge).
#last modified 12.08.2020


use strict;
use warnings;

my $arg_count = @ARGV;
my $corr_args = 1;
if ($arg_count != $corr_args) {
  &print_usage;
  die "Please specify $corr_args arguments!";
}

my $input = shift; #the the gff file

open( INFILE, $input ) or die( "Couldn't open $input: $!.\n" ) ;

while (<INFILE>) {
	chomp;
	my $gene_id = "";
	my $transcript_id = "";
	my @line = split (/\s+/, $_);
	if (@line > 8) {
		if ($line[2] =~ m/mRNA/) { # ID=transcript:Acey_s0009.g673.t1;Parent=gene:Acey_s0009.g673;
			if ($line[8]=~m/ID=transcript:(Acey_s\d{4}.g\d*\.t\d*);Parent=gene:(Acey_s\d{4}.g\d*);/g) {
				$gene_id = $2 ; 
				$transcript_id = $1; 
#				print $2 . "\n";
				print $line[0]."\t".$line[1]."\ttranscript\t".$line[3]."\t".$line[4]."\t".$line[5]."\t".$line[6]."\t".$line[7].
				"\tgene_id \"" . $gene_id ."\"; transcript_id \"" . $transcript_id ."\";  \n";
			}
		}
		if ($line[2] =~ m/exon/) {
#		print $line[2] . "\n";
			if ($line[8]=~m/ID=exon:(Acey_s\d{4}.g\d*)\.t\d*\.\d*;Parent=transcript:(Acey_s\d{4}.g\d*\.t\d*)/g) {
				$gene_id = $1 ; 
				$transcript_id = $2; 
#				print $2 . "\n";
				print $line[0]."\t".$line[1]."\t".$line[2]."\t".$line[3]."\t".$line[4]."\t".$line[5]."\t".$line[6]."\t".$line[7].
				"\tgene_id \"" . $gene_id ."\"; transcript_id \"" . $transcript_id ."\";  \n";
			}			
		}		
		if ($line[2] =~ m/CDS/) { # ID=cds:Acey_s0009.g673.t1;Parent=transcript:Acey_s0009.g673.t1;
			if ($line[8]=~m/ID=cds:(Acey_s\d{4}.g\d*)\.t\d*;Parent=transcript:(Acey_s\d{4}.g\d*\.t\d*)/g) {
				$gene_id = $1 ; 
				$transcript_id = $2; 
#				print $2 . "\n";
				print $line[0]."\t".$line[1]."\t".$line[2]."\t".$line[3]."\t".$line[4]."\t".$line[5]."\t".$line[6]."\t".$line[7].
				"\tgene_id \"" . $gene_id ."\"; transcript_id \"" . $transcript_id ."\";  \n";
			}
		}
	}
}



	sub print_usage {
	  print "Summary.\nYou have to specify exactly 1 argument(s)!\n";
	  print "Usage: perl wbps_gff_2_karyoploteR_gff.pl <file.gff> \n";
	  print "  <file.gff>   = gff file to use for creating gft file\n\n";
	  print "  IMPORTANT !!! this script is highly dependent on the format of the string after \"ID\"\n";
  	  print "                For different species, it will be necessary to change the regular expression in the code\n";
	  
	}


