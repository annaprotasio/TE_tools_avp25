#!/usr/local/bin/perl -w

### this script was written by Anna V.Protasio / date July 2020
### it takes files in EMBL format, collects all features (FT) and makes a GTF 
### DEPENDENCIES: must be run where all the EMBL files are and these, should have the extension *.embl


use strict;
use warnings;
use File::Basename;

my $arg_count = @ARGV;
if ($arg_count != 1) {
  &print_usage;
  die "Please specify 1 argument(s) - the directory path where the *.embl files are located!";
}

my $dir = shift;

my $outfile  =  "all_LTR_annotation.gff";		#name of file where output is to be save, comes from args
my $outfile2  =  "all_LTR_annotation.bed";

open( OUT, ">$outfile");
print OUT "\#\#\#gff-annotation.\n";

open( OUT2, ">$outfile2");

for my $file (glob ($dir . '*.embl')) {
    open my $FH, '<', $file or die "$file: $!";
    open( OUT, ">>$outfile");
    my $chr = "";
    my $family = "";
    my $type = "";
    my $subpart = "";
    my @coords = "";
    if ( basename($file) =~ /(.*)#(.*)\.embl/) {$chr = $1."#".$2; $family = $1};
#    print "$chr\n";
    while (<$FH>) {
     my @line = split /\s+/, $_ ;
     if ($line[1] =~ /^LTR|^CDS_motif/) {
     	@coords = split /\.\./, $line[2];
     	if ($line[1] =~ /LTR/) {
     	$type = "LTR"; if ($coords[0] == 1) {$subpart = "5LTR"} else {$subpart = "3LTR"} 
     	}
     	elsif (($line[1] =~ /CDS_motif/)) {$type = "INT"; $subpart = "INT"} 
     print OUT "$chr\tmanual\t$type\t$coords[0]\t$coords[1]\t.\t+\t.\tFamily \"$family\"; Type \"$type\" Subpart \"$subpart\"\n";
     print OUT2 "$chr\t$coords[0]\t$coords[1]\t$family\_$subpart\t0\t+\n";
} else {} 
}
}
close OUT or die $!;
close OUT2 or die $!;

print "Successfully completed\n";

sub print_usage {
  print "Summary.\n\n";
  print "Usage: perl convert_embls_to_gff_and_bed.pl <path_to_dir_with_embl_files> \n\n";
  print "This script was written by Anna V.Protasio / date July 2020. ";
  print "It takes file(s) in EMBL format, collects all features with LTR and CDS_motif and makes a GTF and BED files\n\n";
  print "DEPENDENCIES: EMBL files should have the extension *.embl\n\n";
}



