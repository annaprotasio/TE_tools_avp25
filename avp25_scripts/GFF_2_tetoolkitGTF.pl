#!/usr/bin/perl -w
use strict;
use warnings;

my $arg_count = @ARGV;
if ($arg_count != 2) {
  &print_usage;
  die "Please specify 2 arguments!";
}

my $infile  =  shift;		#name of file where output is to be save, comes from args
my $outfile = shift;

open( OUT, ">$outfile");
print OUT "\#\#\#gtf-annotation-for-tetoolkit\n";


# open file
open( IN, "$infile" ) or die "Cannot open IN file\n";
my $str = "";
my $transc_id = "";
my $count = 1;

while (<IN>) {
	next if /^\#/ ;
    my @line = split /\s+/, $_ ;
    my $str = substr $line[9], 1, -2;
    my $transc_id = join '', $str, '_dup' , $count; 
print OUT "$line[0]\tte_annot\texon\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\tgene_id $line[9] transcript_id \"$transc_id; family_id $line[11] class_id $line[13]\n";
    $count++;
       
    }

close OUT or die $!;


sub print_usage {
  
  print "Usage: $0 <IN_file_GFF> <OUT_file_GFF>\n\n";
  print "Parser script, turing a GFF into a GTF that can be used with tetoolkit. Turns this:\n\n";
  print "Chr	rpmk	one	1000	2000	7.9	-	.	Family \"HARLEQUIN_LTR_solo\"; Superfamily \"Gypsy\"; Order \"SINE\"; Perc_ref \"0.522\"\n";
  print "Chr	rpmk	one	3000	4000	18.3	-	.	Family \"T2\"; Superfamily \"Unknown\"; Order \"SINE\"; Perc_ref \"0.581\"\n\n";
  print "Into: \n";
  print "Chr	rpmk	one	1000	2000	7.9	+	.	gene_id \"HARLEQUIN_LTR_solo\"; transcript_id \"HARLEQUIN_LTR_solo_dup1; family_id \"Gypsy\"; class_id \"LTR\";\n";
  print "Chr	rpmk	one	3000	4000	18.3	-	.	gene_id \"T2\"; transcript_id \"T2_dup1; family_id \"Unknown\"; class_id \"SINE\";\n\n";
;
}


exit