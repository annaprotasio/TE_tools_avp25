#!/usr/bin/perl -w
use strict;

my $arg_count = @ARGV;
if ($arg_count != 1) {
  &print_usage;
  die "Please specify 1 arguments!";
}
my $input_fasta=$ARGV[0];
open(IN,"<$input_fasta") || die ("Error opening $input_fasta $!");

my $line = <IN>; 
print $line;

while ($line = <IN>)
{
chomp $line;
if ($line=~m/^>/) { print "\n",$line,"\n"; }
else { print $line; }
}

print "\n";

sub print_usage {
  
  print "Usage: fasta_multi_to_one_line.pl <fasta> \n";
  print "  <file1> = file containes fasta file \n";
  print "  Trasnforms this: \n";
  print ">D915_04184\n";
  print "MDALRATLDSTDGQVTTVSRTSSPPPLDILTPPVVAEYWCHTQVRVTKMKYIWTISNFSF\n";
  print "CREEMGEVVKSSFFSCGPNDKLKWCLRINPKGLDEESREYLSLYLLLVNCGTKSEARAKF\n";
  print "KFSLLNAKREETKAMGNFCLFFQKILNSYVLYCIFFFSCIN*\n";
  print ">D915_11077\n";
  print "MLRLGTIHCTQEDYNQCLELVPPSLRFYVEAIYERRQAIVPDFEALHQQRRLFVEARLAA\n";
  print "AAANQAAGQTDENGDPIGPGTVGGRSMDHSGLGRSGLRSLRHSPMDTGAAYTSIFESSTN\n";
  print "ESTLRQLYLSRTADDGDVSRLGSGTSSMGSGGFVPGDTQHYARAAASRGYVVPGGPGGLL\n";
  print "RHRQLPVSVAQSSSGNLGMGPRSFVTSERYIHDRDTKPALICSVAGGSREDMSSTNALVS\n\n";
  print "Into this: \n";
  print ">D915_04184\n";
  print "MDALRATLDSTDGQVTTVSRTSSPPPLDILTPPVVAEYWCHTQVRVTKMKYIWTISNFSFCREEMGEVVKSSFFSCGPNDKLKWCLRINPKGLDEESREYLSLYLLLVNCGTKSEARAKFKFSLLNAKREETKAMGNFCLFFQKILNSYVLYCIFFFSCIN*\n";
  print ">D915_11077\n";
  print "MLRLGTIHCTQEDYNQCLELVPPSLRFYVEAIYERRQAIVPDFEALHQQRRLFVEARLAAAAANQAAGQTDENGDPIGPGTVGGRSMDHSGLGRSGLRSLRHSPMDTGAAYTSIFESSTNESTLRQLYLSRTADDGDVSRLGSGTSSMGSGGFVPGDTQHYARAAASRGYVVPGGPGGLLRHRQLPVSVAQSSSGNLGMGPRSFVTSERYIHDRDTKPALICSVAGGSREDMSSTNALVS\n\n";
  
}

exit
