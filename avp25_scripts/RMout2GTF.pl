#!/usr/bin/perl -w
use strict;
use warnings;

my $arg_count = @ARGV;
if ($arg_count != 2) {
  &print_usage;
  die "Please specify 2 arguments!";
}


my $RM_file =  $ARGV[0]; 		#the *.cat.out file
my $min_len =  $ARGV[1];		#minimum length of repeat match to be considered

open OUT0, ">", "RepeatMasker.$RM_file.min$min_len.gtf" or die "ooooooooooooooops\n" ; 

open (IN, "$RM_file") or die "oops!\n" ;
while (<IN>) {
    next if /^\n/ ;
    chomp ;
    my @line = split /\s+/, $_ ;
    #print $line[1];
    	next if $line[1] =~ /\D/;
    	next if $line[11] =~ /Simple_repeat/ ;
	    next if $line[11] =~ /Low_complexity/ ;
   		next if $line[10] =~ /\(.+\)n/ ;
    	next if $line[11] =~ /Satellite/ ;
	    #next if $line[10] =~ /Unknown/ ;
   		next if $line[11] =~ /rRNA/ ;
   		next if $line[11] =~ /snRNA/ ;
		
		if ($line[7] - $line[6] > $min_len) {
		
		my $strand = "";
    	if ( $line[9] eq 'C') {
		$strand = "-";
    	}
    	else {$strand = "+";
    	}
    
    	my @te = "";
		my $order = "";
    	my $superfam = "";

    	if ( $line[11] =~ m{/}) {
		@te = split "/" , $line[11];
		$order = $te[0] ; $superfam = $te[1]
    	} else { #($line[11] =~ /Unknown/) {
    	$order = $line[11] ; $superfam = "Unknown"
    	}
    print OUT0 "$line[5]\tRepeatMasker\trepeat_region\t$line[6]\t$line[7]\t$line[2]\t$strand\t.\trepeat_id \"$line[10]\"; Order \"$order\"; Superfamily \"$superfam\";\n";	
    #print "$line[5]\tRepeatMasker\trepeat_region\t$line[6]\t$line[7]\t$line[2]\t$strand\t.\trepeat_id \"$line[10]\"; Superfamily \"$superfam\";Order \"$order\";\n";	
    }
    };


close(IN);

sub print_usage {
  
  print "Usage: $0 <RM.out> <min_length> \n\n";
  print "Trasnforms a RepeatMasker <RM.out> file into GTF. Ignores any matches that are less than <min_len> bp long\n";
  print "The'repeat classfamily' column is split into 'Order' and 'Superfamily'\n\n";
  
}

exit
