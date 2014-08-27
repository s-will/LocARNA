#!/usr/bin/env perl

## ATTENTION: absolute path
my $nupackPath = "/usr/local/user/nupack3.0";


#####################################################################
##
##  CALL : ./NuPack_dotplot.pl <SEQID> <SEQUENCE>
##
##  will generate the following files:
##   + SEQID.ppairs : NuPack output
##   + SEQID.pp     : the according conversion to .pp format
##
####################################################################

use warnings;
use strict;

# read sequence parameter
my $seqID = $ARGV[0];
my $seq = $ARGV[1];

# set folding temperature
my $temperature = 37;


my $nupack = "$nupackPath/bin/pairs";

my $verbose = 0;

# run NUPACK

system("export NUPACKHOME=$nupackPath; printf \"$seqID\\n$seq\" | $nupack -pseudo -dangles some -T $temperature -material rna1999 > /dev/null;");

# read/convert NUPACK output

my %dotplot = ();

my $length = length($seq);

open(PAIRSFILE,"$seqID.ppairs");
foreach my $zeile (<PAIRSFILE>){
#		if($verbose){print $zeile;}
    if($zeile =~/\%.+/){
	if($verbose){print "Skipped: $zeile";}
    }elsif($zeile =~ /^(\d+)$/){
	if($verbose){print "Alte L: $length neue L: $zeile";}
	$length = $1;
    }elsif($zeile =~ /(\d+)\s(\d+)\s(\S+)\s(\S+)\s(\S+)/){
	if($verbose){print "Infos extrahiert: $zeile";}
	if($2 != -1 && $2 <= $length){
	    $dotplot{"$1,$2"} = [$1,$2,$3];
	}
    }else{
	if($verbose){print "Skipped $zeile";}
    }
}
close PAIRSFILE;
#	unlink("$seqID.ppairs");

# write dotplot output

open(DPOUT,">$seqID.pp");

print DPOUT "#PP 2.0\n\n";

print DPOUT "$seqID $seq\n\n#END\n\n";

print DPOUT "#SECTION BASEPAIRS\n\n";

foreach my $value (values %dotplot){
    my $prob = sqrt( ${$value}[2]);
    print DPOUT "${$value}[0] ${$value}[1] $prob\n";
}

print DPOUT "\n#END";

close DPOUT;

exit 0;
#### END OF SCRIPT


