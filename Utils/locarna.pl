#!/usr/bin/perl

############################################################
#
#  front-end for locarna
#
#  USAGE: locarna.pl [<seq1> <seq2>] [<options>]
#
#  use this for aligning two RNA given by sequences alone
#
#  calls RNAfold -p for predicting matrices of pair
#  probabilities and hands this to locarna
#  
############################################################

use strict;
use warnings;

### --------------------


my $seq1;
my $seq2;

## get directory of executable "cmd_dir"
my $basename=`basename $0`;
chomp $basename;
my $cmd_dir;
($cmd_dir = $0) =~ s/$basename$//; 
$cmd_dir =~ s/\/$//;


sub usage {
    print STDERR "USAGE: locarna.pl [<seq1> <seq2>] [<options>]\n";
}

## count non-option arguments and search for help
my $non_opt_argc=0;
foreach my $arg (@ARGV) {
    if (substr($arg,0,1) ne "-") {
	$non_opt_argc++;
    } elsif ($arg eq "--help") {
	usage;
	exit -1;
    }
}


if ($non_opt_argc==2) {
    $seq1=$ARGV[0];
    $seq2=$ARGV[1];

    shift @ARGV; shift @ARGV;    
} elsif ($non_opt_argc==0) {
    print "\n";
    print "Input two times name and sequence (upper or lower case)\n";
    print "....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8\n";
    
    (<STDIN> =~ /^>\s*([^\s]*)/) || die "Expected name 1.\n";
    my $line;
    $seq1="";
    while(($line=<STDIN>) =~ /^[^>]/) {chomp($line);$seq1.=$line;}
    ($line =~ /^>\s*([^\s]*)/) || die "Expected name 2.\n";
    $seq2="";
    while($line=<STDIN>) {chomp($line);$seq2.=$line;}
} else {
    usage;
    exit(-1);
}


$seq1=uc $seq1;
$seq2=uc $seq2;

my $seq1name="seqA";
my $seq2name="seqB";

my $tmpbase="/tmp";

my $tmpdir=readpipe("mktemp -d $tmpbase/locarna.XXXXXXXXXX");
chomp $tmpdir;

my $workdir;
chomp($workdir = `pwd`);

chdir("$tmpdir") || die("Cannot change dir to $tmpdir");


system("(echo \\\>$seq1name ; echo $seq1) | RNAfold -p");
system("(echo \\\>$seq2name ; echo $seq2) | RNAfold -p");

chdir("$workdir");

my $dpps="_dp.ps";

my $align_cmd="$cmd_dir/locarna $tmpdir/$seq1name$dpps $tmpdir/$seq2name$dpps";

if ($#ARGV>-1) {
    $align_cmd.=" @ARGV";
}

system($align_cmd);

system("rm -rf $tmpdir");

