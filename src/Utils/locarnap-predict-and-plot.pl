#!/usr/bin/env perl

=head1 NAME

locarnap-predict-and-plot.pl

=head1 SYNOPSIS

locarnap-predict-and-plot.pl [options] [<annotation-file>]

=head1 DESCRIPTION

Performs boundary and reliability prediction and draws all reliability
plots according to annotation (given as file or read from standard
input). The script is usually used after generating alignments with
locarnap-realign-all.pl as third step in a pipeline for refining RNAz
hits with LocARNA-P.

=head1 OPTIONS

=over 4

=item  B<--help>

Brief help message

=item  B<--man>

Full documentation

=item  B<--test>

Test

=item  B<--output-dir>=d

Output directory (def=Relplots)

=item  B<--dont-plot>

Skip plotting, only output

=item  B<--show-sw>

Show the structure weight in the plot

=item  B<--revcompl>

Draw for reverse complement (3'-5')

=item  B<--write-subseq>

Write the subsequence of fit

=item  B<--output-format>=f

Output format (f = pdf or png, def=pdf)

=back

By default plots are written to directory Relplots. The predictions
are written to standard out as a table. A line of the table contains
of the locus name, start,end, and orientation of the RNAz prediction,
the LocARNA prediction and the first annotation, the on and off value
of the fit, and the background and hit reliability.

=cut

use warnings;
use strict;
use FindBin;
my $bindir = "$FindBin::Bin";

## ------------------------------------------------------------
## parameter

my $seqname="DroMel";
my $str_weight=3;
my $delta=0.5;



##------------------------------------------------------------
## options
use Getopt::Long;
use Pod::Usage;

my $help;
my $man;
my $quiet;
my $test;
my $verbose;

my $output_dir="Relplots";

my $dont_plot;
my $show_sw;
my $revcompl;

my $write_subseq;

my $output_format="pdf";

## Getopt::Long::Configure("no_ignore_case");

GetOptions(
    "verbose" => \$verbose,
    "quiet" => \$quiet,
    "test" => \$test,
    "help"=> \$help,
    "man" => \$man,
    "output-dir=s" => \$output_dir,
    "dont-plot" => \$dont_plot,
    "show-sw" => \$show_sw,
    "revcompl" => \$revcompl,
    "write-subseq" => \$write_subseq,
    "output-format=s" => \$output_format
    ) || pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;


## ------------------------------------------------------------
## main part


if ( ! -d $output_dir ) {
    print STDERR "Output directory $output_dir does not exist. Make directory ...";
    mkdir $output_dir;
    if ( ! -d $output_dir ) {
	print STDERR "FAILED.\n";
	exit -1;
    }
    print STDERR "DONE.\n";
}



while(<>) {

    my @l=split /\s+/;

    ## get left offset!

    my $filename="$l[1]:$l[0]";

#    if ($#l>=5 && $l[5] =~ /^MI/) { ## special handling for micro RNA precursor
#	$filename="$l[0]-$l[10]";
#    }

    if ($revcompl) {
	$filename.="-rc";
    }

    my $offset=0;

    if ($#l>=2) {
	my $left_context;
	my $sequences="Realign-Sequences/$filename.mfa"; ## name of sequence file in Realign-Sequences

	open(my $SEQ, "<", "$sequences") || die "Cannot read from $sequences: $!";
	while(<$SEQ>) {if ($_=~/DroMel_CAF1.*left_context=(\d+)/) {$left_context=$1;}}
	close $SEQ;
	if (!defined($left_context)) {
	    print STDERR "Could not find left context in $sequences. Assume left context of 100.";
	    $left_context=100;
	}

	$offset=$l[2]-$left_context;
    }

    my $signals = "";
    for (my $i=2; $i<$#l; $i+=5) {
	$signals .= "$l[$i] $l[$i+1] ";
	if ($l[$i+2] eq "+") { $signals.="+1"; }
	elsif ($l[$i+2] eq "-") { $signals.="-1"; }
	else { $signals.="0"; }
	$signals.="; "
    }
    $signals=~s/; $//;

    my $signal_names="";
    if ($#l>=5) { $signal_names="RNAz "; }
    for (my $i=6; $i<$#l; $i+=5) {
	$signal_names .= "$l[$i] ";
    }

    my $locus_name="$l[1]_$l[0]";
    my $locus_title_name="$l[1]:$l[0]";

    my $outfile="$output_dir/$locus_name".($revcompl?"-rc":"");

    my $command="$bindir/reliability-profile.pl ".
	"--seqname $seqname ".
	"--structure-weight=$str_weight ".
	($show_sw?"--show-sw ":"").
	($dont_plot?"--dont-plot ":"").
	($revcompl?"--revcompl ":"").
	($write_subseq?"--write-subseq ":"").
	"--output-format=$output_format ".
	"--fit-penalty=$delta ".
	"--fit-once-on ".
	"--title=\"$locus_title_name\" ".
	"--signal-names \"$signal_names\" ".
	"Alignment-Results/$filename.dir ".
	"--signals \"$signals\" ".
	"--offset=$offset ".
	"--out=\"$outfile\"";

    my @res;
    if ($test) {
	print STDERR "$command\n";
	exit 0;
    } else {
	@res=readpipe($command);
    }

    my $from;
    my $to;
    my $orientation=($revcompl)?"-":"+";

    my $onoff="";

    my $hit_score="";
    my $bg_score="";

    foreach my $l (@res) {
	if ($l=~/FIT (\d+)\s+(\d+)/) {
	    $from=$1;
	    $to  =$2;
	} elsif ($l=~/ONOFF (\S+\s+\S+)/) {
	    $onoff=$1;
	} elsif ($l=~/SEQ (.*)/) {
	    print "$1\n";
	} elsif ($l=~/SCORE\s+(\S+)\s+(\S+)/) {
	    $hit_score=$1;
	    $bg_score=$2;
	}
    }

    print "$locus_name "; ## the locus name

    if ($#l>=4) {
	print "$l[2] $l[3] $l[4] "; ## rnaz hit
    } else {
	print "? ? ? "; ## no rnaz hit
    }
    print "$from $to $orientation "; ## fit

    if ($#l>=9) {
	## compare to the first annotation
	print "$l[7] $l[8] $l[9] ";
    }

    print "$onoff ";

    print "$bg_score $hit_score";

    print "\n";

    #print "@res";
}



## ------------------------------------------------------------
