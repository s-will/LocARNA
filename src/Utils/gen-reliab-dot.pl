#!/usr/bin/env perl

=head1 NAME

<prg_name>

=head1 SYNOPSIS

<prg_name> [options]

Options:

=over 1

=item  B<--help>                        Brief help message

=item  B<--man>                         Full documentation

=item  B<-v, --verbose>                 Verbose

=item  B<-q, --quiet>                   Quiet


=item  B<-f=<file>>                     reliability file

=item  B<-a=<file>>                     alignment file

=item  B<-s>                            scale probabilities

=back

=head1 DESCRIPTION

Generate a reliability dot plot out of a list of arc match reliabilities
and the alignment in clustalw aln format.

=cut

use warnings;
use strict;

##------------------------------------------------------------
## options
use Getopt::Long;
use Pod::Usage;

my $help;
my $man;
my $quiet;
my $verbose;

my $pair_file = "result.amreliability";
my $align_file = "result.aln";
my $scale_prob;

## Getopt::Long::Configure("no_ignore_case");

GetOptions(
	   "verbose" => \$verbose,
	   "quiet" => \$quiet,
	   "help"=> \$help,
	   "man" => \$man,

           "f=s" => \$pair_file,
           "a=s" => \$align_file,
           "s" => \$scale_prob

	   ) || pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

## ------------------------------------------------------------
## main part


############################################################

# $| = 1;				# flush after every write


if (! -e "$align_file") {
    print STDERR "ERROR: File $align_file does not exist.\n";
    exit -1;
}

### call RNAalifold
my @alires = readpipe("RNAalifold -p $align_file 2>/dev/null");

my $seq=$alires[0];
die "Cannot alifold $align_file\n" if !defined($seq);

chomp $seq;

open(my $PF, "<", $pair_file) ||
	  die "Can't open pair file $pair_file: $!";
my @pairs = <$PF>;
close $PF;


### print start until sequence
###
print "%!PS-Adobe-3.0 EPSF-3.0\n";
print "%%Title: RNA Dot Plot\n";
print "%%Creator: PS_dot.c,v 1.38 2007/02/02 15:18:13 ivo Exp \$, ViennaRNA-1.7.2\n";
print "%%CreationDate: Tue Feb 17 12:00:16 2009\n";
print "%%BoundingBox: 66 211 518 662\n";
print "%%DocumentFonts: Helvetica\n";
print "%%Pages: 1\n";
print "%%EndComments\n";
print "\n";
print "%Options: \n";
print "% \n";
print "%This file contains the square roots of the base pair probabilities in the form\n";
print "% i  j  sqrt(p(i,j)) ubox\n";
print "\n";
print "%%BeginProlog\n";
print "/DPdict 100 dict def\n";
print "DPdict begin\n";
print "/logscale false def\n";
print "/lpmin 1e-05 log def\n";
print "\n";
print "/box { %size x y box - draws box centered on x,y\n";
print "   2 index 0.5 mul sub            % x -= 0.5\n";
print "   exch 2 index 0.5 mul sub exch  % y -= 0.5\n";
print "   3 -1 roll dup rectfill\n";
print "} bind def\n";
print "\n";
print "/ubox {\n";
print "   logscale {\n";
print "      log dup add lpmin div 1 exch sub dup 0 lt { pop 0 } if\n";
print "   } if\n";
print "   3 1 roll\n";
print "   exch len exch sub 1 add box\n";
print "} bind def\n";
print "\n";
print "/lbox {\n";
print "   3 1 roll\n";
print "   len exch sub 1 add box\n";
print "} bind def\n";
print "\n";
print "/drawseq {\n";
print "% print sequence along all 4 sides\n";
print "[ [0.7 -0.3 0 ]\n";
print "  [0.7 0.7 len add 0]\n";
print "  [-0.3 len sub -0.4 -90]\n";
print "  [-0.3 len sub 0.7 len add -90]\n";
print "] {\n";
print "   gsave\n";
print "    aload pop rotate translate\n";
print "    0 1 len 1 sub {\n";
print "     dup 0 moveto\n";
print "     sequence exch 1 getinterval\n";
print "     show\n";
print "    } for\n";
print "   grestore\n";
print "  } forall\n";
print "} bind def\n";
print "\n";
print "/drawgrid{\n";
print "  0.01 setlinewidth\n";
print "  len log 0.9 sub cvi 10 exch exp  % grid spacing\n";
print "  dup 1 gt {\n";
print "     dup dup 20 div dup 2 array astore exch 40 div setdash\n";
print "  } { [0.3 0.7] 0.1 setdash } ifelse\n";
print "  0 exch len {\n";
print "     dup dup\n";
print "     0 moveto\n";
print "     len lineto \n";
print "     dup\n";
print "     len exch sub 0 exch moveto\n";
print "     len exch len exch sub lineto\n";
print "     stroke\n";
print "  } for\n";
print "  [] 0 setdash\n";
print "  0.04 setlinewidth \n";
print "  currentdict /cutpoint known {\n";
print "    cutpoint 1 sub\n";
print "    dup dup -1 moveto len 1 add lineto\n";
print "    len exch sub dup\n";
print "    -1 exch moveto len 1 add exch lineto\n";
print "    stroke\n";
print "  } if\n";
print "  0.5 neg dup translate\n";
print "} bind def\n";
print "\n";
print "end\n";
print "%%EndProlog\n";
print "DPdict begin\n";
print "%delete next line to get rid of title\n";
print "270 665 moveto /Helvetica findfont 14 scalefont setfont (dot.ps) show\n";
print "\n";
print "/sequence { (\\\n";

################
### now
print "$seq\\\n";
################

################
### print until data

print ") } def\n";
print "/len { sequence length } bind def\n";
print "\n";
print "72 216 translate\n";
print "72 6 mul len 1 add div dup scale\n";
print "/Helvetica findfont 0.95 scalefont setfont\n";
print "\n";
print "drawseq\n";
print "0.5 dup translate\n";
print "% draw diagonal\n";
print "0.04 setlinewidth\n";
print "0 len moveto len 0 lineto stroke \n";
print "\n";
print "drawgrid\n";
print "%data starts here\n";


################
### print basepairs (probs in sqrt)

my $max_sqrtprob = 0;

if (defined($scale_prob)) {
    foreach my $line (@pairs) {
	if ($line =~ /([0-9]*) ([0-9]*) ([0-9.]*) *$/ ) {
	    my $sqrtprob=sqrt($3);
	    if ($sqrtprob>$max_sqrtprob) {
		$max_sqrtprob = $sqrtprob;
	    }
	}
    }
}
if ($max_sqrtprob==0) {
    $max_sqrtprob=1;
}


foreach my $line (@pairs) {
    if ($line =~ /([0-9]*) ([0-9]*) ([0-9.]*) *$/ ) {
        my $sqrtprob=sqrt($3);
	if (defined($scale_prob)) {
	    $sqrtprob /= $max_sqrtprob;
	}
        print "$1 $2 $sqrtprob ubox\n";
    }
}

################

################
### rest
print "showpage\n";
print "end\n";
print "%%EOF\n";
