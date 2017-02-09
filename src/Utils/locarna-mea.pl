#!/usr/bin/env perl

=head1 NAME

locarna_mea

=head1 SYNOPSIS

locarna_mea infile.pp

=head1 DESCRIPTION

Reads a base pair probabilities from a pp file and computes the
secondary structure with maximum sum of base pair probabilities.

(Ignores potential annotation including structure constraints.)

=head1 OPTIONS

=over 1

=item  B<--help>                        Brief help message

=item  B<--man>                         Full documentation

=back


=cut

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin/../lib/perl";

use MLocarna;
use MLocarna::Aux;
use MLocarna::MatchProbs;


use Cwd 'abs_path';

##------------------------------------------------------------
## options
use Getopt::Long;
use Pod::Usage;

my $help;
my $man;
my $threshold=0;

## Getopt::Long::Configure("no_ignore_case");

GetOptions(
    "help"=> \$help,
    "man" => \$man,
    "threshold=f" => \$threshold
    ) || pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

my $inputfile = $ARGV[0];
if (!defined($inputfile)) {
    print STDERR "Please specify input file.\n";
    pod2usage(-exitstatus => -1, -verbose => 1);
}


## ------------------------------------------------------------
## main part

my %pairprobs = read_pp_file_pairprobs("$inputfile");
my %aln = read_pp_file_aln_wo_anno("$inputfile");

# compute and print a maximum reliability structure for the consensus reliabilities
my @empty=();
my ($score,$rel_str)=max_weight_structure(aln_length(\%aln),\@empty,\%pairprobs,1,$threshold);
print "$rel_str\n";

