#!/usr/bin/env perl

=head1 NAME

aln2fa.pl

=head1 SYNOPSIS

aln2fa.pl input.aln [options]

Options:

=over 20

=item  B<--help>

Brief help message

=item  B<--man>

Full documentation

=item  B<-v, --verbose>

Verbose

=item  B<-d,--degap>

Remove gaps from sequences

=back

=head1 DESCRIPTION

Convert input.aln to fasta format

=cut

use warnings;
use strict;
use FindBin;
use lib "$FindBin::Bin/../lib/perl";

##------------------------------------------------------------
## options

use Getopt::Long;
use Pod::Usage;

my $help;
my $man;
my $verbose;

my $degap;

GetOptions(
    "verbose" => \$verbose,
    "help"=> \$help,
    "man" => \$man,

    "d|degap" => \$degap

    ) || pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

my $inputfilename;
if ($#ARGV == 0) {
    $inputfilename = $ARGV[0];
}else {
    print STDERR "Input filename missing.\n";
    pod2usage(1);
}

## ------------------------------------------------------------
## main part

use MLocarna;

my $aln = read_clustalw_alnloh("$inputfilename");

if ($degap) {
    for my $i (0..@$aln-1) {
	$aln->[$i]{seq} =~ s/[-~]//g;
    }
}

print sprint_fasta_alnloh($aln,75);

## ------------------------------------------------------------


