#!/usr/bin/env perl

=head1 NAME

alnsel.pl

=head1 SYNOPSIS

alnsel.pl input.aln [options] [names]

=head1 DESCRIPTION

Select subset of sequences from multiple alignment file.


Options:

=over 20

=item  B<--help>

Brief help message

=item  B<--man>

Full documentation

=item  B<-v, --verbose>

Verbose

=back

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

GetOptions(
    "verbose" => \$verbose,
    "help"=> \$help,
    "man" => \$man,
    ) || pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

my $inputfilename;
if (@ARGV > 1) {
    $inputfilename = shift @ARGV;
}else {
    print STDERR "Input filename missing and/or names.\n";
    pod2usage(1);
}

## ------------------------------------------------------------
## main part

use MLocarna;

my $aln = read_clustalw_alnloh("$inputfilename");

my @alnsel=();
for my $i (0..@$aln-1) {
    my $name=$aln->[$i]{name};
    if (grep /^$name$/, @ARGV) {
	push @alnsel,@$aln[$i];
    }
}

project_alnloh(\@alnsel);

print STDOUT "CLUSTAL W\n\n";
write_clustalw_alnloh(*STDOUT, \@alnsel,75,0);

## ------------------------------------------------------------


