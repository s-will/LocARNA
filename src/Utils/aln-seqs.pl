#!/usr/bin/env perl

=head1 NAME

aln-seqs

=head1 SYNOPSIS

aln-seqs files

=head1 DESCRIPTION

List aln files together with the contained sequences, (optionally)
sort output by number of sequences


=head1 OPTIONS

=over 1

=item  B<--help>                        Brief help message

=item  B<--man>                         Full documentation

=item  B<--sort>                        Sort by number of sequences

=item  B<-0,--null>                     Separate file name from sequence names by \0

=item  B<-no-hash>                      Remove names with prefix '#'

=back


=cut

use warnings;
use strict;

##------------------------------------------------------------
## options
use Getopt::Long;
use Pod::Usage;

my $help;
my $man;

my $opt_sort;
my $opt_null;
my $opt_nohash;

## Getopt::Long::Configure("no_ignore_case");

GetOptions(
    "help"=> \$help,
    "man" => \$man,
    "sort" => \$opt_sort,
    "0|null" => \$opt_null,
    "no-hash" => \$opt_nohash
    ) || pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

## ------------------------------------------------------------
## main part

my @rows=();

for my $file (@ARGV) {
    my @row=();
    push @row, "$file";
    open(my $fh, "<", $file) || die "Cannot read file from $file: $!";
    while(my $line=<$fh>) {
	if ($line !~ /^CLUSTAL/) {
	    if ($line =~ /^(\S+)\s+/) {
		my $name=$1;
		if (! member($name,@row)) {
		    push @row, "$name" unless ($opt_nohash && $name=~/^#/);
		}
	    }
	}
    }
    close $fh;
    push @rows, [ @row ];
}

if ($opt_sort) {
    @rows = sort {@$a<=>@$b} @rows;
}

for my $row (@rows) {
    my $filename=shift @$row;
    print $filename;
    if ($opt_null) {
	print "\0";
    } else {
	print " ";
    }
    print join(" ",@$row)."\n";
}

sub member {
    my $x=shift @_;
    return (grep {$x eq $_} @_) > 0;
}
