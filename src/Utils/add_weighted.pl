#!/usr/bin/env perl

=head1 NAME

add_weighted

=head1 SYNOPSIS

<prg_name> <weigth1> <weigth2> ... <weigthN>

Options:

=over 1

=item  B<--help>                        Brief help message

=item  B<--man>                         Full documentation

=item  B<-v, --verbose>                 Verbose

=item  B<-q, --quiet>                   Quiet

=back

=head1 DESCRIPTION

For each row in the stdin, add the columns with weights

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


## Getopt::Long::Configure("no_ignore_case");

GetOptions(
	   "verbose" => \$verbose,
	   "quiet" => \$quiet,
	   "help"=> \$help,
	   "man" => \$man
	   ) || pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

## ------------------------------------------------------------
## main part


while(my $line = <STDIN>) {
    my @nums = split /\s+/,$line;

    my $sum=0;
    for(my $i=0; $i<=$#ARGV; $i++) {
	#print "$ARGV[$i]*$nums[$i] ";
	$sum += $ARGV[$i]*$nums[$i];
    }

    print "$sum\n";
}


## ------------------------------------------------------------


