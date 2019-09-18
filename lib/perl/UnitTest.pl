#!/usr/bin/perl -I. -w

## some unit tests (very incomplete :))

use MLocarna::SparseMatrix;


my $tmpfile="tmp.$$";


########################################
print "Testing MLocarna::SparseMatrix::read_2D\n";

open(OUT,">$tmpfile") || die "Cannot write to $tmpfile\n";
print OUT "1 2 0.5\n10 15 1.75e+2\n";
close OUT;
my %sm = MLocarna::SparseMatrix::read_2D($tmpfile);
if ($sm{1}{2} == 0.5 && $sm{10}{15} == 1.75e+2) {
    print "OK\n";
} else {
    print "FAIL $sm{10}{15}\n";
}
