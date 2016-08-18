#!/usr/bin/env perl

use strict;
use warnings;

my @B=('A','C','U','G');

for (my $i=0; $i<$ARGV[0]; $i++) {
    print "$B[rand 4]"; 
}
print "\n";
