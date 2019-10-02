#!/usr/bin/env perl

use FindBin;
use lib "$FindBin::Bin/../../lib/perl";

use MLocarna::Tree;
use Data::Dumper;

my $newick_ex1="((A,((B1,B2:12.4)I1,(C,D),E:2)):1,F);";

my $tree = new MLocarna::Tree($newick_ex1);

# print Dumper( $tree )."\n";

print $tree->to_newick().";\n";

my $label_num=1;
sub callback {
    my $tree = shift;
    my @child_results=@_;

    my $label = $tree->label();
    if ($tree->is_leaf()) {
        return "$label";
    } else {
        if (!defined $label) {
            $label = $label_num++;
        } else {
            $label="_$label";
        }
        return {label=>$label, children=>\@child_results};
    }
}

my $res = $tree->traverse( \&callback );

#print Dumper($res)."\n";

my $partitions = $tree->traverse_partitions(sub {
    my $partition = shift;
    print "@$partition\n";
});

my $tree2 = new MLocarna::Tree("(A,B,(C,D,E,F)Label:10)Blubb:5");
$tree2 = $tree2->binarize();
print Dumper($tree2)."\n";
