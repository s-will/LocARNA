package MLocarna::Trees;

############################################################
##
## package for constructing trees by PGMA
## and handling tress in NEWICK format
##
############################################################


use 5.008003;
use strict;
use warnings;

use Text::ParseWords qw(parse_line) ; ## for tokenizing
use Scalar::Util qw(looks_like_number);

require Exporter;

# set the version for version checking
our $VERSION     = 1.00;

our @ISA         = qw(Exporter);
our @EXPORT      = qw(
);

our %EXPORT_TAGS = ();

# your exported package globals go here,
# as well as any optionally exported functions
our @EXPORT_OK   = qw(
check_tree_labels
$node_sym
newick_tree_to_postorder
project_tree
quote_newick_label
tree_partitions
unquote_newick_label
upgma_tree
);

our $node_sym="\$\$nodesym";

########################################
## quote_newick_label
##
## @param label unquoted label
## @returns quoted label
##
## quotes special characters, such that the text
## is valid as label in newwick format
########################################
sub quote_newick_label {
    my $label = shift;

    my $add_ticks=0;
    if ( $label =~ /[\':]/ ) {
	$add_ticks = 1;
    }

    $label =~ s/\'/\'\'/g;

    if ($add_ticks) {
	$label = "\'$label\'"
    }

    return $label;
}

########################################
## upgma_tree_dist( $names ref to list, $dist_matrix ref to matrix of distances)
##
## Compute an upgma-tree from distances
##
## @param $names ref to a list of names
## @param $dist_matrix ref to a matrix (array of arrays) of similarity dists
## @param $add_branch_lengths boolean 0/1-flag, control whether to add
##        branch lengths in the tree string [default on]
##
## The output uses the names in list @$names. Consequently, names
## should be unique and indices in $@names and matrix @$dist_matrix
## should correspond
##
## @return tree as string in NEWICK tree format
########################################
sub upgma_tree_dist {
    my ($names, $dist_matrix, $add_branch_lengths) = @_;
    if (!defined($add_branch_lengths)) {
	$add_branch_lengths = 1;
    }

    ## compute tree by upgma applied to dist matrix

    my @clusters; # a list of the clusters
    my @trees; # a list of sub-trees
    my @cluster_sizes;
    my @heights;

    for (my $i=0; $i<@$names; $i++) {
	$clusters[$i] = $i;
	$trees[$i]    = quote_newick_label($names->[$i]);
	$cluster_sizes[$i] = 1;
	$heights[$i]  = 0;
    }

    my $INFINITY = 1e10;

    while ($#clusters>0) {
        ## find the nearest two clusters (with minimal similarity)
	my $min_i;
	my $min_j;
	my $min_dist=$INFINITY;
	for (my $i=0; $i<=$#clusters; $i++) {
	    for (my $j=$i+1; $j<=$#clusters; $j++) {
		my $dist=$dist_matrix->[$clusters[$i]][$clusters[$j]];
		if ($dist < $min_dist) {
		    $min_i=$i;
		    $min_j=$j;
		    $min_dist=$dist;
		}
	    }
	}

	## recompute similarities
	my $cluster_i = $clusters[$min_i];
	my $cluster_j = $clusters[$min_j];

	## update the list clusters
	$clusters[$min_j] = $clusters[$#clusters];
	$clusters[$min_i] = $clusters[0];
	$clusters[0] = $cluster_i;

	for (my $i=1; $i<$#clusters; $i++) {
	    $dist_matrix->[$clusters[0]][$clusters[$i]] =
		($cluster_sizes[$cluster_i] * $dist_matrix->[$cluster_i][$clusters[$i]]
		 + $cluster_sizes[$cluster_j] * $dist_matrix->[$cluster_j][$clusters[$i]])
		/ ($cluster_sizes[$cluster_i]+$cluster_sizes[$cluster_j]);
	    $dist_matrix->[$clusters[$i]][$clusters[0]] =
		$dist_matrix->[$clusters[0]][$clusters[$i]];
	}

	## height of the new node and branch lengths
	my $height = $min_dist / 2.0;

	my $new_tree;
	if ($add_branch_lengths) {
	    my $ilen = sprintf("%.3f",$height - $heights[$cluster_i]);
	    my $jlen = sprintf("%.3f",$height - $heights[$cluster_j]);
	    $new_tree = "(".$trees[$cluster_i].":".$ilen.",".$trees[$cluster_j].":".$jlen.")";
	} else {
	    $new_tree = "(".$trees[$cluster_i].",".$trees[$cluster_j].")";
	}

	$trees[$clusters[0]] = $new_tree;

	$cluster_sizes[$clusters[0]] = $cluster_sizes[$cluster_i]+$cluster_sizes[$cluster_j];
	$heights[$clusters[0]] = $height;

	$#clusters--;
    }

    return $trees[$clusters[0]].";";
}

sub max {
    my ($x,$y) = @_;
    return $x>=$y ? $x : $y;
}

########################################
## scores_to_dists( $score_matrix )
##
## @param score_matrix ref to a matrix of scores
##
## @return distance matrix
########################################
sub scores_to_dists {
    my $score_matrix = shift;

    my $n = scalar( @$score_matrix );
    my $max = $score_matrix->[0][0];
    for(my $i=0;$i<$n;$i++) {
	for (my $j=$i+1; $j<$n; $j++) {
	    $max = max( $score_matrix->[$i][$j], $max );
	}
    }

    my @dist_matrix;
    for(my $i=0;$i<$n;$i++) {
	for (my $j=0; $j<$n; $j++) {
	    $dist_matrix[$i][$j] = $max - $score_matrix->[$i][$j];
	}
    }

    return \@dist_matrix;
}

########################################
## upgma_tree( $names ref to list, $score_matrix ref to matrix of scores)
##
## Compute an upgma-tree from similarities (like locarna scores)
##
## Transforms the scores to distances before performing the upgm-algorithm
##
## @param $names ref to a list of names
## @param $score_matrix ref to a matrix (array of arrays) of similarity scores
##
## The output uses the names in list @$names. Consequently, names
## should be unique and indices in $@names and matrix @$score_matrix
## should correspond
##u
## @return tree as string in NEWICK tree format
##
########################################
sub upgma_tree {
    my ($names, $score_matrix) = @_;

    # convert scores to distances
    my $dist_matrix = scores_to_dists($score_matrix);

    return upgma_tree_dist($names, $dist_matrix);
}

########################################ä
## tree_partitions($tree_postorder ref of list)
##
## generate partitions out of the postorder tree
##
## $tree_postorder ref of list representing tree in postorder (as
## generated by newick_tree_to_postorder)
##
## returns ref of list of refs of partitions of leaves due to the tree; a
## partition is represented as a list/subset of leaves
##
########################################
sub tree_partitions {
    my $tree_postorder = shift;

    my @result;

    my @stack;

    for my $item (@$tree_postorder) {
	if ($item eq $node_sym) {
	    my @op1 = @{ $stack[-2] };
	    my @op2 = @{ $stack[-1] };

	    $#stack-=2;

	    my @op12 = (@op1, @op2);

	    push @stack, \@op12;
	} else {
	    push @stack, [ $item ];
	}
	push @result, $stack[-1];
    }
    $#result-=2; # the last is empty, the one before symmetric
    return \@result;
}

########################################
##
## unquote a newick label
##
## @params string
## @return unquoted string
##
sub unquote_newick_label {
    my ($s)=@_;
    if ($s =~ /^\'(.*)\'$/) {
	$s = "$1";
	$s =~ s/\'\'/\'/g;
    }
    return $s;
}

########################################ä
## newick_tree_to_postorder($tree string)
##
## Translates a newick tree format string into a list of nodes in postorder
##
## $tree string in newick tree format
##
## Ignores branch length and names of inner nodes
##
## returns ref to list of nodes/leaves in postorder, use $node_sym for inner nodes
##
## supports quotation via module Text::ParseWords
##
## dies if tree is not parsable
##
########################################
sub newick_tree_to_postorder {

    my ($tree) = @_;

    $tree =~ s/;$//; # allow that the tree string is terminated by ';'

    # tokenize line and preserve delimiters, also preserve quotes and backslashs
    my @tokens = parse_line('[,()\s:]', "delimiters", $tree);

    my @list=();

    my $brcount=0;
    for (my $i=0; defined($tokens[$i]); $i++) {
	my $tok=$tokens[$i];

	if ($tok eq "(") {
	    $brcount++;
	} elsif ($tok eq ")") {
	    $brcount--;
	    if ($brcount<0) {
		die "Parse error in tree.";
	    }

	    $tok = $tokens[$i+1];
	    if ( defined($tok) && $tok =~ /^[a-zA-Z0-9]/ ) {
		$i++;
		#ignore names of inner nodes
	    }
	    push @list, $node_sym;
	} elsif ($tok eq ",") {
	    ## ignore, although we could do syntax checking
	} elsif ($tok =~ /^\s*$/)  {
	    ## ignore whitespace
	} elsif ($tok eq ":") {
	    $i++;
	    if ( ! looks_like_number($tokens[$i]) ) {
		die "Distance expected in tree.";
	    }
	} else {
	    push @list, unquote_newick_label($tok);
	}
    }

    return \@list;
}

########################################
## old version of newick tree parsing without quoting support
sub newick_tree_to_postorder_old {
    my ($tree) = @_;

    my @list;

    #$tree =~ s/:[\d.e-]+//g;

    $tree =~ s/;$//; # allow that the tree string is terminated by ';'

    my $brcount=0;
    for (my $i=0; $i< length $tree; $i++) {
	my $c=substr $tree,$i,1;

	if ($c eq "(") {
	    $brcount++;
	} elsif ($c eq ")") {
	    $brcount--;
	    if ($brcount<0) {
		die "Parse error in tree.";
	    }

	    push @list, $node_sym;
	} elsif ($c eq ",") {
	    ## ignore, although we could do syntax checking
	} else {
	    my $rest=substr $tree,$i;
	    $rest =~ /^([^(),]+)/;
	    my $token = $1;
	    $i+=(length $token)-1;
	    push @list, $token;
	}
    }

    return \@list;
}


########################################
## check whether a tree contains a set of labels
##
## @param $tree ref to tree in postorder
## @param labels ref to list of labels
sub check_tree_labels {
    my ($tree,$labels) = @_;

    my %existing_labels;

    foreach my $i (@$tree) {
	$existing_labels{$i}=1;
    }

    foreach my $i (@$labels) {
	if (!exists $existing_labels{$i}) {
	    return 0;
	}
    }

    return 1;
}

########################################
## project a tree to a set of labels
##
## @param $tree ref to tree in NEWICK format
## @param labels ref to list of labels
##
## @returns tree in postorder that is a sub-tree of @$tree
## and contains exactly the labels @$labels

sub project_tree {
    my ($tree,$labels) = @_;

    my @stack;

    foreach my $item (@$tree) {
	my @list=();

	if ($item eq $node_sym) {
	    my @x = @{ $stack[-2] };
	    my @y = @{ $stack[-1] };
	    $#stack-=2;

	    if (@x==0 && @y==0) {
		@list = ();
	    } elsif (@x==0) {
		@list = @y;
	    } elsif (@y==0) {
		@list = @x;
	    } else {
		@list=(@x,@y,$node_sym);
	    }
	} else {
	    @list=();
	    if (grep (/^$item$/, @$labels)!=0) {
		@list=($item);
	    }
	}
	push @stack, [ @list ];
    }

    return $stack[0];
}
