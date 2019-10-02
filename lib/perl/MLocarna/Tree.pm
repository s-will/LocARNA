package MLocarna::Tree;

############################################################
##
## Class for representing trees with labeled nodes and branch lengths
##
## Trees can be
## - constructed from and written to NEWICK format
## - constructed by PGMA
## - traversed
## - partitioned
## - projected to leave subsets
##
############################################################

use strict;
use warnings;

use Text::ParseWords qw(parse_line) ; ## for tokenizing
use Scalar::Util qw(looks_like_number reftype);

use Data::Dumper;


# construct an empty Tree object
#
# @note constructor is overloaded
sub new {
    my $class = shift;
    my @args = @_;

    my $self = {};
    bless $self, $class;

    my $reftype0 = reftype($args[0]);

    if (@args == 0) {
        # pass
    } elsif (defined $reftype0) {
        if ($reftype0 eq "MLocarna::Tree") {
            # shallow copy
            my $tree = $args[0];
            for my $k (keys %{$tree}) {
                $self->{$k} = $tree->{$k};
            }
        } elsif ($reftype0 eq "ARRAY") {
            # assume construction from list of children, label, and length
            $self->{children} = $args[0];
            if (defined $args[1]) {
                $self->{label} = $args[1];
            }
            if (defined $args[2]) {
                $self->{length} = $args[2];
            }
        }
    } elsif ( $args[0] eq "LEAVE" ) {
        shift @args;
        $self->{label} = $args[0];
        if (defined $args[1]) {
            $self->{label} = $args[1];
        }
        if (defined $args[2]) {
            $self->{length} = $args[2];
        }
    } elsif ( $args[0] eq "UPGMA_DIST" ) {
        # assume construction by UPGMA from distance matrix
        shift @args;
        $self = _upgma_dist(@args);
    } elsif ( $args[0] eq "UPGMA" ) {
        # assume construction by UPGMA from similarity matrix
        shift @args;
        $self = _upgma(@args);
    } else {
        # assume construction from newick string
        $self = _read_newick($args[0]);
    }

    return $self;
}

# destroy tree object
sub DESTROY {
    my $self = shift;
}

# append child to tree object
sub append_child {
    my $self = shift;
    my $child = shift;
    if(!exists $self->{children}) {
        $self->{children} = [];
    }
    push @{ $self->{children} }, $child;
}

# convert tree object to newick string
sub to_newick {
    my $self = shift;

    my $newick = '';
    if (exists $self->{children}) {
        $newick = '(';
        my $children = $self->{children};
        for(my $i=0; $i<@$children-1;$i++) {
            $newick .= $children->[$i]->to_newick();
            $newick .= ',';
        }
        $newick .= $children->[-1]->to_newick();
        $newick .= ')';
    }
    if (exists $self->{label}) {
        $newick .= _quote_newick_label($self->{label});
    }
    if (exists $self->{length}) {
        $newick .= ':'.$self->{length};
    }

    return $newick;
}


########################################
## quote a newick label
##
## @param label unquoted label
## @returns quoted label
##
## quotes special characters, such that the text
## is valid as label in newick format
########################################
sub _quote_newick_label {
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
## unquote a newick label
##
## @params string
## @return unquoted string
##
sub _unquote_newick_label {
    my ($s)=@_;
    if ($s =~ /^\'(.*)\'$/) {
        $s = "$1";
        $s =~ s/\'\'/\'/g;
    }
    return $s;
}

########################################
## upgma_tree_dist( $names ref to list, $dist_matrix ref to matrix of distances)
##
## Compute an upgma-tree from distances and overwrite/initialize the tree
## object with the result
##
## @poram $self
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
sub _upgma_dist {
    my ($names, $dist_matrix) = @_;

    ## compute tree by upgma applied to dist matrix

    my @clusters; # a list of the clusters
    my @trees; # a list of sub-trees
    my @cluster_sizes;
    my @heights;

    for (my $i=0; $i<@$names; $i++) {
        $clusters[$i] = $i;
        $trees[$i]    = new MLocarna::Tree( $names->[$i] );
        $cluster_sizes[$i] = 1;
        $heights[$i]  = 0;
    }

    my $INFINITY = 1e10;

    my $index=1;

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
        my $ilen = sprintf("%.3f",$height - $heights[$cluster_i]);
        my $jlen = sprintf("%.3f",$height - $heights[$cluster_j]);

        $trees[$cluster_i]->{length} = $ilen;
        $trees[$cluster_j]->{length} = $jlen;

        $trees[$clusters[0]] = new MLocarna::Tree( [ $trees[$cluster_i], $trees[$cluster_j] ], $index++ );

        $cluster_sizes[$clusters[0]] = $cluster_sizes[$cluster_i]+$cluster_sizes[$cluster_j];
        $heights[$clusters[0]] = $height;

        $#clusters--;
    }

    return $trees[$clusters[0]];
}

########################################
## scores_to_dists( $score_matrix )
##
## @param score_matrix ref to a matrix of scores
##
## @return distance matrix
########################################
sub _scores_to_dists {

    sub _max {
        my ($x,$y) = @_;
        return $x>=$y ? $x : $y;
    }

    my $score_matrix = shift;

    my $n = scalar( @$score_matrix );
    my $max = $score_matrix->[0][0];
    for(my $i=0;$i<$n;$i++) {
        for (my $j=$i+1; $j<$n; $j++) {
            $max = _max( $score_matrix->[$i][$j], $max );
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
## Compute an upgma-tree from similarities (like locarna scores) and
## initialize/overwrite the tree object with it
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
sub _upgma {
    my ($names, $score_matrix) = @_;
    # convert scores to distances
    my $dist_matrix = _scores_to_dists($score_matrix);

    return _upgma_dist($names, $dist_matrix);
}

########################################
## Parse a newick tree format string and construct tree
##
## $newick tree in newick format
##
########################################
sub _read_newick {
    my $newick = shift;

    $newick =~ s/;$//; # allow that the tree string is terminated by ';'

    # tokenize line and preserve delimiters, also preserve quotes and backslashs
    my @tokens = parse_line('[,()\s:]', "delimiters", $newick);

    my @stack;

    my $i=0;
    while (1) {
        my $tok = $tokens[$i];
        $i++;

        last if !defined($tok);

        my $top=$stack[-1];

        if ($tok =~ /^\s*$/)  {
            # ignore whitespace
        } elsif ($tok eq "(") {
            # push open tree
            push @stack, new MLocarna::Tree;
        } elsif ($tok eq ")") {
            # close tree
            $#stack--;
            $stack[-1]->append_child($top);

            $tok = $tokens[$i];
            if (defined $tok and $tok =~ /^[a-zA-Z0-9]/ ) {
                $top=$stack[-1];
                $top->{label} = $tok;
                $i++;
            }
        } elsif ($tok eq ",") {
            ## go on to next child
            $#stack--;
            $stack[-1]->append_child($top);

        } elsif ($tok eq ":") {
            # parse distance
            $tok = $tokens[$i];
            if (!defined($tok)) {
                die("Parse error: unexpected end of newick tree string");
            }
            $i++;

            if ( ! looks_like_number($tok) ) {
                die("Distance expected in tree.");
            }
            $top->{length} = $tok;
        } else {
            push @stack, new MLocarna::Tree("LEAVE", _unquote_newick_label($tok));
        }
    }

    return $stack[-1];
}

sub is_leaf {
    my $self = shift;
    return @{$self->children()} == 0;
}

sub children {
    my $self = shift;

    my $children = [];

    if (exists $self->{children}) {
        $children = $self->{children};
    }
    return $children;
}

sub label {
    my $self = shift;
    if (exists $self->{label}) {
        return $self->{label};
    } else {
        return undef;
    }
}

# tree dfs traversal with post-callback
sub traverse {
    my $self = shift;
    my $postcallback = shift;

    my @results = ();
    for my $child (@{$self->children()}) {
        push @results, $child->traverse($postcallback);
    }
    return &$postcallback( $self, \@results );
}

# filter leaves of a tree
# @param self tree object
# @param function boolean function on trees
#
# keep only subtrees where function returns 1 on any leaf
sub filter {
    my $self = shift;
    my $function = shift;


    if ($self->is_leaf() and &$function( $self ) ) {
        return $self;
    } else {
        my @children;

        for my $child (@{$self->children()}) {
            my $filtered_child = $child->filter($function);
            if (defined $filtered_child) {
                push @children, $filtered_child;
            }
        }
        if ( @children != 0 ) {
            return new MLocarna::Tree(\@children,$self->{label},$self->{length});
        }
    }
    return undef;
}


########################################
## collect labels of leafs in a tree
##
## @param $tree tree object
sub leaf_labels {
    my $self = shift;
    my $labels = shift;

    my @labels;

    $self->traverse(sub {
       my $tree=shift;
       if( $tree->is_leaf() ){
            push @labels, $tree->label();
       }
    });

    return \@labels;
}

## check whether the tree contains a set of leaves
sub contains_leaves {
    my $self = shift;
    my $leaves = shift;

    my %labels_hash = map { $_ => undef } @{$self->leaf_labels()};
    for my $l (@$leaves) {
        if (!exists $labels_hash{$l}) {return 0;}
    }
    return 1;
}

## binarize tree
## rearrange a possibly multi-ary tree to be at most binary
## e.g. (A,B,C)L:3 becomes (((A,B)L:0),C:0)L:0),D)L:3
sub binarize {
    my $self = shift;
    return $self->traverse(sub {
        my $tree = shift;
        my $children = shift;
        if (@$children > 2) {
            my $newtree = new MLocarna::Tree( [$children->[0], $children->[1]] );
            if (exists $tree->{label}) {
                $newtree->{label}=$tree->{label};
            }
            if (exists $tree->{length}) {
                $newtree->{length}=0;
            }
            for(my $i=2; $i<@$children; $i++) {
                $newtree = new MLocarna::Tree( [ $newtree, $children->[$i] ] );
                if (exists $tree->{label}) {
                    $newtree->{label}=$tree->{label};
                }
                if (exists $tree->{length}) {
                    $newtree->{length}=0;
                }
            }
            if (exists $tree->{length}) {
                $newtree->{length} = $tree->{length};
            }
            return $newtree;
        } else {
            if (@$children>0) {$tree->{children}=$children;}
            return $tree;
        }
    });
}


########################################
## project a tree to a set of labels
##
## @param $self tree
## @param labels ref to list of labels
##
## @returns tree in postorder that is a sub-tree of @$tree
## and contains only leaves with labels in @$labels
##
sub project {
    my ($self,$labels) = @_;

    # create hash out of list labels for fast lookup
    my %labels_hash = map { $_ => undef } @$labels;

    return $self->filter(sub {
        my $tree = shift;
        return (exists $labels_hash{$tree->label()});
    });
}

########################################
## Traverse all leaf partitions from splits along edges of the tree
##
## @param $self
## @param $callback function on leaf labels list
##
########################################
sub traverse_partitions {
    my $self = shift;
    my $callback = shift;

    my $ignore=0;
    if ($self->children()==2) {
        $ignore = $self->children()->[1];
    }

    $self->traverse(sub {
        my $tree = shift;
        if ($tree!=$self and $tree!=$ignore) {
            &$callback($tree->leaf_labels());
        }
    });
}

return 1;
