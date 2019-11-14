package MLocarna::SparseMatrix;

############################################################
##
## package of functions for handling 2D and 4D sparse matrices that
## are implemented as hashs of hashs
##
############################################################

use 5.008003;
use strict;
use warnings;

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
    add_2D_inplace
    add_4D_inplace
    divide_2D
    divide_2D_inplace
    divide_4D
    filter_2D
    filter_4D
    print_2D
    print_4D
    read_2D
    read_4D
    scale_2D
    scale_4D
    transpose_2D
    transpose_4D
    to_array_2D
    to_sym_array_2D
    to_sym_array_check_2D
    write_2D
    write_4D
);


############################################################
## operations on 2- and 4-dimensional sparse matrices that are
## represented as hashs of hashs
##

########################################
## add_2D_inplace($m1,$m2)
##
## Add a sparse matrices in-place
##
## @post $m1 = $m1 + $m2
##
########################################
sub add_2D_inplace {
    my ($m1,$m2) = @_;

    foreach my $i (keys %$m2) {
        foreach my $j (keys %{ $m2->{$i} } ) {
            $m1->{$i}{$j} += $m2->{$i}->{$j};
        }
    }
}

########################################
## add two 4-dimensional sparse matrices
##
## Add a 4D sparse matrices in-place
##
## @post $m1 = $m1 + $m2
########################################
sub add_4D_inplace {
    my ($m1,$m2) = @_;

    foreach my $i (keys %$m2) {
        foreach my $j (keys %{ $m2->{$i} }) {
            foreach my $k (keys %{ $m2->{$i}{$j} }) {
                foreach my $l (keys %{ $m2->{$i}{$j}{$k} }) {
                    $m1->{$i}{$j}{$k}{$l} += $m2->{$i}{$j}{$k}{$l};
                }
            }
        }
    }
}

sub divide_2D {
    my ($m,$divisor) = @_;
    my %r;

    foreach my $i (keys %$m) {
        foreach my $j (keys %{ $m->{$i} }) {
            $r{$i}{$j} = $m->{$i}{$j} / $divisor;
        }
    }

    return \%r;
}

########################################
## divide_2D_inplace($m,$divisor)
##
## @post $m /= $divisor
##
########################################
sub divide_2D_inplace {
    my ($m,$divisor) = @_;

    foreach my $i (keys %$m) {
        foreach my $j (keys %{ $m->{$i} }) {
            $m->{$i}{$j} = $m->{$i}{$j} / $divisor;
        }
    }
}


sub divide_4D {
    my ($m,$divisor) = @_;
    my %r;

    foreach my $i (keys %$m) {
        foreach my $j (keys %{ $m->{$i} }) {
            foreach my $k (keys %{ $m->{$i}{$j} }) {
                foreach my $l (keys %{ $m->{$i}{$j}{$k} }) {
                    $r{$i}{$j}{$k}{$l} = $m->{$i}{$j}{$k}{$l} / $divisor;
                }
            }
        }
    }

    return \%r;
}


########################################
## filter_2D($m,$threshold)
##
## filter 2-dim matrix by probability threshold
##
## return result sparsematrix
########################################
sub filter_2D {
    my ($m,$threshold) = @_;

    my %r;

    foreach my $i (keys %$m) {
        foreach my $j (keys %{ $m->{$i} }) {
            if ( $m->{$i}{$j} >= $threshold ) {
                $r{$i}{$j} = $m->{$i}{$j};
            }
        }
    }

    return \%r;
}

########################################
## filter_2D($m,$threshold)
##
## filter 4-dim matrix by probability threshold
##
## return result sparsematrix
########################################
sub filter_4D {
    my ($m,$threshold) = @_;

    my %r;

    foreach my $i (keys %$m) {
        foreach my $j (keys %{ $m->{$i} }) {
            foreach my $k (keys %{ $m->{$i}{$j} }) {
                foreach my $l (keys %{ $m->{$i}{$j}{$k} }) {
                    if ( $m->{$i}{$j}{$k}{$l} >= $threshold ) {
                        $r{$i}{$j}{$k}{$l} = $m->{$i}{$j}{$k}{$l};
                    }
                }
            }
        }
    }

    return \%r;
}

sub scale_2D {
    my ($m,$scale) = @_;

    my %r;

    foreach my $i (keys %$m) {
        foreach my $j (keys %{ $m->{$i} }) {
            $r{$i}{$j} = $m->{$i}{$j} * $scale;
        }
    }

    return \%r;
}

sub scale_4D {
    my ($m,$scale) = @_;

    my %r;

    foreach my $i (keys %$m) {
        foreach my $j (keys %{ $m->{$i} }) {
            foreach my $k (keys %{ $m->{$i}{$j} }) {
                foreach my $l (keys %{ $m->{$i}{$j}{$k} }) {
                    $r{$i}{$j}{$k}{$l} = $m->{$i}{$j}{$k}{$l} * $scale;
                }
            }
        }
    }

    return \%r;
}


# read from a sparse matrix file as written by locarna --write-match-probs
sub read_2D {
    my ($file) = @_;

    open(my $SM_IN, "<", $file) || die "Cannot read from $file: $!";

    my %h;

    while( my $line=<$SM_IN> ) {
        if ( $line =~ /(\d+) (\d+) ([\d.e+-]+)/ ) {
            $h{$1}{$2} = $3;
        }
    }

    close $SM_IN;

    return \%h;
}

# read from a sparse matrix file as written by locarna --write-match-probs
sub read_4D {
    my ($file) = @_;

    open(my $SM_IN, "<", $file) || die "Cannot read from $file: $!";

    my %h;

    while( my $line=<$SM_IN> ) {
        if ( $line =~ /(\d+) (\d+) (\d+) (\d+) ([\d.e+-]+)/ ) {
            $h{$1}{$2}{$3}{$4} = $5;
        }
    }

    close $SM_IN;

    return \%h;
}

sub write_2D {
    my ($m,$file)=@_;

    open(my $SM_OUT,">", "$file") || die "Cannot write to $file: $!";

    foreach my $i ((keys %$m)) {
        foreach my $j ((keys %{ $m->{$i} })) {
            print $SM_OUT "$i $j $m->{$i}{$j}\n";
        }
    }

    close $SM_OUT;
}

sub print_2D {
    my ($m)=@_;

    foreach my $i ((keys %$m)) {
        foreach my $j ((keys %{ $m->{$i} })) {
            print "$i $j $m->{$i}{$j}\n";
        }
    }
}

sub write_4D {
    my ($m, $file)=@_;

    open(my $SM_OUT,">", "$file") || die "Cannot write to $file: $!";

    foreach my $i ((keys %$m)) {
        foreach my $j ((keys %{ $m->{$i} })) {
            foreach my $k ((keys %{ $m->{$i}{$j} })) {
                foreach my $l ((keys %{ $m->{$i}{$j}{$k} })) {
                    print $SM_OUT "$i $j $k $l $m->{$i}{$j}{$k}{$l}\n";
                }
            }
        }
    }

    close $SM_OUT;
}

sub print_4D {
    my ($m, $file)=@_;

    foreach my $i ((keys %$m)) {
        foreach my $j ((keys %{ $m->{$i} })) {
            foreach my $k ((keys %{ $m->{$i}{$j} })) {
                foreach my $l ((keys %{ $m->{$i}{$j}{$k} })) {
                    print "$i $j $k $l $m->{$i}{$j}{$k}{$l}\n";
                }
            }
        }
    }
}


# transpose a sparse matrix implemented as hash of hashs
sub transpose_2D {
    my ($m) = @_;

    my %r;

    foreach my $i (keys %$m) {
        foreach my $j (keys %{ $m->{$i} }) {
            $r{$j}{$i} = $m->{$i}{$j};
        }
    }

    return \%r;
}

# transpose a sparse matrix implemented as hash of hashs
sub transpose_4D {
    my ($m) = @_;

    my %r;

    foreach my $i (keys %$m) {
        foreach my $j (keys %{ $m->{$i} }) {
            foreach my $k (keys %{ $m->{$i}{$j} }) {
                foreach my $l (keys %{ $m->{$i}{$j}{$k} }) {
                    $r{$k}{$l}{$i}{$j} = $m->{$i}{$j}{$k}{$l};
                }
            }
        }
    }

    return \%r;
}

## convert a 2D sparse matrix to an 2D array
## @param $hash 2D hash
## @param $matrix if given add to this matrix
## @return 2D matrix
sub to_array_2D {
    my $hash = shift;
    my $matrix = shift;
    for my $i (keys %$hash) {
        for my $j (keys %{ $hash->{$i} }) {
            $matrix->[$i][$j] = $hash->{$i}{$j};
        }
    }

    return $matrix;
}

## convert a 2D sparse matrix to a symmetric 2D array
## @param $hash 2D hash
## @param $matrix if given add to this matrix
## @return 2D matrix
sub to_sym_array_2D {
    my $hash = shift;
    my $matrix = shift;
    for my $i (keys %$hash) {
        for my $j (keys %{ $hash->{$i} }) {
            $matrix->[$i][$j] = $hash->{$i}{$j};
            $matrix->[$j][$i] = $hash->{$i}{$j};
        }
    }

    return $matrix;
}

## convert a 2D sparse matrix to a symmetric 2D array
## @param $hash 2D hash
## @param $matrix if given add to this matrix
## @param $overwrite accumulator for overwrite flag
## @return whether overwrite took place
sub to_sym_array_check_2D {
    my $hash = shift;
    my $matrix = shift;
    my $overwrite = shift;
    for my $i (keys %$hash) {
        for my $j (keys %{ $hash->{$i} }) {
            $overwrite |= defined $matrix->[$i][$j];
            $overwrite |= defined $matrix->[$j][$i];

            $matrix->[$i][$j] = $hash->{$i}{$j};
            $matrix->[$j][$i] = $hash->{$i}{$j};
        }
    }

    return $overwrite;
}

## ------------------------------------------------------------
1;
