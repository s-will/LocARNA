package MLocarna::NameNormalizer;

use strict;
use warnings;

use MLocarna::Aux;

########################################
## Sequence names in fasta and clustalw files can contain characters
## that cannot be written to disk.  Therefore, we introduce name
## normalization that removes special characters and restrict name
## length in order to generate a nice filename.
##

sub new {
    my ($class, $seqs) = @_;

    my $self = {
        ## global hash for storing the association of names to normalized names
        _norm_names => { }
    };

    bless $self, $class;

    $self->register_names($seqs);

    return $self;
}

sub DESTROY {
    my $self = shift;
}

########################################
## _normalize($name, @names list of existing names)
##
## generate a sequence name from $name that
## has at most a length of 16
## and
## does not already exist in @names
##
########################################
sub _normalize {
    my ($self, $name, @names) = @_;

    my $maxlen=16;

    chomp $name;

    # replace all non-alpha-numeric symbols by '_'
    $name =~ s/[^a-zA-Z\d]/_/g;
    # take first 16 characters
    $name = substr $name,0,$maxlen;

    # make $name unique if it already occurs in @names
    # by appending '_' and a number $i to the truncated name
    # ($name is truncated such that maxlen is not exceeded)
    #
    # iterate over numbers $i from 1 until a unique name is generated
    my $i=1;
    while (grep /^$name$/, @names) {
        my $arity = int(log($i)/log(10))+1;
        if ($arity > 10) { # this will never happen ;)
            die "Could not generate unique name";
        }

        $name = substr $name,0,$maxlen-$arity-1;
        $name = sprintf("%s_%0$arity"."d",$name,$i);
        $i++;
    }
    return $name;
}


########################################
## print normalized sequence names hash to stdout
##
sub print {
    my $self = shift;
    foreach my $k (sort keys %{ $self->{_norm_names} }) {
        print "$k => $self->{_norm_names}->{$k}\n";
    }
}

########################################
## register_name( $name )
##
## registers a normalized name for $name
## if the name exists already, do nothing
##
## @param $name a sequence name
##
sub register_name {
    my ( $self, $name ) = @_;

    if ( ! exists $self->{_norm_names}->{$name} ) {
        $self->{_norm_names}->{$name} =
          $self->_normalize( $name, values %{ $self->{_norm_names} } );
    }
}

########################################
## register_names( $loh )
##
## registers all normalized name for names in a loh with name entries
## by calling register_normalized_seqname( $name )
##
## @param $loh list of hashs with entry name
##
sub register_names {
    my ( $self, $loh ) = @_;

    foreach my $h (@$loh) {
        $self->register_name( $h->{name} );
    }
}

########################################
## forget()
##
## forget all normalized names
##
sub forget {
    my $self = shift;
    $self->{_norm_names} = {};
}

########################################
## nname( $name )
## returns normalized name for the given name
##
## @returns normalized name for $name
##
## if normalized name was not registered: prints error message and exits with error code -1.
sub nname {
    my ( $self, $name ) = @_;

    if (exists $self->{_norm_names}->{$name}) {
        return $self->{_norm_names}->{$name};
    } else {
        MLocarna::Aux::printerr "ERROR: No normalized name was registered for the requested name $name.\n";
        $self->print();
        exit(-1);
    }
}

##
## generate a unique string from two sequence names
## that can be used as a filename
##
## normalized names have to be registered
##
## normalized names do not containt '-' (see MLocarna::chp !)
sub nnamepair {
    my ($self, $nameA, $nameB) = @_;

    my $nnameA=$self->nname($nameA);
    my $nnameB=$self->nname($nameB);

    return MLocarna::chp($nnameA,$nnameB);
}


return 1;
