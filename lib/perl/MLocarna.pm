package MLocarna;

use 5.008003;
use strict;
use warnings;


#use FindBin;
#use lib $FindBin::RealBin;
use MLocarna::SparseMatrix;
use MLocarna::MatchProbs;
use MLocarna::Aux;

use Cwd;

require Exporter;

# set the version for version checking
our $VERSION     = 1.00;

our @ISA         = qw(Exporter);
our @EXPORT      =
  qw(
        %loh_associate_nnames_to_names
        %loh_translate_names
        alifold_mfe
        alifold_pf
        alifold_structure
        aln_h2loh
        aln_length_atleastonematch
        aln_to_alnloh
        anchor_constraint_string
        compute_alignment_from_seqs
        compute_alignment_score
        constraint_annotation_is_valid_or_die
        constrain_sequences_from_reliable_structures
        convert_alifold_dp_to_pp
        convert_dp_to_pp_with_constraints
        convert_fix_structure_to_pp
        extract_from_clustal_alignment
        find_in_exec_path
        find_in_exec_path_or_error
        loh_names
        loh_sort
        new_intermediate_name
        parse_bracket_structure
        parse_bracket_structure_single
        project_aln
        project_alnloh
        project_string_to_alignment_sequence
        project_structure_to_alignment_sequence
        read_2D_matrix
        read_aln_wo_anno
        read_anchor_constraints
        read_clustalw_alignment
        read_clustalw_aln
        read_clustalw_alnloh
        read_dp_ps
        read_fasta
        read_pipe
        read_pp_file_aln_w_anno
        read_pp_file_aln_wo_anno
        read_pp_file_pairprob_info
        read_pp_file_pairprobs
        sprint_fasta_alnloh
        system_pipein_list
        write_2D_matrix
        write_aln
        write_clustalw_alnloh
        write_clustalw_loh
        write_pp
   );

our %EXPORT_TAGS = ();

# your exported package globals go here,
# as well as any optionally exported functions
our @EXPORT_OK   =
  qw(
        $PACKAGE_STRING
        parse_mfasta
        parse_mfasta_constraints
   );

our $RNAalifold = "RNAalifold";

## alifold specific options
## turn on ribosum scoring and use "best" factors from the "improved alifold" paper
our @RNAalifold_options = ("--ribosum_scoring", "--cfactor" => "0.6", "--nfactor" => "0.5");

our $PACKAGE_STRING = "MLocarna";



########################################
## new_intermediate_name
##
## returns a new name for intermediate files
##
## These names are used by mlocarna as filenames/identifiers of
## intermediary alignments during the computation of progressive and
## iterative alignments
##
########################################
sub new_intermediate_name {
    my $data = shift;
    my $label = shift;
    my $baseone = shift;

    my $intermediate_name_base="intermediate";

    my $name = "$intermediate_name_base";
    if (defined($label)) {
        $label =~ s/[^A-Za-z0-9]/_/g;
        $name="$name$label";
    }

    my $suf="";
    my $i=0;
    if (defined($baseone)) {
        $i=1;
        $suf="-1";
    }

    while( exists $data->{"$name$suf"} ) {
        $i++;
        $suf="-$i";
    }
    $name="$name$suf";
    $data->{$name}=undef;

    return $name;
}


########################################
## read_fasta($filename)
## ------------------------------
##
## Read an (extended) fasta file. If filename ends with "*.gz" the
## file is automatically uncompressed by gunzip.  Dies if file is not
## readable or gunzippable.  The result is returned as list of the
## sequences.
##
##
## Each sequence is encoded as hash of features
##   name: id of sequence
##   descr: description of the sequence
##   seq:  sequence string
##
## supports special annotation strings as used by locarna
## lines ending in #(\S+) are recognized and
## returned as features "ANNO#$1" of the sequence.
##
## return reference to list of hash representation
##
########################################
sub read_fasta {
    my $filename = shift;

    my $fh;

    if ($filename =~ /\.gz$/) {
	open($fh,"gunzip -c $filename |") || die "Cannot read or uncompress file $filename.";
    } else {
	open($fh,$filename) || die "Cannot read file $filename.";
    }

    my @fasta = ();
    my %seen_names;

    my $line=<$fh>;
    while(defined($line)) {
	if ($line=~/^>\s*(\S+)\s*(.*)/) {
	    my $name=$1;
	    my $description=$2;

	    ## check for duplicate names in fasta
	    if (exists $seen_names{$name}) {
		printerr "Duplicate name \"$name\" in fasta input. ";
		my $bar="some text here";
		if (length($description)>0) {
		    $bar=$description;
		}
		printerr "Note that in \">$name $bar\", only \"$name\" is the name, ";
		printerr "whereas the rest of the line \"$bar\" (after the blank)"
                  ." is interpreted as description.\n";
		exit(-1);
	    }
	    $seen_names{$name}=1;

	    my $seq = { name  => $name,
			descr => $description };

	    while (defined($line=<$fh>) && ($line !~ /^>/)) {
		chomp $line;
		$line =~ s/\s+//g;

		if  ($line =~ /(.+)\s*\#(.+)/) {
		    $seq->{"ANNO#$2"} .= $1;
		} else {
		    $seq->{seq} .= $line;
		}
	    }

	    push @fasta, $seq;
	} else {
	    $line=<$fh>;
	}
    }

    close $fh;

    return \@fasta;
}



########################################
## constraint_annotation_is_valid_or_die($fasta)
##
## Checks validity of anchor and structure constraint description in
## fasta list of hash representation.
##
## $fasta list of hashs representation
##
## require that all anchor constraints have the same number of row (if
## any anchor constraints are given!)
##
## allow mixing of sequences with and without anchor constraint
## annotation (was forbidden <=1.8.8)
##
## die with error message when constraint annotation is invalid
##
########################################
sub constraint_annotation_is_valid_or_die {
    my $fasta=shift;

    my $numseqconstraints = -1;

    for my $seq (@$fasta) {
	if (exists $seq->{"ANNO#S"}) {
	    if (length($seq->{"ANNO#S"}) != length($seq->{seq})) {
		die "Structure constraint length unequals sequence length for "
		    .$seq->{name}.".";
	    }
	}

	my $i=1;
	while (exists $seq->{"ANNO#$i"}) {
	    if (length($seq->{"ANNO#$i"}) != length($seq->{seq})) {
		die "Sequence constraint length unequals sequence length for "
		    .$seq->{name}.".";
	    }
	    $i++;
	}
        $i=$i-1; ## set i to number of rows in anchor constraint
        if ($i>0) { ## ignore sequences without anchor constraints
            if ($numseqconstraints == -1) {$numseqconstraints=$i;}

            if ($i != $numseqconstraints) {
                die "Bad sequence constraints for sequence "
                  .$seq->{name}." -- expected $numseqconstraints anchor constraint row(s)";
            }
        }
    }
}

## print error and exit(-1)
sub error_exit {
    my $msg = shift;
    printerr $msg;
    exit(-1);
}

##maximization
sub max {
    my ($x,$y)=@_;
    return ($x>=$y)?$x:$y;
}

########################################
## read in bed file, check 4-column format, check existence of
## sequence names
sub read_bed4 {
    my $filename = shift;

    open(my $fh, $filename)
      || error_exit "ERROR: Cannot open bed file for reading.\n";

    my @bedentries=(); # entries in the bed file
    my %seqnames=(); # set of sequence names in the bed file
    while(<$fh>) {
        my @F=split /\s+/;
        if (@F != 4) {
            error_exit "ERROR: bed file $filename has wrong format. "
              ."Require four-column bed format (line "
              .(1 + scalar @bedentries).").\n";
        }
        if ($F[1]>$F[2]) {
            error_exit "ERROR: in bed file $filename, start position greater than end position "
              ."in line\n  ".(1 + scalar @bedentries).": $_\n";
        }
        $seqnames{$F[0]}=1;
        push @bedentries, [ @F ];
    }

    return \@bedentries;
}

## generate lookup hash from list of hashs
sub gen_lookup_hash_loh {
    my $loh = shift;
    my $key = shift;

    my %lookup_hash = ();

    my $idx=0;
    for my $entry (@$loh) {
        $lookup_hash{$entry->{$key}}=$idx++;
    }
    return \%lookup_hash;
}



## check whether region name does not occur right of current index in any bin
sub is_notgreater_name {
    my $regname = shift;
    my $rbins = shift;
    my $rbinmaps = shift;
    my $rcur_index = shift;

    for my $bin (keys %$rbins) {
        if ( exists $rbinmaps->{$bin}{$regname}
             and $rbinmaps->{$bin}{$regname} > $rcur_index->{$bin} ) {
            return 0;
        }
    }
    return 1;
}


########################################
## assign indices to region names
##
## Indices are assigned ascending in alignment order.  The alignment
## order is the order of regions enforced by the requirement to occur
## in a common multiple alignment (from smaller to larger positions).
##
## Failure policy:
##   Exit with error
##   * if there is no alignment order
##   * if regions overlap
##
sub assign_region_indices {

    # get bedentries and make copy
    my $bedentries_in = shift;
    my $bedentries = [ @{ $bedentries_in } ];

    # constants: positions in (extended) bed entry
    my $SEQNAME=0;
    my $START=1;
    my $END=2;
    my $REGNAME=3;
    my $LINENO=4;

    # add line numbers to bed entries (to keep track of input order)
    for (my $i=0; $i < @$bedentries; $i++) {
        my @entry = @{ $bedentries->[$i] };
        push @entry, $i;
        $bedentries->[$i] = \@entry;
    }

    for (my $i=0; $i < @$bedentries; $i++) {
        my @entry = @{ $bedentries->[$i] };
    }

    # bin bedentries by sequences
    my %bins = ();
    for my $bedentry (@$bedentries) {
        my @bin = ();
        if (exists $bins{$bedentry->[$SEQNAME]}) {
            @bin = @{ $bins{$bedentry->[$SEQNAME]} };
        }
        push @bin, $bedentry;
        $bins{$bedentry->[$SEQNAME]} = \@bin;
    }


    # sort bins by first position
    for my $bin (keys %bins) {
        $bins{$bin} = [ sort {$a->[$START] <=> $b->[$START]} @{ $bins{$bin} } ];
    }

    # print "Sorted bins:\n ";
    # for my $bin (keys %bins) {
    #     print "BIN $bin ";
    #     for my $x (@{$bins{$bin}}) {print ", @$x";}
    #     print "\n";
    # }

    # detect overlap
    for my $bin (keys %bins) {
        my $last_end=0;
        for my $entry (@{$bins{$bin}}) {
            if ($entry->[$START] < $last_end) {
                error_exit ("ERROR: region \"$entry->[$REGNAME]\" overlaps previous region".
                            " for sequence \"$entry->[$SEQNAME]\" in bed anchor constraints.\n");
            }
            $last_end = $entry->[$END];
        }
    }


    # create name to index map for each bin
    my %binmaps = ();
    for my $bin (keys %bins) {
        my %n2i;
        for (my $i=0; $i < @{ $bins{$bin} }; $i++) {
            $n2i{ $bins{$bin}->[$i]->[$REGNAME] } = $i;
        }
        $binmaps{$bin} = \%n2i;
    }

    # print("Name to index map\n");
    # for my $bin (keys %bins) {
    #     print "  BIN_N2I $bin ";
    #     for my $n (keys %{$binmaps{$bin}}) {print ", $n:".$binmaps{$bin}->{$n};}
    #     print "\n";
    # }

    # merge bins to construct order or detect inconsistency
    my $region_name_index = 0;
    my %region_name_map  = ();

    my %cur_index; #index for each bin
    for my $bin (keys %bins) {$cur_index{$bin}=0;} #initialize

    while (1) {
        # determine smallest region name
        #  * does not occur right of cur index in any bin (is_notgreater_name)
        #  * has smallest line number
        my $min_name = undef;
        my $min_lineno = @{ $bedentries };

        my $done=1;
        for my $bin (keys %bins) {
            if ( $cur_index{$bin} >= @{ $bins{$bin} } ) {next;}
            $done = 0;

            my $cur_entry = $bins{$bin}->[$cur_index{$bin}];
            if ( is_notgreater_name($cur_entry->[$REGNAME], \%bins, \%binmaps, \%cur_index) ) {
                if ( $cur_entry->[$LINENO] < $min_lineno ) {
                    $min_lineno = $cur_entry->[$LINENO];
                    $min_name = $cur_entry->[$REGNAME];
                }
            }
        }
        if ($done) {last;}

        ## detect inconsistencies
        if ( not defined($min_name) ) {
            error_exit "ERROR: Conflicting anchors in bed file.\n";
        }
        if (exists $region_name_map{ $min_name } ) {
            error_exit "ERROR: Inconsistent anchors in bed file.\n";
        }

        ## advance current indices of bins where current name is $min_name
        for my $bin (keys %bins) {
            if ( $cur_index{$bin} >= @{ $bins{$bin} } ) {next;}
            my $cur_entry = $bins{$bin}->[$cur_index{$bin}];
            if ($cur_entry->[$REGNAME] eq $min_name) {
                $cur_index{$bin}++;
            }
        }
        ## register min name
        $region_name_map{ $min_name } = $region_name_index++;
    }

    return \%region_name_map;
}

########################################
## read_anchor_constraints($seqs,$filename)
##
## Read anchor constraints from bed file and set them in $seqs
##
## @param $seqs     sequences (list of hashs)
## @param $filename bed file with anchor constraints
##
## @return ref to list of anchor region names in the order of assigned
## anchor ids
## @result set anchor constraints in $seqs
##
## Anchor constraints in four-column bed format specify positions of
## named anchor regions per sequence. The 'contig' names have to
## correspond to the fasta input sequence names. Anchor names must be
## unique per sequence and regions of the same name for different
## sequences must have the same length. This constrains the alignment
## to align all regions of the same name.
##
## The specification of anchors via this option removes all anchor
## definitions that may be given directly in the fasta input file!
##
## Details: region names receive name numbers in the order the names
## have to appear in an alignment; inconsistencies are detected. As
## well, positions in a region receive position numbers from left to
## right. Numbering starts with 0. Anchor names are than composed of
## region number and position number, both with leading 0s, where the
## number of digits is choosen minimally to encode all occuring anchor
## names with the same length.  Finally the anchor names are encoded
## as anchor constraint lines in the same way LocARNA accepts them in
## the fasta input file.
##
## FAILURE policy:
## exit with error message,
## * if file cannot be read or has the wrong format.
## * if bed file does not contain >=2 sequence names in $seqs
## * if some sequence names in the bed file do not exist in $seqs
## * if some anchor regions is out of range
## * if anchor regions conflict with each other (cross or overlap)
## * if regions of same name have different lengths
##
## @pre names in $seqs are unique
##
########################################
sub read_anchor_constraints {
    my $seqs = shift;
    my $filename = shift;

    # read entries in the bed file
    my $bedentries = read_bed4($filename);

    # generate set of sequence names
    my %seqnames=(); # set of sequence names in the bed file
    foreach my $F (@$bedentries) {$seqnames{$F->[0]}=1;}

    ## count sequences with specification; fill name->idx hash
    my $seq_idx_by_name = gen_lookup_hash_loh($seqs,'name');

    ## --------------------
    ## check input errors
    my $specified_seqs=0;
    my $unknown_seqs="";
    for my $k ( keys %$seq_idx_by_name ) {
        if (exists $seqnames{$k}) {
            $specified_seqs++;
        } else {
            $unknown_seqs="$unknown_seqs $k";
        }
    }
    if ($specified_seqs<2) {
        error_exit
          "ERROR: $specified_seqs sequence name(s) in the anchor "
          ."constraint file match the fasta input (required >=2).\n";
    }
    if ( $unknown_seqs ne "" ) {
        error_exit "ERROR: anchor constraint file contains specifications "
          ."for unknown sequence(s):\n\t$unknown_seqs\n";
    }

    # check equal length of regions
    {
        my %region_lengths=();
        foreach my $F (@$bedentries) {
            my $len = $F->[2]-$F->[1]+1;
            if ( exists $region_lengths{$F->[3]} ) {
                if ( $region_lengths{$F->[3]} != $len ) {
                    error_exit "ERROR: Region ".$F->[3]." occurs with ".
                      "different lengths in anchor constraints.\n";
                }
            } else {
                $region_lengths{$F->[3]} = $len;
            }
        }
    }

    ## end check input errors
    ## --------------------

    ## assign numbers to region names
    my $region_names = assign_region_indices($bedentries);

    my $region_id_width = length((keys %$region_names)-1);

    # determine longest range
    my $position_id_width=0;
    for my $bedentry (@$bedentries) {
        my $region_len=$bedentry->[2]-$bedentry->[1];
        $position_id_width=max($position_id_width,length($region_len-1));
    }

    ########################################
    ## generate anchor constraint lines
    #
    ## 0) remove existing anchor lines
    for my $seq (@{$seqs}) {
        for my $key (keys %{$seq}) {
            if ( $key =~ /ANNO#(\d+)/ ) {
                $seq->{$key}=undef;
            }
        }
    }
    #
    ## 1) generate empty lines
    for my $seq (@{$seqs}) {
        for (my $i=1; $i<=$region_id_width+$position_id_width; $i++) {
            my $key="ANNO#$i";
            my $seq_length=length($seq->{seq});
            $seq->{$key}='.'x$seq_length;
        }
    }
    #
    ## 2) set anchors
    for my $bedentry (@$bedentries) {
        next unless exists $seq_idx_by_name->{$bedentry->[0]};

        my $seq=$seqs->[$seq_idx_by_name->{$bedentry->[0]}];
        my $rid = $region_names->{$bedentry->[3]};
        my @rid = split //, sprintf("%0${region_id_width}d",$rid);

        if ($bedentry->[1]<0 || $bedentry->[2]>length($seq->{'seq'})) {
            error_exit("ERROR: region out of range in anchor constraints".
                       " file ($bedentry->[0], $bedentry->[3]).\n");
        }

        for (my $p=$bedentry->[1]; $p<$bedentry->[2]; $p++) {
            ## set name at pos $p for region $rid
            for ( my $i=1; $i<=$region_id_width; $i++ ) {
                substr($seq->{"ANNO#$i"},$p,1) = $rid[$i-1];
            }
            my @pos =
              split //, sprintf("%0${position_id_width}d",$p-$bedentry->[1]);
            for ( my $i=1; $i<=$position_id_width; $i++ ) {
                substr($seq->{"ANNO#".($i+$region_id_width)},$p,1)
                  = $pos[$i-1];
            }
        }
    }

    ## get list of region names and return it
    my @regions_names_list;
    for my $key (keys %$region_names) {
        $regions_names_list[$region_names->{$key}]=$key;
    }
    return \@regions_names_list;
}

########################################
## anchor_constraint_string($seq hash ref)
##
## Collect sequence constraint string for $seq entry of a fasta list of hashs
##
## returns string "<c1>#<c2>#...<cN>", where ci is the constraints description line i
##
########################################+
sub anchor_constraint_string {
    my $seq=shift;
    my $i=1;

    my $str="";
    while (exists $seq->{"ANNO#$i"}) {
	$str .= $seq->{"ANNO#$i"}."#";
	$i++;
    }
    $str =~ s/\#$//;

    return $str;
}

########################################
# DEPRECATED
# read a (multiple) fasta file and return
# a hash that maps the sequence names to their
# corresponding sequences
#
sub parse_mfasta {
    my ($file) = @_;
    my %mfasta;

    local *PMF_IN;

    printerr "Use of parse_mfasta is deprecated. Use read_fasta instead.";

    open(PMF_IN,$file) || die "Cannot read mfasta file $file\n";

    my $line=<PMF_IN>;
    while(defined($line)) {
	if ($line=~/^>\s*(.+)/) {
	    my $name=$1;
	    chomp $name;
	    $name =~ s/\s+$//g; ## accept and ignore trailing spaces in name
	    $name =~ s/\s+.+$//g; ## ignore everything after first blank
	    $name =~ s/[^a-zA-Z\d]/_/g;
	    $name = substr $name,0,16;

	    if (exists $mfasta{$name}) {
		printerr "Duplicate name in mfasta: $name\n";
		exit -1;
	    }
	    while (defined($line=<PMF_IN>) && ($line !~ /^>/)) {
		chomp $line;
		$line =~ s/\s+//g;
		$mfasta{$name} .= $line;
	    }
	} else {
	    $line=<PMF_IN>;
	}
    }
    return %mfasta;
}

########################################
# DEPRECATED
# read a (multiple) fasta file
# with constraint annotation
# also reads long names (the complete name in the fasta file)
#
# constraints are specified as strings of the same length
# as the sequence. The idea is:
#
#  - a string followed by #S on the same line defines
#    structure constraints (this is passed to RNAfold -p!)
#  - a string followed by #<k> is the k-th line of
#    the sequence constraint (this is passed to locarna)
#
# lines of different "types" can be mixed freely!
#
# return a hash that maps the sequence names to their
# corresponding sequences and constraints:
#
#  for a sequence name $name, we generate hash entries
#    - $mfasta{"$name"}
#    - $mfasta{"$name#S"}
#    - $mfasta{"$name#C"}
#    - $mfasta{"$name#LONG"}
#
sub parse_mfasta_constraints {
    my ($file) = @_;

    printerr "Use of parse_mfasta_constraints is deprecated. Use read_fasta instead.";

    my %mfasta;

    my %seqcons; # sequence constraints hash
    my %strcons; # sequence constraints hash

    local *PMF_IN;

    open(PMF_IN,$file) || die "Cannot read mfasta file $file\n";

    my $line=<PMF_IN>;
    while(defined($line)) {
	if ($line=~/^>\s*(.+)/) {

	    my $longname = $1;
	    $longname =~ s/\s+.*$//; # eat everything after first space in $longname

	    my $name = make_unique_seqname $longname,(keys %mfasta);

	    if (exists $mfasta{$name}) {
		printerr "Duplicate name in mfasta, cannot make unique : $name\n";
		exit -1;
	    }

	    $mfasta{"$name\#LONG"}=$longname;

	    while (defined($line=<PMF_IN>) && ($line !~ /^>/)) {
		chomp $line;
		$line =~ s/\s+//g;
		if ($line =~ /(.+)\#S/) {
		    $strcons{$name} .= $1;
		} elsif ($line =~ /(.+)\#(.+)/) {
		    $seqcons{"$name\#$2"} .= $1;
		} else {
		    $mfasta{$name} .= $line;
		}
	    }
	}
    }

    my @names = keys %mfasta;

    ## ----------------------------------------
    # insert constraints in mfasta hash
    # and check validity
    #
    my $has_seq_constraints=0; # whether there are seq constraints

    for my $name (@names) {
	if (exists $strcons{$name}) {
	    if (length $strcons{$name} != length $mfasta{$name}) {
		printerr "Structure constraint length unequals sequence length for $name.\n";
		exit -1;
	    }

	    $mfasta{"$name\#S"} = $strcons{$name};
	}

	my $i=1;
	while (exists $seqcons{"$name\#$i"}) {
	    if (length $seqcons{"$name\#$i"} != length $mfasta{$name}) {
		printerr "Sequence constraint length unequals sequence length for $name.\n";
		exit -1;
	    }

	    $mfasta{"$name\#C"} .= $seqcons{"$name\#$i"}."\#";
	    $i++;
	}
	if ($i>1) {
	    $mfasta{"$name\#C"} =~ s/\#$//; # remove last '#'
	    $has_seq_constraints=1;
	}

	my $num_seq_constraint_lines=grep /$name\#/, keys %seqcons;
	if ($i-1 != $num_seq_constraint_lines) {
	    printerr "Bad sequence constraints for sequence $name.\n";
	    exit -1;
	}
    }

    for my $name (@names) {
	if (($name!~/#/) &&  $has_seq_constraints) {
	    if (! exists $mfasta{"$name\#C"}) {
		printerr "Sequence constraint string missing for $name.\n";
	    }
	}
    }
    #
    ## ----------------------------------------

    return %mfasta;
}

## convert a dotplot file to a pp file
## thereby insert the sequence constraint strings
sub convert_dp_to_pp_with_constraints {
    my ($dpfile,$ppfile,$name,$sequence,$constraints,$read_condprobs) = @_;

    open(PP_OUT,">$ppfile") || die "Cannot open $ppfile for writing.";
    open(DP_IN,"$dpfile") || die "Cannot open $dpfile for reading.";

    print PP_OUT "#PP 2.0\n\n";
    print PP_OUT "$name $sequence\n";
    if (defined $constraints && $constraints ne "") {
	my @cs = split /\#/, $constraints;
	for (my $i=0; $i<@cs;$i++) {
	    print PP_OUT "#A".($i+1)." ".$cs[$i]."\n";
	}
    }

    print PP_OUT "\n#END\n";

    print PP_OUT "\n#SECTION BASEPAIRS\n\n";

    if ($read_condprobs) {
	print PP_OUT "\n#STACK\n";
    }

    my %pp;

    while (my $line=<DP_IN>) {
	if ($line =~ /^(\d+) (\d+) ([\d\.]+) ubox/) {
	    my $i=$1;
	    my $j=$2;
	    my $p=$3*$3;
	    $pp{$i}{$j}="$p";
	}

	if ($read_condprobs) {
	    # assume ubox always before lbox
	    if ($line =~ /^(\d+) (\d+) ([\d\.]+) lbox/) {
		my $i=$1;
		my $j=$2;
		my $p=$3*$3;
		if (! exists $pp{$i}{$j}) {
		    $pp{$i}{$j}="0";
		}
	    $pp{$i}{$j}.=" $p";
	    }
	}
    }

    close DP_IN;


    for my $i ( keys %pp ) {
	for my $j ( keys %{ $pp{$i} } ) {
	    print PP_OUT "$i $j $pp{$i}{$j}\n";
	}
    }

    print PP_OUT "\n\#END\n";

    close PP_OUT;
}

########################################
## convert alifold dot plot to pp file
##
## @param dot plot file
## @param clustalw alignment file
## @param pp output file
##
sub convert_alifold_dp_to_pp {
    my ($dpfile,$alnfile,$ppfile) = @_;

    my $aln = read_clustalw_alnloh($alnfile);

    my $PP_OUT;
    my $DP_IN;

    open($PP_OUT,">$ppfile") || die "Cannot open $ppfile for writing.";
    open($DP_IN,"$dpfile") || die "Cannot open $dpfile for reading.";

    print $PP_OUT "SCORE: 0\n\n";

    ## copy alignment from alnfile
    write_clustalw_alnloh($PP_OUT, $aln, 75, 0);

    print $PP_OUT "\n\#\n";

    while (my $line=<$DP_IN>) {
	if ($line =~ /(\d+) (\d+) ([\d\.]+) ubox/) {
	    my $i=$1;
	    my $j=$2;
	    my $p=$3*$3;

	    print $PP_OUT "$i $j $p\n";
	}
    }
    close $DP_IN;
    close $PP_OUT;
}


########################################
## extract_from_clustal_alignment($nameA,$nameB,$alignment)
##
## Parse an alignment in CLUSTAL-like format and return the alignment
## strings for sequences of the two names $nameA and $nameB.
##
## $nameA name of sequence A as string
## $nameB name of sequence B as string
## $alignment ref to list of lines of a CLUSTAL-like formatted alignment
##
## returns pair of alignment strings
##
########################################
sub extract_from_clustal_alignment {
    my ($nameA,$nameB,$alignment) = @_;

    my $aliA="";
    my $aliB="";

    for my $line (@$alignment) {
	if ($line =~ /^$nameA\s+(.+)/) {
	    $aliA.=$1;
	} elsif ($line =~ /^$nameB\s+(.+)/) {
	    $aliB.=$1;
	}
    }
    return ($aliA,$aliB);
}

########################################
## read_dp_ps($filename)
##
## Parse a "RNAfold -p"-generated dotplot postscript file
##
## $filename name of dot plot ps file
##
## @returns pair of sequence and reference to list of pair probabilities
##
########################################
sub read_dp_ps {
    my ($filename) = @_;
    local *IN;

    open(IN,$filename) || die "Cannot read $filename for parsing as dp-ps file.\n";

    my $seq="";
    my @pairprobs;

    while (my $line=<IN>) {
	if ($line =~ /^\/sequence \{ \(/) {
	    while (defined($line = <IN>) && ($line !~  /\} def/  ))  {
		chomp $line;
		$line =~ s/\\$//;
		$seq .= $line;
	    }
	    #print "read_dp_ps $filename: $seq\n";
	}

	if ($line =~ /(\d+)\s+(\d+)\s+(\S+)\s+ubox/) {
	    $pairprobs[$1-1][$2-1]=$3*$3;
	}
    }

    close IN;

    $seq ne "" || die "Empty sequence in dp.ps file $filename\n";

    return ($seq,\@pairprobs);
}


## read pp file and return the alignment w/o constraints
## auto-detect version of pp file
sub read_pp_file_aln_wo_anno {
    my ($filename)=@_;
    my %aln = read_pp_file_aln_w_anno($filename);

    for my $k (keys %aln) {
	if ($k =~ /^#/) { ## delete names beginning with '#'
	    delete $aln{$k};
	}
    }
    return %aln;
}

## read pp file and return the alignment w/ constraints
## auto-detect version of pp file
sub read_pp_file_aln_w_anno {
    my ($filename)=@_;
    my $pp_in;

    open($pp_in,$filename) || die "MLocarna::read_pp_file_aln: Cannot read $filename\n";

    my %aln;
    my %pairprobs;

    my $line;
    my $pp_version=1.0;

    if (($line = <$pp_in>) =~ /^#PP ([\d.]+)/ ) {
	$pp_version=$1;
	$line = <$pp_in>
    }

    while ($line) {
	$line = <$pp_in>;

	if ($pp_version>=2 && $line =~ /^\#END/) { last; }
	if ($pp_version<2 && $line =~ /^\#\s*$/) { last; }

	if (($line =~ /^(\S+)\s+(.+)/) && ($line !~ /^SCORE:/) ) {
	    $aln{$1}.=$2;
	}
    }

    close $pp_in;

    return %aln;
}

## read a pp file and return a hash of pair probabilities
##
## @param $filename name of the pp file
## @returns hash of pair probabilities reported in the file
##          the hash has keys "$i $j" and contains the probabilities p_ij (i<j)
##          positions in the hash are in [1..sequence length]
##
## auto-detect version of pp file
sub read_pp_file_pairprobs {
    my ($filename)=@_;
    my $pp_in;

    open($pp_in,$filename) || die "Can not read $filename\n";

    #my %aln;
    my %pairprobs = read_pp_file_pairprob_info($filename);

    foreach my $k (keys %pairprobs) {
	$pairprobs{$k} =~ /^(\S+)/;
	$pairprobs{$k} = $1;
    }

    return %pairprobs;
}

## read a pp file and return a hash of pair probability information
##
## @param $filename name of the pp file
## @returns hash of pair probability information reported in the file
##          the hash has keys "$i $j" and contains the probability information
##          for this pair;
##          this information is the rest of the line starting with "$i $j"
##          in the pp file
##          positions in the hash are in [1..sequence length]
##
## auto-detect version of pp file
sub read_pp_file_pairprob_info {
    my ($filename)=@_;
    my $pp_in;

    open($pp_in,$filename) || die "Can not read $filename\n";

    #my %aln;
    my %pairprobinfo;

    my $line;

    my $pp_version=1.0;
    if (($line = <$pp_in>) =~ /^#PP ([\d.]+)/ ) {
	$pp_version=$1;
	$line = <$pp_in>
    }

    while ($line = <$pp_in>) {
	if ($pp_version>=2 && $line =~ /^\#END/) {
	    while ($line = <$pp_in>) {
		if ($line =~ /^\#SECTION BASEPAIRS/) {
		    last;
		}
	    }
	    last;
	}
	if ($pp_version<2 && $line =~ /^\#\s*$/) { last; }
    }

    while (($line = <$pp_in>)) {
	if ($line =~ /^(\S+)\s+(\S+)\s+(.+)/) {
	    my $i=$1;
	    my $j=$2;
	    my $pi=$3;
	    $pairprobinfo{"$i $j"}=$pi;
	}
	if ($pp_version>=2 && $line =~ /^\#END/) {
	    last;
	}
    }
    close $pp_in;

    return %pairprobinfo;
}

########################################
## minimum of two numbers
##
sub min {
    my ($x,$y) = @_;
    return ($x<=$y)?$x:$y;
}

########################################
## project_string_to_alignment_sequence($str,$seq,$gap_symbols)
##
## Project string to alignment sequence by removing all columns from
## string that have a corresponding gap colunn in the alignment
## sequence
##
## @param $str string
## @param $seq alignment string (sequence with gaps)
## @param $gap_symbols symbols considered as gaps
##
## @pre $str and $seq have same length
##
## @returns projected string
##
########################################
sub project_string_to_alignment_sequence {
    my ($str,$seq,$gap_symbols)=@_;

    die "project_string_to_alignment_sequence: string and alignment sequence must have the same length:\n  $str\n  $seq\n    " if length($str) != length($seq);

    my $res="";
    for (my $i=0; $i<length($seq); $i++) {
	if (substr($seq,$i,1) !~ /[$gap_symbols]/) {
	    $res .= substr $str,$i,1;
	}
    }
    return $res;
}

########################################
## project_structure_to_alignment_sequence($str,$seq,$gap_symbols,$opening_symbols,$closing_symbols,$neutral_symbol)
##
## Project structure string to alignment sequence by removing all
## columns that have a corresponding gap colunn in the alignment
## sequence; in addition we maintain structure: if (only) one end of a base pair is removed,
## the character corresponding to the other base pair is set to $neutral_symbol
##
## @param $str structure string
## @param $seq alignment string (sequence with gaps)
## @param $gap_symbols symbols considered as gaps
## @param $opening_symbols symbols of opening brackets
## @param $closing_symbols symbols of closing brackets
## @param $neutral_symbol neutral symbol
##
## @pre $str and $seq have same length
## @pre opening and closing symbols have same length
##
## @returns projected structure string
##
########################################
sub project_structure_to_alignment_sequence {
    my ($str,
        $seq,
        $gap_symbols,
        $opening_symbols,
        $closing_symbols,
        $neutral_symbol)
      = @_;

    # 0) parse structure
    my $nbs = # number of bracket symbols
      min(length($opening_symbols),
          length($closing_symbols),
         );
    my @str_array = ();
    for (my $i=0; $i<length($str); $i++) {
        push @str_array,-1;
    }
    for (my $i=0; $i<$nbs; $i++) {
        @str_array =
          parse_bracket_structure_single(
                                         $str,
                                         substr($opening_symbols,$i,1),
                                         substr($closing_symbols,$i,1),
                                         \@str_array
                                        );
    }

    # 1) set (potentially kept) ends of deleted base pairs to neutral
    for (my $i=0; $i<length($str); $i++) {
        if ( substr($seq,$i,1) =~ /[$gap_symbols]/ ) {
            if ( $str_array[$i] != -1) {
                substr($str,$str_array[$i],1) = $neutral_symbol;
            }
        }
    }

    # 2) project
    return
      project_string_to_alignment_sequence($str, $seq, $gap_symbols);
}


########################################
## shell_quote(@args)
##
## @param  @args command and argument list
## @return shell-quoted command string
##
########################################

sub shell_quote {
    my @args = @_;
    my $quoted = "";

    foreach my $arg (@args) {
        $arg =~ s/'/'\\''/g;

        $arg = "'$arg'";

        $arg =~ s/^''//;
        $arg =~ s/''$//;

        $quoted .= "$arg ";
    }
    $quoted =~ s/\s+$//;
    return $quoted;
}


########################################
## readpipe_list
##
## safer replacement for readpipe that avoids the shell
##
## @param  @_ command and argument list
## @return reference to list of command output lines
##
########################################
sub readpipe_list {
    my (@cmdlist) = @_;
    my $cmd = shell_quote(@cmdlist);
    $cmd = "$cmd 2>/dev/null";
    # printmsg 1, "$cmd\n";

    open(my $fh, "$cmd |");
    my @result = <$fh>;
    return \@result;
}

########################################
## system_pipein_list
##
## safer replacement for system call with piped in input
## that avoids the shell
##
## @param  @_ command and argument list
##
########################################
sub system_pipein_list {
    my ($input, @cmdlist) = @_;

    my $cmd = shell_quote(@cmdlist);
    $cmd = "echo ".shell_quote($input)." | $cmd >/dev/null";
    printmsg 1, "$cmd\n";

    system($cmd);
}

########################################
## system_pipein_list
##
## safer replacement for system call with piped in input
## that avoids the shell
##
## variant of system_pipein_list using perl 'open'
##
## @param  @_ command and argument list
##
########################################
sub system_pipein_list_open {
    my ($input, @cmd) = @_;

    printmsg 1, "@cmd <<< \"$input\"\n";

    {
        open(my $infh, "|-", @cmd);
        print $infh $input;
        close($infh);
    }
}


########################################
## readpipe_list_open
##
## safer replacement for readpipe that avoids the shell
##
## variant of readpipe_list using perl 'open'
##
## @param  @_ command and argument list
## @return reference to list of command output lines
##
########################################
sub readpipe_list_open {
    open(my $fh, "-|", @_);
    my @result = <$fh>;
    return \@result;
}


########################################
## alifold_structure($filename,@RNAfold_args)
##
## Run RNAalifold on file
##
## $filename string containing name of an clustalw aln file
##
## returns the alifold structure for the given aln file
## as side effect write alirna.ps and aln.ps to current dir
##
########################################
sub alifold_structure {
    my ($filename,@RNAfold_args)=@_;

    my @options = @RNAalifold_options;

    push @options, "--mis";

    if ( not (defined($RNAfold_args[-1]) and $RNAfold_args[-1] eq "--noPS" ) ) {
        push @options, ("--aln", "--color");
    }

    push @options, @RNAfold_args;

    my @aliout = @{
        readpipe_list("$RNAalifold", @options, "$filename")
    };

    if ($#aliout>=1) {
	if ($aliout[1] =~ /([().]+) /) {
	    return $aliout[1];
	}
    }

    return "";
}

########################################
## alifold_mfe($file file name)
##
## compute the alifold mfe for an aln file
##
## returns alifold mfe of the alignment in $file
##
## Does not write ps files.
##
########################################+
sub alifold_mfe {
    my ($file, @RNAfold_args) = @_;

    my @aliout = @{ readpipe_list("$RNAalifold", @RNAalifold_options, "--noPS", @RNAfold_args, "$file") };

    if ($#aliout>=1) {
	if ($aliout[1] =~ /[().]+\s+\(\s*([\d.-]+)\s*=\s*([\d.-]+)\s*\+\s*([\d.-]+)\)/) {
	    my $ali_mfe=$1;
	    my $cons_mfe=$2;
	    my $cov_term=$3;
	    return $ali_mfe;
	}
    }
    return 123456789;
}

########################################
## alifold_pf($file file name)
##
## compute the alifold partition function for an aln file
##
## returns alifold output
##
## side effect: writes files alirna.ps and alidot.ps to current
## directory
##
########################################+
sub alifold_pf {
    my ($file,@RNAfold_args) = @_;

    return
      readpipe_list("$RNAalifold", @RNAalifold_options, "--color", "--mis", @RNAfold_args, "-p", "$file");
}


########################################
## constrain_sequences_from_reliable_structures($seqs,$constraints)
##
## in-place adds (or replaces) structure constraints to a sequences
## data structure, where constraints come from the most reliable
## structure strings.
##
## @param $seqs fasta data structure
## @param $constraints hashs containing constraints (same names as in seqs),
##        all symbols different to '(' or ')' are replaced by '.'
##
########################################
sub constrain_sequences_from_reliable_structures {
    my ($seqs,$constraints) = @_;

    foreach my $i (0..@$seqs-1) {

	if (exists $constraints->{$seqs->[$i]->{name}}) {
	    my $constraint_string = $constraints->{$seqs->[$i]->{name}};
	    $constraint_string =~ s/[^()]/\./g;
	    $seqs->[$i]->{"ANNO#S"} = $constraint_string;
	}
    }
}

########################################
## print_k_dim_hash($h,$k,$s)
##
## Prints a $k dimensional hash of hashs to standard output
##
########################################
sub print_k_dim_hash {
    my ($h,$k,$s) = @_;

    foreach my $i (sort keys %$h) {
	print "$s$i:";
	if ($k==1) {
	    print "$h->{$i}\n";
	} else {
	    print "\n";
	    print_k_dim_hash( $h->{$i},$k-1,"$s\t" );
	}
    }
}




############################################################
## alignment manipulation


########################################
## aln_to_alnloh($seqs,$aln)
##
## Convert hash representation of alignment to list of hash
## representation, such that the resulting alignment representation
## contains an alignment row for each sequence in $seqs and alignment
## rows are in the same order as in $seqs.
##
## The list of hash representation of an alignment is compatible to
## the internal alignment format of RNAz.pm
##
## $seqs ref of list of hash representation of alignment sequences
## $aln ref of alignment hash
##
## The names in $seqs have to occur as keys in $aln
##
## @returns ref of alignment in list of hash representations
########################################
sub aln_to_alnloh {
    my ($seqs,$aln) = @_;

    my @resaln=();

    foreach my $seq (@$seqs) {
       my $id = $seq->{name}; ## the name of sequence as in the alignment

       my %resseq= %{ $seq }; ## copy entries of $seqs
       $resseq{seq} = $aln->{$id};
       push @resaln, \%resseq;
    }

    return \@resaln;
}


########################################
## read_clustalw_alignment($fh)
##
## read multiple alignment in CLUSTALW aln format from filehandle $fh
##
## $fh file handle open for reading
##
## returns tuple of
##     * ref of alignment hash,
##         which  associates names to alignment strings
##     * ref of names list
##         giving information about order (first occurence) of names
##         in the input
##     * flag, whether a clustal header was detected
##
########################################
sub read_clustalw_alignment {
    my ($fh) = @_;

    my %aln;

    my @names=(); ## keep a list of names for order

    my $line;

    my $clustal_header=0;
    if (($line=<$fh>) =~ /^CLUSTAL/) {
	$clustal_header=$line;
	$line=<$fh>;
    }

    do {
        ## regognize lines "name[whitespace]sequence", but allow empty
        ## sequence (even without whitespace after name)!
	if ($line =~ /^(\S+)\s*(\S*)/) {
	    my $name=$1;
	    my $seq=$2;

	    if (!exists $aln{$name}) {
		push @names,$name;
	    }

	    $aln{$name} .= $seq;

	}
    } while ($line = <$fh>);

    return (\%aln,\@names,$clustal_header);
}



########################################
## read_clustalw_alnloh($filename)
##
## read multiple alignment in CLUSTALW aln format from file $filename
##
## @param aln_filename  file name
## @param strict        require header
##
## returns ref of alignment list of hash
##
########################################
sub read_clustalw_alnloh {
    my ($aln_filename, $strict) = @_;

    my $fh;

    open($fh, "$aln_filename")
	|| die "Can not read alignment file $aln_filename\n";
    my ($aln, $names, $header) =
	read_clustalw_alignment($fh);
    close $fh;

    if ($strict and !$header) { die "Missing CLUSTAL header in $aln_filename."; }

    ## generate list of hash alignment
    my @alnloh;

    for my $name (@$names) {
	my $seq = $aln->{$name};
	my $entry = { name=>$name, seq=>$seq};
	push @alnloh,$entry;
    }

    return \@alnloh;
}

########################################
## read_clustalw_aln($filename)
##
## read multiple alignment in CLUSTALW aln format from file $filename
##
## @param aln_filename  file name
##
## the result hash associates names to alignment strings
##
## returns ref of alignment hash and the header line
##
########################################
sub read_clustalw_aln {
    my ($aln_filename) = @_;

    my $fh;

    open($fh, "$aln_filename")
	|| die "Can not read alignment file $aln_filename\n";
    my ($aln, $names, $header) =
	read_clustalw_alignment($fh);
    close $fh;

    if (!$header) { die "Missing CLUSTAL header in $aln_filename."; }

    # the return value was changed form just $aln to $header and $aln
    # in order for stuff that only expects $aln to keep working
    # $aln must be returned last
    return $header, $aln;
}



########################################
## write_clustalw_alnloh($fh, $seqs, $width)
##
## Writes alignment $aln to filehandle $fh, where
## the alignment is given in list of hash format.
## Breaks lines to restrict maximal line width
## Write in CLUSTALW format!
##
## $fh filehandle open for writing
## $seqs loh representation of alignment
## $width line width
##
## @return $namewidth width of sequence names
########################################
sub write_clustalw_alnloh {
    my $fh    = shift;
    my $aln   = shift;
    my $width = shift;
    my $write_header = shift;

    if ($width<=0) {
        $width=undef;
    }

    if (!defined($write_header) || $write_header==1 ) {
	print $fh "CLUSTAL W --- $PACKAGE_STRING\n\n\n"; #  - Local Alignment of RNA
    }

    my $maxlen=0;
    ## determine longest name
    foreach my $seq (@$aln) {
	my $len = length($seq->{name});
	if ($len>$maxlen) {$maxlen=$len;}
    }
    if ($maxlen < 18) {
	$maxlen = 18;
    }

    my $offset=0;
    while(1) {
	my $more=0;

	foreach my $seq (@$aln) {

	    my $s=$seq->{seq};

	    if (defined($width)) {
		if ($offset<length($s)) {
		    $s = substr $s,$offset,$width;
		    $more |= length($seq->{seq})>$offset+$width;
		} else {
		    $s="";
		}
	    }
	    print $fh sprintf("%-".$maxlen."s %s\n",$seq->{name},$s);
	}

	last unless $more;

	print $fh "\n";
	$offset+=$width;
    }
    return $maxlen;
}

########################################
## sprint_fasta_alnloh($fh, $seqs, $width)
##
## Writes alignment $aln to string, where
## the alignment is given in list of hash format.
## Breaks lines to restrict maximal line width.
## Write in FASTA format!
##
## $seqs loh representation of alignment
## $width line width
##
########################################
sub sprint_fasta_alnloh {
    my $aln   = shift;
    my $width = shift;

    my $res="";

    foreach my $seq (@$aln) {

	$res .= ">".$seq->{name};
	if (exists $seq->{desc}) {
	    $res.=" ".$seq->{desc};
	}
	$res.="\n";

	my $offset=0;
	while(1) {
	    my $more=0;

	    my $s=$seq->{seq};
	    if (defined($width)) {
		$s = substr $s,$offset,$width;
		$more |= length($seq->{seq})>$offset+$width;
	    }
	    $res .= "$s\n";

	    last unless $more;

	    $offset+=$width;
	}
    }
    return $res;
}

########################################ä
## loh_names($loh)
##
## returns ref of list of names in list of hash
########################################
sub loh_names {
    my $loh = shift;
    my @names = map { $_->{name}; } @$loh;
    return \@names;
}

########################################ä
## loh_sort($loh,$names)
##
## $loh ref of list of hash alignment or sequences
## $names ref of list of names
##
## returns ref of copy of @$loh sorted by $names
########################################
sub loh_sort {
    my ($loh,$names) = @_;

    my @res;

    my %idx=();
    for my $i (0..(@$names-1)) {
	$idx{$names->[$i]} = $i;
    }

    for my $seq (@$loh) {
        if (exists $idx{$seq->{name}}) {
            $res[$idx{$seq->{name}}] = { %$seq };
        } else {
            push @res, $seq;
        }
    }

    return \@res;
}


# write aln in clustalw format to stdout
sub write_aln {
    my ($aln_ref)=@_;;
    my %aln=%{ $aln_ref };

    foreach my $name (sort(keys %aln)) {
	printf "%-18s %s\n",$name,$aln{$name};
    }
}


########################################
## read_aln_wo_anno($filename)
##
## read multiple alignment in CLUSTALW aln format from file $filename
## ignores annotation starting with #
##
## returns ref of alignment hash
##
########################################
sub read_aln_wo_anno {
    my ($aln_filename) = @_;
    local *ALN_IN;

    open(ALN_IN,$aln_filename) || die "MLocarna::read_aln_wo_anno: Cannot read aln file $aln_filename\n";

    my %aln;

    my %names=();

    my $line;

    if (($line=<ALN_IN>) !~ /CLUSTAL/) {
	printerr "$aln_filename not in clustal aln-format.\n";
	exit(-1);
    }

    while ($line = <ALN_IN>) {
	if ($line =~ /^([^\s]+)\s+(.+)/) {
	    my $name=$1;
	    my $seq=$2;

	    if ($name !~ /^#/) { # ignore annotation
		$aln{$name} .= $seq;
	    }
	}
    }

    close ALN_IN;

    return \%aln;
}


########################################
## project_aln($aln ref of alignment, $names ref of names list)
##
## project multiple alignment to subset of names
## removes only-gap columns
##
sub project_aln {
    my ($aln,$names) = @_;

    my %alnP;

    foreach my $name (@$names) {
	$alnP{$name} = $aln->{$name};
    }

    my $len = aln_length(\%alnP);

    for (my $i=0; $i<$len; ) {
	my $allgap=1;
	foreach my $name (keys %alnP) {
	    my $row = $alnP{$name};
	    # print "$name $row $i\n";
	    if (! is_gap(substr($row,$i,1))) {
		$allgap=0;
		last;
	    }
	}
	if ($allgap) {
	    # remove alignment column i
	    foreach my $name (keys %alnP) {
		substr($alnP{$name},$i,1) = "";
	    }
	    $len--;
	} else {
	    $i++;
	}
    }

    return \%alnP;
}

########################################
## project_alnloh($aln ref of alignment in loh represenation)
##
## project multiple alignment to subset of names
## removes only-gap columns (operate in place!)
##
sub project_alnloh {
    my $aln = shift;

    my $len = length($aln->[0]{seq}); ## number of alignment columns

    for (my $col=0; $col<$len; ) {
	my $allgap=1;
	foreach my $row (0..@$aln-1) {
	    my $rowseq = $aln->[$row]{seq};
	    # print "$name $row $col\n";
	    if (! is_gap(substr($rowseq,$col,1))) {
		$allgap=0;
		last;
	    }
	}
	if ($allgap) {
	    # remove alignment column i
	    foreach my $row (0..@$aln-1) {
		substr($aln->[$row]{seq},$col,1) = "";
	    }
	    $len--;
	} else {
	    $col++;
	}
    }

    return $aln;
}

########################################
## aln_length_atleastonematch($aln ref of alignment)
##
## returns number of alignment columns with at least one match
##
sub aln_length_atleastonematch {
    my ($aln) = @_;

    my $len=aln_length($aln); ## length alignment (and of each alignment string)
    my $mylen=0; ## counts columns with at least one match

    for (my $i=0; $i<$len; $i++) {
	my $non_gap=0;
	foreach my $k (keys %$aln) {
	    if (substr($aln->{$k},$i,1) ne "-") {$non_gap++;}
	}
	if ($non_gap>1) {$mylen++;}
    }

    return $mylen;
}


########################################
## read_2D_matrix($file string, $n uint, $m uint)
##
## read a $nx$m-matrix from file, matrix entries are
## separated by white space, rows are separated by new line
##
########################################
sub read_2D_matrix {
    my ($file,$n,$m) = @_;

    my @matrix=();

    my $MAT;

    open($MAT,"$file") || die "Cannot read matrix $file\n";

    for (my $a=0; $a<$n; $a++) {
	my $line;
	if ($line=<$MAT>) {

	    chomp $line;
	    $line=~s/^\s+//;

	    my @row=split /\s+/,$line;
	    if (int(@row) != $m) {
		die "Expect $m entries per row while reading matrix from $file; found ".int(@row)."\n";
	    }

	    push @matrix, [ @row ];
	}
	else {
	    die "Expect $n rows while reading matrix from $file, found ".($a-1)."\n";
	}
    }
    close $MAT;

    return \@matrix;
}

########################################
## write_2D_matrix($file string, $matrix ref of 2D array)
##
## write a matrix given by 2D-array to file, matrix entries are
## assumed to be integers and output is formatted such that integers
## can have up to 6 digits.
##
########################################
sub write_2D_matrix {
    my ($file,$matrix) = @_;

    open(MAT,">$file") || die "Cannot write matrix $file\n";

    my $size_x=@$matrix;
    for (my $a=0; $a<$size_x; $a++) {
	my @row = @{ $matrix->[$a] };
	my $size_y = $#row+1;
	for (my $b=0; $b<$size_y; $b++) {
	    if ($a==$b) {print MAT "     0 ";} else {
		printf MAT "%6d ",$matrix->[$a][$b];
	    }
	}
	print MAT "\n";
    }
    close MAT;
}

############################################################
## conversion between sequence and alignment positions
##

########################################
sub seqpos_to_alipos {
    my ($ali) = @_;

    my @res;

    my $seqpos=0;
    my $alipos=0;

    for my $c (split //, $ali) {
	if ($c =~ /[ACGUT]/) {
	    push @res,$alipos;
	    $seqpos++;
	}
	$alipos++;
    }

    return @res;
}


########################################
sub alipos_to_seqpos {
    my ($ali) = @_;

    my @res;

    my $seqpos=0;
    my $alipos=0;

    for my $c (split //, $ali) {
	if ($c =~ /[ACGUT]/) {
	    push @res,$seqpos;
	    $seqpos++;
	} else {
	    push @res,-1;
	}
	$alipos++;
    }

    return @res;
}


########################################
## write pp file from aln and consensus dp, write only probs >= min_prob
sub write_pp {
    my ($file,$aln_ref, $cons_ref, $min_prob) = @_;
    my %aln = %{ $aln_ref };
    my @cons = @{ $cons_ref };

    local *OUT;
    open(OUT, ">$file") || die "Cannot open $file for writing\n";

    print OUT "SCORE: 0\n\n";

    foreach my $name (keys %aln) {
	printf OUT "%-18s %s\n", $name, $aln{$name};
    }


    print OUT "\n#\n";

    my $len = aln_length(\%aln);

    for (my $i=1; $i<=$len; $i++) {
	for (my $j=$i+1; $j<=$len; $j++) {
	    if (defined $cons[$i][$j]) {
		my $p = $cons[$i][$j];
		if ($p >= $min_prob) {
		    print OUT "$i $j $p\n";
		}
	    }
	}
    }

    close OUT;
}


########################################
## write_tcoffee_lib_file($filename,$alignments)
##
## takes a filename and an array ref of alignemtns and
## outputs a tcoffe lib file
##
## arg $filename      name of the library file (is created)
## arg \@alignments    ref to array of hash refs where key is the name
##                     of the sequence and value the alignment row
##                     of sequences given by names
##
## Lib format definitions can be found here:
## http://www.tcoffee.org/Documentation/t_coffee/t_coffee_technical.htm#_Toc256781778
########################################
sub write_tcoffee_lib_file {
    my ($filename, $alignments) = @_;

    # open the file
    open(my $FILE, ">$filename");

    # first collect all the sequences
    my %sequences;
    foreach my $aln (@{$alignments}) {
        while( my ($name, $seq) = each(%{$aln->{rows}})){
            if (!defined $sequences{$name}) {
                my $temp = $seq;
                # remove all non alphanumeric characters (i.e. Gap symbols)
                $temp =~ s/\W//g;
                $sequences{$name} = { "sequence" => $temp,
                                      "number"   => undef};
            }
        }
    }

    ## assign sequence numbers in alphanumeric order of sequence names
    my $seqCounter = 1;
    for my $name (sort keys(%sequences)) {
        $sequences{$name}->{"number"} = $seqCounter;
        $seqCounter++;
    }

    # write the header to the file
    print $FILE "! TC_LIB_FORMAT_01\n";
    # print number of sequences
    print $FILE "" . ($seqCounter - 1) . "\n";

    # print sequences
    foreach my $key (sort{ $sequences{$a}->{number} <=> $sequences{$b}->{number}} (keys(%sequences))) {
        print $FILE $key . " " . length($sequences{$key}->{sequence}) . " "
                    . $sequences{$key}->{sequence} . "\n";
    }


    foreach my $aln (@{$alignments}) {
        # get the sequence indexes
        my @tmp = keys(%{$aln->{rows}});
        my $seqAname = $tmp[0];
        my $seqA =  $aln->{rows}->{$tmp[0]};
        my $seqBname = $tmp[1];
        my $seqB =  $aln->{rows}->{$tmp[1]};
        print $FILE "#" . $sequences{$seqAname}->{number} . " " . $sequences{$seqBname}->{number} . "\n";
        # get the alignment edges
        my $edges = get_alignment_edges($seqA, $seqB);
        for(my $i = 0; $i < scalar(@{$edges}); $i++) {
            if ($edges->[$i] != -1) {
                print $FILE ($i + 1) . " " . ($edges->[$i] + 1) . " " . $aln->{score} . "\n";
            }
        }
    }
    print $FILE "! SEQ_1_TO_N\n";
    close($FILE);
}

########################################
## get_alignment_edges($sequenceA,$sequenceB)
##
## takes two rows of an alignment and returns an array ref
## with an element for each position in the first sequence
## where edges[$i] == -1 if the position is unmatched and
##       edges[$i] == $j if the postion is matched to base $j of the other
##                       sequence
##
## arg $sequenceA    first alignment row
## arg $sequenceB    second alignment row
##
########################################
sub get_alignment_edges($$)
{
  my ($sequenceA, $sequenceB) = @_;

  my @gappedA = split //, $sequenceA;
  my @gappedB = split //, $sequenceB;
  my $gapFree = $sequenceA;
  $gapFree =~ s/\W//g;
  my @seq = split //, $gapFree;

  # a counter pair, pointing to the current position in both sequences
  my @currentPos = (0 , 0);

  for(my $i = 0; $i < scalar(@gappedA);$i++)
  {
    # case 1: match / mismatch
    if ($gappedA[$i] ne "-" and $gappedB[$i] ne "-") {
        # set the alignment edge
        $seq[$currentPos[0]] = $currentPos[1];

        # count up the pointers
        $currentPos[0]++;
        $currentPos[1]++;
        next;
    }
    # case 2: gap in first sequence
    if ($gappedA[$i] eq "-" and $gappedB[$i] ne "-") {
        # count up second counter
        $currentPos[1]++;
        next;
    }
    # case 3: gap in second sequence
    if ($gappedA[$i] ne "-" and $gappedB[$i] eq "-") {
        $seq[$currentPos[0]] = "-1";

        # count up first counter
        $currentPos[0]++;
        next;
    }
    # this should never be reached as that would mean - and - are matched
    print "two gaps matched!!!!!";
  }

  return \@seq;
}


########################################
## parse pseudoknotted dot-bracket structures
## and generate pp file for fixed structure
##
sub parse_bracket_structure_single {
    my ($str,$open,$close,$str_array_ref)=@_;

    my @str_array;

    if (defined($str_array_ref)) {
        @str_array = @{ $str_array_ref };
    } else {
        @str_array = (-1)x(length($str));
    }

    my @stack;

    for (my $i=0; $i<length($str); $i++) {
	my $c=substr $str,$i,1;

	if ($c eq $open) {
	    push @stack,$i;
	} elsif ($c eq $close) {
	    my $j=pop @stack;
	    $str_array[$i]=$j;
	    $str_array[$j]=$i;
	}
    }

    return @str_array;
}

########################################
## @returns array for structure
## entry per sequence position: gives paired position or -1
sub parse_bracket_structure {
    my ($str)=@_;
    my @str_array;

    @str_array=parse_bracket_structure_single $str,"(",")";
    @str_array=parse_bracket_structure_single $str,"{","}",\@str_array;
    @str_array=parse_bracket_structure_single $str,"[","]",\@str_array;
    @str_array=parse_bracket_structure_single $str,"<",">",\@str_array;
    @str_array=parse_bracket_structure_single $str,"A","a",\@str_array;
    @str_array=parse_bracket_structure_single $str,"B","b",\@str_array;
    @str_array=parse_bracket_structure_single $str,"C","c",\@str_array;
    @str_array=parse_bracket_structure_single $str,"D","d",\@str_array;

    return @str_array;
}


########################################
## convert_fix_structure_to_pp($ppfilename,$name,$seq,$str,$constraints)
##
## @param $ppfilename output file name
## @param $name sequence name
## @param $seq sequence string
## @param $str structure string
## @param $constraints anchor constraint string
##
sub convert_fix_structure_to_pp {
    my ($ppfilename,$name,$seq,$str,$constraints) = @_;

    my @str=parse_bracket_structure($str);

    local *OUT;

    open(OUT,">$ppfilename") || die "Cannot write $!";

    print OUT "#PP 2.0\n\n";
    print OUT "$name\t$seq\n";

    if (defined $constraints && $constraints ne "") {
    	my @cs = split /\#/, $constraints;
	for (my $i=0; $i<@cs;$i++) {
	    print OUT "#A".($i+1)." ".$cs[$i]."\n";
	}
    }
    print OUT "\n#END\n\n";
    print OUT "#SECTION BASEPAIRS\n\n";

    for (my $i=0; $i<=$#str; $i++) {
	if ($str[$i]>$i) {
	    print OUT ($i+1)." ".($str[$i]+1)." 1.0\n";
	}
    }
    print OUT "\n#END\n";
}


########################################
## find_in_exec_path
##
## search in PATH, unless absolute or relative path is given,
## and return the absolute path
##
## @param prg name or path of executable
##
## @return absolute path of the executable or undef if not found
##
sub find_in_exec_path {
    my $prg = shift;
    my $whichprg = readpipe("which $prg");
    chomp $whichprg;
    if ($? == 0) {
        return Cwd::abs_path($whichprg);
    } else {
        return undef;
    }
}

########################################
## find_in_exec_path_or_error
##
## search in PATH, unless absolute or relative path is given,
## and return the absolute path. On error, exit with error message.
##
## @param prg name or path of executable
##
## @return absolute path of the executable; exit if no success
##
sub find_in_exec_path_or_error {
    my $prg = shift;
    my $res=find_in_exec_path($prg);
    if (! defined($res)) {
        printerr "ERROR: cannot find required executable $prg.\n";
        exit -1;
    }
    return $res;
}



## ------------------------------------------------------------


1;
