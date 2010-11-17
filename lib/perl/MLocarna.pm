package MLocarna;

use 5.008003;
use strict;
use warnings;


#use FindBin;
#use lib $FindBin::RealBin;
use MLocarna::SparseMatrix;
use MLocarna::Trees;
use MLocarna::MatchProbs;
use MLocarna::Aux;

require Exporter;
    
# set the version for version checking
our $VERSION     = 1.00;

our @ISA         = qw(Exporter);
our @EXPORT      = qw(

read_fasta
check_constraints_in_fasta

sequence_constraint_string

compute_alignment_from_seqs
compute_alignment_from_dps
compute_alignment_score
read_dp_ps
read_pp_file_aln
read_pp_file_pairprobs
convert_dp_to_pp_with_constraints
compute_alifold_structure
constrain_sequences


extract_from_clustal_alignment

project_str_to_aln_seq

clone_hash

read_aln
write_aln
project_aln
aln_h2loh
write_clustalw_loh
aln_length_atleastonematch

aln_to_alnloh
read_clustalw_alignment
read_clustalw_alnloh
read_clustalw_aln
write_clustalw_alnloh
loh_translate_names
loh_associate_nnames_to_names
loh_names 
loh_sort


write_pp

write_2D_matrix

new_intermediate_name

extract_score_matrix_from_alignments

aln_reliability_beststruct_fromfile
aln_reliability_fromfile

);

our %EXPORT_TAGS = ();

# your exported package globals go here,
# as well as any optionally exported functions
our @EXPORT_OK   = qw(
$RNAalifold
$PACKAGE_STRING

parse_mfasta
parse_mfasta_constraints
);

our $RNAalifold_params="-r";

our $RNAalifold = "RNAalifold";

our $PACKAGE_STRING = "MLocarna";


########################################
## new_intermediate_name
##
## returns a new name for intermediate files
## uses counter 
my $intermediate_name_counter=0;
## and basename
my $intermediate_name_base="intermediate";
##
## These names are used by mlocarna as filenames/identifiers of
## intermediary alignments during the computation of progressive and
## iterative alignments
##
########################################
sub new_intermediate_name {
    $intermediate_name_counter++;
    return "$intermediate_name_base$intermediate_name_counter";
}


########################################
## threading support
##

# use Config;

# my $THREADS_ENABLED = $Config{useithreads};

# if ($THREADS_ENABLED) {
#     require MLocarna::threaded;
#     import MLocarna::threaded;
# } else {
#     require MLocarna::unthreaded;
#     import MLocarna::unthreaded;
# }



########################################
## make_normalized_seqname($name, @names list of existing names)
##
## generate a sequence name from $name that
## has at most a length of 16
## and
## does not already exist in @names
## 
########################################
sub make_normalized_seqname {
    my ($name, @names) = @_;
    
    my $maxlen=16;

    chomp $name;
    $name =~ s/[^a-zA-Z\d]/_/g;
    $name = substr $name,0,$maxlen;
    
    my $i=1;
    while (grep /^$name$/, @names) {
	my $arity = int(log($i)/log(10))+1;
	if ($arity >= 4) { # arbitrary limit
	    die "Could not generate unique name";
	}
	
	$name = substr $name,0,$maxlen-$arity-1;
	$name = sprintf("%s_%0$arity"."d",$name,$i);
	$i++;
    }
    return $name;
}



########################################
## read_fasta($filename)
## ------------------------------
##
## Read a simplified fasta file. If filename ends with "*.gz" the file
## is automatically uncompressed by gunzip.  Dies if file is not
## readable or gunzippable.  The result is returned as list of the
## sequences.
##
##
## Each sequence is encoded as hash of features
##   name: id of sequence
##   descr: description of the sequence
##   nname: normalized id of the sequence 
##       (generated from id, 
##        can be used as a filename which is unique for the set of
##        sequences in the fasta file)
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
    my @names = ();
    
    my $line=<$fh>;
    while(defined($line)) {
	if ($line=~/^>\s*(\S+)\s*(.*)/) {
	    my $name=$1;
	    my $description=$2;
	    
	    my $nname = make_normalized_seqname $name,@names;
	    
	    push @names, $nname;
	    
	    my $seq = { name  => $name,
			nname => $nname, 
			descr => $description };
	    
	    while (defined($line=<$fh>) && ($line !~ /^>/)) {
		chomp $line;
		$line =~ s/\s+//g;
		
		if  ($line =~ /(.+)\#(.+)/) {
		    $seq->{"ANNO\#$2"} .= $1;
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
## check_constraints_in_fasta($fasta)
##
## Checks validity of constraint description in fasta list of hash
## representation.
##
## $fasta list of hashs representation
##
## die with error message when constraint annotation is invalid
##
########################################
sub check_constraints_in_fasta {
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
	if ($numseqconstraints == -1) {$numseqconstraints=$i-1;}
	
	if ($numseqconstraints != $i-1) {
	    die "Bad sequence constraints for sequence "
		.$seq->{name}.".";
	}
    }
}

########################################
## sequence_constraint_string($seq hash ref)
##
## Collect sequence constraint string for $seq entry of a fasta list of hashs
##
## returns string "<c1>#<c2>#...<cN>", where ci is the constraints description line i
##
########################################+
sub sequence_constraint_string {
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
sub parse_mfasta($) {
    my ($file) = @_;
    my %mfasta;
    
    local *PMF_IN;

    print STDERR "Use of parse_mfasta is deprecated. Use read_fasta instead.";
    
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
		print STDERR "Duplicate name in mfasta: $name\n";
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

    print STDERR "Use of parse_mfasta_constraints is deprecated. Use read_fasta instead.";

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
		print STDERR "Duplicate name in mfasta, cannot make unique : $name\n";
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
		print STDERR "Structure constraint length unequals sequence length for $name.\n";
		exit -1;
	    }

	    $mfasta{"$name\#S"} = $strcons{$name};
	}
	
	my $i=1;
	while (exists $seqcons{"$name\#$i"}) {
	    if (length $seqcons{"$name\#$i"} != length $mfasta{$name}) {
		print STDERR "Sequence constraint length unequals sequence length for $name.\n";
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
	    print STDERR "Bad sequence constraints for sequence $name.\n";
	    exit -1;
	}
    }
    
    for my $name (@names) {
	if (($name!~/#/) &&  $has_seq_constraints) {
	    if (! exists $mfasta{"$name\#C"}) {
		print STDERR "Sequence constraint string missing for $name.\n";
	    }
	}
    }
    #
    ## ----------------------------------------
    
    return %mfasta;
}

## convert a dotplot file to a pp file
## thereby insert the sequence constraint strings
sub convert_dp_to_pp_with_constraints($$$$$$) {
    my ($dpfile,$ppfile,$name,$sequence,$constraints,$read_condprobs) = @_;
    
    open(PP_OUT,">$ppfile") || die "Cannot open $ppfile for writing."; 
    open(DP_IN,"$dpfile") || die "Cannot open $dpfile for reading."; 
    
    print PP_OUT "SCORE: 0\n\n";
    print PP_OUT "$name $sequence\n";
    if (defined $constraints && $constraints ne "") {
	print PP_OUT "#C $constraints\n";
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
    
    print PP_OUT "\n\#\n";
    
    for my $i ( keys %pp ) {
	for my $j ( keys %{ $pp{$i} } ) {
	    print PP_OUT "$i $j $pp{$i}{$j}\n";
	}
    }
    
    close PP_OUT;
}


## compute pairwise alignment for given sequences using RNAfold -p for generating dps
sub compute_alignment($$$$$) {
    my ($bindir,$seqA,$seqB,$locarna_params,$tmpprefix) =  @_;
    
    local *CA_IN;
    
    my $tmpfile="$tmpprefix.clustal";

    system("$bindir/locarna.pl $seqA $seqB $locarna_params --clustal=$tmpfile >/dev/null");
    
    open(CA_IN,"$tmpfile") || die "Cannot read tmp file $tmpfile\n";
    
    my @content=<CA_IN>;
   
    close CA_IN;
    
    unlink $tmpfile;
    
    return @content;
}


## compute the alignment from dp-files or pp-files
sub compute_alignment_from_dps($$$$$) {
    my ($bindir,$seqA_dp,$seqB_dp,$locarna_params,$tmpprefix) =  @_;

    local *CA_IN;

    my $tmpfile="$tmpprefix.clustal";

    my $cmd="$bindir/locarna $seqA_dp $seqB_dp $locarna_params --clustal=$tmpfile >/dev/null";
    system $cmd || die "Command $cmd failed."; 
    
    open(CA_IN,"$tmpfile") || die "compute_alignment: Cannot read temporary file $tmpfile\n";
    
    my @content=<CA_IN>;
   
    close CA_IN;
    
    unlink $tmpfile;
    
    return @content;
}


sub compute_alignment_score($$$$$) {
    my ($bindir,$seqA,$seqB,$locarna_params,$tmpprefix) =  @_;
    
    my @content=compute_alignment($bindir,$seqA,$seqB,$locarna_params,$tmpprefix);
    # print "@content\n";
    $content[0] =~ /Score: ([\d\.\-]+)/ || die "Cannot return score.\n";
    my $score=$1;

    return $score;
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



sub read_pp_file_aln($) {
    my ($filename)=@_;
    local *PP_IN;
   
    open(PP_IN,$filename) || die "Can not read $filename\n";
    
    my %aln;
    my %pairprobs;

    my $line;
    
    while ($line = <PP_IN>) {
	if ($line =~ /^\#\s*$/) { last; }
    
	if (($line =~ /^(\S+)\s+(.+)/) && ($line !~ /^SCORE:/) ) {
	    $aln{$1}.=$2;
	}
    }
    
    close PP_IN;
    
    return %aln;
}

## read a pp file and return a hash of pair probabilities
##
## @param $filename name of the pp file
## @returns hash of pair probabilities reported in the file
##          the hash has keys "$i $j" and contains the probabilities p_ij (i<j)
##          positions in the hash are in [1..sequence length]
##
sub read_pp_file_pairprobs($) {
    my ($filename)=@_;
    local *PP_IN;

    open(PP_IN,$filename) || die "Can not read $filename\n";
    
    #my %aln;
    my %pairprobs;

    my $line;

    while ($line = <PP_IN>) {
	if ($line =~ /^\#\s*$/) { last; }

	if ($line =~ /^(\S+)\s+(.+)/) {
	    #$aln{$1}.=$2;
	}
    }
    
    while (($line = <PP_IN>)) {
	if ($line =~ /(\S+)\s+(\S+)\s+(\S+)/) {
	    my $i=$1;
	    my $j=$2;
	    my $p=$3;
	    $pairprobs{"$i $j"}=$p;
	}
    }
    close PP_IN;
    
    return %pairprobs;
}


########################################
## project_str_to_aln_seq($str,$seq)
##
## Project structure string to alignment sequence
## Remove all columns from structure string that
## have a corresponding gap colunn in alignment sequence
##
## $str structure string (with gaps)
## $seq sequence with gaps
##
## pre: $str and $seq must have same length
## consider everything different from [a-zA-Z] as gap
##
## returns resulting structure string
##
########################################
sub project_str_to_aln_seq($$) {
    my ($str,$seq)=@_;
    my $res="";
    for (my $i=0; $i<length($seq); $i++) {
	if (substr($seq,$i,1) =~ /[a-zA-Z]/) {
	    $res .= substr $str,$i,1;
	}
    }
    return $res;
}



########################################
## clone_hash($levels,$x)
##
## make a non-shared copy of a hash
##
## (one could use ref as "typeof" to make this more generic,
## however its not needed for the moment)
########################################
sub clone_hash {
    my ($levels,$x) = @_;
    my %y;
	
    if ($levels>1) {
	foreach my $k (keys %{ $x }) {
	    $y{$k} = clone_hash( $levels-1, $x->{$k} );
	}
    } else {
	%y = %{ $x };
    }
    return \%y;
}



########################################
## compute_alifold_structure($filename)
##
## Run RNAalifold on file
##
## $filename string containing name of an clustalw aln file
##
## returns the alifold structure for the given aln file
## as side effect write alirna.ps to current dir
##
########################################
sub compute_alifold_structure {
    my ($filename)=@_;
    
    my @aliout =  readpipe("$RNAalifold $RNAalifold_params $filename 2>/dev/null");
    
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
########################################+
sub alifold_mfe {
    my ($file) = @_;
    
    my @aliout =  readpipe("RNAalifold $RNAalifold_params $file 2>/dev/null");

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
## constrain_sequences($seqs,$constraints)
##
## in-place adds (or replaces) structure constraints to a sequences 
## data structure
## 
## @param $seqs fasta data structure
## @param $constraints hashs containing constraints (same names as in seqs),
##        all symbols different to '(' or ')' are replaced by '.'
##
########################################
sub constrain_sequences($$) {
    my ($seqs,$constraints) = @_;
        
    foreach my $i (0..@$seqs-1) {
	
	if (exists $constraints->{$seqs->[$i]->{nname}}) {
	    my $constraint_string = $constraints->{$seqs->[$i]->{nname}}->{"ANNO\#S"};
	    $constraint_string =~ s/[^()]/\./g;	    
	    $seqs->[$i]->{"ANNO\#S"} = $constraint_string;
	    printmsg 3,"ANNO\#S: $constraint_string\n";
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
## rows are in the same order as in $seqs. The normalized names (nname)
## of $seq have to occur as keys in the alignment hash.
##
## The list of hash representation of an alignment is compatible to
## the internal alignment format of RNAz.pm
##
## $seqs ref of list of hash representation of alignment sequences
## $aln ref of alignment hash
##
## @returns ref of alignment in list of hash representations
########################################
sub aln_to_alnloh($$) {
    my ($seqs,$aln) = @_;
    
    my @resaln=();
    
    foreach my $seq (@$seqs) {
       my $id = $seq->{nname}; ## the name of sequence as in the alignment 
       
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
##         which  associates normalized names to 
##         alignment strings
##     * ref of names list
##         giving information about order (first occurence) of names
##         in the input
##     * ref of hash name=>nname
##     * flag, whether a clustal header was detected
##
########################################
sub read_clustalw_alignment {
    my ($fh) = @_;
    
    my %aln;
    
    my %nnames=(); ## keep a hash as association name=>normalized name
    my @names=(); ## keep a list of non-normalized names for order
    
    my $line;
    
    my $clustal_header=0;
    if (($line=<$fh>) =~ /^CLUSTAL/) {
	$clustal_header=1;
	$line=<$fh>;
    }
    
    do {
	if ($line =~ /^([^\s]+)\s+(.+)/) {
	    my $name=$1;
	    my $seq=$2;
	    
	    my $nname;
	    if (exists $nnames{$name}) {
		$nname = $nnames{$name};
	    } else {
		$nname = make_normalized_seqname $name,values %nnames;
		$nnames{$name} = $nname;
		push @names, $nname;
	    }
	    
	    $aln{$nname} .= $seq;
	}
    } while ($line = <$fh>);
    
    return (\%aln,\@names,\%nnames,$clustal_header);
}



########################################
## read_clustalw_alnloh($filename)
##
## read multiple alignment in CLUSTALW aln format from file $filename
##
## the names of sequences are converted to normalized names
## NOTE that normalization of names depends on the order of names
## in the alignment if some names are not unique in the first 16 characters
##
## returns ref of alignment list of hash
##
########################################
sub read_clustalw_alnloh {
    my ($aln_filename) = @_;
		
    my $fh;
    
    open($fh, "$aln_filename") 
	|| die "Can not read alignment file $aln_filename\n";
    my ($aln, $names, $nnames, $header) = 
	read_clustalw_alignment($fh);
    close $fh;
    
    if (!$header) { die "Missing CLUSTAL header in $aln_filename."; }
    
    ## generate list of hash alignment
    my @alnloh;
    
    for my $name (@$names) {
	my $nname = $nnames->{$name};
	my $seq = $aln->{$nname};
	my $entry = { name=>$name, nname=>$nname, seq=>$seq};
	push @alnloh,$entry;
    }
    
    return \@alnloh;
}

########################################
## read_clustalw_aln($filename)
##
## read multiple alignment in CLUSTALW aln format from file $filename
##
## the names of sequences are converted to normalized names
## NOTE that normalization of names depends on the order of names
## in the alignment if some names are not unique in the first 16 characters
##
## the result hash associates normalized names to alignment strings
##
## returns ref of alignment list of hash
##
########################################
sub read_clustalw_aln {
    my ($aln_filename) = @_;
		
    my $fh;
    
    open($fh, "$aln_filename") 
	|| die "Can not read alignment file $aln_filename\n";
    my ($aln, $names, $nnames, $header) = 
	read_clustalw_alignment($fh);
    close $fh;
    
    if (!$header) { die "Missing CLUSTAL header in $aln_filename."; }
		
    return $aln;    
}



########################################
## write_clustalw_alnloh($fh, $seqs, $width)
## 
## 
## Writes alignment $aln to filehandle $fh, where
## the alignment is given in list of hash format.
## Breaks lines to restrict maximal line width
##
## $fh filehandle open for writing
## $seqs loh representation of alignment
## $width line width
##
## Writes non-normalized names
##
########################################
sub write_clustalw_alnloh {
    my $fh    = shift;
    my $aln   = shift;
    my $width = shift;

    print $fh "CLUSTAL W --- $PACKAGE_STRING - Local Alignment of RNA\n\n\n";
    
    my $maxlen=0;
    ## determine longest name
    foreach my $seq (@$aln) {
	my $len = length($seq->{name});
	if ($len>$maxlen) {$maxlen=$len;}
    }
    $maxlen+=1;
    if ($maxlen < 18) {
	$maxlen = 18;
    }
    
    my $offset=0;
    while(1) {	
	my $more=0;

	foreach my $seq (@$aln) {
	    
	    my $s=$seq->{seq};
	    if (defined($width)) {
		$s = substr $s,$offset,$width;
		$more |= length($seq->{seq})>$offset+$width;
	    }
	    print $fh sprintf("%-".$maxlen."s %s\n",$seq->{name},$s);
	}
	
	last unless $more;

	print $fh "\n";
	$offset+=$width;
    }
}

########################################
## loh_translate_names($alnloh,$names)
##
## Translate names in $alnloh according to hash $names
##
## $alnloh alignment as list of hashs
## $names  hash of names
##
## Example: 
##    use 
##       alnloh_translate_names($alnloh,loh_associate_nnames_to_names($seqs))
##    to 'unnormalize' names in $alnloh
##
## returns ref to a copy of $alnloh where names are replaced by their
## associated values in hash $names.
##
########################################
sub loh_translate_names {
    my ($alnloh,$names) = @_;
    
    my @res=();

    for my $seq (@$alnloh) {
	my %entry = %$seq;
	$entry{name}=$names->{$seq->{name}};
	push @res,\%entry;
    }
    
    return \@res;
}

########################################
## loh_associate_nnames_to_names($loh)
##
## Associate normalized names to their long names as in $loh.
##
## $loh sequences or alignment as list of hashs
## 
## returns ref to a hash that associates the names.
##
########################################
sub loh_associate_nnames_to_names {
    my $loh = shift;
    
    my %names=();
    for my $seq (@$loh) {
	$names{$seq->{nname}}=$seq->{name};
    }
    
    return \%names;
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
## $loh ref of list of hash alignemnt or sequences
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
	$res[$idx{$seq->{name}}] = { %$seq };
    }
    
    return \@res;
}


# write aln in clustalw format to stdout
sub write_aln($) {
    my ($aln_ref)=@_;;
    my %aln=%{ $aln_ref };
    
    foreach my $name (sort(keys %aln)) {
	printf "%-18s %s\n",$name,$aln{$name};
    }
}

########################################
## read_aln($filename)
##
## read multiple alignment in CLUSTALW aln format from file $filename
##
## the names of sequences are converted to normalized names,
## potential additional information in long names is lost
##
## NOTE that normalization of names depends on the order of names
## in the alignment if some names are not unique in the first 16 characters
##
## returns ref of alignment hash
##
########################################
sub read_aln($) {
    my ($aln_filename) = @_;
    local *ALN_IN;
    
    open(ALN_IN,$aln_filename) || die "Can not read aln file $aln_filename\n";
    
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
	    
	    my $nname;
	    if (exists $names{$name}) {
		$nname = $names{$name};
	    } else {
		$nname = make_normalized_seqname $name,values %names;
		$names{$name} = $nname;
	    }
	    	    
	    $aln{$nname} .= $seq;
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
    
    return %alnP;
}

########################################
## aln_length_atleastonematch($aln ref of alignment)
##
## returns number of alignment columns with at least one match
##
sub aln_length_atleastonematch($) {
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
sub write_pp($$$$) {
    my ($file,$aln_ref, $cons_ref, $min_prob) = @_;
    my %aln = %{ $aln_ref };
    my @cons = @{ $cons_ref };
    
    local *OUT;
    open(OUT, ">$file") || die "Cannot open $file for writing\n";

    print OUT "SCORE: 0\n\n";
    
    foreach my $name (keys %aln) {
	printf OUT "%-20s %s\n", $name, $aln{$name};
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
## extract_score_matrix_from_alignments($names,$pairwise_aln)
##
## compute the score matrix from all pairwise alinments
##
## arg \@names         ref to list of sequence names
## arg \@pairwise_aln  ref to 2D-array of all pairwise alignments
##                     of sequences given by names
##
## returns ref to 2D-array of scores (symmetric),
##         indices are positions in name string
##
########################################
sub extract_score_matrix_from_alignments($$) {
    my ($names,$pairwise_alns_ref) = @_;
    my @pairwise_alns = @{ $pairwise_alns_ref };
    
    my @score_matrix; ## result
    

    for (my $a=0; $a<@$names; $a++) {
	for (my $b=0; $b<$a; $b++) {
	    
	    my @aln = @{ $pairwise_alns[$a][$b] };
		
	    $aln[0] =~ /Score: ([\d\.\-]+)/ || die "Cannot extract score from $aln[0].\n";
	    my $score=$1;
	    
	    $score_matrix[$a][$b] = $score;
	    $score_matrix[$b][$a] = $score;
	}
    }

    return \@score_matrix;
}


########################################
## compute reliability of an alignment as reliability for the best structure
sub aln_reliability_beststruct_fromfile($$$) {
    my ($file,$bmprobs,$amprobs)=@_;
    my $aln=read_aln($file);
    
    my ($score,$rel_str) = evaluate_alignment($aln,$bmprobs,$amprobs,0);
    
    return $score;
}

########################################
## aln_reliability_fromfile($file,$bmprobs,$amprobs)
##
## compute reliability of an alignment as average of the column reliability
##
sub aln_reliability_fromfile($$$) {
    my ($file,$bmprobs,$amprobs)=@_;
    my $aln=read_aln($file);
    
    return aln_reliability($aln,$bmprobs,$amprobs);
}


## ------------------------------------------------------------


1;