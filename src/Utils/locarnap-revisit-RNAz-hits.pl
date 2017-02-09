#!/usr/bin/env perl

=head1 NAME

locarnap-revisit-RNAz-hits

=head1 SYNOPSIS

locarnap-revisit-RNAz-hits [options]

=head1 DESCRIPTION

Revisits the hits of the fly RNAz screen and prepares input data for
realigning all loci with locarnap-realign-all.pl. The script is
usually used as first step in a pipeline for refining RNAz hits with
LocARNA-P.

=head1 OPTIONS

=over 4

=item  B<--help>

Brief help message

=item  B<--man>

Full documentation

=item  B<--all>

Extract all loci (in contrast, the default mode selects only loci that
overlap with Flybase ncRNA annotation).

=item  B<--random <n>>

Special mode: draw n random instances

=item B<--context <c>>

Extract with maximal context of c columns upstream and downstream.

=back

In default mode, determine the RNAz hits with dmel flybase
annotation. In random mode, draw n random hits. With flag --all select
all hits.  For the selected hits get the annotation and the position
in the pecan alignment. Determine the sequences that are well-aligned
to the dmel sequence, obtain these sequences with genomic context
(option --context, default=100).

Goal: use these sequence for realining by locarna. From the
locarna reliability profile, determine boundaries of the
ncRNA. Compare to the flybase annotation and RNAz boundaries.

=cut


use warnings;
use strict;
use FindBin;
use lib $FindBin::Bin;

use LocARNA_RNAz;

## bioperl packages

use Bio::Seq;
use Bio::SeqIO;

use Bio::AlignIO;

use Bio::LocatableSeq;


### ----------------------------------------
## constants
# that control the study
#


## only look at flybase annotation that is not longer than max_len
my $max_len=400;

## consider that  much context to the left and right
my $context_size=100;

## maximal tolerated content of gaps in the sequence slices
my $max_gap_content=0.25;

## assume that the script is started from a "home-directory" that contains
## all the inut-data in subdirectory "Data"
#
my $homedir = readpipe("pwd");
chomp $homedir;

my $datadir = "$homedir/Data";

# where to put fasta files with sequences for realignment
my $mfa_outdir = "$homedir/Realign-Sequences";

# where to put slices of the alignment (for control)
my $aln_outdir = "$homedir/Alignment-Slices";


## perform study for the chromosomes in the list
#
#my @chromosomes=("2R");
my @chromosomes=("2L", "2R", "3L", "3R", "4", "X");



##------------------------------------------------------------
## options
use Getopt::Long;
use Pod::Usage;

my $help;
my $man;
my $quiet;
my $verbose;

my $random_n;
my $all;

## Getopt::Long::Configure("no_ignore_case");

GetOptions(
    "verbose" => \$verbose,
    "quiet" => \$quiet,
    "help"=> \$help,
    "man" => \$man,
    "random=i" => \$random_n,
    "all" => \$all,
    "context=i" => \$context_size
    ) || pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

### ------------------------------------
## subs

## find a sequence $seq in an alignment, given the directory $dir containing alignments in multiple fasta
## return the file, the start, and the end alignment column
sub find_in_alignments {
    my ($dir,$seq,$start_pos,$end_pos) = @_;

    local *DIR;

    opendir(DIR, $dir) || die "Cannot read directory $dir";

    my @files = grep(/.mfa.gz$/, readdir(DIR));

    closedir(DIR);

    for my $alnfile (@files) {

	my $aln_start_pos;
	my $aln_end_pos;

	if ($alnfile=~/(\d+)_(\d+).mfa.gz$/) {
	    $aln_start_pos=$1;
	    $aln_end_pos=$2;
	} else {
	    ## ignore this strange file
	    print STDERR "WARNING: ignore $alnfile in $dir\n";
	    next;
	}

	## test whether alignment file covers the searched sequence
	if (! ($aln_start_pos<=$start_pos && $aln_end_pos>=$end_pos)) {next;}

	# print "FILE: $alnfile\n";

	my $alignment_in = Bio::SeqIO->new( '-file' => "gunzip -c $dir/$alnfile |",
					    '-format' =>  "fasta");

	## the dmel reference sequence alignment string will be the first one
	my $ref_seq = $alignment_in->next_seq();

	#print $ref_seq->id()."\n";
	#print $ref_seq->seq()."\n";

	## make LocatableSeq out of sequence (actually our sequence is an alignment string)
	#
	my $ref_locseq = Bio::LocatableSeq->new(-seq => $ref_seq->seq(),
						-id  => "DroMel_CAF1",
						-start => $aln_start_pos,
						-end   => $aln_end_pos);
	## remove gaps
        #
	my $ref_locseq_wo_gaps = $ref_locseq->seq();
	$ref_locseq_wo_gaps =~ tr/-//d; ## remove gaps
	$ref_locseq_wo_gaps = uc $ref_locseq_wo_gaps; ## remove gaps

	# print "Search in $ref_locseq_wo_gaps\n";

	#### ============================================================
	#### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	####
        #### ATTENTION: Search is not necessary,
	#### since pos can be computed!
	####
	#### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	#### ============================================================

	## find $seq in $ref_locseq_wo_gaps
	#
	my $pos = index( $ref_locseq_wo_gaps,$seq->seq,0 )+1;



	## the position can be mapped back to the alignment!
	if ($pos >=1) {
	    ## seq is found in alignment

	    my $start_col = $ref_locseq->column_from_residue_number($pos+$aln_start_pos);
	    my $end_col =  $ref_locseq->column_from_residue_number($pos+$aln_start_pos+length($seq->seq)-2);

	    # print "FOUND seq in $alnfile: $start_col-$end_col.\n";

	    return ($alnfile,$start_col,$end_col);
	} else {
	    # print "Not found in $alnfile.\n";
	}
    }

    return ("",-1,-1);
}



# slice_sequences_with_context($aln,$start_col,$end_col,$context)
#
# slice alignment from start_col to end_col
# collect the sequences that have sufficiently many non-gap characters
# in the slice with context of +/- $context positions
#
# write mfasta file $mfa_outfile
# write alignment slices to $slice_outfile
sub slice_sequences_with_context {
    my ($mfa_outfile,$slice_outfile,$aln,$start_col,$end_col,$context_size,$max_gap_content)=@_;

    open(my $MFAOUT, ">", "$mfa_outfile") || die "Cannot write to $mfa_outfile: $!";
    open(my $SLCOUT, ">", "$slice_outfile") || die "Cannot write to $slice_outfile: $!";

    ## write clustalw header
    print $SLCOUT "CLUSTAL W (1.83)\n\n";

    # go through all sequences in the alignment, look at slice
    # sufficiently many non-gap?
    # if yes, get slice with context (involves some mapping from col to residue and back)
    #
    foreach my $seq ( $aln->each_seq() ) {

	my $seq_slice = uc($seq->subseq($start_col,$end_col));

	## get number of gaps in $seq_slice
	my $seq_slice_wo_gaps=$seq_slice;
	$seq_slice_wo_gaps =~ tr/-//d;
	my $no_gaps = length($seq_slice)-length($seq_slice_wo_gaps);


	my $gap_ratio = $no_gaps/length($seq_slice);

	#print $seq->id()." $gap_ratio\n";

	#print "$seq_slice\n";

	if ($no_gaps<length($seq_slice)) {
	    ## write alignment slice
	    print $SLCOUT $seq->id()." $seq_slice\n";
	}

	if ($gap_ratio <= $max_gap_content) {

	    my $seq_w_gaps=$seq->seq(); ## get sequence as string

	    my $seq_wo_gaps = $seq_w_gaps;  ## sequence without gaps
	    $seq_wo_gaps =~ tr/-//d; # remove all gaps


	    my $aln_len=length($seq_w_gaps); ## length of the alingment string (sequence with gaps)
	    my $seq_len=length($seq_wo_gaps); ## length of sequence without gaps

	    # determine start and end index from start_col and end_col
	    my $start=-1;
	    my $end=$seq_len-1;


	    my $idx=0; ## index in sequence (without gaps)
	    for (my $i=1; $i<=$aln_len; $i++) {

		if (substr($seq_w_gaps,$i-1,1) ne "-") {$idx++;}

		if ($i==$start_col) {$start=$idx;}
		if ($i>=$end_col) {$end=$idx;last;} # stop as soon as we get to or over the end column
	    }

	    #print "$start $end\n";

	    my $start_w_context = $start-$context_size;
	    my $end_w_context = $end+$context_size;

	    if ($start_w_context<=0) {$start_w_context=0;}
	    if ($end_w_context>=length($seq_wo_gaps)) {$end_w_context = length($seq_wo_gaps);}

	    #print "$start_w_context $end_w_context\n";


	    my $seq_slice_w_content = uc(substr($seq_wo_gaps,$start_w_context,$end_w_context-$start_w_context+1));

	    print $MFAOUT ">".$seq->id()." left_context=".($start-$start_w_context)."; right_context=".($end_w_context-$end)."; gap_ratio=$gap_ratio\n";
	    print $MFAOUT "$seq_slice_w_content\n";

	}
    }
    close $MFAOUT;
    close $SLCOUT;
}


# print a selected rnaz hit (+ annotation) in list and store the mfa and aln files for the hit
sub store_rnaz_hit {
    #parameter
    my ($locus,$chromosome,$dmel_seq,$start,$end,$orient,$ann_id,$ann_name,$ann_start,$ann_end,$ann_orient) = @_;

    print "$locus $chromosome $start $end $orient";
    if (defined($ann_id)) {
	print " $ann_id $ann_name $ann_start $ann_end $ann_orient";
    }
    print "\n";

    ## .$hit_tab[8]." ".$hit_tab[9]." ".$hit_tab[10]." ".$hit_tab[11]." ".$hit_tab[12]." ".$hit_tab[13]." ".$hit_tab[14]." ".$hit_tab[15]." ".$hit_tab[16]."\n";

    ## ------------------------------
    ## look for the RNAz hit in the pecan alignment, get alignment slice,
    ## determine well-aligned sequences and sequences with context
    ## in the alignment

    ## get drosophila sequence of the RNAz prediction
    my $RNAz_hit_sequence =  $dmel_seq->trunc($start,$end);

    #print "NOW FIND\n";

    my ($alnfile,$start_col,$end_col)
	= find_in_alignments("$datadir/Alignments/$chromosome-pecan-CAF1",$RNAz_hit_sequence,$start,$end);

    #print "PECAN: $alnfile $start_col $end_col\n";

    # load alignment
    my $aln_in = Bio::AlignIO->new( '-file' => "gunzip -c $datadir/Alignments/$chromosome-pecan-CAF1/$alnfile |",
				    '-format' => "fasta" );
    my $aln = $aln_in->next_aln();

    #print "\n";

    ## generate name from locus and ann_id
    my $name = $chromosome.":".$locus;

    my $locus_outfile="$mfa_outdir/$name.mfa";
    my $slice_outfile="$aln_outdir/$name.aln";

    ###
    my $sequences = slice_sequences_with_context($locus_outfile,$slice_outfile,$aln,$start_col,$end_col,$context_size,$max_gap_content);

    #print "\n\n";
}

# draw n random instances from all RNAz hits
sub draw_and_store_random_instances {
    my ($datadir,$mfa_outdir,$aln_outdir,$random_n,$max_gap_content) = @_;

    # count total number of RNAz hits
    my $number_of_hits=0;
    foreach my $chromosome (@chromosomes) {
	## read in RNAz hits for the chromosome
	open(my $RNAZHITS, "-|", "gunzip -c $datadir/Annotation/$chromosome-pecan-CAF1.annotation.dat.gz | grep '^locus'")
	    || die "Cannot read or gunzip $datadir/Annotation/$chromosome-pecan-CAF1.annotation.dat.gz";
	my @rnazhits=<$RNAZHITS>;
	close $RNAZHITS;
	$number_of_hits += $#rnazhits+1;
    }
    # print "Total number of RNAz hits: $number_of_hits\n";

    my @numbers = draw_random_numbers($random_n,$number_of_hits);
    @numbers = sort { $a <=> $b } @numbers;
        # count total number of RNAz hits
    my $hit_idx=0;
    my $num_idx=0;
    foreach my $chromosome (@chromosomes) {
	## read in dmel chromosome
	my $dmel_seqin = Bio::SeqIO->new( '-file' => "gunzip -c $datadir/Dmel-r4.3/dmel-$chromosome-chromosome-r4.3.fasta.gz |",
					  '-format' => "fasta" );
	my $dmel_seq = $dmel_seqin->next_seq;  # we assume there is only one sequence


	## read in RNAz hits for the chromosome
	open(my $RNAZHITS, "-|", "gunzip -c $datadir/Annotation/$chromosome-pecan-CAF1.annotation.dat.gz | grep '^locus'")
	    || die "Cannot read $datadir/Annotation/$chromosome-pecan-CAF1.annotation.dat.gz";
	my @rnazhits=<$RNAZHITS>;
	close $RNAZHITS;

	# while the next random number refers to a hit of the current chromosome
	while ( $num_idx < $random_n && $numbers[$num_idx] <= $hit_idx+$#rnazhits+1 ) {

	    # print STDERR "$num_idx Store $numbers[$num_idx]\n";

	    my @hit_tab = split /\s+/, $rnazhits[$numbers[$num_idx]-$hit_idx];
	    my $locus=$hit_tab[0];
	    my $start=$hit_tab[2];
	    my $end=$hit_tab[3];

	    store_rnaz_hit($locus,$chromosome,$dmel_seq,$start,$end,"?");

	    $num_idx++;
	}

	$hit_idx += $#rnazhits+1;
    }
}

# extract all RNAz hits
sub extract_all_instances {
    my ($datadir,$mfa_outdir,$aln_outdir,$max_gap_content) = @_;

    foreach my $chromosome (@chromosomes) {
	## read in dmel chromosome
	my $dmel_seqin = Bio::SeqIO->new( '-file' => "gunzip -c $datadir/Dmel-r4.3/dmel-$chromosome-chromosome-r4.3.fasta.gz |",
					  '-format' => "fasta" );
	my $dmel_seq = $dmel_seqin->next_seq;  # we assume there is only one sequence


	## read in RNAz hits for the chromosome
	open(my $RNAZHITS, "-|", "gunzip -c $datadir/Annotation/$chromosome-pecan-CAF1.annotation.dat.gz | grep '^locus'")
	    || die "Cannot read $datadir/Annotation/$chromosome-pecan-CAF1.annotation.dat.gz";
	my @rnazhits=<$RNAZHITS>;
	close $RNAZHITS;

	# while the next random number refers to a hit of the current chromosome
	foreach my $rnazhit ( @rnazhits ) {

	    my @hit_tab = split /\s+/, $rnazhit;
	    my $locus=$hit_tab[0];
	    my $start=$hit_tab[2];
	    my $end=$hit_tab[3];

	    store_rnaz_hit($locus,$chromosome,$dmel_seq,$start,$end,"?");

	}
    }
}




### ------------------------------------------------------------
## main part
#
if (defined($random_n)) {
    $mfa_outdir.=".random";
    $aln_outdir.=".random";
}

if (! -e $mfa_outdir) {
    mkdir "$mfa_outdir";
}

if (! -e $aln_outdir) {
    mkdir "$aln_outdir";
}

if (defined($random_n)) {

    # run in special mode
    draw_and_store_random_instances($datadir,$mfa_outdir,$aln_outdir,$random_n,$max_gap_content);
    exit;
}


if (defined($all)) {
    # run in special mode
    extract_all_instances($datadir,$mfa_outdir,$aln_outdir,$max_gap_content);
    exit;
}


# else: default mode --- look for FB-annotated hits

foreach my $chromosome (@chromosomes) {

    ## read in dmel chromosome
    my $dmel_seqin = Bio::SeqIO->new( '-file' => "gunzip -c $datadir/Dmel-r4.3/dmel-$chromosome-chromosome-r4.3.fasta.gz |",
				    '-format' => "fasta" );
    my $dmel_seq = $dmel_seqin->next_seq;  # we assume there is only one sequence


    ## read in RNAz hits for the chromosome
    open(my $RNAZHITS, "-|", "gunzip -c $datadir/Annotation/$chromosome-pecan-CAF1.annotation.dat.gz | grep '^locus'")
	|| die "Cannot read or gunzip $datadir/Annotation/$chromosome-pecan-CAF1.annotation.dat.gz";
    my @rnazhits=<$RNAZHITS>;
    close $RNAZHITS;

    ## read in Flybase ncRNA annotation for the chromosome
    my $flybase_ncRNA_anno =
	Bio::SeqIO->new( '-file' => "gunzip -c $datadir/Flybase-Annotation/dmel-$chromosome-ncRNA-r4.3.fasta.gz |",
			 '-format' => "fasta" );


    ## list RNAz hits that overlap with flybase annotation,
    ## report RNAz location, dmel sequence, Flybase annotation,...

    while (my $annseqobj = $flybase_ncRNA_anno->next_seq) {

	my $ann_id=$annseqobj->id();
	my $seqdesc=$annseqobj->description();
	#print "$ann_id $seqdesc\n";

	my %d = parse_fasta_description($seqdesc);

	my $ann_name = $d{"Name"};

	my %loc= parse_location($d{"loc"});

	my $ann_start=$loc{"start"};
	my $ann_end=$loc{"end"};

	my $ann_len=$d{"len"};

	my $ann_orient=($d{"complement"}==0)?"+":"-";

	if ($ann_len > $max_len) {next;} ## skip annotations that are too long

	#print "$ann_id ".$d{"Name"}." ".$loc{"chr"}." ".$loc{"start"}." ".$loc{"end"}." ".$loc{"complement"}."\n";

	foreach my $hit ( @rnazhits ) {
	    my @hit_tab = split /\s+/, $hit;
	    my $locus=$hit_tab[0];
	    my $start=$hit_tab[2];
	    my $end=$hit_tab[3];

	    # test for overlap of annotation and RNAz hit, if yes: report
	    if (ranges_overlap($ann_start,$ann_end,$start,$end)) {

		store_rnaz_hit($locus,$chromosome,$dmel_seq,$start,$end,"?",$ann_id,$ann_name,$ann_start,$ann_end,$ann_orient);

	    }
	}
    }
}

## ------------------------------------------------------------


