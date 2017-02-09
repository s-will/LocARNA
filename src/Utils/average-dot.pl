#!/usr/bin/env perl

=head1 NAME

average-dot.pl

=head1 SYNOPSIS

average-dot.pl [options] --alignment <pp-file> --sequences <pp-file 1> --sequences <pp-file 2>...

Options:
   --help               brief help message

   --man                full documentation

   --verbose            be verbose

   --quiet              be quiet

   --alignment <file>   pp file representing the (Lo)Carna alignment

   --sequences <file>   pp files for the input sequences

   --outfile <file>     output file (DEFAULT: averagedot)

   --threshold <float>  for each given threshold (this option can be specified more than once)
                        a copy of the orginal dotplots is created where the dots are
                        highlighted that have in the average plot a probability above the
                        threshold.


=head1 EXAMPLE

The following call

./average-dot.pl -t 0.2 -t 0.5 -t 0.75 -a result.pp -s seq1  -s seq2 -s seq3

creates given a result.pp (as found in the folder "results" of
mlocarna) and input pp files seq1 seq2 seq3 (as found in the folder
"input" of mlocarna) the following files

averagedot_comb_20.pp  # average dot plot of the sequences where dots >=0.2 are highlighted
averagedot_comb_50.pp  # average dot plot of the sequences where dots >=0.5 are highlighted
averagedot_comb_75.pp  # average dot plot of the sequences where dots >=0.75 are highlighted

seq1.pp         # dot plot of sequence 1 projected to the alignment with no highlights
seq1_comb_20.pp # dot plot of sequence 1 projected to the alignment and dots highlighted
                # that occur in the average with a probability >=0.2
                # in the lower right triangle the average with colored variance info is given
seq1_comb_50.pp # ... same as previous for probability >=0.5
seq1_comb_75.pp # ... same as previous for probability >=0.5

+ analogous files for seq2 and seq3


All pp files have to be in format "PP 2.0"!

=head1 DESCRIPTION
=cut


############################################################
### Options
###
use warnings;
use strict;

use Getopt::Long;
use Pod::Usage;

use List::Util qw[min max];

use FindBin;
use lib "$FindBin::Bin/../lib/perl";

use MLocarna;

my $verbose;
my $quiet;
my $help;
my $man;
my $ofile="averagedot";
my $alignmentppFile;
my @sequenceppFiles;
my @thresholds;
GetOptions("verbose" => \$verbose,
	   "quiet" => \$quiet,
	   "help"=> \$help,
	   "man" => \$man,
	   "alignment=s" => \$alignmentppFile,
	   "sequences=s" => \@sequenceppFiles,
	   "outfile=s" => \$ofile,
	   "threshold=f" => \@thresholds,
	   ) || pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;


##################################################################################

sub write_pp_annotated {
    my ($filename,$sequence,$probUR,$probLL,$threshold,$colorLL)=@_;

    open(my $OUT, ">", "$filename") || die "Cannot write to $filename: $!";

    print $OUT "#PP 2.0\n\n";

    while ( my ($name, $row) = each %$sequence ) {
	print $OUT "$name\t\t$row\n";
    }
    print $OUT "\n#END\n\n#SECTION BASEPAIRS\n\n";

    foreach my $id (keys %{$probUR}){
	my $highlight = 0;
	if($probLL->{$id} > $threshold){ $highlight = 1;}

	my $fileLine = "$id $probUR->{$id} $probLL->{$id} $highlight $colorLL->{$id}\n";
	print $OUT $fileLine;
	print $fileLine if $verbose;
    }

    print $OUT "\n#END\n";

    close $OUT;
}



############################################################
### main
###

print "read alignment\n" if $verbose;
my %alignment = read_pp_file_aln_wo_anno($alignmentppFile);
my %sequenceData; # a hash from sequence name to sequence string
my %sequenceDotPlots; # dot plots with indices w.r.t. individual sequence
my %alignmentDotPlots; # dot plots with indices w.r.t. alignment

foreach my $filename (@sequenceppFiles){
    print "read sequence $filename\n" if $verbose;
    my %sequenceHash = read_pp_file_aln_wo_anno($filename);
    my $sequenceName = (keys %sequenceHash)[0];
    print "Sequence name = $sequenceName\n" if $verbose;
    my %seqDotPlot = read_pp_file_pairprobs($filename);


    ## convert the indices from sequence space to alignment space
    my $row = $alignment{$sequenceName};
    print "alignment row:$row\n" if $verbose;
    my @seq2aln; # mapping from sequence positions to alignment positions
    my $offset=0;
    my $alIndex =0;
    my $seqIndex =0;
    foreach my $character (split //, $row) {
	if($character ne "-"){
	    #print "$seqIndex --> $alIndex\n" if $verbose;
	    $seq2aln[$seqIndex]=$alIndex;
	    $seqIndex++;
	}
	$alIndex++;
    }
    ## now use the mapping to compute a DotPlot on alignment indices
    my %alnDotPlot=();
    for my $k (keys %seqDotPlot) {
	$k =~ /^(\d+)\s+(\d+)/ or die "invalid index in dot plot:'$k'\n";
	my $sx = $1;
	my $sy = $2;
	my $ax = $seq2aln[$sx-1]+1;
	my $ay = $seq2aln[$sy-1]+1;
	#print "($sx,$sy)-->($ax,$ay)\n" if $verbose;
	$alnDotPlot{"$ax $ay"} = $seqDotPlot{$k};
    }

    $sequenceData{$sequenceName} = \%sequenceHash;
    $sequenceDotPlots{$sequenceName} = \%seqDotPlot;
    $alignmentDotPlots{$sequenceName} = \%alnDotPlot;
}


# compute for each base pair of the alignment the sum of the
# probabilities and the sum of the squares of the probabilities (of
# the individual sequences
my %sumProbs =();
my %sumSquareProbs =();

foreach my $seq (keys %alignmentDotPlots){
    print "average sequence $seq\n" if $verbose;
    foreach my $id (keys %{$alignmentDotPlots{$seq}}){
	if(! exists $sumProbs{$id}){
	    $sumProbs{$id}=0;
	}
	if(! exists $sumSquareProbs{$id}){
	    $sumSquareProbs{$id}=0;
	}
	my $value =  $alignmentDotPlots{$seq}{$id};
	$sumProbs{$id} = $sumProbs{$id}+ $value;
	#print "sumProbs($id)=$sumProbs{$id}\n" if $verbose;
	$sumSquareProbs{$id} = $sumSquareProbs{$id}+ $value*$value;
    }
}

# compute for each base pair of the alignment the average probability and
# the scaled standard deviation

my $numSequences = @sequenceppFiles;
my %averageProbs=();
my %scaledDeviation=();

foreach my $id (keys %sumProbs){
    my $sumProb = $sumProbs{$id};
    my $sumSquareProb = $sumSquareProbs{$id};

    my $average = $sumProb/$numSequences;
    my $squaresAverage =  $sumSquareProb/$numSequences;
    my $deviation = sqrt(max(0,$squaresAverage-$average*$average));

    $averageProbs{$id}=$average;
    $scaledDeviation{$id}= 1 -2* $deviation;
}


#################### OUTPUT


write_pp_annotated("$ofile.pp",\%alignment,\%averageProbs,\%averageProbs,2,\%scaledDeviation);
foreach my $filename (@sequenceppFiles){
    my %sequenceData = read_pp_file_aln_wo_anno($filename);
    my $sequenceName = (keys %sequenceData)[0];
    my %tempHash = ($sequenceName => $alignment{$sequenceName});
    print "write file for: $sequenceName\n" if $verbose;
    write_pp_annotated($filename.".pp",
		       \%tempHash,
		       $alignmentDotPlots{$sequenceName},
		       \%averageProbs,2,\%scaledDeviation);
}

foreach my $threshold (@thresholds){
    print "write average file for threshold $threshold\n" if $verbose;
    write_pp_annotated($ofile."_comb_".($threshold*100).".pp",
		       \%alignment,\%averageProbs,\%averageProbs,
		       $threshold,\%scaledDeviation);

    foreach my $filename (@sequenceppFiles){
	my %sequenceData = read_pp_file_aln_wo_anno($filename);
	my $sequenceName = (keys %sequenceData)[0];
	my %tempHash = ($sequenceName => $alignment{$sequenceName});
	print "write file for: $sequenceName and threshold $threshold\n" if $verbose;
	write_pp_annotated($filename."_comb_".($threshold*100).".pp",
			   \%tempHash,
			   $alignmentDotPlots{$sequenceName},
			   \%averageProbs,$threshold,\%scaledDeviation);
    }
}

print "DONE\n" if $verbose;


## ------------------------------------------------------------

__END__


