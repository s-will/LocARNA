package LocARNA_RNAz;

use 5.008003;
use strict;
use warnings;

use FindBin;
use lib $FindBin::Bin;
use RNAz; # the RNAz module

require Exporter;

our @ISA = qw(Exporter);

our @EXPORT_OK = ();

our $maxGapWindow        = 0.25;
our $maxMaskedWindow     = 0.1;

our $maxGapLocus        = 0.25;
our $maxMaskedLocus     = 0.1;

our @EXPORT = qw(
parse_fasta_description
parse_location
ranges_overlap
find_PECAN_alignment_file
draw_random_numbers
read_alignment
read_alignment_fasta
nice_sequences_win_and_loc
write_alignment
project_alignment
);


## !!!!!!!!!!!!!!!!!!!!
## CONVENTION: the first sequence in an alignment is the reference sequence


########################################
## read_alignment($filename)
##
## Read an alignment from file (containing a single alignment)
## If file ends with .gz, the file is automatically uncompressed
##
## relies on RNAz module routines
##
## returns alignment in RNAz representation
##
########################################
sub read_alignment {
    my $filename = shift;

    my $fh;

    if ($filename =~ /\.gz$/) {
        open($fh, "-|", "gunzip -c $filename") || die "Cannot read or uncompress file $filename: $!";
    } else {
        open($fh, "<", $filename) || die "Cannot read file $filename: $!";
    }

    my $alnFormat = checkFormat($fh);

    my $alnString = getNextAln( $alnFormat, $fh );

    close $fh;

    my $aln = parseAln($alnString,$alnFormat);

    return $aln;
}

########################################
## read_alignment_fasta($filename)
##
## Read an alignment from file containing a single alignment in fasta format.
## If file ends with .gz, the file is automatically uncompressed
##
## returns alignment in RNAz representation
##
########################################
sub read_alignment_fasta {
    my $filename = shift;

    my $fh;

    if ($filename =~ /\.gz$/) {
        open($fh, "-|", "gunzip -c $filename") || die "Cannot read or uncompress file $filename: $!";
    } else {
        open($fh, "<", $filename) || die "Cannot read file $filename: $!";
    }

    my @aln;

    my $line=<$fh>;
    while(defined($line)) {
        if ($line=~/^>\s*(.+)/) {
            my $name=$1;
            $name =~ s/\s+$//g; ## accept and ignore trailing spaces in name

            my $seq="";
            while (defined($line=<$fh>) && ($line !~ /^>/)) {
                chomp $line;
                $line =~ s/\s+//g;
                $seq .= $line;
            }
            push @aln, {name=>$name,
                        seq=>$seq,
                        start=>0,
                        end=>0 # signals that position information is not used
            };
        } else {
            $line=<$fh>;
        }
    }

    close $fh;

    return \@aln;
}

########################################
## project_alignment($aln,$sequence_ids)
## project alignment $aln to alignment that contains
## only the ids in $sequence_ids
## do not remove common gaps
##
## perform shallow copy
##
## @returns ref of loh alignment
########################################
sub project_alignment {
    my ($aln,$sequence_ids) = @_;
    my @res=();
    for my $i (@$sequence_ids) {
        push @res, $aln->[$i];
    }
    return \@res;
}

########################################
## write_alignment($filename,$format,$aln)
## ------------------------------
## Write an alignment to file (in given format)
##
## $filename name of target file
## $aln reference to RNAz alignment representation
## $format ... CLUSTAL,  FASTA, MAF
##
########################################
sub write_alignment {
    my ($filename, $format, $aln) = @_;

    my $fh;

    my $format_string = formatAln($aln, $format);

    open($fh, ">", "$filename") || die "Cannot write file $filename: $!";
    print $fh $format_string;
    close $fh;
}


## parses the description in a fasta file and converts to hash.
## Takes string [<name>=<feature>; ]*
sub parse_fasta_description {
    my ($desc_str)=@_;

    my @desc = split (/\s*;\s*/, $desc_str);

    my %h;

    foreach my $d (@desc) {
        if ($d=~/([^=]+)=(.*)/) {
            $h{$1}=$2;
        }
    }
    return %h;
}

## parses a location of the form <chr>:<start>..<end> or <chr>:complement(<start>..<end>)
## and returns hash
sub parse_location {
    my ($loc_str)=@_;
    my %h;

    if ($loc_str=~/([^:]*):(\d+)\.\.(\d+)/) {
        $h{"chr"}=$1;
        $h{"start"}=$2;
        $h{"end"}=$3;
        $h{"complement"}=0;
    } elsif ($loc_str=~/([^:]*):complement\((\d+)\.\.(\d+)\)/) {
        $h{"chr"}=$1;
        $h{"start"}=$2;
        $h{"end"}=$3;
        $h{"complement"}=1;
    }
    else {
        die "Cannot parse location string $loc_str\n";
    }
    return %h;
}


# test whether two ranges i..j and k..l overlap
sub ranges_overlap {
    my ($i,$j,$k,$l) = @_;

    return
        ($i<=$k && $k<=$j) # k in i..j
        ||
        ($i<=$l && $l<=$j) # l in i..j
        ||
        ($k<=$i && $i<=$l) # i in k..l
        ||
        ($k<=$j && $j<=$l) # j in k..l
        ;
}


########################################
## find_PECAN_alignment_file($dir,$start_pos,$end_pos)
##
## find an alignment file in an PECAN alignment
##
## returns triple of filename, start, and end position
##
########################################
sub find_PECAN_alignment_file {
    my ($dir,$start_pos,$end_pos) = @_;

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
        if (($aln_start_pos<=$start_pos && $aln_end_pos>=$end_pos)) {return ($alnfile,$aln_start_pos,$aln_end_pos);}
    }
    return "";
}


## check whether a sequence is 'nice'
## that is whether it has sufficiently few gaps and masked characters
##
## $seq alignment string that is tested
## max ratio of gaps to sequence length
## max ratio of masked characters to sequence length
sub is_nice_sequence {
    my ($seq, $maxGap,$maxMasked)=@_;

    my $no_gaps = ( $seq =~ tr/-./-/ );
    my $no_masked = ( $seq =~ tr/Na-z/Na-z/ ); # N and small letters count as masked!

    my $gap_ratio = $no_gaps/length($seq);
    my $masked_ratio = $no_masked/length($seq);

    return $gap_ratio <= $maxGap && $masked_ratio <= $maxMasked;
}

## check whether a sequence is 'nice' compared to a reference sequence
## that is whether it has sufficiently few insertions/deletions
## compared to the reference and masked characters
##
## $seq   alignment string that is tested
## $refseq alignment string of reference
## max ratio of gaps to sequence length
## max ratio of masked characters to sequence length

sub is_nice_sequence_compared_to_ref {
    my ($seq, $refseq, $maxGap,$maxMasked)=@_;

    ## determine gap ratio from aligmnent where common gaps are removed
    my @tmpAln = ( { seq => $refseq }, { seq => $seq } );
    removeCommonGaps( \@tmpAln );
    my $numGaps0 = ( $tmpAln[0]->{seq} =~ tr/-./-/ );
    my $numGaps1 = ( $tmpAln[1]->{seq} =~ tr/-./-/ );

    my $tmpLength = length( $tmpAln[0]->{seq} );

    my $gap_ratio = ($numGaps0+$numGaps1)/$tmpLength;

    ## determine masked ratio from sequence directly
    my $no_masked = ( $seq =~ tr/Na-z/Na-z/ ); # N and small letters count as masked!
    my $masked_ratio = $no_masked/length($seq);

    ## both ratios must be sufficiently small
    return $gap_ratio <= $maxGap && $masked_ratio <= $maxMasked;
}


## computes the subset of 'nice' sequence in an alignment that have
## sufficiently few gaps and masked symbols.
## return the sequence indices as a list
sub nice_sequences {
    my ($aln, $maxGap,$maxMasked)=@_;
    my @seqs=();

    if ( is_nice_sequence($aln->[0]{seq},$maxGap,$maxMasked) ) {

        push @seqs,0;

        foreach my $i ( 1 .. @$aln-1 ) {
            push(@seqs, $i) if is_nice_sequence_compared_to_ref($aln->[$i]{seq},$aln->[0]{seq},$maxGap,$maxMasked);
        }

    }

    return @seqs;
}


########################################
## nice_sequences_win_and_loc($aln, $start_pos, $locuslines)
## ------------------------------
##
## Very specialized from of pruning an alignment that removes sequences
## that have too many gaps or masked characters in all windows and the locus
##
## $aln        reference to RNAz aligmnent rep
## $start_pos  location offset of the reference sequence in the alignment
## $locuslines array of strings describing the locus and its windows
##   first line is for the locus
##   next lines describe windows, the 4th and 5th position are window start and end respectively
##
## returns ref of list of indices of nice sequences in $aln
##
########################################
sub nice_sequences_win_and_loc {
    my ($aln, $locuslines) = @_;

    my %niceseqs; ## hash that will contain an entry for each id of a nice sequence

    foreach my $s (nice_sequences($aln, $maxGapLocus, $maxMaskedLocus)) {$niceseqs{$s}=1;}

    my $wstart_pos=-1;
    my $wend_pos=-1;

    for (my $i=1; $i < @$locuslines; $i++) {
        my @window = split(/\s+/, $locuslines->[$i]);

        ## skip check for successing identic window
        next if ($window[3]==$wstart_pos && $window[4]==$wend_pos);

        $wstart_pos = $window[3];
        $wend_pos = $window[4];

        my $aln_window = sliceAlnByPos($aln, 0, $wstart_pos,$wend_pos);

        foreach my $s (nice_sequences($aln_window, $maxGapWindow, $maxMaskedWindow)) {$niceseqs{$s}=1;}
    }

    ## remove all sequences that are not in @niceseqs

    # copy only the nice sequences of aln to the pruned alignment
    #my @pruned_aln = ();
    #foreach my $i (sort {$a<=>$b} (keys %niceseqs)) {
    #   push @pruned_aln, $aln->[$i];
    #}

    #removeCommonGaps(\@pruned_aln);

    #return \@pruned_aln;

    my @nice_sequences_indices = sort {$a<=>$b} (keys %niceseqs);

    return \@nice_sequences_indices;
}



# draw $k random numbers from the set {0,...,N-1}
# return as (unsorted) list
sub draw_random_numbers {
    my ($k, $N) = @_;

    ## randomly draw $n different numbers from 1..$N
    my @bag=();
    for (my $i=0; $i < $N; $i++) {
        push @bag, $i;
    }

    my @numbers=();
    for (my $i=0; $#bag>=0 && $i<$k; $i++) {
        my $r=int(rand ($#bag+1));

        push @numbers,$bag[$r];

        $bag[$r]=$bag[-1];
        $#bag--;
    }
    return @numbers;
}


# ########################################
# ## bioaln_to_rnaz($aln)
# ## ------------------------------
# ##
# ## Translate a bioperl Bio::SimpleAlign object to
# ## the alignment representation of RNAz
# ## (as used in perl module RNAz.pm)
# ##
# ## RNAz alignment representation
# ## ------------------------------
# ## an alignment is a hashlist of hashs
# ## with entries name and seq, start, end, strand
# ##   other hash keys are used in RNAz as dead, chrom, org, fullLength
# ## Sequences can containt gap symbols '-' or '.'.
# ## Small letters are considered as masked.
# ##
# ## The RNAz alignment thus can represent order of sequenes, sequence
# ## strings and ids, and to a limited extent sequence location
# ##
# ## ------------------------------
# ## $aln alignment as Bio::SimpleAlign object
# ##
# ## returns alignment in RNAz internal format
# ##
# ########################################
# sub bioaln_to_rnaz {
#     my ($aln) = @_;

#     my @RNAzaln=();

#     foreach my $seq ( $aln->each_seq ) {
#       push @RNAzaln, { name=> $seq->id,
#                        seq => $seq->seq,
#                        start => $seq->start,
#                        stop => $seq->stop,
#                        strand => $seq->strand };
#     }

#     return \@RNAzaln;
# }


# ########################################
# ## rnazaln_to_bio($aln)
# ## ----------------------------------------
# ## convert alignment in RNAz internal representation
# ## to bioperl Bio::SimpleAlign object
# ##
# ## for each sequence convert information on id,seq,start,stop,strand
# ##
# ########################################
# sub rnazaln_to_bio {
#     my ($aln) = @_;

#     my $bioaln = new Bio::SimpleAlign();

#     foreach my $seq ($aln) {
#       my $seq = new Bio::LocatableSeq(
#           -id => $seq{name},
#           -seq => $seq{seq},
#           -start => $seq{start},
#           -end => $seq{end},
#           -strand => $seq{strand}
#           );
#       $bioaln->addSeq($seq);
#     }
# }
