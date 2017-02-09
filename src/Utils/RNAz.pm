package RNAz;

use 5.008003;
use strict;
use warnings;
require Exporter;

our @ISA = qw(Exporter);

our @EXPORT_OK = ();

our $rnazVersion='1.0';

our @EXPORT = qw(checkFormat
                                 getNextAln
                                 formatAln
                                 readMAF
                                 readClustal
                                 parseAln
                                 sliceAlnByColumn
                                 sliceAlnByPos
                                 alnCol2genomePos
                                 genomePos2alnCol
                                 removeCommonGaps
                                 rangeWarn
                                 revAln
                                 meanPairID
                                 pruneAln
                                 getNextRNAz
                                 parseRNAz
                                 shuffleAln
                                 getSeq
                                 blastSeq
                                 niceNumber);

# Set version of current RNAz package


our $VERSION = '0.1';


sub checkFormat{

  my $fh;

  $fh=shift;
  if (!defined $fh){
        $fh=*STDIN;
  }

  while(<$fh>){
        next if /^\s*\#/;
        next if /^\s*$/;
        if (/CLUSTAL/i){
          return "CLUSTAL";
        } elsif (/^a/){
          return "MAF"
        } else {
          return "UNKNOWN";
        }
  }
}

######################################################################
#
# readMAF($string)
#
# Converts the MAF in string to internal alignment format. Returns
# list of usual array references.
#
######################################################################

sub readMAF{

  my $string=shift;

  return [] if $string eq '';

  my @input=split("\n",$string);

  my @aln=();

  foreach my $i (0..$#input){

        $_=$input[$i];

        next if (/^\s?\#/);
        next if (/^\s?a/);

        if (/^\s?s/) {
          (my $dummy, my $name, my $start, my $length,
           my $strand, my $fullLength, my $seq)=split;

          my $end=$start+$length;

          $seq=~s/\./-/g;

          my $row={name=>$name,
                           start=>$start,
                           end=>$end,
                           fullLength=>$fullLength,
                           seq=>$seq,
                           strand=>$strand};

          if ($name=~/^(.*)\.(.*)$/){
                $row->{org}=$1;
                $row->{chrom}=$2;
          }

          push @aln, $row;
        }
  }
  return \@aln;
}

sub readClustal{

  my $inString=shift;

  my ($order,%order);
  my %input;

  my @lines = split /^/, $inString;

  foreach (@lines){
        next if ( /^\s+$/ );
        next if ( /^\/\/\$/); # ignore block separators '//'
        my ($seqname, $aln_line) = ('', '');
        if ( /^\s*(\S+)\s*\/\s*(\d+)-(\d+)\s+(\S+)\s*\d*\s*$/ ) {
          # clustal 1.4 format
          #($seqname,$aln_line) = ("$1/$2-$3",$4);
          ($seqname,$aln_line) = ($1,$4);
          $seqname.="\_$2\_$3";
        } elsif ( /^(\S+)\s+([A-Z\-]+)\s*\d*\s*$/i ) {

          ($seqname,$aln_line) = ($1,$2);
        } else {
          next;
        }

        if ( !exists $order{$seqname} ) {
          $order{$seqname} = $order++;
        }

        $input{$seqname} .= $aln_line;
  }

  my @aln=();
  my $N=0;

  foreach my $key (keys %input) {
        $input{$key}=~s/\./-/g;
        $aln[$order{$key}] = {name=>$key, seq=>$input{$key}};
        $N++
  }
  return [@aln];
}

######################################################################
#
# formatAln(\@aln ref-to-alignment, $format string)
#
# Formats alignment from list of hash format into various formats. The
# input needs at least the "seq" key with the sequences. The rest is
# optional
#
# \@aln ... alignment in list of hash format
# $format ... CLUSTAL,  FASTA, MAF
#
# Returns formatted alignment as string.
#
######################################################################


sub formatAln{

  my @aln=@{$_[0]};

  my $format=$_[1];

  $format=lc($format);

  my @alnSeqs=();
  my @alnNames=();

  my $counter=1;

  foreach my $row (@aln){

        my $name="seq$counter";
        $counter++;

        $name=$row->{name};

        my $start=$row->{start};
        my $end=$row->{end};
        my $strand=$row->{strand};

        my $pos='';

        if (defined $start and defined $end){
          $pos="/$start-$end";
          if (defined $strand){
                if ($strand eq '+'){
                  $name.=""; # Don't put _fwd
                } else {
                  $name.="_rev";
                }
          }
        }

        push @alnNames, "$name$pos";
        push @alnSeqs, $row->{seq};

  }


  my $output='';

  if ($format eq 'clustal'){

        $output="CLUSTAL W(1.81) multiple sequence alignment\n\n\n";
        my $maxName=0;

        foreach my $name (@alnNames){
          $maxName=($maxName<length($name))?length($name):$maxName;
        }

        for my $i (0..$#alnNames){
          my $buffer=" "x(($maxName+6)-length($alnNames[$i]));
          $alnNames[$i].=$buffer;
        }
        my $columnWidth=60;
        my $currPos=0;
        my $length=length($alnSeqs[0]);

        while ($currPos<$length){
          for my $i (0..$#alnNames){
                $output.=$alnNames[$i];
                $output.=substr($alnSeqs[$i],$currPos,$columnWidth);
                $output.="\n";
          }
          $output.="\n\n";
          $currPos+=$columnWidth;
        }
  } elsif ($format eq 'fasta'){
        foreach my $i (0..$#alnNames){
          my $name=$alnNames[$i];
          my $seq=$alnSeqs[$i];
          $seq=~ s/(.{60})/$1\n/g;
          $output.=">$name\n$seq\n";
        }

  } elsif ($format eq 'maf'){
        $output.="a score=0\n";
        foreach my $row (@aln){
          my $length=$row->{end}-$row->{start};
          $output.="s $row->{name} $row->{start} $length $row->{strand} $row->{fullLength} $row->{seq}\n";
        }
        $output.="\n";
  }
  return $output;

}


sub getNextAln{

  my $format=shift;
  my $out='';

  my $fh=shift;
  if (!defined $fh){
        $fh=*STDIN;
  }

  $format=uc($format);

  while (<$fh>){
        if ($format eq "MAF"){
          last if /^a/;
        }
        if ($format eq "CLUSTAL"){
          last if /CLUSTAL/;
        }

        $out.=$_;

        last if eof; #seems to be necessary if <> does not read from stdin
                 #but from real file given without "<" at the
                 #commandline
  }
  return $out;

}


sub parseAln{

  (my $string,my $format)=@_;

  $format=uc($format);

  return readMAF($string) if $format eq "MAF";

  return readClustal($string) if $format eq "CLUSTAL";

  return [];

}

######################################################################
#
# sliceAlnByColumn(\@aln ref-to-alignment, $start int, $end int)
#
# Returns slice of alignment specified by alignment column.
#
# \@aln ... alignment in list of hash format
# $start, $end ... slice to cut
#
# Returns reference to alignment in list of hash format. This is a new
# alignment, i.e. the input is not sliced in place
#
######################################################################

sub sliceAlnByColumn: prototype($$$) {

  my @aln=@{$_[0]};
  shift;
  (my $start, my $end)=@_;

  # correct ends without warning if outside of valid range
  $start=0 if ($start<0);
  $end=length($aln[0]->{seq}) if ($end > length($aln[0]->{seq}));

  #my @newAln=@aln;

  # make deep copy of list of hash
  my @newAln=();
  foreach (@aln){
        push @newAln,{%{$_}};
  }

  foreach my $i (0..$#newAln){

        if ((defined $newAln[$i]->{start}) and (defined $newAln[$i]->{start})){
          my $oldStart=$newAln[$i]->{start};
          my $oldEnd=$newAln[$i]->{end};
          $newAln[$i]->{start}=alnCol2genomePos($newAln[$i]->{seq},$oldStart,$start);
          $newAln[$i]->{end}=alnCol2genomePos($newAln[$i]->{seq},$oldStart,$end-1)+1;
        }

        $newAln[$i]->{seq}=substr($newAln[$i]->{seq},$start,$end-$start);

  }

  return([@newAln]);

}


######################################################################
#
# sliceAlnByPos(\@aln ref-to-alignment, $index int, $start int, $end int)
#
# Returns slice of alignment relative to the genomic position of one
# given sequence.
#
# \@aln ... alignment in list of hash format
# $index ... Sequence in the alignment that should be used as reference
# $start, $end ... genomic positions of slice to cut
#
# Returns reference to alignment in list of hash format. This is a new
# alignment, i.e. the input is not sliced in place
#
######################################################################

sub sliceAlnByPos{

  my @aln=@{$_[0]};
  shift;
  (my $index, my $start, my $end)=@_;

  #my @newAln=@aln;

  # make deep copy of list of hash
  my @newAln=();
  foreach (@aln){
        push @newAln,{%{$_}};
  }

  # print "Start:$start, end:$end\n";

  my ($colStart, $colEnd);

  if ($start<$aln[$index]->{start}){
        $colStart=0;
  } else {
        $colStart=genomePos2alnCol($aln[$index]->{seq},$aln[$index]->{start},$start);
  }

  if ($end>$aln[$index]->{end}){
        $colEnd=length($aln[$index]->{seq});
  } else {
        $colEnd=genomePos2alnCol($aln[$index]->{seq},$aln[$index]->{start},$end-1)+1;
  }


  foreach my $i (0..$#newAln){

        my $oldStart=$newAln[$i]->{start};
        my $oldEnd=$newAln[$i]->{end};

        if ($newAln[$i]->{end}!=0){
          $newAln[$i]->{start}=alnCol2genomePos($newAln[$i]->{seq},$oldStart,$colStart,'after');
          $newAln[$i]->{end}=alnCol2genomePos($newAln[$i]->{seq},$oldStart,$colEnd-1,'before')+1;
        } else {
          $newAln[$i]->{start}=0;
          $newAln[$i]->{end}=0;
        }

        $newAln[$i]->{seq}=substr($newAln[$i]->{seq},$colStart,$colEnd-$colStart);

  }
  return([@newAln]);

}

######################################################################
#
# alnCol2genomePos($seq string, $start int, $col int)
#
# Calculates the genomic position corresponding to a column in an
# alignment.
#
# $seq ... sequence from alignment (i.e. letters with gaps)
# $start ... Genomic position of first letter in $seq
# $col ... column in the alignment that is to be mapped
# $gapDecision ... 'before' or 'after', if column maps to a gap, a decision
#                  has to be made whether it gets the positon of the
#                  letter before or after.
#
# Returns genomic position. No error handling, so $col must be a valid
# column of the string $seq.
#
#######################################################################


sub alnCol2genomePos{

  (my $seq, my $start, my $col, my $gapDecision)=@_;

  $gapDecision='after' if (!$gapDecision);

  $seq=~s/\./-/g; #Convert all gaps to "-"

  my $newPos=$start;

  # if gap only...
  return $start if ($seq=~/^-+$/);

  #print "$seq\n";
  (my $tmp)=$seq=~/(-*)[^-]/;

  my $leadingGaps=length($tmp);
  #print "$tmp,$leadingGaps\n";
  # if column is in gap before first letter,
  # return position of the first letter
  return $start if ($col<$leadingGaps);

  $newPos=$start-1;


  for my $i ($leadingGaps..$col){
        $newPos++ if ((my $t=substr($seq,$i,1)) ne '-');
  }

  if (substr($seq,$col,1) eq '-'){
        if ($gapDecision eq 'after'){
          $newPos++;
        }
  }
  return $newPos;
}

######################################################################
#
# genomePos2alnCol($seq string, $start int, $pos int)
#
# Calculates the column in an alignment corresponding to a genomic
# position.
#
# $seq ... sequence from alignment (i.e. letters with gaps)
# $start ... Genomic position of first letter in $seq
# $pos ... genomic position that is to be mapped
#
# Returns column in the alignment. No error handling, so $pos must be
# a valid position depending on $start and the length of $seq.
#
#######################################################################

sub genomePos2alnCol{

  (my $seq, my $start, my $pos)=@_;

  $seq=~s/\./-/g; #Convert all gaps to "-"

  my $newPos=$start;

  (my $tmp)=$seq=~/(-*)[^-]/;

  my $leadingGaps=length($tmp);

  $newPos=$start-1;

  for my $i ($leadingGaps..(length($seq)-1)){
        $newPos++ if (substr($seq,$i,1) ne '-');
        return $i if ($newPos==$pos);
  }
  return $newPos;
}


######################################################################
#
# revAln(\@aln ref-to-alignment)
#
# Returns reference to reverse complement of alignment (new alignment,
# i.e. not modified in place)
#
######################################################################

sub revAln{

  my @aln=@{$_[0]};

  # make deep copy of list of hash
  my @newAln=();

  foreach (@aln){
        push @newAln,{%{$_}};
  }

  foreach (@newAln){
        $_->{seq} = reverse $_->{seq};
        $_->{seq}=~tr/AGCTUagctu/TCGAAtcgaa/;

        if (defined $_->{strand}){
          if ($_->{strand} eq '+'){
                $_->{strand}='-';
          } else {
                $_->{strand}='+';
          }
        }
  }

  return [@newAln];

}

######################################################################
#
# rangeWarn(\@aln alnRef)
#
# Checks if alignment is within RNAz definition ranges. Return values:
#
# 0: everything is in range
# 1: sequence is longer than 400 or shorter than 50
# 2: base ratios are above 0.75 or below 0.25
# 3: length and base ratios out of range
#
######################################################################

sub rangeWarn{

  my @aln=@{$_[0]};

  my $lengthWarn=0;
  my $compWarn=0;

  foreach my $row (@aln){

        my $seq=$row->{seq};

        $seq=uc($seq);

        $seq=~s/-//g;
        $seq=~s/\.//g;

        my $length=length($seq);

        return 1 if ($length==0);

        my $g=($seq=~tr/G/G/);
        my $c=($seq=~tr/C/C/);
        my $a=($seq=~tr/A/A/);
        my $t=($seq=~tr/T/T/);
        my $u=($seq=~tr/U/U/);

        $t+=$u;

        my $GC=($g+$c)/$length;


        my $A=0;
        my $C=0;

        if (($t+$a)>0 and ($g+$c)>0){
          $A=$a/($t+$a);
          $C=$c/($g+$c);
        }
        $lengthWarn=1 if ($length<50 or $length>400);

        $compWarn=1 if ($GC<0.25 or $GC>0.75);
        $compWarn=1 if ($A<0.25 or $A>0.75);
        $compWarn=1 if ($C<0.25 or $C>0.75);

  }

  return 3 if ($compWarn==1 and $lengthWarn==1);
  return 1 if ($lengthWarn==1);
  return 2 if ($compWarn==1);

  return 0;

}

######################################################################
#
# removeCommonGaps(\@aln ref-to-alignment)
#
# Removes gap only columns in \@aln (in place!)
#
######################################################################


sub removeCommonGaps{

  my $alnRef=$_[0];

  my @aln=();

  foreach (@$alnRef){
        push @aln, [split(//,$_->{seq})];
  }

  my $currCol=0;

  while ($aln[0][$currCol]){

        my $allGap=1;
        foreach my $currRow (0..$#aln){
          if ($aln[$currRow][$currCol] ne '-'){
                $allGap=0;
                last;
          }
        }

        if ($allGap){
          foreach my $currRow (0..$#aln){
                splice(@{$aln[$currRow]},$currCol,1);
          }
        } else {
          $currCol++;
        }
  }

  foreach my $i (0..$#aln){
        $alnRef->[$i]->{seq}=join('',@{$aln[$i]});
  }

}

######################################################################
#
# meanPairID(\@aln ref-to-alignment)
#
# Returns mean pairwise identity of alignment
#
######################################################################


sub meanPairID{

  my @inputAln=@{$_[0]};

  my @aln=();

  foreach (@inputAln){
        push @aln, [split(//,$_->{seq})];
  }

  my $pairs=0;
  my $matches=0;

  for my $i (0..$#aln){
        for my $j ($i+1..$#aln){
          for my $k (0..(@{$aln[0]}-1)){
                if (($aln[$i][$k] ne '-') or ($aln[$j][$k] ne '-')){
                  if ($aln[$i][$k] eq $aln[$j][$k]){
                        $matches++;
                  }
                  $pairs++;
                }
          }
        }
  }
  return sprintf("%.4f",$matches/$pairs);
}



sub pruneAln{

  my %args = (maxN=>6,
                          minN=>2,
                          optSim => 0.8,
                          keepfirst=> 1,
                          maxID=> 0.95,
                          numAln=>1,
                          verbose=>0,
                          @_);

  my $maxN=$args{maxN};
  my $minN=$args{minN};
  my $optSim=$args{optSim};
  my $keepfirst=$args{keepfirst};
  my $maxID=$args{maxID};
  my $verbose=$args{verbose};
  my $numAln=$args{numAln};

  my @aln=@{$args{alnRef}};
  my $N=@aln;

  my @outAlns=();

  my $remove = sub {
        my $i = shift;
        warn "broken entry $i aln: %{$aln[$i]}" if !exists($aln[$i]->{name});
        $aln[$i]->{dead} = 1;
        print "removing $i:$aln[$i]->{name}\n" if $verbose;
  };

  my $alive =sub {
        my $i = shift;
        return $aln[$i]->{dead}?0:1;
  };

  # Calculate matrix of pairwise identities
  my @idMatrix=();
  foreach my $i (0..$N-1) {
        foreach my $j (0..$N-1) {
          next if $j==$i;
          my $seq1=$aln[$i]->{seq};
          my $seq2=$aln[$j]->{seq};

          $idMatrix[$i][$j]=$idMatrix[$j][$i]=
                meanPairID([{'seq'=>$seq1},{'seq'=>$seq2}]);
        }
  }

  if ($verbose>1) {
        foreach my $i (0..$#idMatrix-1) {
          print $aln[$i]->{name},":\n--------\n";
          foreach my $j ($i+1..@{$idMatrix[$i]}-1) {
                #print $idMatrix[$i][$j], "  ";
                print "$aln[$j]->{name}: $idMatrix[$i][$j]\n";
          }
          print "\n";
        }
  }

  my @used = (0) x @aln;

  for my $alnid (1..$numAln) {
        my $Nalive = $N;

        # Step 1: remove (almost) identical sequences
        foreach my $i (0..$N-1) {
          next unless &$alive($i);
          foreach my $j ($i+1..$N-1) {
                next unless &$alive($j);
                if ($idMatrix[$i][$j]>$maxID) {
                  print "$i,$j:$idMatrix[$i][$j] " if $verbose;
                  # remove one of the 2 seqs, prefer the one that's been used more often
                  my $r = (rand()*($used[$i]+$used[$j])<$used[$i]) ? $i : $j;
                  $r = $j if $keepfirst && ($i==0);
                  $Nalive--;
                  last if $Nalive<$minN;
                  &$remove($r);
                  last if $r==$i;
                }
          }
        }

        if ($alnid>1) {
          # Step 2: pre-select sequence to get several different samples.
          # We remove ($N-$maxN)/2 of the already used sequences.
          # choose sequences to be removed with probability
          # proportional to the number of times they've been used.
          my $s = 0;
          foreach (@used) {$s += $_};
          $s -= $used[0] if $keepfirst;
          print "going to preremove int(($Nalive-$maxN)/2) seqs\n"
                if $verbose;
          for my $i (1..int(($Nalive-$maxN)/2)) {
                my $r = rand($s);
                my $ss=0;
                for my $a ($keepfirst .. $#aln) {
                  next unless &$alive($a);
                  $ss += $used[$a];
                  if ($ss>$r) {
                        print "next aln " if $verbose;
                        &$remove($a); $Nalive--;
                        $s -= $used[$a];
                        last;
                  }
                }
          }
        }
        # Step 3: Optimize mean pairwise similarity (greedily)
        # remove worst sequence until desired number is reached
        while ($Nalive > $maxN) {
          my $maxcost=0;
          my $maxind;
          foreach my $i (0..$N-1) {
                next unless &$alive($i);
                next if $i==0 && $keepfirst; # never delete seq 0
                my $cost = 0;
                foreach my $j (0..$N-1) {
                  next if $i==$j;
                  $cost += ($idMatrix[$i][$j]-$optSim)*($idMatrix[$i][$j]-$optSim);
                }
                ($maxcost,$maxind) = ($cost,$i) if $cost>$maxcost;
          }
          &$remove($maxind); $Nalive--;
        }

        #  my @newaln = grep {!$_->{dead}} @aln;
        my @newaln;
        foreach my $row (@aln) {
          # need to do a deep copy here, else we'll modify the original aln
          push @newaln, {name=>$row->{name},
                                         seq=>$row->{seq},
                                         start=>$row->{start},
                                         end=>$row->{end},
                                         chrom=>$row->{chrom},
                                         org=>$row->{org},
                                         strand=>$row->{strand},
                                         fullLength=>$row->{fullLength}
                                        } unless $row->{dead};
        }

        removeCommonGaps(\@newaln);

        #print formatAln(\@newaln,"CLUSTAL");

        push @outAlns, [@newaln];

        for my $i (0..$#aln) {
          if (&$alive($i)) {
                $used[$i]++;
          } else {
                $aln[$i]->{dead} = 0; # revive
          }
        }
  }

  return [@outAlns];

}

sub getNextRNAz{

  my $fh=shift;
  if (!defined $fh){
        $fh=*STDIN;
  }

  my $out='';

  while (<$fh>){

        next if /^\s?$/;

        last if /^\#.*RNAz.*\#$/ and $out ne '';

        $out.=$_;

        last if eof; #seems to be necessary if <> does not read from stdin
                     #but from real file given without "<" at the
                     #commandline
  }
  return $out;
}

sub parseRNAz{

  my $rnaz=shift;

  my @rnaz=split(/^/, $rnaz);
  my ($N,$identity,$columns,$decValue,$P,$z,$sci,$energy,$strand,
      $covariance,$combPerPair,$meanMFE,$consensusMFE,$consensusSeq,
      $consensusFold, $GCcontent, $ShannonEntropy);

  my @aln=();

  foreach my $i (0..$#rnaz){
        my $line=$rnaz[$i];
        $identity=$1 if ($line=~/Mean pairwise identity:\s*(-?\d+.\d+)/);
        $N=$1 if ($line=~/Sequences:\s*(\d+)/);
        if ($line=~/Reading direction:\s*(forward|reverse)/){
          $strand=($1 eq 'forward')?'+':'-';
        }
        $columns=$1 if ($line=~/Columns:\s*(\d+)/);
        $decValue=$1 if ($line=~/SVM decision value:\s*(-?\d+.\d+)/);
        $P=$1 if ($line=~/SVM RNA-class probability:\s*(-?\d+.\d+)/);
        $z=$1 if ($line=~/Mean z-score:\s*(-?\d+.\d+)/);
        $sci=$1 if ($line=~/Structure conservation index:\s*(-?\d+.\d+)/);
        $energy=$1 if ($line=~/Energy contribution:\s*(-?\d+.\d+)/);
        $covariance=$1 if ($line=~/Covariance contribution:\s*(-?\d+.\d+)/);
        $combPerPair=$1 if ($line=~/Combinations\/Pair:\s*(-?\d+.\d+)/);
        $consensusMFE=$1 if ($line=~/Consensus MFE:\s*(-?\d+.\d+)/);
        $meanMFE=$1 if ($line=~/Mean single sequence MFE:\s*(-?\d+.\d+)/);
        $GCcontent=$1 if ($line=~/G\+C content:\s(\d+.\d+)/);
        $ShannonEntropy=$1 if ($line=~/Shannon entropy:\s*(\d+.\d+)/);

        if ($line=~/^>/){
          chomp($rnaz[$i+1]);
          chomp($rnaz[$i+2]);
          if ($line=~/^>consensus/){
                $consensusSeq=$rnaz[$i+1];
                $consensusFold=substr($rnaz[$i+2],0,length($rnaz[$i+1]));
                last;
          } else {

                if ($line=~/>(.*?) (\d+) (\d+) (\+|\-) (\d+)/){
                  push @aln, {name=>$1,
                                          start=>$2,
                                          end=>$2+$3,
                                          strand=>$4,
                                          fullLength=>$5,
                                          seq=>$rnaz[$i+1],
                                          fold=>substr($rnaz[$i+2],0,length($rnaz[$i+1]))};
                  $i+=2;
                } elsif ($line=~/^(.*)\/(\d+)-(\d+)$/){
                  push @aln, {name=>$1,
                                          start=>$2,
                                          end=>$3,
                                          strand=>$strand,
                                          fullLength=>'',
                                          seq=>$rnaz[$i+1],
                                          fold=>substr($rnaz[$i+2],0,length($rnaz[$i+1]))};
                  $i+=2;
                }
          }
        }
  }

  return {"N"=>$N,
          "identity"=>$identity,
          "decValue"=>$decValue,
          "columns"=>$columns,
          "P"=>$P,
          "z"=>$z,
          "sci"=>$sci,
          "energy"=>$energy,
          "covariance"=>$covariance,
          "combPerPair"=>$combPerPair,
          "consensusMFE"=>$consensusMFE,
          "meanMFE"=>$meanMFE,
          "consensusSeq"=>$consensusSeq,
          "consensusFold"=>$consensusFold,
          "refSeqName"=>$aln[0]->{name},
          "refSeqStart"=>$aln[0]->{start},
          "refSeqEnd"=>$aln[0]->{end},
          "refSeqStrand"=>$aln[0]->{strand},
          "aln"=>[@aln],
          "rawOutput"=>$rnaz,
          "GC"=>$GCcontent,
          "entropy"=>$ShannonEntropy};

}

######################################################################
#
# shuffleAln(\@aln ref-to-alignment, $level string) Shuffles the
#
# alignment. Corresponds to conservative2 mode of shuffle-aln.pl
#
# \@aln ... alignment in list of hash format
#
# Returns reference to alignment (new copy not shuffled in place)
#
######################################################################


sub shuffleAln{

  my @inputAln=@{$_[0]};
  my $level=$_[1];

  if (!defined $level){
        $level=2;
  }

  # convert in list of list format which is use by shuffle-aln.pl
  my @aln=();
  foreach my $line (@inputAln){
        $line->{seq}=~s/\./-/g;
        push @aln,[split(//,$line->{seq})];
  }

  #my @aln=@{$_[0]};
  my $maxRow=$#aln;
  my $maxCol=@{$aln[0]}-1;

  my @list=(); # stores mask for each column
  my %hash=(); # stores columns for each mask as hash of array
               # with mask as key

  # creates characteristic mask for each column and
  # writes @list and %hash
  foreach my $currCol (0..$maxCol){
        my %seen=();
        my $mask='';
        my $counter=0;

        # creates mask for gap-pattern a la: --XX---X
        foreach my $currRow (0..$maxRow){
          my $currNt=$aln[$currRow][$currCol];
          if ($currNt eq '-'){
                $mask.='-';
          } else {
                $mask.='X';
          }
        }
        # Calculates mean pairwise identity for each column
        my $pairs=0;
        my $matches=0;

        for my $i (0..$maxRow){
          for my $j ($i+1..$maxRow){

                my $charI=uc($aln[$i][$currCol]);
                my $charJ=uc($aln[$j][$currCol]);

                if (($charI ne '-') and ($charJ ne '-')){
                  if ($charI eq $charJ){
                        $matches++;
                  }
                  $pairs++;
                }
          }
        }

        # Adds mean pairwise identity to mask. The rounding in printf can
        # be used as a convenient measure for coarse graining.
        my $id;

        if ($pairs>0){
          if ($level==0){
                $id=int(($matches/$pairs)+0.5);
          }
          if ($level==1){
                $id=sprintf("%.1f",$matches/$pairs);
          }
          if ($level==2){
                $id=sprintf("%.2f",$matches/$pairs);
          }

        }
        else {
          $id=sprintf("%.0f",1);
        }

        $mask.=$id;

        push @list, $mask;
        if (!exists $hash{$mask}){
          $hash{$mask}=[$currCol];
        } else {
          push @{$hash{$mask}},$currCol;
        }
  }

  #print "$_\n" foreach (@list);

  #return ();

  # each list of columns with the same mask
  # are shuffled (Fisher-Yates, code from perlfaq4)
  foreach my $arrayRef (values %hash){
        my $i;
        for ($i = @$arrayRef; --$i; ) {
          my $j = int rand ($i+1);
          @$arrayRef[$i,$j] = @$arrayRef[$j,$i];
        }
  }

  # columns are reassembled to a shuffled alignment
  my @shuffledAln;
  foreach my $currCol (0..$maxCol){
        my $randomCol=shift @{$hash{$list[$currCol]}};
        foreach  my $currRow (0..$maxRow){
          $shuffledAln[$currRow][$currCol]=$aln[$currRow][$randomCol];
        }
  }

  my @output=();

  foreach my $i (0..$#shuffledAln){
        my %tmp=();
        $tmp{seq}=join('',@{$shuffledAln[$i]});
        $tmp{name}=$inputAln[$i]->{name};
        $tmp{org}=$inputAln[$i]->{org};
        $tmp{chrom}=$inputAln[$i]->{chrom};
        $tmp{start}=$inputAln[$i]->{start};
        $tmp{end}=$inputAln[$i]->{end};
        $tmp{strand}=$inputAln[$i]->{strand};
        push @output,{%tmp};

  }
  return \@output;
}

sub getSeq{

  (my $file,my $start, my $end, my $strand)=@_;

  $strand='+' unless $strand;

  open(my $SEQFILE, "<", "$file") || die("Could not read $file: $!");

  my $firstline = <$SEQFILE>;
  my $secondline = <$SEQFILE>;

  my $headerLength = length $firstline;
  my $defaultLength = length($secondline)-1;

  my $startAddress=$headerLength+int($start/$defaultLength)+$start;

  my $endAddress=$headerLength+(int($end/$defaultLength))+$end;

  my $seq;
  seek $SEQFILE, $startAddress,0;

  if ($startAddress>$endAddress){
        print "warning!\n";
        return "";
  }

  read $SEQFILE, $seq, ($endAddress-$startAddress);

  $seq =~ s/\n//g;

  close($SEQFILE);
  if ($strand eq "-"){
        $seq=reverse $seq;
        $seq=~tr/AGCTUagctu/TCGAAtcgaa/;
  }
  return $seq;
}
######################################################################
#
# blastSeq($dir string, $db string, $cutoff real,
#          $seq string, $blastExecutable string)
#
# Runs blast search for sequence $seq on the database $db with a
# e-value cutoff of $cutoff. $dir: BLASTDIR
#
# Returns a list of hash. See code for details. Returns empty list
# if no hit was found below the cutoff.
#
######################################################################

sub blastSeq{

  (my $dir, my $db, my $cutoff, my $seq, my $blastExecutable)=@_;

  open(my $TMP,">","/tmp/blast$$.fa");

  print $TMP ">dummy\n$seq\n";

  close($TMP);

  my @results=`$blastExecutable -p blastn -e $cutoff -d $db -m 8 -i /tmp/blast$$.fa`;

  #print @results;

  return () if (!@results);

  my @output=();

  foreach my $line (@results){

        my %tmp=();

        (my $queryID, my $subjectID, my $identity,
         my $length, my $mismatches, my $gaps, my $queryStart,
         my $queryEnd, my $subjectStart, my $subjectEnd,
         my $e, my $bit)=split(/\s/,$line);

        $tmp{queryID}=$queryID;
        $tmp{subjectID}=$subjectID;
        $tmp{identity}=$identity;
        $tmp{length}=$length;
        $tmp{mismatches}=$mismatches;
        $tmp{gaps}=$gaps;
        $tmp{queryStart}=$queryStart;
        $tmp{queryEnd}=$queryEnd;
        $tmp{subjectStart}=$subjectStart;
        $tmp{subjectEnd}=$subjectEnd;
        $tmp{e}=$e;
        $tmp{bit}=$bit;

        push @output,{%tmp};

  }
  return @output;
}

######################################################################
#
# niceNumber(pos int)
#
# Makes 26,456,345 out of 26456345
#
######################################################################

sub niceNumber{
  my $pos=shift;
  my $rev=reverse $pos;
  $rev=~s/(...)/$1,/g;
  $pos=reverse $rev;
  $pos=~s/^,//;
  return $pos;
}


1;
