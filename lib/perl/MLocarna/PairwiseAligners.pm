package MLocarna::PairwiseAligners;


use MLocarna;
use strict;

############################################################
##
## Locarnate supports multiple Pairwise alignment tools
## all of which are defined here
##
############################################################

sub locarna_compute_pairwise_alignments {
    
  my ($options, $results_path, $sequences, $bindir, $pp_dir) = @_;

  my $numSequences = scalar(@{$sequences});

  # get locarna parameters
  my $parameter = $options->{pairwise_aligner_params} or "";

  my $counter = 0;
  # how many pairwise alignments must be calculated?
  my $number = $numSequences * 
      ($numSequences - 1) / 2;

  for (my $i = 0; $i < $numSequences - 1; ++$i) {
    # the path to the .pp file of the this sequence
    my $first = $pp_dir . '/S'.($i + 1) .= '.pp';

    # iterate over all other sequences in the rest of the array and compute
    # pairwise alignment
    for (my $j = $i + 1; $j < $numSequences; ++$j) {
      ++$counter;
      *STDERR->print("\r".(' 'x80));
      *STDERR->print("\r".'Calculate pairwise alignment with locarna '.
                     $counter.'/'.$number.' (Seq '.($i + 1).' <-> Seq '.
                     ($j + 1).')...');
      # path to the dot plot of the other sequence
      my $second = $pp_dir . '/S'.($j + 1) . '.pp';

      # path to the output file
      my $result = $results_path.'/pair/S'.($i + 1).'_S'.($j + 1).'.aln';

      # assemble the call to locarna
      my $call = $bindir .'/locarna '.$parameter.
          ' --clustal='.$result.' --write-structure '
          .$first.' '.$second.' 1>/dev/null 2>/dev/null';

      my $code = system($call);
      if ($code) {
        print('Can\'t perform locARNA calculation'.
                            " (developed under locarna 1.7.8, call:\n".$call.")\n");
      }
    }
  }
  
  *STDERR->print(' done!'."\n");

  *STDERR->print('parse pairwise alignments...'."\r");

  my @alignments = ();

  for (my $i = 0; $i < $numSequences - 1; ++$i) {
    for (my $j = $i + 1; $j < $numSequences; ++$j) {
      my $file_name = $results_path.'/pair/S'.($i + 1).
          '_S'.($j + 1).'.aln';
          
      my ($header, $aln) = MLocarna::read_clustalw_aln($file_name);
      my $score;
      if ($header =~ m/Score: ([\d|-]+)/) {
        $score = $1;
      }
      push(@alignments, {"rows" => $aln, "score" => $score});
    }
  }

  *STDERR->print('Parse pairwise alignments... done!'."\n");

  return(\@alignments);
}

sub exparna_compute_pairwise_alignments {
    my ($options, $results_path, $sequences, $bindir, $pp_dir) = @_;
    my $numSequences = scalar(@{$sequences});

    # get exparna parameters
    my $parameter = $options->{pairwise_aligner_parameter} or "";

    my $counter = 0;
    # how many pairwise alignments must be calculated?
    my $number = $numSequences * 
      ($numSequences - 1) / 2;

    # this has saves an alignment of 2 sequences with the key
    # "<firstSequence>-<secondSequence>"
    my %result = ();

    for (my $i = 0; $i < $numSequences - 1; ++$i) {
        # the path to the .pp file of the this sequence
        my $first = $pp_dir . '/S'.($i + 1) .= '.pp';

        # iterate over all other sequences in the rest of the array and compute
        # pairwise alignment
        for (my $j = $i + 1; $j < $numSequences; ++$j) {
            ++$counter;
            *STDERR->print("\r".(' 'x80));
            *STDERR->print("\r".'Calculate pairwise alignment with exparna '.
                         $counter.'/'.$number.' (Seq '.($i + 1).' <-> Seq '.
                         ($j + 1).')...');
            # path to the dot plot of the other sequence
            my $second = $pp_dir . '/S'.($j + 1) . '.pp';

            # path to the output file
            my $result = $results_path.'/pair/S'.($i + 1).'_S'.($j + 1).'.aln';

            # assemble the call to exparna
            my $call = $bindir ."/locarna_X $first $second " .
                "--number-of-EPMs 20 --suboptimal --inexact-struct-match --output-epm-list=" .
                $result . " 1>/dev/null 2>/dev/null";

            my $code = system($call);
            if ($code) {
                print('Can\'t perform Exparna calculation'.
                                " (developed under locarna 1.7.8, call:\n".$call.")\n");
            }

            $result{"$i-$j"} = parse_exparna_output($result);
        }
    }
    print "\n";
    return \%result;
}

# Takes an exparna output file and returns a list of baspairs with score
sub parse_exparna_output {
    my ($file) = @_;

    open(my $FILE, "<$file");

    my $header = <$FILE>;

    # this array holds hashes of the form {posA => int, posB => int, score => int}
    my  @result = ();

    # iterate over all EPMs
    while (my $line = <$FILE>) {
        # expects output like this: "<epm_id>\t<score>\t<structure>\t<list_of_pairs separated with :>"
        if ($line =~ /^\d+\t(\d+)\t.+\t(.+)$/) {
            my $score = $1;
            my @pairs = split /\s+/, $2;
            # iterate over all basepairs of this EPM
            foreach my $p (@pairs) {
                my ($posA, $posB) = split /:/, $p;

                push(@result, {"score" => $score, "posA" => $posA, "posB" => $posB});
            }
        }
    }
    print scalar(@result);
    close($FILE);
    return \@result;
}

sub exparna_result_to_tcoffee_lib_file {
    my ($filename, $alignments, $sequences) = @_;

    open(my $FILE, ">$filename");
    # print headers
    
    print $FILE "! TC_LIB_FORMAT_01\n";
    # print number of sequences
    print $FILE "" . (scalar(@{$sequences})) . "\n";

    # print sequences
    for(my $i = 0; $i < scalar(@{$sequences}); $i++) {
        print $FILE "S" . ($i + 1) . " " . length($sequences->[$i]->{seq}) . " " . $sequences->[$i]->{seq} . "\n";
    }

    # iterate over alignments
    while (my ($key, $value) = each(%{$alignments})) {
        my ($s1, $s2) = split /-/, $key;
        $s1++;
        $s2++;
        print "$key $value " . scalar(@{$value}) ."\n";
        print $FILE "#$s1 $s2\n";

        # iterate over all basepairs
        for(my $i = 0; $i < scalar(@{$value}); $i++) {
            print $FILE $value->[$i]->{posA} . " " . $value->[$i]->{posB} . " " . $value->[$i]->{score} . "\n";
        }
    }
    
    print $FILE "! SEQ_1_TO_N\n";
    close($FILE);
}

1;
