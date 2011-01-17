#!/usr/bin/perl -w

=head1 NAME

compare_dotplots

=head1 SYNOPSIS

compare_dotplots [options] dotplot1 dotplot2

Options:

=over 1

=item  B<--help>                        Brief help message

=item  B<--man>                         Full documentation

=item  B<-v, --verbose>                 Verbose

=item  B<-q, --quiet>                   Quiet

=item  B<--compare-by>                  Measure for comparing column distributions (default=cor)

=back

=head1 DESCRIPTION

Compare two dot plots by average divergence/similarity of column probability distributions.

=cut


use strict;

##------------------------------------------------------------
## options
use Getopt::Long;
use Pod::Usage;

my $help;
my $man;
my $quiet;
my $verbose;

my $compare_by = "cor";

## Getopt::Long::Configure("no_ignore_case");

GetOptions(	   
    "verbose" => \$verbose,
    "quiet" => \$quiet,   
    "help"=> \$help,
    "man" => \$man,
    "compare-by=s" => \$compare_by
    ) || pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

sub parse_dp_ps {
    my ($filename) = @_;
    local *IN;
    
    open(IN,$filename) || die "Cannot read $filename for parsing as dp-ps file.\n";
    
    my $seq="";
    my %pairprobs;
    
    while (my $line=<IN>) {
	if ($line =~ /^\/sequence \{ \(/) {
	    while (defined($line = <IN>) && ($line !~  /\} def/  ))  {
		chomp $line;
		$line =~ s/\\$//;
		$seq .= $line;
	    }
	    #print "parse_dp_ps $filename: $seq\n";
	}
	
	if ($line =~ /(\d+)\s+(\d+)\s+(\S+)\s+ubox$/) {
	    $pairprobs{$1}{$2}=$3*$3;
	}
    }
    
    close IN;
    
    $seq ne "" || die "Empty sequence in dp.ps file $filename\n";
    
    return ($seq,\%pairprobs);
}

sub print_matrix($$$) {
    my ($file,$len,$matrix_ref) = @_;
    my %matrix = %{ $matrix_ref };

    open(OUT,">$file") || die "Cannot write tmp file $file.";
    
    for (my $i=1; $i<=$len; $i++) {
	for (my $j=1; $j<=$len; $j++) {
	    my $p=0;
	    if (exists ($matrix{$i}{$j})) {
		$p=$matrix{$i}{$j};
	    } elsif  (exists ($matrix{$j}{$i})) {
		$p=$matrix{$j}{$i};
	    }
	    print OUT "$p ";
	}
	print OUT "\n";
    }

    close OUT;
}


## ------------------------------------------------------------
## main part

($#ARGV == 1) || pod2usage(1);

my $dpfile1=$ARGV[0];
my $dpfile2=$ARGV[1];


my ($seq1,$bpprobs1_ref) = parse_dp_ps($dpfile1);
my ($seq2,$bpprobs2_ref) = parse_dp_ps($dpfile2);

my %bpprobs1=%{$bpprobs1_ref};
my %bpprobs2=%{$bpprobs2_ref};


my $len1=length($seq1);
my $len2=length($seq2);

my $tmp1="matrix1.$$.tmp";
my $tmp2="matrix2.$$.tmp";

print_matrix($tmp1,$len1,\%bpprobs1);
print_matrix($tmp2,$len2,\%bpprobs2);


my $rscript="gargs <- commandArgs()
t1<-read.table(gargs[5])
t2<-read.table(gargs[6])

measures=numeric(length(t1))

for (i in 1:length(t1)) { 
   x1 <- t1[[i]];
   x2 <- t2[[i]];
   x1[i] <- 1-sum(x1);
   x2[i] <- 1-sum(x2);
   measures[i]<-$compare_by(x1,x2);
}

#print(measures)
print(mean(measures))
";

system("printf '$rscript' | R -q --slave --args $tmp1 $tmp2|cut -f2 -d' '");

unlink $tmp1;
unlink $tmp2;

## ------------------------------------------------------------
