#!/usr/bin/perl -w

=head1 NAME

realing-all

=head1 SYNOPSIS

realign-all [options] <annotation-file>

Options:

=over 1

=item  B<--help>                        Brief help message

=item  B<--man>                         Full documentation

=item  B<-v, --verbose>                 Verbose

=item  B<-q, --quiet>                   Quiet

=item  B<--test>                        Test only. Jobs are not submitted!

=item  B<--revcompl>                    Realign reverse complement

=back

=head1 DESCRIPTION

calls mlocarna on sequence sets in Realign-Sequences, writes result
files to Alignment-Results, takes alignment jobs from annotation file

Distribute jobs to SGE-cluster. Assume script is run on a submission host!

=cut


use strict;


## ----------------------------------------
##
#constants

my $homedir = readpipe("pwd");
chomp $homedir;

my $source_dir=$homedir."/Realign-Sequences";
my $tgt_dir=$homedir."/Alignment-Results";
my $tmp_dir=$tgt_dir; # specify a directory with fast access from the
		      # compute server here,
		      # e.g. "/scratch/0/$HOME/Alignment-Results";

my $mlocarna="/home/will/Soft/locarna-1.4.8/bin/mlocarna";
my $locarna_options="--probabilistic --consistency-transformation --max-diff-match=100 --temperature=150 --struct-weight=200 --mea-beta 400 --moreverb";

##------------------------------------------------------------
## options
use Getopt::Long;
use Pod::Usage;

my $help;
my $man;
my $quiet;
my $verbose;
my $test;

my $revcompl;

## Getopt::Long::Configure("no_ignore_case");

GetOptions(	   
    "verbose" => \$verbose,
    "quiet" => \$quiet,   
    "help"=> \$help,
    "man" => \$man,
    "test" => \$test,
    "revcompl" => \$revcompl
    ) || pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

if ($#ARGV!=0) {print STDERR "Need annotation file.\n"; pod2usage(-exitstatus => -1);}
 
my $joblist="$ARGV[0]";

## ------------------------------------------------------------
## main part

if (!-d $tgt_dir) {
    print STDERR "$tgt_dir does not exist. Exit.\n";
    exit -1;
}
if (!-d $tmp_dir) {
    print STDERR "$tmp_dir does not exist. Exit.\n";
    exit -1;
}

## write sge job-script

my $tmpjoblist="$tmp_dir/realign.joblist";
my $jobscript="$tmp_dir/realign.sge";

system "cp",$joblist,$tmpjoblist || die "Cannot read and copy joblist $joblist";

my $num_tasks=`wc -l $tmpjoblist | cut -f1 -d' '`;
chomp $num_tasks;


my $rcsuf=""; ## suffix string in case of reverse complement alignment
if ($revcompl) {$rcsuf="-rc";}


open(OUT,">$jobscript") || die "Cannot write jobscript $jobscript.";
print OUT "#!/bin/bash
#\$ -e $tmp_dir/stderr
#\$ -o $tmp_dir/stdout
#\$ -l h_vmem=4G

line=`cat $tmpjoblist | head -\$SGE_TASK_ID | tail -1`;
locus=`echo \$line | cut -f1 -d' '`;
id=`echo \$line | cut -f5 -d' '`;
name=\$locus-\$id$rcsuf

$mlocarna $source_dir/\$name.mfa $locarna_options --tgtdir $tgt_dir/\$name.dir  > $tmp_dir/\$name.output 2>&1;

if [ $tmp_dir != $tgt_dir ] ; then cp $tmp_dir/\$name.output $tgt_dir; fi
";


close OUT;


my $submission_cmd = "qsub -t 1-$num_tasks $jobscript";
print "EXEC: $submission_cmd\n";
system($submission_cmd) unless $test;

## ------------------------------------------------------------
