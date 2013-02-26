#!/bin/bash

## USAGE: mlocarna-svmgdnspdk.sh fastafilename [options to mlocarna]

## call mlocarna, where similarity matrix is computed by svmsgdnspdk
## (using fasta2shrep_gspan)
## 

FASTA2SHREP=/usr/local/user/RNAtools/fasta2shrep_gspan.pl
SVMSGDNSPDK=/home/costa/Projects/svmsgdnspdk/svmsgdnspdk
MLOCARNA=/usr/local/user/locarna-1.7.2/bin/mlocarna

FASTA_FILENAME=$1
FASTA2SHREP_PARAMS="-wins 200 -shift 100 -stack -t 3 -M 3" ## -stack: extra stack vertices -t grammar abstraction level (of RNAshape) ## -M number of shapes

echo "Run $FASTA2SHREP $FASTA2SHREP_PARAMS -fasta $FASTA_FILENAME"
rm -rf GSPAN/
time $FASTA2SHREP $FASTA2SHREP_PARAMS -fasta $FASTA_FILENAME ## writes output directory GSPAN
bzcat GSPAN/*bz2 > $FASTA_FILENAME.gspan

#gspan format defines graphs

SVMSGDNSPDK_RADIUS=2
SVMSGDNSPDK_DISTANCE=4

echo "Run $SVMSGDNSPDK -d $FASTA_FILENAME.gspan -R $SVMSGDNSPDK_RADIUS -D $SVMSGDNSPDK_DISTANCE -gt DIRECTED -a MATRIX"
$SVMSGDNSPDK -d $FASTA_FILENAME.gspan -R $SVMSGDNSPDK_RADIUS -D $SVMSGDNSPDK_DISTANCE -gt DIRECTED -a MATRIX ##-R = radius # -D = distance # -gt = graph type # -a = action

MATRIX_FILENAME="matrix_R${SVMSGDNSPDK_RADIUS}D${SVMSGDNSPDK_DISTANCE}.mtx" # changed in new version to ~  $FASTA_FILENAME_R${SVMSGDNSPDK_RADIUS}D${SVMSGDNSPDK_DISTANCE}.mtx 

## the produced matrix file has entries in [0,1] => convert it
cat $MATRIX_FILENAME | awk -v FACTOR=10000 '{for(i=1;i<=NF;i++){if (i==NR){printf("0 ")} else{ printf("%d ",$(i)*FACTOR)}} print ""}' > $MATRIX_FILENAME.converted

echo "Run $MLOCARNA --similarity-matrix $MATRIX_FILENAME.converted $*"
$MLOCARNA --similarity-matrix $MATRIX_FILENAME.converted $*
