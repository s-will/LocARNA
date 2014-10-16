#!/bin/bash

cd Tests

bin=bin
## perform pseudo installation in subdir $bin

if [ -e $bin ] ; then rm -rf $bin ; fi
mkdir $bin

topdir=../$srcdir/.. # top level directory of the locarna source

cp -l $topdir/Utils/mlocarna $bin
cp -l ../locarna.bin $bin/locarna
cp -l ../locarna_p $bin
cp -l ../sparse $bin
cp -l ../locarna_rnafold_pp $bin
ln -sf $topdir/lib . # also put perl lib modules in the right place

function calltest {
    echo "============================================================"
    echo CALL $*
    echo
    if $* ; then
      echo "======================================== OK"
    else
      echo "======================================== FAILED"
      exit -1
    fi
}

## run mlocarna
calltest $bin/mlocarna $topdir/Examples/archaea.fa --tgtdir test.out --alifold-cons -p 0.05 --max-diff 10
rm -rf test.out

calltest $bin/mlocarna $topdir/Examples/haca.snoRNA.fa --tgtdir test.out -p 0.05 --max-diff 20
rm -rf test.out

calltest $bin/mlocarna $topdir/Examples/archaea.fa --tgtdir test.out --pw-aligner bin/sparse -p 0.05  --max-diff 10
rm -rf test.out

calltest $bin/mlocarna $topdir/Examples/archaea.fa --tgtdir test.out --probabilistic --consistency-transformation -p 0.05  --max-diff 10
rm -rf test.out

## cleanup
rm -rf $bin
rm -f lib
