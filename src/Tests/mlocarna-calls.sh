#!/bin/bash

## Test mlocarna and locarnate

cd Tests

bin=bin
## perform pseudo installation in subdir $bin

if [ -e $bin ] ; then rm -rf $bin ; fi
mkdir $bin

topdir=../$srcdir/.. # top level directory of the locarna source

ln -sf ../$topdir/Utils/mlocarna $bin/mlocarna
ln -sf ../$topdir/Utils/locarnate $bin/locarnate
ln -sf ../../locarna.bin $bin/locarna
ln -sf ../../locarna_p $bin/locarna_p
ln -sf ../../sparse $bin/sparse
ln -sf ../../locarna_rnafold_pp $bin/locarna_rnafold_pp
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

## ========================================
## test mlocarna
##

calltest $bin/mlocarna $topdir/Examples/archaea.fa --tgtdir test.out --alifold-cons -p 0.05 --max-diff 10
rm -rf test.out

calltest $bin/mlocarna $topdir/Examples/haca.snoRNA.fa --tgtdir test.out -p 0.05 --max-diff 20
rm -rf test.out

calltest $bin/mlocarna $topdir/Examples/archaea.fa --tgtdir test.out --sparse -p 0.05  --max-diff 10
rm -rf test.out

calltest $bin/mlocarna $topdir/Examples/archaea.fa --tgtdir test.out --probabilistic --consistency-transformation -p 0.05  --max-diff 10
rm -rf test.out


## ========================================
## test locarnate (if t_coffee is available)
##
calltest $bin/locarnate $topdir/Examples/archaea.fa
if [ -d "test_results" ] ; then
    rm -rf test_results
fi
if [ -e "input.dnd" ] ; then
    rm -f input.dnd
fi

## cleanup
rm -rf $bin
rm -f lib
