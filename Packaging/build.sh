#!/usr/bin/env bash

##
## build a new source package of locarna and test it
##
##

## release and package version have to be used consistently
## in distribution tar file and debian/changelog
RELEASE=1.8.7
PACKAGE_VERSION=1

## distribution and architecture are relevant only for the creation of
## the binary package
DISTRIBUTION=wily
ARCHITECTURE=amd64

## pre: make distcheck -j -C _build

## copy the distribution tar file from _build 
cp ../_build/locarna-$RELEASE.tar.gz locarna_$RELEASE.orig.tar.gz

## extract the original tar file
tar xzf ../_build/locarna-$RELEASE.tar.gz

## call editor to edit changelog
emacs debian/changelog 

## copy debian directory to release directory
cp -r debian locarna-$RELEASE

## create source package (signed).
## 
cd locarna-$RELEASE
dpkg-buildpackage -S

## test package building with pbuilder-dist;
## pbuilder creates 
cd ..
pbuilder-dist $DISTRIBUTION $ARCHITECTURE build locarna_$RELEASE-$PACKAGE_VERSION.dsc
