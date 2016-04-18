#!/usr/bin/env bash

##
## build a new source package of locarna and test it
##
##

if [ "$1" = "--dput" ] ; then
  shift
  INSTALL_IT=true
fi

## release and package version have to be used consistently
## in distribution tar file and debian/changelog
PACKAGE_NAME=locarna
RELEASE=1.8.8
PACKAGE_VERSION=1

## distribution and architecture are relevant only for the creation of
## the binary package
DISTRIBUTION=${1:-wily}
ARCHITECTURE=amd64

PACKAGE_FULLVERSION="$RELEASE-$PACKAGE_VERSION~${DISTRIBUTION}"
PACKAGE_BASENAME="${PACKAGE_NAME}_${PACKAGE_FULLVERSION}"

CHANGELOG_HEADER="${PACKAGE_NAME} (${PACKAGE_FULLVERSION}) ${DISTRIBUTION}; urgency=low"


function runIt {
  echo $*
  if ! $* ; then
      echo "FAILED"
      exit -1
  fi
}

## pre: make distcheck -j -C _build

## copy the distribution tar file from _build 
runIt cp ../_build/locarna-$RELEASE.tar.gz locarna_$RELEASE.orig.tar.gz

## extract the original tar file
runIt tar xzf ../_build/locarna-$RELEASE.tar.gz

## modify changelog header
runIt cp debian/changelog debian/changelog.orig

echo "Set changelog header to ${CHANGELOG_HEADER}" 
(echo $CHANGELOG_HEADER; tail -n +2 debian/changelog.orig ) > debian/changelog

## call editor to edit changelog
## runIt emacs debian/changelog 

## copy debian directory to release directory
runIt cp -Lr debian locarna-$RELEASE

## create source package (signed).
## 
runIt cd locarna-$RELEASE
runIt dpkg-buildpackage -S

## test package building with pbuilder-dist;
## pbuilder creates 
runIt cd ..
pbuilder_basefile=$HOME/pbuilder/${DISTRIBUTION}-base.tgz
if [ -e "$pbuilder_basefile" ] ; then
  runIt pbuilder-dist $DISTRIBUTION $ARCHITECTURE build ${PACKAGE_BASENAME}.dsc
else
  echo "Do not test build package for distribution ${DISTRIBUTION}, since base file ${pbuilder_basefile} is missing."
fi

# upload to launchpad
echo "The package was build successfully."
echo ""
dput_cmd=dput ppa:swill/locarna ${PACKAGE_BASENAME}_source.changes

if [ "$INSTALL_IT" = "true" ] ; then
    runIt $dput_cmd
else
    echo "Upload to launchpad by running"
    echo "  $dput_cmd"
fi

