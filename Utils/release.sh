#!/usr/bin/env bash

############################################################
##
## Script for automatic new release of LocARNA
##
## Runs various checks, copies release to web server,
## updated web page and online resources (release, doc).
##
##
## WARNING: ALPHA; still immature: perform non-test runs with care
##
## reads configuration from release.cfg
############################################################

BINDIR=$(dirname "$(readlink -fn "$0")")

. "$BINDIR/release.cfg"

function runcmd {
    if [ "$TEST" = "false" ] ; then
	echo $*
	$*
    else
	echo "DRY: $*"
    fi
}

function runcmdto {
    to=$1
    shift
    if [ "$TEST" = "false" ] ; then
	echo "$* \> $to"
	$* > $to
    else
	echo "DRY: $* \> $to"
    fi
}

cd $LOCARNA_HOME

########################################
## read and check arguments
RELEASE=$1

if [ -z "$RELEASE" ] ;then
    echo "USAGE: $0 <release>"
    echo "Please specify the release number."
    
    exit -1
fi

########################################
## check that release fits with configure.ac (AC_INIT)
if grep AC_INIT configure.ac | grep -F "AC_INIT([LocARNA], [$RELEASE]" ; then
    echo "Release number matches configure.ac"
else
    echo "ERROR: Release number differs from configure.ac"
    echo -n "  line ";
    grep -n AC_INIT $LOCARNA_HOME/configure.ac
    exit -1
fi


########################################
## check that ChangeLog contains release version number
if ! grep -w $RELEASE ChangeLog > /dev/null ; then
    echo "ERROR: no release entry in ChangeLog."
    exit -1
fi

########################################
## check that all files are checked in (clean status)
#
if [ "`hg -q status | wc -l`" != 0 ] ; then
    hg -q status
    echo "ERROR: work directory differs from repository. Please checkin."
    exit -1
fi

########################################
## check distribution and run tests
#
if !  make -j -C $BUILDDIR distcheck ; then
    echo "ERROR: Make dist check failed."
    exit -1
fi

########################################
## check documentation
#
if cd $BUILDDIR && tar xzf $PACKAGE-$RELEASE.tar.gz && cd $PACKAGE-$RELEASE &&\
	./configure && make doxygen-doc
then
    echo "Doxygen documentation successfully created."
else
    echo "ERROR: Doxygen documentation cannot be created."
    exit -1
fi

########################################
## tag and publish to central repository
#
runcmd hg tag -f "$PACKAGE-$RELEASE"
runcmd hg push

########################################
## copy to web directory, modify web page, checkin
#

runcmd cp $BUILDDIR/$PACKAGE-$RELEASE.tar.gz $WEBWORKDIR/Releases
if cd $WEBWORKDIR/Releases ; then
    runcmd  cvs add $PACKAGE-$RELEASE.tar.gz
    cd ..
else
    print "ERROR: could not copy the release to the web directory."
fi

## copy ChangeLog
runcmd \cp $BUILDDIR/$PACKAGE-$RELEASE/ChangeLog $WEBWORKDIR

if cd $WEBWORKDIR ; then
    runcmdto index.html.new perl $BINDIR/update_index.pl $RELEASE index.html
    runcmd \cp index.html index.html.bkp
    runcmd \mv index.html.new index.html 
    runcmd cvs commit -m \"Release $PACKAGE-$RELEASE\"
fi
       
## copy documentation to webserver
runcmd rsync -av $BUILDDIR/$PACKAGE-$RELEASE/Doc/ $WEBSERVER:$WEBSERVER_BASEDIR/$WEB_LOCARNA_HOME/Doc

runcmd ssh $WEBSERVER cd $WEBSERVER_BASEDIR/$WEB_LOCARNA_HOME \; cvs update

runcmd ssh $WEBSERVER chmod -R g+w $WEBSERVER_BASEDIR/$WEB_LOCARNA_HOME
