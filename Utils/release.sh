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
    if [ "$LOCAL" != true ] ; then
        echo "ERROR: work directory differs from repository. Please checkin."
        exit -1
    else
        echo "WARNING: work directory differs from repository. Please checkin."
    fi
fi

########################################
## check distribution and run tests
#
if [ "$CHECK" = "true" ] ; then
    echo "**************************************************"
    echo "**                                              **"
    echo "** Run make distcheck with g++/stdlibc++        **"
    echo "**                                              **"
    echo "**************************************************"
    if ! make -j -C $BUILDDIR distcheck ; then
        echo "ERROR: Make dist check failed (g++/libstdc++)."
        exit -1
    fi

    if [ "$CHECK_CLANG" = "true" ] ; then
        ##repeat dist check with clang / libc++
        echo "**************************************************"
        echo "**                                              **"
        echo "** Run make distcheck with clang++/libc++       **"
        echo "**                                              **"
        echo "**************************************************"
        if ! make -j -C $BUILDDIR distcheck \
             CC=/usr/bin/clang \
             CXX=/usr/bin/clang++ \
             CXXFLAGS="-g -O2 --stdlib=libc++" ; then
            echo "ERROR: Make dist check failed (only clang/libc++)."
            exit -1
        fi
    fi
else
    if ! make -j -C $BUILDDIR dist-gzip ; then
        echo "ERROR: Make dist-gzip failed."
        exit -1
    fi
fi

echo "Extract tar-ball."
if ! ( cd $BUILDDIR && tar xzf $PACKAGE-$RELEASE.tar.gz ) ; then
    echo "ERROR: Tar-ball $PACKAGE-$RELEASE.tar.gz cannot be extracted."
    exit -1
fi

########################################
## check documentation
#
if [ "$DOXYGEN" = "true" ] ; then
    if cd $BUILDDIR/$PACKAGE-$RELEASE &&\
	    ./configure && make doxygen-doc
    then
        echo "Doxygen documentation successfully created."
    else
        echo "ERROR: Doxygen documentation cannot be created."
        exit -1
    fi
fi

########################################
## tag and publish to central repository
#
if [ "$LOCAL" != "true" ] ; then
    runcmd hg tag -f "$PACKAGE-$RELEASE"
    runcmd hg push
else
    echo "Local mode: don't tag and push"
fi

########################################
## copy to web directory, modify web page, checkin
#

runcmd cp $BUILDDIR/$PACKAGE-$RELEASE.tar.gz $WEBWORKDIR/Releases
if cd $WEBWORKDIR/Releases ; then
    if [ "$LOCAL" != "true" ] ; then
        runcmd  cvs add $PACKAGE-$RELEASE.tar.gz
    else
        echo "Local mode: don't cvs add released tar.gz"
    fi
    cd ..
else
    print "ERROR: could not copy the release to the web directory."
    exit -1
fi

## copy ChangeLog
runcmd \cp $BUILDDIR/$PACKAGE-$RELEASE/ChangeLog $WEBWORKDIR

if cd $WEBWORKDIR ; then
    runcmdto index.html.new perl $BINDIR/update_index.pl $RELEASE index.html
    runcmd \cp index.html index.html.bkp
    runcmd \mv index.html.new index.html 
    if [ "$LOCAL" != "true" ] ; then
        runcmd cvs commit -m "\"Release $PACKAGE-$RELEASE\""
    else
        echo "Local mode: don't cvs commit web page"
    fi
else
    echo "ERROR: cannot change to web working directory $WEBWORKDIR"
    exit -1;
fi

if [ "$DOXYGEN" = "true" ] ; then
    if [ "$LOCAL" != "true" ] ; then
        ## copy documentation to webserver
        runcmd rsync -av $BUILDDIR/$PACKAGE-$RELEASE/Doc/ $WEBSERVER:$WEBSERVER_BASEDIR/$WEB_LOCARNA_HOME/Doc
    
    else
        echo "Local mode: rsync to local webdir"
        runcmd rsync -av $BUILDDIR/$PACKAGE-$RELEASE/Doc/ $WEBWORKDIR/Doc
    fi
fi

if [ "$LOCAL" != "true" ] ; then
    runcmd ssh $WEBSERVER cd $WEBSERVER_BASEDIR/$WEB_LOCARNA_HOME \; cvs update
    runcmd ssh $WEBSERVER chmod -R g+w $WEBSERVER_BASEDIR/$WEB_LOCARNA_HOME
else
    echo "Local mode: no cvs update on web server"
fi


echo "--------------------------------------------------"
echo "SUMMARY:"
echo "--------------------------------------------------"
if [ "$TEST" = "true" ] ; then
    echo "Dry run: nothing (outside of builddir) was changed."
else
    if [ "$DOXYGEN" != "true" ] ; then
        echo "Doxygen documentation was neither checked nor installed."
    fi
    if [ "$LOCAL" = "true" ] ; then
        echo "The release was published locally. NO web release; no commits to vc"
    else
        echo "The new release is published."
    fi
fi
echo
