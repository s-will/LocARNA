DISTRIBUTIONS="precise trusty xenial"
#wily 

# assume that locarna is installed to ../_inst
LOCARNA_BIN="../_inst/bin/locarna"
# so that we can get the version by 
RELEASE_VERSION=`$LOCARNA_BIN --version | cut -d' ' -f2`

PACKAGE_VERSION=${1:-0}
ARCHITECTURE=${2:-amd64}

for distrib in $DISTRIBUTIONS ; do
    ./build.sh --dput $RELEASE_VERSION $distrib $PACKAGE_VERSION $ARCHITECTURE
done
