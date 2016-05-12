DISTRIBUTIONS="precise trusty xenial"
#wily 

for distrib in $DISTRIBUTIONS ; do
    ./build.sh --dput $distrib
done
