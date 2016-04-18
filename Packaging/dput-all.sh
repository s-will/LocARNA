DISTRIBUTIONS=wily precise trusty xenial  

for distrib in $DISTRIBUTIONS ; do
  ./build.sh --dput $distrib
fi
