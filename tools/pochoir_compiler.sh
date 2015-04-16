#!/bin/sh

dir="."
for arg in $* ; do
  if test -f "$arg" ; then
    cp "$arg" .
    dir=`dirname "$arg"`
    gen=`basename "$arg" .cpp`.i
  fi
done
tmp=`mktemp`
pochoir $* &> $tmp

if test $? != 0 ; then
  cp $gen $dir/
  pochoir $* &> $tmp
  cat $tmp
fi
rm $tmp

