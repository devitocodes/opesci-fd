#!/bin/sh

pochoir $* &> /dev/null

# If it fails then run the compiler again and let the stdout/stderr go to screen.
if test $? != 0 ; then
  pochoir $*
fi

