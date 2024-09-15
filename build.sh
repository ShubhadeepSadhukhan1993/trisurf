#!/bin/bash
#This is a build script that helps with the initial build as described in README file
libtoolize --force
aclocal
autoheader
automake -ac
autoconf
./configure --prefix=/home/shubsad/shubsad_user
make clean
make -j8
make install -j8
echo "Warning!"
echo "Python scripts for running the trisurf-ng have been moved to separate package"

