#! /bin/bash

# Minato version
export SPOTVER=2.11.6.dev
# restrict/relax version
# export SPOTVER=2.10.4.dev


#wget --progress=dot:mega --no-check-certificate http://www.lrde.epita.fr/dload/spot/spot-$SPOTVER.tar.gz
wget --progress=dot:mega --no-check-certificate https://www.lrde.epita.fr/~adl/spot-$SPOTVER.tar.gz

tar zxf spot-$SPOTVER.tar.gz
rm spot-$SPOTVER.tar.gz

export IFOLDER=$(pwd)
cd spot-*

# From Etienne Renault : allow more acceptance
./configure -C VALGRIND=false --without-included-lbtt --disable-devel --disable-shared --prefix=$IFOLDER --enable-max-accsets=64

make -j2

# make check -j2
make install 



