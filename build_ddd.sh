#! /bin/bash

git clone https://github.com/lip6/libDDD.git

cd libDDD

./configure --prefix=$PWD/../

make -j2 install




