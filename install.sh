#!/bin/sh

set -e

cd external/gmp
gmp_path="$PWD"
./configure --prefix="$gmp_path"
make -j
make install

cd ../ntl/src
./configure GMP_PREFIX="$gmp_path"
make -j

