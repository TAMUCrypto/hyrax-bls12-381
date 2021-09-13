#!/bin/bash

set -x

rm ../output/poly_time.txt

for (( c=15; c<=21; c++ ))
do
   /usr/bin/time -v ../cmake-build-release/src/hyrax_time $c >> ../output/poly_time.txt
done