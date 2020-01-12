#!/bin/bash
make clean
make

size=128
thdNum="8"
args="$thdNum 30 0.25 $size"
file=sample_"$size"x"$size"
echo file is $file.in $file.out

# for i in {0..5}
# do
    bash -c "./rainfall_pt $args $file.in > result.txt"
    bash -c "./check.py $size $file.out result.txt"
# done

