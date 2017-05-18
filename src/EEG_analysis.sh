#!/bin/bash

CHANNEL=20
echo $CHANNEL

rm *.log script.gplot
./hellomake $CHANNEL /Users/vinay/Desktop/A10_C1.bdf
gnuplot script.gplot


