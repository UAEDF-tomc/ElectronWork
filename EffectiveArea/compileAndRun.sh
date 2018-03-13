#!/bin/bash
eval `scram runtime -sh`
echo "Compiling..."
echo
echo
g++ -O3 -I `root-config --incdir` -o $1 $1.C `root-config --libs` -lMinuit $ROOTSYS/lib/libTMVA.so -std=c++1y 2>&1 >/dev/null | head -n 35
./$1
