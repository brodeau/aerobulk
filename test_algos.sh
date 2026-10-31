#!/usr/bin/env bash

# FORT_FMT_RECL=512: wide standard output for the summary tables of the toy
# (needed with Intel and LLVM flang, which cut output lines at about 80 characters; harmless with gfortran)
FORT_FMT_RECL=512 ./bin/aerobulk_toy.x <<EOF
10.
2.
22.
20.
12.
5
EOF
