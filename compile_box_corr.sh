#!/bin/bash
set -e
GSL_INC="${GSL_INC:-/usr/include}"
GSL_LIB="${GSL_LIB:-/usr/lib64}"

gcc -O2 -fopenmp \
    src/box_corr_main.c \
    src/interface_state.c \
    src/spin_init.c \
    -o BoxCorr \
    -I src/ -I"${GSL_INC}" \
    -L"${GSL_LIB}" -lgsl -lgslcblas -lm \
    -Wl,-rpath,"${GSL_LIB}"

echo "Build succeeded: ./BoxCorr"
