#!/bin/bash
set -e
GSL_INC="${GSL_INC:-/usr/include}"
GSL_LIB="${GSL_LIB:-/usr/lib64}"

gcc -O2 -fopenmp \
    src/corr2pt_main.c \
    src/interface_state.c \
    src/spin_init.c \
    src/configuration.c \
    src/configuration_statistics.c \
    -o Corr2pt \
    -I src/ -I"${GSL_INC}" \
    -L"${GSL_LIB}" -lgsl -lgslcblas -lm \
    -Wl,-rpath,"${GSL_LIB}"

echo "Build succeeded: ./Corr2pt"
