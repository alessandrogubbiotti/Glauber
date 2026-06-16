#!/bin/bash
set -e

# Adjust GSL paths if installed locally (e.g. ~/.local)
GSL_INC="${GSL_INC:-/usr/include}"
GSL_LIB="${GSL_LIB:-/usr/lib64}"

echo "Compiling Gillespie interface simulation..."

gcc -O2 -g -fopenmp \
    src/gillespie_main.c \
    src/interface_state.c \
    src/interface_stats.c \
    src/spin_init.c \
    src/configuration.c \
    src/configuration_statistics.c \
    -o Gillespie \
    -I src/ \
    -I"${GSL_INC}" \
    -L"${GSL_LIB}" \
    -lgsl -lgslcblas -lm \
    -Wl,-rpath,"${GSL_LIB}"

echo "Build succeeded: ./Gillespie"
