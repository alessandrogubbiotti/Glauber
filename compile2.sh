#!/bin/bash
set -e

# Glauber_tidy — the first-generation, discrete-time parity-update engine.
# Kept as the INDEPENDENT REFERENCE implementation the Gillespie engine was
# validated against (different algorithm, different error class), and as the
# kymograph generator.  For production runs use ./compile_gillespie.sh.
#
# Reproducibility: export GLAUBER_SEED=<n> to repeat a run exactly; otherwise
# each run draws a fresh seed (printed at startup).

# Adjust GSL paths if installed locally (e.g. ~/.local, or Homebrew on macOS:
# GSL_INC=/opt/homebrew/include GSL_LIB=/opt/homebrew/lib ./compile2.sh)
GSL_INC="${GSL_INC:-/usr/include}"
GSL_LIB="${GSL_LIB:-/usr/lib64}"

CC="${CC:-gcc}"

echo "Compiling Glauber_tidy (reference engine)..."

# Glauber_tidy.c defines its own initialize_spins_* copies, so spin_init.c is
# deliberately NOT linked here (it would be a duplicate-symbol error).
"${CC}" -O2 -g \
    src/Glauber_tidy.c \
    src/configuration.c \
    src/configuration_statistics.c \
    -o Glauber \
    -I src/ \
    -I"${GSL_INC}" \
    -L"${GSL_LIB}" \
    -lgsl -lgslcblas -lm \
    -Wl,-rpath,"${GSL_LIB}"

echo "Build succeeded: ./Glauber"
