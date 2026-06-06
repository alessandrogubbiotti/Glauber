#include "configuration.h"
#include <stdlib.h>

/* Initialization functions shared between Glauber_tidy and Gillespie binaries.
 * Declared in configuration.h (forward declarations). */

void initialize_spins_1(int *spins, int l, void *params) {
    spins[0] = (rand() % 2) * 2 - 1;
    double *p = (double *)params;
    for (int i = 1; i < l; i++)
        spins[i] = ((rand() / (double)RAND_MAX) < *p) ? -spins[i - 1] : spins[i - 1];
}

void initialize_spins_2(int *spins, int l, void *params) {
    (void)params;
    for (int i = 0;     i < l;     i++) spins[i] =  1;
    for (int i = l / 2; i < l;    i++) spins[i] = -1;
}

void initialize_spins_3(int *spins, int l, void *params) {
    (void)params;
    for (int i = 0; i < l; i++) spins[i] = -1;
}

void initialize_spins_4(int *spins, int l, void *params) {
    double *p = (double *)params;
    for (int i = 0; i < l; i++)
        spins[i] = ((rand() / (double)RAND_MAX) < *p) ? 1 : -1;
}

void initialize_spins_5(int *spins, int l, void *params) {
    (void)params;
    for (int i = 0; i < l; i++)
        spins[i] = (i % 2 == 0) ? 1 : -1;
}
