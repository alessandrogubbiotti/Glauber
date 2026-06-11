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

/* Exact invariant measure of the reversible interface dynamics:
 * all l bonds carry an interface independently with probability p,
 * conditioned on an even total number of interfaces (ring constraint).
 * Sampled by rejection (acceptance ~ (1 + (1-2p)^l)/2 ~ 1/2), which is
 * exact and translation invariant, unlike the sequential construction
 * of initialize_spins_1 whose seam bond is forced by parity. */
void initialize_spins_6(int *spins, int l, void *params) {
    double *p = (double *)params;
    int parity;
    do {
        parity = 0;
        spins[0] = (rand() % 2) * 2 - 1;
        int bond;
        for (int i = 1; i < l; i++) {
            bond = (rand() / (double)RAND_MAX) < *p;
            spins[i] = bond ? -spins[i - 1] : spins[i - 1];
            parity ^= bond;
        }
        /* seam bond between site l-1 and site 0, sampled like the others */
        bond = (rand() / (double)RAND_MAX) < *p;
        parity ^= bond;
        /* accept iff the independently sampled seam bond agrees with the
         * parity closure imposed by the ring (spins[l-1] vs spins[0]) */
    } while (parity); /* odd total => seam contradicts ring => resample */
}
