#include "configuration.h"
#include "spin_init.h"
#include <stdlib.h>
#include <gsl/gsl_randist.h>

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

/* Equilibrium (reversible) measure of the interface dynamics: each of the l
 * bonds carries an interface independently with probability p (pass
 * p = nu/(N+nu)), CONDITIONED on an even total number of interfaces (the ring
 * constraint).  Sampled exactly by rejection (acceptance ~ 1/2), then integrated
 * to spins with sigma_0 = +/-1 uniform.
 *
 * Rejection is required: this conditioned-even Bernoulli(p) field is the exact
 * Ising-ring measure (pair correlator eta^|n|).  Merely flipping one uniformly
 * chosen bond when the parity is odd produces a DIFFERENT measure, biased at
 * O(1/l) (~5% of the correlator at l=48), which is not stationary for the
 * dynamics and does not match the analytic formulas.
 *
 * rand()-based variant (used by Glauber_tidy via the init_function table); it
 * uses full-config rejection.  GSL-driven binaries should call
 * sample_equilibrium_interfaces(), which samples the same measure the faster
 * "particle" way (Binomial count + uniform placement). */
void initialize_spins_6(int *spins, int l, void *params) {
    double p = *(double *)params;
    int *bond = malloc(l * sizeof(int));
    int parity;
    do {
        parity = 0;
        for (int b = 0; b < l; b++) {
            bond[b] = (rand() / (double)RAND_MAX) < p;
            parity ^= bond[b];
        }
    } while (parity);                   /* reject odd: condition on even parity */
    spins[0] = (rand() % 2) * 2 - 1;
    for (int i = 1; i < l; i++)
        spins[i] = bond[i - 1] ? -spins[i - 1] : spins[i - 1];
    free(bond);
}

/* GSL-driven equilibrium sampler — see spin_init.h.
 *
 * "Particle" sampling (Prop 1.2), exact and faster than full-config rejection
 * for small p:
 *   1. number of interfaces k ~ Binomial(l, p) conditioned on k even;
 *   2. place k interfaces on a uniformly random k-subset of the l bonds
 *      (Floyd's algorithm — distinct bonds = the exclusion rule, O(k));
 *   3. integrate bonds -> spins, sigma_0 = +/-1 uniform.
 * This is the same conditioned-even Bernoulli(p) (Ising-ring) measure. */
void sample_equilibrium_interfaces(int *spins, int l, double p, gsl_rng *r) {
    unsigned int k;
    do { k = gsl_ran_binomial(r, p, (unsigned int)l); } while (k & 1u);

    int *bond = calloc((size_t)l, sizeof(int));
    if (!bond) { perror("sample_equilibrium_interfaces"); exit(EXIT_FAILURE); }
    /* Floyd's algorithm: uniform k-subset of {0,...,l-1} in O(k). */
    for (unsigned int i = (unsigned int)l - k; i < (unsigned int)l; i++) {
        unsigned long j = gsl_rng_uniform_int(r, i + 1);   /* 0..i */
        bond[bond[j] ? (int)i : (int)j] = 1;
    }

    spins[0] = (gsl_rng_uniform(r) < 0.5) ? 1 : -1;
    for (int i = 1; i < l; i++)
        spins[i] = bond[i - 1] ? -spins[i - 1] : spins[i - 1];
    free(bond);
}

/* GSL-driven disordered sampler — see spin_init.h. */
void sample_disordered_spins(int *spins, int l, gsl_rng *r) {
    for (int i = 0; i < l; i++)
        spins[i] = (gsl_rng_uniform(r) < 0.5) ? 1 : -1;
}
