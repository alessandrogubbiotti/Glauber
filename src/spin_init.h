#ifndef SPIN_INIT_H
#define SPIN_INIT_H

#include <gsl/gsl_rng.h>

/* Initial configurations drawn with a GSL generator (defined in spin_init.c).
 * The rand()-based ones, initialize_spins_1 ... initialize_spins_6, are the
 * ones chosen by init_function in config.txt and are declared in
 * configuration.h. */

/* The equilibrium (reversible) measure of the interface dynamics: an interface
 * on each of the l bonds independently with probability p, conditioned on an
 * even number of interfaces.  Pass p = nu/(N+nu). */
void sample_equilibrium_interfaces(int *spins, int l, double p, gsl_rng *r);

/* Independent spins, +1 or -1 with probability 1/2. */
void sample_disordered_spins(int *spins, int l, gsl_rng *r);

#endif /* SPIN_INIT_H */
