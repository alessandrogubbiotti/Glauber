#ifndef INTERFACE_STATE_H
#define INTERFACE_STATE_H

#include <gsl/gsl_rng.h>

/*
 * Interface-picture Gillespie simulation for 1D Glauber dynamics.
 *
 * An interface is a bond b where spin[b] != spin[(b+1)%l].
 * Each interface carries the spin value to its left (left_color).
 *
 * Domain j lies between interface j and interface (j+1)%k.
 * Width: w_j = (pos[j+1] - pos[j] + l) % l  (number of sites).
 *
 * Domain-centric event rates:
 *   w_j == 1  ->  annihilation at rate ann_rate
 *   w_j >= 2  ->  movement (left or right interface steps in) at rate 2
 * Global creation event at rate n_bulk * create_rate,
 * where n_bulk = sum_j max(0, w_j - 2).
 *
 * Macroscopic time = microscopic Gillespie time / N^2.
 */

typedef struct {
    int pos;        /* bond index on torus: 0 .. l-1         */
    int left_color; /* spin on the left of this bond: +1/-1  */
} Interface;

typedef struct {
    Interface *ifaces;    /* sorted array of interfaces, length k_alloc */
    int        k;         /* current interface count                     */
    int        k_alloc;   /* allocated capacity                          */
    int        l;         /* torus length = N * L                        */
    int        N;         /* sites per unit length (for N^2 rescaling)   */
    double     ann_rate;  /* annihilation rate per width-1 domain        */
    double     create_rate; /* creation rate per bulk spin               */
    int        n_bulk;    /* sum_j max(0, w_j - 2)                       */
    int        uniform_color; /* spin color of the ring when k == 0      */
    double     micro_time;  /* accumulated microscopic time              */
    gsl_rng   *rng;

    /* Event queue (k domain events + 1 global creation event).
     * domain_tau[j] is the ABSOLUTE microscopic time of the next event
     * for domain j.  create_tau is the absolute time for global creation.
     * INF (= 1e300) means the event does not exist in the current state. */
    double    *domain_tau;
    double     create_tau;

    long n_events;   /* total Gillespie steps taken since alloc */
} InterfaceState;

/* Build an InterfaceState from an existing spin array.
 * spins[0..l-1] must already be populated (+1 or -1).
 * seed: RNG seed for this simulation instance. */
InterfaceState *istate_alloc(const int *spins, int l, int N,
                             double ann_rate, double create_rate,
                             unsigned long seed);

/* Free all memory owned by s (does not free the original spin array). */
void istate_free(InterfaceState *s);

/* Reconstruct spin array from interface positions and left_colors.
 * spins_out must be pre-allocated with at least l ints. O(l). */
void istate_to_spins(const InterfaceState *s, int *spins_out);

/* Advance the simulation until macro_time >= target_macro_time.
 * Modifies the interface list in place. */
void istate_advance_to(InterfaceState *s, double target_macro_time);

/* Convenience accessor. */
static inline int istate_n_interfaces(const InterfaceState *s) { return s->k; }

#endif /* INTERFACE_STATE_H */
