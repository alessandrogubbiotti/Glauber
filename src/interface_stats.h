#ifndef INTERFACE_STATS_H
#define INTERFACE_STATS_H

#include "interface_state.h"

/*
 * Running accumulators for time- and ensemble-averaged statistics.
 *
 * At each snapshot the caller passes the current InterfaceState and the
 * reconstructed spin array; both correlation functions are accumulated.
 *
 * At the end, stats_write() normalises by n_samples and writes two files:
 *   spin_corr.txt   — equal-time spin-spin correlation C(r)
 *   iface_corr.txt  — interface pair correlation g(r)
 */

typedef struct {
    double *spin_corr;   /* spin_corr[r]  += (1/l) sum_x sigma(x)*sigma(x+r) */
    double *iface_corr;  /* iface_corr[r] += pairwise interface distance histogram */
    double  k_sum;       /* sum of interface counts across snapshots */
    long    n_samples;
    int     l;
} StatsAccum;

void stats_alloc      (StatsAccum *a, int l);
void stats_free       (StatsAccum *a);
void stats_accumulate (StatsAccum *a, const InterfaceState *s, const int *spins);
void stats_write      (const StatsAccum *a, const char *outdir);

#endif /* INTERFACE_STATS_H */
