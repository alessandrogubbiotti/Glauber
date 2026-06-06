#include "interface_stats.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void stats_alloc(StatsAccum *a, int l)
{
    a->l         = l;
    a->n_samples = 0;
    a->k_sum     = 0.0;
    a->spin_corr  = calloc(l / 2 + 1, sizeof(double));
    a->iface_corr = calloc(l / 2 + 1, sizeof(double));
    if (!a->spin_corr || !a->iface_corr) {
        perror("stats_alloc"); exit(EXIT_FAILURE);
    }
}

void stats_free(StatsAccum *a)
{
    free(a->spin_corr);
    free(a->iface_corr);
    a->spin_corr  = NULL;
    a->iface_corr = NULL;
}

void stats_accumulate(StatsAccum *a, const InterfaceState *s, const int *spins)
{
    int l = a->l;

    /* --- Spin-spin correlation C(r) = (1/l) sum_x sigma(x) sigma(x+r) --- */
    for (int r = 0; r <= l / 2; r++) {
        double sum = 0.0;
        for (int x = 0; x < l; x++)
            sum += spins[x] * spins[(x + r) % l];
        a->spin_corr[r] += sum / (double)l;
    }

    /* --- Interface pair correlation g(r): histogram of pairwise distances --- */
    int k = s->k;
    a->k_sum += k;
    if (k >= 2) {
        for (int i = 0; i < k; i++) {
            for (int j = i + 1; j < k; j++) {
                int diff = abs(s->ifaces[i].pos - s->ifaces[j].pos);
                int d    = (diff <= l / 2) ? diff : l - diff; /* torus distance */
                if (d >= 1 && d <= l / 2)
                    a->iface_corr[d] += 1.0;
            }
        }
    }

    a->n_samples++;
}

void stats_write(const StatsAccum *a, const char *outdir)
{
    if (a->n_samples == 0) return;
    int    l  = a->l;
    double ns = (double)a->n_samples;

    /* spin_corr.txt: columns  r_over_l   C_spin(r) */
    {
        char path[512];
        snprintf(path, sizeof(path), "%s/spin_corr.txt", outdir);
        FILE *f = fopen(path, "w");
        if (!f) { perror("spin_corr.txt"); return; }
        fprintf(f, "# r/l\tC_spin(r)\n");
        for (int r = 0; r <= l / 2; r++)
            fprintf(f, "%.6f\t%.8f\n",
                    (double)r / (double)l,
                    a->spin_corr[r] / ns);
        fclose(f);
    }

    /* iface_corr.txt: normalise to pair density.
     * g(r) = [count(r) / (n_samples * k_avg*(k_avg-1)/2)] * (l / 2)
     * where k_avg = k_sum / n_samples.
     * For a Poisson process with intensity rho = k_avg/l, the expected
     * count is rho^2 * l/2 per unit r; g(r)=0 in that case. */
    {
        char path[512];
        snprintf(path, sizeof(path), "%s/iface_corr.txt", outdir);
        FILE *f = fopen(path, "w");
        if (!f) { perror("iface_corr.txt"); return; }

        double k_avg   = a->k_sum / ns;
        double rho     = (l > 0) ? k_avg / (double)l : 0.0;
        double rho_sq  = rho * rho;

        fprintf(f, "# r/l\tg(r)\t(Poisson baseline: %.6f)\n", rho_sq);
        for (int r = 1; r <= l / 2; r++) {
            /* Empirical pair density: count pairs at distance r per snapshot,
             * normalised per unit length. */
            double pair_density = a->iface_corr[r] / ns / (double)l;
            fprintf(f, "%.6f\t%.8f\n",
                    (double)r / (double)l,
                    pair_density);
        }
        fclose(f);
    }
}
