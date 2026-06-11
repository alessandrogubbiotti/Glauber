/*
 * corr2pt_main.c — Time-delayed spin correlation G_n(T) = <sigma_0(T) sigma_n(0)>
 *
 * Starts from the Poissonian initial condition (each bond has an interface
 * independently with probability p = 1/N) — the same IC as the main
 * Gillespie simulations.  No burn-in.
 *
 * For each simulation:
 *   1. Sample Poissonian IC.
 *   2. Record full spin array at t=0.
 *   3. Advance sequentially to each T in T_VALUES and record spin array.
 *   4. Accumulate cross-correlation (translation-averaged over ring).
 *
 * Output: corr2pt_multi.txt — n, n/l, G_n(T1), stderr, G_n(T2), stderr, ...
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <sys/stat.h>

#include "configuration.h"
#include "interface_state.h"

#include <gsl/gsl_rng.h>

/* Exact conditioned-even Bernoulli(p) bond sampler (the invariant measure),
 * driven by GSL MT19937 instead of C rand(): glibc rand() is an additive-
 * feedback generator with known structural correlations (lag-31 recurrence)
 * that could imprint spurious spatial structure on the initial conditions. */
static void sample_invariant_ic(int *spins, int l, double p, gsl_rng *r)
{
    int parity;
    do {
        parity = 0;
        spins[0] = (gsl_rng_uniform(r) < 0.5) ? 1 : -1;
        int bond;
        for (int i = 1; i < l; i++) {
            bond = gsl_rng_uniform(r) < p;
            spins[i] = bond ? -spins[i - 1] : spins[i - 1];
            parity ^= bond;
        }
        bond = gsl_rng_uniform(r) < p;     /* seam bond, sampled like the rest */
        parity ^= bond;
    } while (parity);                      /* accept iff even total parity */
}

#define N_SIMS    2000   /* independent simulations */

static double T_VALUES[] = { 0.25, 0.5, 1.0, 2.0 };
#define N_T (int)(sizeof(T_VALUES)/sizeof(T_VALUES[0]))

static unsigned long read_base_seed(void)
{
    unsigned long seed = 0;
    FILE *f = fopen("/dev/urandom", "rb");
    if (f) { fread(&seed, sizeof(seed), 1, f); fclose(f); }
    if (!seed)
        seed = (unsigned long)time(NULL) ^ ((unsigned long)getpid() << 16);
    return seed;
}

int main(int argc, char **argv)
{
    const char *config_path = (argc > 1) ? argv[1] : "config.txt";
    int n_sims = (argc > 2) ? atoi(argv[2]) : N_SIMS;
    if (n_sims <= 0) n_sims = N_SIMS;

    ModelConfig *model = calloc(1, sizeof(ModelConfig));
    read_model_config_file(model, config_path);

    int    l = model->N * model->L;
    double p = model->initialization_param ?
                   *(double *)model->initialization_param : 1.0 / model->N;

    printf("=== Time-delayed spin correlation G_n(T) (Poissonian IC) ===\n");
    printf("  N=%d  L=%d  l=%d  p=%.5f  config=%s\n",
           model->N, model->L, l, p, config_path);
    printf("  N_sims=%d  ann=%.10g  create=%.10g\n",
           n_sims, model->annihilation, model->creation);
    printf("  T values:");
    for (int t = 0; t < N_T; t++) printf(" %.2f", T_VALUES[t]);
    printf("\n\n");

    /* Output directory */
    char outdir[512];
    snprintf(outdir, sizeof(outdir),
             "results_N%d_T%.1f_L%d_Ann%.2f_create%.2e_corr2pt",
             model->N, model->T, model->L,
             model->annihilation, model->creation);
    if (mkdir(outdir, 0755) && errno != EEXIST) { perror("mkdir"); exit(1); }

    /* Accumulators for G_n(T): sum and sum-of-squares over sims */
    int n_out = l / 2 + 1;
    double **G_sum = calloc(N_T, sizeof(double *));
    double **G_sq  = calloc(N_T, sizeof(double *));
    for (int t = 0; t < N_T; t++) {
        G_sum[t] = calloc(n_out, sizeof(double));
        G_sq[t]  = calloc(n_out, sizeof(double));
    }

    /* Full second-moment accumulator sum_sims xc_n * xc_m per T, for a
     * covariance-aware chi^2 test of MC vs exact formula in Python. */
    double **G_xx = calloc(N_T, sizeof(double *));
    for (int t = 0; t < N_T; t++)
        G_xx[t] = calloc((size_t)n_out * n_out, sizeof(double));
    double *xc_vec = malloc(n_out * sizeof(double));

    /* Accumulators for the equal-time correlator C(r) at the N_T+1
     * snapshot times {0, T_1, ..., T_NT}.  Under the invariant IC these
     * must all coincide with the static correlator (stationarity test). */
    double **C_sum = calloc(N_T + 1, sizeof(double *));
    double **C_sq  = calloc(N_T + 1, sizeof(double *));
    for (int t = 0; t <= N_T; t++) {
        C_sum[t] = calloc(n_out, sizeof(double));
        C_sq[t]  = calloc(n_out, sizeof(double));
    }

    /* Two-time correlators between all pairs of positive snapshot times:
     * P[a][b] ~ <sigma_0(T_a) sigma_n(T_b)>, 0 <= a < b < N_T.
     * With the invariant IC these depend on T_b - T_a only (stationarity);
     * after a T=0 quench they exhibit ageing (Henkel, arXiv:2501.05912). */
    int n_pairs = N_T * (N_T - 1) / 2;
    double **P_sum = calloc(n_pairs, sizeof(double *));
    double **P_sq  = calloc(n_pairs, sizeof(double *));
    for (int q = 0; q < n_pairs; q++) {
        P_sum[q] = calloc(n_out, sizeof(double));
        P_sq[q]  = calloc(n_out, sizeof(double));
    }

    int *spins_init = malloc(l * sizeof(int));
    int **snaps     = malloc(N_T * sizeof(int *));
    for (int t = 0; t < N_T; t++)
        snaps[t] = malloc(l * sizeof(int));

    unsigned long base_seed = read_base_seed();
    /* Seed C rand() used by the generic spin initialisers: without this every
     * run replays the identical IC sequence (default seed 1), so the finite-
     * sample fluctuation of the initial ensemble is frozen across runs. */
    srand((unsigned int)(base_seed & 0xFFFFFFFFu));
    /* Dedicated GSL stream for the invariant-measure initialiser. */
    gsl_rng *ic_rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(ic_rng, base_seed ^ 0x9E3779B97F4A7C15UL);

    for (int sim = 0; sim < n_sims; sim++) {
        if (sim % 200 == 0)
            printf("  sim %d/%d\n", sim, n_sims), fflush(stdout);

        unsigned long seed = base_seed ^ ((unsigned long)(sim + 1) * 6364136223846793005UL);

        /* Sample the IC.  For the invariant-measure initialiser (function 6)
         * use the GSL-driven sampler; otherwise the configured init_function. */
        if (model->initialize == initialize_spins_6)
            sample_invariant_ic(spins_init, l, p, ic_rng);
        else
            model->initialize(spins_init, l, model->initialization_param);

        /* Build interface state */
        InterfaceState *gs = istate_alloc(spins_init, l, model->N,
                                          model->annihilation, model->creation, seed);

        /* Record spin array at t=0 */
        istate_to_spins(gs, spins_init);

        /* Advance sequentially to each T and record */
        for (int ti = 0; ti < N_T; ti++) {
            istate_advance_to(gs, T_VALUES[ti]);
            istate_to_spins(gs, snaps[ti]);
        }

        istate_free(gs);

        /* Translation-averaged cross-correlation.
         * xcorr[n] = (1/l) * sum_x spins_init[x] * snaps[ti][(x+n) % l]
         * Estimator for G_n(T) = <sigma_0(0) sigma_n(T)>. */
        for (int ti = 0; ti < N_T; ti++) {
            for (int n = 0; n < n_out; n++) {
                double xc = 0.0;
                for (int x = 0; x < l; x++) {
                    int y = x + n; if (y >= l) y -= l;
                    xc += (double)spins_init[x] * (double)snaps[ti][y];
                }
                xc /= l;
                xc_vec[n] = xc;
                G_sum[ti][n] += xc;
                G_sq[ti][n]  += xc * xc;
            }
            /* outer product for the covariance matrix */
            for (int n = 0; n < n_out; n++) {
                double *row = G_xx[ti] + (size_t)n * n_out;
                double xn = xc_vec[n];
                for (int m = 0; m < n_out; m++)
                    row[m] += xn * xc_vec[m];
            }
        }

        /* Translation-averaged equal-time correlator C(r) at every
         * snapshot time (index 0 = t=0, index ti+1 = T_VALUES[ti]). */
        for (int ti = 0; ti <= N_T; ti++) {
            const int *snap = (ti == 0) ? spins_init : snaps[ti - 1];
            for (int r = 0; r < n_out; r++) {
                double c = 0.0;
                for (int x = 0; x < l; x++) {
                    int y = x + r; if (y >= l) y -= l;
                    c += (double)snap[x] * (double)snap[y];
                }
                c /= l;
                C_sum[ti][r] += c;
                C_sq[ti][r]  += c * c;
            }
        }

        /* Two-time cross-correlations between snapshot pairs (T_a, T_b). */
        {
            int q = 0;
            for (int a = 0; a < N_T; a++) {
                for (int b = a + 1; b < N_T; b++, q++) {
                    for (int n = 0; n < n_out; n++) {
                        double xc = 0.0;
                        for (int x = 0; x < l; x++) {
                            int y = x + n; if (y >= l) y -= l;
                            xc += (double)snaps[a][x] * (double)snaps[b][y];
                        }
                        xc /= l;
                        P_sum[q][n] += xc;
                        P_sq[q][n]  += xc * xc;
                    }
                }
            }
        }
    }

    /* Write output */
    char out_path[600];
    snprintf(out_path, sizeof(out_path), "%s/corr2pt_multi.txt", outdir);
    FILE *fout = fopen(out_path, "w");
    fprintf(fout, "# n\tn/l");
    for (int ti = 0; ti < N_T; ti++)
        fprintf(fout, "\tG_T%.2f\tstderr_T%.2f", T_VALUES[ti], T_VALUES[ti]);
    fprintf(fout, "\n");
    for (int n = 0; n < n_out; n++) {
        fprintf(fout, "%d\t%.6f", n, (double)n / l);
        for (int ti = 0; ti < N_T; ti++) {
            double mean  = G_sum[ti][n] / n_sims;
            double var   = G_sq[ti][n]  / n_sims - mean * mean;
            double sterr = sqrt(fmax(0.0, var) / n_sims);
            fprintf(fout, "\t%.8f\t%.8f", mean, sterr);
        }
        fprintf(fout, "\n");
    }
    fclose(fout);
    printf("\nWrote %s\n", out_path);

    /* Covariance raw moments (binary): header [N_T, n_out, n_sims] as int32,
     * then per T the n_out*n_out sums of xc_n*xc_m as float64. */
    snprintf(out_path, sizeof(out_path), "%s/corr2pt_cov.bin", outdir);
    {
        FILE *fb = fopen(out_path, "wb");
        int hdr[3] = { N_T, n_out, n_sims };
        fwrite(hdr, sizeof(int), 3, fb);
        for (int ti = 0; ti < N_T; ti++)
            fwrite(G_xx[ti], sizeof(double), (size_t)n_out * n_out, fb);
        fclose(fb);
        printf("Wrote %s\n", out_path);
    }

    /* Equal-time correlator at all snapshot times */
    snprintf(out_path, sizeof(out_path), "%s/corr_eqtime_multi.txt", outdir);
    fout = fopen(out_path, "w");
    fprintf(fout, "# r\tr/l\tC_t0.00\tstderr_t0.00");
    for (int ti = 0; ti < N_T; ti++)
        fprintf(fout, "\tC_T%.2f\tstderr_T%.2f", T_VALUES[ti], T_VALUES[ti]);
    fprintf(fout, "\n");
    for (int r = 0; r < n_out; r++) {
        fprintf(fout, "%d\t%.6f", r, (double)r / l);
        for (int ti = 0; ti <= N_T; ti++) {
            double mean  = C_sum[ti][r] / n_sims;
            double var   = C_sq[ti][r]  / n_sims - mean * mean;
            double sterr = sqrt(fmax(0.0, var) / n_sims);
            fprintf(fout, "\t%.8f\t%.8f", mean, sterr);
        }
        fprintf(fout, "\n");
    }
    fclose(fout);
    printf("Wrote %s\n", out_path);

    /* Two-time pair correlators */
    snprintf(out_path, sizeof(out_path), "%s/corr2pt_pairs.txt", outdir);
    fout = fopen(out_path, "w");
    fprintf(fout, "# n\tn/l");
    for (int a = 0; a < N_T; a++)
        for (int b = a + 1; b < N_T; b++)
            fprintf(fout, "\tG_s%.2f_t%.2f\tstderr_s%.2f_t%.2f",
                    T_VALUES[a], T_VALUES[b], T_VALUES[a], T_VALUES[b]);
    fprintf(fout, "\n");
    for (int n = 0; n < n_out; n++) {
        fprintf(fout, "%d\t%.6f", n, (double)n / l);
        for (int q = 0; q < n_pairs; q++) {
            double mean  = P_sum[q][n] / n_sims;
            double var   = P_sq[q][n]  / n_sims - mean * mean;
            double sterr = sqrt(fmax(0.0, var) / n_sims);
            fprintf(fout, "\t%.8f\t%.8f", mean, sterr);
        }
        fprintf(fout, "\n");
    }
    fclose(fout);
    printf("Wrote %s\n", out_path);

    /* Config JSON */
    {
        char path[600]; snprintf(path, sizeof(path), "%s/configuration.json", outdir);
        FILE *f = fopen(path, "w");
        fprintf(f, "{\n  \"N\":%d,\n  \"L\":%d,\n  \"p_init\":%.10f,\n"
                   "  \"N_sims\":%d,\n  \"annihilation\":%.10f,\n  \"creation\":%.10e\n}\n",
                model->N, model->L, p, n_sims, model->annihilation, model->creation);
        fclose(f);
    }

    /* Plotter */
    char cmd[800];
    snprintf(cmd, sizeof(cmd),
             "python3 python_scripts/plot_time_delayed_corr.py \"%s\"", outdir);
    system(cmd);

    for (int t = 0; t < N_T; t++) {
        free(G_sum[t]); free(G_sq[t]); free(snaps[t]); free(G_xx[t]);
    }
    free(G_xx); free(xc_vec);
    for (int t = 0; t <= N_T; t++) { free(C_sum[t]); free(C_sq[t]); }
    free(C_sum); free(C_sq);
    for (int q = 0; q < n_pairs; q++) { free(P_sum[q]); free(P_sq[q]); }
    free(P_sum); free(P_sq);
    free(G_sum); free(G_sq); free(snaps);
    free(spins_init);
    if (model->initialization_param) free(model->initialization_param);
    free(model);
    return 0;
}
