/*
 * corr2pt_main.c — finite-N correlators of the interface (Glauber) dynamics,
 * for comparison with the exact theory of docs/section1.tex.
 *
 * Run modes (argv[3], default "equilibrium"):
 *   equilibrium : start from the reversible measure, interface on each bond
 *                 independently with p = nu/(N+nu) (even parity enforced).
 *                 The natural observables are the time-delayed correlators
 *                 G_n(T) = <sigma_0(0) sigma_n(T)>  and
 *                 H_n(T) = <eta_0(0)  eta_n(T)>,  eta_b = (1 - s_b s_{b+1})/2.
 *   disordered  : start from i.i.d. uniform spins (interfaces Bernoulli(1/2)).
 *                 The natural observable is the single-time quench correlator
 *                 C_n(T) = <sigma_0(T) sigma_n(T)>  (the equal-time column).
 *
 * Both modes also record the equal-time correlator at every snapshot and the
 * two-time pair correlators, for stationarity / ageing diagnostics.
 *
 * Usage:  ./Corr2pt [config.txt] [n_sims] [equilibrium|disordered]
 *
 * Output (one tidy text file per observable, columns n  n/l  mean  stderr ... ):
 *   corr2pt_multi.txt        G_n(T)
 *   iface_corr2pt_multi.txt  H_n(T)
 *   corr_eqtime_multi.txt    C_n at snapshot times {0, T_1, ..., T_NT}
 *   corr2pt_pairs.txt        two-time pair correlators
 *   corr2pt_cov.bin          raw second moments of G for a covariance chi^2
 *   configuration.json       N, nu, L, ann, create, derived block, corr_start
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
#include "spin_init.h"

#include <gsl/gsl_rng.h>

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

/* Interface (bond) occupation eta_b = (1 - s_b s_{b+1})/2 in {0,1}. */
static void spins_to_eta(const int *spins, int l, int *eta)
{
    for (int b = 0; b < l; b++) {
        int next = (b + 1 == l) ? 0 : b + 1;
        eta[b] = (1 - spins[b] * spins[next]) / 2;
    }
}

int main(int argc, char **argv)
{
    const char *config_path = (argc > 1) ? argv[1] : "config.txt";
    int n_sims = (argc > 2) ? atoi(argv[2]) : N_SIMS;
    if (n_sims <= 0) n_sims = N_SIMS;
    const char *start_mode = (argc > 3) ? argv[3] : "equilibrium";
    int disordered = (strcmp(start_mode, "disordered") == 0);
    if (!disordered && strcmp(start_mode, "equilibrium") != 0) {
        fprintf(stderr, "Unknown start mode '%s' (use equilibrium|disordered)\n",
                start_mode);
        return 1;
    }

    ModelConfig *model = calloc(1, sizeof(ModelConfig));
    read_model_config_file(model, config_path);

    int    l = model->N * model->L;
    double p = model->p_equilibrium;   /* = nu/(N+nu); the equilibrium density */

    printf("=== Finite-N correlators (%s start) ===\n", start_mode);
    printf("  N=%d  nu=%.6g  L=%d  l=%d  p=%.6f  config=%s\n",
           model->N, model->nu, model->L, l, p, config_path);
    printf("  N_sims=%d  ann=%.10g  create=%.10g\n",
           n_sims, model->annihilation, model->creation);
    print_derived_banner(model);
    printf("  T values:");
    for (int t = 0; t < N_T; t++) printf(" %.2f", T_VALUES[t]);
    printf("\n\n");

    /* Output directory */
    char outdir[512];
    snprintf(outdir, sizeof(outdir),
             "results_N%d_nu%.2f_L%d_%s_Ann%.2f_create%.2e_corr2pt",
             model->N, model->nu, model->L, start_mode,
             model->annihilation, model->creation);
    if (mkdir(outdir, 0755) && errno != EEXIST) { perror("mkdir"); exit(1); }

    int n_out = l / 2 + 1;

    /* Accumulators for G_n(T) = <sigma_0(0) sigma_n(T)>: sum and sum-of-squares. */
    double **G_sum = calloc(N_T, sizeof(double *));
    double **G_sq  = calloc(N_T, sizeof(double *));
    /* Accumulators for H_n(T) = <eta_0(0) eta_n(T)>. */
    double **H_sum = calloc(N_T, sizeof(double *));
    double **H_sq  = calloc(N_T, sizeof(double *));
    for (int t = 0; t < N_T; t++) {
        G_sum[t] = calloc(n_out, sizeof(double));
        G_sq[t]  = calloc(n_out, sizeof(double));
        H_sum[t] = calloc(n_out, sizeof(double));
        H_sq[t]  = calloc(n_out, sizeof(double));
    }

    /* Full second-moment accumulator sum_sims xc_n * xc_m per T, for a
     * covariance-aware chi^2 test of MC vs exact formula in Python. */
    double **G_xx = calloc(N_T, sizeof(double *));
    for (int t = 0; t < N_T; t++)
        G_xx[t] = calloc((size_t)n_out * n_out, sizeof(double));
    double *xc_vec = malloc(n_out * sizeof(double));

    /* Equal-time correlator C(r) at the N_T+1 snapshot times {0, T_1, ..., T_NT}.
     * In the equilibrium start these must all coincide with the static
     * correlator eta^|r| (stationarity); after a disordered quench the columns
     * are the single-time C_n(T). */
    double **C_sum = calloc(N_T + 1, sizeof(double *));
    double **C_sq  = calloc(N_T + 1, sizeof(double *));
    for (int t = 0; t <= N_T; t++) {
        C_sum[t] = calloc(n_out, sizeof(double));
        C_sq[t]  = calloc(n_out, sizeof(double));
    }

    /* Two-time correlators between all pairs of positive snapshot times. */
    int n_pairs = N_T * (N_T - 1) / 2;
    double **P_sum = calloc(n_pairs, sizeof(double *));
    double **P_sq  = calloc(n_pairs, sizeof(double *));
    for (int q = 0; q < n_pairs; q++) {
        P_sum[q] = calloc(n_out, sizeof(double));
        P_sq[q]  = calloc(n_out, sizeof(double));
    }

    int *spins_init = malloc(l * sizeof(int));
    int *eta_init   = malloc(l * sizeof(int));
    int *eta_snap   = malloc(l * sizeof(int));
    int **snaps     = malloc(N_T * sizeof(int *));
    for (int t = 0; t < N_T; t++)
        snaps[t] = malloc(l * sizeof(int));

    unsigned long base_seed = read_base_seed();
    /* Dedicated GSL stream for the initial condition. */
    gsl_rng *ic_rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(ic_rng, base_seed ^ 0x9E3779B97F4A7C15UL);

    for (int sim = 0; sim < n_sims; sim++) {
        if (sim % 200 == 0)
            printf("  sim %d/%d\n", sim, n_sims), fflush(stdout);

        unsigned long seed = base_seed ^ ((unsigned long)(sim + 1) * 6364136223846793005UL);

        /* Sample the initial condition. */
        if (disordered)
            sample_disordered_spins(spins_init, l, ic_rng);
        else
            sample_equilibrium_interfaces(spins_init, l, p, ic_rng);

        /* Build interface state from the IC. */
        InterfaceState *gs = istate_alloc(spins_init, l, model->N,
                                          model->annihilation, model->creation, seed);

        /* Record spin array at t=0 and its interface field. */
        istate_to_spins(gs, spins_init);
        spins_to_eta(spins_init, l, eta_init);

        /* Advance sequentially to each T and record. */
        for (int ti = 0; ti < N_T; ti++) {
            istate_advance_to(gs, T_VALUES[ti]);
            istate_to_spins(gs, snaps[ti]);
        }
        istate_free(gs);

        /* Translation-averaged time-delayed correlators.
         * G_n(T): xcorr[n] = (1/l) sum_x spins_init[x] * snaps[ti][(x+n)%l].
         * H_n(T): same with the interface fields eta_b. */
        for (int ti = 0; ti < N_T; ti++) {
            spins_to_eta(snaps[ti], l, eta_snap);
            for (int n = 0; n < n_out; n++) {
                double xc = 0.0, hc = 0.0;
                for (int x = 0; x < l; x++) {
                    int y = x + n; if (y >= l) y -= l;
                    xc += (double)spins_init[x] * (double)snaps[ti][y];
                    hc += (double)eta_init[x]   * (double)eta_snap[y];
                }
                xc /= l; hc /= l;
                xc_vec[n] = xc;
                G_sum[ti][n] += xc;  G_sq[ti][n] += xc * xc;
                H_sum[ti][n] += hc;  H_sq[ti][n] += hc * hc;
            }
            /* outer product for the covariance matrix of G */
            for (int n = 0; n < n_out; n++) {
                double *row = G_xx[ti] + (size_t)n * n_out;
                double xn = xc_vec[n];
                for (int m = 0; m < n_out; m++)
                    row[m] += xn * xc_vec[m];
            }
        }

        /* Translation-averaged equal-time correlator C(r) at every snapshot. */
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

    char out_path[600];

    /* G_n(T) */
    snprintf(out_path, sizeof(out_path), "%s/corr2pt_multi.txt", outdir);
    {
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
    }

    /* H_n(T) — time-delayed interface correlator */
    snprintf(out_path, sizeof(out_path), "%s/iface_corr2pt_multi.txt", outdir);
    {
        FILE *fout = fopen(out_path, "w");
        fprintf(fout, "# n\tn/l");
        for (int ti = 0; ti < N_T; ti++)
            fprintf(fout, "\tH_T%.2f\tstderr_T%.2f", T_VALUES[ti], T_VALUES[ti]);
        fprintf(fout, "\n");
        for (int n = 0; n < n_out; n++) {
            fprintf(fout, "%d\t%.6f", n, (double)n / l);
            for (int ti = 0; ti < N_T; ti++) {
                double mean  = H_sum[ti][n] / n_sims;
                double var   = H_sq[ti][n]  / n_sims - mean * mean;
                double sterr = sqrt(fmax(0.0, var) / n_sims);
                fprintf(fout, "\t%.8f\t%.8f", mean, sterr);
            }
            fprintf(fout, "\n");
        }
        fclose(fout);
        printf("Wrote %s\n", out_path);
    }

    /* Covariance raw moments (binary) */
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
    {
        FILE *fout = fopen(out_path, "w");
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
    }

    /* Two-time pair correlators */
    snprintf(out_path, sizeof(out_path), "%s/corr2pt_pairs.txt", outdir);
    {
        FILE *fout = fopen(out_path, "w");
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
    }

    /* Config JSON (full derived block + run-specific fields). */
    {
        char path[600]; snprintf(path, sizeof(path), "%s/configuration.json", outdir);
        FILE *f = fopen(path, "w");
        fprintf(f, "{\n");
        fprintf(f, "  \"N\": %d,\n",            model->N);
        fprintf(f, "  \"L\": %d,\n",            model->L);
        fprintf(f, "  \"T\": %lf,\n",           model->T);
        fprintf(f, "  \"N_sims\": %d,\n",       n_sims);
        fprintf(f, "  \"corr_start\": \"%s\",\n", start_mode);
        fprintf(f, "  \"p_init\": %.10g,\n",    disordered ? 0.5 : p);
        fprintf(f, "  \"T_values\": [");
        for (int t = 0; t < N_T; t++)
            fprintf(f, "%s%.4g", t ? ", " : "", T_VALUES[t]);
        fprintf(f, "],\n");
        fprint_derived_block(f, model);
        fprintf(f, "}\n");
        fclose(f);
    }

    /* Theory-vs-MC plotter (non-fatal if matplotlib/scipy are missing). */
    char cmd[800];
    snprintf(cmd, sizeof(cmd),
             "python3 python_scripts/plot_theory_vs_mc.py \"%s\"", outdir);
    if (system(cmd) != 0)
        fprintf(stderr, "Warning: plot_theory_vs_mc.py failed.\n");

    gsl_rng_free(ic_rng);
    for (int t = 0; t < N_T; t++) {
        free(G_sum[t]); free(G_sq[t]); free(H_sum[t]); free(H_sq[t]);
        free(snaps[t]); free(G_xx[t]);
    }
    free(G_sum); free(G_sq); free(H_sum); free(H_sq);
    free(G_xx); free(xc_vec);
    for (int t = 0; t <= N_T; t++) { free(C_sum[t]); free(C_sq[t]); }
    free(C_sum); free(C_sq);
    for (int q = 0; q < n_pairs; q++) { free(P_sum[q]); free(P_sq[q]); }
    free(P_sum); free(P_sq);
    free(snaps); free(spins_init); free(eta_init); free(eta_snap);
    if (model->initialization_param) free(model->initialization_param);
    free(model);
    return 0;
}
