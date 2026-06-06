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

int main(void)
{
    ModelConfig *model = calloc(1, sizeof(ModelConfig));
    read_model_config_file(model, "config.txt");

    int    l = model->N * model->L;
    double p = model->initialization_param ?
                   *(double *)model->initialization_param : 1.0 / model->N;

    printf("=== Time-delayed spin correlation G_n(T) (Poissonian IC) ===\n");
    printf("  N=%d  L=%d  l=%d  p=%.5f\n", model->N, model->L, l, p);
    printf("  N_sims=%d  ann=%.5g  create=%.5g\n",
           N_SIMS, model->annihilation, model->creation);
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

    int *spins_init = malloc(l * sizeof(int));
    int **snaps     = malloc(N_T * sizeof(int *));
    for (int t = 0; t < N_T; t++)
        snaps[t] = malloc(l * sizeof(int));

    unsigned long base_seed = read_base_seed();

    for (int sim = 0; sim < N_SIMS; sim++) {
        if (sim % 200 == 0)
            printf("  sim %d/%d\n", sim, N_SIMS), fflush(stdout);

        unsigned long seed = base_seed ^ ((unsigned long)(sim + 1) * 6364136223846793005UL);

        /* Sample Poissonian IC using the configured init_function / init_param */
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
                G_sum[ti][n] += xc;
                G_sq[ti][n]  += xc * xc;
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
            double mean  = G_sum[ti][n] / N_SIMS;
            double var   = G_sq[ti][n]  / N_SIMS - mean * mean;
            double sterr = sqrt(fmax(0.0, var) / N_SIMS);
            fprintf(fout, "\t%.8f\t%.8f", mean, sterr);
        }
        fprintf(fout, "\n");
    }
    fclose(fout);
    printf("\nWrote %s\n", out_path);

    /* Config JSON */
    {
        char path[600]; snprintf(path, sizeof(path), "%s/configuration.json", outdir);
        FILE *f = fopen(path, "w");
        fprintf(f, "{\n  \"N\":%d,\n  \"L\":%d,\n  \"p_init\":%.6f,\n"
                   "  \"N_sims\":%d,\n  \"annihilation\":%.6f,\n  \"creation\":%.5e\n}\n",
                model->N, model->L, p, N_SIMS, model->annihilation, model->creation);
        fclose(f);
    }

    /* Plotter */
    char cmd[800];
    snprintf(cmd, sizeof(cmd),
             "python3 python_scripts/plot_time_delayed_corr.py \"%s\"", outdir);
    system(cmd);

    for (int t = 0; t < N_T; t++) {
        free(G_sum[t]); free(G_sq[t]); free(snaps[t]);
    }
    free(G_sum); free(G_sq); free(snaps);
    free(spins_init);
    if (model->initialization_param) free(model->initialization_param);
    free(model);
    return 0;
}
