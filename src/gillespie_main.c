#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "configuration.h"
#include "configuration_statistics.h"
#include "interface_state.h"
#include "interface_stats.h"

/* -------------------------------------------------------------------------
 * Timing (wall-clock, monotonic)
 * ------------------------------------------------------------------------- */
static double now_seconds(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* -------------------------------------------------------------------------
 * Seeding
 * ------------------------------------------------------------------------- */
static unsigned long read_base_seed(void)
{
    /* Reproducibility hook: a fixed seed via $GILLESPIE_SEED (used by the
     * serial-vs-parallel validation) overrides the random source. */
    const char *env = getenv("GILLESPIE_SEED");
    if (env && *env) {
        unsigned long s = strtoul(env, NULL, 10);
        if (s) return s;
    }
    unsigned long seed = 0;
    FILE *f = fopen("/dev/urandom", "rb");
    if (f) { fread(&seed, sizeof(seed), 1, f); fclose(f); }
    if (!seed)
        seed = (unsigned long)time(NULL) ^ ((unsigned long)getpid() << 16);
    return seed;
}

/* -------------------------------------------------------------------------
 * Output directory
 * ------------------------------------------------------------------------- */
static void create_outdir(const ModelConfig *conf, char *buf, size_t sz)
{
    snprintf(buf, sz,
             "results_N%d_T%.1f_L%d_Ann%.2f_create%.2e_gillespie",
             conf->N, conf->T, conf->L, conf->annihilation, conf->creation);
    if (mkdir(buf, 0755) && errno != EEXIST) { perror("mkdir"); exit(EXIT_FAILURE); }
}

/* -------------------------------------------------------------------------
 * configuration.json
 * ------------------------------------------------------------------------- */
static void write_conf_json(const ModelConfig *conf, const char *dir)
{
    char p[600]; snprintf(p, sizeof(p), "%s/configuration.json", dir);
    FILE *f = fopen(p, "w");
    if (!f) { perror("configuration.json"); exit(EXIT_FAILURE); }
    fprintf(f, "{\n");
    fprintf(f, "  \"N\": %d,\n",               conf->N);
    fprintf(f, "  \"T\": %f,\n",               conf->T);
    fprintf(f, "  \"L\": %d,\n",               conf->L);
    fprintf(f, "  \"annihilation\": %.5f,\n",  conf->annihilation);
    fprintf(f, "  \"creation\": %.5e,\n",      conf->creation);
    fprintf(f, "  \"N_simulations\": %d,\n",   conf->N_simulations);
    fprintf(f, "  \"resolution\": %d,\n",      conf->resolution);
    fprintf(f, "  \"Micro_n_steps\": %d,\n",   conf->Micro_n_steps);
    fprintf(f, "  \"algorithm\": \"gillespie\",\n");
    fprint_derived_block(f, conf);
    fprintf(f, "}\n");
    fclose(f);
}

/* -------------------------------------------------------------------------
 * Statistics.txt header and row  (k, magnetization per sim, as before)
 * ------------------------------------------------------------------------- */
static void write_stats_header(FILE *f, int N_sim)
{
    fprintf(f, "time");
    for (int i = 0; i < N_sim; i++)
        fprintf(f, "\tsim%d_k\tsim%d_magnetization", i, i);
    fprintf(f, "\n");
}
static void write_stats_row(FILE *f, double t,
                            const int **snap, const InterfaceState **gs,
                            int N_sim, int l)
{
    fprintf(f, "%.6f", t);
    for (int i = 0; i < N_sim; i++) {
        double mag = 0;
        for (int x = 0; x < l; x++) mag += snap[i][x];
        fprintf(f, "\t%d\t%.6f", gs[i]->k, mag / l);
    }
    fprintf(f, "\n");
}

/* -------------------------------------------------------------------------
 * Binary spin snapshot (same format as the parity-scheme binary)
 * ------------------------------------------------------------------------- */
static void write_binary_snapshot(FILE *f, const int **snap, int N_sim, int l)
{
    static const int marker = -9999;
    for (int i = 0; i < N_sim; i++)
        fwrite(snap[i], sizeof(int), l, f);
    fwrite(&marker, sizeof(int), 1, f);
}

/* -------------------------------------------------------------------------
 * Equal-time spin-spin correlation  C(r) = (1/l) Σ_x σ(x)σ(x+r)
 * Accumulated into corr_acc[0..l/2]  (caller must pre-zero).
 * ------------------------------------------------------------------------- */
static void accumulate_C_r(const int *spins, int l, double *corr_acc)
{
    for (int r = 0; r <= l / 2; r++) {
        double c = 0;
        for (int x = 0; x < l; x++)
            c += spins[x] * spins[(x + r) % l];
        corr_acc[r] += c / (double)l;
    }
}

/* =========================================================================
 * main
 * ========================================================================= */
int main(void)
{
    /* --- Config ---------------------------------------------------------- */
    ModelConfig *model = calloc(1, sizeof(ModelConfig));
    if (!model) { perror("model calloc"); exit(EXIT_FAILURE); }
    read_model_config_file(model, "config.txt");

    if (model->Micro_n_steps > 0)
        fprintf(stderr, "INFO: Micro_n_steps ignored by Gillespie engine.\n");

    if (model->N <= 0 || model->T <= 0.0 || model->L <= 0 ||
        model->resolution <= 0 || model->N_simulations <= 0 ||
        !model->initialize) {
        fprintf(stderr, "Invalid config.\n"); exit(EXIT_FAILURE);
    }

    int    l       = model->N * model->L;
    int    N_sim   = model->N_simulations;
    int    Mac_T   = (int)(model->T * model->resolution);
    double inv_res = 1.0 / model->resolution;
    int    half_l  = l / 2;

    printf("=== Gillespie Interface Simulation ===\n");
    printf("  N=%d  T=%.2f  L=%d  l=%d\n", model->N, model->T, model->L, l);
    printf("  annihilation=%.5f  creation=%.2e\n",
           model->annihilation, model->creation);
    printf("  N_simulations=%d  resolution=%d  Mac_T=%d\n",
           N_sim, model->resolution, Mac_T);
    print_derived_banner(model);

    /* --- Observation sites ----------------------------------------------- */
    int x0 = l / 2;          /* site for mean-spin evolution            */
    int x1 = l / 4;          /* first site for two-point correlation    */
    int x2 = 3 * l / 4;      /* second site                            */
    printf("  Observation sites: x0=%d  x1=%d  x2=%d\n", x0, x1, x2);

    /* --- Selected time indices for full C(r,t) snapshots ---------------- */
    /* ~20 equally spaced snapshots across the run. */
    int n_tsel = 20;
    if (n_tsel > Mac_T) n_tsel = Mac_T;
    int *t_sel = malloc(n_tsel * sizeof(int));
    for (int i = 0; i < n_tsel; i++)
        t_sel[i] = (int)((double)i / (n_tsel - 1) * (Mac_T - 1) + 0.5);

    /* --- Output directory ----------------------------------------------- */
    char outdir[600];
    create_outdir(model, outdir, sizeof(outdir));
    write_conf_json(model, outdir);
    printf("  Output: %s\n\n", outdir);

    /* --- Time-resolved accumulators (summed over simulations) ------------ */
    /* mean_spin_sum[t]   = Σ_i σ_i(x0, t)                                  */
    /* mean_spin_sq[t]    = Σ_i σ_i(x0, t)^2  (for std dev)                 */
    /* corr_pair_sum[t]   = Σ_i σ_i(x1,t)⋅σ_i(x2,t)                        */
    /* corr_pair_sq[t]    = Σ_i [σ_i(x1,t)⋅σ_i(x2,t)]^2                    */
    /* C_r_sel[s*(half_l+1)+r] = full C(r) at selected snapshot index s      */
    double *mean_spin_sum  = calloc(Mac_T, sizeof(double));
    double *mean_spin_sq   = calloc(Mac_T, sizeof(double));
    double *corr_pair_sum  = calloc(Mac_T, sizeof(double));
    double *corr_pair_sq   = calloc(Mac_T, sizeof(double));
    double *C_r_sel        = calloc((size_t)n_tsel * (half_l + 1), sizeof(double));
    int    *t_sel_next     = malloc(sizeof(int));  /* current index into t_sel */

    /* --- Per-simulation scalars ----------------------------------------- */
    double *sim_wall_sec  = malloc(N_sim * sizeof(double));
    long   *sim_n_events  = malloc(N_sim * sizeof(long));

    /* --- Shared-output stats (interface correlations etc.) -------------- */
    StatsAccum *sa = malloc(sizeof(StatsAccum));
    stats_alloc(sa, l);

    /* --- Per-simulation working arrays ---------------------------------- */
    unsigned long base_seed = read_base_seed();
    /* Seed C rand() so that spin initialisation is not deterministic. */
    srand((unsigned int)(base_seed & 0xFFFFFFFFu));

    int   **spins_init = malloc(N_sim * sizeof(int *));
    int   **spins_snap = malloc(N_sim * sizeof(int *));
    InterfaceState **gs = malloc(N_sim * sizeof(InterfaceState *));
    for (int i = 0; i < N_sim; i++) {
        spins_init[i] = malloc(l * sizeof(int));
        spins_snap[i] = malloc(l * sizeof(int));
        model->initialize(spins_init[i], l, model->initialization_param);
        unsigned long seed = base_seed ^ ((unsigned long)(i + 1) * 2654435761UL);
        gs[i] = istate_alloc(spins_init[i], l, model->N,
                              model->annihilation, model->creation, seed);
    }

    /* --- Time-autocorrelation accumulators ------------------------------ */
    /* A(lag) = Σ_{sim,t} σ(x0,t)⋅σ(x0,t+lag),  count(lag) = N_sim*(Mac_T-lag) */
    double *autocorr_sum   = calloc(Mac_T, sizeof(double));
    long   *autocorr_count = calloc(Mac_T, sizeof(long));

    /* --- Trajectory subfolder (sim 0 only, compatible with existing plotter) */
    char traj_dir[700];
    snprintf(traj_dir, sizeof(traj_dir), "%s/trajectory", outdir);
    if (mkdir(traj_dir, 0755) && errno != EEXIST) { perror("mkdir traj"); exit(EXIT_FAILURE); }
    /* Write a configuration.json with N_simulations=1 so the plotter parses correctly. */
    {
        char tp[800]; snprintf(tp, sizeof(tp), "%s/configuration.json", traj_dir);
        FILE *tf = fopen(tp, "w");
        if (!tf) { perror("traj config.json"); exit(EXIT_FAILURE); }
        fprintf(tf, "{\n");
        fprintf(tf, "  \"N\": %d,\n",              model->N);
        fprintf(tf, "  \"T\": %f,\n",              model->T);
        fprintf(tf, "  \"L\": %d,\n",              model->L);
        fprintf(tf, "  \"annihilation\": %.5f,\n", model->annihilation);
        fprintf(tf, "  \"creation\": %.5e,\n",     model->creation);
        fprintf(tf, "  \"N_simulations\": 1,\n");
        fprintf(tf, "  \"resolution\": %d,\n",     model->resolution);
        fprintf(tf, "  \"Micro_n_steps\": %d,\n",  model->Micro_n_steps);
        fprintf(tf, "  \"algorithm\": \"gillespie\",\n");
        fprint_derived_block(tf, model);
        fprintf(tf, "}\n");
        fclose(tf);
    }
    char path[700];
    snprintf(path, sizeof(path), "%s/spins_output.bin", traj_dir);
    FILE *fbin = fopen(path, "wb");
    if (!fbin) { perror("trajectory/spins_output.bin"); exit(EXIT_FAILURE); }
    printf("  Trajectory (sim 0): %s\n", path);

    snprintf(path, sizeof(path), "%s/Statistics.txt", outdir);
    FILE *fstat = fopen(path, "w");
    if (!fstat) { perror("Statistics.txt"); exit(EXIT_FAILURE); }
    write_stats_header(fstat, N_sim);

    /* =====================================================================
     * Main simulation loop, parallelized over the N_sim independent
     * simulations.  Each gs[i] carries its own deterministic RNG stream
     * (seeded per i in the setup loop above), so the result is independent of
     * the thread schedule up to floating-point summation order in the
     * reductions.  Each thread keeps private scratch + private accumulators
     * and sums them into the shared arrays under one critical section.
     *
     * NB: sim_wall_sec[] now measures per-sim wall time UNDER CONTENTION
     * (several sims run at once), so timing_stats.txt is no longer a clean
     * per-sim cost; sim_n_events and all physics outputs are unaffected.
     * ===================================================================== */
    printf("Running %d simulations...\n", N_sim);
#ifdef _OPENMP
    printf("  OpenMP: %d threads\n", omp_get_max_threads());
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        /* Per-thread scratch. */
        double *spin_hist   = malloc(Mac_T * sizeof(double));
        double *C_r_scratch = calloc(half_l + 1, sizeof(double));

        /* Per-thread accumulators (reduced into the shared arrays below). */
        double *mean_spin_sum_t  = calloc(Mac_T, sizeof(double));
        double *mean_spin_sq_t   = calloc(Mac_T, sizeof(double));
        double *corr_pair_sum_t  = calloc(Mac_T, sizeof(double));
        double *corr_pair_sq_t   = calloc(Mac_T, sizeof(double));
        double *C_r_sel_t        = calloc((size_t)n_tsel * (half_l + 1), sizeof(double));
        double *autocorr_sum_t   = calloc(Mac_T, sizeof(double));
        long   *autocorr_count_t = calloc(Mac_T, sizeof(long));
        StatsAccum sa_t; stats_alloc(&sa_t, l);

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
        for (int i = 0; i < N_sim; i++) {
            double t_start = now_seconds();
            int t_sel_cursor = 0; /* index into t_sel for this simulation */

            for (int t = 0; t < Mac_T; t++) {
                double t_macro = (double)(t + 1) * inv_res;
                istate_advance_to(gs[i], t_macro);
                istate_to_spins(gs[i], spins_snap[i]);

                /* ---- Scalar observables at every snapshot ---- */
                double s0 = spins_snap[i][x0];
                spin_hist[t]        = s0;
                mean_spin_sum_t[t] += s0;
                mean_spin_sq_t [t] += s0 * s0;

                double sp = (double)spins_snap[i][x1] * spins_snap[i][x2];
                corr_pair_sum_t[t] += sp;
                corr_pair_sq_t [t] += sp * sp;

                /* ---- Full C(r,t) at selected snapshots ---- */
                if (t_sel_cursor < n_tsel && t == t_sel[t_sel_cursor]) {
                    memset(C_r_scratch, 0, (half_l + 1) * sizeof(double));
                    accumulate_C_r(spins_snap[i], l, C_r_scratch);
                    double *dest = C_r_sel_t + (size_t)t_sel_cursor * (half_l + 1);
                    for (int r = 0; r <= half_l; r++) dest[r] += C_r_scratch[r];
                    t_sel_cursor++;
                }

                /* ---- Save trajectory for sim 0 (only one thread runs i==0) ---- */
                if (i == 0) {
                    static const int marker = -9999;
                    fwrite(spins_snap[0], sizeof(int), l, fbin);
                    fwrite(&marker, sizeof(int), 1, fbin);
                }

                /* ---- Combined iface+spin stats (time-averaged) ---- */
                stats_accumulate(&sa_t, gs[i], spins_snap[i]);
            }

            sim_wall_sec[i] = now_seconds() - t_start;
            sim_n_events[i] = gs[i]->n_events;

            /* ---- Time autocorrelation from this simulation's history ---- */
            for (int lag = 0; lag < Mac_T; lag++) {
                int pairs = Mac_T - lag;
                for (int t = 0; t < pairs; t++)
                    autocorr_sum_t[lag] += spin_hist[t] * spin_hist[t + lag];
                autocorr_count_t[lag] += pairs;
            }

            if (i % 50 == 0) {
#ifdef _OPENMP
#pragma omp critical
#endif
                printf("  sim %4d/%d  wall=%.2fs  events=%ld\n",
                       i + 1, N_sim, sim_wall_sec[i], sim_n_events[i]);
            }
        }

        /* ---- Reduce per-thread accumulators into the shared arrays ---- */
#ifdef _OPENMP
#pragma omp critical
#endif
        {
            for (int t = 0; t < Mac_T; t++) {
                mean_spin_sum[t]   += mean_spin_sum_t[t];
                mean_spin_sq [t]   += mean_spin_sq_t [t];
                corr_pair_sum[t]   += corr_pair_sum_t[t];
                corr_pair_sq [t]   += corr_pair_sq_t [t];
                autocorr_sum [t]   += autocorr_sum_t [t];
                autocorr_count[t]  += autocorr_count_t[t];
            }
            for (size_t k = 0; k < (size_t)n_tsel * (half_l + 1); k++)
                C_r_sel[k] += C_r_sel_t[k];
            for (int r = 0; r <= l / 2; r++) {
                sa->spin_corr[r]  += sa_t.spin_corr[r];
                sa->iface_corr[r] += sa_t.iface_corr[r];
            }
            sa->k_sum     += sa_t.k_sum;
            sa->n_samples += sa_t.n_samples;
        }

        stats_free(&sa_t);
        free(spin_hist); free(C_r_scratch);
        free(mean_spin_sum_t); free(mean_spin_sq_t);
        free(corr_pair_sum_t); free(corr_pair_sq_t);
        free(C_r_sel_t); free(autocorr_sum_t); free(autocorr_count_t);
    }

    /* Write the single combined binary: re-run would be needed to write
       interleaved format — instead write one representative binary for sim 0. */
    /* For the visualisation script, write sim 0's spins snapshot separately. */
    /* Actually write the Statistics.txt row-by-row: requires rerunning or
       storing.  For simplicity, write a compact summary. */

    /* =====================================================================
     * Write time-resolved outputs
     * ===================================================================== */
    double inv_Nsim = 1.0 / N_sim;

    /* mean_spin_t.txt */
    snprintf(path, sizeof(path), "%s/mean_spin_t.txt", outdir);
    {
        FILE *f = fopen(path, "w");
        fprintf(f, "# t\tmean_sigma_x%d\tstd_sigma_x%d\n", x0, x0);
        for (int t = 0; t < Mac_T; t++) {
            double mean = mean_spin_sum[t] * inv_Nsim;
            double var  = mean_spin_sq[t]  * inv_Nsim - mean * mean;
            double std  = (var > 0) ? sqrt(var) : 0.0;
            fprintf(f, "%.6f\t%.8f\t%.8f\n",
                    (double)(t + 1) * inv_res, mean, std);
        }
        fclose(f);
        printf("Wrote %s\n", path);
    }

    /* corr_pair_t.txt */
    snprintf(path, sizeof(path), "%s/corr_pair_t.txt", outdir);
    {
        FILE *f = fopen(path, "w");
        fprintf(f, "# t\tC(x%d,x%d,t)\tstd\n", x1, x2);
        for (int t = 0; t < Mac_T; t++) {
            double mean = corr_pair_sum[t] * inv_Nsim;
            double var  = corr_pair_sq[t]  * inv_Nsim - mean * mean;
            double std  = (var > 0) ? sqrt(var) : 0.0;
            fprintf(f, "%.6f\t%.8f\t%.8f\n",
                    (double)(t + 1) * inv_res, mean, std);
        }
        fclose(f);
        printf("Wrote %s\n", path);
    }

    /* time_autocorr.txt — A(tau) = <sigma(x0,t) sigma(x0,t+tau)>_{t,sim} */
    snprintf(path, sizeof(path), "%s/time_autocorr.txt", outdir);
    {
        /* Overall mean spin (time- and ensemble-averaged). */
        double grand_mean = 0;
        long   grand_n    = 0;
        for (int t = 0; t < Mac_T; t++) {
            grand_mean += mean_spin_sum[t];
            grand_n    += N_sim;
        }
        grand_mean /= grand_n;

        FILE *f = fopen(path, "w");
        fprintf(f, "# tau\tA(tau)=<s(t)s(t+tau)>\tG(tau)=A(tau)-<s>^2\t<s>=%.6f\n",
                grand_mean);
        for (int lag = 0; lag < Mac_T; lag++) {
            double A = autocorr_sum[lag] / (double)autocorr_count[lag];
            double G = A - grand_mean * grand_mean;   /* connected */
            fprintf(f, "%.6f\t%.8f\t%.8f\n",
                    (double)lag * inv_res, A, G);
        }
        fclose(f);
        printf("Wrote %s\n", path);
    }

    /* spin_corr_sel_t.txt — full C(r) at selected snapshots */
    snprintf(path, sizeof(path), "%s/spin_corr_sel_t.txt", outdir);
    {
        FILE *f = fopen(path, "w");
        /* header */
        fprintf(f, "# r/l");
        for (int s = 0; s < n_tsel; s++)
            fprintf(f, "\tC_t=%.3f", (double)(t_sel[s] + 1) * inv_res);
        fprintf(f, "\n");
        for (int r = 0; r <= half_l; r++) {
            fprintf(f, "%.6f", (double)r / l);
            for (int s = 0; s < n_tsel; s++)
                fprintf(f, "\t%.8f",
                        C_r_sel[(size_t)s * (half_l + 1) + r] * inv_Nsim);
            fprintf(f, "\n");
        }
        fclose(f);
        printf("Wrote %s\n", path);
    }

    /* timing_stats.txt */
    snprintf(path, sizeof(path), "%s/timing_stats.txt", outdir);
    {
        double tw_sum = 0, tw_sq = 0;
        double te_sum = 0, te_sq = 0;
        for (int i = 0; i < N_sim; i++) {
            double w = sim_wall_sec[i], e = (double)sim_n_events[i];
            tw_sum += w; tw_sq += w * w;
            te_sum += e; te_sq += e * e;
        }
        double mean_w = tw_sum * inv_Nsim;
        double std_w  = sqrt(fmax(0, tw_sq * inv_Nsim - mean_w * mean_w));
        double mean_e = te_sum * inv_Nsim;
        double std_e  = sqrt(fmax(0, te_sq * inv_Nsim - mean_e * mean_e));

        FILE *f = fopen(path, "w");
        fprintf(f, "N_simulations\t%d\n",     N_sim);
        fprintf(f, "mean_wall_time_s\t%.6f\n", mean_w);
        fprintf(f, "std_wall_time_s\t%.6f\n",  std_w);
        fprintf(f, "min_wall_time_s\t%.6f\n",  sim_wall_sec[0]);
        fprintf(f, "max_wall_time_s\t%.6f\n",  sim_wall_sec[N_sim-1]);
        fprintf(f, "mean_n_events\t%.1f\n",    mean_e);
        fprintf(f, "std_n_events\t%.1f\n",     std_e);
        fprintf(f, "events_per_macrotime\t%.2f\n",
                   mean_e / (model->T * (double)(model->N * model->N)));
        fclose(f);

        printf("\n=== Timing & Events ===\n");
        printf("  Wall time per sim: %.4f ± %.4f s\n", mean_w, std_w);
        printf("  Gillespie events:  %.0f ± %.0f\n",   mean_e, std_e);
        printf("  Events per N²⋅T:   %.2f (expected ~l)\n",
               mean_e / (model->T * (double)(model->N * model->N)));
        printf("Wrote %s\n", path);
    }

    /* ---- Combined spin and interface correlation (time-averaged) ---- */
    stats_write(sa, outdir);
    printf("Wrote spin_corr.txt and iface_corr.txt\n");

    fclose(fbin);
    fclose(fstat);

    /* =====================================================================
     * Clean up
     * ===================================================================== */
    for (int i = 0; i < N_sim; i++) {
        istate_free(gs[i]);
        free(spins_init[i]);
        free(spins_snap[i]);
    }
    free(gs); free(spins_init); free(spins_snap);
    free(mean_spin_sum); free(mean_spin_sq);
    free(corr_pair_sum); free(corr_pair_sq);
    free(C_r_sel);
    free(sim_wall_sec); free(sim_n_events);
    free(autocorr_sum); free(autocorr_count);
    free(t_sel); free(t_sel_next);
    stats_free(sa); free(sa);
    if (model->initialization_param) free(model->initialization_param);
    free(model);

    /* =====================================================================
     * Python scripts
     * ===================================================================== */
    char cmd[800];
    snprintf(cmd, sizeof(cmd),
             "python3 python_scripts/plot_gillespie.py \"%s\"", outdir);
    if (system(cmd) != 0)
        fprintf(stderr, "Warning: plot_gillespie.py failed.\n");

    printf("\nDone.\n");
    return 0;
}
