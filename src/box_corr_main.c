/*
 * box_corr_main.c — box-magnetization correlation, exact-theory-comparable.
 *
 *   C(x, T) = < m_box(L/2, 0) * m_box(x, T) >,
 *   m_box(x, t) = average spin over the 2w+1 sites centred at x, w = round(h*N/2).
 *
 * Sampling: one INDEPENDENT ensemble of nsims realizations per (x_j, T_k)
 * grid point — fresh equilibrium IC, fresh trajectory, independent RNG
 * stream — so the grid-point estimates are statistically independent and the
 * per-time chi^2/dof against the exact theory is meaningful.
 *
 * IC: the exact reversible measure via sample_equilibrium_interfaces()
 * (interface count ~ even-conditioned Binomial(l, p), p = nu/(N+nu), placed
 * uniformly with exclusion).  Rates: the closing "glauber" choice in code
 * units (hop rate 1 per direction): ann = 1+gamma, create = 1-gamma with
 * gamma = (N^2-nu^2)/(N^2+nu^2).  Macroscopic time tau = micro/N^2; the
 * analytic formulas are evaluated at t = 2*N^2*tau by the Python plotter:
 *
 *   python3 python_scripts/mc_box_correlation.py --plot-only <outdir>
 *
 * Usage:
 *   ./BoxCorr [--N 12] [--L 4] [--nu 1.0] [--h 0.5] [--M 15]
 *             [--T 0.05,0.1,0.2,0.4] [--nsims 1000] [--seed 0]
 *             [--out box_corr_c]
 *
 * OpenMP: compiled with -fopenmp, grid points run in parallel (each point's
 * ensemble is a single task; RNG streams are derived per (point, sim) by
 * splitmix64, so results are independent of the thread schedule).
 */

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "interface_state.h"
#include "spin_init.h"

#include <gsl/gsl_rng.h>

#define MAX_T 32

/* splitmix64: derive independent 64-bit seeds from (base, point, sim). */
static unsigned long long mix64(unsigned long long z)
{
    z += 0x9E3779B97F4A7C15ULL;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

/* =========================================================================
 * Checkpoint / resume.
 *
 * Each grid point pt writes a disjoint slot Cmc[pt]/Cerr[pt] computed purely
 * from (seed, pt), so checkpointing completed points and skipping them on
 * resume is byte-identical to an uninterrupted run (the seed must be pinned).
 * File <out>/checkpoint.bin: a magic+version+params header, then T_values[K],
 * then done[K*M], Cmc[K*M], Cerr[K*M].  Written field-group-wise (no struct
 * blob) to avoid padding; written atomically via a .tmp file + rename().
 * ========================================================================= */

#define CKPT_MAGIC   0x42434B31  /* "BCK1" */
#define CKPT_VERSION 1

typedef struct {
    int                N, L, M, nsims, K, iface;
    double             nu, h;
    unsigned long long seed;     /* the PINNED seed actually used */
} CkptParams;

static void ckpt_write_atomic(const char *out, const CkptParams *p,
                              const double *T_values,
                              const int *done, const double *Cmc, const double *Cerr)
{
    int  KM = p->K * p->M;
    char tmp[600], dst[600];
    snprintf(tmp, sizeof(tmp), "%s/checkpoint.bin.tmp", out);
    snprintf(dst, sizeof(dst), "%s/checkpoint.bin", out);
    FILE *f = fopen(tmp, "wb");
    if (!f) { perror("checkpoint.bin.tmp"); return; }   /* non-fatal */
    int magic = CKPT_MAGIC, version = CKPT_VERSION;
    fwrite(&magic,    sizeof(int), 1, f);
    fwrite(&version,  sizeof(int), 1, f);
    fwrite(&p->N,     sizeof(int), 1, f);
    fwrite(&p->L,     sizeof(int), 1, f);
    fwrite(&p->M,     sizeof(int), 1, f);
    fwrite(&p->nsims, sizeof(int), 1, f);
    fwrite(&p->K,     sizeof(int), 1, f);
    fwrite(&p->iface, sizeof(int), 1, f);
    fwrite(&p->nu,    sizeof(double), 1, f);
    fwrite(&p->h,     sizeof(double), 1, f);
    fwrite(&p->seed,  sizeof(unsigned long long), 1, f);
    fwrite(T_values,  sizeof(double), p->K, f);
    fwrite(done,      sizeof(int),    KM, f);
    fwrite(Cmc,       sizeof(double), KM, f);
    fwrite(Cerr,      sizeof(double), KM, f);
    fflush(f);
    fsync(fileno(f));
    fclose(f);
    if (rename(tmp, dst) != 0) {           /* fallback if rename won't replace */
        unlink(dst);
        if (rename(tmp, dst) != 0) perror("rename checkpoint.bin");
    }
}

/* Peek only the stored seed from a checkpoint header.  Returns 1 if a
 * valid-magic checkpoint exists (seed_out set), else 0. */
static int ckpt_peek_seed(const char *out, unsigned long long *seed_out)
{
    char path[600];
    snprintf(path, sizeof(path), "%s/checkpoint.bin", out);
    FILE *f = fopen(path, "rb");
    if (!f) return 0;
    int magic = 0, version = 0, hi[6];
    double hd[2];
    unsigned long long s = 0;
    int ok = fread(&magic,   sizeof(int), 1, f) == 1 && magic   == CKPT_MAGIC
          && fread(&version, sizeof(int), 1, f) == 1 && version == CKPT_VERSION
          && fread(hi, sizeof(int), 6, f) == 6        /* N,L,M,nsims,K,iface */
          && fread(hd, sizeof(double), 2, f) == 2     /* nu,h                */
          && fread(&s, sizeof(unsigned long long), 1, f) == 1;
    fclose(f);
    if (!ok) return 0;
    *seed_out = s;
    return 1;
}

/* Load a checkpoint from <out>/checkpoint.bin.
 *   1  = loaded, all params match (done/Cmc/Cerr filled)
 *   0  = no usable file (missing / corrupt / short read) -> start fresh
 *  -1  = valid file but parameters MISMATCH -> caller must refuse to resume */
static int ckpt_load(const char *out, const CkptParams *p,
                     const double *T_values,
                     int *done, double *Cmc, double *Cerr)
{
    int  KM = p->K * p->M;
    char path[600];
    snprintf(path, sizeof(path), "%s/checkpoint.bin", out);
    FILE *f = fopen(path, "rb");
    if (!f) return 0;

    int magic = 0, version = 0;
    CkptParams q;
    int ok = fread(&magic,   sizeof(int), 1, f) == 1 && magic   == CKPT_MAGIC
          && fread(&version, sizeof(int), 1, f) == 1 && version == CKPT_VERSION
          && fread(&q.N,     sizeof(int), 1, f) == 1
          && fread(&q.L,     sizeof(int), 1, f) == 1
          && fread(&q.M,     sizeof(int), 1, f) == 1
          && fread(&q.nsims, sizeof(int), 1, f) == 1
          && fread(&q.K,     sizeof(int), 1, f) == 1
          && fread(&q.iface, sizeof(int), 1, f) == 1
          && fread(&q.nu,    sizeof(double), 1, f) == 1
          && fread(&q.h,     sizeof(double), 1, f) == 1
          && fread(&q.seed,  sizeof(unsigned long long), 1, f) == 1;
    if (!ok) { fclose(f); return 0; }
    if (q.K < 1 || q.K > MAX_T || q.M < 1) { fclose(f); return 0; } /* implausible */

    double Tv[MAX_T];
    if ((int)fread(Tv, sizeof(double), q.K, f) != q.K) { fclose(f); return 0; }

    int match = (q.N == p->N && q.L == p->L && q.M == p->M && q.nsims == p->nsims &&
                 q.K == p->K && q.iface == p->iface && q.nu == p->nu && q.h == p->h &&
                 q.seed == p->seed);
    if (match) for (int i = 0; i < q.K; i++) if (Tv[i] != T_values[i]) match = 0;
    if (!match) { fclose(f); return -1; }

    if ((int)fread(done, sizeof(int),    KM, f) != KM ||
        (int)fread(Cmc,  sizeof(double), KM, f) != KM ||
        (int)fread(Cerr, sizeof(double), KM, f) != KM) {
        /* short / corrupt payload: discard and start fresh */
        memset(done, 0, (size_t)KM * sizeof(int));
        memset(Cmc,  0, (size_t)KM * sizeof(double));
        memset(Cerr, 0, (size_t)KM * sizeof(double));
        fclose(f);
        return 0;
    }
    fclose(f);
    return 1;
}

static double box_mean(const int *spins, int l, int center, int w)
{
    long sum = 0;
    for (int o = -w; o <= w; o++)
        sum += spins[(center + o + l) % l];
    return (double)sum / (double)(2 * w + 1);
}

/* Box-averaged interface (wall) density: eta_b = (1 - s_b s_{b+1})/2 in {0,1},
 * averaged over the 2w+1 bonds centred at `center`. */
static double box_eta_mean(const int *spins, int l, int center, int w)
{
    long sum = 0;
    for (int o = -w; o <= w; o++) {
        int b = (center + o + l) % l;
        sum += (1 - spins[b] * spins[(b + 1) % l]) / 2;
    }
    return (double)sum / (double)(2 * w + 1);
}

int main(int argc, char **argv)
{
    /* ---- defaults ---- */
    int    N = 12, L = 4, M = 15, nsims = 1000;
    double nu = 1.0, h = 0.5;
    unsigned long long seed = 0;
    const char *out = "box_corr_c";
    double T_values[MAX_T] = { 0.05, 0.1, 0.2, 0.4 };
    int    K = 4;
    int    iface = 0;   /* 0 = spin box correlator, 1 = interface (eta) memory f */
    double ckpt_interval = 5.0; /* seconds between checkpoint writes            */
    int    allow_resume  = 1;   /* auto-resume from <out>/checkpoint.bin        */

    /* ---- flag parsing ---- */
    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--N")     && i + 1 < argc) N     = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--L")     && i + 1 < argc) L     = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--nu")    && i + 1 < argc) nu    = atof(argv[++i]);
        else if (!strcmp(argv[i], "--h")     && i + 1 < argc) h     = atof(argv[++i]);
        else if (!strcmp(argv[i], "--M")     && i + 1 < argc) M     = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--nsims") && i + 1 < argc) nsims = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seed")  && i + 1 < argc) seed  = strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--out")   && i + 1 < argc) out   = argv[++i];
        else if (!strcmp(argv[i], "--observable") && i + 1 < argc)
            iface = !strcmp(argv[++i], "iface");
        else if (!strcmp(argv[i], "--T")     && i + 1 < argc) {
            K = 0;
            char *tok = strtok(argv[++i], ",");
            while (tok && K < MAX_T) { T_values[K++] = atof(tok); tok = strtok(NULL, ","); }
        }
        else if (!strcmp(argv[i], "--checkpoint-interval") && i + 1 < argc)
            ckpt_interval = atof(argv[++i]);
        else if (!strcmp(argv[i], "--no-resume")) allow_resume = 0;
        else { fprintf(stderr, "Unknown/incomplete flag: %s\n", argv[i]); return 1; }
    }
    if (N <= 0 || L <= 0 || nu <= 0 || M <= 0 || nsims <= 0 || K <= 0) {
        fprintf(stderr, "Invalid parameters.\n"); return 1;
    }
    /* Pin the seed (required for deterministic resume).  With no explicit
     * --seed, adopt a resumable checkpoint's seed if one exists; otherwise
     * generate one and freeze it. */
    int seed_explicit = (seed != 0);
    if (allow_resume && !seed_explicit) {
        unsigned long long stored;
        if (ckpt_peek_seed(out, &stored)) seed = stored;
    }
    if (!seed) seed = (unsigned long long)time(NULL) ^ ((unsigned long long)getpid() << 20);

    /* ---- derived quantities (the (N, nu) dictionary) ---- */
    int    l     = N * L;
    double gamma = ((double)N * N - nu * nu) / ((double)N * N + nu * nu);
    double eta   = ((double)N - nu) / ((double)N + nu);
    double p     = nu / ((double)N + nu);
    double ann    = 1.0 + gamma;   /* closing rates, code units (hop 1/dir) */
    double create = 1.0 - gamma;

    int w = (int)floor(h * N / 2.0 + 0.5);
    if (w < 0) w = 0;
    int boxw = 2 * w + 1;
    int c0 = l / 2;
    if (M > l) M = l;

    /* Box centres: round(j*l/M), j = 0..M-1 (same as the Python tool). */
    int *centers = malloc(M * sizeof(int));
    for (int j = 0; j < M; j++)
        centers[j] = (int)floor((double)j * l / M + 0.5) % l;

    printf("=== Box-magnetization correlation (C engine, independent ensembles) ===\n");
    printf("  N=%d L=%d l=%d nu=%g  gamma=%.6f eta=%.6f p=%.6f\n",
           N, L, l, nu, gamma, eta, p);
    printf("  ann=%.10g create=%.10g  (closing rates, hop 1/dir)\n", ann, create);
    printf("  box h=%g -> w=%d (%d sites)  reference site %d (x=L/2)\n", h, w, boxw, c0);
    printf("  observable: %s\n",
           iface ? "interface memory f  (N^2 * connected eta_box correlator)"
                 : "spin box correlator");
    printf("  M=%d space pts  K=%d times  nsims=%d per point  seed=%llu\n",
           M, K, nsims, seed);
    printf("  => %d independent ensembles (%ld trajectories total)\n",
           M * K, (long)M * K * nsims);
#ifdef _OPENMP
    printf("  OpenMP: %d threads\n", omp_get_max_threads());
#endif

    /* ---- accumulators ---- */
    double *Cmc  = calloc((size_t)K * M, sizeof(double));
    double *Cerr = calloc((size_t)K * M, sizeof(double));

    /* ---- output directory (created early so checkpoints can be written) ---- */
    if (mkdir(out, 0755) && errno != EEXIST) { perror("mkdir"); return 1; }

    /* ---- checkpoint / resume ---- */
    int  KM = K * M;
    int *done = calloc((size_t)KM, sizeof(int));
    if (!done) { perror("calloc done"); return 1; }
    CkptParams cp = { N, L, M, nsims, K, iface, nu, h, seed };
    int done_count = 0;
    if (allow_resume) {
        int rc = ckpt_load(out, &cp, T_values, done, Cmc, Cerr);
        if (rc < 0) {
            fprintf(stderr,
                "checkpoint.bin in '%s' has parameters that differ from this run; "
                "refusing to resume.\n  Use a different --out, or pass --no-resume "
                "to start over.\n", out);
            return 1;
        }
        for (int pt = 0; pt < KM; pt++) done_count += done[pt];
        if (rc == 1)
            printf("  resuming: %d/%d grid points already done\n", done_count, KM);
        else
            /* fresh start: persist the pinned seed now, so a crash before the
             * first point completes still resumes with the same seed. */
            ckpt_write_atomic(out, &cp, T_values, done, Cmc, Cerr);
        fflush(stdout);
    } else {
        /* --no-resume: drop any stale checkpoint so a later crash cannot resume
         * onto an old seed. */
        char p[600]; snprintf(p, sizeof(p), "%s/checkpoint.bin", out); unlink(p);
    }
    time_t last_ckpt = time(NULL);

    /* ---- one independent ensemble per grid point ---- */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int pt = 0; pt < K * M; pt++) {
        if (done[pt]) continue;   /* already computed in a previous run */
        int ti = pt / M, j = pt % M;
        double T = T_values[ti];

        int *spins = malloc(l * sizeof(int));
        gsl_rng *ic_rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(ic_rng, (unsigned long)mix64(seed ^ mix64((unsigned long long)pt * 2 + 1)));

        double acc = 0.0, acc2 = 0.0;
        for (int s = 0; s < nsims; s++) {
            unsigned long eseed =
                (unsigned long)mix64(seed ^ mix64(((unsigned long long)pt << 32) | (unsigned long long)s));
            sample_equilibrium_interfaces(spins, l, p, ic_rng);
            /* reference at t=0; for the interface estimator subtract p, since
             * <eta_box> = p exactly at equilibrium, so (eta_box-p) gives the
             * connected correlator and a far smaller variance. */
            double m0 = iface ? box_eta_mean(spins, l, c0, w) - p
                              : box_mean(spins, l, c0, w);

            InterfaceState *gs = istate_alloc(spins, l, N, ann, create, eseed);
            istate_advance_to(gs, T);
            istate_to_spins(gs, spins);
            double mt = iface ? box_eta_mean(spins, l, centers[j], w) - p
                              : box_mean(spins, l, centers[j], w);
            double pr = m0 * mt;
            istate_free(gs);

            acc  += pr;
            acc2 += pr * pr;
        }
        double mean = acc / nsims;
        double var  = acc2 / nsims - mean * mean;
        /* interface: scale the connected correlator by N^2 -> box-average of f. */
        double sc   = iface ? (double)N * (double)N : 1.0;

        gsl_rng_free(ic_rng);
        free(spins);

#ifdef _OPENMP
#pragma omp critical
#endif
        {
            /* store the disjoint slot + mark done + throttled checkpoint, all in
             * one critical region so each snapshot is internally consistent. */
            Cmc [(size_t)ti * M + j] = sc * mean;
            Cerr[(size_t)ti * M + j] = sc * sqrt(fmax(var, 0.0) / nsims);
            done[pt] = 1;
            done_count++;
            time_t ckpt_now = time(NULL);
            if (difftime(ckpt_now, last_ckpt) >= ckpt_interval || done_count == KM) {
                ckpt_write_atomic(out, &cp, T_values, done, Cmc, Cerr);
                last_ckpt = ckpt_now;
            }
            if (done_count % M == 0)
                printf("  %d/%d grid points done\n", done_count, K * M), fflush(stdout);
        }
    }

    /* ---- output (directory already created before the compute loop) ---- */
    char path[600];

    snprintf(path, sizeof(path), "%s/box_corr.txt", out);
    FILE *f = fopen(path, "w");
    if (!f) { perror("box_corr.txt"); return 1; }
    fprintf(f, "# x");
    for (int ti = 0; ti < K; ti++)
        fprintf(f, "\tMC_T%g\tstderr_T%g", T_values[ti], T_values[ti]);
    fprintf(f, "\n");
    for (int j = 0; j < M; j++) {
        fprintf(f, "%.6f", (double)centers[j] / N);
        for (int ti = 0; ti < K; ti++)
            fprintf(f, "\t%.8f\t%.8f",
                    Cmc[(size_t)ti * M + j], Cerr[(size_t)ti * M + j]);
        fprintf(f, "\n");
    }
    fclose(f);
    printf("Wrote %s\n", path);

    snprintf(path, sizeof(path), "%s/configuration.json", out);
    f = fopen(path, "w");
    if (!f) { perror("configuration.json"); return 1; }
    fprintf(f, "{\n");
    fprintf(f, "  \"N\": %d,\n  \"L\": %d,\n  \"nu\": %.10g,\n  \"l\": %d,\n", N, L, nu, l);
    fprintf(f, "  \"gamma\": %.12f,\n  \"eta\": %.12f,\n  \"p_equilibrium\": %.12f,\n",
            gamma, eta, p);
    fprintf(f, "  \"ann\": %.12f,\n  \"create\": %.12f,\n", ann, create);
    fprintf(f, "  \"box_width_h\": %.10g,\n  \"box_sites\": %d,\n", h, boxw);
    fprintf(f, "  \"reference_site\": %d,\n  \"reference_x\": %.10g,\n", c0, L / 2.0);
    fprintf(f, "  \"M\": %d,\n  \"T_values\": [", M);
    for (int ti = 0; ti < K; ti++)
        fprintf(f, "%s%.10g", ti ? ", " : "", T_values[ti]);
    fprintf(f, "],\n");
    fprintf(f, "  \"nsims\": %d,\n  \"mode\": \"independent\",\n", nsims);
    fprintf(f, "  \"observable\": \"%s\",\n", iface ? "iface" : "spin");
    fprintf(f, "  \"engine\": \"C\",\n");
    fprintf(f, "  \"time_scale_formula\": %.1f,\n", 2.0 * N * N);
    fprintf(f, "  \"seed\": %llu\n", seed);
    fprintf(f, "}\n");
    fclose(f);
    printf("Wrote %s\n", path);

    printf("\nOverlay the exact theory with:\n"
           "  python3 python_scripts/mc_box_correlation.py --plot-only %s\n", out);

    free(centers); free(Cmc); free(Cerr); free(done);
    return 0;
}
