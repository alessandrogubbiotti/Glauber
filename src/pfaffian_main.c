/*
 * pfaffian_main.c — the Pfaffian kernel and the correlation functions of the
 * interfaces, estimated from the simulation.
 *
 * The interfaces start independent, one on each bond with probability
 * nu'/(N+nu') (conditioned to an even number on the ring), and evolve with the
 * Glauber rates of intensity nu: hop 1 per direction, annihilation 1+gamma,
 * creation 1-gamma, gamma = (N^2-nu^2)/(N^2+nu^2).  The time is the macroscopic
 * one, tau = micro/N^2.  nu' = nu is the invariant measure and nu' = N is the
 * start from independent spins.
 *
 * With s_x the spins and eta_x = (1 - s_x s_{x+1})/2 the interface on the bond
 * x, the program estimates at every time of --T, for the distances r >= 1,
 *
 *   the four entries of the rescaled kernel of Garrod
 *       K_11 = -1/2   E[ s_x s_{x+r} ]
 *       K_12 =  N     E[ s_x s_{x+r} eta_{x+r} ]
 *       K_21 =  N     E[ eta_x s_x s_{x+r} ]
 *       K_22 = -2 N^2 E[ eta_x eta_{x+r} s_x s_{x+r} ]
 *   the two-point function    rho_2 = N^2 E[ eta_x eta_{x+r} ]
 *   the three-point function  rho_3 = N^3 E[ eta_x eta_{x+r} eta_{x+r+r'} ],
 *
 * for every pair of times sigma < tau of --T the two-time function
 *       rho(sigma, 0; tau, r) = N^2 E[ eta_x(sigma) eta_{x+r}(tau) ],
 * and the density rho_1 = N E[ eta_x ] on a grid of times.  The expectations
 * are averages over the translations x and over the simulations; the standard
 * errors are computed over the simulations, which are independent.
 *
 * The distances are grouped in bins of --bin bonds: the bin j contains
 * r = j*bin+1, ..., (j+1)*bin, and the estimate is the average over the bin
 * (for the two-time function over +r and -r, and r = 0 is given apart).
 * The exact values and the limit are computed by python_scripts/pfaffian_kernel.py,
 * and the comparison is drawn by
 *
 *   python3 python_scripts/plot_pfaffian_comparison.py <out>
 *
 * Usage:
 *   ./Pfaffian [--N 100] [--L 8] [--nu 1.0] [--nuprime 3.0]
 *              [--T 0.01,0.03,0.1] [--R 1.5] [--bin 5] [--ndens 40]
 *              [--nsims 2000] [--seed 0] [--out pfaffian_mc]
 *
 * --R is the largest distance, in macroscopic units (smaller than L/2); the
 * two gaps of the three-point function go up to R/2.  --ndens is the number of
 * steps of the grid of the density, from 0 to the largest time of --T.  The
 * times of --T are at most 16, and they enter the names of the files written
 * with %g, so they have to differ in their first six digits.
 *
 * rho_2, rho_3, K_22 and the two-time function are counts of pairs and of
 * triples of interfaces.  With few of them in a bin (some tens) the standard
 * error of the bin is not reliable: take wider bins or more simulations.
 *
 * Output, in the directory --out:
 *   density.txt                  tau, rho_1, stderr
 *   kernel_T<tau>.txt            for each bin: rmin, rmax, d, then K_11, K_12,
 *                                K_21, K_22, rho_2, each with its stderr
 *   rho3_T<tau>.txt              for each pair of bins: the two gaps, rho_3, stderr
 *   twotime_S<sigma>_T<tau>.txt  for r = 0 and for each bin: rho, stderr
 *   configuration.json           the parameters
 *
 * OpenMP: with -fopenmp the simulations run in parallel.  The seeds of a
 * simulation depend only on --seed and on its index, and the sums over the
 * simulations are sums of integers, so the files do not depend on the number
 * of threads.  Without --seed one is drawn and written in configuration.json.
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

#define MAX_T 16

/* The observables at one time, as functions of the bin of the distance. */
enum { Q_K11, Q_K12, Q_K21, Q_K22, Q_RHO2, NQ };

/* splitmix64: a hash of 64 bits. */
static unsigned long long mix64(unsigned long long z)
{
    z += 0x9E3779B97F4A7C15ULL;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

/* The seeds of the simulations.  The generator of the GSL takes 32 bits of
 * the seed.  If the seed of a simulation were a hash of (seed, simulation) cut
 * to 32 bits, two simulations of a long run would get the same one (with
 * 100000 simulations, about one pair) and would start from the same
 * configuration.  So the seeds of a run are the images of 0, 1, 2, ... under a
 * permutation of the integers of 32 bits, which depends on --seed: they are
 * all different.  The permutation is a Feistel network of four rounds on the
 * two halves of 16 bits, and it is a permutation whatever the function f is. */
static unsigned long permute32(unsigned long long key, unsigned long x)
{
    unsigned long left = (x >> 16) & 0xffffUL, right = x & 0xffffUL;
    for (int round = 1; round <= 4; round++) {
        unsigned long f = (unsigned long)(mix64(key ^ ((unsigned long long)round << 32) ^ right) & 0xffffULL);
        unsigned long next = left ^ f;
        left  = right;
        right = next;
    }
    return (left << 16) | right;
}

/* The numbers of the command line: the whole word must be a number, so that
 * "1e5" for an integer, or "0.01;0.03", are refused and not cut. */
static int read_int(const char *flag, const char *word, int *value)
{
    char *end;
    errno = 0;
    long v = strtol(word, &end, 10);
    if (end == word || *end != '\0' || errno || v > 2000000000L || v < -2000000000L) {
        fprintf(stderr, "%s: \"%s\" is not an integer.\n", flag, word);
        return 0;
    }
    *value = (int)v;
    return 1;
}

static int read_double(const char *flag, const char *word, double *value)
{
    char *end;
    *value = strtod(word, &end);
    if (end == word || *end != '\0') {
        fprintf(stderr, "%s: \"%s\" is not a number.\n", flag, word);
        return 0;
    }
    return 1;
}

static int read_seed(const char *flag, const char *word, unsigned long long *value)
{
    char *end;
    errno = 0;
    *value = strtoull(word, &end, 10);
    if (end == word || *end != '\0' || errno || word[0] == '-') {
        fprintf(stderr, "%s: \"%s\" is not a non-negative integer.\n", flag, word);
        return 0;
    }
    return 1;
}

/* The seed of the generator number x of the run.  The GSL takes the seed 0
 * for 4357, which another generator of the run may have: instead of 0 we give
 * the image of the last integer, that no generator of the run uses. */
static unsigned long seed_of(unsigned long long key, unsigned long x)
{
    unsigned long seed = permute32(key, x);
    return seed ? seed : permute32(key, 0xffffffffUL);
}

/* The name of a file contains its time written with %g.  Returns 1 if the two
 * times have the same name. */
static int same_name(double a, double b)
{
    char sa[64], sb[64];
    snprintf(sa, sizeof(sa), "%g", a);
    snprintf(sb, sizeof(sb), "%g", b);
    return strcmp(sa, sb) == 0;
}

/* All the estimates of one simulation are the entries of one vector, so that
 * the sums and the sums of the squares over the simulations are taken in one
 * place.  The layout says where each observable is.
 *
 * Every entry is an integer, a count or a sum of +1 and -1, and the estimate
 * is this integer times a factor that depends only on the entry (set_scale).
 * So the sums over the simulations are sums of integers: they are exact, and
 * they do not depend on the order in which the threads add them. */
typedef struct {
    int NT;        /* number of times of --T                      */
    int nb;        /* number of bins of the distance              */
    int nb3;       /* number of bins of each gap of rho_3         */
    int npairs;    /* number of pairs of times                    */
    int ndens;     /* number of steps of the grid of the density  */
    int one, three, two, dens, densT, total;   /* offsets and total length */
} Layout;

static int at_one(const Layout *y, int t, int q, int j)   { return y->one   + (t * NQ + q) * y->nb + j; }
static int at_three(const Layout *y, int t, int i, int j) { return y->three + (t * y->nb3 + i) * y->nb3 + j; }
static int at_two(const Layout *y, int q, int j)          { return y->two   + q * (y->nb + 1) + j; }

/* The factor of each entry: the estimate is the integer times its factor. */
static void set_scale(const Layout *y, int l, int N, int bin, double *scale)
{
    const double norm = 1.0 / ((double)l * bin);    /* average over x and over the bin */
    const double dN   = (double)N;
    static const double one[NQ] = { -0.5, 1.0, 1.0, -2.0, 1.0 };   /* K_11, K_12, K_21, K_22, rho_2 */
    static const int power[NQ]  = { 0, 1, 1, 2, 2 };                /* the power of N of each     */
    for (int t = 0; t < y->NT; t++) {
        for (int q = 0; q < NQ; q++)
            for (int j = 0; j < y->nb; j++)
                scale[at_one(y, t, q, j)] = one[q] * pow(dN, power[q]) * norm;
        for (int i = 0; i < y->nb3 * y->nb3; i++)
            scale[at_three(y, t, 0, 0) + i] = dN * dN * dN * norm / bin;
        scale[y->densT + t] = dN / l;
    }
    for (int q = 0; q < y->npairs; q++) {
        scale[at_two(y, q, 0)] = dN * dN / l;                       /* r = 0             */
        for (int j = 1; j <= y->nb; j++)
            scale[at_two(y, q, j)] = dN * dN / ((double)l * 2.0 * bin);   /* +r and -r  */
    }
    for (int g = 0; g <= y->ndens; g++)
        scale[y->dens + g] = dN / l;
}

static void layout_init(Layout *y, int NT, int nb, int ndens)
{
    y->NT = NT;  y->nb = nb;  y->nb3 = nb / 2;  y->ndens = ndens;
    y->npairs = NT * (NT - 1) / 2;
    y->one   = 0;
    y->three = y->one   + NT * NQ * nb;
    y->two   = y->three + NT * y->nb3 * y->nb3;
    y->dens  = y->two   + y->npairs * (nb + 1);
    y->densT = y->dens  + ndens + 1;
    y->total = y->densT + NT;
}

/* The interfaces of a configuration: the bonds x with s_x != s_{x+1}, in
 * increasing order.  Returns their number. */
static int find_walls(const int *s, int l, int *walls)
{
    int k = 0;
    for (int x = 0; x < l; x++)
        if (s[x] != s[x + 1 == l ? 0 : x + 1]) walls[k++] = x;
    return k;
}

/* The observables at one time: the contributions of the configuration s, with
 * its k interfaces, are added to val.  They are integers; the factors are in
 * set_scale. */
static void measure_one_time(const Layout *y, int t, const int *s, int l, int bin,
                             const int *walls, int k, long long *val)
{
    const int rmax  = y->nb * bin;        /* the largest distance          */
    const int rmax3 = y->nb3 * bin;       /* the largest gap of rho_3      */

    /* K_11 = -1/2 E[s_x s_{x+r}] */
    for (int r = 1; r <= rmax; r++) {
        long long sum = 0;
        for (int x = 0; x < l; x++) {
            int z = x + r; if (z >= l) z -= l;
            sum += s[x] * s[z];
        }
        val[at_one(y, t, Q_K11, (r - 1) / bin)] += sum;
    }

    /* K_12 and K_21: an interface on the bond x and the spin at distance r,
     * to the left for K_12 (the interface is the second point) and to the
     * right for K_21 (the interface is the first point). */
    for (int i = 0; i < k; i++) {
        int x = walls[i];
        for (int r = 1; r <= rmax; r++) {
            int right = x + r; if (right >= l) right -= l;
            int left  = x - r; if (left < 0)   left  += l;
            int j = (r - 1) / bin;
            val[at_one(y, t, Q_K12, j)] += s[left] * s[x];
            val[at_one(y, t, Q_K21, j)] += s[x] * s[right];
        }
    }

    /* K_22, rho_2 and rho_3: the pairs and the triples of interfaces.  For
     * each interface the following ones are visited in order, until the
     * distance exceeds rmax; rmax < l/2, so each pair is met once. */
    for (int i = 0; i < k; i++) {
        for (int m = 1; m < k; m++) {
            int j = i + m; if (j >= k) j -= k;
            int r = walls[j] - walls[i]; if (r < 0) r += l;
            if (r > rmax) break;
            int b1 = (r - 1) / bin;
            val[at_one(y, t, Q_RHO2, b1)] += 1;
            val[at_one(y, t, Q_K22, b1)]  += s[walls[i]] * s[walls[j]];

            if (r > rmax3) continue;
            for (int m2 = 1; m + m2 < k; m2++) {
                int jj = j + m2; if (jj >= k) jj -= k;
                int r2 = walls[jj] - walls[j]; if (r2 < 0) r2 += l;
                if (r2 > rmax3) break;
                val[at_three(y, t, b1, (r2 - 1) / bin)] += 1;
            }
        }
    }
}

/* The two-time function: the interfaces wa at the time sigma and wb at the
 * time tau > sigma.  The entry 0 is the displacement r = 0, the entry j+1 is
 * the bin j of |r|. */
static void measure_two_times(const Layout *y, int q, int l, int bin,
                              const int *wa, int ka, const int *wb, int kb, long long *val)
{
    const int rmax = y->nb * bin;
    for (int i = 0; i < ka; i++) {
        for (int j = 0; j < kb; j++) {
            int r = wb[j] - wa[i];
            if (r < 0) r = -r;
            if (r > l - r) r = l - r;               /* the distance on the ring */
            if (r == 0)
                val[at_two(y, q, 0)] += 1;
            else if (r <= rmax)
                val[at_two(y, q, 1 + (r - 1) / bin)] += 1;
        }
    }
}

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, z = *(const double *)b;
    return (x > z) - (x < z);
}

/* A double written with the fewest digits that give it back when it is read:
 * the scripts rebuild from configuration.json the names of the files and the
 * exact parameters. */
static void print_double(FILE *f, double x)
{
    char s[40];
    for (int digits = 15; digits <= 17; digits++) {
        snprintf(s, sizeof(s), "%.*g", digits, x);
        if (strtod(s, NULL) == x) break;
    }
    fputs(s, f);
}

/* The integer sums, their squares and the factors, to write the output. */
typedef struct {
    const long long          *sum;
    const unsigned long long *sq;
    const double             *scale;
    int                       n;       /* number of simulations */
} Totals;

/* mean and standard error of an entry, from the sum and the sum of the squares */
static void mean_err(const Totals *T, int i, double *mean, double *err)
{
    double m   = (double)T->sum[i] / T->n;
    double var = ((double)T->sq[i] - T->n * m * m) / (T->n - 1);     /* n >= 2 */
    *mean = T->scale[i] * m;
    *err  = fabs(T->scale[i]) * sqrt(fmax(var, 0.0) / T->n);
}

int main(int argc, char **argv)
{
    /* ---- defaults ---- */
    int    N = 100, L = 8, nsims = 2000, bin = 5, ndens = 40;
    double nu = 1.0, nuprime = 3.0, R = 1.5;
    double T_values[MAX_T] = { 0.01, 0.03, 0.1 };
    int    NT = 3;
    unsigned long long seed = 0;
    const char *out = "pfaffian_mc";

    /* ---- flag parsing ---- */
    for (int i = 1; i < argc; i++) {
        const char *flag = argv[i];
        int ok = 1;
        if (i + 1 >= argc) { fprintf(stderr, "Unknown/incomplete flag: %s\n", flag); return 1; }
        char *word = argv[++i];
        if      (!strcmp(flag, "--N"))       ok = read_int(flag, word, &N);
        else if (!strcmp(flag, "--L"))       ok = read_int(flag, word, &L);
        else if (!strcmp(flag, "--nu"))      ok = read_double(flag, word, &nu);
        else if (!strcmp(flag, "--nuprime")) ok = read_double(flag, word, &nuprime);
        else if (!strcmp(flag, "--R"))       ok = read_double(flag, word, &R);
        else if (!strcmp(flag, "--bin"))     ok = read_int(flag, word, &bin);
        else if (!strcmp(flag, "--ndens"))   ok = read_int(flag, word, &ndens);
        else if (!strcmp(flag, "--nsims"))   ok = read_int(flag, word, &nsims);
        else if (!strcmp(flag, "--seed"))    ok = read_seed(flag, word, &seed);
        else if (!strcmp(flag, "--out"))     out = word;
        else if (!strcmp(flag, "--T")) {
            NT = 0;
            for (char *tok = strtok(word, ","); tok && ok; tok = strtok(NULL, ",")) {
                if (NT == MAX_T) { fprintf(stderr, "--T takes at most %d times.\n", MAX_T); return 1; }
                ok = read_double(flag, tok, &T_values[NT++]);
            }
        }
        else { fprintf(stderr, "Unknown/incomplete flag: %s\n", flag); return 1; }
        if (!ok) return 1;
    }
    /* The comparisons are written so that a NaN is refused too. */
    if (N <= 0 || L <= 0 || nsims <= 0 || bin <= 0 || ndens <= 0 || NT <= 0 ||
        !(nu > 0 && isfinite(nu)) || !(nuprime > 0 && isfinite(nuprime)) || !(R > 0 && isfinite(R))) {
        fprintf(stderr, "Invalid parameters.\n"); return 1;
    }
    if (nsims < 2) {
        fprintf(stderr, "--nsims must be at least 2: the standard errors are computed over the simulations.\n");
        return 1;
    }
    if ((double)N * (double)L > 1e9) { fprintf(stderr, "N*L is too large.\n"); return 1; }
    qsort(T_values, (size_t)NT, sizeof(double), cmp_double);
    for (int t = 0; t < NT; t++) {
        if (!(T_values[t] > 0 && isfinite(T_values[t]))) {
            fprintf(stderr, "The times of --T must be positive and finite.\n"); return 1;
        }
        if (t > 0 && same_name(T_values[t - 1], T_values[t])) {
            fprintf(stderr, "Two times of --T are written in the same way by %%g (%g): "
                            "their files would have the same name.\n", T_values[t]);
            return 1;
        }
    }
    if (!seed) seed = (unsigned long long)time(NULL) ^ ((unsigned long long)getpid() << 20);

    /* ---- derived quantities ---- */
    int    l       = N * L;
    double gamma   = ((double)N * N - nu * nu) / ((double)N * N + nu * nu);
    double ann     = 1.0 + gamma;              /* the Glauber rates, hop 1 per direction */
    double create  = 1.0 - gamma;
    double p       = nu / ((double)N + nu);            /* invariant measure */
    double p_init  = nuprime / ((double)N + nuprime);  /* initial measure   */
    if (!isfinite(gamma) || !(p < 1.0) || !(p_init < 1.0)) {      /* nu or nu' of the order of 1e16 N or more */
        fprintf(stderr, "--nu or --nuprime is too large with respect to N.\n"); return 1;
    }
    if (!(R * N / bin < (double)l)) {
        fprintf(stderr, "--R must be at least bin/N and smaller than L/2.\n"); return 1;
    }
    int    nb      = (int)floor(R * N / bin + 1e-9);
    if (nb < 1 || 2 * nb * bin >= l) {
        fprintf(stderr, "--R must be at least bin/N and smaller than L/2.\n"); return 1;
    }
    double Tmax = T_values[NT - 1];

    Layout y;
    layout_init(&y, NT, nb, ndens);

    printf("=== The kernel and the correlation functions of the interfaces ===\n");
    printf("  N=%d L=%d l=%d  nu=%g (p=%.6f)  nu'=%g (p=%.6f)\n", N, L, l, nu, p, nuprime, p_init);
    printf("  ann=%.10g create=%.10g  (Glauber rates, hop 1 per direction)\n", ann, create);
    printf("  times:");
    for (int t = 0; t < NT; t++) printf(" %g", T_values[t]);
    printf("   distances up to %g in %d bins of %d bonds\n", (double)nb * bin / N, nb, bin);
    printf("  nsims=%d  seed=%llu\n", nsims, seed);
#ifdef _OPENMP
    printf("  OpenMP: %d threads\n", omp_get_max_threads());
#endif
    fflush(stdout);

    /* The output directory is checked before the simulations, not after:
     * configuration.json is written now. */
    struct stat st;
    if (mkdir(out, 0755) && errno != EEXIST) { perror("mkdir"); return 1; }
    if (stat(out, &st) || !S_ISDIR(st.st_mode)) {
        fprintf(stderr, "%s is not a directory.\n", out); return 1;
    }
    const size_t path_len = strlen(out) + 96;        /* the directory and the longest name of a file */
    char *path = malloc(path_len);
    if (!path) { perror("malloc"); return 1; }
    {
        snprintf(path, path_len, "%s/configuration.json", out);
        FILE *f = fopen(path, "w");
        if (!f) { perror("configuration.json"); return 1; }
        fprintf(f, "{\n");
        fprintf(f, "  \"N\": %d,\n  \"L\": %d,\n  \"l\": %d,\n", N, L, l);
        fprintf(f, "  \"nu\": ");         print_double(f, nu);
        fprintf(f, ",\n  \"nu_prime\": "); print_double(f, nuprime);
        fprintf(f, ",\n");
        fprintf(f, "  \"gamma\": %.12f,\n  \"p_equilibrium\": %.12f,\n  \"p_initial\": %.12f,\n", gamma, p, p_init);
        fprintf(f, "  \"ann\": %.12f,\n  \"create\": %.12f,\n", ann, create);
        fprintf(f, "  \"T_values\": [");
        for (int t = 0; t < NT; t++) { fprintf(f, "%s", t ? ", " : ""); print_double(f, T_values[t]); }
        fprintf(f, "],\n");
        fprintf(f, "  \"bin\": %d,\n  \"n_bins\": %d,\n  \"n_bins_rho3\": %d,\n  \"ndens\": %d,\n", bin, nb, y.nb3, ndens);
        fprintf(f, "  \"nsims\": %d,\n  \"engine\": \"C\",\n", nsims);
        fprintf(f, "  \"time_scale_formula\": %.1f,\n", 2.0 * N * N);
        fprintf(f, "  \"seed\": %llu\n", seed);
        fprintf(f, "}\n");
        fclose(f);
    }

    const unsigned long long key = mix64(seed);     /* the key of the permutation of the seeds */

    /* The sums of the squares must stay below 2^64: the largest entry of a
     * simulation is l*bin, for K_11. */
    if ((double)l * bin * (double)l * bin * (double)nsims > 1.8e19) {
        fprintf(stderr, "Too many simulations for this size: (N L bin)^2 nsims must be below 1.8e19.\n"); return 1;
    }
    long long          *sum   = calloc((size_t)y.total, sizeof(long long));
    unsigned long long *sq    = calloc((size_t)y.total, sizeof(unsigned long long));
    double             *scale = calloc((size_t)y.total, sizeof(double));
    if (!sum || !sq || !scale) { perror("calloc"); return 1; }
    set_scale(&y, l, N, bin, scale);

    int done = 0;
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        /* Per-thread buffers and sums. */
        int       *spins  = malloc((size_t)l * sizeof(int));
        int       *walls  = malloc((size_t)NT * l * sizeof(int));   /* the interfaces at each time */
        int       *nwalls = malloc((size_t)NT * sizeof(int));
        long long *val    = malloc((size_t)y.total * sizeof(long long));
        long long          *sum_t = calloc((size_t)y.total, sizeof(long long));
        unsigned long long *sq_t  = calloc((size_t)y.total, sizeof(unsigned long long));
        gsl_rng *ic_rng = gsl_rng_alloc(gsl_rng_mt19937);
        if (!spins || !walls || !nwalls || !val || !sum_t || !sq_t) { perror("malloc"); exit(1); }

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
        for (int sim = 0; sim < nsims; sim++) {
            memset(val, 0, (size_t)y.total * sizeof(long long));

            /* The initial configuration and the trajectory have their own seeds,
             * different from the ones of all the other simulations of the run. */
            gsl_rng_set(ic_rng, seed_of(key, 2UL * (unsigned long)sim));
            unsigned long traj_seed = seed_of(key, 2UL * (unsigned long)sim + 1UL);

            sample_equilibrium_interfaces(spins, l, p_init, ic_rng);
            InterfaceState *gs = istate_alloc(spins, l, N, ann, create, traj_seed);
            val[y.dens] = istate_n_interfaces(gs);

            /* The times of the grid of the density and the times of --T are
             * visited in order. */
            int g = 1, t = 0;
            while (g <= ndens || t < NT) {
                double tg = (g <= ndens) ? Tmax * g / ndens : INFINITY;
                double tt = (t < NT) ? T_values[t] : INFINITY;
                if (tg <= tt) {
                    istate_advance_to(gs, tg);
                    val[y.dens + g] = istate_n_interfaces(gs);
                    g++;
                } else {
                    istate_advance_to(gs, tt);
                    istate_to_spins(gs, spins);
                    int *w = walls + (size_t)t * l;
                    nwalls[t] = find_walls(spins, l, w);
                    val[y.densT + t] = nwalls[t];
                    measure_one_time(&y, t, spins, l, bin, w, nwalls[t], val);
                    t++;
                }
            }
            istate_free(gs);

            for (int a = 0, q = 0; a < NT; a++)
                for (int b = a + 1; b < NT; b++, q++)
                    measure_two_times(&y, q, l, bin, walls + (size_t)a * l, nwalls[a],
                                      walls + (size_t)b * l, nwalls[b], val);

            for (int i = 0; i < y.total; i++) {
                sum_t[i] += val[i];
                sq_t[i]  += (unsigned long long)(val[i] * val[i]);
            }

#ifdef _OPENMP
#pragma omp critical
#endif
            {
                done++;
                if (done % (nsims >= 10 ? nsims / 10 : 1) == 0)
                    printf("  sim %d/%d\n", done, nsims), fflush(stdout);
            }
        }

        /* Sums of integers: the order in which the threads arrive does not matter. */
#ifdef _OPENMP
#pragma omp critical
#endif
        for (int i = 0; i < y.total; i++) { sum[i] += sum_t[i]; sq[i] += sq_t[i]; }

        gsl_rng_free(ic_rng);
        free(spins); free(walls); free(nwalls); free(val); free(sum_t); free(sq_t);
    }

    /* ---- output ---- */
    double mean, err;
    FILE *f;
    const Totals tot = { sum, sq, scale, nsims };

    snprintf(path, path_len, "%s/density.txt", out);
    f = fopen(path, "w");
    if (!f) { perror("density.txt"); return 1; }
    fprintf(f, "# tau\trho_1\tstderr\n");
    for (int g = 0; g <= ndens; g++) {
        mean_err(&tot, y.dens + g, &mean, &err);
        fprintf(f, "%.10g\t%.8f\t%.8f\n", Tmax * g / ndens, mean, err);
    }
    fclose(f);

    for (int t = 0; t < NT; t++) {
        static const char *names[NQ] = { "K_11", "K_12", "K_21", "K_22", "rho_2" };
        snprintf(path, path_len, "%s/kernel_T%g.txt", out, T_values[t]);
        f = fopen(path, "w");
        if (!f) { perror("kernel"); return 1; }
        mean_err(&tot, y.densT + t, &mean, &err);
        fprintf(f, "# tau = %g   rho_1 = %.8f +- %.8f\n", T_values[t], mean, err);
        fprintf(f, "# rmin\trmax\td");
        for (int q = 0; q < NQ; q++) fprintf(f, "\t%s\tstderr", names[q]);
        fprintf(f, "\n");
        for (int j = 0; j < nb; j++) {
            fprintf(f, "%d\t%d\t%.6f", j * bin + 1, (j + 1) * bin, (j * bin + 0.5 * (bin + 1)) / N);
            for (int q = 0; q < NQ; q++) {
                mean_err(&tot, at_one(&y, t, q, j), &mean, &err);
                fprintf(f, "\t%.8f\t%.8f", mean, err);
            }
            fprintf(f, "\n");
        }
        fclose(f);

        snprintf(path, path_len, "%s/rho3_T%g.txt", out, T_values[t]);
        f = fopen(path, "w");
        if (!f) { perror("rho3"); return 1; }
        fprintf(f, "# tau = %g   the bin j is r = j*%d+1, ..., (j+1)*%d\n", T_values[t], bin, bin);
        fprintf(f, "# bin1\tbin2\td1\td2\trho_3\tstderr\n");
        for (int i = 0; i < y.nb3; i++)
            for (int j = 0; j < y.nb3; j++) {
                mean_err(&tot, at_three(&y, t, i, j), &mean, &err);
                fprintf(f, "%d\t%d\t%.6f\t%.6f\t%.8f\t%.8f\n", i, j,
                        (i * bin + 0.5 * (bin + 1)) / N, (j * bin + 0.5 * (bin + 1)) / N, mean, err);
            }
        fclose(f);
    }

    for (int a = 0, q = 0; a < NT; a++)
        for (int b = a + 1; b < NT; b++, q++) {
            snprintf(path, path_len, "%s/twotime_S%g_T%g.txt", out, T_values[a], T_values[b]);
            f = fopen(path, "w");
            if (!f) { perror("twotime"); return 1; }
            fprintf(f, "# sigma = %g   tau = %g\n", T_values[a], T_values[b]);
            fprintf(f, "# rmin\trmax\td\trho\tstderr\n");
            mean_err(&tot, at_two(&y, q, 0), &mean, &err);
            fprintf(f, "0\t0\t%.6f\t%.8f\t%.8f\n", 0.0, mean, err);
            for (int j = 0; j < nb; j++) {
                mean_err(&tot, at_two(&y, q, j + 1), &mean, &err);
                fprintf(f, "%d\t%d\t%.6f\t%.8f\t%.8f\n", j * bin + 1, (j + 1) * bin,
                        (j * bin + 0.5 * (bin + 1)) / N, mean, err);
            }
            fclose(f);
        }

    printf("Wrote %s/{density.txt,kernel_T*.txt,rho3_T*.txt,twotime_S*_T*.txt,configuration.json}\n", out);
    free(sum); free(sq); free(scale); free(path);
    return 0;
}
