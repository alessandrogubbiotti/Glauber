#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <ctype.h>

#include "configuration.h"
#include "configuration_statistics.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Continuous-time Glauber Dynamics on the torus with parity update.
// Rates: interface spins flip at rate 1; bulk spins flip at annihilation/creation rate.
// The parity split (odd/even sites) avoids spurious interface creation at synchronous updates.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Parameters passed to the Glauber update function.
typedef struct {
    int l;
    double prob_annihilation;
    double prob_creation;
    double prob_movement;
    int n_micro_evolution;
    gsl_rng *rng;
} Parameters_update;

void free_parameters_update(Parameters_update *par) {
    gsl_rng_free(par->rng);
    free(par);
}

/* Seed source.  A fixed seed via $GLAUBER_SEED makes a run reproducible (used
 * when this engine is the reference implementation in a comparison); otherwise
 * draw a fresh one, so that repeated runs give INDEPENDENT statistics.  Same
 * convention as $GILLESPIE_SEED in gillespie_main.c. */
static unsigned long read_base_seed(void) {
    const char *env = getenv("GLAUBER_SEED");
    if (env && *env) {
        unsigned long s = strtoul(env, NULL, 10);
        if (s) return s;
    }
    unsigned long seed = 0;
    FILE *f = fopen("/dev/urandom", "rb");
    if (f) {
        if (fread(&seed, sizeof(seed), 1, f) != 1) seed = 0;
        fclose(f);
    }
    if (!seed) seed = (unsigned long)time(NULL);
    return seed;
}

void set_parameters(const ModelConfig *conf, Parameters_update *parameters) {
    parameters->l = conf->L * conf->N;

    /* The parity update splits the ring into the even and the odd sublattice,
     * so that no two ADJACENT sites are ever updated in the same half-pass.
     * On an odd ring sites l-1 and 0 share a parity while being neighbours,
     * which re-opens the simultaneous-adjacent-flip hole the scheme exists to
     * close (it would create/destroy interface pairs spuriously). */
    if (parameters->l % 2 != 0) {
        fprintf(stderr,
                "set_parameters: L*N = %d is ODD; the parity update requires an "
                "even ring (choose L or N even).\n", parameters->l);
        exit(EXIT_FAILURE);
    }

    parameters->prob_movement    = 1 - exp(-1.0 / conf->Micro_n_steps);
    parameters->prob_annihilation = 1 - exp(-(double)conf->annihilation / conf->Micro_n_steps);
    parameters->prob_creation     = 1 - exp(-(double)conf->creation     / conf->Micro_n_steps);
    // Factor 2: only half the sites are updated per parity step.
    parameters->n_micro_evolution = 2 * conf->Micro_n_steps * conf->N * conf->N / conf->resolution;

    printf("DEBUG: creation rate=%.6f, Micro_n_steps=%d, prob_creation=%.6f\n",
           conf->creation, conf->Micro_n_steps, parameters->prob_creation);
    printf("DEBUG: annihilation rate=%lf, Micro_n_steps=%d, prob_annihilation=%.6f\n",
           conf->annihilation, conf->Micro_n_steps, parameters->prob_annihilation);
    printf("DEBUG: N^2 * prob_creation=%.6f, prob_annihilation=%.6f\n",
           (double)conf->N * (double)conf->N * parameters->prob_creation, parameters->prob_annihilation);

    const gsl_rng_type *G;
    gsl_rng_env_setup();
    G = gsl_rng_mt19937;
    parameters->rng = gsl_rng_alloc(G);

    /* Seed BOTH generators explicitly: gsl_rng drives the dynamics, and C
     * rand() drives the spin initialisation in initialize_spins_*.  Without
     * this both start from their fixed default seeds, so every run of this
     * binary reproduces the previous one bit for bit -- repeated runs would
     * add no statistics at all.  The seed is printed so a run can be repeated. */
    unsigned long seed = read_base_seed();
    gsl_rng_set(parameters->rng, seed);
    srand((unsigned int)(seed & 0xFFFFFFFFu));
    printf("Seed: %lu   (set $GLAUBER_SEED to reproduce this run)\n", seed);
}


// Output helpers

void write_stats_header(FILE *f, int N_sim, const char **stat_names, int N_stats) {
    fprintf(f, "time");
    for (int i = 0; i < N_sim; i++)
        for (int j = 0; j < N_stats; j++)
            fprintf(f, "\tsim%d_%s", i, stat_names[j]);
    fprintf(f, "\n");
}

void write_stats_row(FILE *f, double *corr, double *magnetization, double *time_delayed_correlations, int *interfaces, int N_sim) {
    for (int i = 0; i < N_sim; i++) {
        fprintf(f, "%.6f\t", corr[i]);
        fprintf(f, "%.6f\t", magnetization[i]);
        fprintf(f, "%.6f\t", time_delayed_correlations[i]);
        fprintf(f, "%d\t",   interfaces[i]);
    }
    fprintf(f, "\n");
}

void print_spins_to_file_binary2(int *spins, FILE *f, int length) {
    fwrite(spins, sizeof(int), length, f);
}


// Statistics

double correlation(int *spins1, int *spins2, Two_Intervals *intervals) {
    assert(intervals != NULL);
    double corr = 0;
    for (int i = intervals->a1; i < intervals->b1; i++)
        for (int j = intervals->a2; j < intervals->b2; j++)
            corr += spins1[i] * spins2[j];
    return corr / (double)((intervals->b1 - intervals->a1) * (intervals->b2 - intervals->a2));
}

double compute_magnetization(int *spin, Interval *interval) {
    double sum = 0;
    for (int i = interval->a1; i < interval->b1; i++)
        sum += spin[i];
    return sum / (double)(interval->b1 - interval->a1);
}

int n_interfaces(int *spins, int l) {
    int sum = (spins[0] != spins[l - 1]) ? 1 : 0;
    for (int i = 0; i < l - 1; i++)
        if (spins[i] != spins[i + 1]) sum++;
    return sum;
}

int creation_of_interfaces(int a1, int b1, int l, int *spin1, int *spin2) {
    int count = 0;
    for (int i = a1; i < b1; i++) {
        if (spin1[i] != spin1[(i + 1 + l) % l]) count++;
        if (spin2[i] != spin2[(i + 1 + l) % l]) count--;
    }
    return count;
}

double F(double x, int N) {
    return sin(4 * M_PI * x / N);
}

double GF(int l, int N, int *spins) {
    double integral = 0;
    for (int i = 0; i < l; i++)
        integral += spins[i] * F(i / N, N);
    return integral / (double)N;
}


// Spin initialization

void initialize_spins_1(int *spins, int l, void *params) {
    spins[0] = (rand() % 2) * 2 - 1;
    double *p = (double *)params;
    for (int i = 1; i < l; i++)
        spins[i] = ((rand() / (double)RAND_MAX) < *p) ? -spins[i - 1] : spins[i - 1];
}

void initialize_spins_2(int *spins, int l, void *params) {
    (void)params;
    for (int i = 0; i < l;     i++) spins[i] =  1;
    for (int i = l / 2; i < l; i++) spins[i] = -1;
}

void initialize_spins_3(int *spins, int l, void *params) {
    (void)params;
    for (int i = 0; i < l; i++) spins[i] = -1;
}

void initialize_spins_4(int *spins, int l, void *params) {
    double *p = (double *)params;
    for (int i = 0; i < l; i++)
        spins[i] = ((rand() / (double)RAND_MAX) < *p) ? 1 : -1;
}

void initialize_spins_5(int *spins, int l, void *params) {
    (void)params;
    for (int i = 0; i < l; i++)
        spins[i] = (i % 2 == 0) ? 1 : -1;
}


// Glauber update (parity). a=0: update even sites; a=1: update odd sites.
void Glauber_update_parity(int *spins_aux, int *new_spins, Parameters_update *params, int a) {
    int left, right, si;
    double p;
    for (int s = 0; s < params->n_micro_evolution; s++) {
        memcpy(spins_aux, new_spins, params->l * sizeof(int));

        // Boundary site: index 0 (a==1) or index l-1 (a==0).
        if (a == 1) {
            p     = gsl_rng_uniform(params->rng);
            left  = spins_aux[params->l - 1];
            right = spins_aux[1];
            si    = spins_aux[0];
            if (left == right) {
                if (left == si  && p < params->prob_creation)    new_spins[0] = -new_spins[0];
                if (left != si  && p < params->prob_annihilation) new_spins[0] = -new_spins[0];
            } else {
                if (p < params->prob_movement) new_spins[0] = -new_spins[0];
            }
        } else {
            p     = gsl_rng_uniform(params->rng);
            left  = spins_aux[params->l - 2];
            right = spins_aux[0];
            si    = spins_aux[params->l - 1];
            if (left == right) {
                if (left == si  && p < params->prob_creation)    new_spins[params->l - 1] = -new_spins[params->l - 1];
                if (left != si  && p < params->prob_annihilation) new_spins[params->l - 1] = -new_spins[params->l - 1];
            } else {
                if (p < params->prob_movement) new_spins[params->l - 1] = -new_spins[params->l - 1];
            }
        }

        // Interior sites (stride 2 for parity).
        for (int i = 1 + a; i < params->l - 1; i += 2) {
            p     = gsl_rng_uniform(params->rng);
            left  = spins_aux[i - 1];
            right = spins_aux[i + 1];
            si    = spins_aux[i];
            if (left == right) {
                if (left == si  && p < params->prob_creation)    new_spins[i] = -new_spins[i];
                if (left != si  && p < params->prob_annihilation) new_spins[i] = -new_spins[i];
            } else {
                if (p < params->prob_movement) new_spins[i] = -new_spins[i];
            }
        }
        a = 1 - a;
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// MAIN
//////////////////////////////////////////////////////////////////////////////////////

int main(void) {
    ModelConfig *model = malloc(sizeof(ModelConfig));
    if (!model) { perror("model malloc failed"); exit(EXIT_FAILURE); }

    StatsConfig *stats = malloc(sizeof(StatsConfig));
    if (!stats) { perror("stats malloc failed"); exit(EXIT_FAILURE); }

    read_model_config_file(model, "config.txt");
    read_stats_config_file(stats, "config.txt");
    validate_model_config(model);
    validate_stats_config(stats, model->N, model->L);
    print_model_config(model);

    char dirname[256];
    create_output_directory(model, dirname, sizeof(dirname));
    print_conf_json(model, dirname);

    int timestep_marker = -9999;

    Parameters_update *parameters = malloc(sizeof(Parameters_update));
    if (!parameters) { perror("parameters malloc failed"); exit(EXIT_FAILURE); }

    Two_Intervals *corr_par              = malloc(sizeof(Two_Intervals));
    Two_Intervals *time_delayed_corr_par = malloc(sizeof(Two_Intervals));
    Interval      *par_av_magnetization  = malloc(sizeof(Interval));
    if (!corr_par || !time_delayed_corr_par || !par_av_magnetization) {
        perror("interval malloc failed"); exit(EXIT_FAILURE);
    }

    set_parameters(model, parameters);
    set_interval_magnetization(stats, model->N, par_av_magnetization);
    set_correlation_params(stats,     model->N, corr_par);
    set_time_delayed_params(stats,    model->N, time_delayed_corr_par);

    double *corr             = malloc(model->N_simulations * sizeof(double));
    double *magnetization    = malloc(model->N_simulations * sizeof(double));
    double *time_delayed_corr = malloc(model->N_simulations * sizeof(double));
    int    *interfaces        = malloc(model->N_simulations * sizeof(int));

    int length = model->N * model->L * model->N_simulations;

    printf("Flipping probabilities per microscopic step:\n");
    printf("  interface spin:   %lf\n", parameters->prob_movement);
    printf("  creation:         %lf\n", parameters->prob_creation);
    printf("  annihilation:     %lf\n", parameters->prob_annihilation);
    printf("DEBUG: N=%d T=%lf corr_par=(%d,%d,%d,%d)\n",
           model->N, model->T, corr_par->a1, corr_par->a2, corr_par->b1, corr_par->b2);

    int *data      = malloc(length * sizeof(int));
    int *new_data  = malloc(length * sizeof(int));
    int *past_data = malloc(length * sizeof(int));
    int *data_aux  = malloc(length * sizeof(int));

    int **spins      = malloc(model->N_simulations * sizeof(int *));
    int **first_spins = malloc(model->N_simulations * sizeof(int *));
    int **spins_aux  = malloc(model->N_simulations * sizeof(int *));
    int **new_spins  = malloc(model->N_simulations * sizeof(int *));

    for (int i = 0; i < model->N_simulations; i++) {
        spins[i]       = &data[i       * parameters->l];
        first_spins[i] = &past_data[i  * parameters->l];
        spins_aux[i]   = &data_aux[i   * parameters->l];
        new_spins[i]   = &new_data[i   * parameters->l];
    }

    for (int i = 0; i < model->N_simulations; i++)
        model->initialize(spins[i], parameters->l, model->initialization_param);
    memcpy(past_data, data, length * sizeof(int));
    memcpy(new_data,  data, length * sizeof(int));

    char filepath[512];
    snprintf(filepath, sizeof(filepath), "%s/spins_output.bin", dirname);
    FILE *output = fopen(filepath, "wb");
    if (!output) { perror("File open failed"); exit(1); }

    char filepath2[512];
    snprintf(filepath2, sizeof(filepath2), "%s/Statistics.txt", dirname);
    FILE *Statistics = fopen(filepath2, "w");
    if (!Statistics) { perror("File open failed"); exit(1); }

    const char *stat_names[] = {"correlations", "magnetization", "time_delayed_corr", "interfaces"};
    write_stats_header(Statistics, model->N_simulations, stat_names, 4);

    int a = 0;
    int Macro_Times = (int)(model->T * model->resolution);

    for (int t = 0; t < Macro_Times; t++) {
        fprintf(Statistics, "%lf\t", (double)t / (double)model->resolution);
        for (int i = 0; i < model->N_simulations; i++) {
            Glauber_update_parity(spins_aux[i], new_spins[i], parameters, a);
            /* All four observables are read from spins[i] = the state at the
             * time written in this row's first column.  new_spins[i] already
             * holds the NEXT frame (the update above ran first), so measuring
             * the ageing correlator there would shift it by one frame,
             * i.e. report C(t + 1/resolution, 0) under the label t. */
            corr[i]              = correlation(spins[i], spins[i], corr_par);
            magnetization[i]     = compute_magnetization(spins[i], par_av_magnetization);
            time_delayed_corr[i] = correlation(first_spins[i], spins[i], time_delayed_corr_par);
            interfaces[i]        = n_interfaces(spins[i], parameters->l);
        }
        print_spins_to_file_binary2(data, output, length);
        fwrite(&timestep_marker, sizeof(int), 1, output);
        write_stats_row(Statistics, corr, magnetization, time_delayed_corr, interfaces, model->N_simulations);
        memcpy(data, new_data, length * sizeof(int));
        a = (a + 1) % 2;
    }

    fclose(output);
    fclose(Statistics);

    free(data); free(new_data); free(past_data); free(data_aux);
    free(spins); free(new_spins); free(spins_aux); free(first_spins);
    free(corr); free(magnetization); free(time_delayed_corr); free(interfaces);
    free_parameters_update(parameters);
    free(corr_par); free(time_delayed_corr_par); free(par_av_magnetization);
    if (model->initialization_param) free(model->initialization_param);
    free(model);
    free(stats);

    char command[600];
    snprintf(command, sizeof(command), "python3 python_scripts/Plot_Glauber_tidy_Binary.py \"%s\"", dirname);
    if (system(command) != 0)
        fprintf(stderr, "Failed to run Plot_Glauber_tidy_Binary.py on %s\n", dirname);

    char command2[600];
    snprintf(command2, sizeof(command2), "python3 python_scripts/Plot_stats.py ./%s/Statistics.txt --outdir ./%s/", dirname, dirname);
    if (system(command2) != 0)
        fprintf(stderr, "Failed to run Plot_stats.py on %s\n", dirname);

    return 0;
}


// Commented-out interactive configuration functions (replaced by file-based parsing):
//
// void read_config(Configuration *conf) { ... }          // stdin-based full config read
// void choose_rates(double *, double *, int) { ... }     // interactive rate selection menu
// void choose_initialization_function(Configuration *) { ... }  // interactive init choice
