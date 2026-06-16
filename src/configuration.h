#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <stdlib.h>
#include <stdio.h>

#define MAX_LINE 256

typedef void (*Initialization_function)(int *spin, int l, void *variables);

/* Rate parametrization mode (conventions: docs/section1.tex, Section 1).
 *
 *   glauber  : ann = 1 + gamma, create = 1 - gamma  (closes the hierarchy,
 *              ann + create = 2 = 2*hop_rate).            [default]
 *   scaled   : glauber rates multiplied by pow(N, a) (slowdown_exponent a);
 *              same invariant measure, slower reactions.
 *   explicit : ann, create read verbatim from the config file; the implied
 *              nu and gamma are reported.
 *
 * In every mode gamma = (N^2 - nu^2)/(N^2 + nu^2),  beta = log(N/nu)/2,
 * eta = (N - nu)/(N + nu),  p = nu/(N + nu),  and the engine's macroscopic
 * time is tau = micro_time / N^2, so the analytic (Glauber-time) formulas
 * are evaluated at t_formula = 2 * N^2 * tau  (time_scale_formula = 2*N^2). */
typedef enum {
    RATE_GLAUBER  = 0,
    RATE_SCALED   = 1,
    RATE_EXPLICIT = 2
} RateMode;

typedef struct {
    int N;
    double T;
    int L;
    double annihilation;
    double creation;
    int N_simulations;
    int resolution;
    int Micro_n_steps;
    Initialization_function initialize;
    void *initialization_param;

    /* --- Re-parametrization (Section 1 conventions) ------------------------ */
    double   nu;                 /* target Poisson intensity (default 1.0)     */
    RateMode rate_mode;          /* glauber | scaled | explicit (default glauber) */
    double   slowdown_exponent;  /* a, used by scaled mode (default 0.0)       */
    int      ann_set;            /* annihilation key present in config?        */
    int      create_set;         /* creation key present in config?            */

    /* --- Derived block (filled by model_apply_rate_mode) ------------------- */
    double beta;                 /* log(N/nu)/2                                */
    double gamma;                /* (N^2 - nu^2)/(N^2 + nu^2)                  */
    double eta;                  /* (N - nu)/(N + nu)                          */
    double p_equilibrium;        /* nu/(N + nu)  (= 1/(1+sqrt(ann/create)))    */
    double time_scale_formula;   /* 2*N^2  (t_formula = time_scale_formula*tau) */
} ModelConfig;

typedef int (*ModelConfigSetter)(ModelConfig *conf, const char *value);

typedef struct {
    const char *key;
    ModelConfigSetter setter;
} ModelConfigEntry;
/* Function names automatically decay to function pointers.
 * Equivalent to: conf->initialize = &initialize_spins_5;
 */
// Forward declarations for initialization functions defined in Glauber_tidy.c
void initialize_spins_1(int *spins, int l, void *params);
void initialize_spins_2(int *spins, int l, void *params);
void initialize_spins_3(int *spins, int l, void *params);
void initialize_spins_4(int *spins, int l, void *params);
void initialize_spins_5(int *spins, int l, void *params);
void initialize_spins_6(int *spins, int l, void *params);

// Setters
int set_N(ModelConfig *conf, const char *value);
int set_T(ModelConfig *conf, const char *value);
int set_L(ModelConfig *conf, const char *value);
int set_annihilation(ModelConfig *conf, const char *value);
int set_creation(ModelConfig *conf, const char *value);
int set_N_simulations(ModelConfig *conf, const char *value);
int set_resolution(ModelConfig *conf, const char *value);
int set_Micro_n_steps(ModelConfig *conf, const char *value);
int set_init_function(ModelConfig *conf, const char *value);
int set_init_param(ModelConfig *conf, const char *value);
int set_nu(ModelConfig *conf, const char *value);
int set_rate_mode(ModelConfig *conf, const char *value);
int set_slowdown_exponent(ModelConfig *conf, const char *value);

// Config operations
void model_config_defaults(ModelConfig *conf);
void read_model_config_file(ModelConfig *conf, const char *filename);
void model_apply_rate_mode(ModelConfig *conf);
const char *rate_mode_name(RateMode m);
void validate_model_config(const ModelConfig *conf);
void print_model_config(const ModelConfig *conf);
void print_derived_banner(const ModelConfig *conf);
void fprint_derived_block(FILE *f, const ModelConfig *conf);
void print_conf_json(const ModelConfig *conf, const char *dirname);
void create_output_directory(const ModelConfig *conf, char *dirname_out, size_t len);

#endif
