#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <stdlib.h>

#define MAX_LINE 256

typedef void (*Initialization_function)(int *spin, int l, void *variables);

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

// Config operations
void read_model_config_file(ModelConfig *conf, const char *filename);
void validate_model_config(const ModelConfig *conf);
void print_model_config(const ModelConfig *conf);
void print_conf_json(const ModelConfig *conf, const char *dirname);
void create_output_directory(const ModelConfig *conf, char *dirname_out, size_t len);

#endif
