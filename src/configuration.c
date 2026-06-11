#include "configuration.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>

static int is_comment_or_blank(const char *line) {
    while (isspace(*line)) line++;
    return (*line == '#' || *line == '/' || *line == '\0' || *line == '\n');
}

int set_N(ModelConfig *conf, const char *value) {
    int v = atoi(value);
    if (v <= 0) return -1;
    conf->N = v;
    return 0;
}

int set_T(ModelConfig *conf, const char *value) {
    double v = atof(value);
    if (v <= 0.0) return -1;
    conf->T = v;
    return 0;
}

int set_L(ModelConfig *conf, const char *value) {
    int v = atoi(value);
    if (v <= 0) return -1;
    conf->L = v;
    return 0;
}

int set_annihilation(ModelConfig *conf, const char *value) {
    conf->annihilation = atof(value);
    return 0;
}

int set_creation(ModelConfig *conf, const char *value) {
    conf->creation = atof(value);
    return 0;
}

int set_N_simulations(ModelConfig *conf, const char *value) {
    int v = atoi(value);
    if (v <= 0) return -1;
    conf->N_simulations = v;
    return 0;
}

int set_resolution(ModelConfig *conf, const char *value) {
    int v = atoi(value);
    if (v <= 0) return -1;
    conf->resolution = v;
    return 0;
}

int set_Micro_n_steps(ModelConfig *conf, const char *value) {
    int v = atoi(value);
    if (v <= 0) return -1;
    conf->Micro_n_steps = v;
    return 0;
}

// Maps the integer key in the config file to the corresponding init function pointer.
int set_init_function(ModelConfig *conf, const char *value) {
    int v = atoi(value);
    switch (v) {
        case 1: conf->initialize = initialize_spins_1; break;
        case 2: conf->initialize = initialize_spins_2; break;
        case 3: conf->initialize = initialize_spins_3; break;
        case 4: conf->initialize = initialize_spins_4; break;
        case 5: conf->initialize = initialize_spins_5; break;
        case 6: conf->initialize = initialize_spins_6; break;
        default: return -1;
    }
    return 0;
}

// For init functions that need a scalar parameter (e.g. lambda for geometric, p for Bernoulli).
int set_init_param(ModelConfig *conf, const char *value) {
    double *p = malloc(sizeof(double));
    if (!p) return -1;
    *p = atof(value);
    conf->initialization_param = p;
    return 0;
}

static ModelConfigEntry model_config_table[] = {
    {"N",              set_N},
    {"T",              set_T},
    {"L",              set_L},
    {"annihilation",   set_annihilation},
    {"creation",       set_creation},
    {"N_simulations",  set_N_simulations},
    {"resolution",     set_resolution},
    {"Micro_n_steps",  set_Micro_n_steps},
    {"init_function",  set_init_function},
    {"init_param",     set_init_param},
    {NULL, NULL}
};

void read_model_config_file(ModelConfig *conf, const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening config file");
        exit(EXIT_FAILURE);
    }

    char line[MAX_LINE], key[64], value[128];
    while (fgets(line, sizeof(line), file)) {
        if (is_comment_or_blank(line)) continue;
        if (sscanf(line, "%63[^=]=%127s", key, value) != 2) continue;

        for (int i = 0; model_config_table[i].key; i++) {
            if (strcmp(model_config_table[i].key, key) == 0) {
                if (model_config_table[i].setter(conf, value) != 0) {
                    fprintf(stderr, "Invalid value for key: %s\n", key);
                    exit(EXIT_FAILURE);
                }
                break;
            }
        }
        // Keys not in this table belong to StatsConfig — silently skipped.
    }

    fclose(file);
}

void validate_model_config(const ModelConfig *conf) {
    if (conf->N <= 0 || conf->T <= 0.0 || conf->L <= 0) {
        fprintf(stderr, "N, T, L must be positive\n");
        exit(EXIT_FAILURE);
    }
    if ((conf->N * conf->N * conf->Micro_n_steps) % conf->resolution != 0) {
        fprintf(stderr, "N*N*Micro_n_steps must be divisible by resolution\n");
        exit(EXIT_FAILURE);
    }
    if (conf->initialize == NULL) {
        fprintf(stderr, "init_function not set in config file\n");
        exit(EXIT_FAILURE);
    }
}

void print_model_config(const ModelConfig *conf) {
    printf("Model configuration:\n");
    printf("  N:              %d\n",   conf->N);
    printf("  T:              %lf\n",  conf->T);
    printf("  L:              %d\n",   conf->L);
    printf("  annihilation:   %lf\n",  conf->annihilation);
    printf("  creation:       %lf\n",  conf->creation);
    printf("  N_simulations:  %d\n",   conf->N_simulations);
    printf("  resolution:     %d\n",   conf->resolution);
    printf("  Micro_n_steps:  %d\n",   conf->Micro_n_steps);
}

void print_conf_json(const ModelConfig *conf, const char *dirname) {
    char filepath[512];
    snprintf(filepath, sizeof(filepath), "%s/configuration.json", dirname);
    FILE *f = fopen(filepath, "w");
    if (!f) { perror("Failed to write configuration.json"); exit(EXIT_FAILURE); }
    fprintf(f, "{\n");
    fprintf(f, "  \"N\": %d,\n",              conf->N);
    fprintf(f, "  \"T\": %lf,\n",             conf->T);
    fprintf(f, "  \"L\": %d,\n",              conf->L);
    fprintf(f, "  \"annihilation\": %.5f,\n", conf->annihilation);
    fprintf(f, "  \"creation\": %.5f,\n",     conf->creation);
    fprintf(f, "  \"N_simulations\": %d,\n",  conf->N_simulations);
    fprintf(f, "  \"resolution\": %d,\n",     conf->resolution);
    fprintf(f, "  \"Micro_n_steps\": %d\n",   conf->Micro_n_steps);
    fprintf(f, "}\n");
    fclose(f);
}

void create_output_directory(const ModelConfig *conf, char *dirname_out, size_t len) {
    snprintf(dirname_out, len, "results_N%d_T%.1lf_L%d_Ann%.2f_micro%d",
             conf->N, conf->T, conf->L, conf->annihilation, conf->Micro_n_steps);
    if (mkdir(dirname_out, 0755) && errno != EEXIST) {
        perror("mkdir failed");
        exit(EXIT_FAILURE);
    }
}
