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
    conf->ann_set = 1;
    return 0;
}

int set_creation(ModelConfig *conf, const char *value) {
    conf->creation = atof(value);
    conf->create_set = 1;
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

int set_nu(ModelConfig *conf, const char *value) {
    double v = atof(value);
    if (v <= 0.0) return -1;
    conf->nu = v;
    return 0;
}

/* Accepts the strings glauber/scaled/explicit (case-insensitive) or 0/1/2. */
int set_rate_mode(ModelConfig *conf, const char *value) {
    char buf[32];
    size_t i = 0;
    while (value[i] && i < sizeof(buf) - 1) { buf[i] = (char)tolower((unsigned char)value[i]); i++; }
    buf[i] = '\0';
    if      (!strcmp(buf, "glauber")  || !strcmp(buf, "0")) conf->rate_mode = RATE_GLAUBER;
    else if (!strcmp(buf, "scaled")   || !strcmp(buf, "1")) conf->rate_mode = RATE_SCALED;
    else if (!strcmp(buf, "explicit") || !strcmp(buf, "2")) conf->rate_mode = RATE_EXPLICIT;
    else return -1;
    return 0;
}

int set_slowdown_exponent(ModelConfig *conf, const char *value) {
    conf->slowdown_exponent = atof(value);
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
    {"N",                 set_N},
    {"T",                 set_T},
    {"L",                 set_L},
    {"annihilation",      set_annihilation},
    {"creation",          set_creation},
    {"N_simulations",     set_N_simulations},
    {"resolution",        set_resolution},
    {"Micro_n_steps",     set_Micro_n_steps},
    {"init_function",     set_init_function},
    {"init_param",        set_init_param},
    {"nu",                set_nu},
    {"rate_mode",         set_rate_mode},
    {"slowdown_exponent", set_slowdown_exponent},
    {NULL, NULL}
};

const char *rate_mode_name(RateMode m) {
    switch (m) {
        case RATE_GLAUBER:  return "glauber";
        case RATE_SCALED:   return "scaled";
        case RATE_EXPLICIT: return "explicit";
        default:            return "unknown";
    }
}

/* Sensible defaults, independent of how the struct was allocated. */
void model_config_defaults(ModelConfig *conf) {
    memset(conf, 0, sizeof(*conf));
    conf->nu                = 1.0;
    conf->rate_mode         = RATE_GLAUBER;
    conf->slowdown_exponent = 0.0;
}

/* Derive gamma/beta/eta/p and the (ann, create) rates from (N, nu, rate_mode),
 * following docs/section1.tex.  Must be called after the config file is read. */
void model_apply_rate_mode(ModelConfig *conf) {
    if (conf->N <= 0 || conf->nu <= 0.0) {
        /* validate_model_config / the per-binary checks will report this. */
        return;
    }

    double N  = (double)conf->N;
    double nu = conf->nu;

    conf->time_scale_formula = 2.0 * N * N;

    if (conf->rate_mode == RATE_EXPLICIT) {
        /* Rates come straight from the file; report what they imply. */
        double ann = conf->annihilation, create = conf->creation;
        if (!(ann > 0.0) || !(create > 0.0)) {
            fprintf(stderr,
                "WARNING: rate_mode=explicit requires positive annihilation and "
                "creation keys (got ann=%g, create=%g).\n", ann, create);
            return;
        }
        double nu_implied    = N * sqrt(create / ann);
        double gamma_implied = (ann - create) / (ann + create);
        conf->nu            = nu_implied;
        conf->gamma         = gamma_implied;
        conf->beta          = 0.5 * log(N / nu_implied);
        conf->eta           = (N - nu_implied) / (N + nu_implied);
        conf->p_equilibrium = 1.0 / (1.0 + sqrt(ann / create));
        if (fabs(ann + create - 2.0) > 1e-9)
            fprintf(stderr,
                "WARNING: ann + create = %.10g != 2; the correlation hierarchy "
                "does not close, so the analytic formulas are only asymptotic.\n",
                ann + create);
        return;
    }

    /* glauber / scaled: derive everything from (N, nu). */
    if ((conf->ann_set || conf->create_set))
        fprintf(stderr,
            "WARNING: rate_mode=%s derives the rates from (N, nu); the "
            "annihilation/creation keys in the config are ignored.\n",
            rate_mode_name(conf->rate_mode));

    double gamma = (N * N - nu * nu) / (N * N + nu * nu);
    conf->gamma         = gamma;
    conf->beta          = 0.5 * log(N / nu);
    conf->eta           = (N - nu) / (N + nu);
    conf->p_equilibrium = nu / (N + nu);

    double ann    = 1.0 + gamma;
    double create = 1.0 - gamma;
    if (conf->rate_mode == RATE_SCALED) {
        double s = pow(N, conf->slowdown_exponent);
        ann    *= s;
        create *= s;
    }
    conf->annihilation = ann;
    conf->creation     = create;
}

void read_model_config_file(ModelConfig *conf, const char *filename) {
    model_config_defaults(conf);

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

    /* Fix (gamma, beta, eta, p) and (ann, create) from (N, nu, rate_mode). */
    model_apply_rate_mode(conf);
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
    print_derived_banner(conf);
}

/* Human-readable derived block for the run banner. */
void print_derived_banner(const ModelConfig *conf) {
    printf("  --- derived (section1.tex) ---\n");
    printf("  rate_mode:      %s", rate_mode_name(conf->rate_mode));
    if (conf->rate_mode == RATE_SCALED)
        printf("  (a=%.4g)", conf->slowdown_exponent);
    printf("\n");
    printf("  nu:             %.10g\n", conf->nu);
    printf("  beta:           %.10g\n", conf->beta);
    printf("  gamma:          %.10g\n", conf->gamma);
    printf("  eta:            %.10g\n", conf->eta);
    printf("  p_equilibrium:  %.10g\n", conf->p_equilibrium);
    printf("  ann:            %.10g\n", conf->annihilation);
    printf("  create:         %.10g\n", conf->creation);
    printf("  time_scale:     2*N^2 = %.10g\n", conf->time_scale_formula);
}

/* Common derived JSON fields, shared by every binary's configuration.json. */
void fprint_derived_block(FILE *f, const ModelConfig *conf) {
    fprintf(f, "  \"nu\": %.10g,\n",                  conf->nu);
    fprintf(f, "  \"beta\": %.10g,\n",                conf->beta);
    fprintf(f, "  \"gamma\": %.10g,\n",               conf->gamma);
    fprintf(f, "  \"eta\": %.10g,\n",                 conf->eta);
    fprintf(f, "  \"p_equilibrium\": %.10g,\n",       conf->p_equilibrium);
    fprintf(f, "  \"rate_mode\": \"%s\",\n",          rate_mode_name(conf->rate_mode));
    fprintf(f, "  \"slowdown_exponent\": %.10g,\n",   conf->slowdown_exponent);
    fprintf(f, "  \"ann\": %.10g,\n",                 conf->annihilation);
    fprintf(f, "  \"create\": %.10g,\n",              conf->creation);
    fprintf(f, "  \"time_scale_formula\": %.10g\n",   conf->time_scale_formula);
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
    fprintf(f, "  \"annihilation\": %.10g,\n", conf->annihilation);
    fprintf(f, "  \"creation\": %.10g,\n",     conf->creation);
    fprintf(f, "  \"N_simulations\": %d,\n",  conf->N_simulations);
    fprintf(f, "  \"resolution\": %d,\n",     conf->resolution);
    fprintf(f, "  \"Micro_n_steps\": %d,\n",  conf->Micro_n_steps);
    fprint_derived_block(f, conf);
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
