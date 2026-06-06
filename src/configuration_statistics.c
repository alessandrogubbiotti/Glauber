#include "configuration_statistics.h"
#include <string.h>
#include <math.h>
#include <ctype.h>

#define MAX_LINE 256

static int is_comment_or_blank(const char *line) {
    while (isspace(*line)) line++;
    return (*line == '#' || *line == '/' || *line == '\0' || *line == '\n');
}

int set_A1_magnetization(StatsConfig *conf, const char *value) { conf->A1_magnetization = atof(value); return 0; }
int set_B1_magnetization(StatsConfig *conf, const char *value) { conf->B1_magnetization = atof(value); return 0; }
int set_A1_corr(StatsConfig *conf, const char *value)          { conf->A1_corr = atof(value);          return 0; }
int set_A2_corr(StatsConfig *conf, const char *value)          { conf->A2_corr = atof(value);          return 0; }
int set_B1_corr(StatsConfig *conf, const char *value)          { conf->B1_corr = atof(value);          return 0; }
int set_B2_corr(StatsConfig *conf, const char *value)          { conf->B2_corr = atof(value);          return 0; }

int set_A1_time_delayed_corr(StatsConfig *conf, const char *value) { conf->A1_time_delayed_corr = atof(value); return 0; }
int set_A2_time_delayed_corr(StatsConfig *conf, const char *value) { conf->A2_time_delayed_corr = atof(value); return 0; }
int set_B1_time_delayed_corr(StatsConfig *conf, const char *value) { conf->B1_time_delayed_corr = atof(value); return 0; }
int set_B2_time_delayed_corr(StatsConfig *conf, const char *value) { conf->B2_time_delayed_corr = atof(value); return 0; }

static StatsConfigEntry stats_config_table[] = {
    {"A1_magnetization",      set_A1_magnetization},
    {"B1_magnetization",      set_B1_magnetization},
    {"A1_corr",               set_A1_corr},
    {"A2_corr",               set_A2_corr},
    {"B1_corr",               set_B1_corr},
    {"B2_corr",               set_B2_corr},
    {"A1_time_delayed_corr",  set_A1_time_delayed_corr},
    {"A2_time_delayed_corr",  set_A2_time_delayed_corr},
    {"B1_time_delayed_corr",  set_B1_time_delayed_corr},
    {"B2_time_delayed_corr",  set_B2_time_delayed_corr},
    {NULL, NULL}
};

void read_stats_config_file(StatsConfig *conf, const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening config file");
        exit(EXIT_FAILURE);
    }

    char line[MAX_LINE], key[64], value[128];
    while (fgets(line, sizeof(line), file)) {
        if (is_comment_or_blank(line)) continue;
        if (sscanf(line, "%63[^=]=%127s", key, value) != 2) continue;

        for (int i = 0; stats_config_table[i].key; i++) {
            if (strcmp(stats_config_table[i].key, key) == 0) {
                if (stats_config_table[i].setter(conf, value) != 0) {
                    fprintf(stderr, "Invalid value for key: %s\n", key);
                    exit(EXIT_FAILURE);
                }
                break;
            }
        }
        // Keys not in this table belong to ModelConfig — silently skipped.
    }

    fclose(file);
}

static int is_scaled_integer(float val, int N) {
    return fabsf(val * N - roundf(val * N)) < 1e-5;
}

void validate_stats_config(const StatsConfig *conf, int N, int L) {
    float coords[] = {
        conf->A1_magnetization, conf->B1_magnetization,
        conf->A1_corr,          conf->A2_corr,
        conf->B1_corr,          conf->B2_corr,
        conf->A1_time_delayed_corr, conf->A2_time_delayed_corr,
        conf->B1_time_delayed_corr, conf->B2_time_delayed_corr
    };
    for (int i = 0; i < (int)(sizeof(coords) / sizeof(float)); i++) {
        if (coords[i] <= 0 || coords[i] >= L || !is_scaled_integer(coords[i], N)) {
            fprintf(stderr, "Statistic parameter %.2f is invalid: must be in (0, L) and N*value must be an integer\n", coords[i]);
            exit(EXIT_FAILURE);
        }
    }
}

void set_interval_magnetization(const StatsConfig *conf, int N, Interval *interval) {
    interval->a1 = (int)(N * conf->A1_magnetization);
    interval->b1 = (int)(N * conf->B1_magnetization);
}

void set_correlation_params(const StatsConfig *conf, int N, Two_Intervals *intervals) {
    intervals->a1 = (int)(N * conf->A1_corr);
    intervals->a2 = (int)(N * conf->A2_corr);
    intervals->b1 = (int)(N * conf->B1_corr);
    intervals->b2 = (int)(N * conf->B2_corr);
}

void set_time_delayed_params(const StatsConfig *conf, int N, Two_Intervals *intervals) {
    intervals->a1 = (int)(N * conf->A1_time_delayed_corr);
    intervals->a2 = (int)(N * conf->A2_time_delayed_corr);
    intervals->b1 = (int)(N * conf->B1_time_delayed_corr);
    intervals->b2 = (int)(N * conf->B2_time_delayed_corr);
}
