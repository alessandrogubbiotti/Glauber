#ifndef CONFIGURATION_STATISTICS_H
#define CONFIGURATION_STATISTICS_H

#include <stdio.h>
#include <stdlib.h>

// Fractional macroscopic coordinates — multiplied by N to get microscopic indices.
typedef struct {
    float A1_magnetization;
    float B1_magnetization;
    float A1_corr;
    float A2_corr;
    float B1_corr;
    float B2_corr;
    float A1_time_delayed_corr;
    float A2_time_delayed_corr;
    float B1_time_delayed_corr;
    float B2_time_delayed_corr;
} StatsConfig;

// Microscopic integer intervals derived from StatsConfig.
typedef struct {
    int a1;
    int b1;
} Interval;

typedef struct {
    int a1;
    int a2;
    int b1;
    int b2;
} Two_Intervals;

typedef int (*StatsConfigSetter)(StatsConfig *conf, const char *value);

typedef struct {
    const char *key;
    StatsConfigSetter setter;
} StatsConfigEntry;

// Setters
int set_A1_magnetization(StatsConfig *conf, const char *value);
int set_B1_magnetization(StatsConfig *conf, const char *value);
int set_A1_corr(StatsConfig *conf, const char *value);
int set_A2_corr(StatsConfig *conf, const char *value);
int set_B1_corr(StatsConfig *conf, const char *value);
int set_B2_corr(StatsConfig *conf, const char *value);
int set_A1_time_delayed_corr(StatsConfig *conf, const char *value);
int set_A2_time_delayed_corr(StatsConfig *conf, const char *value);
int set_B1_time_delayed_corr(StatsConfig *conf, const char *value);
int set_B2_time_delayed_corr(StatsConfig *conf, const char *value);

// Config operations
void read_stats_config_file(StatsConfig *conf, const char *filename);
void validate_stats_config(const StatsConfig *conf, int N, int L);

// Translate fractional StatsConfig coordinates to integer microscopic intervals.
void set_interval_magnetization(const StatsConfig *conf, int N, Interval *interval);
void set_correlation_params(const StatsConfig *conf, int N, Two_Intervals *intervals);
void set_time_delayed_params(const StatsConfig *conf, int N, Two_Intervals *intervals);

#endif
