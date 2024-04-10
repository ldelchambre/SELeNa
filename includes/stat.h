#ifndef _STAT_H_
#define _STAT_H_

#include <stddef.h>

double stat_mean(const double *x, const size_t n);
double stat_mean2(double *var, const double *x, const size_t n);
double stat_std(const double *x, const size_t n);
double stat_var(const double *x, const size_t n);
double stat_cov(const double *x, const double *y, const size_t n);
double stat_corr(const double *x, const double *y, const size_t n);
double stat_abserr(const double *x, const double *y, const size_t n);
void stat_covm(double *covm, const double *x, const double *y, const size_t n);

#endif
