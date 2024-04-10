#include <stat.h>
#include <math.h>

double _stat_cov(const double *x, const double ex, const double *y, const double ey, const size_t n) {
  double cov = 0;

  for(size_t i = 0; i < n; i++)
    cov += x[i]*y[i];

  return cov/n-ex*ey;
}

double stat_mean(const double *x, const size_t n) {
  double mean = 0;

  for(size_t i = 0; i < n; i++)
    mean += x[i];

  return mean / n;
}

double stat_mean2(double *var, const double *x, const size_t n) {
  const double ex = stat_mean(x,n);
  if(var) *var = _stat_cov(x,ex,x,ex,n);
  return ex;
}

double stat_std(const double *x, const size_t n) {
  return sqrt(stat_var(x,n));
}

double stat_var(const double *x, const size_t n) {
  const double ex = stat_mean(x,n);
  return _stat_cov(x,ex,x,ex,n);
}

double stat_cov(const double *x, const double *y, const size_t n) {
  const double ex = stat_mean(x,n);
  const double ey = stat_mean(y,n);
  return _stat_cov(x,ex,y,ey,n);
}

double stat_corr(const double *x, const double *y, const size_t n) {
  double covm[3];
  stat_covm(covm, x, y, n);
  return covm[2]/sqrt(covm[0]*covm[1]);
}

double stat_abserr(const double *x, const double *y, const size_t n) {
  double abserr = 0;

  for(size_t i = 0; i < n; i++)
    abserr += fabs(x[i]-y[i]);

  return abserr / n;
}

void stat_covm(double *covm, const double *x, const double *y, const size_t n) {
  const double ex = stat_mean(x,n);
  const double ey = stat_mean(y,n);

  covm[0] = _stat_cov(x,ex,x,ex,n);
  covm[1] = _stat_cov(y,ey,y,ey,n);
  covm[2] = _stat_cov(x,ex,y,ey,n);
}
