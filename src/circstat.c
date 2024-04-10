#include <circstat.h>

#include <float.h>
#include <math.h>

double circstat_mean(const double *x, const size_t n) {
  double ssin, scos;
  size_t i;

  ssin = scos = 0;
  for(i = 0; i < n; i++) {
    ssin += sin(x[i]);
    scos += cos(x[i]);
  }

  return atan2(ssin,scos);
}

double circstat_r(const double *x, const size_t n) {
  double ssin, scos;
  size_t i;

  ssin = scos = 0;
  for(i = 0; i < n; i++) {
    ssin += sin(x[i]);
    scos += cos(x[i]);
  }

  return sqrt(scos*scos+ssin*ssin)/n;
}

double circstat_var(const double *x, const size_t n) {
  return 1.-circstat_r(x, n);
}

double circstat_diff(const double a, const double b) {
  const double d = a - b;
  return atan2(sin(d),cos(d));
}

double circstat_corr(const double *x, const double *y, const size_t n) {
  double xmean, ymean;
  double dx, dy;
  double xvar, yvar, covar;
  size_t i;

  // Compute the mean angle of x and y
  xmean = circstat_mean(x,n);
  ymean = circstat_mean(y,n);

  // Compute the covariance matrix
  xvar = yvar = covar = 0;
  for(i = 0; i < n; i++) {
    dx = sin(x[i]-xmean);
    dy = sin(y[i]-ymean);
    covar += dx*dy;
    xvar += dx*dx;
    yvar += dy*dy;
  }

  return covar/sqrt(xvar*yvar);
}

double circstat_abserr(const double *x, const double *y, const size_t n) {
  double adiff, abserr = 0;
  size_t i;

  for(i = 0; i < n; i++) {
    adiff = x[i]-y[i];
    abserr += fabs(atan2(sin(adiff),cos(adiff)));
  }

  return abserr/n;
}

double circstat_confmean(const double *x, const size_t n, const double sigma) {
  const double r = circstat_r(x, n);
  const double R = r*n;
  const double c2 = sigma*sigma;
  double t;

  if(r < 0.9 && r > sqrt(c2/2/n))
    t = sqrt((2.*n*(2.*R*R-n*c2))/(4.*n-c2));
  else if(r >= 0.9) {
    const double n2 = (double) n*n;
    t = sqrt(n2-(n2-R*R)*exp(c2/n));
  } else
    return -1;
  return acos(t/R);
}

double circstat_kappa(const double *x, const size_t n) {
  const double r = circstat_r(x, n);
  const double r2 = r*r;
  double k;

  // XXX A MLE of kappa can be acquired by finding I1(kappa)/I0(kappa)-r=0; where I0, I1 are the Bessel functions of order 0 and 1, respectively.

  if(r < 0.53)
    k = r*(2.+r2*(1.+5./6*r2));
  else if(r < 0.85)
    k = -0.4+1.39*r+0.43/(1.-r);
  else
    k = 1./(3.-4.*r+r2)/r;

  if(n <= 15) {
    if(k < 2) {
      k -= 2./(k*n);
      k = (k < 0) ? 0 : k;
    } else {
      const double _n = 1.-n;
      k *= _n*_n*_n/(1.+n*n)/n;
    }
  }

  return k;
}
