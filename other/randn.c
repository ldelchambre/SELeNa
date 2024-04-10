#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define NDRAW (1U << 24)

double randn(void) {
  static char iset = 0;
  static double next = 0;
  double x1, x2, r2, v;
  if(iset == 0) {
    do {
      x1 = 2. / RAND_MAX * rand() - 1;
      x2 = 2. / RAND_MAX * rand() - 1;
      r2 = x1 * x1 + x2 * x2;
    } while(r2 >= 1.0 || r2 == 0.0);
    v = sqrt( -2. * log(r2) / r2);
    next = x2 * v;
    iset = 1;
    return x1 * v;
  } else {
    iset = 0;
    return next;
  }
}

int main(int argc, char **argv) {
  unsigned int i, nn = 0, n1s = 0, n2s = 0, n3s = 0;
  double r, fr;

  for(i = 0; i < NDRAW; i++) {
    r = randn();
    fr = fabs(r);

    nn += (r < 0.);
    n1s += (fr < 1.);
    n2s += (fr < 2.);
    n3s += (fr < 3.);
  }

  fprintf(stdout, "Ratio of negative values: %g%%\n", 100.*nn/NDRAW);
  fprintf(stdout, "Ratio of draw between -1 and 1: %g%% (theory: %g%%)\n", 100.*n1s/NDRAW, 100.*erf(M_SQRT1_2));
  fprintf(stdout, "Ratio of draw between -2 and 2: %g%% (theory: %g%%)\n", 100.*n2s/NDRAW, 100.*erf(2.*M_SQRT1_2));
  fprintf(stdout, "Ratio of draw between -3 and 3: %g%% (theory: %g%%)\n", 100.*n3s/NDRAW, 100.*erf(3.*M_SQRT1_2));

  return 0;
}
