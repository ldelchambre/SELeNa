#include <graddlsqr.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <utils.h>

// Problem size
#define NATTR 8     // The number of attributes in the problem
#define NINST 65536 // The number of instances
#define RAND_SEED 0 // The random seed to use

// Gradient descent parameters
#define GRADD_EPS        1e-15      // The required precision in the change of the cost function
#define GRADD_BATCH_SIZE 0.5*NINST  // The size of the chunck to treat
#define GRADD_NTHREAD    4          // The number of threads to run
#define GRADD_MAXITER    10000      // The maximal number of iterations to perform
#define GRADD_VERBOSE    1          // Should we run in verbose mode

// y = a[0] * x[0] + ... + a[nattr-1] * x[nattr-1] + a[nattr]
double f(const double *a, const double *x, const unsigned long nattr, void *arg) {
  const unsigned long na = *((unsigned long *) arg);
  double tmp;
  unsigned long i;

  tmp = a[na - 1];
  for (i = 0; i < nattr; i++)
    tmp += x[i] * a[i];

  return tmp;
}

// dy/da = ( x[0] ... x[nattr-1] 1 )
void df(double *grad, const double *a, const double *x,
    const unsigned long nattr, void *arg) {
  const unsigned long na = *((unsigned long *) arg);
  unsigned long i;

  grad[na - 1] = 1.;
  for (i = 0; i < nattr; i++)
    grad[i] = x[i];
}

int main(void) {
  // The number of parameters we will fit
  const unsigned long na = NATTR + 1;
  // Allocate vector and matrices
  double *a = (double *) calloc(na, sizeof(double));
  double *a0 = (double *) calloc(na, sizeof(double));
  double *areal = (double *) calloc(na, sizeof(double));
  double *X = (double *) calloc(NATTR*NINST, sizeof(double));
  double *y = (double *) calloc(NINST, sizeof(double));
  double *w = (double *) calloc(NINST, sizeof(double));
  // The random seed to use
  unsigned int seed = RAND_SEED;
  // The parameters of the gradient descent algorithm
  graddlsqr_opt_t opt = GRADDLSQR_DEFAULT_OPTS;
  // General usage counters
  unsigned long i, j;

  // Check that all allocations succeed
  if(a == NULL || a0 == NULL || areal == NULL || X == NULL || y == NULL || w == NULL) {
    fprintf(stderr, "Arrays allocation fails");
    free(a);free(a0);free(areal);free(X);free(y);free(w);
    return -1;
  }

  // Randomly initialize a and areal with numbers uniform drawn in [-1,1],
  // then copy the initial estimate, a, into a0
  for (i = 0; i < na; i++) {
    a[i] = 2. * mrand_r(&seed) / MRAND_MAX - 1;
    areal[i] = 2. * mrand_r(&seed) / MRAND_MAX - 1;
    a0[i] = a[i];
  }

  // Randomly initialize X with numbers uniform drawn in [-1,1]
  // Compute the associated outcome value, y
  // Assign weight, w, to each observation
  for (i = 0; i < NINST; i++) {
    double *x = X + i * NATTR;
    for (j = 0; j < NATTR; j++)
      x[j] = 2. * mrand_r(&seed) / MRAND_MAX - 1;
    y[i] = f(areal, x, NATTR, (void *) &na);
    w[i] = 1.;
  }

  // Initialize the options of the gradient descent algorithm
  opt.eps = GRADD_EPS;
  opt.batch_size = GRADD_BATCH_SIZE;
  opt.nthread = GRADD_NTHREAD;
  opt.maxiter = GRADD_MAXITER;
  opt.verbose = GRADD_VERBOSE;

  // Optimize the cost function with respect to a
  graddlsqr(a, na, X, y, w, NINST, NATTR, f, df, (void *) &na, &opt);

  // Print result
  fprintf(stdout, "   Real solution: (");
  for (i = 0; i < na; i++)
    fprintf(stdout, " %+g", areal[i]);
  fprintf(stdout, " )\n  Found solution: (");
  for (i = 0; i < na; i++)
    fprintf(stdout, " %+g", a[i]);
  fprintf(stdout, " )\nInitial solution: (");
  for (i = 0; i < na; i++)
    fprintf(stdout, " %+g", a0[i]);
  fprintf(stdout, " )\n");

  // Free vectors and matrices
  free(a);free(a0);free(areal);free(X);free(y);free(w);

  return 0;
}
