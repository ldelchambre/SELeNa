/**
 * @file nelder_test.c
 * 
 * Nelder-Mead optimization algorithm, test file.
 * 
 * Copyright 2018 Ludovic Delchambre <ldelchambre@uliege.be>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */
#include <math.h>
#include <nelder.h>
#include <stdio.h>
#include <stdlib.h>

// Perform 65536 tests
#define NTEST (1U << 16)

// The alpha parameter to use
#define NELDER_ALPHA NELDER_DEFAULT_ALPHA

// The beta parameter to use
#define NELDER_BETA  NELDER_DEFAULT_BETA

// The gamma parameter to use
#define NELDER_GAMMA NELDER_DEFAULT_GAMMA

// The delta parameter to use
#define NELDER_DELTA NELDER_DEFAULT_DELTA

// The eps parameter to use
#define NELDER_EPS   NELDER_DEFAULT_EPS

#ifdef __GNUC__
#define UNUSED(param) __attribute__((unused)) param
#else
#define UNUSED(param) param
#endif

/*******************************************************************************
 * Typedefs
 ******************************************************************************/
typedef struct point_range {
  double min;
  double max;
} point_range_t;

typedef struct nelder_test_data {
  double (*fct)(const double *, void *); // The function to minimize
  int (*fct_check)(const double *, const double, void *); // The function designed to check the result of the minimization
  int (*fct_end)(const double *, const double *, void *); // The ending function to use or NULL in order to use the default one

  unsigned int n; // The problem size
  unsigned int narg; // The number of argument to the function

  point_range_t *x0lim; // The limits of the values of each point in the initial simplex
  point_range_t *args; // The limits of the arguments to pass to the function

} nelder_test_data_t;

/*******************************************************************************
 * Rosenbrock's parabolic valley
 ******************************************************************************/
double fct1(const double *x, void *arg) {
  const double *a = (const double *) arg;

  // Rosenbrock's parabolic valley (Rosenbrock, 1960)
  // Global minimum stand at (a[0], a[0]^2), where f(x) = 0
  // No local minimum
  const double a0 = (x[1] - x[0] * x[0]);
  const double a1 = (a[0] - x[0]);

  return a[1] * a0 * a0 + a1 * a1;
}

int fct1_check(const double *x, const double y, void *arg) {
  const double xtolerance = 1e-6;
  const double ytolerance = 1e-10;
  const double *a = (const double *) arg;

  if(fabs(a[0]-x[0]) > xtolerance || fabs(a[0]*a[0]-x[1]) > xtolerance
     || y > ytolerance) {
    fprintf(stderr, "Rosenbrock's parabolic valley\n");
    fprintf(stderr, "\tArguments: (%g, %g)\n", a[0], a[1]);
    fprintf(stderr, "\tFound optimum: f(%g, %g) = %g\n", x[0], x[1], y);
    fprintf(stderr, "\tEffective optimum: f(%g, %g) = 0\n", a[0], a[0]*a[0]);
    return -1;
  }

  return 0;
}

point_range_t fct1_x0lim[] = {{-1, 1}, {-1, 1}};
point_range_t fct1_args[]  = {{-1, 1}, {1, 100}};

/*******************************************************************************
 * Powell's quadratic function
 ******************************************************************************/
double fct2(const double *x, void *arg) {
  const double *a = (const double *) arg;

  // Powell's quadratic function (Powell, 1962)
  // Global optimum stand at (0, 0, 0, 0), where f(x) = 0
  // No local minimum
  const double a0 = x[0] + a[0] * x[1];
  const double a1 = x[1] + a[1] * x[2];
  const double a2 = x[2] + a[2] * x[3];
  const double a3 = x[3] + a[3] * x[0];
  return a[4] * a0 * a0 + a[5] * a1 * a1 + a[6] * a2 * a2 * a2 * a2 + a[7] * a3 * a3 * a3 * a3;
}

int fct2_check(const double *x, const double y, void *arg) {
  const double xtolerance = 1e-3;
  const double ytolerance = 1e-15;
  const double *a = (const double *) arg;

  if(   fabs(x[0]) > xtolerance || fabs(x[1]) > xtolerance
     || fabs(x[2]) > xtolerance || fabs(x[3]) > xtolerance
     || y > ytolerance) {
    fprintf(stderr, "Powell's quadratic function\n");
    fprintf(stderr, "\tArguments: (%g, %g, %g, %g, %g, %g, %g, %g)\n",
                    a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]);
    fprintf(stderr, "\tFound optimum: f(%g, %g, %g, %g) = %g\n", x[0], x[1], x[2], x[3], y);
    fprintf(stderr, "\tEffective optimum: f(0, 0, 0, 0) = 0\n");
    return -1;
  }

  return 0;
}

int fct2_end(const double UNUSED(*X), const double *y, void *arg) {
  const double ytolerance = 1e-15;
  const nelder_test_data_t *test = (nelder_test_data_t *) arg;
  const unsigned int n = test->n;
  unsigned int i;

  for(i = 0; i <= n; i++) {
    if(y[i] < ytolerance)
      return 1;
  }

  return 0;
}

point_range_t fct2_x0lim[] = {
  {-1,1}, {-1,1}, {-1,1}, {-1,1}
};
point_range_t fct2_args[]  = {
  {-1,1}, {-1,1}, {-1,1}, {-1,1}, {1,2}, {1,2}, {1,2}, {1,2}
};

/*******************************************************************************
 * Test cases execution
 ******************************************************************************/
int perform_test(nelder_test_data_t *test, unsigned int seed) {
  double X[test->n * ( test->n + 1 )];
  double  y;
  double arg[test->narg];
  unsigned long i, j;

  // Initialize the random number generator
  srand(seed);

  // Pick up random arguments
  for(i = 0; i < test->narg; i++)
    arg[i] = (test->args[i].max - test->args[i].min) * rand() / (1. + RAND_MAX) + test->args[i].min;

  {
  double (*_X)[test->n] = (double (*)[test->n]) X;

  // Pick up random points
  for(i = 0; i <= test->n; i++)
    for(j = 0; j < test->n; j++)
      _X[i][j] = (test->x0lim[j].max - test->x0lim[j].min) * rand() / RAND_MAX + test->x0lim[j].min;
  }

  // Execute the Nelder-Meader algorithm
  if(test->fct_end)
    y = nelder2(X, test->n, test->fct, arg, NELDER_ALPHA, NELDER_BETA, NELDER_GAMMA, NELDER_DELTA, test->fct_end, test);
  else
    y = nelder(X, test->n, test->fct, arg, NELDER_ALPHA, NELDER_BETA, NELDER_GAMMA, NELDER_DELTA, NELDER_EPS);

  // Check the result of the optimization
  return test->fct_check(X, y, arg);
}

nelder_test_data_t tests[] = {
  {fct1, fct1_check, NULL,     2, 2, fct1_x0lim, fct1_args},
  {fct2, fct2_check, fct2_end, 4, 8, fct2_x0lim, fct2_args},
  {NULL, NULL, NULL, 0, 0, NULL, NULL}
};

int main(void) {
  unsigned int i, j;

  for(i = 0; i < NTEST; i++) {
    for(j = 0; tests[j].fct; j++) {
      if(perform_test(&tests[j], 1 + i)) {
        fprintf(stderr, "\tSeed: %u\n", 1+i);
      }
    }
  }

  exit(EXIT_SUCCESS);
}
