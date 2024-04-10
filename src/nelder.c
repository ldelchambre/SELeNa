/**
 * @file nelder.c
 * 
 * Nelder-Mead optimization algorithm, source file.
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
#include <nelder.h>

// #define _NELDER_DEBUG // Uncomment in order to print debug messages

#ifdef _NELDER_DEBUG
#include <stdio.h> // In debug mode, stdio.h is needed for the fprintf prototype
#endif

#ifdef __GNUC__
/**
 * Mark attributes as unused in order to prevent warning flags.
 */
#define UNUSED(param) __attribute__((unused)) param
#else
#define UNUSED(param) param
#endif

/**
 * The arguments that are needed by the default end function of the Nelder-Mead
 * algorithm.
 */
typedef struct _nelder_default_end_fct_arg {
  double eps2; ///< The square of the eps parameter that was passed to the nelder function
  unsigned int n; ///< The size of the input parameters
} _nelder_default_end_fct_arg_t;

/**
 * The default end function of the Nelder-Mead algorithm.
 * 
 *  This function return 1 if the variance of the values of the function 
 * evaluated in each of the vertices of the simplex falls below the limit given
 * by arg->eps2.
 * 
 * @param X   The last computed simplex (not used here)
 * @param y   The values of the function evaluated in each vertices of the simplex
 * @param arg Pointer to a structure of type _nelder_default_end_fct_arg_t, containing the arguments to pass to the function
 * 
 * @return 1 if the Nelder-Mead algorithm should stop iterating, 0 otherwise
 */
int _nelder_default_end_fct(const double UNUSED(*X), const double *y, void *arg) {
  const double eps2    = ((_nelder_default_end_fct_arg_t *) arg)->eps2;
  const unsigned int n = ((_nelder_default_end_fct_arg_t *) arg)->n;
  double mean, var;
  unsigned int i;

  // Compute the mean and variance of y
  mean = var = 0;
  for(i = 0; i <= n; i++) {
    mean += y[i];
    var  += y[i] * y[i];
  }
  mean /= (1. + n);
  var  = var / (1. + n) - mean * mean;

  return (var < eps2);
}

double nelder(double *X, const unsigned int n,
              double (*fct)(const double *, void *), void *fct_arg,
              const double alpha, const double beta, const double gamma, const double delta,
              const double eps) {
  _nelder_default_end_fct_arg_t arg = {eps*eps, n};
  return nelder2(X, n, fct, fct_arg, alpha, beta, gamma, delta, _nelder_default_end_fct, &arg);
}

double nelder2(double *_X, const unsigned int n,
            double (*fct)(const double *, void *), void *fct_arg,
            const double alpha, const double beta, const double gamma, const double delta,
            int (*end_fct)(const double *, const double *, void *), void *end_fct_arg) {
  // Re-interpret X for clarity purpose
  double (*X)[n] = (double (*)[n]) _X;

  // The Y values associated with each X[i], i in [0, n]
  double y[n+1];

  // The centroid of X
  double P[n];

  // Candidates points and associated values
  double ys, yss, Ps[n], Pss[n];

  // Indices of the minimal and maximal yvalue
  unsigned int il, ih;

  // General usage counters
  unsigned int i, j;

#ifdef _NELDER_DEBUG
  // In debug mode, the number of iteration that were performed
  unsigned int iter = 0;
#endif

  // Evaluate the function in each point
  for(i = 0; i <= n; i++)
    y[i] = fct(X[i], fct_arg);

#ifdef _NELDER_DEBUG
  fprintf(stderr, "Nelder-Mead initial conditions:\n");
  fprintf(stderr, "\talpha = %f, beta = %f, gamma = %f, delta = %f\n", alpha, beta, gamma, delta);
  for(i = 0; i <= n; i++) {
    fprintf(stderr, "\tx[%u] = (", i);
    for(j = 0; j < n; j++)
      fprintf(stderr, "%s%g", (j > 0) ? ", " : "", X[i][j]);
    fprintf(stderr, "), y[%u] = %g\n", i, y[i]);
  }
#endif

  // Stop the Nelder-Mead algorithm when the end function return something different than 0
  while(end_fct(_X,y,end_fct_arg) == 0) {
    // Get the indices of the minimal and maximal points
    il = ih = 0;
    for(i = 1; i <= n; i++) {
      if(y[i] < y[il]) il = i;
      if(y[i] > y[ih]) ih = i;
    }

    // Compute the centroid of X[i], i != ih
    for(j = 0; j < n; j++)
      P[j] = 0;
    for(i = 0; i <= n; i++) {
      if(i != ih) {
        for(j = 0; j < n; j++)
          P[j] += X[i][j];
      }
    }
    for(j = 0; j < n; j++)
      P[j] /= n;

#ifdef _NELDER_DEBUG
    fprintf(stderr, "%u) Centroid -> P = (", iter);
    for(j = 0; j < n; j++) fprintf(stderr, "%s%g", (j > 0) ? ", " : "", P[j]);
    fprintf(stderr, ")\n");
#endif

    // Compute the reflexion of X[ih] over P
    for(j = 0; j < n; j++)
      Ps[j] = (1. + alpha) * P[j] - alpha * X[ih][j];

    // Evaluate the function in Ps
    ys = fct(Ps, fct_arg);

#ifdef _NELDER_DEBUG
    fprintf(stderr, "%u) Reflexion of x[%u] over P -> Ps = (", iter, ih);
    for(j = 0; j < n; j++) fprintf(stderr, "%s%g", (j > 0) ? ", " : "", Ps[j]);
    fprintf(stderr, "), ys = %g\n", ys);
#endif

    // Check whether we should expand Ps (i.e. ys is a new minimum)
    if(ys < y[il]) {
      // Compute the expansion of Ps over P
      for(j = 0; j < n; j++)
        Pss[j] = gamma * Ps[j] + (1. - gamma) * P[j];

      // Evaluate the function in Pss
      yss = fct(Pss, fct_arg);

#ifdef _NELDER_DEBUG
      fprintf(stderr, "%u) Expansion of Ps over P -> Pss = (", iter);
      for(j = 0; j < n; j++) fprintf(stderr, "%s%g", (j > 0) ? ", " : "", Pss[j]);
      fprintf(stderr, "), yss = %g\n", yss);
#endif
      // If yss < y[il], then replace X[ih] by Pss
      if(yss < y[il]) {
        for(j = 0; j < n; j++) X[ih][j] = Pss[j];
        y[ih] = yss;
#ifdef _NELDER_DEBUG
        fprintf(stderr, "%u) x[%u] = Pss\n", iter, ih);
#endif
      } else { // Otherwise replace keep Ps as the best solution (replace X[ih] by Ps)
        for(j = 0; j < n; j++) X[ih][j] = Ps[j];
        y[ih] = ys;
#ifdef _NELDER_DEBUG
        fprintf(stderr, "%u) x[%u] = Ps\n", iter, ih);
#endif
      }
    } else { // ys is not a new minimum
      // Check whether Ps is a maximum over y[i], i != ih
      for(i = 0; i <= n && (i == ih || y[i] < ys); i++)
        ;
      // If ys is not a maximum, then take this new solution
      if(i <= n) {
        for(j = 0; j < n; j++) X[ih][j] = Ps[j];
        y[ih] = ys;
#ifdef _NELDER_DEBUG
        fprintf(stderr, "%u) x[%u] = Ps\n", iter, ih);
#endif
      } else { // Otherwise, ys is a maximum
        // If y[ih] is lower than ys
        if(y[ih] < ys) {
          // Compute the contraction of X[ih] towards P
          for(j = 0; j < n; j++)
            Pss[j] = beta * X[ih][j] + (1. - beta) * P[j];
          yss = fct(Pss, fct_arg);
#ifdef _NELDER_DEBUG
          fprintf(stderr, "%u) Contraction of x[%u] over P -> Pss = (", iter, ih);
          for(j = 0; j < n; j++) fprintf(stderr, "%s%g", (j > 0) ? ", " : "", Pss[j]);
          fprintf(stderr, "), yss = %g\n", yss);
#endif
        } else { // If ys is lower or equal to y[ih]
          // Compute the contraction of Ps towards P
          for(j = 0; j < n; j++)
            Pss[j] = beta * Ps[j] + (1. - beta) * P[j];
          yss = fct(Pss, fct_arg);
#ifdef _NELDER_DEBUG
          fprintf(stderr, "%u) Contraction of Ps over P -> Pss = (", iter);
          for(j = 0; j < n; j++) fprintf(stderr, "%s%g", (j > 0) ? ", " : "", Pss[j]);
          fprintf(stderr, "), yss = %g\n", yss);
#endif
        }

        // If yss is now lower than y[ih]
        if(yss < y[ih]) {
          // Replace X[ih] by Pss
          for(j = 0; j < n; j++) X[ih][j] = Pss[j];
          y[ih] = yss;
#ifdef _NELDER_DEBUG
          fprintf(stderr, "%u) x[%u] = Pss\n", iter, ih);
#endif
        } else {
          // Shrink all X[i] towards X[il]
          for(i = 0; i <= n; i++) {
            if(i != il) {
              for(j = 0; j < n; j++)
                X[i][j] = (1. - delta) * X[il][j] + delta * X[i][j];
              y[i] = fct(X[i], fct_arg);
            }
          }
#ifdef _NELDER_DEBUG
          fprintf(stderr, "%u) Shrinkage towards x[%u]\n", iter, il);
#endif
        }
      }
    }
#ifdef _NELDER_DEBUG
  fprintf(stderr, "%u) Intermediate simplex:\n", iter);
  for(i = 0; i <= n; i++) {
    fprintf(stderr, "%u)\tx[%u] = (", iter, i);
    for(j = 0; j < n; j++)
      fprintf(stderr, "%s%g", (j > 0) ? ", " : "", X[i][j]);
    fprintf(stderr, "), y[%u] = %g\n", i, y[i]);
  }
  iter++;
#endif
  }

  // Get the index of the minimal point
  il = 0;
  for(i = 1; i <= n; i++) {
    if(y[i] < y[il]) il = i;
  }

  // If il != 0
  if(il != 0) {
    // Swap X[il] with X[0]
    for(j = 0; j < n; j++) {
      ys = X[0][j]; // Use ys as temporary storage
      X[0][j] = X[il][j];
      X[il][j] = ys;
    }
  }

  return y[il];
}
