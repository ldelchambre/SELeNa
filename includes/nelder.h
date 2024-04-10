/**
 * @file nelder.h
 * 
 * Nelder-Mead optimization algorithm, header file.
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
#ifndef _NELDER_H_
#define _NELDER_H_

/**
 * The default reflexion coefficient of the Nelder-Mead algorithm.
 * 
 * @see Nelder, J. A. and Mead, R., "A simplex method for Function minimization", 1965, The Computer Journal, 7(4), 308, doi:10.1093/comjnl/7.4.308
 */
#define NELDER_DEFAULT_ALPHA 1.0

/**
 * The default contraction coefficient of the Nelder-Mead algorithm.
 * 
 * @see Nelder, J. A. and Mead, R., "A simplex method for Function minimization", 1965, The Computer Journal, 7(4), 308, doi:10.1093/comjnl/7.4.308
 */
#define NELDER_DEFAULT_BETA  0.5

/**
 * The default expansion coefficient of the Nelder-Mead algorithm.
 * 
 * @see Nelder, J. A. and Mead, R., "A simplex method for Function minimization", 1965, The Computer Journal, 7(4), 308, doi:10.1093/comjnl/7.4.308
 */
#define NELDER_DEFAULT_GAMMA 2.0

/**
 * The default shrinkage coefficient of the Nelder-Mead algorithm.
 * 
 * @see Nelder, J. A. and Mead, R., "A simplex method for Function minimization", 1965, The Computer Journal, 7(4), 308, doi:10.1093/comjnl/7.4.308
 */
#define NELDER_DEFAULT_DELTA 0.5

/**
 * The default standard deviation of the values of the function evaluated at the
 * vertices of the simplex below which the Nelder-Mead algorithm is considered 
 * to have reached the convergence.
 * 
 * @see Nelder, J. A. and Mead, R., "A simplex method for Function minimization", 1965, The Computer Journal, 7(4), 308, doi:10.1093/comjnl/7.4.308
 */
#define NELDER_DEFAULT_EPS   1e-15

/**
 * Minimize a function of n input parameters using a Nelder-Mead algorithm 
 * along with a starting simplex given by X.
 * 
 *   On input, X should contains a set of (n+1) vertices of dimension n taken in
 * the input space of parameters (i.e. X[0:n-1] contains the first point, 
 * X[n:2*n-1] contain the second point and so on). Each point of X must be 
 * distinct in order for the algorithm to work. On output, X will contain the 
 * vertices of the simplex used in the last iteration of the algorithm, more 
 * particularly X[0:n-1] will contain the vertex for which the function fct is 
 * minimized, other vertices are unordered.
 * 
 * The function to minimize must have a prototype of the form
 * > double fct(const double *x, void *arg);
 * where x is an input array of size n and where arg is an user input argument
 * given by fct_arg. This function should return the value against which the
 * minimization should be performed.
 * 
 * The Nelder-Mead algorithm stops iterating when the standard deviation of the
 * values of the function taken in each vertices of the simplex X falls below 
 * eps.
 * 
 * @param X       The input simplex of size n x (n+1), in column-major order and final simplex on output..
 * @param n       The number of input parameters.
 * @param fct     The function to minimize.
 * @param fct_arg The arguments to pass to fct
 * @param alpha   The reflexion coefficient of the Nelder-Mead algorithm (default NELDER_DEFAULT_ALPHA)
 * @param beta    The contraction coefficient of the Nelder-Mead algorithm (default NELDER_DEFAULT_BETA)
 * @param gamma   The expansion coefficient of the Nelder-Mead algorithm (default NELDER_DEFAULT_GAMMA)
 * @param delta   The shrinkage coefficient of the Nelder-Mead algorithm (default NELDER_DEFAULT_DELTA)
 * @param eps     The standard deviation of the values of the function below which the algorithm stops
 * 
 * @return The minimal value of the function we found (evaluated in X[0:n-1])
 * 
 * @see Nelder, J. A. and Mead, R., "A simplex method for Function minimization", 1965, The Computer Journal, 7(4), 308, doi:10.1093/comjnl/7.4.308
 */
double nelder(double *X, const unsigned int n,
              double (*fct)(const double *, void *), void *fct_arg,
              const double alpha, const double beta, const double gamma, const double delta,
              const double eps);

/**
 * Minimize a function of n input parameters using a Nelder-Mead algorithm 
 * along with a starting simplex given by X.
 * 
 *   On input, X should contains a set of (n+1) vertices of dimension n taken in
 * the input space of parameters (i.e. X[0:n-1] contains the first point, 
 * X[n:2*n-1] contain the second point and so on). Each point of X must be 
 * distinct in order for the algorithm to work. On output, X will contain the 
 * vertices of the simplex used in the last iteration of the algorithm, more 
 * particularly X[0:n-1] will contain the vertex for which the function fct is 
 * minimized, other vertices are unordered.
 * 
 * The function to minimize must have a prototype of the form
 * > double fct(const double *x, void *arg);
 * where x is an input array of size n and where arg is an user input argument
 * given by fct_arg. This function should return the value against which the
 * minimization should be performed.
 * 
 * The Nelder-Mead algorithm stops iterating when end_fct return something 
 * different than zero. The end_fct function has a prototype given by
 * > int end_fct(const double *X, const double *y, void *arg);
 * where X is the last computed simplex, y are the values of the function 
 * evaluated in each vertices of the simplex (it thus have a size of n+1 
 * elements) and arg is a user-provided arguments corresponding to end_fct_arg.
 * 
 * @param X       The input simplex of size n x (n+1), in column-major order and final simplex on output..
 * @param n       The number of input parameters.
 * @param fct     The function to minimize.
 * @param fct_arg The arguments to pass to fct
 * @param alpha   The reflexion coefficient of the Nelder-Mead algorithm (default NELDER_DEFAULT_ALPHA)
 * @param beta    The contraction coefficient of the Nelder-Mead algorithm (default NELDER_DEFAULT_BETA)
 * @param gamma   The expansion coefficient of the Nelder-Mead algorithm (default NELDER_DEFAULT_GAMMA)
 * @param delta   The shrinkage coefficient of the Nelder-Mead algorithm (default NELDER_DEFAULT_DELTA)
 * @param end_fct The function which determine whether the Nelder-Mead algorithm should stop iterating or not
 * @param end_fct_arg The argument to pass to the end_fct function
 * 
 * @return The minimal value of the function we found (evaluated in X[0:n-1])
 * 
 * @see Nelder, J. A. and Mead, R., "A simplex method for Function minimization", 1965, The Computer Journal, 7(4), 308, doi:10.1093/comjnl/7.4.308
 */
double nelder2(double *X, const unsigned int n,
              double (*fct)(const double *, void *), void *fct_arg,
              const double alpha, const double beta, const double gamma, const double delta,
              int (*end_fct)(const double *, const double *, void *), void *end_fct_arg);

#endif
