#ifndef _GRADDLSQR_H_
#define _GRADDLSQR_H_

/*******************************************************************************
 * Macros
 ******************************************************************************/

/**
 * Default gradient descent options
 */
#define GRADDLSQR_DEFAULT_OPTS { \
  /*        eps = */ 1e-6, \
  /*      alpha = */ 0.001, \
  /*         B1 = */ 0.9, \
  /*         B2 = */ 0.999, \
  /*    epsilon = */ 1e-8, \
  /* batch_size = */ -1UL, \
  /*    nthread = */ 1,\
  /*    maxiter = */ 10000, \
  /*    verbose = */ 0 \
}

/*******************************************************************************
 * Typedefs
 ******************************************************************************/

/**
 * Structure containing the options of the stochastic gradient descent algorithm
 * with momentum.
 */
typedef struct graddlsqr_opt {
         double eps;        // The required precision
         double alpha;      // The step size/learning rate
         double B1;         // Exponential decay rate for the first moment estimate
         double B2;         // Exponential decay rate for the second moment estimate
         double epsilon;    // Small number designed to prevent division by 0
  unsigned long batch_size; // The batch size
  unsigned long nthread;    // The number of threads to run
  unsigned long maxiter;    // The maximal number of iterations to perform
           char verbose;    // Should we have to be verbose?
} graddlsqr_opt_t;

/**
 * Prototype of the function to optimize
 */
typedef double (*graddlsqr_fct_t)(     const double * /* a */,
                                       const double * /* x */,
                                  const unsigned long /* nattr */,
                                               void * /* arg */);

/**
 * Prototype of the derivative of the function to optimize
 */
typedef void   (*graddlsqr_deriv_t)(           double * /* grad */,
                                         const double * /* a */,
                                         const double * /* x */,
                                    const unsigned long /* nattr */,
                                                 void * /* arg */);

/*******************************************************************************
 * Function prototypes
 ******************************************************************************/

/**
 * Minimize the cost function
 *
 * C = 0.5 * | diag(w) * ( f(X;a) - y ) |^2                                  (1)
 *
 * through stochastic gradient descent.
 *
 * In the latter equation:
 *   w is the vector of weights associated with each observations (size: ninst).
 *   a is the vector of parameters we have to optimize (size: nattr).
 *   X is the matrix of input attributes associated with each observation in
 *     column-major order (size: nattr x ninst).
 *   f is the function we would like to fit.
 *   y is the vector of outcome values (size: ninst).
 *
 *   The stochastic gradient descent uses the Adam algorithm described in
 * Kingma & Ba, "ADAM: A METHOD  FOR STOCHASTIC OPTIMIZATION", ICLR 2015,
 * arXiv:1412.6980.
 *
 * ---- ALGORITHM --------------------------------------------------------------
 * > C(-1) = 0;
 * > m(0) = 0; { Initialize the first momentum vector }
 * > v(0) = 0; { Initialize the second momentum vector }
 * > t = 0;
 * > while | C(t) - C(t-1) | >= eps && t < maxiter
 * >  t = t + 1;
 * >  { Compute the gradient of C wrt. a}
 * >  g(t) = dC(t) / da;
 * >  { Update the biased first moment estimate }
 * >  m(t) = B1 * m(t-1) + (1-B1) * g(t);
 * >  { Update the biased second moment estimate }
 * >  v(t) = B2 * v(t-1) + (1-B2) * g(t)^2;
 * >  { Update parameters }
 * >                          sqrt( 1 - B2^t )             m(t)
 * >  a(t) = a(t-1) - alpha * ---------------- * ------------------------
 * >                              1 - B1^t        sqrt( v(t) ) + epsilon
 * -----------------------------------------------------------------------------
 *
 *   At each step of the algorithm, the equation (1) is evaluated over a random
 * subset of observations, batch_size.
 *
 * @param a     The parameters to be fitted (size: na).
 *              On input a contain the initial solution, a(0).
 *              On output a contains the final solution, a(t).
 * @param na    The size of the vector a.
 * @param X     The matrix of observations in column-major order (size: ninst x nattr).
 * @param y     The vector of outcome values (size: ninst).
 * @param w     The vector of weights associated with each observations (size: ninst).
 * @param ninst The number of instances.
 * @param nattr The number of attributes.
 * @param f     The function that should be fitted to the observations, f(X;a).
 * @param dfda  The gradient function of f with respect to a, df(X;a) / da.
 * @param arg   User-defined parameters to pass to the f and dfda functions.
 * @param opt   The options of the stochastich gradient descent algorithm.
 *              If opt == NULL, then the default options will be used, though these
 *              might be inappropriate.
 *              opt contains the following fields:
 *                - eps        The tolerance in the change of the cost function between
 *                             each iteration of the algorithm.
 *                - alpha      The step size or learning rate of the algorithm.
 *                - B1         The exponential decay rate for the first moment estimates
 *                - B2         The exponential decay rate for the seconf moment estimates
 *                - epsilon    A very small number to prevent any division by zero
 *                - batch_size The number of random observation to use at each iteration
 *                             of the algorithm. If batch_size >= ninst, then all
 *                             observations will be used.
 *                - nthread    The number of threads to use.
 *                - maxiter    The maximum number of iterations to perform.
 *                - verbose    Print intermediate result to stderr if set.
 *
 *  @return -1 in case of error, 0 otherwise
 */
int graddlsqr(               double * a,
                  const unsigned long na,
                       const double * X,
                       const double * y,
                       const double * w,
                  const unsigned long ninst,
                  const unsigned long nattr,
                      graddlsqr_fct_t f,
                    graddlsqr_deriv_t dfda,
                               void * arg,
              const graddlsqr_opt_t * opt);

#endif
