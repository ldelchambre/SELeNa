#include <graddlsqr.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <pthread.h>
#include <semaphore.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <utils.h>

// #define _GRADDLSQR_DEBUG_MODE // Uncomment to enable verbose debugging

/**
 * Print a debug message if we are in debug mode
 */
#ifdef _GRADDLSQR_DEBUG_MODE
#define _GRADDLSQR_DEBUG(format,...) do{fprintf(stderr, " graddlsqr - " format, ##__VA_ARGS__);fflush(stderr);}while(0)
#else
#define _GRADDLSQR_DEBUG(format,...)
#endif

/*******************************************************************************
 * Threads arguments
 ******************************************************************************/
/**
 * The kind of task to execute
 */
typedef enum _graddlsqr_task {
  _GRADDLSQR_TASK_NONE, // Nothing has to be done
  _GRADDLSQR_TASK_COST, // Compute the cost function
  _GRADDLSQR_TASK_GRAD, // Compute gradient
  _GRADDLSQR_TASK_QUIT  // Exit
} _graddlsqr_task_t;

/**
 * The structure to pass to the threads
 */
typedef struct _graddlsqr_thread_arg {
  // graddlsqr function parameter
  double *a;
  unsigned long na;
  const double *X;
  const double *y;
  const double *w;
  unsigned long ninst;
  unsigned long nattr;
  graddlsqr_fct_t f;
  graddlsqr_deriv_t dfda;
  void *arg;
  const graddlsqr_opt_t *opt;

  // The array of indices of the observations
  unsigned long *idx;

  // The starting and ending indices from idx defining the chunk that is currently being processed
  unsigned long ifrom;
  unsigned long ito;

  // The task to be executed
  _graddlsqr_task_t task;

  // Synchronization variables
  unsigned long ijob; // The next job that should be executed
  unsigned long njob; // The total number of jobs

  // Synchronization variables
  pthread_mutex_t lock;
  sem_t *job_done;
  sem_t *job_waiting;

  // The arrays where the outcome values shall be stored
  double *grads; // size: njob * na
  double *C; // size: njob
} _graddlsqr_thread_arg_t;

/**
 * Initializer of the structure to pass to the threads
 */
#define _GRADDLSQR_THREAD_ARG_INIT { \
    /*           a = */ NULL, \
    /*          na = */ 0, \
    /*           X = */ NULL, \
    /*           y = */ NULL, \
    /*           w = */ NULL, \
    /*       ninst = */ 0, \
    /*       nattr = */ 0, \
    /*           f = */ NULL, \
    /*        dfda = */ NULL, \
    /*         arg = */ NULL, \
    /*         opt = */ NULL, \
    /*         idx = */ NULL, \
    /*       ifrom = */ 0, \
    /*         ito = */ 0, \
    /*        task = */ _GRADDLSQR_TASK_NONE, \
    /*        ijob = */ 0, \
    /*        njob = */ 0, \
    /*        lock = */ PTHREAD_MUTEX_INITIALIZER, \
    /*    job_done = */ NULL, \
    /* job_waiting = */ NULL, \
    /*       grads = */ NULL, \
    /*           C = */ NULL \
}

/**
 * Check the options passed to the graddlsqr function
 */
int _graddlsqr_check_opt(const graddlsqr_opt_t *opt) {
  // Check eps (if eps < 0 then the algorithm stops when maxiter is attained)
  if (isfinite(opt->eps) == 0) {
    fprintf(stderr, "The value of eps is not finite in graddlsqr function.");
    return -1;
  }

  // Check B1
  if (isfinite(opt->B1) == 0 || opt->B1 < 0 || opt->B1 >= 1) {
    fprintf(stderr,
        "The value of B1 is not finite, is lower than zero or is greater or equal to one in graddlsqr function.");
    return -1;
  }

  // Check B2
  if (isfinite(opt->B2) == 0 || opt->B2 < 0 || opt->B2 >= 1) {
    fprintf(stderr,
        "The value of B2 is not finite, is lower than zero or is greater or equal to one in graddlsqr function.");
    return -1;
  }

  // Check alpha
  if (isfinite(opt->alpha) == 0) {
    fprintf(stderr, "The value of alpha is not finite in graddlsqr function.");
    return -1;
  }

  // Check batch_size
  if (opt->batch_size == 0) {
    fprintf(stderr, "The batch size is set to zero in graddlsqr function.");
    return -1;
  }

  // Check nthread
  if (opt->nthread == 0) {
    fprintf(stderr,
        "The number of threads to run is set to zero in graddlsqr function.");
    return -1;
  }

  return 0;
}

/**
 * Evaluate the cost function for the observations targ->idx[ijob:njob:end]
 */
double _graddlsqr_thread_cost(_graddlsqr_thread_arg_t *targ,
    const unsigned long ijob) {
  // Obtain cost parameters from targ
  const unsigned long ifrom = targ->ifrom;
  const unsigned long ito = targ->ito;
  const unsigned long *idx = targ->idx;
  const unsigned long nattr = targ->nattr;
  const unsigned long njob = targ->njob;
  const double *w = targ->w;
  const double *a = targ->a;
  const double *X = targ->X;
  const double *y = targ->y;
  graddlsqr_fct_t f = targ->f;
  void *arg = targ->arg;
  // Cost variables
  double err, C = 0;
  unsigned long i, j;

  // Evaluate the cost function for each instance (0.5 factor is discarded here)
  // C = | W * ( f(a;X) - y ) |^2
  for (i = ifrom + ijob; i <= ito; i += njob) {
    j = idx[i];
    err = w[j] * (f(a, X + j * nattr, nattr, arg) - y[j]);
    C += err * err;
  }

  return C;
}

/**
 * Evaluate the gradient for the observations targ->idx[ijob:njob:end]
 */
void _graddlsqr_thread_grad(double *grad, _graddlsqr_thread_arg_t *targ,
    const unsigned long ijob) {
  // Obtain grad parameters from targ
  const unsigned long ifrom = targ->ifrom;
  const unsigned long ito = targ->ito;
  const unsigned long *idx = targ->idx;
  const unsigned long nattr = targ->nattr;
  const unsigned long na = targ->na;
  const unsigned long njob = targ->njob;
  const double *w = targ->w;
  const double *a = targ->a;
  const double *X = targ->X;
  const double *y = targ->y;
  graddlsqr_fct_t f = targ->f;
  graddlsqr_deriv_t dfda = targ->dfda;
  void *arg = targ->arg;
  // Gradient variables
  double tmp, gradf[na];
  unsigned long i, j;

  // Initialize the gradient to 0
  for (i = 0; i < na; i++)
    grad[i] = 0;

  // Compute the gradient for each instance
  // dC/da = W^2 * ( f(a;X) - y ) * df(a;X) / da
  for (i = ifrom + ijob; i <= ito; i += njob) {
    j = idx[i];
    tmp = w[j] * w[j] * (f(a, X + j * nattr, nattr, arg) - y[j]);
    dfda(gradf, a, X + j * nattr, nattr, arg);
    for (j = 0; j < na; j++)
      grad[j] += tmp * gradf[j];
  }
}

/**
 * The main thread function
 */
void *_graddlsqr_thread(void *arg) {
  // Obtain the thread arguments
  _graddlsqr_thread_arg_t *targ = (_graddlsqr_thread_arg_t *) arg;

  // The job index that is currently executed
  unsigned long ijob;

  // Print debug message
  _GRADDLSQR_DEBUG("Thread %012X: Launched\n", (unsigned ) pthread_self());

  // Loop indefinitely
  while (1) {
    // Wait for some jobs to be available
    _GRADDLSQR_DEBUG("Thread %012X: Waiting for a job\n",
        (unsigned ) pthread_self());
    sem_wait(targ->job_waiting);

    // Take in charge the job numbered ijob
    pthread_mutex_lock(&targ->lock);
    ijob = targ->ijob++;
    pthread_mutex_unlock(&targ->lock);

    // Execute the job
    switch (targ->task) {
      case _GRADDLSQR_TASK_NONE:
        // Nothing has to be done
        _GRADDLSQR_DEBUG("Thread %012X: Running an empty job (job = %lu)\n",
            (unsigned ) pthread_self(), ijob);
        break;
      case _GRADDLSQR_TASK_COST:
        // Compute cost
        _GRADDLSQR_DEBUG("Thread %012X: Computing cost (job = %lu)\n",
            (unsigned ) pthread_self(), ijob);
        targ->C[ijob] = _graddlsqr_thread_cost(targ, ijob);
        break;
      case _GRADDLSQR_TASK_GRAD:
        // Compute gradient
        _GRADDLSQR_DEBUG("Thread %012X: Computing gradient (job = %lu)\n",
            (unsigned ) pthread_self(), ijob);
        _graddlsqr_thread_grad(targ->grads + ijob * targ->na, targ, ijob);
        break;
      case _GRADDLSQR_TASK_QUIT:
        // Exit thread
        _GRADDLSQR_DEBUG("Thread %012X: Exiting (job = %lu)\n",
            (unsigned ) pthread_self(), ijob);
        pthread_exit(NULL);
        break;
    }

    // Signal that we finished the job
    _GRADDLSQR_DEBUG("Thread %012X: Job finished\n",
        (unsigned ) pthread_self());
    sem_post(targ->job_done);
  }
}

/**
 * Launch all threads that are waiting for jobs, then wait for all the jobs to be finished
 */
void _graddlsqr_launch_thread(_graddlsqr_thread_arg_t *targ) {
  unsigned long i;

  _GRADDLSQR_DEBUG("Adding %lu job(s)\n", targ->njob);
  targ->ijob = 0;
  for (i = 0; i < targ->njob; i++)
    sem_post(targ->job_waiting); // Add nthread jobs

  _GRADDLSQR_DEBUG("Waiting for %lu job(s) to finish\n", targ->njob);
  for (i = 0; i < targ->njob; i++)
    sem_wait(targ->job_done); // Wait for all jobs to be done

  _GRADDLSQR_DEBUG("All %lu job(s) finished\n", targ->njob);
}

/**
 * Compute the cost function for observations targ->idx[ifrom:ito] using thread(s)
 */
double _graddlsqr_build_cost(_graddlsqr_thread_arg_t *targ,
    const unsigned long ifrom, const unsigned long ito) {
  double C = 0;
  unsigned long i;

  // Compute cost using threads
  targ->task = _GRADDLSQR_TASK_COST;
  targ->ifrom = ifrom;
  targ->ito = ito;
  _graddlsqr_launch_thread(targ);

  // Sum up the cost computed by each thread
  for (i = 0; i < targ->njob; i++)
    C += targ->C[i];

  // Don't forget the 0.5 term
  return 0.5 * C;
}

/**
 * Compute the gradient for observations targ->idx[ifrom:ito] using thread(s)
 */
void _graddlsqr_build_grad(double *grad, _graddlsqr_thread_arg_t *targ,
    const unsigned long ifrom, const unsigned long ito) {
  const unsigned long na = targ->na;
  unsigned long i, j;

  // Compute gradient using threads
  targ->task = _GRADDLSQR_TASK_GRAD;
  targ->ifrom = ifrom;
  targ->ito = ito;
  _graddlsqr_launch_thread(targ);

  // Initialize the gradient to zero
  for (i = 0; i < na; i++)
    grad[i] = 0;

  // Sump up the gradient computed by each thread
  for (i = 0; i < targ->njob; i++) {
    const double *g = targ->grads + i * na;
    for (j = 0; j < na; j++)
      grad[j] += g[j];
  }

  // FIXME Scale the gradient according to the weight of the processed instances and total weights
}

/**
 * Main function of the gradient descent algorithm
 */
int _graddlsqr_build(_graddlsqr_thread_arg_t *targ) {
  // The options of the gradient descent algorithm
  const graddlsqr_opt_t *opt = targ->opt;

  // The index of the last instance we have to process
  const unsigned long ito =
      (opt->batch_size < targ->ninst) ? opt->batch_size - 1 : targ->ninst - 1;

  // The number of attributes we have
  const unsigned long na = targ->na;

  // The number of iteration performed so far
  unsigned long iter = 0;

  // The random seed to use in order to shuffle the list of indices
  unsigned int seed = 0;

  // The current and previous value of the cost function
  double C, Cprev = DBL_MAX;

  // The current gradient, the first and second momentum vector
  double g[na], m[na], v[na];

  // The t power of B1 and B2
  double B1t = opt->B1, B2t = opt->B2;

  // General usage counter
  unsigned long i;

  // Initialize the momentum vectors
  for (i = 0; i < na; i++) {
    m[i] = 0;
    v[i] = 0;
  }

  // Compute the initial cost function
  C = _graddlsqr_build_cost(targ, 0, targ->ninst - 1);

  // In verbose mode, print the initial solution
  if (opt->verbose) {
    fprintf(stderr, "Initial solution: (");
    for (i = 0; i < na; i++)
      fprintf(stderr, " %g", targ->a[i]);
    fprintf(stderr, " ), Cost: %g\n", C);
  }

  // Main loop
  while (iter < opt->maxiter && fabs(C - Cprev) > opt->eps) {
    // In online mode, select ito + 1 random indices and store them
    // as first elements of targ->idx
    if (opt->batch_size < targ->ninst) {
      unsigned long j, itmp;
      for (i = 0; i <= ito; i++) {
        j = rand_r(&seed) % (targ->ninst - i) + i; // j = [i, ninst[
        itmp = targ->idx[i];
        targ->idx[i] = targ->idx[j];
        targ->idx[j] = itmp;
      }
    }

    // Compute the gradient based on instances [0, ito]
    _graddlsqr_build_grad(g, targ, 0, ito);

    // Update the biased momentum vectors
    for (i = 0; i < na; i++) {
      m[i] = opt->B1 * m[i] + (1. - opt->B1) * g[i];
      v[i] = opt->B2 * v[i] + (1. - opt->B2) * g[i] * g[i];
    }

    // Update parameters
    for (i = 0; i < na; i++)
      targ->a[i] -= opt->alpha * sqrt(1. - B2t) / (1. - B1t) * m[i]
          / (sqrt(v[i]) + opt->epsilon);

    // Print the computed gradient, update vector and new solution in debug mode
#ifdef _GRADDLSQR_DEBUG_MODE
    fprintf(stderr, " graddlsqr -         Gradient : (");
    for (i = 0; i < na; i++)
    fprintf(stderr, " %+g", g[i]);
    fprintf(stderr, " )\n graddlsqr - Updated solution : (");
    for (i = 0; i < na; i++)
    fprintf(stderr, " %+g", targ->a[i]);
    fprintf(stderr, " )\n");
    fflush(stderr);
#endif

    // Compute the new cost function and save the previous one
    Cprev = C;
    C = _graddlsqr_build_cost(targ, 0, targ->ninst - 1);

    // In verbose mode print progress informations
    if (opt->verbose) {
      fprintf(stderr, "%lu) Solution: (", iter);
      for (i = 0; i < na; i++)
        fprintf(stderr, " %+g", targ->a[i]);
      fprintf(stderr, " ), Cost: %g, Previous cost: %g, Difference: %g\n", C,
          Cprev, Cprev - C);
    }

    // Go to the next iteration
    iter++;

    // Compute the iter power of B1, B2
    B1t *= opt->B1;
    B2t *= opt->B2;
  }

  return 0;
}

/**
 * Dstroy the thread arguments
 */
void _graddlsqr_thread_arg_destroy(_graddlsqr_thread_arg_t *targ) {
  free(targ->idx);
  free(targ->C);
  pthread_mutex_destroy(&targ->lock);
  if (targ->job_done)
    sem_destroy(targ->job_done);
  if (targ->job_waiting)
    sem_destroy(targ->job_waiting);
}

/**
 * Initialize the thread arguments, launch threas then destroy the threads arguments
 */
int _graddlsqr(double *a, const unsigned long na, const double *X,
    const double *y, const double *w, const unsigned long ninst,
    const unsigned long nattr, graddlsqr_fct_t f, graddlsqr_deriv_t dfda,
    void *arg, const graddlsqr_opt_t *opt) {
  // The threads identifiers
  pthread_t tid[opt->nthread];

  // Thread arguments
  _graddlsqr_thread_arg_t targ = _GRADDLSQR_THREAD_ARG_INIT;

  // The array containing all gradients from the threads + the values of cost
      double *gradsC;

  // The semaphores we will use
      sem_t job_done,
  job_waiting;

  // General usage counters
  unsigned long i, j;

  // Initialize the common arguments with the given parameters
  targ.a = a;
  targ.na = na;
  targ.X = X;
  targ.y = y;
  targ.w = w;
  targ.ninst = ninst;
  targ.nattr = nattr;
  targ.f = f;
  targ.dfda = dfda;
  targ.arg = arg;
  targ.opt = opt;
  targ.njob = opt->nthread;

  // Allocate the array of indices and the array containing the gradients from all threads
  targ.idx = (unsigned long *) calloc(ninst, sizeof(unsigned long));
  if (targ.idx == NULL) {
    fprintf(stderr,
        "Unable to allocate the array of indices in graddlsqr function\n");
    _graddlsqr_thread_arg_destroy(&targ);
    return -1;
  }

  // Initialize the array of indices
  for (i = 0; i < ninst; i++)
    targ.idx[i] = i;

  // Allocate the array of gradients and cost value
  gradsC = (double *) calloc((na + 1) * targ.njob, sizeof(double));
  if (gradsC == NULL) {
    fprintf(stderr,
        "Unable to allocate the array of gradients and cose in graddlsqr function\n");
    _graddlsqr_thread_arg_destroy(&targ);
    return -1;
  }
  targ.grads = gradsC + targ.njob;
  targ.C = gradsC;

  // Initialize the number of jobs currently done and the number of waiting jobs
  if (sem_init(&job_done, 0,
      0) || (targ.job_done = &job_done) == NULL
      || sem_init(&job_waiting, 0, 0) || (targ.job_waiting = &job_waiting) == NULL) {
    fprintf(stderr, "Unable to initialize semaphores in graddlsqr function\n");
    _graddlsqr_thread_arg_destroy(&targ);
    return -1;
  }

  // Launch all threads
  _GRADDLSQR_DEBUG("Launching %lu threads\n", opt->nthread);
  for (i = 0; i < opt->nthread; i++) {
    if (pthread_create(&tid[i], NULL, _graddlsqr_thread, &targ)) {
      fprintf(stderr, "Unable to create thread in graddlsqr function\n");
      for (j = 0; j < i; j++)
        pthread_cancel(tid[j]);
      _graddlsqr_thread_arg_destroy(&targ);
      return -1;
    }
  }

  // Launch the gradient descent algorithm
  if (_graddlsqr_build(&targ)) {
    for (i = 0; i < opt->nthread; i++)
      pthread_cancel(tid[i]);
    _graddlsqr_thread_arg_destroy(&targ);
    return -1;
  }

  // Cancel all threads
  targ.task = _GRADDLSQR_TASK_QUIT;
  for (i = 0; i < opt->nthread; i++)
    sem_post(targ.job_waiting);
  for (i = 0; i < opt->nthread; i++)
    pthread_join(tid[i], NULL);

  // Destroy synchronization variables and arrays
  _graddlsqr_thread_arg_destroy(&targ);

  return 0;
}

int graddlsqr(double *a, const unsigned long na, const double *X,
    const double *y, const double *w, const unsigned long ninst,
    const unsigned long nattr, graddlsqr_fct_t f, graddlsqr_deriv_t dfda,
    void *arg, const graddlsqr_opt_t *opt) {
  // If opt is NULL, then use the default options
  if (opt == NULL) {
    const graddlsqr_opt_t default_opt = GRADDLSQR_DEFAULT_OPTS;
    return graddlsqr(a, na, X, y, w, ninst, nattr, f, dfda, arg, &default_opt);
  } else if (_graddlsqr_check_opt(opt)) // Otherwise, check options
    return -1;

  // Launch the gradient descent algorithm
  return _graddlsqr(a, na, X, y, w, ninst, nattr, f, dfda, arg, opt);
}
