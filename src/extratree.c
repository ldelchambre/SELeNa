#include <extratree.h>
#include <utils.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Comment in order to disable the debug mode
 */
// #define _EXTRATREE_DEBUG_MODE

/**
 * In debug mode, _THREADPOOL_DEBUG will print a debug message
 */
#ifdef _EXTRATREE_DEBUG_MODE
#define _EXTRATREE_DEBUG(format,...) do{fprintf(stdout, format, ##__VA_ARGS__);fflush(stdout);}while(0)
#else
#define _EXTRATREE_DEBUG(format,...)
#endif

/**
 * Common error handling within the thread pool
 */
#define _EXTRATREE_ERROR(err,message, ...)\
do {\
  int errsv = err; \
  fprintf(stderr, "extratree - "message, ##__VA_ARGS__); \
  if(errsv) \
    fprintf(stderr, ": %s", strerror(errsv)); \
  fprintf(stderr,"\n");\
  fflush(stderr);\
} while(0)

typedef struct _extratree_algorithm_arg_t {
  unsigned int k;
  unsigned int nmin;
  unsigned int seed;

  split_quality_fct split_quality;
  void *split_quality_arg;

  int (*ycmp)(const void *, const void *, void *);
  size_t ysize;
  void *ycmp_arg;
} _extratree_algorithm_arg;

void extratree_default_parameters(unsigned int *k, unsigned int *nmin, const size_t nattr, const int regression) {
  if(k) *k = (regression) ? nattr-1 : rint(sqrt(nattr-1));
  if(nmin) *nmin = (regression) ? 5 : 2;
}

int _extratree_init(const double *X, const void *y, const size_t ninst, const size_t nattr, void *arg) {
#ifdef _EXTRATREE_DEBUG_MODE
  _extratree_algorithm_arg *ert = (_extratree_algorithm_arg *) arg;
  _EXTRATREE_DEBUG("extratree initialization (k=%u, nmin=%u, seed= %u)\n", ert->k, ert->nmin, ert->seed);
#endif

  return 0;
}

size_t _extratree_getsplit(double *cutval, const double *_X, const void *_y, const size_t *inst, const size_t ninst, const size_t nattr, void *arg) {
  _extratree_algorithm_arg *ert = (_extratree_algorithm_arg *) arg; // The algorithm parameters
  const double (*X)[nattr] = (const double (*)[nattr]) _X; // The input dataset
  const char (*y)[ert->ysize] = (const char (*)[ert->ysize]) _y;
  size_t eattr[nattr]; // The set of elligible attributes
  size_t neattr; // The number of elligible attributes

  // Create a leaf node if we reached the minimal number of node necessary in order to split the node
  if(ninst < ert->nmin)
    return nattr;

  // Get the set of elligible attributes
  neattr = tree_elligible_attributes(eattr, _X, _y, inst, ninst, nattr, ert->ysize, ert->ycmp, ert->ycmp_arg);

  // If we don't have any elligible attribute, we should create a leaf node
  if(neattr == 0)
    return nattr;

  {
    char (*ysplit)[ert->ysize]; // The splitted classes
    double rcutval; // The randomly chosen cut val
    double score; // The score of a cut
    double best_score; // The best score encountered
    double best_iattr; // The attribute associated with the best score
    double x, minx, maxx; // The current, minimal and maximal values of the K randomly selected attributes
    unsigned int seed = hashcode(inst, ninst*sizeof(size_t))+ert->seed; // The random seed we will use (use the hash code of inst array as 'unique' identifier)
    size_t i, k, kattr, head, tail;

    // Allocate the array of splitted classes
    ysplit = (char (*)[ert->ysize]) calloc(ninst, ert->ysize);
    if(ysplit == NULL)
      return -1UL;

    // How many attributes should we test (i.e. the actual K value)?
    k = (neattr < ert->k) ? neattr : ert->k;

    // For k randomly selected attributes
    best_score = -DBL_MAX;
    for(; k; k--) {
      // Select a random attribute
      i = mrand_r(&seed) % neattr;
      kattr = eattr[i];
      eattr[i] = eattr[--neattr];

      _EXTRATREE_DEBUG("Selected attribute: %zu\n", kattr);

      // Get the minimal and maximal value of this attribute
      minx = DBL_MAX;
      maxx = -DBL_MAX;
      for(i = 0; i < ninst; i++) {
        x = X[inst[i]][kattr];
        if(x < minx)
          minx = x;
        if(maxx < x)
          maxx = x;
      }

      // Get a random cut-point for this attribute
      rcutval = (maxx-minx)*(1.+mrand_r(&seed))/(1.+MRAND_MAX)+minx; // r = (1.+mrand_r(&seed))/(1.+MRAND_MAX) produce a random number that is such that 0 < r <= 1

      _EXTRATREE_DEBUG("Selected cut value: %f (min=%f, max=%f)\n", rcutval, minx, maxx);

      // Split the input classes according to the found split criterion
      head = 0;
      tail = ninst;
      for(i = 0; i < ninst; i++) {
        x = X[inst[i]][kattr];
        if(x < rcutval)
          memcpy(&ysplit[head++], &y[inst[i]], ert->ysize);
        else
          memcpy(&ysplit[--tail], &y[inst[i]], ert->ysize);
      }
      // Here head == tail

      // Get the score of this split
      score = ert->split_quality(&ysplit[0], head, &ysplit[tail], ninst-tail, ert->split_quality_arg);
      if(best_score < score) {
        best_score = score;
        best_iattr = kattr;
        *cutval = rcutval;
      }

      _EXTRATREE_DEBUG("Score of split X[%zu] < %f: %f\n", kattr, rcutval, score);
    }

    // Delete the array of splitted classes
    free(ysplit);

    return best_iattr;
  }
}

void *_extratree_getdata(size_t *data_size, const void *y, const size_t *inst,  const size_t ninst, void *arg) {
  _extratree_algorithm_arg *ert = (_extratree_algorithm_arg *) arg;
  *data_size = ninst*ert->ysize;
  return tree_flatten_y(y, inst, ninst, ert->ysize);
}

void extratree_algorithm_destroy(tree_algorithm *algorithm) {
  if(algorithm) {
    free(algorithm->arg);
    free(algorithm);
  }
}

tree_algorithm *extratree_algorithm(const unsigned int k,
                                    const unsigned int nmin,
                                    const unsigned int seed,
                                    split_quality_fct split_quality,
                                    void *split_quality_arg,
                                    int (*ycmp)(const void *, const void *, void *),
                                    size_t ysize,
                                    void *ycmp_arg) {
  tree_algorithm *algorithm;
  _extratree_algorithm_arg *arg;

  if(k == 0) {
    _EXTRATREE_ERROR(0, "k parameter is zero\n");
    return NULL;
  }

  if(nmin == 0) {
    _EXTRATREE_ERROR(0, "nmin parameter is zero\n");
    return NULL;
  }

  // Allocated the extratree algorithm and parameter
  algorithm = (tree_algorithm *) malloc(sizeof(tree_algorithm));
  arg = (_extratree_algorithm_arg *) malloc(sizeof(_extratree_algorithm_arg));
  if(algorithm == NULL || arg == NULL) {
    _EXTRATREE_ERROR(ENOMEM, "Could not allocate the extratree algorithm and parameter");
    free(algorithm);free(arg);
    return NULL;
  }

  // Initialize the extratree parameters
  arg->k = k;
  arg->nmin = nmin;
  arg->seed = seed;
  arg->split_quality = split_quality;
  arg->split_quality_arg = split_quality_arg;
  arg->ycmp = ycmp;
  arg->ysize = ysize;
  arg->ycmp_arg = ycmp_arg;

  // Initialize the extratree algorithm
  algorithm->init = _extratree_init;
  algorithm->getsplit = _extratree_getsplit;
  algorithm->getdata = _extratree_getdata;
  algorithm->arg = arg;

  return algorithm;
}
