#include <treealgorithm.h>
#include <circstat.h>
#include <information.h>
#include <stat.h>
#include <utils.h>

#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/**
 * Comment in order to disable the debug mode
 */
// #define _TREEALGORITHM_DEBUG_MODE

/**
 * In debug mode, _TREEALGORITHM_DEBUG will print a debug message
 */
#ifdef _TREEALGORITHM_DEBUG_MODE
#define _TREEALGORITHM_DEBUG(format,...) do{fprintf(stdout, format, ##__VA_ARGS__);fflush(stdout);}while(0)
#else
#define _TREEALGORITHM_DEBUG(format,...)
#endif

/**
 * Common error handling
 */
#define _TREEALGORITHM_ERROR(format, ...)\
do {\
  int errsv = errno; \
  fprintf(stderr, "Tree - Error while "action": ", ##__VA_ARGS__); \
  if(errsv) \
    fprintf(stderr, "%s", strerror(errsv)); \
  fprintf(stderr, "\n");\
  fflush(stderr);\
} while(0)

/******************************************************************************
 * GENERAL PURPOSE UTILITY FUNCTIONS                                          *
 ******************************************************************************/
size_t tree_elligible_attributes(size_t *eattr, const double *_X, const void *_y, const size_t *inst,
                                                const size_t ninst, const size_t nattr, const size_t ysize,
                                                int (*ycmp)(const void *, const void *, void *), void *ycmparg) {
  const double (*X)[nattr] = (const double (*)[nattr]) _X;
  const char (*y)[ysize] = (const char (*)[ysize]) _y;
  size_t neattr, i, j;

  // Check if all y are equals
  for(i = 1; i < ninst && ycmp(&y[inst[0]], &y[inst[i]],ycmparg) == 0; i++)
    ;
  // If all y are equals, then no attribute is elligible
  if(i == ninst)
    return 0;

  // Get the elligible attributes
  neattr = 0; // Intially, no attribute is elligible
  for(i = 0; i < nattr; i++) {
    // Find non-constant attributes
    for(j = 1; j < ninst && X[inst[0]][i] == X[inst[j]][i]; j++)
      ;
    if(j < ninst)
      eattr[neattr++] = i;
  }

  return neattr;
}

size_t tree_best_split(double *cutval, const double *_X, const void *_y, const size_t *inst,
                                       const size_t ninst, const size_t nattr, const size_t ysize,
                                       split_quality_fct split_quality, int (*ycmp)(const void *, const void *, void *),
                                       void *ycmparg, void *arg) {
  const double (*X)[nattr] = (const double (*)[nattr]) _X; // The input dataset
  const char (*y)[ysize] = (const char (*)[ysize]) _y; // The input classes
  size_t eattr[nattr]; // The elligible attributes
  size_t neattr; // The number of elligible attributes

  double *v; // The values of the subset of attributes
  size_t *is; // The set of indices that would sort v
  char (*ys)[ysize]; // The set of classes sorted according to their associated attributes

  size_t best_iattr;
  double best_cutval;
  double best_score;
  double s;
  size_t i,j;

  // Get the set of elligible attributes
  neattr = tree_elligible_attributes(eattr, _X, _y, inst, ninst, nattr, ysize, ycmp, ycmparg);

  // If we don't have any elligible attribute, we should create a leaf node
  if(neattr == 0)
    return nattr;

  // Allocate the various arrays
  v = (double *) calloc(ninst, sizeof(double));
  is = (size_t *) calloc(ninst, sizeof(size_t));
  ys = (char (*)[ysize]) calloc(ninst, ysize);
  if(v == NULL || is == NULL || ys == NULL) {
    free(v); free(is); free(ys);
    return -1UL;
  }

  // For each elligible attribute
  best_score = -DBL_MAX;
  for(i = 0; i < neattr; i++) { // XXX We might use multiprocessing in order to fasten this part (1 thread per attribute)
    // Copy the subset of attribute values within vattr
    for(j = 0; j < ninst; j++)
      v[j] = X[inst[j]][eattr[i]];

    // Find the set of indices that would sort v (i.e. v[is[0...ninst]] is sorted)
    argsort(is, v, ninst, sizeof(double), dbl_compar);

    // Copy the classes sorted according to their associated current attribute
    for(j = 0; j < ninst; j++)
      memcpy(&ys[j], &y[inst[is[j]]], ysize);

    // Find the initial split value
    for(j = 1; j < ninst && v[is[j-1]] == v[is[j]]; j++)
      ;

    // Browse all cut point within this attribute
    while(j < ninst) {
      // Compute the score associated with this cut point
      s = split_quality(&ys[0], j, &ys[j], ninst-j, arg);

      // If the returned score is -DBL_MAX, then an error occurs
      if(s == -DBL_MAX) {
        free(v); free(is); free(ys);
        return -1UL;
      }

      // Check whether this split is the best one
      if(best_score < s) {
        best_score = s;
        best_iattr = eattr[i];
        best_cutval = 0.5 * (v[is[j-1]]+v[is[j]]);;
      }

      // Go to the next split value
      for(j++; j < ninst && v[is[j-1]] == v[is[j]]; j++)
        ;
    }
  }

  // Delete arrays
  free(v);free(is);free(ys);

  // Return the best split attribute and its value
  *cutval = best_cutval;
  return best_iattr;
}

void *tree_flatten_y(const void *y, const size_t *inst, const size_t ninst, const size_t ysize) {
  const char (*yin)[ysize] = (const char (*)[ysize]) y;
  char (*yout)[ysize] = (char (*)[ysize]) calloc(ninst, ysize);

  if(y == NULL)
    return NULL;
  for(size_t i = 0; i < ninst; i++)
      memcpy(&yout[i], &yin[inst[i]], ysize);
  return yout;
}

/******************************************************************************
 * REGRESSION TREE ALGORITHM AND UTILITY FUNCTIONS                            *
 ******************************************************************************/
double regressiontree_split_quality(const void *S1, const size_t n1, const void *S2, const size_t n2, void *arg) {
  const double var1 = stat_var((const double *) S1, n1);
  const double var2 = stat_var((const double *) S2, n2);
  return -(var1*n1 + var2*n2) / (n1 + n2);
}

int _regressiontree_algorithm_init(const double *X, const void *y, const size_t ninst, const size_t nattr, void *arg) {
  _TREEALGORITHM_DEBUG("Regression tree initialization\n");
  return 0;
}

/**
 * Find the best split amongst all attributes and all instances regarding the
 * numerical attributes contained within _y.
 */
size_t _regressiontree_algorithm_getsplit(double *cutval, const double *X, const void *y, const size_t *inst, const size_t ninst, const size_t nattr, void *arg) {
  return tree_best_split(cutval, X, y, inst, ninst, nattr, sizeof(double), regressiontree_split_quality, dbl_compar_r, NULL, arg);
}

void *_regressiontree_algorithm_getdata(size_t *data_size, const void *_y, const size_t *inst,  const size_t ninst, void *arg) {
  double *y;
  double *res;
  double var;

  // Flatten the y array
  y = (double *) tree_flatten_y(_y, inst, ninst, sizeof(double));
  if(y == NULL)
    return NULL;

  // Allocate the result array
  res = (double *) calloc(2, sizeof(double));
  if(res == NULL) {
    free(y);
    return NULL;
  }

  // Associated data are the mean value along with its associated standard deviation
  *data_size = 2*sizeof(double);
  res[0] = stat_mean2(&var, y, ninst);
  res[1] = sqrt(var);

  _TREEALGORITHM_DEBUG("data = [ninstances=%zu, mean = %f, std = %f]\n", ninst, res[0], res[1]);

  // Delete the flattenedd y array
  free(y);

  return res;
}

const tree_algorithm regressiontree_algorithm = {
  /* init = */       _regressiontree_algorithm_init,
  /* getsplit = */   _regressiontree_algorithm_getsplit,
  /* getdata =  */   _regressiontree_algorithm_getdata,
  /* arg = */        NULL
};

/******************************************************************************
 * PERIODIC REGRESSION TREE ALGORITHM AND UTILITY FUNCTIONS                   *
 ******************************************************************************/
double pregressiontree_split_quality(const void *S1, const size_t n1, const void *S2, const size_t n2, void *arg) {
  const double var1 = circstat_var((const double *) S1, n1);
  const double var2 = circstat_var((const double *) S2, n2);
  return -(var1*n1 + var2*n2) / (n1 + n2);
}

int _pregressiontree_algorithm_init(const double *X, const void *y, const size_t ninst, const size_t nattr, void *arg) {
  _TREEALGORITHM_DEBUG("Periodic regression tree initialization\n");
  return 0;
}

size_t _pregressiontree_algorithm_getsplit(double *cutval, const double *X, const void *y, const size_t *inst, const size_t ninst, const size_t nattr, void *arg) {
  return tree_best_split(cutval, X, y, inst, ninst, nattr, sizeof(double), pregressiontree_split_quality, dbl_compar_r, NULL, arg);
}

void *_pregressiontree_algorithm_getdata(size_t *data_size, const void *_y, const size_t *inst,  const size_t ninst, void *arg) {
  const double sigmaconf = (arg == NULL) ? 2. : *((const double *) arg);
  double *y;
  double *res;

  // Flatten the y array
  y = (double *) tree_flatten_y(_y, inst, ninst, sizeof(double));
  if(y == NULL)
    return NULL;

  // Allocate the result array
  res = (double *) calloc(2, sizeof(double));
  if(res == NULL) {
    free(y);
    return NULL;
  }

  // Associated data are the mean value along with its associated standard deviation
  *data_size = 2*sizeof(double);
  res[0] = circstat_mean(y,ninst);
  res[1] = circstat_confmean(y,ninst,sigmaconf);

  _TREEALGORITHM_DEBUG("data = [ninstances=%zu, mean = %f, %.2f%% confidence interval = %f]\n", ninst, res[0], 100.*erf(sigmaconf/M_SQRT2), res[1]);

  // Delete the flattenedd y array
  free(y);

  return res;
}

double pregression_sigmaconf = 2.;

const tree_algorithm pregressiontree_algorithm = {
  /* init = */       _pregressiontree_algorithm_init,
  /* getsplit = */   _pregressiontree_algorithm_getsplit,
  /* getdata =  */   _pregressiontree_algorithm_getdata,
  /* arg = */        &pregression_sigmaconf
};

/******************************************************************************
 * CLASSIFICATION TREE ALGORITHM AND UTILITY FUNCTIONS                        *
 ******************************************************************************/
int classificationtree_measure = CLASSIFICATION_MEASURE_IG;

double classificationtree_split_quality(const void *S1, const size_t n1, const void *S2, const size_t n2, void *arg) {
  const int *p = (const int *) arg;
  double i1, i2;
  int e1, e2;

  // Compute the impurity os each subset, based on arg
  if(p == NULL || *p == CLASSIFICATION_MEASURE_IG) {
    e1 = information_entropy(&i1, (const int *) S1, n1);
    e2 = information_entropy(&i2, (const int *) S2, n2);
  } else if(*p == CLASSIFICATION_MEASURE_GINI) {
    e1 = information_gini(&i1, (const int *) S1, n1);
    e2 = information_gini(&i2, (const int *) S2, n2);
  } else {
    e1 = information_entropy(&i1, (const int *) S1, n1);
    e2 = information_entropy(&i2, (const int *) S2, n2);
  }

  // If one of the impurity function fails, signal error
  if(e1 || e2)
    return -DBL_MAX;

  return -(i1*n1+i2*n2)/(n1+n2);
}

int _classificationtree_algorithm_init(const double *X, const void *y, const size_t ninst, const size_t nattr, void *arg) {
  _TREEALGORITHM_DEBUG("Classification tree initialization\n");
  return 0;
}

/**
 * Find the best split amongst all attributes and all instances regarding the
 * categorical attributes contained within _y.
 */
size_t _classificationtree_algorithm_getsplit(double *cutval, const double *X, const void *y, const size_t *inst, const size_t ninst, const size_t nattr, void *arg) {
  return tree_best_split(cutval, X, y, inst, ninst, nattr, sizeof(int), classificationtree_split_quality, int_compar_r, NULL, arg);
}

void *_classificationtree_algorithm_getdata(size_t *data_size, const void *_y, const size_t *inst,  const size_t ninst, void *arg) {
  int *y;
  classificationtree_data *data;
  int mode;
  double p;

  // Flatten the y array
  y = (int *) tree_flatten_y(_y, inst, ninst, sizeof(int));
  if(y == NULL)
    return NULL;

  // Get the mode and associated probability
  if(information_mode(&mode, &p, y, ninst)) {
    free(y);
    return NULL;
  }

  // Allocate the resulting data
  data = calloc(1,sizeof(classificationtree_data));
  if(data == NULL) {
    free(y);
    return NULL;
  }

  // Fill the resulting data
  *data_size = sizeof(classificationtree_data);
  data->c = mode;
  data->pc = p;

  _TREEALGORITHM_DEBUG("data = [ninstances=%zu, class = %d, P(%d) = %f]\n", ninst, data->c, data->c, data->pc);

  // Delete the flatened y array
  free(y);

  return data;
}

const tree_algorithm classificationtree_algorithm = {
  /* init = */       _classificationtree_algorithm_init,
  /* getsplit = */   _classificationtree_algorithm_getsplit,
  /* getdata =  */   _classificationtree_algorithm_getdata,
  /* arg = */        &classificationtree_measure
};
