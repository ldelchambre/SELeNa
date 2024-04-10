#include <float.h>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <circstat.h>
#include <stat.h>
#include <information.h>

#include <tree.h>
#include <treealgorithm.h>
#include <extratree.h>
#include <treecommittee.h>
#include <ExtraTrees.h>

#include <utils.h>

#define NINSTANCES 2500
#define NATTRIBUTES 5
#define NTHREAD 12

/**
 * Function type designed to produced data based on an input dataset
 */
typedef void *(*dataproducer)(const double *X, const size_t ninst, const size_t nattr, unsigned int *seed);
/**
 * Function type designed to make prediction based on a input data set
 */
typedef int (*predictor)(threadpool_t *pool, void *LSypred, void *TSypred, const double *LS, const void *LSy, const double *TS, const size_t ninst, const size_t nattr, void *arg);
/**
 * Function designed to report prediction statistics
 */
typedef void (*reporter)(FILE *fout, const void *y, const void *ypred, const size_t ninst);

/**
 * Produce regression data based on an input dataset
 */
void *regression_data(const double *_X, const size_t ninst, const size_t nattr, unsigned int *seed) {
  const double noise = 0.1; // The (uniform) noise to add to the resulting regression data
  const double (*X)[nattr] = (const double (*)[nattr]) _X;
  double *y = (double *) calloc(ninst, sizeof(double));
  double a[nattr];
  double ymin, ymax;
  size_t i, j;

  // If allocation fails, return NULL
  if(y == NULL)
    return NULL;

  // Initialize the linear coefficient
  for(i = 0; i < nattr; i++)
    a[i] = 2.*mrand_r(seed)/MRAND_MAX-1.; // a[i] in [-1, 1]

  // For each instance
  ymin = DBL_MAX;
  ymax = -DBL_MAX;
  for(i = 0; i < ninst; i++) {
    // Compute the linear combination of the attributes as input
    y[i] = 0;
    for(j = 0; j < nattr; j++)
      y[i] += a[j]*X[i][j]; // y[i] in [-nattr, nattr]
    // Compute the minimal/maximal value of y[i]
    if(y[i] < ymin)
      ymin = y[i];
    if(ymax < y[i])
      ymax = y[i];
  }

  // Normalize the data in [0,1] then add noise
  for(i = 0; i < ninst; i++)
    y[i] = (y[i] - ymin) / (ymax - ymin) + noise*(2.*mrand_r(seed)/MRAND_MAX-1.);

  return y;
}

/**
 * Produce periodic regression data based on an input dataset
 */
void *pregression_data(const double *X, const size_t ninst, const size_t nattr, unsigned int *seed) {
  // Compute classical regression data
  double *y = regression_data(X, ninst, nattr, seed);

  // If allocation fails, return NULL
  if(y == NULL)
    return NULL;

  // Normalize all prediction such as to stands within [0,2 pi] +- 2*pi
  for(size_t i = 0; i < ninst; i++)
    y[i] = 2.*M_PI*(y[i]+(mrand_r(seed)%3-1));

  return y;
}

/**
 * Produce classification data based on an input dataset
 */
void *classification_data(const double *_X, const size_t ninst, const size_t nattr, unsigned int *seed) {
  const size_t nclasses = 3; // The number of class to produce (must be at least 2)
  const double misclassified = 0.1; // The ratio of misclassified instances
  const double (*X)[nattr] = (const double (*)[nattr]) _X;
  int *y = (int *) calloc(ninst, sizeof(int));
  double cclass[nclasses-1][nattr];
  size_t i, j, k;

  // If allocation fails, return NULL
  if(y == NULL)
    return NULL;

  // Initialize the center of each class
  for(i = 0; i < nclasses-1; i++)
    for(j = 0; j < nattr; j++)
      cclass[i][j] = 2.*mrand_r(seed)/MRAND_MAX-1.;

  // For each instance
  for(i = 0; i < ninst; i++) {

    // Seach the class of this instance
    y[i] = -1; // By default, class is unknown
    for(j = 0; j < nclasses-1 && y[i] == -1; j++) {
      // Check if this instance is contained within a box of length 2 centered around cclass
      for(k = 0; k < nattr && cclass[j][k]-1 <= X[i][k] && X[i][k] <= cclass[j][k]+1; k++)
        ;
      if(k == nattr)
        y[i] = j;
    }

    // If no class was found, set it to the last one by default
    if(y[i] == -1)
      y[i] = nclasses-1;

    // Add some misclassified instances
    if((double) mrand_r(seed)/MRAND_MAX < misclassified)
      y[i] = mrand_r(seed) % nclasses;
  }

  // Uncomment the line below in order to be sure that the only split is X[i][0] < 0.666
  // for(i = 0; i < NINSTANCES; i++) y[i] = (X[i][0] < 0.666);

  return y;
}

void regression_report(FILE *fout, const void *_y, const void *_ypred, const size_t ninst) {
  const double *y = (const double *) _y;
  const double *ypred = (const double *) _ypred;

  fprintf(fout, "\t\tCorrelation: %f%%\n", 100.*stat_corr(y,ypred,ninst));
  fprintf(fout, "\t\tMean absolute error: %f\n", stat_abserr(y,ypred,ninst));
}

void pregression_report(FILE *fout, const void *_y, const void *_ypred, const size_t ninst) {
  const double *y = (const double *) _y;
  const double *ypred = (const double *) _ypred;

  fprintf(fout, "\t\tPeriodic correlation: %f%%\n", 100.*circstat_corr(y,ypred,ninst));
  fprintf(fout, "\t\tMean absolute error: %f\n", circstat_abserr(y,ypred,ninst));
}

void classification_report(FILE *fout, const void *_y, const void *_ypred, const size_t ninst) {
  const int *y = (const int *) _y;
  const int *ypred = (const int *) _ypred;
  size_t ngood, nbad, ntot, nclasses;
  size_t i;
  int c;

  // Get the number of classes
  c = -1;
  for(i = 0; i < NINSTANCES; i++) {
    if(c < y[i])
      c = y[i];
    if(c < ypred[i])
      c = ypred[i];
  }
  nclasses = c + 1;


  // Compute the global number of correctly/incorrectly classified instances
  ngood = nbad = 0;
  for(i = 0; i < NINSTANCES; i++) {
    if(y[i] == ypred[i])
      ngood++;
  }
  nbad = ninst-ngood;

  fprintf(fout, "\t\tNumber of correctly classified instances: %zu/%zu (%.3f%%)\n", ngood, ninst, 100.*ngood/ninst);
  fprintf(fout, "\t\tNumber of incorrectly classified instances: %zu/%zu (%.3f%%)\n", nbad, ninst, 100.*nbad/ninst);

  // For each class
  for(c = 0; (size_t) c < nclasses; c++) {
    // Compute its True positive rate and false positive rate
    ngood = nbad = ntot = 0;
    for(i = 0; i < ninst; i++) {
      if(y[i] == c) {
        if(y[i] == ypred[i])
          ngood++;
        ntot++;
      } else if(ypred[i] == c)
        nbad++;
    }

    if(ntot > 0) {
      fprintf(fout, "\t\tClass %d (%zu instances):\n", c, ntot);
      fprintf(fout, "\t\t\tTrue positive rate: %.2f%%\n", 100.*ngood/ntot);
      fprintf(fout, "\t\t\tFalse positive rate: %.2f%%\n", 100.*nbad/(ninst-ntot));
    }
  }
}

int tree_predictor(threadpool_t *pool, void *LSypred, void *TSypred,
                                     const double *LS, const void *LSy,
                                     const double *TS,
                                     const size_t ninst, const size_t nattr,
                                     const size_t ysize,
                                     const tree_algorithm *algorithm) {
  tree_t *tree = NULL;
  struct timespec t0, t1;
  clock_t t0cpu, t1cpu;
  const char *tmpfilename = "tree_model.tmp";
  tree_t *tree_reload = NULL;
  FILE *f;
  size_t i;

  // Build the tree
  clock_gettime(CLOCK_MONOTONIC, &t0);t0cpu=clock();
  tree = tree_build(LS, LSy, ninst, nattr, pool, algorithm);

  if(tree == NULL) {
    fprintf(stderr, "Tree building failed");
    return -1;
  }
  clock_gettime(CLOCK_MONOTONIC, &t1);t1cpu=clock();
  fprintf(stdout, "Time taken to build the model: %fs (CPU time=%fs)\n", (double)(t1.tv_sec-t0.tv_sec) + 1e-9*(t1.tv_nsec-t0.tv_nsec), (double) (t1cpu-t0cpu)/CLOCKS_PER_SEC);

  // Compute the predictions on LS and on TS
  for(i = 0 ; i < ninst; i++) {
    memcpy(LSypred+i*ysize, tree_get(NULL,tree, LS+i*nattr), ysize);
    memcpy(TSypred+i*ysize, tree_get(NULL,tree, TS+i*nattr), ysize);
  }

  // Save the tree into a temporary file
  f = fopen(tmpfilename,"w");
  if(tree_dump(f, tree)) {
    fprintf(stderr, "Unable to dump tree");
    fclose(f);tree_free(tree);
    return -1;
  }
  fclose(f);

  // Load the tree
  f = fopen(tmpfilename,"r");
  tree_reload = tree_load(f);
  if(tree_reload == NULL) {
    fprintf(stderr, "Unable to load the tree");
    fclose(f);tree_free(tree);
    return -1;
  }
  fclose(f);

  // Delete the temporary file
  remove(tmpfilename);

  // Check that both the original tree and the reloaded tree gives the same results
  for(i = 0; i < ninst; i++) {
    size_t n1,n2;
    void *d1 = tree_get(&n1,tree, TS+i*nattr);
    void *d2 = tree_get(&n2,tree, TS+i*nattr);
    if(n1 != n2 || memcmp(d1,d2,n1)) {
      fprintf(stderr, "Reloaded tree don't give the same result as the original tree\n");
      tree_free(tree_reload);tree_free(tree);
      return -1;
    }
  }

  // Delete the reloaded tree
  tree_free(tree_reload);

  // Delete the tree
  tree_free(tree);

  return 0;
}

int regression_predictor(threadpool_t *pool, void *LSypred, void *TSypred, const double *LS, const void *LSy, const double *TS, const size_t ninst, const size_t nattr, void *arg) {
  return tree_predictor(pool, LSypred, TSypred, LS, LSy, TS, ninst, nattr,
                        sizeof(double), &regressiontree_algorithm);
}

int pregression_predictor(threadpool_t *pool, void *LSypred, void *TSypred, const double *LS, const void *LSy, const double *TS, const size_t ninst, const size_t nattr, void *arg) {
  return tree_predictor(pool, LSypred, TSypred, LS, LSy, TS, ninst, nattr,
                        sizeof(double), &pregressiontree_algorithm);
}

int classification_predictor(threadpool_t *pool, void *LSypred, void *TSypred, const double *LS, const void *LSy, const double *TS, const size_t ninst, const size_t nattr, void *arg) {
  return tree_predictor(pool, LSypred, TSypred, LS, LSy, TS, ninst, nattr,
                        sizeof(int), &classificationtree_algorithm);
}

int extratree(threadpool_t *pool, void *LSypred, void *TSypred, const double *LS, const void *LSy, const double *TS, const size_t ninst, const size_t nattr, void *arg) {
  tree_algorithm *algorithm;
  unsigned int k, nmin;
  int res = -1;

  // Build the extratree algorithm in regression mode
  extratree_default_parameters(&k,&nmin,nattr,1);
  algorithm = extratree_algorithm(k, nmin, 0, regressiontree_split_quality, NULL, dbl_compar_r, sizeof(double), NULL);
  if(algorithm == NULL)
    return res;

  // Replace the original getdata method with the one of the regression tree
  algorithm->getdata = regressiontree_algorithm.getdata;

  // Make the tree predictions with the provided algorithm
  res = tree_predictor(pool, LSypred, TSypred, LS, LSy, TS, ninst, nattr,
                       sizeof(double), algorithm);

  // Destroy the extratree algorithm
  extratree_algorithm_destroy(algorithm);

  return res;
}

int ExtraTrees_predictor(threadpool_t *pool, void *LSypred, void *TSypred, const double *LS, const void *LSy, const double *TS, const size_t ninst, const size_t nattr, void *arg) {
  ExtraTrees *ert;
  ExtraTreesPrediction *lspred, *tspred;
  const ExtraTreesMode mode = *((ExtraTreesMode *)arg);
  const unsigned int N = 250;
  const unsigned int k = 0;
  const unsigned int nmin = 0;
  const unsigned int seed = 0;
  struct timespec t0, t1;
  clock_t t0cpu, t1cpu;

  // Build the tree committee
  clock_gettime(CLOCK_MONOTONIC, &t0);t0cpu=clock();
  ert = ExtraTrees_build(LS,LSy,ninst,nattr,mode,N,k,nmin,seed,pool);
  if(ert == NULL)
    return -1;
  clock_gettime(CLOCK_MONOTONIC, &t1);t1cpu=clock();
  fprintf(stdout, "Time taken to build the model: %fs (CPU time=%fs)\n", (double)(t1.tv_sec-t0.tv_sec) + 1e-9*(t1.tv_nsec-t0.tv_nsec), (double) (t1cpu-t0cpu)/CLOCKS_PER_SEC);


  // Make the predictions on LS and TS
  for(size_t i = 0; i < ninst; i++) {
    lspred = ExtraTrees_get(ert, LS+i*nattr);
    tspred = ExtraTrees_get(ert, TS+i*nattr);
    
    switch(mode) {
      case ExtraTreesClassificationMode:
        ((int *) LSypred)[i] = lspred->classification_val;
        ((int *) TSypred)[i] = tspred->classification_val;
        break;
      case ExtraTreesRegressionMode:
        ((double *) LSypred)[i] = lspred->regression_val;
        ((double *) TSypred)[i] = tspred->regression_val;
        break;
      case ExtraTreesPeriodicRegressionMode:
        ((double *) LSypred)[i] = lspred->pregression_val;
        ((double *) TSypred)[i] = tspred->pregression_val;
        break;
    }
    free(lspred);free(tspred);
  }

  // Destroy the ExtraTrees
  ExtraTrees_free(ert);

  return 0;
}

void tree_test(threadpool_t *pool, dataproducer dp, predictor predic, reporter report,
               const double *LS, const double *TS, const size_t ninst,
               const size_t nattr, const size_t ysize, void *arg) {
  void *LSy, *TSy, *LSypred, *TSypred;
  unsigned int seed;

  // Get the predictions from LS and TS
  seed = 0; LSy = dp(LS, ninst, nattr, &seed);
  seed = 0; TSy = dp(TS, ninst, nattr, &seed);
  if(LSy == NULL || TSy == NULL) {
    fprintf(stderr, "Error while producing the LS and TS predictions\n");
    free(LSy);free(TSy);
    return;
  }

  // Allocate the result arrays
  LSypred = calloc(ninst, ysize);
  TSypred = calloc(ninst, ysize);
  if(LSypred == NULL || TSypred == NULL) {
    fprintf(stderr, "Error while allocating the resulting array\n");
    free(LSy);free(TSy);free(TSypred);free(LSypred);
    return;
  }

  // Get the prediction on LS and TS
  if(predic(pool,LSypred,TSypred,LS,LSy,TS,ninst,nattr, arg)) {
    fprintf(stderr, "Error while retrieving the LS and TS predictions\n");
    free(LSy);free(TSy);free(TSypred);free(LSypred);
    return;
  }
  fprintf(stdout, "\tLS Prediction\n");
  report(stdout, LSy, LSypred, ninst);
  fprintf(stdout, "\tTS Prediction\n");
  report(stdout, TSy, TSypred, ninst);
  free(LSy);free(TSy);free(TSypred);free(LSypred);
}

int main(const int argc, const char **argv) {
  threadpool_t *pool = NULL;
  double *LS, *TS;
  unsigned int seed = 0;

  // Allocate the thread pool
  pool = threadpool_create(NTHREAD);
  if(pool == NULL) {
    perror("Thread pool creation failed");
    return -1;
  }

  // Allocate the LS and TS array
  LS = (double *) calloc(NINSTANCES*NATTRIBUTES,sizeof(double));
  TS = (double *) calloc(NINSTANCES*NATTRIBUTES,sizeof(double));

  // Check allocations
  if(LS == NULL || TS == NULL) {
    perror("Dataset allocation fails");
    threadpool_destroy(pool);free(LS);free(TS);
    return -1;
  }

  // Fill datasets with random values
  srand(0);
  for(size_t i = 0; i < NINSTANCES*NATTRIBUTES; i++) {
    LS[i] = 2.*mrand_r(&seed)/MRAND_MAX-1.;
    TS[i] = 2.*mrand_r(&seed)/MRAND_MAX-1.;
  }

  fprintf(stdout, "CART classification tree test\n");
  fprintf(stdout, "=============================\n");
  tree_test(pool,classification_data,classification_predictor,classification_report,
            LS, TS, NINSTANCES, NATTRIBUTES,sizeof(int),NULL);
  fprintf(stdout, "CART regression tree test\n");
  fprintf(stdout, "=========================\n");
  tree_test(pool,regression_data,regression_predictor,regression_report,
            LS, TS, NINSTANCES, NATTRIBUTES,sizeof(double),NULL);

  fprintf(stdout, "CART regression tree test on periodic data\n");
  fprintf(stdout, "==========================================\n");
  tree_test(pool,pregression_data,regression_predictor,regression_report,
            LS, TS, NINSTANCES, NATTRIBUTES,sizeof(double),NULL);

  fprintf(stdout, "CART periodic regression tree test\n");
  fprintf(stdout, "==================================\n");
  tree_test(pool,pregression_data,pregression_predictor,pregression_report,
            LS, TS, NINSTANCES, NATTRIBUTES,sizeof(double),NULL);

  fprintf(stdout, "extratree regression test\n");
  fprintf(stdout, "=========================\n");
  tree_test(pool,regression_data,extratree,regression_report,
            LS, TS, NINSTANCES, NATTRIBUTES,sizeof(double),NULL);

  {
    ExtraTreesMode mode;

    fprintf(stdout, "ExtraTrees classification test\n");
    fprintf(stdout, "==============================\n");
    mode = ExtraTreesClassificationMode;
    tree_test(pool,classification_data,ExtraTrees_predictor,classification_report,
              LS, TS, NINSTANCES, NATTRIBUTES,sizeof(double),&mode);

    fprintf(stdout, "ExtraTrees regression test\n");
    fprintf(stdout, "==========================\n");
    mode = ExtraTreesRegressionMode;
    tree_test(pool,regression_data,ExtraTrees_predictor,regression_report,
              LS, TS, NINSTANCES, NATTRIBUTES,sizeof(double),&mode);

    fprintf(stdout, "ExtraTrees periodic regression test\n");
    fprintf(stdout, "===================================\n");
    mode = ExtraTreesPeriodicRegressionMode;
    tree_test(pool,pregression_data,ExtraTrees_predictor,pregression_report,
              LS, TS, NINSTANCES, NATTRIBUTES,sizeof(double),&mode);
  }

  threadpool_destroy(pool);free(LS);free(TS);

  sleep(1); // This line allow all memory ressources to be correctly freed (e.g. thread system ressources)

  return 0;
}
