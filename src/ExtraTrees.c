#include <ExtraTrees.h>

#include <stdlib.h>

#include <circstat.h>
#include <information.h>
#include <stat.h>

#include <extratree.h>
#include <treealgorithm.h>
#include <treecommittee.h>

#include <utils.h>

#include <math.h>

/******************************************************************************
 * The ExtraTrees structure                                                   *
 ******************************************************************************/
struct ExtraTrees_t {
  ExtraTreesMode mode;
  treecommittee *tcommittee;
};

/******************************************************************************
 * ExtraTrees prediction management functions                                 *
 ******************************************************************************/
void *_ExtraTreesClassification_getdata(size_t *data_size, const void *_y, const size_t *inst,  const size_t ninst, void *arg) {
  information_event_count *data;
  int *y;
  size_t n;

  // Flatten the y array
  y = (int *) tree_flatten_y(_y, inst, ninst, sizeof(int));
  if(y == NULL)
    return NULL;

  // Get the count number of each class
  data = information_count(&n, y, ninst);
  free(y);
  if(data == NULL)
    return NULL;

  *data_size = n*sizeof(information_event_count);
  return data;
}

int _ExtraTreesClassification_combine_sort(const void *a, const void *b) {
  return ((const classificationtree_data *) a)->c-((const classificationtree_data *) b)->c;
}

void *_ExtraTreesClassification_combine(const void **datas, const size_t *data_sizes, const size_t N, void *arg) {
  const information_event_count **counts = (const information_event_count **) datas;
  ExtraTreesPrediction *pred = (ExtraTreesPrediction *) malloc(sizeof(ExtraTreesPrediction));
  classificationtree_data *preds;
  size_t n;
  
  // If pred is NULL, then exit
  if(!pred) return NULL;

  // Build the array of combined predictions
  {
    size_t ninst[N];
    size_t ipred = 0;

    // Compute the number of instances within each prediction and the number of
    // (not necesseraly distinct) classes
    n = 0;
    for(size_t i = 0; i < N; i++) {
      const size_t ni = data_sizes[i]/sizeof(information_event_count);
      ninst[i] = 0;
      for(size_t j = 0; j < ni; j++)
        ninst[i] += counts[i][j].count;
      n += ni;
    }

    // Allocate the array of predictions
    preds = (classificationtree_data *) calloc(n, sizeof(classificationtree_data));
    if(preds == NULL) {
      free(pred);
      return NULL;
    }

    // Fill the array of predictions
    for(size_t i = 0; i < N; i++) {
      const size_t ni = data_sizes[i]/sizeof(information_event_count);
      for(size_t j = 0; j < ni; j++) {
        preds[ipred].c = counts[i][j].event;
        preds[ipred].pc = 1./N*counts[i][j].count/ninst[i];
        ipred++;
      }
    }
  }

  // Sort the array of predictions
  // XXX counts arrays are sorted according to their event, it may hence be possible to have a quicker sort
  qsort(preds, n, sizeof(classificationtree_data), _ExtraTreesClassification_combine_sort);

  // Combine predictions
  {
    double p, pbest = -1.;
    size_t i = 0;

    while(i < n) {
      const int c = preds[i].c;
      p = 0;
      while(i < n && preds[i].c == c) {
        p += preds[i].pc;
        i++;
      }
      if(pbest < p) {
        pbest = p;
        pred->classification_val = c;
        pred->uncertainty = 1.-p;
      }
    }
  }

  // Due to round-off errors, uncertainty may get negative (i.e. -1e-16)
  if(pred->uncertainty < 0.)
    pred->uncertainty = 0.;

  // Delete the predictions array
  free(preds);

  return pred;
}

void *_ExtraTreesRegression_getdata(size_t *data_size, const void *_y, const size_t *inst,  const size_t ninst, void *arg) {
  double *data = (double *) malloc(sizeof(double));
  double *y;

  // Check that data was correctly allocated
  if(data == NULL)
    return NULL;

  // Flatten the y array
  y = (double *) tree_flatten_y(_y, inst, ninst, sizeof(double));
  if(y == NULL) {
    free(data);
    return NULL;
  }

  // Get the data associated with this leaf
  *data = stat_mean(y, ninst);

  // Delete the y array
  free(y);

  *data_size = sizeof(double);
  return data;
}

void *_ExtraTreesRegression_combine(const void **datas, const size_t *data_sizes, const size_t N, void *arg) {
  ExtraTreesPrediction *pred = (ExtraTreesPrediction *) malloc(sizeof(ExtraTreesPrediction));
  double y[N];
  double var;

  for(size_t i = 0; i < N; i++)
    y[i] = ((const double **) datas)[i][0];

  pred->regression_val = stat_mean2(&var, y, N);
  pred->uncertainty = sqrt(var);

  return pred;
}

void *_ExtraTreesPeriodicRegression_getdata(size_t *data_size, const void *_y, const size_t *inst,  const size_t ninst, void *arg) {
  double *data = (double *) malloc(sizeof(double));
  double *y;

  // Check that data was correctly allocated
  if(data == NULL)
    return NULL;

  // Flatten the y array
  y = (double *) tree_flatten_y(_y, inst, ninst, sizeof(double));
  if(y == NULL) {
    free(data);
    return NULL;
  }

  // Get the data associated with this leaf
  *data = circstat_mean(y,ninst);

  // Delete the y array
  free(y);

  *data_size = sizeof(double);
  return data;
}

void *_ExtraTreesPeriodicRegression_combine(const void **datas, const size_t *data_sizes, const size_t N, void *arg) {
  ExtraTreesPrediction *pred = (ExtraTreesPrediction *) malloc(sizeof(ExtraTreesPrediction));
  double y[N];
  double var = 0;

  // Compute the circular mean of the provided angles
  for(size_t i = 0; i < N; i++)
    y[i] = ((const double **) datas)[i][0];
  pred->pregression_val = circstat_mean(y,N);

  // Compute the standard deviation of the values around the mean.
  // Note: We compute the standard deviation as if the samples were drawn from
  //       a Gaussian distribution while they comes from a Von Mises distribution.
  //       This is justified by the fact that the Von Mises distribution tends to
  //       a normal distribution when its kappa parameter (i.e. the inverse
  //       spread of the distribution) tends to infinity.
  //       We can reasonably suppose that each tree will give quite similar
  //       results (i.e. large kappa; small standard deviation) such that the
  //       sample distribution can be considered as Gaussian-profiled.
  for(size_t i = 0; i < N; i++) {
    const double d = circstat_diff(y[i],pred->pregression_val);
    var += d*d;
  }
  pred->uncertainty = sqrt(var/N);

  return pred;
}

/******************************************************************************
 * ExtraTrees methods                                                         *
 ******************************************************************************/
ExtraTrees *ExtraTrees_build(const double *X,
                             const void *y,
                             const size_t ninst,
                             const size_t nattr,
                             const ExtraTreesMode mode,
                             const size_t N,
                             const unsigned int k,
                             const unsigned int nmin,
                             const unsigned int seed,
                             threadpool_t *pool) {
  ExtraTrees *ert;

  // We can not create a RandomCommittee with 0 voters
  if(N == 0) {
    fprintf(stderr, "Can not create an empty ExtraTrees\n");
    return NULL;
  }

  // Allocate the ExtraTrees
  ert = (ExtraTrees *) malloc(sizeof(ExtraTrees));
  if(ert == NULL)
    return NULL;

  // Initialize the ExtraTrees
  ert->mode = mode;

  // Build the ExtraTrees
  {
    tree_algorithm *algorithms[N];
    unsigned int kp, nminp;

    // Initialize the k and nmin parameters according to the current problem size
    extratree_default_parameters(&kp, &nminp, nattr, mode != ExtraTreesClassificationMode);
    if(k > 0) kp = k;
    if(nmin > 0) nminp = nmin;

    // Initialize the N extratree algorithms
    for(size_t i = 0; i < N; i++) {
      switch(mode) {
        case ExtraTreesClassificationMode:
          algorithms[i] = extratree_algorithm(kp, nminp, seed+i,
                                              classificationtree_split_quality,
                                              NULL, int_compar_r, sizeof(int), NULL);
          if(algorithms[i]) algorithms[i]->getdata = _ExtraTreesClassification_getdata;
          break;
        case ExtraTreesRegressionMode:
          algorithms[i] = extratree_algorithm(kp, nminp, seed+i,
                                            regressiontree_split_quality,
                                            NULL,dbl_compar_r, sizeof(double),NULL);
          if(algorithms[i]) algorithms[i]->getdata = _ExtraTreesRegression_getdata;
          break;
        case ExtraTreesPeriodicRegressionMode:
          algorithms[i] = extratree_algorithm(kp, nminp, seed+i,
                                            pregressiontree_split_quality,
                                            NULL,dbl_compar_r, sizeof(double),NULL);
          if(algorithms[i]) algorithms[i]->getdata = _ExtraTreesPeriodicRegression_getdata;
          break;
      }
      // If the last agorithm allocation fails, delete all previously allocated algorithms, then exit
      if(algorithms[i] == NULL) {
        for(size_t j = 0; j < i; j++) extratree_algorithm_destroy(algorithms[j]);
        free(ert);
        return NULL;
      }
    }

    // Build the tree committee
    switch(ert->mode) {
      case ExtraTreesClassificationMode:
        ert->tcommittee = treecommittee_build(X,y,ninst,nattr,N,algorithms,pool,_ExtraTreesClassification_combine, ert);
        break;
      case ExtraTreesRegressionMode:
        ert->tcommittee = treecommittee_build(X,y,ninst,nattr,N,algorithms,pool,_ExtraTreesRegression_combine, ert);
        break;
      case ExtraTreesPeriodicRegressionMode:
        ert->tcommittee = treecommittee_build(X,y,ninst,nattr,N,algorithms,pool,_ExtraTreesPeriodicRegression_combine, ert);
        break;
    }

    // If the tree committee Building fails
    if(ert->tcommittee == NULL) {
      free(ert);
      for(size_t i = 0; i < N; i++) extratree_algorithm_destroy(algorithms[i]);
      return NULL;
    }

    // Destroy the extratree algorithms
    for(size_t i = 0; i < N; i++)
      extratree_algorithm_destroy(algorithms[i]);
  }

  return ert;
}

ExtraTreesPrediction *ExtraTrees_get(const ExtraTrees *ert, const double *x) {
  return (ExtraTreesPrediction *) treecommittee_get(ert->tcommittee, x);
}

ExtraTreesMode ExtraTrees_mode(const ExtraTrees *ert) {
  return ert->mode;
}

int ExtraTrees_dump(FILE *fout, const ExtraTrees *ert) {
  if(fwrite(&ert->mode, sizeof(ExtraTreesMode), 1, fout) == 0) return -1;
  return treecommittee_dump(fout, ert->tcommittee);
}

ExtraTrees *ExtraTrees_load(FILE *fin) {
  ExtraTrees *ert = (ExtraTrees *) malloc(sizeof(ExtraTrees));

  if(ert == NULL) return NULL;

  // Read the classification mode of the ExtraTrees
  if(fread(&ert->mode,sizeof(ExtraTreesMode), 1, fin) == 0) {
    fprintf(stderr, "Unable to read the classification scheme used within the ExtraTree building\n");
    free(ert);
    return NULL;
  }

  // Read the tree committee from file
  switch(ert->mode) {
    case ExtraTreesClassificationMode:
      ert->tcommittee = treecommittee_load(fin, _ExtraTreesClassification_combine, ert);
      break;
    case ExtraTreesRegressionMode:
      ert->tcommittee = treecommittee_load(fin, _ExtraTreesRegression_combine, ert);
      break;
    case ExtraTreesPeriodicRegressionMode:
      ert->tcommittee = treecommittee_load(fin, _ExtraTreesPeriodicRegression_combine, ert);
      break;
  }

  // If the tree committee loading fails
  if(ert->tcommittee == NULL) {
    free(ert);
    return NULL;
  }

  return ert;
}

void ExtraTrees_free(ExtraTrees *ert) {
  if(ert) {
    treecommittee_free(ert->tcommittee);
    free(ert);
  }
}
