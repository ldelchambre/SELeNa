#ifndef _EXTRATREES_WRAPPER_H_
#define _EXTRATREES_WRAPPER_H_

#include <tree.h>

typedef struct ExtraTrees_t ExtraTrees;

typedef enum ExtraTreesMode_t {
  ExtraTreesClassificationMode = 0,
  ExtraTreesRegressionMode = 1,
  ExtraTreesPeriodicRegressionMode = 2
} ExtraTreesMode;

typedef struct ExtraTreesPrediction_t {
  union {
    int classification_val;
    double regression_val;
    double pregression_val;
  };
  double uncertainty;
} ExtraTreesPrediction;

ExtraTrees *ExtraTrees_build(const double *X,
                             const void *y,
                             const size_t ninst,
                             const size_t nattr,
                             const ExtraTreesMode mode,
                             const size_t N,
                             const unsigned int k,
                             const unsigned int nmin,
                             const unsigned int seed,
                             threadpool_t *pool);

ExtraTreesPrediction *ExtraTrees_get(const ExtraTrees *ert, const double *x);

ExtraTreesMode ExtraTrees_mode(const ExtraTrees *ert);

int ExtraTrees_dump(FILE *fout, const ExtraTrees *ert);

ExtraTrees *ExtraTrees_load(FILE *fin);

void ExtraTrees_free(ExtraTrees *ert);

#endif
