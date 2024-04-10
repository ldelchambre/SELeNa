#ifndef _TREECOMMITTEE_H_
#define _TREECOMMITTEE_H_

#include <tree.h>

typedef struct treecommittee_t treecommittee;

typedef void *(*treecommittee_combinefct)(const void **datas,
                                          const size_t *data_sizes,
                                          const size_t N,
                                          void *arg);

treecommittee *treecommittee_build(const double *X,
                                   const void *y,
                                   const size_t ninst,
                                   const size_t nattr,
                                   const size_t N,
                                   tree_algorithm * const * algorithms,
                                   threadpool_t *pool,
                                   treecommittee_combinefct combinefct,
                                   void *arg);

void *treecommittee_get(const treecommittee *tc, const double *x);

int treecommittee_dump(FILE *fout, const treecommittee *tc);

treecommittee *treecommittee_load(FILE *fin,
                                  treecommittee_combinefct combinefct,
                                  void *arg);

void treecommittee_free(treecommittee *tc);

#endif
