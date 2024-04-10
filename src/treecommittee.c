#include <treecommittee.h>

#include <stdlib.h>

struct treecommittee_t {
  size_t N;
  tree_t **trees;
  treecommittee_combinefct combinefct;
  void *combinefctarg;
};

treecommittee *treecommittee_build(const double *X,
                                   const void *y,
                                   const size_t ninst,
                                   const size_t nattr,
                                   const size_t N,
                                   tree_algorithm * const * algorithms,
                                   threadpool_t *pool,
                                   treecommittee_combinefct combinefct,
                                   void *arg) {
  treecommittee *rc;

  // We can not create a RandomCommittee with 0 voters
  if(N == 0) {
    fprintf(stderr, "Can not create an empty treecommittee\n");
    return NULL;
  }

  // Allocate the Random committee
  rc = (treecommittee *) malloc(sizeof(treecommittee));
  if(rc == NULL)
    return NULL;

  // Initialize the Random Committee
  rc->N = N;
  rc->trees = (tree_t **) calloc(N, sizeof(tree_t *));
  rc->combinefct = combinefct;
  rc->combinefctarg = arg;
  if(rc->trees == NULL) {
    free(rc);
    return NULL;
  }

  // Build the N trees, one after another so as to spare memory consumption
  for(size_t i = 0; i < N; i++) {
    rc->trees[i] = tree_build(X,y,ninst,nattr,pool,algorithms[i]);
    fprintf(stderr, "Tree %zd built!\n", i);
  }

  // Check that all trees were correctly built
  for(size_t i = 0; i < N; i++) {
    if(rc->trees[i] == NULL) {
      treecommittee_free(rc);
      return NULL;
    }
  }

  return rc;
}

void *treecommittee_get(const treecommittee *tc, const double *x) {
  const void *datas[tc->N];
  size_t data_sizes[tc->N];

  // Get the predictions from all trees
  for(size_t i = 0; i < tc->N; i++)
    datas[i] = tree_get(&data_sizes[i], tc->trees[i], x);

  // Combine the returned predictions
  return tc->combinefct(datas, data_sizes, tc->N, tc->combinefctarg);
}

typedef struct _treecommittee_get_thread_arg_t {
  tree_t *tree;
  const double *x;
  const void **data;
  size_t *data_size;
} _treecommittee_get_thread_arg;

void _treecommittee_get2(void *arg) {
  _treecommittee_get_thread_arg *thread_arg = (_treecommittee_get_thread_arg *) arg;
  
  *(thread_arg->data) = tree_get(thread_arg->data_size, thread_arg->tree, thread_arg->x);
}

int treecommittee_dump(FILE *fout, const treecommittee *rc) {
  // Write the number of tree within the file
  if(fwrite(&rc->N, sizeof(size_t), 1, fout) == 0) {
    fprintf(stderr, "Unable to write the RandomCommittee header\n");
    return -1;
  }

  // Write all trees
  for(size_t i = 0; i < rc->N; i++) {
    if(tree_dump(fout,rc->trees[i])) {
      fprintf(stderr, "Unable to write the %zuth tree\n", i);
      return -1;
    }
  }

  return 0;
}

treecommittee *treecommittee_load(FILE *fin, treecommittee_combinefct combinefct, void *arg) {
  treecommittee *rc;

  // Allocate the Random Committee
  rc = (treecommittee *) malloc(sizeof(treecommittee));
  if(rc == NULL)
    return NULL;

  // Read the number of tree present within the file
  if(fread(&rc->N,sizeof(size_t), 1, fin) == 0) {
    fprintf(stderr, "Unable to read the number of trees within the file\n");
    free(rc);
    return NULL;
  }

  // Initialize the Random Committee
  rc->trees = (tree_t **) calloc(rc->N, sizeof(tree_t *));
  rc->combinefct = combinefct;
  rc->combinefctarg = arg;
  if(rc->trees == NULL) {
    fprintf(stderr, "Unable to allocate trees in the Random committee\n");
    free(rc);
    return NULL;
  }

  // Read all the N trees
  for(size_t i = 0; i < rc->N; i++) {
    rc->trees[i] = tree_load(fin);
    if(rc->trees[i] == NULL) {
      fprintf(stderr, "Unable to load the %zuth tree from the Random committee\n", i);
      for(size_t j = 0; j < i; j++) tree_free(rc->trees[j]);
      free(rc);free(rc->trees);
      return NULL;
    }
  }

  return rc;
}

void treecommittee_free(treecommittee *rc) {
  if(rc) {
    for(size_t i = 0; i < rc->N; i++)
      tree_free(rc->trees[i]);
    free(rc->trees);
    free(rc);
  }
}
