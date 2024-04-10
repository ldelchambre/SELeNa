#ifndef _TREE_H_
#define _TREE_H_

#include <threadpool.h>
#include <stddef.h>
#include <stdio.h>

/**
 * The version of the serialization protocol
 */
#define TREE_SERIALIZATION_VERSION 0x01

/**
 * Structure designed to contain the representation of a decision tree
 */
typedef struct tree tree_t;

/**
 * Structure containing the functions that are necessary in order to build a tree
 */
typedef struct tree_algorithm_t {
  /**
   * Initialize the algorithm according to the input dataset and associated classes (can be NULL).
   *
   * @param[in] X The input dataset (size = ninst x nattr), row-major order
   * @param[in] y The array of classes to predict (size = ninst)
   * @param[in] ninst The number of instances within X and y
   * @param[in] nattr The number of attributes within X
   * @param[in] arg User provided pointer
   *
   * @return 0 if we should process this dataset, something else otherwise
   */
  int (*init)(const double *X, const void *y, const size_t ninst, const size_t nattr, void *arg);
  /**
   * Given a subset of instances whose attributes are in X and classes in y,
   * this function have to provide the attribute index and associated value
   * according to which this subset should be splitted.
   *
   * @param[out] cutval The attribute value associated with this split
   * @param[in] X The input dataset
   * @param[in] y The array of classes to predict
   * @param[in] inst The indices of the instances composing this subset
   * @param[in] ninst The number of instances within this subset (i.e. the size of inst)
   * @param[in] nattr The number of attributes within the dataset
   * @param[in] arg User provided pointer
   *
   * @return The index of the attribute to use for splitting; nattr if we should create a leaf node; -1UL upon error
   */
  size_t (*getsplit)(double *cutval, const double *X, const void *y, const size_t *inst, const size_t ninst, const size_t nattr, void *arg);

  /**
   * Return a pointer to a newly allocated memory area containing the data that
   * should be associated with a leaf node (i.e. its predicted value).
   *
   * @param[out] data_size The size in bytes of the allocated memory area
   * @param[in] y The input set of classes
   * @param[in] inst The index of the instances composing this subset
   * @param[in] ninst The number of instances within this subset (i.e. the size of inst)
   * @param[in] arg User provided pointer
   *
   * Note that the tree data must be contiguous in memory (i.e. it should not rely
   * on pointers) and will be freed with freedata or with free if freedata is NULL.
   *
   * @return The newly allocated data associated with this subset; NULL upon error
   */
  void *(*getdata)(size_t *data_size, const void *y, const size_t *inst,  const size_t ninst, void *arg);

  /**
   * A user provided pointer that will be passed to the init, getsplit, getdata, serialize, unserialize and freedata functions
   */
  void *arg;
} tree_algorithm;

/**
 * Build a decision tree for the given dataset and classes using the provided
 * algorithm.
 *
 * @param[in] X The input dataset
 * @param[in] y The array of classes to predict
 * @param[in] ninst The number of instances within X and y
 * @param[in] nattr The number of attributes within X
 * @param[in] pool The thread pool to use in order to build the tree
 * @param[in] algorithm The algorithm that will be used in order to build the tree
 *
 * @return tree upon success; NULL upon failure
 */
tree_t *tree_build(const double *X, const void *y, const size_t ninst, const size_t nattr, threadpool_t *pool, const tree_algorithm *algorithm);

/**
 * Structure containing the parameters that are necessary in order for a thread
 * pool to develop a tree based on an input dataset
 */
typedef struct tree_build_job_t tree_build_job;

/**
 * Build a decision tree for the given dataset and classes using the provided
 * algorithm without waiting for the threadpool jobs to end.
 *
 * @param[in] X The input dataset
 * @param[in] y The array of classes to predict
 * @param[in] ninst The number of instances within X and y
 * @param[in] nattr The number of attributes within X
 * @param[in] pool The thread pool to use in order to build the tree
 * @param[in] algorithm The algorithm that will be used in order to build the tree
 *
 * @return A pointer to the tree building job that should be passed to the tree_build_end function, or NULL upon error
 */
tree_build_job *tree_build_nowait(const double *X, const void *y, const size_t ninst, const size_t nattr, threadpool_t *pool, const tree_algorithm *algorithm);

/**
 * Return the current build progress.
 *
 * @param[in] The job returned by the tree_buid_nowait function
 *
 * @return The current progress of building [0, 1]
 */
int tree_build_progess(size_t *n, size_t *ntot, size_t *nnode, tree_build_job *job);

/**
 * Get the tree developed during a job execution.
 *
 * @param[in,out] The job returned by the tree_buid_nowait function.
 *
 *  Note that the user should ensure that the pool used within tree_build_nowait
 * function finished all the jobs related to the tree building before calling this
 * function.
 *
 * @return The developed tree upon success, NULL otherwise
 */
tree_t *tree_build_end(tree_build_job *job);

/**
 * Get the predicted value associatedd with the instance stored in x
 *
 * @param[out] data_size If not NULL, set to the size of the returned data [bytes]
 * @param[in] tree The decision tree to search in
 * @param[in] x The attributes associated of the given instance (must have the same length as the number of attributes used in order to build the tree)
 *
 * @return A pointer to the predicted value for this instance or NULL upon failure
 */
void *tree_get(size_t *data_size, const tree_t *tree, const double *x);

/**
 * Dump tree into file.
 *
 * @param[in,out] fout The file where the tree should be written
 * @param[in] tree The tree to write
 *
 * @return 0 upon success
 */
int tree_dump(FILE *fout, const tree_t *tree);

/**
 * Load tree from file.
 *
 * @param[in] fin The file from which to read the tree
 *
 * @return The loaded tree upon success, NULL otherwise
 */
tree_t *tree_load(FILE *fin);

/**
 * Destroy a decision tree.
 *
 * @param[in,out] tree Pointer to the decision tree to destroy
 */
void tree_free(tree_t *tree);

#endif
