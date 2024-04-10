#ifndef _TREEALGORITHM_H_
#define _TREEALGORITHM_H_

#include <tree.h>

/******************************************************************************
 * GENERAL PURPOSE UTILITY FUNCTIONS                                          *
 ******************************************************************************/

/**
 * Prototype of the split quality functions.
 * Functions of this kind should provide the score of a split regarding its classes.
 *
 * The given subset is split into 2 subsets, S1 and S2 in the following way
 *   S1 = { y[inst[i]] | X[inst[i]][attr] < vattr }; i = 0 ... ninst
 *   S2 = { y[inst[i]] | X[inst[i]][attr] >= vattr }; i = 0 ... ninst
 *
 * @param[in] S1 The input set of classes for the first subset
 * @param[in] n1 The size of the first subset
 * @param[in] S2 The input set of classes for the second subset
 * @param[in] n2 The size of the second subset
 * @param[in] arg A user provided pointer
 *
 * @return The quality of this split (the higher is the better) or -DBL_MAX upon error
 */
typedef double (*split_quality_fct)(const void *, const size_t, const void *, const size_t, void *);

/**
 * Compute the best split over all attributes and attribute values based on the
 * provided split quality function.
 *
 * @param[out] cutval The attribute value associated with this split
 * @param[in] X The input set of attributes
 * @param[in] y The input set of classes
 * @param[in] inst The indices of the instances composing this subset
 * @param[in] ninst The number of instances within this subset (i.e. the size of inst)
 * @param[in] nattr The number of attributes
 * @param[in] ysize The size of the elements contained within y (e.g. sizeof(double) if y contains doubles)
 * @param[in] split_quality The split quality function to use
 * @param[in] ycmp A function designed to compare the elements within y (must return -1,0,1 respectively if the first element is lower, equal or greater than the second one)
 * @param[in] ycmparg A user provided pointer to pass to the ycmp function
 * @param[in] arg User provided pointer to pass to the split quality function
 *
 * @return The index of the attribute to use for splitting; nattr if we should create a leaf node or -1UL upon error
 */
size_t tree_best_split(double *cutval, const double *X, const void *y, const size_t *inst,
                                       const size_t ninst, const size_t nattr, const size_t ysize,
                                       split_quality_fct split_quality, int (*ycmp)(const void *, const void *, void *),
                                       void *ycmparg, void *arg);


/**
 * Function flattening the given y prediction.
 *
 * The flattening is performed as follow for the given result array, res:
 * res[i] = y[inst[i]] for i=0...ninst
 *
 * @param[in] y The input set of classes
 * @param[in] inst The index of the instances composing this subset
 * @param[in] ninst The number of instances within this subset (i.e. the size of inst)
 * @param[in] ysize The size of the elements contained within y
 *
 * @return A pointer to the flattened y prediction or NULL upon error (MUST be freed with free)
 */
void *tree_flatten_y(const void *y, const size_t *inst, const size_t ninst, const size_t ysize);

/**
 * Find the elligible attributes for a given subset of the dataset.
 * Elligible attributes are those for which the class to predict is not constant
 * and for which the set of attribute's values is not constant.
 *
 * @param[out] eattr The set of attributes that are elligible (must be at least of size nattr)
 * @param[in] X The input dataset
 * @param[in] y The input classes to predict
 * @param[in] inst The indices of the instances composing this subset
 * @param[in] ninst The number of instances within this subset (i.e. the size of inst)
 * @param[in] nattr The number of attributes within X
 * @param[in] ysize The size of the elements contained within y (e.g. sizeof(double) if y contains doubles)
 * @param[in] ycmp A function designed to compare the elements within y (must return -1,0,1 respectively if the first element is lower, equal or greater than the second one)
 * @param[in] ycmparg A user provided pointer to pass to the ycmp function
 *
 * @return The number of elligible attributes
 */
size_t tree_elligible_attributes(size_t *eattr, const double *X, const void *y, const size_t *inst,
                                                const size_t ninst, const size_t nattr, const size_t ysize,
                                                int (*ycmp)(const void *, const void *, void *), void *ycmparg);

/******************************************************************************
 * REGRESSION TREE ALGORITHM AND UTILITY FUNCTIONS                            *
 ******************************************************************************/
/**
 * Regression tree split quality function.
 *
 * @param[in] S1 The input set of classes for the first subset
 * @param[in] n1 The size of the first subset
 * @param[in] S2 The input set of classes for the second subset
 * @param[in] n2 The size of the second subset
 * @param[in] arg A user provided pointer
 *
 * @return The variance reduction minus the variance of the subset.
 */
double regressiontree_split_quality(const void *S1, const size_t n1, const void *S2, const size_t n2, void *arg);

/**
 * CART regression tree algorithm
 *
 * Used in the cases where y are of type double within the tree_build function
 *
 * Data returned by the tree_get function are arrays of 2 doubles:
 *  - tree_get[0] : The predicted valued (mean value amongst the leaf nodes)
 *  - tree_get[1] : The variance around the mean value (standard deviation amongst the leaf nodes)
 */
extern const tree_algorithm regressiontree_algorithm;

/******************************************************************************
 * PERIODIC REGRESSION TREE ALGORITHM AND UTILITY FUNCTIONS                   *
 ******************************************************************************/
/**
 * Periodic regression tree split quality function.
 *
 * @param[in] S1 The input set of classes for the first subset
 * @param[in] n1 The size of the first subset
 * @param[in] S2 The input set of classes for the second subset
 * @param[in] n2 The size of the second subset
 * @param[in] arg A user provided pointer
 *
 * @return The cyclic variance reduction minus the cyclic variance of the subset.
 */
double pregressiontree_split_quality(const void *S1, const size_t n1, const void *S2, const size_t n2, void *arg);

/**
 * The default sigma confidence interval to use for the periodic tree prediction
 * uncertainties. By default, the latter is set to 2-sigma = 95.45% confidence.
 */
extern double pregression_sigmaconf;

/**
 * CART periodic regression tree algorithm.
 *
 * Used in the cases where y are of type double within the tree_build function
 * and that y[i] ± k*2*pi are equivalent solutions.
 *
 * Data returned by the tree_get function are arrays of 2 doubles:
 *  - tree_get[0] : The predicted valued (cyclic mean value amongst the leaf nodes)
 *  - tree_get[1] : The 95% uncertainty interval of the mean value (can be changed by setting the pregression_sigmaconf to a value different that 2.)
 *
 * arg is set by default to a pointer to pregression_sigmaconf but it
 * can be set to a user-defined double variable or can be NULL in which case
 * a confidence interval of 2 sigma will be assumed.
 */
extern const tree_algorithm pregressiontree_algorithm;

/******************************************************************************
 * CLASSIFICATION TREE ALGORITHM AND UTILITY FUNCTIONS                        *
 ******************************************************************************/
 /**
  * The data type produced by the classificationtree_algorithm;
  */
typedef struct classificationtree_data_t {
  int c; // The class at the leaf node
  double pc; // The probability of having this class at the leaf node
} classificationtree_data;
/**
 * Constant identifying the information gain as measure of impurity.
 */
#define CLASSIFICATION_MEASURE_IG 0

/**
 * Constant identifying the Gini impurity measure
 */
#define CLASSIFICATION_MEASURE_GINI 1

/**
 * The default impurity measure to use within the CART classification tree development.
 * By default, set to CLASSIFICATION_MEASURE_IG.
 */
extern int classificationtree_measure;

/**
 * Classification tree split quality function.
 *
 * @param[in] S1 The input set of classes for the first subset
 * @param[in] n1 The size of the first subset
 * @param[in] S2 The input set of classes for the second subset
 * @param[in] n2 The size of the second subset
 * @param[in] arg Pointer to a const int that select the impurity function to use
 *
 * If arg point to a value containing:
 *  - CLASSIFICATION_MEASURE_IG: This function return the information gain minus the entropy of the subsets
 *  - CLASSIFICATION_MEASURE_GINI : This function return minus the weighted sum of the Gini impurity of the subsets
 * Note that if arg is NULL, then CLASSIFICATION_MEASURE_IG is assumed.
 *
 * @return The score of the provided split (see the arg parameter for more details)
 */
double classificationtree_split_quality(const void *S1, const size_t n1, const void *S2, const size_t n2, void *arg);

/**
 * CART classification tree algorithm.
 *
 * Used in the cases where y are of type int ** within the tree_build function
 *
 * Data returned by the tree_get function are of type classificationtree_data.
 *
 * arg value is set by default to a pointer to classificationtree_measure but it
 * can be set to a user-defined int variable or can be NULL in which case
 * CLASSIFICATION_MEASURE_IG is assumed as an impurity measure function.
 */
extern const tree_algorithm classificationtree_algorithm;

#endif
