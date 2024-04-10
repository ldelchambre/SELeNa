#ifndef _EXTRATREE_H_
#define _EXTRATREE_H_

#include <tree.h>
#include <treealgorithm.h>

void extratree_default_parameters(unsigned int *k, unsigned int *nmin, const size_t nattr, const int regression);

/**
 * @param arg The argument to pass to the split quality function
 */
tree_algorithm *extratree_algorithm(const unsigned int k,
                                    const unsigned int nmin,
                                    const unsigned int seed,
                                    split_quality_fct split_quality,
                                    void *split_quality_arg,
                                    int (*ycmp)(const void *, const void *, void *),
                                    size_t ysize,
                                    void *ycmp_arg);

void extratree_algorithm_destroy(tree_algorithm *algorithm);

#endif
