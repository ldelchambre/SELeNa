#ifndef _COMBINATION_H_
#define _COMBINATION_H_

/**
 * Return the unique combination of k numbers taken in [0, n-1]
 * 
 * Arrays are returned through the callback function fct taking as argument
 * void fct(const unsigned int *idx, const unsigned int n, void *arg)
 * 
 * Example: combine(4,2) = { [0,1], [0,2], [0, 3], [1, 2], [1, 3], [2, 3]}
 * 
 * Return 0 upon success, -1 if the function failed to allocated the array of indices
 */
int combine(const unsigned int n, const unsigned int k, void (*fct)(const unsigned int *, const unsigned int, void *), void *arg);

#endif
