#ifndef _UTILS_H_
#define _UTILS_H_

#include <stddef.h>

# define MRAND_MAX 0x7FFFFFFFU

void *memdup(const void *src, const size_t n);
int int_compar(const void *a, const void *b);
int int_compar_r(const void *a, const void *b, void *arg);
int dbl_compar(const void *a, const void *b);
int dbl_compar_r(const void *a, const void *b, void *arg);
const void *min(const void *x, const size_t nmemb, const size_t size, int (*compar)(const void *, const void *));
const void *max(const void *x, const size_t nmemb, const size_t size, int (*compar)(const void *, const void *));
const void *min_r(const void *x, const size_t nmemb, const size_t size, int (*compar)(const void *, const void *, void *), void *arg);
const void *max_r(const void *x, const size_t nmemb, const size_t size, int (*compar)(const void *, const void *, void *), void *arg);
void argsort(size_t *idx, const void *base, const size_t nmemb, const size_t size, int (*compar)(const void *, const void *));
void argsort_r(size_t *idx, const void *base, const size_t nmemb, const size_t size, int (*compar)(const void *, const void *, void *), void *arg);
int mrand_r(unsigned int *seed); // This random number generator is specified is POSIX.1-2001 and allow to have the same sequence of random numbers over all machines
double mrandn(unsigned int *seed);
unsigned int hashcode(const void *data, const size_t n);
void swap(void *a, void *b, size_t size);
void shuffle(void *base, size_t nmemb, size_t size, unsigned int *seed);
int parse_double(double *d, const char *s, char **endptr, int last);

#endif
