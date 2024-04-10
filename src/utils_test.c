#include <utils.h>
#include <stdio.h>
#include <stdlib.h>

#define SIZE 16536

int _utils_test_compar(const void *a, const void *b) {
  return *((const int *) a) - *((const int *) b);
}

int _utils_test_compar_r(const void *a, const void *b, void *arg) {
  if(*((const int *) arg))
    return  *((const int *) a) - *((const int *) b);
  else
    return  *((const int *) b) - *((const int *) a);
}

int main(const int argc, const char **argv) {
  size_t idx[SIZE];
  int base[SIZE];
  int ascend;
  unsigned int seed = 0;
  size_t i;

  // Initialize the array to sort
  for(i = 0; i < SIZE; i++)
    base[i] = mrand_r(&seed);

  // Perform the argument sort
  argsort(idx, base, SIZE, sizeof(int), _utils_test_compar);

  // Check that the sort is correct
  for(i = 1; i < SIZE; i++) {
    if(base[idx[i-1]] > base[idx[i]])
      fprintf(stderr, "Error base[idx[%zu]] > base[idx[%zu]] (base[%zu]=%d, base[%zu]=%d)\n", i-1, i, idx[i-1], base[idx[i-1]], idx[i], base[idx[i]]);
  }

  // Perform the argument sort with 3 arguments (sort in ascending order)
  ascend = 1;
  argsort_r(idx, base, SIZE, sizeof(int), _utils_test_compar_r, &ascend);

  // Check that the sort is correct
  for(i = 1; i < SIZE; i++) {
    if(base[idx[i-1]] > base[idx[i]])
      fprintf(stderr, "Error base[idx[%zu]] > base[idx[%zu]] (base[%zu]=%d, base[%zu]=%d)\n", i-1, i, idx[i-1], base[idx[i-1]], idx[i], base[idx[i]]);
  }

  // Perform the argument sort with 3 arguments (sort in descending order)
  ascend = 0;
  argsort_r(idx, base, SIZE, sizeof(int), _utils_test_compar_r, &ascend);

  // Check that the sort is correct
  for(i = 1; i < SIZE; i++) {
    if(base[idx[i-1]] < base[idx[i]])
      fprintf(stderr, "Error base[idx[%zu]] < base[idx[%zu]] (base[%zu]=%d, base[%zu]=%d)\n", i-1, i, idx[i-1], base[idx[i-1]], idx[i], base[idx[i]]);
  }
}
