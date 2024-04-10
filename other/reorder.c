#include <stdio.h>
#include <stdlib.h>

#define ARRAY_SIZE 32

void reorder(double *a, const unsigned int *idx, const unsigned int n) {
  unsigned int i, j;
  double tmp;

  fprintf(stdout, "Initial arrays\n");
  for(i = 0; i < n; i++)
    fprintf(stdout, "\t%u %u %f\n", i, idx[i], a[i]);

  for(i = 0; i < n; i++) {
    fprintf(stdout, "a[%u] <-  a[%u]\n", i, idx[i]);
    for(j = idx[i]; j < i; j = idx[j])
      fprintf(stdout, "a[%u]  -> a[%u]\n", j, idx[j]);
    if(i != j) {
      fprintf(stdout, "a[%u] <-> a[%u]\n", i, j);
      tmp = a[i];
      a[i] = a[j];
      a[j] = tmp;
    }
  }

  fprintf(stdout, "Final array\n");
  for(i = 0; i < n; i++)
    fprintf(stdout, "\t%u %f\n", i, a[i]);
}

int shuffle(const void *a, const void *b) {
  return (rand() & 0x01) ? -1 : 1;
}

int main(int argc, char **argv) {
  unsigned int idx[ARRAY_SIZE];
  double a[ARRAY_SIZE];
  unsigned int i;

  for(i = 0; i < ARRAY_SIZE; i++)
    idx[i] = i;

  qsort(idx, ARRAY_SIZE, sizeof(int), shuffle);

  for(i = 0; i < ARRAY_SIZE; i++)
    a[idx[i]] = i;

  reorder(a, idx, ARRAY_SIZE);

  return 0;
}
