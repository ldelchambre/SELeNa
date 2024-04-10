#include <combination.h>
#include <stdio.h>

#define TEST_N 5
#define TEST_K 3

const unsigned int test[10][TEST_K] = {
  {0, 1, 2},
  {0, 1, 3},
  {0, 1, 4},
  {0, 2, 3},
  {0, 2, 4},
  {0, 3, 4},
  {1, 2, 3},
  {1, 2, 4},
  {1, 3, 4},
  {2, 3, 4}
};

void fct(const unsigned int *idx, const unsigned int k, void *arg) {
  unsigned int *i = (unsigned int *) arg;
  unsigned int j;
  
  for(j = 0; j < k; j++) {
    if(idx[j] != test[*i][j]) {
      fprintf(stderr, "Wrong combination: combination %u is ( ", *i);
      for(j = 0; j < k; j++) fprintf(stderr, "%u ", idx[j]);
      fprintf(stderr, ") instead of ( ");
      for(j = 0; j < k; j++) fprintf(stderr, "%u ", test[*i][j]);
      fprintf(stderr, ")\n");
    }
  }
  
  (*i)++;
}

int main(void) {
  unsigned int i = 0;
  
  combine(TEST_N, TEST_K, fct, &i);
}
