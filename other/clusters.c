#include <stdio.h>
#include <stdlib.h>

#define NELEM 16

unsigned long get_subcluster(unsigned long *icluster, const char *_distm, const unsigned long iobj, const unsigned long nobj) {
  char (*distm)[nobj] = (char (*)[nobj]) _distm;
  unsigned long i, istart, iinsert, iend;
  
  istart = iend = 0;
  for(i = iobj; i < nobj; i++)
    icluster[iend++] = i;
  
  for(istart = 0; istart < iend; istart++) {
    for(i = iinsert = istart+1; i < iend; i++) {
      if(distm[icluster[istart]][icluster[i]])
        icluster[iinsert++] = icluster[i];
    }
    iend = iinsert;
  }
  
  return istart;
}

void fill_distm(char *_distm) {
  char (*distm)[NELEM] = (char (*)[NELEM]) _distm;
  unsigned int i,j;
  
  for(i = 0; i < NELEM; i++) {
    distm[i][i] = 1;
    for(j = 0; j < i; j++)
      distm[i][j] = distm[j][i] = (rand() < (RAND_MAX >> 1));
  }
}

void print_distm(const char *_distm, const unsigned long *icluster, const unsigned long n) {
  char (*distm)[NELEM] = (char (*)[NELEM]) _distm;
  unsigned long i,j;
  
  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++)
      fprintf(stdout, "%2d", distm[icluster[i]][icluster[j]]);
    fprintf(stdout, "\n");
  }
}

int main(void) {
  char distm[NELEM*NELEM];
  unsigned long icluster[NELEM];
  unsigned int i, j, n;
  
  fill_distm(distm);
    
  for(i = 0; i < NELEM; i++)
    icluster[i] = i;
  
  fprintf(stdout, "Initial matrix\n");
  print_distm(distm, icluster, NELEM);
  
  for(i = 0; i < NELEM; i++) {
    n = get_subcluster(icluster, distm, i, NELEM);
    fprintf(stdout, "Indices: ");
    for(j = 0; j < n; j++)
      fprintf(stdout, "%5lu", icluster[j]);
    fprintf(stdout, "\n");
    print_distm(distm, icluster, n);
  }
}
