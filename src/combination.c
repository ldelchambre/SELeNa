#include <stdlib.h>

int combine(const unsigned int n, const unsigned int k, void (*fct)(const unsigned int *, const unsigned int, void *), void *arg) {
  unsigned int *idx;
  unsigned int i, j;
  
  // Exit with success if nothing has to be produced
  if(k == 0 || n < k)
    return 0;
  
  // Allocate the array of indices
  idx = (unsigned int *) calloc(k, sizeof(unsigned int));
  if(idx == NULL)
    return -1;
  
  // Fill the initial array of indices  
  // idx = [0, ..., k-1]
  for(i = 0; i < k; i++)
    idx[i] = i;
    
  // Browse each combination
  while(idx[k-1] != n) {
    // Callback function
    fct(idx, k, arg);
    
    // Increase the last index if allowed
    if(idx[k-1] < n - 1)
      idx[k-1]++;
    else { // Otherwise
      // Search for the last index satisfying idx[j-1] < idx[j] - 1
      for(j = k - 1; j && idx[j-1] == idx[j] - 1; j--)
        ;
      // If we do not find it, say idx[0] is sufficient
      if(j == 0) j++;
      
      // Update the array of indices
      idx[j-1]++;
      for(; j < k; j++)
        idx[j] = idx[j-1] + 1; 
    }
  }
  
  // Delete the array of indices
  free(idx);
  
  return 0;
}
