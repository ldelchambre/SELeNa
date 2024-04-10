#ifndef _FLOAD_H_
#define _FLOAD_H_

#include <stdio.h>
#include <stdlib.h>

struct fload_data;

typedef struct fload_data fload_data_t;

fload_data_t *fload(FILE *in, const char sep);
size_t fload_nobs(const fload_data_t *data);
size_t fload_nattr(const fload_data_t *data);
const char *fload_get(const fload_data_t *data, const size_t iobs, const size_t iattr);
void fload_destroy(fload_data_t *data);

#endif
