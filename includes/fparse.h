#ifndef _FPARSE_H_
#define _FPARSE_H_

#include <stdio.h>

typedef int (*fparse_fct_t)(const unsigned int iline,
                            const char **fields,
                            const unsigned int nfield,
                            void *arg);

int fparse(FILE *in, const char sep, const fparse_fct_t fct, void *arg);

#endif
