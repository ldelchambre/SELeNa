#ifndef _RFPRINTF_H_
#define _RFPRINTF_H_

#include <stdio.h>

void rfprintf_reset();
int rfprintf(FILE *stream, const char *format, ...);

#endif
