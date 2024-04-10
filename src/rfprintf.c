/*
 * rfprintf.c
 *
 *  Created on: 21 août 2018
 *      Author: ludovic
 */
#include <rfprintf.h>
#include <stdarg.h>
#include <stdlib.h>

int _rfprintf_length = 0;

void rfprintf_reset() {
  _rfprintf_length = 0;
}

int rfprintf(FILE *stream, const char *format, ...) {
  int r;
  va_list ap;
  va_start(ap, format);
  fprintf(stream, "\r");
  r = vfprintf(stream, format, ap);
  if(r >= 0) {
    if(r < _rfprintf_length)
      fprintf(stream, "%*s", _rfprintf_length - r, "");
    else
      _rfprintf_length = r;
  }
  va_end(ap);
  return r;
}
