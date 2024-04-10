#include <fparse.h>
#include <stdlib.h> // size_t ssize_t

int fparse(FILE *in, const char sep, const fparse_fct_t fct, void *arg) {
  char *line = NULL, *s;
  size_t len = 0;
  ssize_t read = 0;

  unsigned int iline = 0;
  char **fields = NULL;
  unsigned int nfield = 0;

  unsigned int i;


  // Read each line from file
  while ((read = getline(&line, &len, in)) != -1) {
    // Remove the trailing '\n' character (if any)
    if(line[read-1] == '\n')
      line[read-1] = '\0';

    // Count the number of field(s) in this string
    for(nfield = 1, s = line; *s; s++)
      nfield += (*s == sep);

    // (Re-)allocate the fields array
    fields = (char **) realloc(fields, nfield * sizeof(char *));
    if(fields == NULL) {
      free(line);
      return -1;
    }

    // Initialize the fields array
    for(i = 0, s = line; *s; s++) {
      fields[i++] = s;
      while(*s && *s != sep)
        s++;
      *s = '\0';
    }

    // Callback function
    if(fct(iline++, (const char **) fields, nfield, arg)) {
      free(fields);
      free(line);
      return -1;
    }
  }

  free(fields);
  free(line);

  return 0;
}
