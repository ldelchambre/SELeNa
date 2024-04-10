/**
 * @file xmatch.c
 *
 * Perform a cross match of objects on the celestial sphere.
 *
 * Copyright 2018 Ludovic Delchambre <ldelchambre@uliege.be>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */
#include <astrostr.h>
#include <ctype.h>
#include <errno.h>
#include <getopt.h>
#include <htmdb.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#ifdef __GNUC__
#define UNUSED(param) __attribute__((unused)) param
#else
#define UNUSED(param) param
#endif

/*************************************************************************************************
 * Units
 ************************************************************************************************/
#define ARCSEC 4.84813681109536e-6

/*************************************************************************************************
 * Default program parameters
 ************************************************************************************************/
#define DEFAULT_RADIUS (5.*ARCSEC)
#define DEFAULT_LEVEL 10U

/*************************************************************************************************
 * Variable containing program parameters
 ************************************************************************************************/
char all = 0;
unsigned int level = DEFAULT_LEVEL;
double r = DEFAULT_RADIUS;
char *dmag_arg = NULL;
char *pos0 = NULL;
char *pos1 = NULL;
FILE *buffer = NULL; // fopen(buffer_path, "w+");
char printd = 0;

double *dmag = NULL;
unsigned int nmag = 0;

char *buffer_path = NULL;
double cosr; // cos(r)

char verbose = 0;

/*************************************************************************************************
 * Help message
 ************************************************************************************************/
static const char
    *help_message =
        "USAGE:\n"
        "  xmatch [options] input0 input1\n"
        "\n"
        "DESCRIPTION:\n"
        "  Perform a cross match between celestial objects from input0 and celestial objects from input1.\n"
        "\n"
        "  Objects from input1 are matched to objects of input0 based on their spherical coordinates and absolute\n"
        "  differences in magnitudes. If multiple objects are found within the specified radius (i.e. -r option)\n"
        "  and within the given magnitude range (i.e. -m option) then either the one having a minimal\n"
        "  great-circle distance will be returned or all objects will be returned if the -a option is provided.\n"
        "\n"
        "  If input0 or input1 (not both) is set to \"-\", then coordinates are read from standard input.\n"
        "\n"
        "OPTIONS:\n"
        "  -r, --radius,  The minimal distance between two celestial object in order to consider them as the same object [default = %garcsec]\n"
        "  -m, --dmag,    The maximal absolute difference in magnitude in order to consider two object as similar [default = none]\n"
        "                 Multiple value can be provided as comma seprarated values (e.g. -m \"1.5,0.25,6.25\").\n"
        "                 See notes for more informations\n"
        "  -a, --all,     Should we return all objects falling wihtin the given radius and magnitude range?\n"
        "  -d, --printd,  Print the great circle distance in the first output column.\n"
        "  -l, --level,   HTM level to use [default = %u].\n"
        "                 A high HTM level ensure fast execution at the expense of a high memory consumption.\n"
        "                 Low HTM level typically degrade performances but consumes less memory.\n"
        "  -b, --buffer,  Path to the buffer file to use. If specified, input files will written to the specified file, hence allowing the\n"
        "                 program to deal with very large files large; otherwise files are loaded into memory, improving the execution time.\n"
        "  -v, --verbose, Should we have to be verbose?\n"
        "  -?, --help,    Print this help message then exit.\n"
        "\n"
        "NOTES:\n"
        "  The search radius can be expressed in degree, arcmin, arcsec, mas, hour, minute or second (e.g. \"-r 100mas\").\n"
        "  Shorthand to these units are d, ', \", h, m, s respectively for degree, arcmin, arcsec, hour, minute and second.\n"
        "\n"
        "  Right ascension and declination should be contained as the two first columns of both input0 and input1\n"
        "  By default, right ascension and declination are considered to be expressed in radian, if followed by 'd' or 'D', then\n"
        "  these are considered to be in degree whereas if right ascension is in the form \"hh:mm:ss\", sexagesimal coordinates\n"
        "  are supposed corresponding to hour, minute and second. Declination in the form \"dd:mm:ss\" are understood as\n"
        "  sexagesimal coordinates in degree, arcmin and arcsec.\n"
        "\n"
        "  If -m parameter is set, then subsequent columns will be interpreted as magnitudes (e.g. if -m \"2,3,5\") then the input\n"
        "  files must be in the form \"ra dec mag0 mag1 mag2 ...\" for which the absolute difference in magnitude are restricted\n"
        "  to be below 2, 3 and 5, respectively for mag0, mag1 and mag2.\n"
        "\n";

void printHelp(FILE *out) {
  fprintf(out, help_message, DEFAULT_RADIUS/ARCSEC, DEFAULT_LEVEL);
}

/*************************************************************************************************
 * Command line parameter
 ************************************************************************************************/
static struct option options_long[] = {
/*-r*/{ "radius", required_argument, 0, 'r' },
/*-m*/{ "dmag", required_argument, 0, 'm' },
/*-a*/{ "all", no_argument, 0, 'a' },
/*-d*/{ "printd", no_argument, 0, 'd' },
/*-l*/{ "level", required_argument, 0, 'l' },
/*-b*/{ "buffer", required_argument, 0, 'b' },
/*-v*/{ "verbose", no_argument, 0, 'v' },
/*-?*/{ "help", no_argument, 0, '?' },
/*  */{ 0, 0, 0, 0 }
};

static const char options_short[] = "r:m:adl:b:v?";

/*******************************************************************************
 * Progression report functions
 ******************************************************************************/
typedef struct progress_data {
  int finished;
  double t0;
  unsigned long long n0;
  unsigned long long n1;
  unsigned long long nfound;
  void (*print_fct)(void);
  pthread_t thread_id;
} progress_data_t;

progress_data_t progress = {0, 0, 0, 0, 0, NULL, 0};

void *progress_thread(void UNUSED(*arg)) {
  while(!progress.finished) {
    sleep(1); // Print progression report every second
    if(progress.print_fct) progress.print_fct();
  }

  if(progress.print_fct) progress.print_fct();

  pthread_exit(NULL);
}

double cluster_gettime() {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return (double) t.tv_sec + 1e-9 * t.tv_nsec;
}

void progress_print_add(void) {
  const double dt = cluster_gettime() - progress.t0;
  fprintf(stderr, "\r%lld objects added to the database in %fs",
                  progress.n0, dt);
}

void progress_print_browse(void) {
  const double dt = cluster_gettime() - progress.t0;
  fprintf(stderr, "\r%llu matches found over %llu objects, run time: %.3fs",
                  progress.nfound, progress.n1, dt);
}

/*************************************************************************************************
 * Program functions
 ************************************************************************************************/
typedef struct xmatch_data {
  double *mags;
  char *line;
} xmatch_data_t;

void xmatch_free_data(xmatch_data_t *data) {
  if(data) {
    free(data->line);
    free(data->mags);
    free(data);
  }
}

void *xmatch_set_data(xmatch_data_t *data) {
  // If we have to use a buffer
  if(buffer) {
    const size_t line_size = strlen(data->line);
    long buffer_pos;
    
    // Set the file pointer to the end of the file
    if(fseek(buffer, 0, SEEK_END)) {
      fprintf(stderr, "Unable to reach the end of file\n");
      return NULL;
    }
    
    // Get the position within buffer
    buffer_pos = ftell(buffer);
    if(buffer_pos < 0) {
      fprintf(stderr, "Unable get the buffer position\n");
      return NULL;
    }
    
    // Write the xmatch_data_t structure in buffer
    if(fwrite(&line_size, sizeof(size_t), 1, buffer) == 0 ||
       ( line_size > 0 && fwrite(data->line, sizeof(char), line_size, buffer) == 0 ) ||
       ( nmag > 0 && fwrite(data->mags, sizeof(double), nmag, buffer) == 0) ) {
      fprintf(stderr, "Unable to write xmatch data structure to buffer\n");
      return NULL;   
    }
    
// If long values can fit into an intptr_t type  
#if LONG_MAX <= INTPTR_MAX
    // Directly return buffer_pos, cast as a pointer
    return (void *) ((intptr_t) buffer_pos);
#else
    // Otherwise, allocate a long to store the buffer position, then return it
    long *ret = (long *) malloc(sizeof(long));
    if(ret == NULL) {
      fprintf(stderr, "Unable to allocate the buffer position\n");
      return NULL;
    }
    *ret = buffer_pos;
    return ret;
#endif
  }
  return data;
}

xmatch_data_t *xmatch_get_data(void *arg) {
  // If we have to use a buffer
  if(buffer) {
    xmatch_data_t *data;
    size_t line_size;
    
  // If long values can fit into an intptr_t type
#if LONG_MAX <= INTPTR_MAX
    // Get the position within buffer
    const long buffer_pos = (const intptr_t) arg;
#else
    // Otherwise, the position is allocated as a long value
    const long buffer_pos = *((const long *) arg);
#endif
    // Seek position within buffer
    if(fseek(buffer, buffer_pos, SEEK_SET) != 0) {
      fprintf(stderr, "Unable to seek position within buffer\n");
      return NULL;
    }
    
    // Read the line size
    if(fread(&line_size, sizeof(size_t), 1, buffer) == 0) {
      fprintf(stderr, "Unable to read line size from buffer\n");
      return NULL;
    }
    
    // Allocate the xmatch_data_t structure
    data = malloc(sizeof(xmatch_data_t));
    if(data) {
      data->line = (char *) calloc(line_size+1,sizeof(char));
      data->mags = (double *) calloc(nmag, sizeof(double));
    }
    if(!data || !data->line || !data->mags ) {
      fprintf(stderr, "Unable to allocate the xmatch data structure");
      xmatch_free_data(data);
      return NULL;
    }
    
    // Read line and magnitudes
    if( (line_size > 0 && fread(data->line, sizeof(char), line_size, buffer) == 0) ||
        (nmag > 0 && fread(data->mags, sizeof(double), nmag, buffer) == 0) ) {
      fprintf(stderr, "Unable to read the xmatch data structure from buffer\n");
      xmatch_free_data(data);
      return NULL;
    }
    
    return data;
  }
  // If no buffer has to be used, arg is a pointer to the xmatch_data_t structure
  return (xmatch_data_t *) arg;
}

void xmatch_release_data(xmatch_data_t *data) {
  // If we have to use a buffer, then the xmatch_data_t structure was allocated
  // and read from buffer such that we have to delete it
  if(buffer)
    xmatch_free_data(data);
}

void xmatch_destroy_data(void *arg) {
#if LONG_MAX > INTPTR_MAX
  if(buffer)
    free(arg);
  else
#else
  if(!buffer)
#endif
    xmatch_free_data(arg);
}

int xmatch_parse_file(const char *filename, int (*process)(const double, const double, const double *, const char *, void *), void *arg) {
  // Open file in read mode
  FILE *in = (strcmp(filename,"-") == 0) ? stdin : fopen(filename, "r");

  // Object position
  double ra, dec;

  // Object magnitudes
  double *mags;

  // Variables of readline
  char *line = NULL;
  size_t len = 0;
  ssize_t read;

  // If the file opening failed
  if(!in) {
    fprintf(stderr, "Unable to open file %s in read mode (%s)\n", filename, strerror(errno));
    return -1;
  }

  // Allocate the array of magnitudes (if nmag == 0, then mags == NULL)
  mags = calloc(nmag, sizeof(double));
  if(!mags) {
    fprintf(stderr, "Unable to allocate to array of magnitudes\n");
    fclose(in);
    return -1;
  }

  // Browse each line of file
  while((read = getline(&line, &len, in)) != -1) {
    char *s, *endptr;
    unsigned int i;

    // Remove the '\n' end character (if any)
    if(line[read-1] == '\n')
      line[read-1] = '\0';

    // Extract RA
    if(str2ra(&ra,line,&endptr)) {
      fprintf(stderr, "Unable to get right ascension from line \"%s\"\n", line);
      free(line);free(mags);fclose(in);
      return -1;
    }

    // Extract DEC
    if(str2dec(&dec,endptr,&endptr)) {
      fprintf(stderr, "Unable to get declination from line \"%s\"\n", line);
      free(line);free(mags);fclose(in);
      return -1;
    }

    // Extract magnitudes
    s = endptr;
    for(i = 0; i < nmag; i++) {
      mags[i] = strtod(s, &endptr);
      if(s == endptr) {
        fprintf(stderr, "Unable to extract magnitudes from line \"%s\"\n", line);
        free(line);free(mags);fclose(in);
      }
      s = endptr;
    }

    // Process this object
    if(process(ra, dec, mags, line, arg)) {
      free(mags);free(mags);free(line);
      return -1;
    }
  }

  free(line);
  free(mags);
  fclose(in);

  return 0;
}

int xmatch_add_to_db(const double ra, const double dec, const double *mags, const char *line, void *arg) {
  htmdb_t *db = (htmdb_t *) arg;
  double v[3];
  xmatch_data_t *data;
  void *htmdb_data;

  // Allocate the data associated with the cross-match
  data = (xmatch_data_t *) malloc(sizeof(xmatch_data_t));
  if(!data) {
    fprintf(stderr, "Unable to allocate data instance for line \"%s\"\n", line);
    return -1;
  }

  // Allocate line and magnitude array
  data->line = calloc(strlen(line)+1, sizeof(char));
  data->mags = calloc(nmag, sizeof(double));
  if(!data->line || !data->mags) {
    fprintf(stderr, "Unable to allocate data content of line \"%s\"\n", line);
    xmatch_free_data(data);
    return -1;
  }

  // Copy line and magnitudes in the data structure
  strcpy(data->line, line);
  {
    unsigned int i;
    for(i = 0; i < nmag; i++) data->mags[i] = mags[i];
  }

  // Convert spherical coordinates into cartesian coordinates
  htm_cart_coord(v, ra, dec);
  
  // Convert the xmatch data structure as htmdb data
  htmdb_data = xmatch_set_data(data);

  // Add this entry to the database
  if(!htmdb_data || htmdb_add(db, v, htmdb_data)) {
    fprintf(stderr, "Unable to add line \"%s\" to the database\n", line);
    xmatch_free_data(data);
    return -1;
  }
  
  // Release data
  xmatch_release_data(data);

  // Increase the number of objects in the database
  progress.n0++;

  return 0;
}

int xmatch_dmag_match(const double *m0, const double *m1) {
  unsigned int i;

  for(i = 0; i < nmag; i++) {
    if(fabs(m0[i]-m1[i]) > dmag[i])
      return 0;
  }

  return 1;
}

typedef struct xmatch_best_match {
  double v[3];
  const double *mags;
  char *line;
  double cost;
} xmatch_best_match_t;

int _xmatch_find_best_match(const double *v, void *data, void *arg) {
  xmatch_data_t *obj = xmatch_get_data(data);
  xmatch_best_match_t *bmatch = (xmatch_best_match_t *) arg;
  const double cost = v[0]*bmatch->v[0]+v[1]*bmatch->v[1]+v[2]*bmatch->v[2];
  
  if(!obj) {
    fprintf(stderr, "Unable to retrieve object\n");
    return -1;
  }

  if(bmatch->cost < cost && xmatch_dmag_match(bmatch->mags, obj->mags)) {
    bmatch->cost = cost;
    // Copy line
    free(bmatch->line);
    bmatch->line = (char *) calloc(strlen(obj->line)+1,sizeof(char));
    if(bmatch->line == NULL) {
      fprintf(stderr, "Unable to copy line \"%s\"\n", obj->line);
      xmatch_release_data(obj);
      return -1;
    }
    strcpy(bmatch->line, obj->line);
  }
  
  xmatch_release_data(obj);

  return 0;
}

int xmatch_find_best_match(const double ra, const double dec, const double *mags, const char *line, void *arg) {
  htmdb_t *db = (htmdb_t *) arg;
  xmatch_best_match_t bmatch;

  // Convert spherical coordinates into cartesian coordinates
  htm_cart_coord(bmatch.v, ra, dec);

  // Initialize the best match structure
  bmatch.line = NULL;
  bmatch.mags = mags;
  bmatch.cost = -1;

  // Increase the number of objects we browsed
  progress.n1++;

  // Perform constraint on the DB
  if(htmdb_constraint(db, bmatch.v, cosr, _xmatch_find_best_match, &bmatch)) {
    fprintf(stderr, "Unable to perform a constraint around point (%g, %g), line \"%s\"\n", ra, dec, line);
    return -1;
  }

  // If a match was found, print it
  if(bmatch.line) {
    // Increase the number of matches we found
    progress.nfound++;

    // If we have to print the great circle distance between the two objects
    if(printd)
      fprintf(stdout, "%g ", (bmatch.cost > 1.) ? 0 : acos(bmatch.cost)); // Print the great circle dist
    fprintf(stdout, "%s %s\n", bmatch.line, line);
    free(bmatch.line);
  }

  return 0;
}

typedef struct xmatch_all_match {
  double v[3];
  const double *mags;
  const char *line;
} xmatch_all_match_t;

int _xmatch_find_all_match(const double *v, void *data, void  *arg) {
  const xmatch_all_match_t *amatch = (xmatch_all_match_t *) arg;
  xmatch_data_t *obj = xmatch_get_data(data);
  
  if(!obj) {
    fprintf(stderr, "Unable to retrieve object\n");
    return -1;
  }

  // If the absolute difference in magnitude match
  if(xmatch_dmag_match(amatch->mags, obj->mags)) {
    // Increase the number of matches we found
    progress.nfound++;

    // If we have to print the great circle distance between the two objects
    if(printd) {
      double cost = v[0]*amatch->v[0]+v[1]*amatch->v[1]+v[2]*amatch->v[2];
      fprintf(stdout, "%g ", (cost > 1.) ? 0. : acos(cost));
    }

    fprintf(stdout, "%s %s\n", obj->line, amatch->line);
  }
  
  xmatch_release_data(obj);

  return 0;
}

int xmatch_find_all_match(const double ra, const double dec, const double *mags, const char *line, void *arg) {
  htmdb_t *db = (htmdb_t *) arg;
  xmatch_all_match_t amatch;

  // Convert spherical coordinates into cartesian coordinates
  htm_cart_coord(amatch.v, ra, dec);

  // Initialize the all match structure
  amatch.mags = mags;
  amatch.line = line;

  // Increase the number of objects we browsed
  progress.n1++;

  // Perform constraint on the DB
  if(htmdb_constraint(db, amatch.v, cosr, _xmatch_find_all_match, &amatch)) {
    fprintf(stderr, "Unable to perform a constraint around point (%g, %g), line \"%s\"\n", ra, dec, line);
    return -1;
  }

  return 0;
}

/*************************************************************************************************
 * Main program
 ************************************************************************************************/
int main(int argc, char **argv) {
  htmdb_t *db;
  int opt;

  // Get the command line options
  while ((opt = getopt_long(argc, argv, options_short, options_long, NULL)) != -1) {
    switch (opt) {
      case 'r':
        if(str2angle(&r, optarg, NULL)) {
          fprintf(stderr, "Unrecognized angle format: \"%s\"\n", optarg);
          printHelp(stderr);
          return -1;
        }
        break;
      case 'm':
        dmag_arg = optarg;
        break;
      case 'a':
        all = 1;
        break;
      case 'd':
        printd = 1;
        break;
      case 'l':
        level = strtoul(optarg, NULL, 10);
        break;
      case 'b':
        buffer_path = optarg;
        break;
      case 'v':
        verbose = 1;
        break;
      case '?':
        printHelp(stdout);
        exit(EXIT_SUCCESS);
        break;
      default:
        printHelp(stderr);
        exit(EXIT_FAILURE);
        break;
    }
  }

  // Compute cos(r)
  cosr = cos(r);

  // Get the absolute difference in magnitude
  if(dmag_arg){
    const char *s = dmag_arg;
    char *endptr;
    unsigned int i;

    // Count the number of magnitudes in the --dmag option
    nmag = 1;
    while(*s)
      nmag += (*(s++) == ',');

    // Allocate the array of absolute difference in magnitudes
    dmag = calloc(nmag, sizeof(double));
    if(!dmag) {
      fprintf(stderr, "Unable to allocate the array of absolute difference in magnitudes\n");
      exit(EXIT_FAILURE);
    }

    // Fill the array of absolute difference in magnitude
    s = dmag_arg;
    for(i = 0; i < nmag; i++) {
      // Try to convert the ith absolute difference in magnitude
      dmag[i] = strtod(s, &endptr);
      // If nothing was converted
      if(s == endptr) {
        fprintf(stderr, "Empty absolute difference in magnitude at index %u\n", i+1);
        free(dmag);
        exit(EXIT_FAILURE);
      }
      // Skip trailing white spaces
      while(isspace(*endptr))
        endptr++;
      // If the next character to convert is not ','
      if((i + 1 < nmag && *endptr != ',') || (i + 1 == nmag && *endptr != '\0')) {
        fprintf(stderr, "Unrecognized absolute difference in magnitude at index %d \"%s\"\n", i+1, s);
        free(dmag);
        exit(EXIT_FAILURE);
      }
      // All is fine so go to next absolute difference in magnitude
      s = endptr + 1;
    }
  }

  // Get pos0 file name
  if(optind < argc)
    pos0 = argv[optind++];
  else {
    fprintf(stderr, "No input file found\n");
    printHelp(stderr);
    free(dmag);
    return -1;
  }

  // Get pos1 file name
  if(optind < argc)
    pos1 = argv[optind++];
  else {
    fprintf(stderr, "No file to cross match\n");
    printHelp(stdout);
    free(dmag);
    exit(EXIT_FAILURE);
  }

  // Check that pos0 and pos1 are not set to "-" together
  if(strcmp(pos0, "-") == 0 && strcmp(pos1, "-") == 0) {
    fprintf(stderr, "input0 and intput1 cannot be set together to stdin\n");
    printHelp(stdout);
    free(dmag);
    exit(EXIT_FAILURE);
  }

  // In verbose mode, create a watching thread, designed to print progression
  // reports
  if(verbose && pthread_create(&progress.thread_id, NULL, progress_thread, NULL)) {
    fprintf(stderr, "Unable to create the progression thread\n");
    free(dmag);
    exit(EXIT_FAILURE);
  }
  
  // Open the buffer file if asked to do so, then write nmag
  if(buffer_path) {
    buffer = fopen(buffer_path, "w+");
    if(buffer_path && (!buffer || fwrite(&nmag, sizeof(unsigned int), 1, buffer) == 0)) {
      fprintf(stderr, "Unable to open the buffer file\n");
      free(dmag);
      exit(EXIT_FAILURE);
    }
  }

  // Initialize the HTM db
  db = htmdb_init(level);
  if(!db) {
    fprintf(stderr, "Unable to initialize HTM db with level %u\n", level);
    free(dmag);
    if(buffer){fclose(buffer);unlink(buffer_path);}
    exit(EXIT_FAILURE);
  }

  // Fill the HTM db using pos0
  if(verbose) {
    fprintf(stderr, "Building the HTM database...\n");
    progress.print_fct = progress_print_add;
    progress.t0 = cluster_gettime();
  }
  if(xmatch_parse_file(pos0, xmatch_add_to_db, db)) {
    htmdb_destroy(db,xmatch_destroy_data);free(dmag);
    if(buffer){fclose(buffer);unlink(buffer_path);}
    exit(EXIT_FAILURE);
  }
  if(verbose) {
    progress.print_fct = NULL;
    progress_print_add();
    fprintf(stderr, "\n");
  }

  // Perform the cross match from pos1
  {
    int (*xmatch_find_match)(const double, const double, const double *, const char *, void *);

    if(all)
      xmatch_find_match = xmatch_find_all_match;
    else
      xmatch_find_match = xmatch_find_best_match;

    if(verbose) {
      fprintf(stderr, "Browsing the HTM database...\n");
      progress.print_fct = progress_print_browse;
      progress.t0 = cluster_gettime();
    }
    if(xmatch_parse_file(pos1, xmatch_find_match, db)) {
      htmdb_destroy(db,xmatch_destroy_data);free(dmag);
      if(buffer){fclose(buffer);unlink(buffer_path);}
      exit(EXIT_FAILURE);
    }
    if(verbose) {
      progress.print_fct = NULL;
      progress_print_browse();
      fprintf(stderr, "\n");
    }
  }

  // In verbose mode, terminate the thread
  if(verbose) {
    progress.finished = 1;
    pthread_join(progress.thread_id, NULL);
  }

  // Destroy the database, then exit
  htmdb_destroy(db,xmatch_destroy_data);
  free(dmag);
  if(buffer){fclose(buffer);unlink(buffer_path);}
  exit(EXIT_SUCCESS);
}
