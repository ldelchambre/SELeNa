/**
 * @file cluster.c
 *
 * Search for clusters of celestial objects
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
#include <combination.h>
#include <ctype.h>
#include <errno.h>
#include <getopt.h>
#include <htmdb.h>
#include <math.h>
#include <pthread.h>
#include <semaphore.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <threadpool.h>
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

/*******************************************************************************
 * Program parameters
 ******************************************************************************/
#define DEFAULT_SEP (5.*ARCSEC)
#define DEFAULT_LEVEL 10
#define DEFAULT_NMIN 3
#define DEFAULT_NMAX -1U
#define DEFAULT_RDENSITY 0.
#define DEFAULT_DMAX 1e8
#define DEFAULT_MAX_PARENT_SIZE -1U
#define DEFAULT_NTHREAD 1

double sep = DEFAULT_SEP;
double sep_cos /* = cos(sep)*/;
char *dmag_arg = NULL;
double *dmag = NULL;
double dmax = DEFAULT_DMAX;
unsigned int nmag = 0;
unsigned int level = DEFAULT_LEVEL;
unsigned int nmin = DEFAULT_NMIN;
unsigned int nmax = DEFAULT_NMAX;
double rdensity = DEFAULT_RDENSITY;
double rdensity_cos /* = cos(rdensity) */;
double rdensity_area /* = 4*M_PI*sin(rdensity/2)^2 */;
unsigned int max_parent_size = DEFAULT_MAX_PARENT_SIZE;
unsigned int nthread = DEFAULT_NTHREAD;
char header = 1;
char verbose = 0;

#define NOHEADER_PARAM 65540
#define NTHREAD_PARAM 65541

static struct option options_long[] = {
/*-s*/{ "sep", required_argument, 0, 's' },
/*-m*/{ "dmag", required_argument, 0, 'm' },
/*-n*/{ "nmin", required_argument, 0, 'n' },
/*-N*/{ "nmax", required_argument, 0, 'N' },
/*-d*/{ "rdensity", required_argument, 0, 'd' },
/*-D*/{ "dmax", required_argument, 0, 'D' },
/*-P*/{ "max-parent-size", required_argument, 0, 'P' },
/*  */{ "nthread", required_argument, 0, NTHREAD_PARAM },
/*  */{ "noheader", no_argument, 0, NOHEADER_PARAM },
/*-l*/{ "level", required_argument, 0, 'l' },
/*-v*/{ "verbose", no_argument, 0, 'v' },
/*-?*/{ "help", no_argument, 0, '?' },
/*  */{ 0, 0, 0, 0 } };

static const char options_short[] = "s:m:n:N:d:D:P:l:v?";

/*************************************************************************************************
 * Help message
 ************************************************************************************************/
static const char help_message[] = "USAGE:\n"
    "  cluster [options]\n"
    "\n"
    "DESCRIPTION:\n"
    "  Search for clusters of celestial objects.\n"
    "\n"
    "OPTIONS:\n"
    "  -s, --sep,      The maximal great circle distance between any two objects of the cluster [default = %garcsec]\n"
    "  -m, --dmag,     The maximal absolute difference in magnitude between cluster objects [default = none]\n"
    "                  Multiple value can be provided as comma seprarated values (e.g. -m \"1.5,0.25,6.25\").\n"
    "                  See notes for more informations\n"
    "  -n, --nmin,     The minimal number of objects in the cluster [default = %u]\n"
    "  -N, --nmax,     The maximal number of objects in the cluster [default = %u]\n"
    "  -d, --rdensity, The radius to use in order to estimate the density of objects in the vicinity of the candidate.\n"
    "                  If rdensity <= 0, then no density will be outputed [default = %g]\n"
    "  -D, --dmax,     The maximal density of objects sourrounding candidates [default = %g objects/sr]\n"
    "                  Use in conjuction with the --rdensity option.\n"
    "  -P, --max-parent-size, Do not consider combinations of clusters larger than max-parent-size in order to prevent\n"
    "                  too many clusters to be produced (e.g. the decomposition of a cluster composed of 30 objects leads\n"
    "                  to 27,405 clusters composed of 4 images) [default = %u].\n"
    "      --noheader, Don't print header line\n"
    "  -l, --level,    HTM level to use [default = %u].\n"
    "                  A high HTM level ensure fast execution at the expense of a high memory consumption.\n"
    "                  Low HTM level typically degrade performances but consumes less memory.\n"
    "      --nthread,  The number of thread to use [default=%u]\n"
    "  -v, --verbose,  Should we have to be verbose?\n"
    "  -?, --help,     Print this help message\n"
    "\n"
    "NOTE:\n"
    "  Right ascension and declination should be contained as the two first columns of the standard input"
    "\n"
    "  If -m parameter is set, then subsequent columns will be interpreted as magnitudes (e.g. if -m \"2,3,5\") then the input\n"
    "  file must be in the form \"ra dec mag0 mag1 mag2 ...\" for which the absolute difference in magnitude are restricted\n"
    "  to be below 2, 3 and 5, respectively for mag0, mag1 and mag2.\n"
    "\n"
    "  Input lines starting with a '*' will be ignored for building clusters but will be taken into account for computing the\n"
    "  density surrounding each cluster.\n"
    "\n"
    "  At the opposite, lines starting with a '!' will be used for building clusters but will be ignored for density estimation.\n"
    "  By default (lines starting with coordinates), input lines will be used both for estimating the field density and for \n"
    "  building clusters.\n"
    "\n"
    "ANGLE FORMAT:\n"
    "  Angles are expressed in radian by default. Other possibilities are degree,\n"
    "  arcmin, arcsec, mas, hour, minute or second (e.g. \"30arcsec\" is a valid input)\n"
    "  Shorthand to these units are d, ', \", h, m, s respectively for degree, \n"
    "  arcmin, arcsec, hour, minute and second (e.g. \"30d\").\n"
    "  Right ascension and declination can be expressed in sexagecimal notation\n"
    "  respectively through the \"hh:mm:ss\" and \"dd:mm:ss\" notations (e.g. \"10:45:57.14\")\n";

void printHelp(FILE *out) {
  fprintf(out, help_message, DEFAULT_SEP / ARCSEC,
  DEFAULT_NMIN, DEFAULT_NMAX,
  DEFAULT_RDENSITY, DEFAULT_DMAX, DEFAULT_MAX_PARENT_SIZE, DEFAULT_LEVEL,
  DEFAULT_NTHREAD);
}

/*******************************************************************************
 * Multi-threading attributes
 ******************************************************************************/
threadpool_t *thread_pool = NULL;
pthread_mutex_t thread_mutex;
sem_t thread_sem;

/*******************************************************************************
 * Progression report functions
 ******************************************************************************/
typedef enum progress_mode {
  PROGRESS_MODE_NONE, PROGRESS_MODE_ADD, PROGRESS_MODE_BROWSE
} progress_mode_t;

typedef struct progress_data {
  int finished;
  progress_mode_t mode;
  double t0;
  unsigned long long iobject;
  unsigned long long nobjects;
  unsigned long long ncluster;
  pthread_t thread_id;
} progress_data_t;

progress_data_t progress;

double progress_gettime() {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return (double) t.tv_sec + 1e-9 * t.tv_nsec;
}

void progress_print_add(void) {
  const double dt = progress_gettime() - progress.t0;
  fprintf(stderr, "\r%lld objects added to the database in %fs", progress.nobjects, dt);
}

void progress_print_browse(void) {
  const double dt = progress_gettime() - progress.t0;
  const double p = 100. * progress.iobject / progress.nobjects;
  fprintf(stderr, "\r%7.3f%%) %llu clusters found, run time: %.3fs", p, progress.ncluster, dt);
}

void progress_print(void) {
  switch (progress.mode) {
    case PROGRESS_MODE_NONE:
      break;
    case PROGRESS_MODE_ADD:
      progress_print_add();
      break;
    case PROGRESS_MODE_BROWSE:
      progress_print_browse();
      break;
  }
}

void *progress_thread(void UNUSED(*arg)) {
  while (!progress.finished) {
    sleep(1); // Print progression report every second
    progress_print();
  }

  pthread_exit(NULL);
}

int progress_start() {
  progress.finished = 0;
  progress.mode = PROGRESS_MODE_NONE;
  progress.t0 = progress_gettime();
  progress.iobject = progress.nobjects = progress.ncluster = 0;
  return pthread_create(&progress.thread_id, NULL, progress_thread, NULL);
}

void progress_set_mode(progress_mode_t mode) {
  progress_print(); // Print the last status of the previous mode
  progress.mode = mode;
  progress.t0 = progress_gettime();
}

int progress_stop() {
  progress.finished = 1;
  return pthread_join(progress.thread_id, NULL);
}

void progress_new_object() {
  progress.nobjects++;
}

void progress_next_object() {
  progress.iobject++;
}

void progress_new_cluster() {
  progress.ncluster++;
}

/*******************************************************************************
 * Program functions
 ******************************************************************************/
#define OBJ_TYPE_EMPTY   0x00
#define OBJ_TYPE_DENSITY 0x01
#define OBJ_TYPE_CLUSTER 0x02

typedef struct cluster_data {
  char obj_type;
  double *mags;
  char *line;
} cluster_data_t;

typedef struct object_list {
  const double *v;
  const cluster_data_t *data;
  struct object_list *next;
} object_list_t;

typedef struct cluster_add_arg {
  object_list_t **cluster;
  const cluster_data_t *data;
  unsigned long ncluster, nparent;
} cluster_add_arg_t;

void cluster_free_data(void *arg) {
  cluster_data_t *data = (cluster_data_t *) arg;
  free(data->line);
  free(data->mags);
  free(data);
}

htmdb_t *cluster_parse_file(FILE *in) {
  // The HTM databse
  htmdb_t *db;

  // Object position
  double ra, dec;
  double v[3];

  // Whether this object should be considered only for the computation of the density
  char obj_type = OBJ_TYPE_EMPTY;

  // The HTM data
  cluster_data_t *data;

  // Variables of readline
  char *line = NULL;
  size_t len = 0;
  ssize_t read;

  // Initialize the HTM data base
  db = htmdb_init(level);
  if (!db) {
    fprintf(stderr, "Unable to allocate the HTM data base (level=%u)\n", level);
    return NULL;
  }

  // Browse each line of file
  while ((read = getline(&line, &len, in)) != -1) {
    char *s, *scoords, *endptr = line;
    unsigned int i;

    // Remove the '\n' end character (if any)
    if (line[read - 1] == '\n')
      line[read - 1] = '\0';

    // Check whether this object is only needed for the computation of the density,
    // (lines start with a '*'), only needed for the computation of the clusters (
    // (lines start with a '!') or both (by default)
    while (*endptr && isspace(*endptr))
      endptr++;
    if (*endptr == '*') {
      obj_type = OBJ_TYPE_DENSITY;
      endptr++;
    } else if (*endptr == '!') {
      obj_type = OBJ_TYPE_CLUSTER;
      endptr++;
    } else
      obj_type = OBJ_TYPE_CLUSTER + OBJ_TYPE_DENSITY;

    // Only add this object if it is not designed to compute the density or
    // if we have to compute the density
    if ((obj_type & OBJ_TYPE_CLUSTER) != 0 || rdensity > 0) {
      // Save the position of the coordinate string
      while (*endptr && isspace(*endptr))
        endptr++;
      scoords = endptr;

      // Try to extract RA
      if (str2ra(&ra, endptr, &endptr)) {
        fprintf(stderr, "Unable to get right ascension from line \"%s\"\n", scoords);
        free(line);
        htmdb_destroy(db, cluster_free_data);
        return NULL;
      }

      // Try to extract DEC
      if (str2dec(&dec, endptr, &endptr)) {
        fprintf(stderr, "Unable to get declination from line \"%s\"\n", scoords);
        free(line);
        htmdb_destroy(db, cluster_free_data);
        return NULL;
      }

      // Convert spherical coordinates into cartesian coordinates
      htm_cart_coord(v, ra, dec);

      // Allocate the data associated with the cluster
      data = (cluster_data_t *) malloc(sizeof(cluster_data_t));
      if (!data) {
        fprintf(stderr, "Unable to allocate data instance for line \"%s\"\n", scoords);
        free(line);
        htmdb_destroy(db, cluster_free_data);
        return NULL;
      }

      // Fill the object's informations
      data->obj_type = obj_type;
      if ((obj_type & OBJ_TYPE_CLUSTER) == 0) { // If this object will not be used to compute clusters, don't fill its fields
        data->mags = NULL;
        data->line = NULL;
      } else { // If it is a normal object
        data->mags = calloc(nmag, sizeof(double));
        data->line = calloc(strlen(scoords) + 1, sizeof(char));
        if (!data->mags || !data->line) {
          fprintf(stderr, "Unable to allocate data content for line \"%s\"\n", line);
          cluster_free_data(data);
          free(line);
          htmdb_destroy(db, cluster_free_data);
          return NULL;
        }

        // Extract magnitudes
        s = endptr;
        for (i = 0; i < nmag; i++) {
          data->mags[i] = strtod(s, &endptr);
          if (s == endptr) {
            fprintf(stderr, "Unable to extract magnitudes from line \"%s\"\n", line);
            cluster_free_data(data);
            free(line);
            htmdb_destroy(db, cluster_free_data);
          }
          s = endptr;
        }

        // Copy line
        strcpy(data->line, scoords);
      }

      // Add this entry to the database
      if (htmdb_add(db, v, data)) {
        fprintf(stderr, "Unable to add line \"%s\" to the database\n", line);
        cluster_free_data(data);
        free(line);
        htmdb_destroy(db, cluster_free_data);
        return NULL;
      }

      // Update the progress monitor
      progress_new_object();
    }
  }

  // Delete the input line
  free(line);

  // Return the built database
  return db;
}

int cluster_dmag_match(const double *m0, const double *m1) {
  unsigned int i;

  for (i = 0; i < nmag; i++) {
    if (fabs(m0[i] - m1[i]) > dmag[i])
      return 0;
  }

  return 1;
}

void cluster_free(object_list_t *cluster) {
  object_list_t *next;

  while (cluster) {
    next = cluster->next;
    free(cluster);
    cluster = next;
  }
}

int cluster_add(const double *v, void *data, void *arg) {
  // The cluster building arguments
  cluster_add_arg_t *caa = (cluster_add_arg_t *) arg;

  // The reference object around which the constraint was done
  const cluster_data_t *rdata = caa->data;

  // The current matching object
  const cluster_data_t *cdata = (const cluster_data_t *) data;

  // Only consider cluster objects for which the absolute difference in magnitudes match
  if ((cdata->obj_type & OBJ_TYPE_CLUSTER) != 0 && cluster_dmag_match(rdata->mags, cdata->mags)) {
    // Only add clusters for which the data pointer is greater than the one
    // around which this constraint was built so as to avoid duplicate
    // cluster to appear.
    if (rdata < cdata) {
      // Get the list of already found objects
      object_list_t **cluster = caa->cluster;

      // Allocate the new object
      object_list_t *obj = malloc(sizeof(object_list_t));
      if (!obj) {
        fprintf(stderr, "Unable to allocate memory for cluster object \"%s\"\n", cdata->line);
        cluster_free(*cluster);
        *cluster = NULL;
        return -1;
      }

      // Fill the object information
      obj->v = v;
      obj->data = cdata;

      // Add this object to the list
      obj->next = *cluster;
      *cluster = obj;
      caa->ncluster++;
    }
    // This cluster is part of a larger cluster (i.e. parent cluster) where all
    // absolute difference in magnitudes match
    caa->nparent++;
  }

  // Return 1 if the parent cluster has a size that is too large wrt. max_parent_size;
  // otherwise return 0
  return (max_parent_size < caa->nparent) ? 1 : 0;
}

int cluster_count(const double UNUSED(*v), void *data, void *arg) {
  const cluster_data_t *cdata = (const cluster_data_t *) data;

  // Only consider density objects
  if ((cdata->obj_type & OBJ_TYPE_DENSITY) != 0) {
    unsigned long *nnear = (unsigned long *) arg;
    (*nnear)++;
  }

  return 0;
}

char *cluster_distm(const object_list_t *cluster, const unsigned long ncluster) {
  char (*distm)[ncluster] = (char (*)[ncluster]) calloc(ncluster * ncluster, sizeof(char));
  const object_list_t *o1, *o2;
  unsigned long i, j;

  if (!distm)
    return NULL;

  for (i = 0, o1 = cluster; o1; i++, o1 = o1->next) {
    distm[i][i] = 1; // Distance bewteen o1 and o1 is always fine
    for (j = 1, o2 = cluster->next; o2; j++, o2 = o2->next) {
      // Set the distance to 1 if the angular distance between o1 and o2
      // is less than sep and that the absolute differences in magnitudes match
      const double cost = o1->v[0] * o2->v[0] + o1->v[1] * o2->v[1] + o1->v[2] * o2->v[2];
      distm[i][j] = distm[j][i] = (cost >= sep_cos && cluster_dmag_match(o1->data->mags, o2->data->mags));
    }
  }

  return (char *) distm;
}

void cluster_identifier(char *s, const double ra, const double dec) {
  char stmp[14];

  // Build RA in the form "hh:mm:ss.sss"
  ra2str(stmp, ra);

  // Convert RA in the form "hhmmsssss"
  s[0] = stmp[0];  // h
  s[1] = stmp[1];  // h
  s[2] = stmp[3];  // m
  s[3] = stmp[4];  // m
  s[4] = stmp[6];  // s
  s[5] = stmp[7];  // s
  s[6] = stmp[9];  // s
  s[7] = stmp[10]; // s
  s[8] = stmp[11]; // s

  // Build DEC in the form "+dd:mm:ss.sss"
  dec2str(stmp, dec);

  s[9] = stmp[0]; // +
  s[10] = stmp[1]; // d
  s[11] = stmp[2]; // d
  s[12] = stmp[4]; // m
  s[13] = stmp[5]; // m
  s[14] = stmp[7]; // s
  s[15] = stmp[8]; // s
  s[16] = stmp[10]; // s
  s[17] = stmp[11]; // s
  s[18] = stmp[12]; // s

  s[19] = '\0';
}

void subcluster_print(const double *v, const cluster_data_t *data, const object_list_t *cluster,
                      const unsigned long *icluster, const unsigned long ncluster, const double density) {
  const object_list_t *obj;
  double ra, dec;
  char cluster_id[20];
  unsigned long i, j;

  { // Get the sub-cluster identifier
    double n, vmean[3];

    // Compute the mean position of the cluster
    vmean[0] = v[0];
    vmean[1] = v[1];
    vmean[2] = v[2];
    for (i = j = 0, obj = cluster; i < ncluster; j++, obj = obj->next) {
      if (icluster[i] == j) {
        vmean[0] += obj->v[0];
        vmean[1] += obj->v[1];
        vmean[2] += obj->v[2];
        i++;
      }
    }
    n = sqrt(vmean[0] * vmean[0] + vmean[1] * vmean[1] + vmean[2] * vmean[2]);
    vmean[0] /= n;
    vmean[1] /= n;
    vmean[2] /= n;

    // Convert cartesian coordinates into spherical coordinates
    htm_sphere_coord(&ra, &dec, vmean);

    // Get the cluster identifier
    cluster_identifier(cluster_id, ra, dec);
  }

  // Acquire the lock
  pthread_mutex_lock(&thread_mutex);

  // Print the object around which the search was performed
  if (rdensity > 0)
    fprintf(stdout, "%19s %4lu %12.6e %s\n", cluster_id, ncluster + 1, density, data->line);
  else
    fprintf(stdout, "%19s %4lu %s\n", cluster_id, ncluster + 1, data->line);

  // Browse all other objects from the sub-cluster
  for (i = j = 0; i < ncluster; j++, cluster = cluster->next) {
    if (icluster[i] == j) {
      if (rdensity > 0)
        fprintf(stdout, "%19s %4lu %12.6e %s\n", cluster_id, ncluster + 1, density, cluster->data->line);
      else
        fprintf(stdout, "%19s %4lu %s\n", cluster_id, ncluster + 1, cluster->data->line);
      i++;
    }
  }

  // Increase the number of cluster we found so far
  progress_new_cluster();

  // Release the lock
  pthread_mutex_unlock(&thread_mutex);
}

typedef struct subcluster_arg {
  const double *v;
  const cluster_data_t *data;
  const object_list_t *cluster;
  const char *distm;
  unsigned long ncluster;
  double density;
} subcluster_arg_t;

void subcluster_process(const unsigned int *cidx, const unsigned int n, void *arg) {
  const subcluster_arg_t *sarg = (const subcluster_arg_t *) arg;
  char (*distm)[sarg->ncluster] = (char (*)[sarg->ncluster]) sarg->distm;
  unsigned long icluster[n];
  char pass = 1;
  unsigned long i, j;

  // Check that this sub-cluster match the maximal separation and maximal
  // magnitude difference criteria based on the distance matrix
  for (i = 0; i < n && pass; i++) {
    for (j = i + 1; j < n && pass; j++)
      pass &= (distm[cidx[i]][cidx[j]]);
    icluster[i] = cidx[i];
  }

  // If this subcluster pass the tests, then print it
  if (pass)
    subcluster_print(sarg->v, sarg->data, sarg->cluster, icluster, n, sarg->density);
}

typedef struct cluster_browse_db_thread_arg {
  const double *v;
  const cluster_data_t *cdata;
  const htmdb_t *db;
  int *error;
} cluster_browse_db_thread_arg_t;

void cluster_browse_db_thread(void *arg) {
  const double *v = ((cluster_browse_db_thread_arg_t *) arg)->v;
  const cluster_data_t *cdata = ((cluster_browse_db_thread_arg_t *) arg)->cdata;
  const htmdb_t *db = ((cluster_browse_db_thread_arg_t *) arg)->db;
  int *error = ((cluster_browse_db_thread_arg_t *) arg)->error;
  object_list_t *cluster = NULL;
  double density;
  unsigned long ncluster, nparent;
  char *distm;

  // Perform a constraint around this specific object so as to get cluster
  {
    cluster_add_arg_t caa;
    caa.cluster = &cluster;
    caa.data = cdata;
    caa.ncluster = caa.nparent = 0;
    if (htmdb_constraint(db, v, sep_cos, cluster_add, &caa) && caa.nparent <= max_parent_size) {
      fprintf(stderr, "Unable to set constraint around line \"%s\"\n", cdata->line);
      *error = -1;
      return;
    }
    ncluster = caa.ncluster;
    nparent = caa.nparent;
  }

  // If not enough object are present in this cluster or if the parent
  // cluster has a size larger than max_parent_size, skip it
  if (ncluster + 1 < nmin || max_parent_size < nparent) {
    cluster_free(cluster);
    return;
  }

  // Compute the density of objects in the vicinity of this object
  if (rdensity > 0) {
    unsigned long nnear = 0;
    if (htmdb_constraint(db, v, rdensity_cos, cluster_count, &nnear)) {
      fprintf(stderr, "Unable to estimate the density of objects around line \"%s\"\n", cdata->line);
      cluster_free(cluster);
      *error = -1;
      return;
    }
    density = 1. / rdensity_area * nnear;
  } else
    density = 0;

  // If the density of objects around this candidate is below dmax, then process
  // this candidate
  if (density < dmax) {
    // Compute the distance matrix between all objects from the cluster
    distm = cluster_distm(cluster, ncluster);
    if (!distm) {
      fprintf(stderr, "Unable to allocate distance matrix of the cluster around line \"%s\"\n", cdata->line);
      cluster_free(cluster);
      *error = -1;
      return;
    }

    // Browse all combinations of {nmin-1, ..., nmax-1} images from the cluster
    {
      const subcluster_arg_t sarg = { v, cdata, cluster, distm, ncluster, density };
      const unsigned int nlow = nmin - 1;
      const unsigned int nhigh = (ncluster < nmax) ? ncluster : nmax - 1;
      unsigned int k;

      for (k = nlow; k <= nhigh; k++)
        combine(ncluster, k, subcluster_process, (void *) &sarg);
    }

    // Destroy the distance matrix
    free(distm);
  }

  // Destroy cluster
  cluster_free(cluster);
}

void cluster_browse_db_thread_finish(void *arg) {
  free(arg);
  sem_post(&thread_sem);
}

int cluster_browse_db(const double *v, void *data, void *arg) {
  int error = 0;

  // Update the progress monitor
  progress_next_object();

  // Only consider clusters objects
  if ((((const cluster_data_t *) data)->obj_type & OBJ_TYPE_CLUSTER) == 0)
    return 0;

  // Launch a thread in order to process this entry
  cluster_browse_db_thread_arg_t *targ = (cluster_browse_db_thread_arg_t *) malloc(sizeof(cluster_browse_db_thread_arg_t));
  if (targ == NULL)
    return -1;
  else {
    targ->v = v;
    targ->cdata = (const cluster_data_t *) data;
    targ->db = (const htmdb_t *) arg;
    targ->error = &error;
    sem_wait(&thread_sem); // Ensure that at most nthread are running
    if (threadpool_addjob(thread_pool, cluster_browse_db_thread, targ, cluster_browse_db_thread_finish) != 0)
      return -1;
  }

  return 0;
}

/*******************************************************************************
 * Main function
 ******************************************************************************/
int main(int argc, char **argv) {
  htmdb_t *db; // The HTM database to use
  int opt;

  // Get the program parameters
  while ((opt = getopt_long(argc, argv, options_short, options_long, NULL)) != -1) {
    switch (opt) {
      case 's':
        str2angle(&sep, optarg, NULL);
        break;
      case 'm':
        dmag_arg = optarg;
        break;
      case 'n':
        nmin = strtoul(optarg, NULL, 10);
        break;
      case 'N':
        nmax = strtoul(optarg, NULL, 10);
        break;
      case 'd':
        str2angle(&rdensity, optarg, NULL);
        break;
      case 'D':
        dmax = strtod(optarg, NULL);
        break;
      case 'P':
        max_parent_size = strtoul(optarg, NULL, 10);
        break;
      case NOHEADER_PARAM:
        header = 0;
        break;
      case NTHREAD_PARAM:
        nthread = strtoul(optarg, NULL, 10);
        break;
      case 'l':
        level = strtoul(optarg, NULL, 10);
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

  // Check options
  if (nmin < 2) {
    fprintf(stderr, "nmin option must be greater than 1");
    exit(EXIT_FAILURE);
  }

  if (nmax < nmin) {
    fprintf(stderr, "nmax option must be greater or equal to nmin");
    exit(EXIT_FAILURE);
  }

  if (nthread == 0) {
    fprintf(stderr, "nthread option must be greater than zero");
    exit(EXIT_FAILURE);
  }

  // Compute the cosine of the search radius
  sep_cos = cos(sep);

  // Compute the cosine of the density radius and its area
  rdensity_cos = cos(rdensity);
  rdensity_area = 4. * M_PI * pow(sin(0.5 * rdensity), 2);

  // Get the absolute difference in magnitude
  if (dmag_arg) {
    const char *s = dmag_arg;
    char *endptr;
    unsigned int i;

    // Count the number of magnitudes in the --dmag option
    nmag = 1;
    while (*s)
      nmag += (*(s++) == ',');

    // Allocate the array of absolute difference in magnitudes
    dmag = calloc(nmag, sizeof(double));
    if (!dmag) {
      fprintf(stderr, "Unable to allocate the array of absolute difference in magnitudes\n");
      exit(EXIT_FAILURE);
    }

    // Fill the array of absolute difference in magnitude
    s = dmag_arg;
    for (i = 0; i < nmag; i++) {
      // Try to convert the ith absolute difference in magnitude
      dmag[i] = strtod(s, &endptr);
      // If nothing was converted
      if (s == endptr) {
        fprintf(stderr, "Empty absolute difference in magnitude at index %u\n", i + 1);
        free(dmag);
        exit(EXIT_FAILURE);
      }
      // Skip trailing white spaces
      while (isspace(*endptr))
        endptr++;
      // If the next character to convert is not ','
      if ((i + 1 < nmag && *endptr != ',') || (i + 1 == nmag && *endptr != '\0')) {
        fprintf(stderr, "Unrecognized absolute difference in magnitude at index %d \"%s\"\n", i + 1, s);
        free(dmag);
        exit(EXIT_FAILURE);
      }
      // All is fine so go to next absolute difference in magnitude
      s = endptr + 1;
    }
  }

  // Check for extra arguments
  if (optind < argc) {
    fprintf(stderr, "Warning: extra argument(s) found (");
    while (optind < argc)
      fprintf(stderr, " %s", argv[optind++]);
    fprintf(stderr, ")\n");
  }

  // Create a progress monitor
  if (verbose && progress_start() != 0) {
    fprintf(stderr, "Unable to create the progression thread\n");
    free(dmag);
    exit(EXIT_FAILURE);
  }

  // Initialize the thread pool, semaphore and mutex
  thread_pool = threadpool_create(nthread);
  if (thread_pool == NULL) {
    fprintf(stderr, "Unable to initialize the thread pool\n");
    free(dmag);
    exit(EXIT_FAILURE);
  }
  if (pthread_mutex_init(&thread_mutex, NULL)) {
    fprintf(stderr, "Unable to initialize the thread mutex\n");
    threadpool_destroy(thread_pool);
    free(dmag);
    exit(EXIT_FAILURE);
  }
  if (sem_init(&thread_sem, 0, nthread)) {
    fprintf(stderr, "Unable to initialize the thread semaphore\n");
    threadpool_destroy(thread_pool);
    pthread_mutex_destroy(&thread_mutex);
    free(dmag);
    exit(EXIT_FAILURE);
  }

  // Build the HTM data base
  if (verbose) {
    fprintf(stderr, "Building the HTM database...\n");
    progress_set_mode(PROGRESS_MODE_ADD);
  }
  db = cluster_parse_file(stdin);
  if (!db) {
    fprintf(stderr, "Unable to build the HTM database\n");
    free(dmag);
    threadpool_destroy(thread_pool);
    pthread_mutex_destroy(&thread_mutex);
    sem_destroy(&thread_sem);
    exit(EXIT_FAILURE);
  }

  if (verbose)
    progress_set_mode(PROGRESS_MODE_NONE);

  // Print the header line if asked to do so
  if (header) {
    if (rdensity > 0)
      fprintf(stdout, "%-19s %-4s %-12s input\n", "candidate_id", "nimg", "density");
    else
      fprintf(stdout, "%-19s %-4s input\n", "candidate_id", "nimg");
  }

  if (verbose) {
    fprintf(stderr, "\nBrowsing the HTM database...\n");
    progress_set_mode(PROGRESS_MODE_BROWSE);
  }

  { // Browse the whole data base
    const double p[] = { 1., 0., 0. };
    if (htmdb_constraint(db, p, -1., cluster_browse_db, db)) {
      fprintf(stderr, "Unable to browse the entire data base\n");
      htmdb_destroy(db, cluster_free_data);
      free(dmag);
      threadpool_destroy(thread_pool);
      pthread_mutex_destroy(&thread_mutex);
      sem_destroy(&thread_sem);
      exit(EXIT_FAILURE);
    }
    threadpool_wait(thread_pool);
  }

  if (verbose) {
    progress_set_mode(PROGRESS_MODE_NONE);
    fprintf(stderr, "\n");
  }

  // Destroy multi-threading variables
  threadpool_destroy(thread_pool);
  pthread_mutex_destroy(&thread_mutex);
  sem_destroy(&thread_sem);

  // In verbose mode, terminate progression monitoring
  if (verbose)
    progress_stop();

  // Destroy the database
  htmdb_destroy(db, cluster_free_data);
  free(dmag);

  // Exit with success
  exit(EXIT_SUCCESS);
}
