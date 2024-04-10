#include <astrostr.h>
#include <fitsio.h>
#include <getopt.h>
#include <htmdb.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/**
 * Limit to set on the proper motions and parallaxes of the objects so as to
 * consider them as non moving [mas]
 */
#define _PMP_EPSILON 1e-10

#define ARCMIN (M_PI/10800.)
#define ARCSEC (M_PI/648000.)

#ifdef __GNUC__
#define UNUSED(param) __attribute__((unused)) param
#else
#define UNUSED(param) param
#endif

/*******************************************************************************
 * Program parameters
 ******************************************************************************/
#define DEFAULT_SEP (5.*ARCSEC)
#define DEFAULT_DMAG 2.5 // [mag]
#define DEFAULT_LEVEL 10
#define DEFAULT_NMIN 3
#define DEFAULT_NMAX -1U
#define DEFAULT_ABSBMIN 0.
#define DEFAULT_RDENSITY (0.5*ARCMIN)
#define DEFAULT_RAMIN 0.
#define DEFAULT_RAMAX (2.*M_PI)
#define DEFAULT_DECMIN (-0.5*M_PI)
#define DEFAULT_DECMAX (0.5*M_PI)

double sep = DEFAULT_SEP;
double sep_cos /* = cos(sep)*/;
double dmag = DEFAULT_DMAG;
unsigned int level = DEFAULT_LEVEL;
unsigned int nmin = DEFAULT_NMIN;
unsigned int nmax = DEFAULT_NMAX;
double absb = DEFAULT_ABSBMIN;
double rdensity = DEFAULT_RDENSITY;
double rdensity_cos /* = cos(rdensity) */;
double rdensity_area /* = 4*M_PI*sin(rdensity/2)^2 */;
double ramin = DEFAULT_RAMIN;
double ramax = DEFAULT_RAMAX;
double decmin = DEFAULT_DECMIN;
double decmax = DEFAULT_DECMAX;
char header = 1;
char verbose = 0;

double raminext /* = ramin - max(sep,rdensity) */;
double ramaxext /* = ramax + max(sep,rdensity) */;
double decminext /* = decmin - max(sep,rdensity) */;
double decmaxext /* = decmax + max(sep,rdensity) */;
double absbext /* = absb - max(sep,rdensity) */;

#define RAMIN_PARAM    65536
#define RAMAX_PARAM    65537
#define DECMIN_PARAM   65538
#define DECMAX_PARAM   65539
#define NOHEADER_PARAM 65540

static struct option options_long[] = {
/*-s*/{ "sep", required_argument, 0, 's' },
/*-d*/{ "dmag", required_argument, 0, 'd' },
/*-l*/{ "level", required_argument, 0, 'l' },
/*-n*/{ "nmin", required_argument, 0, 'n' },
/*-N*/{ "nmax", required_argument, 0, 'N' },
/*-b*/{ "absb", required_argument, 0, 'b' },
/*-D*/{ "rdensity", required_argument, 0, 'D' },
/*  */{ "ramin", required_argument, 0, RAMIN_PARAM },
/*  */{ "ramax", required_argument, 0, RAMAX_PARAM },
/*  */{ "decmin", required_argument, 0, DECMIN_PARAM },
/*  */{ "decmax", required_argument, 0, DECMAX_PARAM },
/*  */{ "noheader", no_argument, 0, NOHEADER_PARAM },
/*-v*/{ "verbose", no_argument, 0, 'v' },
/*-?*/{ "help", no_argument, 0, '?' },
/*  */{ 0, 0, 0, 0 } };

static const char options_short[] = "s:d:l:n:N:b:D:v?";

/*************************************************************************************************
 * Help message
 ************************************************************************************************/
static const char help_message[] =
  "USAGE:\n"
  "  candidates [options] file\n"
  "\n"
  "DESCRIPTION:\n"
  "  Search for clusters of celestial objects in catalogue.\n"
  "\n"
  "OPTIONS:\n"
  "  -s, --sep,       The maximal great circle distance between any two objects of the cluster [default = %garcsec]\n"
  "  -d, --dmag,      The maximal magnitude difference between any two objects of the cluster [default = %gmag]\n"
  "  -l, --level,     The HTM recursion level to use in order to build database [default = %u]\n"
  "  -n, --nmin,      The minimal number of objects in the cluster [default = %u]\n"
  "  -N, --nmax,      The maximal number of objects in the cluster [default = %u]\n"
  "  -b, --absb,      The minimal absolute galactic latitude to consider [default = %g]\n"
  "  -D, --rdensity,  The radius to use in order to estimate the density of objects in the vicinity of the candidate [default = %garcmin]\n"
  "      --ramin,     The minimal right ascension to consider [default = %g]\n"
  "      --ramax,     The maximal right ascension to consider [default = %g]\n"
  "      --decmin,    The minimal declination to consider [default = %g]\n"
  "      --decmax,    The minimal declination to consider [default = %g]\n"
  "      --noheader,  Don't print header line\n"
  "  -v, --verbose,   Should we have to be verbose\n"
  "  -?, --help,      Print this help message\n"
  "\n"
  "FILE FORMAT:\n"
  "  file must be a FITS file containing the following fields\n"
  "    'source_id',       The integer source identifier\n"
  "    'ra',              The right ascension [degree]\n"
  "    'dec',             The declination [degree]\n"
  "    'parallax',        The parallax of the source [mas]\n"
  "    'pmra',            The proper motion in right ascension [mas/year]\n"
  "    'pmdec',           The proper motion in declination [mas/year]\n"
  "    'phot_g_mean_mag', The G-band magnitude of the source [mag]\n"
  "    'b',               The galactic latitude of the source [degree]\n"
  "\n"
  "ANGLE FORMAT:\n"
  "  Angles are expressed in radian by default. Other possibilities are degree,\n"
  "  arcmin, arcsec, mas, hour, minute or second (e.g. \"30arcsec\" is a valid input)\n"
  "  Shorthand to these units are d, ', \", h, m, s respectively for degree, \n"
  "  arcmin, arcsec, hour, minute and second (e.g. \"30d\").\n"
  "  Right ascension and declination can be expressed in sexagecimal notation\n"
  "  respectively through the \"hh:mm:ss\" and \"dd:mm:ss\" notations (e.g. \"10:45:57.14\")\n";

void printHelp(FILE *out) {
  fprintf(out, help_message, DEFAULT_SEP/ARCSEC, DEFAULT_DMAG, DEFAULT_LEVEL,
                             DEFAULT_NMIN, DEFAULT_NMAX, DEFAULT_ABSBMIN,
                             DEFAULT_RDENSITY/ARCMIN, DEFAULT_RAMIN,
                             DEFAULT_RAMAX, DEFAULT_DECMIN,
                             DEFAULT_DECMAX);
}

/*******************************************************************************
 * Utility functions
 ******************************************************************************/
void icrs2gal(double *vgal, const double *vicrs) {
  // Matrix taken from the SOFA release 2018-01-30, icrs2g function
  // Copyright (C) 2018 IAU SOFA Board.
  const double r[3][3] = {
    {-0.054875560416215368492398900454,-0.873437090234885048760383168409,-0.483835015548713226831774175116 },
    {+0.494109427875583673525222371358,-0.444829629960011178146614061616,+0.746982244497218890527388004556 },
    {-0.867666149019004701181616534570,-0.198076373431201528180486091412,+0.455983776175066922272100478348 }
  };
  unsigned char i, j;

  for (i = 0; i < 3; i++) {
    vgal[i] = 0.;
    for (j = 0; j < 3; j++)
      vgal[i] += r[i][j] * vicrs[j];
  }
}
/*******************************************************************************
 * Type defintions
 ******************************************************************************/
typedef struct gaia_data {
  long long source_id;
  double gmag;
} gaia_data_t;

/*******************************************************************************
 * Column to export out of the Gaia catalogue
 ******************************************************************************/
#define NEXPORT 8

typedef struct export_data {
  long long source_id;
  double ra;
  double dec;
  double parallax;
  double pmra;
  double pmdec;
  double phot_g_mean_mag;
  double b;
} export_data_t;

static const struct {
  char *name;
  int type; // TBIT,TBYTE,TSBYTE,TLOGICAL,STRING,TUSHORT,TSHORT,TUINT,TINT,
            // TULONG,TLONG,TFLOAT,TLONGLONG,TDOUBLE,TCOMPLEX,TDBLCOMPLEX
} export[NEXPORT] = {
  {"source_id", TLONGLONG},
  {"ra", TDOUBLE},
  {"dec", TDOUBLE},
  {"parallax", TDOUBLE},
  {"pmra", TDOUBLE},
  {"pmdec", TDOUBLE},
  {"phot_g_mean_mag", TDOUBLE},
  {"b", TDOUBLE}
};

void fill_export_data(export_data_t *data, const iteratorCol *column_infos, const long irow) {
  data->source_id =       ((const LONGLONG *) column_infos[0].array)[irow];
  data->ra =              ((const double *)   column_infos[1].array)[irow] * M_PI/180.;
  data->dec =             ((const double *)   column_infos[2].array)[irow] * M_PI/180.;
  data->parallax =        ((const double *)   column_infos[3].array)[irow];
  data->pmra =            ((const double *)   column_infos[4].array)[irow];
  data->pmdec =           ((const double *)   column_infos[5].array)[irow];
  data->phot_g_mean_mag = ((const double *)   column_infos[6].array)[irow];
  data->b               = ((const double *)   column_infos[7].array)[irow] * M_PI/180.;
}

/*******************************************************************************
 * Progression report functions
 ******************************************************************************/

/**
 * Structure used for printing progression report in verbose mode
 */
struct {
  time_t t0;
  LONGLONG n; ///< The current number of object
  LONGLONG ncluster; ///< The number of cluster we found
  LONGLONG ndb; ///< The total number of objects in the HTM database
  LONGLONG ntot; ///< The total number of objects in the input file
} progress;

void print_progress(FILE *out, const LONGLONG n, const LONGLONG ntot) {
  const double p0 = floor(20. * n / ntot);
  const double p1 = floor(20. * (1.+n) / ntot);
  if(p0 != p1) fprintf(out, "%3.0f%%, run time: %.0fs\n", 5*p1, difftime(time(NULL),progress.t0));
}

/*******************************************************************************
 * HTM database construction functions
 ******************************************************************************/
int fill_db(long UNUSED(nperiter), long UNUSED(offset), long UNUSED(irow), long nrow,
            int UNUSED(ncol), iteratorCol *column_infos, void *arg) {
  htmdb_t *db = (htmdb_t *) arg;
  export_data_t exp;
  gaia_data_t *data;
  double v[3];
  long i;

  // Browse each row of this iteration
  for(i = 0; i < nrow; i++) {
    // Fill the exported data
    fill_export_data(&exp, column_infos, i+1);

    // If this object is not moving (i.e. no parallax, no proper motions)
    // and corresponds to our selection criterion
    if(fabs(exp.parallax) < _PMP_EPSILON
       && fabs(exp.pmra) < _PMP_EPSILON
       && fabs(exp.pmdec) < _PMP_EPSILON
       && raminext   <= exp.ra
       && exp.ra     <= ramaxext
       && decminext  <= exp.dec
       && exp.dec    <= decmaxext
       && absbext    <= fabs(exp.b)) {

      // Get the cartesian coordinates of this point
      htm_cart_coord(v, exp.ra, exp.dec);

      // Allocate the data that are associated with this point
      data = malloc(sizeof(gaia_data_t));
      if(!data) {
        fprintf(stderr, "Unable to allocate source informations (source_id=%lld)\n", exp.source_id);
        return 1001; // FITSIO uses return value below 1000 for personal purpose
      }

      // Set data associated to this point
      data->source_id = exp.source_id;
      data->gmag = exp.phot_g_mean_mag;

      // Add this object to the HTM data base
      if(htmdb_add(db, v, data)) {
        fprintf(stderr, "Unable to add observation into the HTM database (source_id=%lld)\n", exp.source_id);
        return 1002;
      }

      // We added one more object in the database
      progress.ndb++;
    }

    // In verbose mode, print progress
    if(verbose) print_progress(stderr, ++progress.n, progress.ntot);
  }

  return 0;
}

htmdb_t *read_catalog(const char *filename) {
  htmdb_t *db;
  iteratorCol column_infos[NEXPORT];
  int status = 0;
  fitsfile *fptr;
  unsigned int i;

  // Initialize the HTM data base
  db = htmdb_init(level);
  if(!db) {
    fprintf(stderr, "Unable to allocate the HTM data base (level=%u)\n", level);
    return NULL;
  }

  // Open the FITS table
  if(fits_open_table(&fptr, filename, READONLY, &status)) {
    fits_report_error(stderr, status);
    htmdb_destroy(db, free);
    return NULL;
  }

  // In verbose mode, get the total number of rows
  if(verbose) {
    progress.n = progress.ndb = 0;
    if(fits_get_num_rowsll(fptr, &progress.ntot, &status)) {
      fits_report_error(stderr, status);
      fits_close_file(fptr, &status);
      htmdb_destroy(db, free);
      return NULL;
    }
    fprintf(stderr, "%lld objects found in file %s\n", progress.ntot, filename);
  }

  // Fill columns information
  for(i = 0; i < NEXPORT; i++) {
    if(fits_iter_set_by_name(&column_infos[i], fptr, export[i].name, export[i].type, InputCol)) {
      fits_report_error(stderr, status);
      fits_close_file(fptr, &status);
      htmdb_destroy(db, free);
      return NULL;
    }
  }

  // Iterate over the data
  if (fits_iterate_data(NEXPORT, column_infos, 0, 0, fill_db, db, &status) > 0) {
    fits_report_error(stderr, status);
    fits_close_file(fptr, &status);
    htmdb_destroy(db, free);
    return NULL;
  }

  // Close the FITS file
  fits_close_file(fptr, &status);

  return db;
}

/*******************************************************************************
 * Cluster identification functions
 ******************************************************************************/
typedef struct cluster_obj {
  const double *v;
  const gaia_data_t *data;
  struct cluster_obj *next;
} cluster_obj_t;

typedef struct cluster_with_constraint {
  cluster_obj_t **cluster;
  const gaia_data_t *constraints;
} cluster_with_constraint_t;

void free_cluster(cluster_obj_t *cluster) {
  cluster_obj_t *next;

  while(cluster) {
    next = cluster->next;
    free(cluster);
    cluster = next;
  }
}

int get_cluster(const double *v, void *_data, void *arg) {
  // Get the current cluster and associated constraints
  cluster_with_constraint_t *clusterc = (cluster_with_constraint_t *) arg;
  // Get the data
  const gaia_data_t *data = (const gaia_data_t *) _data;

  // Only add this object if it has a source_id that is greater than the one
  // around which the constraint was performed and if the absolute difference
  // in magnitude is less than dmag
  if(clusterc->constraints->source_id < data->source_id
     && fabs(clusterc->constraints->gmag-data->gmag) < dmag) {
    // Get the list of already found objects
    cluster_obj_t **cluster = clusterc->cluster;

    // Allocate the new object
    cluster_obj_t *obj = malloc(sizeof(cluster_obj_t));
    if(!obj) {
      fprintf(stderr, "Unable to allocate memory for cluster object (source_id=%lld)\n", data->source_id);
      free_cluster(*cluster); *cluster = NULL;
      return -1;
    }

    // Fill the object information
    obj->v = v;
    obj->data = data;

    // Add this object to the list
    obj->next = *cluster;
    *cluster = obj;
  }

  return 0;
}

char *distm_cluster(const cluster_obj_t *cluster, const unsigned long nobj) {
  char (*distm)[nobj] = (char (*)[nobj]) calloc(nobj*nobj,sizeof(char));
  unsigned long i,j;
  const cluster_obj_t *o1, *o2;

  if(!distm)
    return NULL;

  for(i = 0, o1 = cluster; o1; i++, o1 = o1->next) {
    distm[i][i] = 1; // Distance bewteen o1 and o1 is always fine
    for(j = 1, o2 = cluster->next; o2; j++, o2 = o2->next) {
      // Set the distance to 1 if angular distance match and that the difference
      // in magnitude do also match
      const double cost = o1->v[0]*o2->v[0]+o1->v[1]*o2->v[1]+o1->v[2]*o2->v[2];
      distm[i][j] = distm[j][i] = (cost >= sep_cos && fabs(o1->data->gmag-o2->data->gmag) <= dmag);
    }
  }

  return (char *) distm;
}

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

void get_identifier(char *s, const double ra, const double dec) {
  char stmp[14];

  // Build RA in the form "hh:mm:ss.sss"
  ra2str(stmp, ra);

  // Convert RA in the form "hhmmsssss"
  s[0]=stmp[0];  // h
  s[1]=stmp[1];  // h
  s[2]=stmp[3];  // m
  s[3]=stmp[4];  // m
  s[4]=stmp[6];  // s
  s[5]=stmp[7];  // s
  s[6]=stmp[9];  // s
  s[7]=stmp[10]; // s
  s[8]=stmp[11]; // s

  // Build DEC in the form "+dd:mm:ss.sss"
  dec2str(stmp, dec);

  s[ 9]=stmp[0]; // +
  s[10]=stmp[1]; // d
  s[11]=stmp[2]; // d
  s[12]=stmp[4]; // m
  s[13]=stmp[5]; // m
  s[14]=stmp[7]; // s
  s[15]=stmp[8]; // s
  s[16]=stmp[10]; // s
  s[17]=stmp[11]; // s
  s[18]=stmp[12]; // s

  s[19]='\0';
}

void print_cluster(const double *v, const gaia_data_t *data,
                   const cluster_obj_t *cluster, const unsigned long *icluster, const unsigned long ncluster,
                   const double density) {
  const cluster_obj_t *obj;
  double ra, dec;
  char candidate_id[20];
  unsigned long i,j;


  { // Get the candidate identifier
    double n, b, vmean[3];

    // Compute the mean position of the candidate
    vmean[0] = v[0]; vmean[1] = v[1]; vmean[2] = v[2];
    for(i = j = 0, obj = cluster; i < ncluster; j++, obj = obj->next) {
      if(icluster[i] == j) {
        vmean[0] += obj->v[0]; vmean[1] += obj->v[1]; vmean[2] += obj->v[2];
        i++;
      }
    }
    n = sqrt(vmean[0]*vmean[0]+vmean[1]*vmean[1]+vmean[2]*vmean[2]);
    vmean[0] /= n; vmean[1] /= n; vmean[2] /= n;

    // Convert cartesian coordinates into spherical coordinates
    htm_sphere_coord(&ra,&dec,vmean);

    // Compute the galactic coordinate of the mean position
    {
      double l, vmeang[3];
      icrs2gal(vmeang, vmean);
      htm_sphere_coord(&l, &b, vmeang);
    }

    // If this cluster is not contained in our region of interest, skip it
    if(ra < ramin || ramax < ra || dec < decmin || decmax < dec || fabs(b) < absb)
      return;

    // Get the candidate identifier
    get_identifier(candidate_id, ra, dec);
  }

  // Print the object around which the search was performed
  htm_sphere_coord(&ra,&dec,v);
  fprintf(stdout, "%19s %4lu %20lld %+18.15e %+18.15e %+18.15e %12.6e\n",
                  candidate_id, ncluster+1, data->source_id, ra, dec, data->gmag, density);

  // Browse all other objects from the cluster
  for(i = j = 0; i < ncluster; j++, cluster = cluster->next) {
    if(icluster[i] == j) {
      htm_sphere_coord(&ra,&dec,cluster->v);
      fprintf(stdout, "%19s %4lu %20lld %+18.15e %+18.15e %+18.15e %12.6e\n",
                      candidate_id, ncluster+1, cluster->data->source_id, ra, dec, cluster->data->gmag, density);
      i++;
    }
  }
}

int count_nearby_objects(const double UNUSED(*v), void UNUSED(*data), void *arg) {
  unsigned long *nnear = (unsigned long *) arg;
  (*nnear)++;
  return 0;
}

int validate_cluster(const htmdb_t *db, const double *v, const gaia_data_t *data, const cluster_obj_t *cluster) {
  const cluster_obj_t *obj;
  unsigned long nobj;
  char *distm;
  unsigned long *icluster;
  unsigned long iobj, ncluster;
  unsigned long nnear;
  double density;

  // Count the number of objects (not including the object data->source_id)
  for(obj = cluster, nobj = 0; obj; obj = obj->next, nobj++)
    ;

  // If not enough object are present in this cluster
  if(nobj + 1 < nmin)
    return 0;

  // Compute the distance matrix between all objects from the cluster
  distm = distm_cluster(cluster, nobj);
  if(!distm) {
    fprintf(stderr, "Unable to allocate distance matrix (source_id=%lld)\n", data->source_id);
    return -1;
  }

  // Compute the density of objects in the vicinity of this point
  nnear = 0;
  if(rdensity_cos < 1. && htmdb_constraint(db, v, rdensity_cos, count_nearby_objects, &nnear)) {
    fprintf(stderr, "Unable to estimate the density of objects (source_id=%lld)\n", data->source_id);
    free(distm);
    return -1;
  }
  density = (rdensity_cos < 1.) ? 1./rdensity_area*nnear : 0.;

  // Allocate the subcluster indices
  icluster = calloc(nobj,sizeof(unsigned long));
  if(!icluster) {
    fprintf(stderr, "Unable to allocate sub-cluster indices (source_id=%lld)\n", data->source_id);
    free(distm);
    return -1;
  }

  // Get all sub-clusters
  for(iobj = 0; iobj <= nobj+1-nmin; iobj++) {
    // Get the largest sub-cluster containing the iobj object
    ncluster = get_subcluster(icluster, distm, iobj, nobj);
    // If we have enough images in the sub-cluster, print the solution
    if(nmin <= ncluster + 1 && ncluster + 1 <= nmax) {
      print_cluster(v,data,cluster,icluster,ncluster,density);
      progress.ncluster++;
    }
  }

  // Free the distance matrix and the sub-cluster indices
  free(distm);
  free(icluster);

  return 0;
}

int browse_db(const double *v, void *_data, void *arg) {
  const htmdb_t *db = (const htmdb_t *) arg;
  const gaia_data_t *data = (const gaia_data_t *) _data;
  cluster_obj_t *cluster = NULL;

  // Perform a constraint around this specific object so as to get cluster with
  // the constraints that
  //   data->source_id < source_id (avoid to duplicate the identification of the clusters)
  // and that
  //   |data->gmag - gmag| <= dmag (only relevant objects are identified)
  {
    cluster_with_constraint_t clusterc;
    clusterc.cluster = &cluster;
    clusterc.constraints = data;
    if(htmdb_constraint(db,v,sep_cos,get_cluster,&clusterc)) {
      fprintf(stderr, "Unable to set constraint around object (source_id=%lld)\n", data->source_id);
      return -1;
    }
  }

  // Validate this cluster and print matching objects
  if(validate_cluster(db, v,data,cluster)) {
    fprintf(stderr, "Unable to validate cluster around object (source_id=%lld)\n", data->source_id);
    free_cluster(cluster);
    return -1;
  }

  // Delete this cluster
  free_cluster(cluster);

  // Print progress in verbose mode
  if(verbose) print_progress(stderr, ++progress.n, progress.ndb);

  return 0;
}

/*******************************************************************************
 * Main function
 ******************************************************************************/
int main(int argc, char **argv) {
  const char *filename;
  htmdb_t *db; // The HTM database to use
  int opt;

  // Get the program parameters
  while ((opt = getopt_long(argc, argv, options_short, options_long, NULL)) != -1) {
    switch (opt) {
      case 's':
        str2angle(&sep, optarg, NULL);
        break;
      case 'd':
        dmag = strtod(optarg, NULL);
        break;
      case 'l':
        level = strtoul(optarg, NULL, 10);
        break;
      case 'n':
        nmin = strtoul(optarg, NULL, 10);
        break;
      case 'N':
        nmax = strtoul(optarg, NULL, 10);
        break;
      case 'b':
        str2dec(&absb, optarg, NULL);
        break;
      case 'D':
        str2angle(&rdensity, optarg, NULL);
        break;
      case RAMIN_PARAM:
        str2ra(&ramin, optarg, NULL);
        break;
      case RAMAX_PARAM:
        str2ra(&ramax, optarg, NULL);
        break;
      case DECMIN_PARAM:
        str2ra(&decmin, optarg, NULL);
        break;
      case DECMAX_PARAM:
        str2ra(&decmax, optarg, NULL);
        break;
      case NOHEADER_PARAM:
        header = 0;
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

  // Compute the cosine of the search radius
  sep_cos = cos(sep);

  // Compute the cosine of the density radius and its area
  rdensity_cos = cos(rdensity);
  rdensity_area = 4.*M_PI*pow(sin(0.5*rdensity),2);

  // Adjust the ramin, ramax, decmin, decmax, absb parameters by taking into
  // account the extended region that is needed is order to compute the density
  // and/or to extract candidates having a mean position standing in the specified
  // range of positions.
  {
    const double ext = (sep < rdensity)    ? rdensity     : sep;
    raminext =  (ramin  - ext  > 0.)       ? ramin  - ext : 0.;
    ramaxext =  (ramax  + ext  < 2.*M_PI)  ? ramax  + ext : 2.*M_PI;
    decminext = (decmin - ext > -0.5*M_PI) ? decmin - ext : -0.5*M_PI;;
    decmaxext = (decmax + ext < 0.5*M_PI)  ? decmax + ext : 0.5*M_PI;
    absbext   = (absb   - ext > 0.)        ? absb   - ext : 0.;
  }

  // Get the input file name
  if(optind < argc)
    filename = argv[optind++];
  else {
    fprintf(stderr, "No input file found\n");
    printHelp(stderr);
    exit(EXIT_FAILURE);
  }

  // Check for extra arguments
  if(optind < argc) {
    fprintf(stderr, "Warning: extra argument(s) found (");
    while(optind < argc)
      fprintf(stderr, " %s", argv[optind++]);
    fprintf(stderr, ")\n");
  }

  // Print the header line if asked to do so
  if(header)
    fprintf(stdout, "%-19s %-4s %-20s %-22s %-22s %-22s %-12s\n",
            "candidate_id", "nimg", "source_id", "ra", "dec", "gmag", "density");

  // Build the HTM data base
  if(verbose) {
    fprintf(stderr, "Building the HTM database...\n");
    progress.t0 = time(NULL);
  }
  db = read_catalog(filename);
  if(!db) {
    fprintf(stderr, "Unable to build the HTM database\n");
    exit(EXIT_FAILURE);
  }
  if(verbose)
    fprintf(stderr, "%lld objects added to the database in %fs\n", progress.ndb, difftime(time(NULL), progress.t0));

  if(verbose) {
    fprintf(stderr, "Browsing the HTM database...\n");
    progress.t0 = time(NULL);
  }
  { // Browse the whole data base
    const double p[] = {1.,0.,0.};
    progress.n = progress.ncluster = 0;
    if(htmdb_constraint(db, p, -1., browse_db, db)) {
      fprintf(stderr, "Unable to browse the entire data base\n");
      htmdb_destroy(db, free);
      exit(EXIT_FAILURE);
    }
  }
  if(verbose)
    fprintf(stderr, "%lld clusters found in %fs\n", progress.ncluster, difftime(time(NULL), progress.t0));

  // Destroy the database
  htmdb_destroy(db, free);

  // Exit with success
  exit(EXIT_SUCCESS);
}
