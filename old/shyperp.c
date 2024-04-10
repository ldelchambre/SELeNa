#include <fload.h>
#include <float.h>
#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <pthread.h>
#include <rfprintf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <utils.h>

#define HYPERPLANE_STR_LENGTH 256

/*************************************************************************************************
 * Typdefs
 ************************************************************************************************/
/**
 * Structure designed to contain the best splitting hyperplane we found so far.
 *
 *  The best splitting hyperplane is defined by two point (ip and in) and a
 * distance from ip, d.
 *
 *  Associated with this hyperplane, we have a score measure that should be
 *  maximized according to a user-defined function.
 */
typedef struct shyperp_solution {
  // Index of the positive instance defining the hyperplane
  size_t ip;

  // Index of the negative instance defining the hyperplane
  size_t in;

  // The distance from ip
  double d;

  // The score measure associated with this solution
  double score;
} shyperp_solution_t;

/**
 * Initializer for the shyperp_solution_t structure
 */
#define SHYPERP_SOLUTION_INIT {0, 0, 0., -DBL_MAX}

/**
 *  The function to use in order to define a score
 *
 *  Example:
 *    double my_score(double *bestd, double *pd, double *nd)
 */
typedef double (*score_function_t)(double *, double *, double *);

/*************************************************************************************************
 * Functions prototypes
 ************************************************************************************************/
double default_score_fct(double *d, double *Dp, double *Dn);
double mean_score_fct(double *d, double *Dp, double *Dn);

/*************************************************************************************************
 * Program parameters
 ************************************************************************************************/
#define DEFAULT_NTHREAD 1
#define DEFAULT_UPDATE_INTERVAL 2 /* s */

struct {
  char *name;
  score_function_t fct;
  char *description;
} score_measure[] =
    { { "default", default_score_fct,
        "The default score function (see \"exhaustive\")" },
        { "exhaustive", default_score_fct,
            "Search in all the projected distances for the threshold maximizing the difference between TPR and FPR" },
        { "mean", mean_score_fct,
            "Compute the difference between TPR and FPR based on the threshold coming from the mid-point between the mean distances of the positive and negative instances" },
        { NULL, NULL, NULL } };

unsigned int nthread = DEFAULT_NTHREAD;
unsigned int update_interval = DEFAULT_UPDATE_INTERVAL;

#define PARAM_NTHREAD   65538
#define PARAM_UPDATE_INTERVAL 65539

static struct option options_long[] = {
/*-s*/{ "score", required_argument, 0, 's' },
/*-l*/{ "list", no_argument, 0, 'l' },
/*  */{ "nthread", required_argument, 0, PARAM_NTHREAD },
/*  */{ "update-interval", required_argument, 0, PARAM_UPDATE_INTERVAL },
/*-?*/{ "help", no_argument, 0, '?' },
/*  */{ 0, 0, 0, 0 } };

static const char options_short[] = "s:l?";

/*************************************************************************************************
 * Help message
 ************************************************************************************************/
static const char *help_message =
    "USAGE:\n"
        "  shyperp [OPTIONS] LS\n"
        "\n"
        "DESCRIPTION:\n"
        "  Find a nearly optimal splitting hyperplane based on LS.\n"
        "\n"
        "BUILD_OPTIONS:\n"
        "  -s, --score,     Specify the score measure to use (see -l option for a list)\n"
        "  -l, --list,      List the score measures available\n"
        "      --nthread,   The number of thread to use [default=%u]\n"
        "      --update-interval, The update interval for printing the progress report [default=%u]\n"
        "  -?, --help,      Print this help message\n"
        "\n"
        "FILE FORMAT:\n"
        "  Each line of LS is assumed to be in the form \"y x0 x1 x2 x3 ...\n"
        "    where\n"
        "      y  is either 0 or 1\n"
        "      xi are double values of attributes.\n";

void print_help(FILE *out) {
  fprintf(out, help_message, DEFAULT_NTHREAD, DEFAULT_UPDATE_INTERVAL);
}

/*************************************************************************************************
 * Global variables
 ************************************************************************************************/

// The matrix of input attributes. First np columns are positive instances
double *X = NULL;

// The number of instances
size_t ninst = 0;

// The number of attributes
size_t nattr = 0;

// The number of positive instances (having y[i] = 1)
size_t np = 0;

// The number of negative instances (having y[i] = 0)
size_t nn = 0;

// Global lock
pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;

// The index of the next positive instance to consider
size_t nip = 0;

// The index of the next negative instance to consider
size_t nin = 0;

// The best solution found so far
shyperp_solution_t solution = SHYPERP_SOLUTION_INIT;

// The score function to use in order to compute solution.score and solution.d
score_function_t score_fct = NULL;

/*************************************************************************************************
 * File parsing functions
 ************************************************************************************************/
int load_learning_set(double **X, char **y, size_t *nobs, size_t *nattr,
    const char *filename) {
  FILE *in;
  fload_data_t *input_data;

  // Open the input file in read mode
  in = fopen(filename, "r");
  if (!in) {
    fprintf(stderr, "Unable to open input file %s (%s)\n", filename,
        strerror(errno));
    return -1;
  }

  // Load the input file
  input_data = fload(in, ' ');
  if (input_data == NULL) {
    fprintf(stderr, "Unable to load input file %s\n", filename);
    fclose(in);
    return -1;
  }

  // Close the input file
  fclose(in);

  // Set the number of observations/attributes
  *nobs = fload_nobs(input_data);
  *nattr = fload_nattr(input_data) - 1;

  // Allocate the X matrix and y array
  *X = (double *) calloc(*nobs * *nattr, sizeof(double));
  *y = (char *) calloc(*nobs, sizeof(char));
  if (*X == NULL || *y == NULL) {
    fload_destroy(input_data);
    free(*X);
    free(*y);
    return -1;
  }

  { // Convert the fload_data_t structure into an array
    double (*_X)[*nattr] = (double (*)[*nattr]) *X;
    char *_y = *y;
    size_t i, j;
    const char *s;
    char *endptr;

    // Copy the outcome value into y
    for (i = 0; i < *nobs; i++) {
      s = fload_get(input_data, i, 0);
      if (*s < '0' || '1' < *s || *(s + 1) != '\0') {
        fprintf(stderr,
            "Unable to parse \"%s\" as a boolean (outcome value at line %zu)\n",
            s, i + 1);
        fload_destroy(input_data);
        free(*X);
        free(*y);
        return -1;
      }
      _y[i] = *s - '0';
    }

    // Copy the input attributes into X
    for (i = 0; i < *nobs; i++) {
      for (j = 0; j < *nattr; j++) {
        s = fload_get(input_data, i, j + 1);
        _X[i][j] = strtod(s, &endptr);
        if (s == endptr || *endptr != '\0') {
          fprintf(stderr,
              "Unable to parse \"%s\" as a double (%zuth attribute of line %zu)\n",
              s, j + 1, i + 1);
          fload_destroy(input_data);
          free(*X);
          free(*y);
          return -1;
        }
      }
    }
  }

  // Destroy the input data
  fload_destroy(input_data);

  return 0;
}

/*************************************************************************************************
 * Score functions
 ************************************************************************************************/
double score_function_gain(const size_t ip, const size_t in) {
  const size_t n = ip + in;
  const double ptp = (double) ip / n;
  const double pfp = (double) in / n;
  const double ptn = (double) (nn - in) / (ninst - n);
  const double pfn = (double) (nn - ip) / (ninst - n);
  return ( ptp * log(ptp) + pfp * log(pfp) ) * n / ninst + ( ptn * log(ptn) + pfn * log(pfn) ) * (1. - n / ninst);
}

double default_score_fct(double *d, double *Dp, double *Dn) {
  // Score and threshold on the distance
  double sres = -DBL_MAX;
  double dres;

  // Temporary score and threshold on the distance
  double stmp;
  double dtmp;

  // Indices of the positive and negative instances
  size_t ip = 0;
  size_t in = 0;

  // Sort the arrays of distances
  qsort(Dp, np, sizeof(double), dbl_compar);
  qsort(Dn, nn, sizeof(double), dbl_compar);

  // Find a threshold on the distance such that the difference between the
  // true positive rate (TPR) and the false positive rate (FPR) is maximized
  while (ip < np && in < nn) {
    // The current threshold on the distance
    dtmp = 0.5 * (Dp[ip] + Dn[in]);

    // Update the number of positive/negative instances according
    // to the current threshold
    if (Dp[ip] < Dn[in])
      ip++;
    else
      in++;

    // Compute the difference between TPR and FPR
    // stmp = (double) ip / np - (double) in / nn;
    stmp = score_function_gain(ip, in);

    // Update score
    if (sres < stmp) {
      sres = stmp;
      dres = dtmp;
    }
  }

  // Return results
  *d = dres;
  return sres;
}

double mean_score_fct(double *d, double *Dp, double *Dn) {
  double meanp, meann;
  size_t i, ip, in;

  // Compute the mean distance of the positive instances
  meanp = 0;
  for (i = 0; i < np; i++)
    meanp += Dp[i];
  meanp /= np;

  // Compute the mean distance of the positive instances
  meann = 0;
  for (i = 0; i < nn; i++)
    meann += Dn[i];
  meann /= nn;

  // Compute the threshold on the distance
  *d = 0.5 * (meanp + meann);

  // Count the number of true positive instances
  for (i = ip = 0; i < np; i++)
    ip += (Dp[i] < *d);

  // Count the number of false positive instances
  for (i = in = 0; i < nn; i++)
    in += (Dn[i] < *d);

  return (double) ip / np - (double) in / nn;
}

/*************************************************************************************************
 * Build functions
 ************************************************************************************************/
void shyperp_build_thread_cleanup(void *arg) {
  free(*((void **) arg));
}

void *shyperp_build_thread(void *_arg) {
  // The set of positive instances
  double (*P)[nattr] = (double (*)[nattr]) X;

  // The set of negative instances
  double (*N)[nattr] = (double (*)[nattr]) (X + np * nattr);

  // Distances between the hyper-plane and all observations
  double *D = NULL;

  // The current projection vector
  double v[nattr];

  // Pointer to the positive and negative instances
  size_t ip, in;

  // The best distance and score associated with the combination of (ip, in)
  // (computed based on score_fct)
  double d, score;

  // General usage counter
  size_t i, j;

  // Add the cleanup procedure
  pthread_cleanup_push(shyperp_build_thread_cleanup, &D);

  // Allocate the distances array
  D = (double *) calloc(ninst, sizeof(double));
  if (D == NULL) {
    free(D);
    D = NULL;
    pthread_exit(NULL);
  }

  // Main thread loop, browse all positive instances
  pthread_mutex_lock(&lock);
  while (nip < np) {
    // Use the hyper-planes perpendicular to the line joining
    // P[ip] and N[in] while containing P[ip]
    ip = nip;
    in = nin++;
    if (in == nn) {
      nip++;
      nin = 0;
    }
    pthread_mutex_unlock(&lock);

    // Compute the current projection vector
    for (i = 0; i < nattr; i++)
      v[i] = N[in][i] - P[ip][i];

    { // Compute the distance between all X[i] and the hyper-plane
      double (*_X)[nattr] = P;

      for (i = 0; i < ninst; i++) {
        D[i] = 0;
        for (j = 0; j < nattr; j++)
          D[i] += v[j] * (_X[i][j] - P[ip][j]);
      }
    }

    // Compute the score measure
    score = score_fct(&d, D, D + np);

    // Update the best solution find so far
    pthread_mutex_lock(&lock);
    if (solution.score < score) {
      solution.score = score;
      solution.d = d;
      solution.ip = ip;
      solution.in = in;
    }
  }
  pthread_mutex_unlock(&lock);

  // Execute the cleanup procedure
  pthread_cleanup_pop(1);

  // Exit thread
  pthread_exit(NULL);
}

void shyperp_hyperplane_str(char *s, const size_t n) {
  const double (*P)[nattr] = (const double (*)[nattr]) X;
  const double (*N)[nattr] = (const double (*)[nattr]) (X + np * nattr);
  const size_t ip = solution.ip;
  const size_t in = solution.in;
  double b = solution.d;
  double v;
  size_t i, l;

  // Build the formula of the hyperplane
  for (i = l = 0; i < nattr; i++) {
    v = N[in][i] - P[ip][i];
    b += v * P[ip][i];
    if (l < n)
      l += snprintf(&s[l], n - l, "%+g*x%zu", v, i);
  }

  if (l < n)
    snprintf(&s[l], n - l, "<%+g", b);
}

int shyperp_build(const char *filename) {
  pthread_t threads[nthread];
  char formula[HYPERPLANE_STR_LENGTH];
  char *y;
  size_t i, j, k;

  // Load the learning set of observations
  if (load_learning_set(&X, &y, &ninst, &nattr, filename)) {
    fprintf(stderr,
        "Unable to load the learning set of observation from file %s\n",
        filename);
    return -1;
  }

  // Set all positive instance in the first np columns of X
  {
    double (*_X)[nattr] = (double (*)[nattr]) X;
    double tmp;

    // Loop invariant: y[k] == 1 for all k < i
    //                 y[k] == 0 for all j < k
    i = 0;
    j = ninst - 1;
    while (i <= j) {
      // Skip all positive instances
      while (y[i] == 1 && i <= j)
        i++;
      // Skip all negative instances
      while (y[j] == 0 && i <= j)
        j--;
      // Swap X[i][:] with X[j][:]
      if (i <= j) {
        for (k = 0; k < nattr; k++) {
          tmp = _X[i][k];
          _X[i][k] = _X[j][k];
          _X[j][k] = tmp;
        }
        // y[i] is now positive and y[j] negative
        y[i] = 1;
        y[j] = 0;
      }
    }

    // Set the number of positive and negative instances
    np = i;
    nn = ninst - np;
  }

  // We don't need y anymore
  free(y);

  // Check that 2 classes are present in the dataset
  if (np == 0 || nn == 0) {
    fprintf(stderr,
        "Only one outcome value is present, nothing has to be done\n");
    free(X);
    return 0;
  }

  // Launch all threads
  for (i = 0; i < nthread; i++) {
    if (pthread_create(&threads[i], NULL, shyperp_build_thread, NULL)) {
      fprintf(stderr, "Unable to create thread %zu\n", i);
      for (j = 0; j < i; j++)
        pthread_cancel(threads[j]);
      free(X);
      return -1;
    }
  }

  // Print progress report
  while (nip < np) {
    pthread_mutex_lock(&lock);
    shyperp_hyperplane_str(formula, HYPERPLANE_STR_LENGTH);
    rfprintf(stderr, "Finished: %f%%, Score: %g, Hyperplane: %s",
        100. * (nip * nn + nin) / (np * nn), solution.score, formula);
    pthread_mutex_unlock(&lock);
    sleep(update_interval);
  }

  // Wait for all threads to finish
  for (i = 0; i < nthread; i++) {
    if (pthread_join(threads[i], NULL)) {
      fprintf(stderr, "Unable to join thread number %zu\n", i);
      for (j = i + 1; j <= nthread; j++)
        pthread_cancel(threads[j]);
      free(X);
      return -1;
    }
  }

  shyperp_hyperplane_str(formula, HYPERPLANE_STR_LENGTH);
  rfprintf(stderr, "Finished: %f%%, Score: %g, Hyperplane: %s\n",
      100. * (nip * nn + nin) / (np * nn), solution.score, formula);

  return 0;
}

void shyperp_result(void) {
  // The set of positive instances
  double (*P)[nattr] = (double (*)[nattr]) X;

  // The set of negative instances
  double (*N)[nattr] = (double (*)[nattr]) (X + np * nattr);

  // String containing the formula describing the hyperplane
  char formula[HYPERPLANE_STR_LENGTH];

  size_t ip, in, i, j;
  double tmp;

  // Compute the number of positive instances standing on the negative part
  // of the hyperplane (i.e. the true positive)
  ip = 0;
  for (i = 0; i < np; i++) {
    tmp = 0;
    for (j = 0; j < nattr; j++)
      tmp += (N[solution.in][j] - P[solution.ip][j])
          * (P[i][j] - P[solution.ip][j]);
    ip += (tmp < solution.d);
  }

  // Compute the number of negative instances standing on the negative part
  // of the hyperplane (i.e. the false positive)
  in = 0;
  for (i = 0; i < nn; i++) {
    tmp = 0;
    for (j = 0; j < nattr; j++)
      tmp += (N[solution.in][j] - P[solution.ip][j])
          * (N[i][j] - P[solution.ip][j]);
    in += (tmp < solution.d);
  }

  // Build the hyperplane formula
  shyperp_hyperplane_str(formula, HYPERPLANE_STR_LENGTH);

  // Print the hyperplane formula + confusion matrix
  fprintf(stdout, "%s\n", formula);
  fprintf(stderr, "         |  Positive |  Negative | <= Predicted class\n");
  fprintf(stderr, "---------+-----------+-----------|\n");
  fprintf(stderr, "Positive | %9zu | %9zu | TPR = %f\n", ip, np - ip,
      (double) ip / np);
  fprintf(stderr, "Negative | %9zu | %9zu | FPR = %f\n", in, nn - in,
      (double) in / nn);
  fprintf(stderr, "---------+-----------+-----------|\n");
}

/*************************************************************************************************
 * Main program
 ************************************************************************************************/
int main(int argc, char **argv) {
  const char *filename;
  int opt, i;

  // Get the program parameters
  while ((opt = getopt_long(argc, argv, options_short, options_long, NULL))
      != -1) {
    switch (opt) {
    case 's':
      for (i = 0;
          strcmp(score_measure[i].name, optarg) != 0
              && score_measure[i].fct != NULL; i++)
        ;
      if (score_measure[i].fct != NULL)
        score_fct = score_measure[i].fct;
      else
        fprintf(stderr,
            "Unrecognized score function \"%s\", using the default one\n",
            optarg);
      break;
    case 'l':
      for (i = 0; score_measure[i].fct != NULL; i++)
        fprintf(stdout, "%s : %s\n", score_measure[i].name, score_measure[i].description);
      exit(EXIT_SUCCESS);
      break;
    case PARAM_NTHREAD:
      nthread = strtoul(optarg, NULL, 10);
      break;
    case PARAM_UPDATE_INTERVAL:
      update_interval = strtoul(optarg, NULL, 10);
      break;
    case '?':
      print_help(stdout);
      exit(EXIT_SUCCESS);
      break;
    default:
      print_help(stderr);
      exit(EXIT_FAILURE);
      break;
    }
  }

  // Get the input file name
  if (optind < argc)
    filename = argv[optind++];
  else {
    fprintf(stderr, "No input file found\n");
    print_help(stderr);
    exit(EXIT_FAILURE);
  }

  // Check options
  if (score_fct == NULL)
    score_fct = default_score_fct;

  // Build the splitting hyperplane
  if (shyperp_build(filename)) {
    fprintf(stderr, "An error occurs while building the hyperplane\n");
    exit(EXIT_FAILURE);
  }

  // Print result
  shyperp_result();

  // Exit
  exit(EXIT_SUCCESS);
}
