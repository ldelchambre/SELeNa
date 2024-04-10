#include <extratree.h>
#include <ExtraTrees.h>
#include <errno.h>
#include <fload.h>
#include <fparse.h>
#include <getopt.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*************************************************************************************************
 * Program parameters
 ************************************************************************************************/
#define DEFAULT_MODE 'R'
#define DEFAULT_N 10
#define DEFAULT_K 0
#define DEFAULT_NMIN 0
#define DEFAULT_NTHREAD 1
#define DEFAULT_SEED 0

enum {
  RunModeNone, RunModeBuild, RunModePredict
} RunMode = RunModeNone;
char mode = DEFAULT_MODE;
unsigned int N = DEFAULT_N;
unsigned int k = DEFAULT_K;
unsigned int nmin = DEFAULT_NMIN;
unsigned int nthread = DEFAULT_NTHREAD;
unsigned int seed = DEFAULT_SEED;
char noheader = 0;
char headercol = 0;

#define PARAM_BUILD     65536
#define PARAM_PREDICT   65537
#define PARAM_NTHREAD   65538
#define PARAM_NOHEADER  65539
#define PARAM_HEADERCOL 65540

static struct option options_long[] = {
/*  */{ "build", no_argument, 0, PARAM_BUILD },
/*  */{ "predict", no_argument, 0, PARAM_PREDICT },
/*-m*/{ "mode", required_argument, 0, 'm' },
/*-N*/
/*-k*/
/*-n*/{ "nmin", required_argument, 0, 'n' },
/*  */{ "nthread", required_argument, 0, PARAM_NTHREAD },
/*-s*/{ "seed", required_argument, 0, 's' },
/*  */{ "noheader", no_argument, 0, PARAM_NOHEADER },
/*  */{ "headercol", no_argument, 0, PARAM_HEADERCOL },
/*-?*/{ "help", no_argument, 0, '?' },
/*  */{ 0, 0, 0, 0 } };

static const char options_short[] = "m:N:k:n:s:?";

/*************************************************************************************************
 * Help message
 ************************************************************************************************/
static const char *help_message =
    "USAGE:\n"
        "  ert --build [OPTIONS] LS\n"
        "  ert --predict [OPTIONS] TS\n"
        "\n"
        "DESCRIPTION:\n"
        "  First command build a set of Extremely Randomized Trees (Geurts, 2006) from LS and output result to standard output.\n"
        "  Second command read a set of Extremely Randomized Trees from standard input and makes prediction based on TS.\n"
        "\n"
        "BUILD_OPTIONS:\n"
        "  -m, --mode,      'R'/'r': regression mode, 'C'/'c': classification mode, 'P'/'p': periodic regression mode [default='%c']\n"
        "  -N,              The number if extremely randomized trees to build [default=%u]\n"
        "  -k,              The number of attributes to consider at each split (see Geurt 2006, for more informations) [default=%u]\n"
        "  -n, --nmin,      The minimal number of instances in a node in order to split it [default=%u]\n"
        "      --nthread,   The number of thread to use [default=%u]\n"
        "  -s, --seed,      The random seed to use in order to build the trees [default=%u]\n"
        "      --noheader,  If set, LS or TS are not supposed to contain a header line\n"
        "      --headercol, Consider the first column of LS/TS as a header and skip it\n"
        "  -?, --help,      Print this help message\n"
        "\n"
        "FILE FORMAT:\n"
        "  Each line of LS is assumed to be in the form \"y x0 x1 x2 x3 ...\n"
        "    where\n"
        "      y  is an integer value to predict in classification mode or a floating point value in (periodic) regression\n"
        "      xi are double values of attributes\n"
        "  Each line of TS is assumed to be in the form \"x0 x1 x2 x3 ...\n"
        "    where xi are double values of attributes\n"
        "  Output of \"ert --predict\" is in the form \"y, yuncertainty\"\n"
        "    where y is the predicted class/value and where yuncertainty is the uncertainty that is associated with y\n"
        "\n"
        "EXAMPLES:\n"
        "  ert --build LS > model.ert : Build a set of Extremely randomized tree and store them in file \"model.ert\"\n"
        "  cat model.ert | ert --predict TS : Make predictions based on TS and on the previously build model\n"
        "  ert --build LS | ert --predict TS: Combine both previous commands\n"
        "\n"
        "NOTES:\n"
        "  If k == 0, then k = nattr-1 in (periodic) regression mode or rint(sqrt(nattr-1)) in classification mode where nattr is the number of attributes in LS\n"
        "  If nmin == 0, then nmin = 5 in (periodic) regression mode or 2 in classification mode\n"
        "\n"
        "REFERENCE:\n"
        "P. Geurts, D. Ernst., and L. Wehenkel, \"Extremely randomized trees\", Machine Learning, 63(1), 3-42, 2006\n";

void print_help(FILE *out) {
  fprintf(out, help_message, DEFAULT_MODE, DEFAULT_N, DEFAULT_K, DEFAULT_NMIN,
  DEFAULT_NTHREAD, DEFAULT_SEED);
}

/*************************************************************************************************
 * Extremely randomized trees building functions
 ************************************************************************************************/
int build_dataset(double **X, void **y, size_t *nobs, size_t *nattr, const fload_data_t *input_data) {
  const size_t ifrom = (noheader) ? 0 : 1;
  const size_t jfrom = (headercol) ? 1 : 0;

  *nobs = fload_nobs(input_data) - ifrom;
  *nattr = fload_nattr(input_data) - jfrom - 1;

  double (*_X)[*nattr];
  void *_y;

  size_t ix, jx;
  size_t iin, jin;

  char *s, *endptr;

  // Allocate X and y
  _X = (double (*)[*nattr]) calloc(*nattr * *nobs, sizeof(double));
  _y = (mode == 'c') ? calloc(*nobs, sizeof(int)) : calloc(*nobs, sizeof(double));
  if (_X == NULL || _y == NULL) {
    free(_X);
    free(_y);
    return -1;
  }

  // Copy input data into X and y
  for (ix = 0, iin = ifrom; ix < *nobs; ix++, iin++) {
    s = endptr = (char *) fload_get(input_data, iin, jfrom);
    switch (mode) {
    case 'c':
      ((int *) _y)[ix] = strtol(s, &endptr, 10);
      break;
    case 'p':
    case 'r':
      ((double *) _y)[ix] = strtod(s, &endptr);
      break;
    }

    if (s == endptr || *endptr != '\0') {
      fprintf(stderr,
          "Unable to parse \"%s\" as a double (outcome value of line %zu)\n", s,
          iin + 1 - ifrom);
      free(_X);
      free(_y);
      return -1;
    }

    for (jx = 0, jin = jfrom + 1; jx < *nattr; jx++, jin++) {
      s = endptr = (char *) fload_get(input_data, iin, jin);
      _X[ix][jx] = strtod(s, &endptr);
      if (s == endptr || *endptr != '\0') {
        fprintf(stderr,
            "Unable to parse \"%s\" as a double (%zuth attribute of line %zu)\n",
            s, jin + 1 - jfrom, iin + 1 - ifrom);
        free(_X);
        free(_y);
        return -1;
      }
    }
  }

  // Copy pointers
  *X = (double *) _X;
  *y = _y;

  return 0;
}

int build(const char *filename) {
  FILE *in;
  fload_data_t *input_data;
  ExtraTrees *ert;
  double *X;
  void *y;
  size_t ninst, nattr;
  ExtraTreesMode ert_mode = ExtraTreesRegressionMode;
  threadpool_t *pool;

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

  // Build the input dataset
  if (build_dataset(&X, &y, &ninst, &nattr, input_data)) {
    fprintf(stderr, "Unable to build the learning set of observations\n");
    fload_destroy(input_data);
    return -1;
  }

  // Destroy the input data
  fload_destroy(input_data);

  // Set the prediction mode to use
  switch (mode) {
  case 'c':
    ert_mode = ExtraTreesClassificationMode;
    break;
  case 'r':
    ert_mode = ExtraTreesRegressionMode;
    break;
  case 'p':
    ert_mode = ExtraTreesPeriodicRegressionMode;
    break;
  }

  // Set default parameters to k and/or nmin in case these are set to 0
  if (k == 0 || nmin == 0) {
    unsigned int kdefault, nmindefault;
    extratree_default_parameters(&kdefault, &nmindefault, nattr,
        (mode == 'p' || mode == 'r'));
    if (k == 0)
      k = kdefault;
    if (nmin == 0)
      nmin = nmindefault;
  }

  // Allocate the thread pool
  pool = threadpool_create(nthread);
  if (!pool) {
    fprintf(stderr, "Unable to allocate a thread pool with %u threads\n",
        nthread);
    free(X);
    free(y);
    return -1;
  }

  // Build the ExtraTrees
  ert = ExtraTrees_build(X, y, ninst, nattr, ert_mode, N, k, nmin, seed, pool);
  if (!ert) {
    fprintf(stderr, "Error while building the Extremely randomized trees\n");
    free(X);
    free(y);
    threadpool_destroy(pool);
    return -1;
  }

  // Delete the dataset and the thread pool
  free(X);
  free(y);
  threadpool_destroy(pool);

  // Dump the extra trees to stdout
  if (ExtraTrees_dump(stdout, ert)) {
    fprintf(stderr, "Error while dumping the Extremely randomized trees\n");
    ExtraTrees_free(ert);
    return -1;
  }

  // Delete the extra trees
  ExtraTrees_free(ert);

  return 0;
}

/*************************************************************************************************
 * Prediction functions
 ************************************************************************************************/
typedef struct predict_input {
  double *x;
  size_t iline;
} predict_input_t;

typedef struct predict_wait_node {
  ExtraTreesPrediction *pred;
  size_t iline;
  struct predict_wait_node *next;
} predict_wait_node_t;

struct {
  pthread_mutex_t lock;
  ExtraTrees *ert;
  threadpool_t *pool;
  predict_wait_node_t *wait_list;
  size_t nextout;
  int status;
} predict_args = {
PTHREAD_MUTEX_INITIALIZER,
NULL,
NULL,
NULL, 0, 0 };

void predict_destroy_wait_list(void) {
  predict_wait_node_t *next;

  while (predict_args.wait_list) {
    next = predict_args.wait_list->next;
    free(predict_args.wait_list->pred);
    free(predict_args.wait_list);
    predict_args.wait_list = next;
  }
}

void predict_input_free(void *arg) {
  predict_input_t *input = (predict_input_t *) arg;
  free(input->x);
  free(input);
}

void predict_output(const ExtraTreesPrediction *pred) {
  switch (ExtraTrees_mode(predict_args.ert)) {
  case ExtraTreesClassificationMode:
    fprintf(stdout, "%d %f\n", pred->classification_val, pred->uncertainty);
    break;
  case ExtraTreesRegressionMode:
    fprintf(stdout, "%g %g\n", pred->regression_val, pred->uncertainty);
    break;
  case ExtraTreesPeriodicRegressionMode:
    fprintf(stdout, "%g %g\n", pred->pregression_val, pred->uncertainty);
    break;
  }
}

void predict_thread(void *arg) {
  predict_input_t *input = (predict_input_t *) arg;
  ExtraTreesPrediction *pred;

  // If we are in an error state, return
  if (predict_args.status)
    return;

  // Get the ERT prediction
  pred = ExtraTrees_get(predict_args.ert, input->x);
  if (!pred) {
    predict_args.status = 1;
    return;
  }

  // Try to acquire the lock for outputting the prediction line
  if (pthread_mutex_lock(&predict_args.lock)) {
    fprintf(stderr, "Unable to acquire output lock (%s)\n", strerror(errno));
    predict_args.status = 1;
    return;
  }

  // If we have to output the prediction, do so
  if (predict_args.nextout == input->iline) {
    // Output this prediction, then free it
    predict_output(pred);
    predict_args.nextout++;
    free(pred);

    // Output all previously stored prediction
    while (predict_args.wait_list
        && predict_args.nextout == predict_args.wait_list->iline) {
      // Remove this node from the head of the waiting list
      predict_wait_node_t *node = predict_args.wait_list;
      predict_args.wait_list = node->next;
      // Output this prediction
      predict_output(node->pred);
      // Free the prediction and the node
      free(node->pred);
      free(node);
      // Get the next line to output
      predict_args.nextout++;
    }
  } else {
    // If it is too soon in order to output this line, add it to the waiting list
    predict_wait_node_t *prev, *current;
    predict_wait_node_t *node = (predict_wait_node_t *) malloc(
        sizeof(predict_wait_node_t));
    if (!node) {
      predict_args.status = 1;
      pthread_mutex_unlock(&predict_args.lock);
      return;
    }
    // Fill the node's informations
    node->pred = pred;
    node->iline = input->iline;
    // Add this node to the (sorted) waiting list
    prev = NULL;
    current = predict_args.wait_list;
    while (current && current->iline < node->iline) {
      prev = current;
      current = current->next;
    }
    // If there is no previous, then we are in the head of the list
    if (!prev)
      predict_args.wait_list = node;
    else
      prev->next = node;
    node->next = current;
  }

  // Release the lock
  pthread_mutex_unlock(&predict_args.lock);
}

int predict_launch_thread(const unsigned int iline, const char **fields,
    const unsigned int nfield, void *arg) {
  char *endptr;
  size_t i, j;

  if (predict_args.status)
    return -1;

  if (noheader || iline > 0) {
    predict_input_t *input = (predict_input_t *) malloc(
        sizeof(predict_input_t));
    if (input == NULL)
      return -1;

    input->iline = iline - (noheader == 0);
    input->x = (double *) calloc(nfield, sizeof(double));
    if (input->x == NULL) {
      free(input);
      return -1;
    }

    for (i = (headercol) ? 1 : 0, j = 0; i < nfield; i++, j++) {
      input->x[j] = strtod(fields[i], &endptr);
      if (fields[i] == endptr || *endptr != '\0') {
        fprintf(stderr,
            "Unable to parse \"%s\" as a double (%zuth attribute of line %u)\n",
            fields[i], j + 1, iline + 1);
        free(input->x);
        free(input);
        return -1;
      }
    }

    return threadpool_addjob(predict_args.pool, predict_thread, input,
        predict_input_free);
  }

  return 0;
}

int predict(const char *filename) {
  ExtraTrees *ert;
  FILE *in;
  threadpool_t *pool;

  // Load the ERT from stdin
  ert = ExtraTrees_load(stdin);
  if (!ert) {
    fprintf(stderr,
        "Unable to load the Extremely randomized trees from the standard input\n");
    return -1;
  }

  // Open the input file in read mode
  in = fopen(filename, "r");
  if (!in) {
    fprintf(stderr, "Unable to open input file %s (%s)\n", filename,
        strerror(errno));
    return -1;
  }

  // Allocate the thread pool
  pool = threadpool_create(nthread);
  if (!pool) {
    fprintf(stderr, "Unable to allocate a thread pool with %u threads\n",
        nthread);
    fclose(in);
    ExtraTrees_free(ert);
    return -1;
  }

  // Assign the ERT to the parameters we will pass to the thread
  predict_args.ert = ert;

  // Assign the thread pool to the parameters we will pass to the thread
  predict_args.pool = pool;

  // Parse the input file
  if (fparse(in, ' ', predict_launch_thread, NULL)) {
    // if(parse_file(filename, &ninst, &nattr, predict_launch_thread, NULL)) {
    fprintf(stderr, "Unable to parse file %s\n", filename);
    ExtraTrees_free(ert);
    threadpool_destroy(pool);
    predict_destroy_wait_list();
    return -1;
  }

  // Wait for all threads to finish
  threadpool_wait(pool);

  // Delete the thread pool
  threadpool_destroy(pool);

  // Close the input file
  fclose(in);

  // Delete the extra trees
  ExtraTrees_free(ert);

  // In case of error
  if (predict_args.status) {
    predict_destroy_wait_list();
    return -1;
  }

  return 0;
}

/*************************************************************************************************
 * Main program
 ************************************************************************************************/
int main(int argc, char **argv) {
  const char *filename;
  int opt;

  // Get the program parameters
  while ((opt = getopt_long(argc, argv, options_short, options_long, NULL))
      != -1) {
    switch (opt) {
    case PARAM_BUILD:
      RunMode = RunModeBuild;
      break;
    case PARAM_PREDICT:
      RunMode = RunModePredict;
      break;
    case 'm':
      mode = optarg[0];
      break;
    case 'N':
      N = strtoul(optarg, NULL, 10);
      break;
    case 'k':
      k = strtoul(optarg, NULL, 10);
      break;
    case 'n':
      nmin = strtoul(optarg, NULL, 10);
      break;
    case PARAM_NTHREAD:
      nthread = strtoul(optarg, NULL, 10);
      break;
    case 's':
      seed = strtoul(optarg, NULL, 10);
      break;
    case PARAM_NOHEADER:
      noheader = 1;
      break;
    case PARAM_HEADERCOL:
      headercol = 1;
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

  // Check command line arguments
  if (RunMode == RunModeNone) {
    fprintf(stderr,
        "No run mode selected (choices are \"--build\" or \"--predict\")\n");
    exit(EXIT_FAILURE);
  }
  switch (mode) {
  case 'c':
  case 'C':
    mode = 'c';
    break;
  case 'r':
  case 'R':
    mode = 'r';
    break;
  case 'p':
  case 'P':
    mode = 'p';
    break;
  default:
    fprintf(stderr, "Unrecognized mode '%c' (choices are 'p'/'r'/'c')\n", mode);
    exit(EXIT_FAILURE);
  }
  if (N == 0) {
    fprintf(stderr, "N parameter must be greater than zero\n");
    exit(EXIT_FAILURE);
  }
  if (nthread == 0) {
    fprintf(stderr, "nthread parameter must be greater than zero\n");
    exit(EXIT_FAILURE);
  }

  if (RunMode == RunModeBuild) {
    if (build(filename))
      exit(EXIT_FAILURE);
  } else {
    if (predict(filename))
      exit(EXIT_FAILURE);
  }

  // Exit
  exit(EXIT_SUCCESS);
}
