#include <fload.h>
#include <float.h>
#include <errno.h>
#include <getopt.h>
#include <graddlsqr.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <utils.h>

/*************************************************************************************************
 * Program parameters
 ************************************************************************************************/
// The default options of the gradient descent algorithm
graddlsqr_opt_t opts = GRADDLSQR_DEFAULT_OPTS;

// Parameter values
char *param_winit = NULL;
char *param_balance = NULL;
char param_wobs[FILENAME_MAX] = { '\0' };
char param_verbose = 0;

// Long-only option identifier
#define PARAM_NTHREAD 65536
#define PARAM_B1      65537
#define PARAM_B2      65538
#define PARAM_EPSILON 65539

// Parameter list
static struct option options_long[] = {

// Options of the stochastic gradient descent algorithm
    /*-e*/{ "eps", required_argument, 0, 'e' },
    /*-a*/{ "alpha", required_argument, 0, 'a' },
    /*  */{ "B1", required_argument, 0, PARAM_B1 },
    /*  */{ "B2", required_argument, 0, PARAM_B2 },
    /*  */{ "epsilon", required_argument, 0, PARAM_EPSILON },
    /*-b*/{ "batch-size", required_argument, 0, 'b' },
    /*  */{ "nthread", required_argument, 0, PARAM_NTHREAD },
    /*-N*/{ "maxiter", required_argument, 0, 'N' },
    /*-v*/{ "verbose", no_argument, 0, 'v' },

// Weight initialization procedures
    /*-w*/{ "winit", required_argument, 0, 'w' },
    /*-W*/{ "wobs", required_argument, 0, 'W' },
    /*-B*/{ "balance", required_argument, 0, 'B' },

// Help option
    /*-?*/{ "help", no_argument, 0, '?' },
    /*  */{ 0, 0, 0, 0 }

};

static const char options_short[] = "e:a:b:N:vw:W:B:?";

/*************************************************************************************************
 * Help message
 ************************************************************************************************/
static const char *help_message =
    "USAGE:\n"
        "  split [OPTIONS] input\n"
        "\n"
        "DESCRIPTION:\n"
        "  Find an hyperplane that split at best a learning set of "
        "observations regarding a given binary class.\n"
        "\n"
        "  The hyperplane is obtained by minimizing the cost function\n"
        "\n"
        "    C = 0.5 * | diag(W) * ( sig( X' * a + b ) - y ) |^2\n"
        "\n"
        " through stochastic gradient descent.\n"
        "\n"
        " In the latter equation:\n"
        "   W         is the vector of weights associated with each observations\n"
        "             (see --wobs option).\n"
        "   w = (a b) is the vector of weights defining the hyperplane.\n"
        "   X         is the matrix of input observation.\n"
        "   sig(x)    is the sigmoid function defined as sig(x) = 1 / (1 + exp(-x)).\n"
        "   y         is the vector of binary outcome values.\n"
        "\n"
        "   The stochastic gradient descent uses the Adam algorithm described in\n"
        " Kingma & Ba, \"ADAM: A METHOD  FOR STOCHASTIC OPTIMIZATION\", ICLR 2015,\n"
        " arXiv:1412.6980.\n"
        "\n"
        "ALGORITHM:\n"
        "   C(-1) = +Infinity;\n"
        "   m(0) = 0; { Initialize the first momentum vector }\n"
        "   v(0) = 0; { Initialize the second momentum vector }\n"
        "   t = 0;\n"
        "   while | C(t) - C(t-1) | >= eps && t < maxiter\n"
        "     t = t + 1;\n"
        "     { Compute the gradient of C wrt. w}\n"
        "     g(t) = dC(t) / dw;\n"
        "     { Update the biased first moment estimate }\n"
        "     m(t) = B1 * m(t-1) + ( 1 - B1 ) * g(t);\n"
        "     { Update the biased second moment estimate }\n"
        "     v(t) = B2 * v(t-1) + ( 1 - B2 ) * g(t)^2;\n"
        "     { Update parameters }\n"
        "                             sqrt( 1 - B2^t )             m(t)\n"
        "     w(t) = w(t-1) - alpha * ---------------- * ------------------------\n"
        "                                 1 - B1^t        sqrt( v(t) ) + epsilon\n"
        "\n"
        "BUILD_OPTIONS:\n"
        "  -e, --eps,        The tolerance in the change of the cost function between two\n"
        "                    iterations [default = %g]\n"
        "  -a, --alpha,      The learning rate of the gradient descent algorithm [default = %g]\n"
        "      --B1,         The exponential decay rate for the first moment estimates [default = %g]\n"
        "      --B2,         The exponential decay rate for the second moment estimates [default = %g]\n"
        "      --epsilon,    A very small number that prevents any division by zero [default = %g]\n"
        "  -b, --batch-size, The size of the chunck to use in online mode [default = %lu]\n"
        "      --nthread,    The number of thread to use [default=%u]\n"
        "  -N, --maxiter,    The maximal number of iterations to perform [default = %lu]\n"
        "  -v, --verbose,    Print various progression message during execution.\n"
        "  -w, --winit,      Initialize w with the specified weights instead of producing an\n"
        "                    heuristic initial estimate (e.g. --winit \"0.2 0.4 0.3 1\")\n"
        "  -W, --wobs,       Path to the the file containing the weights associated with each\n"
        "                    observation, W. By default all weights are set to 1.\n"
        "                    Note that the --balance option can still modify these weights.\n"
        "  -B, --balance,    Scale the W vector according to the provided weights, Sp and Sn:\n"
        "                      W[i] = Sp * W[i] / Wp, if y[i] == 1\n"
        "                      W[i] = Sn * W[i] / Wn, if y[i] == 0\n"
        "                    where Wp and Wn are respectively the sum of the weights of the\n"
        "                    positive and negative instances and where Sp, Sn are provided as\n"
        "                    parameters of the option (i.e. --balance \"Sp Sn\").\n"
        "                    If Sp = Sn, then positive and negative instances will have the same\n"
        "                    influence on the final output.\n"
        "  -?, --help,     Print this help message\n"
        "\n"
        "FILE FORMAT:\n"
        "  Each line of input is assumed to be in the form \"y x0 x1 x2 x3 ...\n"
        "    where\n"
        "      y  is either 0 or 1\n"
        "      xi are double values of attributes.\n";

void print_help(FILE *out) {
  fprintf(out, help_message, opts.eps, opts.alpha, opts.B1, opts.B2, opts.epsilon, opts.batch_size,
      opts.nthread, opts.maxiter);
}

/*************************************************************************************************
 * File loading function
 ************************************************************************************************/
int load_learning_set(double **X, double **y, size_t *ninst, size_t *nattr,
    const char *filename) {
  FILE *in;
  fload_data_t *input_data;

  // Print debug message
  if (param_verbose)
    fprintf(stderr, "Reading learning set from file %s...\n", filename);

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
  *ninst = fload_nobs(input_data);
  *nattr = fload_nattr(input_data) - 1;
  if (param_verbose)
    fprintf(stderr, "Found %zu instances with %zu attributes\n", *ninst,
        *nattr);

  // Allocate the X matrix and y array
  *X = (double *) calloc(*ninst * *nattr, sizeof(double));
  *y = (double *) calloc(*ninst, sizeof(double));
  if (*X == NULL || *y == NULL) {
    fload_destroy(input_data);
    free(*X);
    free(*y);
    return -1;
  }

  { // Convert the fload_data_t structure into an array
    double (*_X)[*nattr] = (double (*)[*nattr]) *X;
    double *_y = *y;
    size_t i, j;

    // Copy the outcome value into y
    for (i = 0; i < *ninst; i++) {
      if (parse_double(&_y[i], fload_get(input_data, i, 0), NULL, 1)
          || (_y[i] != 0. && _y[i] != 1.)) {
        fprintf(stderr,
            "Unable to parse \"%s\" as a boolean (outcome value at line %zu)\n",
            fload_get(input_data, i, 0), i + 1);
        fload_destroy(input_data);
        free(*X);
        free(*y);
        return -1;
      }
    }

    // Copy the input attributes into X
    for (i = 0; i < *ninst; i++) {
      for (j = 0; j < *nattr; j++) {
        if (parse_double(&_X[i][j], fload_get(input_data, i, j + 1), NULL, 1)) {
          fprintf(stderr,
              "Unable to parse \"%s\" as a double (%zuth attribute of line %zu)\n",
              fload_get(input_data, i, j + 1), j + 1, i + 1);
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
 * Global variables
 ************************************************************************************************/

// The matrix of input attributes (size: ninst x nattr)
double *X = NULL;

// The vector of outcome values (size: ninst, y in {0, 1})
double *y = NULL;

// The vector of weights associated with each instance
double *W = NULL;

// The current set of weights defining the hyperplane (size: nattr+1)
double *w = NULL;

// The number of instances
size_t ninst = 0;

// The number of attributes
size_t nattr = 0;

/*************************************************************************************************
 * Basic functions
 ************************************************************************************************/

/**
 * The activation function
 */
double act_fct(const double x) {
  if (x < 0) {
    double z = exp(x);
    return z / (1. + z);
  } else
    return 1. / (1. + exp(-x));
}

/**
 * The derivative of the activation function
 */
double act_fct_deriv(const double x) {
  const double y = act_fct(x);
  return y * (1. - y);
}

/**
 * The input of the activation function
 */
double act_fct_input(const double *w, const double *x,
    const unsigned long nattr) {
  double res;
  size_t i;

  res = w[nattr];
  for (i = 0; i < nattr; i++)
    res += x[i] * w[i];

  return res;
}

double f(const double *w, const double *x, const unsigned long nattr,
    void * arg) {
  return act_fct(act_fct_input(w, x, nattr));
}

void df(double *grad, const double *w, const double *x,
    const unsigned long nattr, void *arg) {
  size_t i;

  grad[nattr] = act_fct_deriv(act_fct_input(w, x, nattr));
  for (i = 0; i < nattr; i++)
    grad[i] = grad[nattr] * x[i];
}

/*************************************************************************************************
 * Weights initialization procedures
 ************************************************************************************************/

/**
 * Initialize the weights vector defining the hyperplane
 */
int w_init(void) {
  size_t i;

  // If a set of weights was given as argument
  if (param_winit != NULL) {
    char *s, *endptr;

    // Parse the string containing the weights
    for (i = 0, s = param_winit; i <= nattr && *s; i++, s = endptr) {
      if (parse_double(&w[i], s, &endptr, 0)) {
        fprintf(stderr,
            "Unable to parse the initial hyperplane weights \"%s\" near \"%s\"\n",
            param_winit, s);
        return -1;
      }
    }

    // Check that we have enough weights
    if (i <= nattr) {
      fprintf(stderr,
          "Not enough weights given for describing the hyperplane (%zu given, %zu expected)\n",
          i, nattr + 1);
      return -1;
    }

    // Check that we parse the whole string (but don't stop if we don't)
    if (*s)
      fprintf(stderr,
          "Too many weights given for describing the hyperplane, ignoring \"%s\"\n",
          s);
  } else {
    /*
     // Otherwise, initialize the weights vector to the one defining the hyperplane
     // that is perpendicular to the line joining the mean of the positive instances and the
     // mean of the negative instances
     double *meanp = w, *dmean = w; // w, meanp and dmean are never used together
     double meann[nattr];
     double *midpoint = meann;
     double sumwp = 0, sumwn = 0, tmp;
     size_t i, j;

     // Initialize meanp and meann
     for (i = 0; i < nattr; i++)
     meanp[i] = meann[i] = 0;

     // Compute the mean value of the positive and negative instances
     for (i = 0; i < ninst; i++) {
     const double *x = X + i * nattr;

     // If this is a positive instance
     if (y[i] == 1.) {
     for (j = 0; j < nattr; j++)
     meanp[j] += W[i] * x[j];
     sumwp += W[i];
     } else { // If this is a negative instance
     for (j = 0; j < nattr; j++)
     meann[j] += W[i] * x[j];
     sumwn += W[i];
     }
     }
     for (i = 0; i < nattr; i++) {
     meanp[i] /= sumwp;
     meann[i] /= sumwn;
     }

     // Compute the mid-point between both mean vectors and
     // the difference between the mean positive vector and the mid-point
     for (i = 0; i < nattr; i++) {
     // Equiv. meann[i] = sumwp / ( sumwp + sumwn ) * ( meanp[i] + meann[i] )
     midpoint[i] = 0.5 * ( meanp[i] + meann[i] );
     // Equiv. w[i] = w[i] - meann[i]
     dmean[i] = meanp[i] - midpoint[i];
     }

     // Compute the initial w vector
     w[nattr] = 0;
     for (i = 0; i < nattr; i++) {
     // w[i] = dmean[i]; // Equiv. w[i] = w[i];
     w[nattr] -= midpoint[i] * dmean[i];
     }

     // Normalize the w vector such that | w[nattr] | = 1e-10
     tmp = fabs(w[nattr]);
     if(tmp != 0.) { // This may happens... rarely...
     for (i = 0; i <= nattr; i++)
     w[i] /= tmp;
     }
     */

    for (i = 0; i <= nattr; i++)
      w[i] = 0;
  }

  return 0;
}

/**
 * Initialize the weight vector of observations
 */
int wobs_init(double *W, const double *y, const size_t ninst) {
  FILE *in;
  fload_data_t *input_data;
  double wp = 0, wn = 0;
  size_t i;

  // If no file should be loaded, initialize wobs with ones
  if (strlen(param_wobs) == 0) {
    if (param_verbose)
      fprintf(stderr, "Initializing all observation weights to 1\n");
    for (i = 0; i < ninst; i++)
      W[i] = 1;
  } else { // Otherwise, read the provided file
    // Print a debug message
    if (param_verbose)
      fprintf(stderr, "Reading observation weights from file \"%s\"...\n",
          param_wobs);

    // Open the input file in read mode
    in = fopen(param_wobs, "r");
    if (!in) {
      fprintf(stderr, "Unable to open wobs file %s (%s)\n", param_wobs,
          strerror(errno));
      return -1;
    }

    // Load the input file
    input_data = fload(in, ' ');
    if (input_data == NULL) {
      fprintf(stderr, "Unable to load input file %s\n", param_wobs);
      fclose(in);
      return -1;
    }

    // Close the input file
    fclose(in);

    // Check that the file has the right number of columns and lines
    if (fload_nattr(input_data) != 1 || fload_nobs(input_data) != ninst) {
      fprintf(stderr,
          "File %s has %zu rows and %zu columns, expected %zu row and 1 columns\n",
          param_wobs, fload_nobs(input_data), fload_nattr(input_data), ninst);
      fload_destroy(input_data);
      return -1;
    }

    // Parse weights
    for (i = 0; i < ninst; i++) {
      if (parse_double(&W[i], fload_get(input_data, i, 0), NULL, 1)) {
        fprintf(stderr,
            "Unable to parse \"%s\" as a double value from line %zu of file %s\n",
            fload_get(input_data, i, 0), i + 1, param_wobs);
        fload_destroy(input_data);
        return -1;
      }
    }

    // Destroy the input data
    fload_destroy(input_data);
  }

  // Count the sum of the weights of the positive and negative instances
  for (i = 0; i < ninst; i++) {
    if (y[i] == 1)
      wp += W[i];
    else
      wn += W[i];
  }

  // Check that we have two class
  if (wp == 0 || wn == 0) {
    fprintf(stderr, "Only one outcome value affect the cost function\n");
    return -1;
  }

  // If we have to balance observations
  if (param_balance) {
    double S[2];
    const char *s = param_balance;
    char *endptr;

    // Parse Sp and Sn
    if (parse_double(&S[0], s, &endptr, 0)
        || parse_double(&S[1], endptr, NULL, 1)) {
      fprintf(stderr, "Unable to parse \"%s\" as two double values\n", s);
      return -1;
    }

    if (param_verbose)
      fprintf(stderr,
          "Balancing the observation weights. Scaling of the positive/negative instances: %g/%g\n",
          S[0] / wp, S[1] / wn);

    // Adjust weight in order to mimic an equal number of positive and negative instances"
    for (i = 0; i < ninst; i++) {
      if (y[i] == 1)
        W[i] *= S[0] / wp; // This formulation ensure that the total sum of weights
                           // taken over all instances do not change.
      else
        W[i] *= S[1] / wn;
    }
  }

  return 0;
}

/*************************************************************************************************
 * Core functions
 ************************************************************************************************/

void split_result(void) {
  size_t tp, fp, np, nn;
  size_t i;

  // Print the solution to stdout
  for (i = 0; i <= nattr; i++)
    fprintf(stdout, " %.10g", w[i]);
  fprintf(stdout, "\n");

  // In debug mode, also print the true/false positive rates
  if (param_verbose) {
    // Compute the number of (true) positive and (false) positive instances
    tp = fp = np = nn = 0;
    for (i = 0; i < ninst; i++) {
      if (y[i] == 1.) {
        tp += (act_fct_input(w, X + i * nattr, nattr) > 0);
        np++;
      } else {
        fp += (act_fct_input(w, X + i * nattr, nattr) > 0);
        nn++;
      }
    }

    // Print the hyperplane formula + confusion matrix
    fprintf(stderr, "         |  Positive |  Negative | <= Predicted class\n");
    fprintf(stderr, "---------+-----------+-----------|\n");
    fprintf(stderr, "Positive | %9zu | %9zu | TPR = %f\n", tp, np - tp,
        (double) tp / np);
    fprintf(stderr, "Negative | %9zu | %9zu | FPR = %f\n", fp, nn - fp,
        (double) fp / nn);
    fprintf(stderr,
        "---------+-----------+-----------| Precision = %f, Accuracy = %f\n",
        (double) tp / (tp + fp), (double) (tp + nn - fp) / ninst);
  }
}

/**
 * Optimize the weight vector
 */
int split_build(const char *filename) {
  // Load the learning set of observations
  if (load_learning_set(&X, &y, &ninst, &nattr, filename)) {
    fprintf(stderr,
        "Unable to load the learning set of observation from file %s\n",
        filename);
    return -1;
  }

  // Allocate the vector of weights defining the hyperplane
  w = (double *) calloc(nattr + 1, sizeof(double));
  if (w == NULL) {
    fprintf(stderr, "Unable to allocate the weights defining the hyperplane\n");
    free(X);
    free(y);
    return -1;
  }

  // Allocate the vector of weights associated with each observation
  W = (double *) calloc(ninst, sizeof(double));
  if (W == NULL) {
    fprintf(stderr,
        "Unable to allocate the weights associated with each observations\n");
    free(X);
    free(y);
    free(w);
    return -1;
  }

  // Initialize the vectors of weights
  if (wobs_init(W, y, ninst) || w_init()) {
    free(X);
    free(y);
    free(w);
    return -1;
  }

  // Optimize the weights of the hyperplane through gradient descent
  if (graddlsqr(w, nattr + 1, X, y, W, ninst, nattr, f, df, NULL, &opts)) {
    free(X);
    free(y);
    free(w);
    free(W);
    return -1;
  }

  // Print result
  split_result();

  // Delete X, y, w and W
  free(X);
  free(y);
  free(w);
  free(W);

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
      case 'e': /* --eps */
        opts.eps = strtod(optarg, NULL);
        break;
      case 'a': /* --alpha */
        opts.alpha = strtod(optarg, NULL);
        break;
      case PARAM_B1: /* --B1 */
        opts.B1 = strtod(optarg, NULL);
        break;
      case PARAM_B2: /* --B2 */
        opts.B2 = strtod(optarg, NULL);
        break;
      case PARAM_EPSILON: /* --epsilon */
        opts.epsilon = strtod(optarg, NULL);
        break;
      case 'b': /* --batch-size */
        opts.batch_size = strtoul(optarg, NULL, 10);
        break;
      case PARAM_NTHREAD: /* --nthread */
        opts.nthread = strtoul(optarg, NULL, 10);
        break;
      case 'N': /* --maxiter */
        opts.maxiter = strtoul(optarg, NULL, 10);
        break;
      case 'v': /* --verbose */
        opts.verbose = param_verbose = 1;
        break;
      case 'w': /* --winit */
        param_winit = (char *) calloc(strlen(optarg) + 1, sizeof(char));
        strcpy(param_winit, optarg);
        break;
      case 'W': /* --wobs */
        strncpy(param_wobs, optarg, FILENAME_MAX);
        break;
      case 'B': /* --balance */
        param_balance = (char *) calloc(strlen(optarg) + 1, sizeof(char));
        strcpy(param_balance, optarg);
        break;
      case '?': /* --help */
        print_help(stdout);
        exit(EXIT_SUCCESS);
        break;
      default:
        print_help(stderr);
        free(param_winit);
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
    free(param_winit);
    exit(EXIT_FAILURE);
  }

  // Build the splitting hyperplane
  if (split_build(filename)) {
    fprintf(stderr, "An error occurs while building the hyperplane\n");
    free(param_winit);
    exit(EXIT_FAILURE);
  }

  // Delete winit (if we allocated it)
  free(param_winit);

  // Exit
  exit(EXIT_SUCCESS);
}
