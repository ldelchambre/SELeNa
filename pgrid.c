#include <errno.h>
#include <fload.h>
#include <float.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <utils.h>

/*************************************************************************************************
 * Program parameters
 ************************************************************************************************/
#define DEFAULT_STEP 1e-2

// Parameter values
double param_xstep = DEFAULT_STEP;
double param_ystep = DEFAULT_STEP;
char param_verbose = 0;
enum {
  PGRID_PARAM_BUILD, PGRID_PARAM_PREDICT
} param_mode = PGRID_PARAM_BUILD;

// Long-only option identifier
#define PARAM_BUILD 65536
#define PARAM_PREDICT 65537

// Parameter list
static struct option options_long[] = {
/*  */{ "build", no_argument, 0, PARAM_BUILD },
/*  */{ "predict", no_argument, 0, PARAM_PREDICT },
/*-x*/{ "xstep", required_argument, 0, 'x' },
/*-y*/{ "ystep", required_argument, 0, 'y' },
/*-s*/{ "step", required_argument, 0, 's' },
/*-v*/{ "verbose", no_argument, 0, 'v' },
/*-?*/{ "help", no_argument, 0, '?' },
/*  */{ 0, 0, 0, 0 }

};

static const char options_short[] = "x:y:s:v?";

/*************************************************************************************************
 * Help message
 ************************************************************************************************/
static const char *help_message =
    "USAGE:\n"
        "  pgrid [ --build | --predict ] [OPTIONS] input\n"
        "\n"
        "DESCRIPTION:\n"
        "  Build a probability grid based on the two dimensional histograms of the input values\n"
        "  or predict probabilities based on a built grid.\n"
        "\n"
        "BUILD_OPTIONS:\n"
        "  -x, --xstep,    The step size in the x drection [default = %g]\n"
        "  -y, --ystep,    The step size in the y drection [default = %g]\n"
        "  -s, --step,     Combination of --xstep and --ystep [default = %g]\n"
        "\n"
        "GENERAL OPTIONS:\n"
        "  -v, --verbose,  Should we be verbose.\n"
        "  -?, --help,     Print this help message\n"
        "\n"
        "FILE FORMAT:\n"
        "  Each line of input is assumed to be in the form \"class x y\" in build mode or \"x y\" in prediction mode\n"
        "    where\n"
        "      class  is either 0 or 1\n"
        "      x and y are double values of attributes."
        "\n"
        "INPUT:\n"
        "  In prediction mode, standard input must be the the binary grid as outputted from the build mode.\n"
        "\n"
        "OUTPUT:\n"
        "  In build mode, standard output is the probability grid, in binary form.\n"
        "  In prediction mode, standard output is the found probability.\n";

void print_help(FILE *out) {
  fprintf(out, help_message, param_xstep, param_ystep);
}

/*************************************************************************************************
 * Core functions
 ************************************************************************************************/
int pgrid_build(const char *filename) {
  FILE *in;
  fload_data_t *input_data;
  size_t nx = 0, ny = 0, np = 0, nn = 0;
  size_t ninst = 0;
  double x, y, xmin = DBL_MAX, xmax = -DBL_MAX, ymin = DBL_MAX, ymax = -DBL_MAX;
  typedef struct _pgrid_counts {
    size_t np;
    size_t nn;
  } _pgrid_counts_t;
  _pgrid_counts_t *pgrid = NULL;
  size_t i, j;

  // Print debug message
  if (param_verbose)
    fprintf(stderr, "Reading input from file %s...\n", filename);

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

  // Check the number of observations/attributes
  ninst = fload_nobs(input_data);
  if (fload_nattr(input_data) != 3) {
    fprintf(stderr, "Wrong number of attributes. Expected \"class x y\".");
    fload_destroy(input_data);
    return -1;
  }

  if (param_verbose) {
    fprintf(stderr, "Found %zu instances\n", ninst);
    fprintf(stderr, "Computing the probability grid...\n");
  }

  // Find the minimal/maximal x,y values
  for (i = 0; i < ninst; i++) {
    if (strcmp(fload_get(input_data, i, 0), "0") != 0
        && strcmp(fload_get(input_data, i, 0), "1") != 0) {
      fprintf(stderr,
          "Line %zu of file %s: \"%s\" is not a valid binary value (0 or 1)\n",
          i + 1, filename, fload_get(input_data, i, 0));
      fload_destroy(input_data);
      return -1;
    }
    if (parse_double(&x, fload_get(input_data, i, 1), NULL, 1)
        || parse_double(&y, fload_get(input_data, i, 2), NULL, 1)) {
      fprintf(stderr, "Unable to parse line %zu from file \"%s\"\n", i + 1,
          filename);
      fload_destroy(input_data);
      return -1;
    }
    if (x < xmin)
      xmin = x;
    if (xmax < x)
      xmax = x;
    if (y < ymin)
      ymin = y;
    if (ymax < y)
      ymax = y;
  }

  // Compute the size of the probability grid
  nx = ceil((xmax - xmin) / param_xstep);
  ny = ceil((ymax - ymin) / param_ystep);

  if(param_verbose)
    fprintf(stderr, "Grid size: (%zu x %zu), x = [%g, %g], y = [%g, %g]\n", nx, ny, xmin, xmax, ymin, ymax);

  // Allocate the probability grid
  pgrid = (_pgrid_counts_t *) calloc(nx * ny, sizeof(_pgrid_counts_t));
  if (pgrid == NULL) {
    fprintf(stderr, "Unable to allocate the probability grid\n");
    fload_destroy(input_data);
    return -1;
  }

  { // Fill the probability grid
    _pgrid_counts_t (*g)[ny] = (_pgrid_counts_t (*)[ny]) pgrid; // Re-interpret grid in order to ease understanding

    // Initialize the grid
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
        g[i][j].nn = g[i][j].np = 0;

    // Browse all instances
    for (i = 0; i < ninst; i++) {
      const char class = fload_get(input_data, i, 0)[0] - '0';
      const size_t ix = (strtod(fload_get(input_data, i, 1), NULL) - xmin)
          / param_xstep;
      const size_t iy = (strtod(fload_get(input_data, i, 2), NULL) - ymin)
          / param_ystep;
      if (class) {
        g[ix][iy].np++;
        np++;
      } else {
        g[ix][iy].nn++;
        nn++;
      }
    }

    // Output grid
    fwrite(&nx, sizeof(size_t), 1, stdout);
    fwrite(&ny, sizeof(size_t), 1, stdout);
    fwrite(&xmin, sizeof(double), 1, stdout);
    fwrite(&ymin, sizeof(double), 1, stdout);
    fwrite(&param_xstep, sizeof(double), 1, stdout);
    fwrite(&param_ystep, sizeof(double), 1, stdout);
    for (i = 0; i < nx; i++) {
      for (j = 0; j < ny; j++) {
        const double p =
            (g[i][j].np + g[i][j].nn > 0) ?
                ((double) g[i][j].np / np)
                    / ((double) g[i][j].np / np + (double) g[i][j].nn / nn) :
                -1;
        fwrite(&p, sizeof(double), 1, stdout);
      }
    }
  }

  // Destroy input file and pgrid
  free(pgrid);
  fload_destroy(input_data);

  return 0;
}

int pgrid_predict(const char *filename) {
  double *pgrid = NULL;
  size_t ix, iy, nx = 0, ny = 0;
  double x, y, xmin = DBL_MAX, ymin = DBL_MAX, xstep = 0, ystep = 0;
  // getline variables
  FILE *in;
  char *line = NULL, *endptr;
  size_t len = 0, iline = 0;
  ssize_t read;

  if (param_verbose)
    fprintf(stderr, "Reading the probability grid\n");

  // Read the informations about the grid
  if (fread(&nx, sizeof(size_t), 1, stdin) == 0)
    return -1;
  if (fread(&ny, sizeof(size_t), 1, stdin) == 0)
    return -1;
  if (fread(&xmin, sizeof(double), 1, stdin) == 0)
    return -1;
  if (fread(&ymin, sizeof(double), 1, stdin) == 0)
    return -1;
  if (fread(&xstep, sizeof(double), 1, stdin) == 0)
    return -1;
  if (fread(&ystep, sizeof(double), 1, stdin) == 0)
    return -1;

  // Allocate the grid
  pgrid = (double *) calloc(nx * ny, sizeof(double));

  // Read the probability grid
  if (fread(pgrid, sizeof(double), nx * ny, stdin) == 0) {
    free(pgrid);
    return -1;
  }

  // Open file in read mode
  in = fopen(filename, "r");
  if (!in) {
    fprintf(stderr, "Unable to open input file %s (%s)\n", filename,
        strerror(errno));
    return -1;
  }

  // Parse the entire input file
  while ((read = getline(&line, &len, in)) != -1) {
    if(line[read - 1] == '\n') line[read - 1] = '\0';
    if(parse_double(&x, line, &endptr, 0) || parse_double(&y, endptr, NULL, 1)) {
      fprintf(stderr,
          "Unable to parse line %zu from file \"%s\": Can not convert \"%s\" into a pair of double",
          iline + 1, filename, line);
      fclose(in);
      free(pgrid);
      free(line);
      return -1;
    }

    // Output probabilities
    ix = (x - xmin) / xstep;
    iy = (y - ymin) / ystep;
    if (xmin <= x && ymin <= y && ix < nx && iy < ny)
      fprintf(stdout, "%g\n", *(pgrid + ix * ny + iy));
    else
      fprintf(stdout, "-1\n");

    // Go to next line
    iline++;
  }

  fclose(in);
  free(pgrid);
  free(line);

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
      param_mode = PGRID_PARAM_BUILD;
      break;
    case PARAM_PREDICT:
      param_mode = PGRID_PARAM_PREDICT;
      break;
    case 'x': /* --xstep */
      param_xstep = strtod(optarg, NULL);
      break;
    case 'y': /* --ystep */
      param_ystep = strtod(optarg, NULL);
      break;
    case 's': /* --step */
      param_xstep = param_ystep = strtod(optarg, NULL);
      break;
    case 'v': /* --verbose */
      param_verbose = 1;
      break;
    case '?': /* --help */
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

  // Build the splitting hyperplane
  switch (param_mode) {
  case PGRID_PARAM_BUILD:
    if (pgrid_build(filename)) {
      fprintf(stderr,
          "An error occurs during the building of the probability grid\n");
      exit(EXIT_FAILURE);
    }
    break;
  case PGRID_PARAM_PREDICT:
    if (pgrid_predict(filename)) {
      fprintf(stderr, "An error occurs during predictions\n");
      exit(EXIT_FAILURE);
    }
    break;
  }

  // Exit
  exit(EXIT_SUCCESS);
}
