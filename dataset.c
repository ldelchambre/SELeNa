/**
 * @file dataset.c
 *
 * Build a learning data set of observations for gravitational lenses
 * identification based either on simulations or on observations.
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
#include <float.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// If defined, the normalization of the candidates will use the faintest image
// instead of the second brightest image for normalization
#define _NORMALIZATION_USE_LAST

/*************************************************************************************************
 * Program parameters
 ************************************************************************************************/
// Default input parameters
#define DEFAULT_SEED     1
#define DEFAULT_INPUT    "1,2,3,4"
#define DEFAULT_OUTPUT   "class,b,e,te,g,tg,s,xs,ys,xd,yd,x1,y1,m1,x2,y2,m2,x3,y3,m3,x4,y4,m4"
#define DEFAULT_XYNOISE   0.0
#define DEFAULT_MAGNOISE 0.0
#define DEFAULT_NNOISE   1

// Variables containing the program parameters
char        *input    = DEFAULT_INPUT;
char        *output   = DEFAULT_OUTPUT;
char         obs      = 0;
char         sph      = 0;
char         norm     = 0;
unsigned int seed     = DEFAULT_SEED;
double       xynoise   = DEFAULT_XYNOISE;
char         nrel     = 0;
double       magnoise = DEFAULT_MAGNOISE;
char         noheader = 0;
unsigned int nnoise   = DEFAULT_NNOISE;

// Result of the parsing of the --input option
struct {
  unsigned int *idx;
  unsigned int n;
} pinput;

// Result of the parsing of the --output option

// The output type
typedef enum poutput_type {
  OUTPUT_CLASS = 0x0001,
  OUTPUT_B     = 0x0002,
  OUTPUT_E     = 0x0004,
  OUTPUT_TE    = 0x0008,
  OUTPUT_G     = 0x0010,
  OUTPUT_TG    = 0x0020,
  OUTPUT_S     = 0x0040,
  OUTPUT_XS    = 0x0080,
  OUTPUT_YS    = 0x0100,
  OUTPUT_XD    = 0x0200,
  OUTPUT_YD    = 0x0400,
  OUTPUT_X     = 0x0800,
  OUTPUT_Y     = 0x1000,
  OUTPUT_MAG   = 0x2000,
  OUTPUT_NIMG  = 0x4000,
} poutput_type_t;

// Matching between output type and string
struct {
  const char *str;
  const poutput_type_t type;
} poutput_str[] = {
    { "class", OUTPUT_CLASS},
    { "b",     OUTPUT_B    },
    { "e",     OUTPUT_E    },
    { "te",    OUTPUT_TE   },
    { "g",     OUTPUT_G    },
    { "tg",    OUTPUT_TG   },
    { "s",     OUTPUT_S    },
    { "xs",    OUTPUT_XS   },
    { "ys",    OUTPUT_YS   },
    { "xd",    OUTPUT_XD   },
    { "yd",    OUTPUT_YD   },
    { "x",     OUTPUT_X    },
    { "y",     OUTPUT_Y    },
    { "m",     OUTPUT_MAG  },
    { "nimg",  OUTPUT_NIMG },
    { NULL,    0           }
};

// Output node (type + optional index)
typedef struct poutput_node {
  poutput_type_t type;
  unsigned int idx;
  struct poutput_node *next;
} poutput_node_t;

// The output parameters
struct {
  poutput_type_t types;
  poutput_node_t *head;
  poutput_node_t *tail;
} poutput;

// Long option return value
#define PARAM_XYNOISE   65538
#define PARAM_MAGNOISE  65539
#define PARAM_NOHEADER  65540
#define PARAM_NNOISE    65541
#define PARAM_SPH       65542

// Long options definitions
static struct option options_long[] = {
/*-i*/{ "input",     required_argument, 0, 'i' },
/*-o*/{ "output",    required_argument, 0, 'o' },
/*-O*/{ "obs",       no_argument,       0, 'O'},
/*  */{ "sph",       no_argument,       0, PARAM_SPH },
/*-n*/{ "norm",      no_argument,       0, 'n'},
/*-s*/{ "seed",      required_argument, 0, 's' },
/*  */{ "xynoise",   required_argument, 0, PARAM_XYNOISE},
/*-r*/{ "nrel",      no_argument,       0, 'r'},
/*  */{ "magnoise",  required_argument, 0, PARAM_MAGNOISE},
/*  */{ "noheader",  no_argument,       0, PARAM_NOHEADER},
/*  */{ "nnoise",    required_argument, 0, PARAM_NNOISE},
/*-?*/{ "help",      no_argument,       0, '?' },
/*  */{ 0,           0,                 0, 0 } };

// Short options definition
static const char options_short[] = "i:o:Oc:ns:rN:?";

/*************************************************************************************************
 * Help message
 ************************************************************************************************/
static const char
    *help_message =
        "USAGE:\n"
        "  dataset [options]\n"
        "\n"
        "DESCRIPTION:\n"
        "  Build a data set of observations for gravitational lenses identification and \n"
        "  characterization based either on simulations of Non-singular Isothermal\n"
        "  ellipsoid lens models in presence of an external shear (NSIEg lens models) or\n"
        "  based on observations.\n"
        "\n"
        "  The procedure is as follow:\n"
        "    1) Read the lensed imaged based on simulation or on observation (--obs option) from\n"
        "       standard input\n"
        "    2) If this configuration do not have enough images for the --input option, then\n"
        "       discard it.\n"
        "    3) Compute the maximal separation between any two images of the system (used by the\n"
        "       --nrel option).\n"
        "    4) Add noise according to the --xnoise, --ynoise, --magnoise and --nrel options\n"
        "    5) Keep only those images given by the --input option (beware that the noise added\n"
        "       on the magnitude may alter the ordering of the lensed images in brightness).\n"
        "    6) Normalize the lens configuration (i.e. --norm option).\n"
        "    7) Output the lens configuration according to the --output option.\n"
        "\n"
        "OPTIONS:\n"
        "  -i, --input,    Comma-separated indices of the lensed images to keep (e.g. -i 1,3 will\n"
        "                  only keep the brightest image, plus the third brightest image)\n"
        "                  [default = \"%s\"].\n"
        "  -o, --output,   Comma-separated values of the columns to export\n"
        "                  Possible values are:\n"
        "                        class : Set to 1 if the current line is associated with a simulated \n"
        "                                lens configuration, 0 if it is associated with an observation\n"
        "                     b, e, te : Parameters on the NSIEg lens model\n"
        "                     g, tg, s\n"
        "                       xs, ys\n"
        "                       xd, yd : The x and y position of the deflector \n"
        "                       xn, yn : The x, y position and magnitude of the nth lensed images (e.g.\n"
        "                           mn   x1 is the x position of the brightest lensed image)\n"
        "                         nimg : The number of images on which this configuration is based\n"
        "                  [default = \"%s\"].\n"
        "  -O, --obs,      Input file is in the form \"nimg x y mag\", rather than in the form described\n"
        "                  in the \"input format\" section (i.e. NSIEg lens model parameters are unknown).\n"
        "                  This is equivalent to have as an input \"0 0 0 0 0 0 0 0 nimg x y mag 0\".\n"
        "      --sph,      Input coordinates are spherical coordinates expressed in radian.\n"
        "                  These coorinates will be converted to ( x * cos(y1), y ) before processing.\n"
#ifdef _NORMALIZATION_USE_LAST
        "  -n, --norm,     Normalize the output configuration so as to have the brightest image standing at\n"
        "                  (0,0) with a magnitude of 0 and the faintest image standing at (1,0). Note that\n"
        "                  the output lens model parameters are tuned accordingly.\n"
#else
        "  -n, --norm,     Normalize the output configuration so as to have the brightest image standing at\n"
        "                  (0,0) with a magnitude of 0 and the second brightest image standing at (1,0). Note that\n"
        "                  the output lens model parameters are tuned accordingly.\n"
#endif
        "  -s, --seed,     The random seed to use [default = %u].\n"
        "      --nnoise,   The number of noise realization to perform [default = %lu].\n"
        "      --xynoise,  The Gaussian noise to add to the x and y positions of the lensed images [default = %g].\n"
        "  -r, --nrel,     If set, then --xynoise will be considered as being relative to the maximal separation\n"
        "                  between the lensed images (not only those specified by the --input option) rather than\n"
        "                  being absolute.\n"
        "      --magnoise, The Gaussian noise to add to the magnitudes of the lensed images [default = %g].\n"
        "      --noheader, Don't print an header line\n"
        "  -?, --help,     Print this help message.\n"
        "\n"
        "INPUT FORMAT:\n"
        "  The mandatory input format is given by \"b e te g tg s xs ys nimg x y mag sdf\", as outputed from\n"
        "  the simulation program. See however the --obs option for a minimal input format.\n"
        "  These input lines are taken from standard input.\n"
        "\n"
        "SEE:\n"
        "  The documentation of the simulation program in order to obtain a detailed description of the b, e, te,\n"
        "  g, tg, s, xs and ys NSIEg lens model parameters.\n";

void printHelp(FILE *out) {
  fprintf(out, help_message, DEFAULT_INPUT, DEFAULT_OUTPUT, DEFAULT_SEED, DEFAULT_NNOISE, DEFAULT_XYNOISE, DEFAULT_MAGNOISE);
}
/*************************************************************************************************
 * Options parsing functions
 ************************************************************************************************/

// Compare unsigned int
int sortuint(const void *a, const void *b) {
  return (*((const unsigned int *) a) < *((const unsigned int *) b)) ? -1 : 1;
}

// parsing of the --input option
int parse_input_option(void) {
  const char *s;
  char *send;
  unsigned int i, j;

  // Count the number of comma we have in the input string
  i = 0;
  for(s = input; *s; s++)
    i += (*s == ',');

  // Set the initial number of indices to i + 1 (i.e. if --input is "1,2,3",
  // then we have 3 indices but two comma
  pinput.n = i + 1;

  // Allocate the array of indices
  pinput.idx = (unsigned int *) calloc(pinput.n, sizeof(unsigned int));
  if(!pinput.idx)
    return -1;

  // Convert indices from string into unsigned int
  s = input;
  for(i = 0; i < pinput.n; i++) {
    pinput.idx[i] = strtoul(s, &send, 10);
    // If conversion failed, that we found a delimiter that differ from ',' or '\0'
    // or that the provided index is equal to zero
    if(s == send || (*send != ',' && *send != '\0') || pinput.idx[i] == 0) {
      free(pinput.idx);
      return -1;
    }
    s = send;
    if(*s == ',') s++;
  }

  // Skip the trailing spaces
  while(isspace(*s))
    s++;

  // If some trailing characters remain, raise an error
  if(*s != '\0') {
    free(pinput.idx);
    return -1;
  }

  // Sort the array of indices
  qsort(pinput.idx, pinput.n, sizeof(unsigned int), sortuint);

  // Remove duplicates from the list of indices and set zero-based indices
  i = j = 0;
  while(j < pinput.n) {
    pinput.idx[i] = pinput.idx[j++] - 1;
    while(pinput.idx[i] == pinput.idx[j])
      j++;
    i++;
  }
  pinput.n = i;

  return 0;
}

// Return string from poutput_type_t
const char *poutput_get_str(const poutput_type_t type) {
  unsigned int i = 0;
  while(poutput_str[i].type != type)
    i++;
  return poutput_str[i].str;
}

// Destroy the parsing of the --output option
void poutput_destroy(void) {
  poutput_node_t *node;

  while(poutput.head) {
    node = poutput.head->next;
    free(poutput.head);
    poutput.head = node;
  }
  poutput.tail = NULL;
}

// Parsing of the --output option
int parse_output_option(void) {
  const char *s;
  char *send;
  unsigned int bmatch, nmatch;
  unsigned int i, imatch;
  poutput_node_t *node;

  // Initialize the poutput structure
  poutput.head  = NULL;
  poutput.tail  = NULL;
  poutput.types = 0;

  // Parse the output string
  s = output;
  while(*s) {
    // Find the best correspondence of this string with one of the poutput_type_t
    bmatch = nmatch = 0;
    for(i = 0; poutput_str[i].str; i++) {
      imatch = strlen(poutput_str[i].str);
      if(imatch > nmatch && strncmp(poutput_str[i].str, s, imatch) == 0) {
        bmatch = i;
        nmatch = imatch;
      }
    }
    // If no match was found
    if(nmatch == 0) {
      poutput_destroy();
      return -1;
    }

    // Allocate the output node
    node = (poutput_node_t *) malloc(sizeof(poutput_node_t));
    if(!node) {
      poutput_destroy();
      return -1;
    }

    // Initialize the output node
    node->type = poutput_str[bmatch].type;
    node->next = NULL;
    poutput.types |= node->type;

    // Add this node to the list
    if(poutput.tail) {
      poutput.tail->next = node;
      poutput.tail = node;
    } else {
      poutput.head = poutput.tail = node;
    }

    // Skip the best matching string
    s += nmatch;

    // In case of "x", "y", "m", also read the associated index
    switch(node->type){
      case OUTPUT_X:
      case OUTPUT_Y:
      case OUTPUT_MAG:
        node->idx = strtoul(s, &send, 10);
        if(s == send || node->idx == 0) {
          poutput_destroy();
          return -1;
        }
        node->idx--;
        s = send;
        break;
      default:
        node->idx = -1U;
        break;
    }

    // If some trailing character remains
    if(*s != ',' && *s != '\0') {
      poutput_destroy();
      return -1;
    }

    // Skip the ',' character
    if(*s == ',') s++;
  }

  return 0;
}

// Check that all options are consistent
int check_options() {
  poutput_node_t *node = poutput.head;

  // Check that the --output options x, y and m have
  // indices that are within the range of the --input option
  while(node) {
    if(node->idx < -1U) {
      if(pinput.n <= node->idx) {
        fprintf(stderr, "Unable to output %s%u whereas only %u images are kept on input\n",
            poutput_get_str(node->type), node->idx+1, pinput.n);
        return -1;
      }
    }
    node = node->next;
  }

  // Check that the --input option is not NULL
  if(pinput.n == 0) {
    fprintf(stderr, "No input image provided\n");
    return -1;
  }

  // If normalization is used, check that we have at least two lensed images on input
  if(norm && pinput.n < 2) {
    fprintf(stderr, "Two images are required in order to normalize the lens configuration, but only one is given by the --input option\n");
    return -1;
  }

  // Check that the --xynoise parameter is not negative
  if(xynoise < 0.) {
    fprintf(stderr, "Error: negative --xynoise parameter (value: %g)\n", xynoise);
    return -1;
  }

  // Check that the magnoise parameter is not negative
  if(magnoise < 0.) {
    fprintf(stderr, "Error: negative --magnoise parameter (value: %g)\n", magnoise);
    return -1;
  }

  // Warn if we perform many realization of noise with no noise
  if(nnoise > 1 && xynoise == 0 && magnoise == 0)
    fprintf(stderr, "Many realization of the noise will be performed by no noise is given");

  // Print a warning message if the list of output parameters is empty
  if(poutput.types == 0)
    fprintf(stderr, "Warning: no output provided\n");

  // Print a warning message if we want to add noise the the observations
  if(obs && (xynoise > 0 || magnoise > 0))
    fprintf(stderr, "Warning: noise will be added to the observations (option --obs combined with --xynoise or --magnoise)\n");

  return 0;
}

/*************************************************************************************************
 * Program functions
 ************************************************************************************************/

// Print a header line to stdout
void print_header(void) {
  poutput_node_t *node = poutput.head;

  while(node) {
    fprintf(stdout, "%s", poutput_get_str(node->type));
    if(node->idx < -1U)
      fprintf(stdout, "%u", node->idx+1);
    node = node->next;
    if(node)
      fprintf(stdout, " ");
  }
  fprintf(stdout, "\n");
}

/**
 * Rotate and scale two dimensional points
 *
 * Given a rotation angle, theta, and a scaling factor, a, we will have that
 *   xt = a * cos(theta)
 *   yt = a * sin(theta)
 * where the resulting rotated points, xp, yp will be given by
 *   xp = a *  ( x * cos(theta) + y * sin(theta) )
 *   yp = a *  ( y * cos(theta) - x * sin(theta) )
 */
void rotate(double *x, double *y, const double xt, const double yt) {
  const double _x = *x, _y = *y;
  *x = _x * xt + _y * yt;
  *y = _y * xt - _x * yt;
}

/**
 * Normalize the NSIEg lens parameters
 *
 * Example for testing the data set normalization 
 * (only valid if _NORMALIZATION_USE_LAST is not defined)
 * ------------------------------------------------------
 * > ./simulation -s 0 1 0.37 54 0.1 162 0.2 0.07 0.05 [-5:1000:5] [-5:1000:5] \
 *     | ./dataset -i "1,2,3,4,5" -o "b,e,te,g,tg,s,xs,ys,xd,yd,x1,y1,m1,x2,y2,m2,x3,y3,m3,x4,y4,m4,x5,y5,m5" --norm \
 *     | awk '(NR > 1){printf("%s %s %s %s %s %s %s %s\n",$1,$2,$3,$4,$5,$6,$7-$9,$8-$10);for(i=11;i<=NF;i+=3) printf("\t%s %s %s\n", $i-$9, $(i+1)-$10, $(i+2));}'
 * $ +6.6740172295589928e-01 +3.7000000000000000e-01 +8.3399853823873826e+01 +1.0000000000000001e-01 +1.9139985382387383e+02 +1.3348034459117986e-01 0.0243201 0.0520065
 * $    -0.459859 0.121092 +0.0000000000000000e+00
 * $    0.540141 0.121092 +2.7874800000000000e-01
 * $    -0.0415128 0.364032 +2.9394200000000015e-01
 * $    -0.0285532 -0.197461 +4.3694940000000004e-01
 * $    -0.016074 -0.0757006 +1.2123800000000000e+00
 * > ./simulation -s 0 +6.6740172295589928e-01 +3.7000000000000000e-01 +8.3399853823873826e+01 +1.0000000000000001e-01 +1.9139985382387383e+02 +1.3348034459117986e-01 0.0243201 0.0520065 [-5:1000:5] [-5:1000:5]
 * $ 0.667402 0.37 83.3999 0.1 191.4 0.13348 0.0243201 0.0520065 5 -4.598594e-01 +1.210917e-01 -1.406727e+00 +3.009266e-34
 * $ 0.667402 0.37 83.3999 0.1 191.4 0.13348 0.0243201 0.0520065 5 +5.401406e-01 +1.210917e-01 -1.127979e+00 +2.034264e-33
 * $ 0.667402 0.37 83.3999 0.1 191.4 0.13348 0.0243201 0.0520065 5 -4.151277e-02 +3.640322e-01 -1.112785e+00 +2.046301e-34
 * $ 0.667402 0.37 83.3999 0.1 191.4 0.13348 0.0243201 0.0520065 5 -2.855319e-02 -1.974611e-01 -9.697772e-01 +2.407412e-34
 * $ 0.667402 0.37 83.3999 0.1 191.4 0.13348 0.0243201 0.0520065 5 -1.607399e-02 -7.570055e-02 -1.943461e-01 +1.203706e-35
 */
void normalizep(double *b,
                double *te,
                double *tg,
                double *s,
                const double *x,
                const double *y,
                const unsigned int nimg) {
#ifdef _NORMALIZATION_USE_LAST                
  const unsigned int inorm = nimg - 1;
#else
  const unsigned int inorm = 1;
#endif
  const double dx    = x[inorm] - x[0];
  const double dy    = y[inorm] - y[0];
  const double r     = sqrt( dx * dx + dy * dy );
  const double theta = atan2(dy, dx)*180./M_PI;

  *b /= r;
  *s /= r;
  *te -= theta;
  *tg -= theta;
}

/**
 * Normalize the lensed images positions
 */
void normalize(double *xs,
               double *ys,
               double *xd,
               double *yd,
               double *x,
               double *y,
               double *mag,
               const unsigned int nimg) {
#ifdef _NORMALIZATION_USE_LAST                
  const unsigned int inorm = nimg - 1;
#else
  const unsigned int inorm = 1;
#endif
  unsigned int i;

  // Re-center the lensed images, source position and deflector position on
  // the brightest image (i.e. x[0], y[0])
  {
    const double x0 = x[0];
    const double y0 = y[0];
    
    *xs -= x0; *ys -= y0;
    *xd -= x0; *yd -= y0;
    for(i = 0; i < nimg; i++) {
      x[i] -= x0; y[i] -= y0;
    }
  }

  // Rotate all images so as to have x[inorm] = ( x[inorm]^2 + y[inorm]^2 ) and y[inorm] = 0
  // Beware that these rotations comes along with a scaling of sqrt( x[inorm]^2 + y[inorm]^2 )
  {
    const double xrot = x[inorm];
    const double yrot = y[inorm];
    
    rotate(xs, ys, xrot, yrot);
    rotate(xd, yd, xrot, yrot); 
    for(i = 1; i < nimg; i++)
      rotate(&x[i], &y[i], xrot, yrot);
  }

  // Scale the system by x[inorm] so as to have x[inorm] = 1, y[inorm] = 0
  {
    const double r = x[inorm];
    
    *xs /= r; *ys /= r;
    *xd /= r; *yd /= r;
    for(i = 1; i < nimg; i++) {
      x[i] /= r; y[i] /= r;
    }
  }

  // Set the magnitude of the brightest image to zero
  {
  const double m0 = mag[0];
  
  for(i = 0; i < nimg; i++)
    mag[i] -= m0;
  }
}

// Compute the maximal separation between any two lensed images
double max_dist(const double *x, const double *y, const unsigned int nimg) {
  double maxd2 = 0, d2, dx, dy;
  unsigned int i,j;

  for(i = 0; i < nimg; i++) {
    for(j = i + 1; j < nimg; j++) {
      dx = x[i] - x[j];
      dy = y[i] - y[j];
      d2 = dx * dx + dy * dy;
      if(maxd2 < d2)
        maxd2 = d2;
    }
  }

  return sqrt(maxd2);
}

// Output a lens configuration to the standard output
void output_lens( char class,
                  double b,
                  double e,
                  double te,
                  double g,
                  double tg,
                  double s,
                  double xs,
                  double ys,
                  double xd,
                  double yd,
                  double *x,
                  double *y,
                  double *mag,
                  const unsigned int nimg) {
  poutput_node_t *node = poutput.head;

  while(node) {
    switch(node->type) {
      case OUTPUT_CLASS:
        fprintf(stdout, "%d", (int) class);
        break;
      case OUTPUT_B:
        fprintf(stdout, "%+.16e", b);
        break;
      case OUTPUT_E:
        fprintf(stdout, "%+.16e", e);
        break;
      case OUTPUT_TE:
        fprintf(stdout, "%+.16e", te);
        break;
      case OUTPUT_G:
        fprintf(stdout, "%+.16e", g);
        break;
      case OUTPUT_TG:
        fprintf(stdout, "%+.16e", tg);
        break;
      case OUTPUT_S:
        fprintf(stdout, "%+.16e", s);
        break;
      case OUTPUT_XS:
        fprintf(stdout, "%+.16e", xs);
        break;
      case OUTPUT_YS:
        fprintf(stdout, "%+.16e", ys);
        break;
      case OUTPUT_XD:
        fprintf(stdout, "%+.16e", xd);
        break;
      case OUTPUT_YD:
        fprintf(stdout, "%+.16e", yd);
        break;
      case OUTPUT_X:
        fprintf(stdout, "%+.16e", x[node->idx]);
        break;
      case OUTPUT_Y:
        fprintf(stdout, "%+.16e", y[node->idx]);
        break;
      case OUTPUT_MAG:
        fprintf(stdout, "%+.16e", mag[node->idx]);
        break;
      case OUTPUT_NIMG:
        fprintf(stdout, "%u", nimg);
        break;
    }
    node = node->next;
    if(node)
      fprintf(stdout, " ");
  }
  fprintf(stdout, "\n");
}

// Read and convert the next double contained in s
int nextdouble(double *d, const char **s) {
  char *endptr;
  *d = strtod(*s, &endptr);
  if(*s == endptr) {
    fprintf(stderr, "Unable to read double from line \"%s\"\n", *s);
    return -1;
  }
  *s = endptr;
  return 0;
}

// Read and convert the next unsigned int contained in s
int nextuint(unsigned int *ui, const char **s) {
  char *endptr;
  *ui = strtoul(*s, &endptr, 10);
  if(*s == endptr) {
    fprintf(stderr, "Unable to read unsigned int from line \"%s\"\n", *s);
    return -1;
  }
  *s = endptr;
  return 0;
}

int nextline(double *b,  double *e, double *te, double *g, double *tg, double *s,
    double *xs, double *ys, unsigned int *nimg, double *x, double *y, double *mag,
    double *sdf, const char *line) {
  if(obs) { // This is an observation in the form "nimg x y mag"
    *b = 1;
    *e = *te = *g = *tg = *s = *xs = *ys = *sdf = 0;
    if(   nextuint(nimg,  &line)
       || nextdouble(x, &line)
       || nextdouble(y, &line)
       || nextdouble(mag, &line))
      return -1;
  } else if(nextdouble(b, &line) // This is a simulation in the form "b e te g tg s xs ys nimg x y mag sdf"
     || nextdouble(e,   &line)
     || nextdouble(te,  &line)
     || nextdouble(g,   &line)
     || nextdouble(tg,  &line)
     || nextdouble(s,   &line)
     || nextdouble(xs,  &line)
     || nextdouble(ys,  &line)
     || nextuint(  nimg,&line)
     || nextdouble(x, &line)
     || nextdouble(y, &line)
     || nextdouble(mag, &line)
     ||nextdouble( sdf, &line))
    return -1;
  if(*line != '\0')
    return -1;
  return 0;
}

// Return the next lens configuration having a sufficient number of images
unsigned int nextconf(FILE *in, char **line, size_t *nline, size_t *nxymag,
                      double *b,  double *e, double *te, double *g, double *tg, double *s,
                      double *xs, double *ys, double **x, double **y, double **mag,
                      double *sdf) {
  double br, er, ter, gr, tgr, sr, xsr, ysr, xr, yr, magr, sdfr;
  unsigned int nimgr, nimg;
  unsigned int i;
  size_t read;

  // Parse the current input line
  if(nextline(&br, &er, &ter, &gr, &tgr, &sr, &xsr, &ysr, &nimgr, &xr, &yr, &magr, &sdfr, *line)) {
    fprintf(stderr, "Unable to parse line \"%s\"\n", *line);
    return -1U;
  }

  // Skip configurations having an insufficient number of images
  while(nimgr <= pinput.idx[pinput.n-1]) {
    // Read the next nimgr lines
    while(nimgr--) {
      if((read = getline(line, nline, in)) == (size_t) -1)
        return 0; // Return 0 if we reached the EOF
    }
    // Parse the last read line
    if((*line)[read-1] == '\n') (*line)[read-1] = '\0';
    if(nextline(&br, &er, &ter, &gr, &tgr, &sr, &xsr, &ysr, &nimgr, &xr, &yr, &magr, &sdfr, *line)) {
      fprintf(stderr, "Unable to parse line \"%s\"\n", *line);
      return -1U;
    }
  }

  // Set the b, e, te, g, tg, s, xs, ys, sdf and nimg parameters
  *b   = br;
  *e   = er;
  *te  = ter;
  *g   = gr;
  *tg  = tgr;
  *s   = sr;
  *xs  = xsr;
  *ys  = ysr;
  *sdf = sdfr;
  nimg = nimgr;

  // Only reallocate the x, y and mag array if their size is lower than nimg
  if(*nxymag < nimg) {
    // Allocate the new arrays
    *x   = (double *) realloc(*x,   nimg * sizeof(double));
    *y   = (double *) realloc(*y,   nimg * sizeof(double));
    *mag = (double *) realloc(*mag, nimg * sizeof(double));
    if(!x || !y || !mag) {
      fprintf(stderr, "Unable to reallocate the x, y and mag arrays\n");
      return -1U;
    }
    // The x, y and mag array now have a size of nimg elements
    *nxymag = nimg;
  }

  // Fill the x, y and mag arrays
  (*x)[0] = xr; (*y)[0] = yr; (*mag)[0] = magr;
  for(i = 1; i < nimg; i++) {
    // Get the next line from stdin
    if((read = getline(line, nline, in)) == (size_t) -1) {
      fprintf(stderr, "Configuration b=%g e=%g te=%g g=%g tg=%g s=%g xs=%g ys=%g have %u missing lensed images\n",
                       *b, *e, *te, *g, *tg, *s, *xs, *ys, nimg - i);
      return -1U;
    }
    // Parse this line
    if((*line)[read-1] == '\n') (*line)[read-1] = '\0';
    if(nextline(&br, &er, &ter, &gr, &tgr, &sr, &xsr, &ysr, &nimgr, &xr, &yr, &magr, &sdfr, *line)) {
      fprintf(stderr, "Unable to parse line \"%s\"\n", *line);
      return -1U;
    }
    // Check that the last parsed line matches the current lens configuration
    if(*b != br || *e != er || *te != ter || *g != gr || *tg != tgr || *s != sr || *xs != xsr || *ys != ysr || nimg != nimgr) {
      fprintf(stderr, "Configuration b=%g e=%g te=%g g=%g tg=%g s=%g xs=%g ys=%g have %u missing lensed images\n",
                      *b, *e, *te, *g, *tg, *s, *xs, *ys, nimg - i);
      return -1U;
    }
    // Add this line to the x, y and mag array
    (*x)[i] = xr; (*y)[i] = yr; (*mag)[i] = magr;
  }
  
  // If x, y coordinates are spherical coordinates expressed in radian then scale
  // x by cos(y1)
  if(sph) {
    for(i = 0; i < nimg; i++)
      (*x)[i] *= cos((*y)[0]);
  }

  return nimg;
}

// Compare two double contained at arg[a] and arg[b].
// This is used in order to obtain the indices that allow to sort the arg array
int argsortd(const void *a, const void *b, void *arg) {
  const double *mag = (const double *) arg;
  const int ia = *((const int *) a);
  const int ib = *((const int *) b);
  return (mag[ia] < mag[ib]) ? -1 : 1;
}

// Reorder the array a based on the set of (unique) indices contained in idx
void reorder(double *a, const unsigned int *idx, const unsigned int n) {
  unsigned int i, j;
  double tmp;

  for(i = 0; i < n; i++) {
    for(j = idx[i]; j < i; j = idx[j])
      ;
    if(i != j) {
      tmp = a[i];
      a[i] = a[j];
      a[j] = tmp;
    }
  }
}

// Draw a normally distributed variable with mean 0 and variance 1
double randn(void) {
  static char iset = 0;
  static double next = 0;
  double x1, x2, r2, v;
  if(iset == 0) {
    do {
      x1 = 2. / RAND_MAX * rand() - 1;
      x2 = 2. / RAND_MAX * rand() - 1;
      r2 = x1 * x1 + x2 * x2;
    } while(r2 >= 1.0 || r2 == 0.0);
    v = sqrt( -2. * log(r2) / r2);
    next = x2 * v;
    iset = 1;
    return x1 * v;
  } else {
    iset = 0;
    return next;
  }
}

// Core function for producing the dataset
int build(void) {
  // Variables of the lensed system
  double b, e, te, g, tg, s, xs, ys, xd, yd, *x = NULL, *y = NULL, *mag = NULL, sdf;
  unsigned int nimg;

  // Variables of nextconf and getline
  char *line     = NULL;
  size_t nline   = 0;
  size_t nxymag  = 0;
  size_t read;

  // General usage counters
  unsigned int i, j;

  // Initialize the random number generator
  srand(seed);

  // Read the initial line
  if((read = getline(&line, &nline, stdin)) == (size_t) -1)
    return -1;

  // If this is a header line, skip it
  if(line[read-1] == '\n') line[read-1] = '\0';
  if(strncmp(line, "b e te g tg s xs ys nimg x y mag", strlen("b e te g tg s xs ys nimg x y mag")) == 0
     || strncmp(line, "nimg x y mag", strlen("nimg x y mag")) == 0) {
    if((read = getline(&line, &nline, stdin)) != (size_t) -1) {
      if(line[read-1] == '\n') line[read-1] = '\0';
    }
  }

  // The deflector stand at position (0,0)
  xd = yd = 0;

  // While some valid configuration remain
  while (read != (size_t) -1 && (nimg = nextconf(stdin, &line, &nline, &nxymag, &b, &e, &te, &g, &tg, &s, &xs, &ys, &x, &y, &mag, &sdf)) > 0) {
    // If the number of images we have is equal to -1U, then nextconf encountered an error
    if(nimg == -1U) {
      free(x); free(y); free(mag);free(line);
      return -1;
    }

    // If some relative noise will be added to the lensed image positions,
    // then compute the maximal separation between the whole lensed images
    double snoise = 1.; // The relative noise scaling that should be applied (--nrel option)
    if(nrel && (xynoise > 0))
      snoise = max_dist(x,y,nimg);

    // Perform nnoise realization
    for(i = 0; i < nnoise; i++) {
      // The noisy version of the lens parameters
      double bn, en, ten, gn, tgn, sn, xsn, ysn, xdn, ydn;
      double xn[nimg];
      double yn[nimg];
      double magn[nimg];

      // Copy the configuration parameters
      bn  = b;
      en  = e;
      ten = te;
      gn  = g;
      tgn = tg;
      sn  = s;
      xsn = xs;
      ysn = ys;
      xdn = xd;
      ydn = yd;

      // Add noise to this configuration
      for(j = 0; j < nimg; j++)  {
        xn[j]   = x[j]   + xynoise   * snoise * randn();
        yn[j]   = y[j]   + xynoise   * snoise * randn();
        magn[j] = mag[j] + magnoise  * randn();
      }

      { // Sort images according to their magnitudes (brightest first)
        unsigned int is[nimg];
        for(j = 0; j < nimg; j++) is[j] = j;
        qsort_r(&is, nimg, sizeof(int), argsortd, magn);
        reorder(magn, is, nimg);
        reorder(xn, is, nimg);
        reorder(yn, is, nimg);
      }

      // Only keep those images specified by the --input option
      for(j = 0; j < pinput.n; j++) {
        xn[j]   = xn[pinput.idx[j]];
        yn[j]   = yn[pinput.idx[j]];
        magn[j] = magn[pinput.idx[j]];
      }

      // Normalize this configuration
      if(norm) {
        // Normalize the lens parameters, if needed
        if(   (poutput.types & OUTPUT_B)
           || (poutput.types & OUTPUT_S)
           || (poutput.types & OUTPUT_TE)
           || (poutput.types & OUTPUT_TG))
          normalizep(&bn, &ten, &tgn, &sn, xn, yn, pinput.n);
        // Normalize the source, deflector and lensed images positions
        normalize(&xsn, &ysn, &xdn, &ydn, xn, yn, magn, pinput.n);
      }

      // Output this configuration
      output_lens((obs == 0), bn, en, ten, gn, tgn, sn, xsn, ysn, xdn, ydn, xn, yn, magn, pinput.n);
    }

    // Read next line
    if((read = getline(&line, &nline, stdin)) != (size_t) -1) {
      if(line[read-1] == '\n') line[read-1] = '\0';
    }
  }

  // Delete arrays and input line
  free(x); free(y); free(mag);free(line);

  return 0;
}

/*************************************************************************************************
 * Main program
 ************************************************************************************************/
int main(int argc, char **argv) {
  int opt;

  // Get the program parameters
  while ((opt = getopt_long(argc, argv, options_short, options_long, NULL)) != -1) {
    switch (opt) {
      case 'i':            /* --input   */
        input  = optarg;
        break;
      case 'o':            /* --output   */
        output = optarg;
        break;
      case 'O':            /* --obs      */
        obs = 1;
        break;
      case PARAM_SPH:      /* --sph      */
        sph = 1;
        break;
      case 'n':            /* --norm     */
        norm = 1;
        break;
      case 's':            /* --seed     */
        seed = strtoul(optarg, NULL, 10);
        break;
      case PARAM_NNOISE:   /* --nnoise   */
        nnoise = strtoul(optarg, NULL, 10);
        break;
      case PARAM_XYNOISE:  /* --xynoise  */
        xynoise = strtod(optarg, NULL);
        break;
      case 'r':            /* --nrel     */
        nrel = 1;
        break;
      case PARAM_MAGNOISE: /* --magnoise */
        magnoise = strtod(optarg, NULL);
        break;
      case PARAM_NOHEADER: /* --noheader */
        noheader = 1;
        break;
      case '?':            /* --help     */
        printHelp(stdout);
        exit(EXIT_SUCCESS);
        break;
      default:
        printHelp(stderr);
        exit(EXIT_FAILURE);
        break;
    }
  }

  // Parse the --input and --output options
  parse_input_option();
  parse_output_option();
  if(check_options()) {
    free(pinput.idx);
    poutput_destroy();
    exit(EXIT_FAILURE);
  }

  // Print the header line if asked to do so
  if(noheader == 0)
    print_header();

  // Build the learning set of observations and output it to standard output
  if(build()) {
    fprintf(stderr, "An error occurred during the building of the learning set of observations\n");
    free(pinput.idx);
    poutput_destroy();
    exit(EXIT_FAILURE);
  }

  // Delete the parsed input and output options
  free(pinput.idx);
  poutput_destroy();

  // Exit
  exit(EXIT_SUCCESS);
}
