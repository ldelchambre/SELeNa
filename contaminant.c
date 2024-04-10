/**
 * @file contaminant.c
 *
 * Simulate gravitational lens configurations based on a non-singular isothermal
 * ellipsoid lens model in presence of an external shear.
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
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*************************************************************************************************
 * Program parameters
 ************************************************************************************************/
#define DEFAULT_NOBS 1000000
#define DEFAULT_NIMG 4
#define DEFAULT_MAG  3.
#define DEFAULT_PREC 6
#define DEFAULT_SEED 1U

unsigned long nobs = DEFAULT_NOBS;
unsigned int  nimg = DEFAULT_NIMG;
double        mag  = DEFAULT_MAG;
unsigned int  prec = DEFAULT_PREC;
unsigned int  seed = DEFAULT_SEED;
char noheader = 0;

#define PARAM_NOHEADER 65536

static struct option options_long[] = {
/*-N*/{ "nobs", required_argument, 0, 'N'},
/*-n*/{ "nimg", required_argument, 0, 'n'},
/*-m*/{ "mag", required_argument, 0, 'm'},
/*-p*/{ "prec", required_argument, 0, 'p'},
/*-s*/{ "seed", required_argument, 0, 's'},
/*  */{ "noheader", no_argument, 0, PARAM_NOHEADER},
/*-?*/{ "help", no_argument, 0, '?' },
/*  */{ 0, 0, 0, 0 } };

static const char options_short[] = "N:n:m:p:s:?";

/*************************************************************************************************
 * Help message
 ************************************************************************************************/
static const char
    *help_message =
        "USAGE:\n"
        "  contaminant [options]\n"
        "\n"
        "DESCRIPTION:\n"
        "  Produce a set of fortuitous clusters of celestial objects in a given magnitude range.\n"
        "\n"
        "  Clusters are uniformly drawn in a circle of radius 1\n"
        "  and within [0, m] regarding their magnitudes where m is specified by the --mag option.\n"
        "\n"
        "OPTIONS:\n"
        "  -N, --nobs,     The number of contaminants observations to produce [default = %lu]\n"
        "  -n, --nimg,     The number of images in each cluster [default = %u]\n"
        "  -m, --mag,      The upper range of magnitude to produce [default = %g]\n"
        "  -p, --prec,     The desired output precision [default=%u]\n"
        "  -s, --seed,     The random seed to use [default = %u].\n"
        "      --noheader, Don't print an header line\n"
        "  -?, --help,     Print this help message\n";

void printHelp(FILE *out) {
  fprintf(out, help_message, DEFAULT_NOBS, DEFAULT_NIMG, DEFAULT_MAG, DEFAULT_PREC, DEFAULT_SEED);
}
/*************************************************************************************************
 * Main function
 ************************************************************************************************/
int build(void) {
  double m, a, r;
  unsigned int i, j;

  // Initialize the random number generator for contaminants
  srand(seed+1);

  // Print the header line
  if(noheader == 0)
    fprintf(stdout, "nimg x y mag\n");

  // Simulate nobs contaminants
  for(i = 0; i < nobs; i++)
    for(j = 0; j < nimg; j++) {
      m = mag * rand() / RAND_MAX;
      a = 2. * M_PI * rand() / ( 1. + RAND_MAX );
      r = sqrt((double) rand() / ( 1. + RAND_MAX ) );
      fprintf(stdout, "%u %+*.*e %+*.*e %+*.*e\n",
                      nimg,
                      (prec == 0) ? 6 : prec+7, prec, r*cos(a),
                      (prec == 0) ? 6 : prec+7, prec, r*sin(a),
                      (prec == 0) ? 6 : prec+7, prec, m);
    }
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
      case 'N':
        nobs = strtoul(optarg, NULL, 10);
        break;
      case 'n':
        nimg = strtoul(optarg, NULL, 10);
        break;
      case 'm':
        mag = strtod(optarg, NULL);
        break;
      case 'p':
        prec = strtoul(optarg, NULL, 10);
        break;
      case 's':
        seed = strtoul(optarg, NULL, 10);
        break;
      case PARAM_NOHEADER:
        noheader = 1;
        break;
      case '?':
        printHelp(stdout);
        exit(EXIT_SUCCESS);
        break;
      default:
        fprintf(stderr, "Unrecognized options \"%s\"\n", argv[optind]);
        printHelp(stderr);
        exit(EXIT_FAILURE);
        break;
    }
  }

  // Build the contaminants
  build();

  exit(EXIT_SUCCESS);
}
