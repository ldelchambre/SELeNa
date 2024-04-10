/**
 * @file pconvert.c
 * 
 * Astronomical coordinates converter.
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
#include <getopt.h>
#include <htm.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*************************************************************************************************
 * Program parameters
 ************************************************************************************************/
#define DEFAULT_PREC 6

int prec = DEFAULT_PREC;

static struct option options_long[] = {
/*-p*/{ "prec", required_argument, 0, 'p' },
/*  */{ 0, 0, 0, 0 }
};

static const char options_short[] = "p:?";

/*************************************************************************************************
 * Help message
 ************************************************************************************************/
static const char
    *help_message =
        "USAGE:\n"
        "  convert [options] [conversion]\n"
        "\n"
        "DESCRIPTION:\n"
        "  Astronomical coordinates converter.\n"
        "\n"
        "  Available conversions are:\n"
        "    \"sph\"     : Spĥerical coordinates expressed in radian [default]\n"
        "    \"sphd\"    : Spĥerical coordinates expressed in degree\n"
        "    \"sex\"     : Sexagesimal coordinates\n"
        "    \"cart\"    : Cartesian coordinates\n"
        "    \"icrs2gal\": Galactic coordinates from ICRS coordinates (\"sph\" format)\n"
        "    \"gal2icrs\": ICRS coordintes from galactic coordinates (\"sph\" format)\n"
        "\n"
        "OPTIONS:\n"
        "  -p, --prec, The decimal precision to used [default = %d]\n"
        "  -?, --help, Print this help message\n"
        "\n"
        "NOTES:\n"
        "  Input coordinates are read from standard input and outputed to standard output\n"
        "\n"
        "  Right ascension and declination should be contained as the two first columns of standard input\n"
        "  By default, right ascension and declination are considered to be expressed in radian, if followed by 'd' or 'D', then\n"
        "  these are considered to be in degree whereas if right ascension is in the form \"hh:mm:ss\", sexagesimal coordinates\n"
        "  are supposed corresponding to hour, minute and second. Declination in the form \"dd:mm:ss\" are understood as\n"
        "  sexagecimal coordinates in degree, arcmin and arcsec.\n"
        "\n";

void printHelp(FILE *out) {
  fprintf(out, help_message, DEFAULT_PREC);
}

/*******************************************************************************
 * Program functions
 ******************************************************************************/
void sph_process(const double ra, const double dec) {
  fprintf(stdout, "%.*f %.*f", prec, ra, prec, dec);
}

void sphd_process(const double ra, const double dec) {
  fprintf(stdout, "%.*f %.*f", prec, 180.*ra/M_PI, prec, 180.*dec/M_PI);
}

void sex_process(const double ra, const double dec) {
  char sra[10+prec];
  char sdec[11+prec];

  ra2strn(sra, ra, prec);
  dec2strn(sdec, dec, prec);

  fprintf(stdout, "%s %s", sra, sdec);
}

void cart_process(const double ra, const double dec) {
  double v[3];
  htm_cart_coord(v, ra, dec);
  fprintf(stdout, "%.*f %.*f %.*f", prec, v[0], prec, v[1], prec, v[2]);
}

void process_transform(const double ra0, const double dec0, const double r[3][3]) {
  // const double (*r)[3] = (const double (*)[3]) t;
  unsigned char i, j;
  double v0[3], v1[3];
  double ra1, dec1;

  htm_cart_coord(v0, ra0, dec0);
  for (i = 0; i < 3; i++) {
    v1[i] = 0.;
    for (j = 0; j < 3; j++)
      v1[i] += r[i][j] * v0[j];
  }
  htm_sphere_coord(&ra1, &dec1, v1);
  fprintf(stdout, "%.*f %.*f", prec, ra1, prec, dec1);
}

void icrs2gal_process(const double ra, const double dec) {
  // Matrix taken from the SOFA release 2018-01-30, icrs2g function
  // Copyright (C) 2018 IAU SOFA Board.
  const double r[3][3] = {
    {-0.054875560416215368492398900454,-0.873437090234885048760383168409,-0.483835015548713226831774175116},
    {+0.494109427875583673525222371358,-0.444829629960011178146614061616,+0.746982244497218890527388004556},
    {-0.867666149019004701181616534570,-0.198076373431201528180486091412,+0.455983776175066922272100478348}
  };
  process_transform(ra, dec, r);
}

void gal2icrs_process(const double ra, const double dec) {
  // Matrix taken from the SOFA release 2018-01-30, icrs2g function
  // Copyright (C) 2018 IAU SOFA Board.
  const double r[3][3] = {
    {-0.054875560416215368492398900454,+0.494109427875583673525222371358,-0.867666149019004701181616534570},
    {-0.873437090234885048760383168409,-0.444829629960011178146614061616,-0.198076373431201528180486091412},
    {-0.483835015548713226831774175116,+0.746982244497218890527388004556,+0.455983776175066922272100478348 }
  };
  process_transform(ra, dec, r);
}

/*************************************************************************************************
 * Main program
 ************************************************************************************************/
int main(int argc, char **argv) {
  // Object position
  double ra, dec;

  // The processing function to call
  void (*process)(const double, const double) = sph_process;

  // End of converted string
  char *endptr;

  // Variables of readline
  char *line = NULL;
  size_t len = 0;
  ssize_t read;

  // Current option
  int opt;

  // Get the program parameters
  while ((opt = getopt_long(argc, argv, options_short, options_long, NULL)) != -1) {
    switch (opt) {
      case 'p':
        prec = strtoul(optarg, NULL, 10);
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

  // Get the conversion to perform
  if(optind < argc) {
    if     (strcmp(argv[optind], "sph")      == 0) process = sph_process;
    else if(strcmp(argv[optind], "sphd")     == 0) process = sphd_process;
    else if(strcmp(argv[optind], "sex")      == 0) process = sex_process;
    else if(strcmp(argv[optind], "cart")     == 0) process = cart_process;
    else if(strcmp(argv[optind], "icrs2gal") == 0) process = icrs2gal_process;
    else if(strcmp(argv[optind], "gal2icrs") == 0) process = gal2icrs_process;
    else {
      fprintf(stderr, "Unrecognized conversion \"%s\"\n", argv[optind]);
      exit(EXIT_FAILURE);
    }
    optind++;
  }

  // Check for extra arguments
  if(optind < argc) {
    fprintf(stderr, "Warning: extra argument(s) found (");
    while(optind < argc)
      fprintf(stderr, " %s", argv[optind++]);
    fprintf(stderr, ")\n");
  }

  // Browse each line of file
  while((read = getline(&line, &len, stdin)) != -1) {
    // Try to extract RA
    if(str2ra(&ra,line,&endptr)) {
      fprintf(stderr, "Unable to get right ascension from line \"%s\"\n", line);
      free(line);
      exit(EXIT_FAILURE);
    }

    // Try to extract DEC
    if(str2dec(&dec,endptr,&endptr)) {
      fprintf(stderr, "Unable to get declination from line \"%s\"\n", line);
      free(line);
      exit(EXIT_FAILURE);
    }

    // Process this entry
    process(ra,dec);

    // Print the rest of the line
    fprintf(stdout, "%s", endptr);
  }

  // Delete the input line
  free(line);

  return 0;
}
