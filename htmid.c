/**
 * @file htmid.c
 * 
 * Hierarchical Triangular Mesh identifier computation program.
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
#include <errno.h>
#include <getopt.h>
#include <htm.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*************************************************************************************************
 * Program parameters
 ************************************************************************************************/
#define DEFAULT_LEVEL 10U

unsigned int level = DEFAULT_LEVEL;

static struct option options_long[] = {
/*-l*/{ "level", required_argument, 0, 'l' },
/*-?*/{ "help", no_argument, 0, '?' },
/*  */{ 0, 0, 0, 0 } };

static const char options_short[] = "l:?";

/*************************************************************************************************
 * Help message
 ************************************************************************************************/
static const char
    *help_message =
        "USAGE:\n"
        "  htmid [options]\n"
        "\n"
        "DESCRIPTION:\n"
        "  Compute the HTM identifiers of object on the celestial sphere\n"
        "\n"
        "OPTIONS:\n"
        "  -l, --level,   HTM level to use [default = %u]\n"
        "  -?, --help,    Print this help message\n"
        "\n"
        "OUTPUT FORMAT:\n"
        "  All lines from standard input will be prepend by the HTM identifier derived from the two first columns\n"
        "\n"
        "NOTES:\n"
        "  Right ascension and declination should be contained as the two first columns of file\n"
        "  By default, right ascension and declination are considered to be expressed in radian, if followed by 'd' or 'D', then\n"
        "  these are considered to be in degree whereas if right ascension is in the form \"hh:mm:ss\", sexagecimal coordinates\n"
        "  are supposed corresponding to hour, minute and second. Declination in the form \"dd:mm:ss\" are understood as\n"
        "  sexagecimal coordinates in degree, arcmin and arcsec.\n";

void printHelp(FILE *out) {
  fprintf(out, help_message, DEFAULT_LEVEL);
}

/*************************************************************************************************
 * Functions
 ************************************************************************************************/
int process(FILE *in) {
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  double ra,dec;
  double v[3];
  char *endptr;
  char sid[level+3];

  // Browse all entries from the input file
  while((read = getline(&line, &len, in)) != -1) {

    // Try to extract RA
    if(str2ra(&ra,line,&endptr)) {
      fprintf(stderr, "Unable to get right ascension from line \"%s\"\n", line);
      free(line);
      return -1;
    }

    // Try to extract DEC
    if(str2dec(&dec,endptr,NULL)) {
      fprintf(stderr, "Unable to get declination from line \"%s\"\n", line);
      free(line);
      return -1;
    }

    // Process this object
    htm_cart_coord(v, ra, dec);
    fprintf(stdout, "%s %s", htm_string(sid, htm_get(v,level)), line);
  }

  // Free the array used in getline
  free(line);
  
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
      case 'l':
        level = strtoul(optarg, NULL, 10);
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

  // Process file
  if(process(stdin))
    exit(EXIT_FAILURE);
  
  // Exit
  exit(EXIT_SUCCESS);
}
