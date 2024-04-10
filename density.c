/**
 * @file density.c
 * 
 * Compute the density of object on the celestial sphere based on triangles
 * from the Hirarchical Triangular Mesh.
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
#include <math.h>
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
        "  density [options]\n"
        "\n"
        "DESCRIPTION:\n"
        "  Compute the density of object on the celestial sphere based on triangles\n"
        "  from the Hirarchical Triangular Mesh\n"
        "\n"
        "OPTIONS:\n"
        "  -l, --level,   HTM level to use [default = %u]\n"
        "  -?, --help,    Print this help message\n"
        "\n"
        "OUTPUT FORMAT:\n"
        "  \"htmd_id ra dec area density\"\n"
        "  where \n"
        "    htm_id is the identifier for the returned HTM triangle\n"
        "    ra, dec are the mean position in the HTM triangle [rad]\n"
        "    area is the solid angle subtended by the HTM triangle [sr]\n"
        "    density in number of object per steradian found in the triangle\n"
        "NOTES:\n"
        "  Right ascension and declination should be contained as the two first columns of both pos\n"
        "  By default, right ascension and declination are considered to be expressed in radian, if followed by 'd' or 'D', then\n"
        "  these are considered to be in degree whereas if right ascension is in the form \"hh:mm:ss\", sexagecimal coordinates\n"
        "  are supposed corresponding to hour, minute and second. Declination in the form \"dd:mm:ss\" are understood as\n"
        "  sexagecimal coordinates in degree, arcmin and arcsec.\n"
        "\n"
        "  Input is taken from standard input\n";

void printHelp(FILE *out) {
  fprintf(out, help_message, DEFAULT_LEVEL);
}

/*************************************************************************************************
 * Functions
 ************************************************************************************************/
unsigned long *count_object(unsigned long long *n, FILE *in) {
  const unsigned long long mask = ( 0x08ULL << ( 2 * level ) ) - 1ULL;
  unsigned long *nobj;
  double ra,dec;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  double v[3];
  char *endptr;

  // Allocate the array containing the number of objects
  *n = 0x08UL << ( 2 * level);
  nobj = (unsigned long *) calloc(*n, sizeof(unsigned long));
  if(!nobj) {
    fprintf(stderr, "Unable to allocate the array of objects\n");
    return NULL;
  }

  // Browse all entries from the input file
  while((read = getline(&line, &len, in)) != -1) {

    // Try to extract RA
    if(str2ra(&ra,line,&endptr)) {
      fprintf(stderr, "Unable to get right ascension from line \"%s\"\n", line);
      free(line);free(nobj);
      return NULL;
    }

    // Try to extract DEC
    if(str2dec(&dec,endptr,NULL)) {
      fprintf(stderr, "Unable to get declination from line \"%s\"\n", line);
      free(line);free(nobj);
      return NULL;
    }

    // Remove the '\n' end character
    if(line[read-1] == '\n')
      line[read-1] = '\0';

    // Process this object
    htm_cart_coord(v, ra, dec);
    nobj[htm_get(v,level) & mask]++;
  }

  // Free the array used in getline
  free(line);

  return nobj;
}

#define dot_prod(u,v) (u[0]*v[0]+u[1]*v[1]+u[2]*v[2])

void normalv(double *r, const double *u, const double *v) {
  double n;
  // Cross product u X v = |u| |v| sin(theta) n
  r[0]=u[1]*v[2]-u[2]*v[1];
  r[1]=u[2]*v[0]-u[0]*v[2];
  r[2]=u[0]*v[1]-u[1]*v[0];
  // Normalize cross product r = n
  n = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
  r[0] /= n;
  r[1] /= n;
  r[2] /= n;
}

double triangle_area(const double *v0, const double *v1, const double *v2) {
  double n01[3], n02[3], n12[3];
  double A, B, C;
  normalv(n01,v0,v1);
  normalv(n02,v0,v2);
  normalv(n12,v1,v2);
  A = acos(dot_prod(n01,n02));
  B = acos(dot_prod(n02,n12));
  C = acos(-dot_prod(n01,n12));
  return A + B + C - M_PI;
}

void report_densities(const unsigned long *nobj, const unsigned long long n) {
  char s[level+3], sra[13], sdec[14];
  double ra, dec, S, v0[3], v1[3], v2[3];
  unsigned long long i;

  for(i = 0; i < n; i++) {
    // Get the HTM vertices for this id
    htm_vertices(v0,v1,v2, n + i);

    // If this HTM id contains objects
    if(nobj[i]) {
      double vm[3];
      double nvm;

      // Compute the mean vertex
      vm[0]=v0[0]+v1[0]+v2[0];
      vm[1]=v0[1]+v1[1]+v2[1];
      vm[2]=v0[2]+v1[2]+v2[2];
      nvm = sqrt(dot_prod(vm,vm));
      vm[0]/=nvm; vm[1]/=nvm; vm[2]/=nvm;
      htm_sphere_coord(&ra,&dec,vm);

      // Compute the triangle area
      S = triangle_area(v0,v1,v2);

      // Print result
      fprintf(stdout, "%s %s %s %g %g\n", htm_string(s, n + i), ra2str(sra, ra), dec2str(sdec,dec), S, 1./S*nobj[i]);
    }
  }
}

/*************************************************************************************************
 * Main program
 ************************************************************************************************/
int main(int argc, char **argv) {
  unsigned long *nobj;
  unsigned long long n;
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

  // Count objects
  nobj = count_object(&n, stdin);
  if(!nobj)
    exit(EXIT_FAILURE);
  
  // Browse counts
  report_densities(nobj, n);

  // Free the array of objects
  free(nobj);

  // Exit
  exit(EXIT_SUCCESS);
}
