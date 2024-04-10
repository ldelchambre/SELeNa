/**
 * @file astrostr_test.c
 * 
 * Astronomical positions and angle string manipulation, test file.
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NTEST (1U << 16)

#define ARCSEC (M_PI/648000)
#define SECOND (M_PI/43200)

int main(void) {
  const char *units[] = {
    /*  0 */ "degree",
    /*  1 */ "arcmin",
    /*  2 */ "arcsec",
    /*  3 */ "mas",
    /*  4 */ "hour",
    /*  5 */ "minute",
    /*  6 */ "second",
    /*  7 */ "d",
    /*  8 */ "h",
    /*  9 */ "m",
    /* 10 */ "s",
    /* 11 */ "'",
    /* 12 */ "\""
  };
  double ra, dec, theta;
  double ra2, dec2, theta2;
  char s[256], *send;
  unsigned int ndec, iunit;
  unsigned int i; 
  
  for(i = 0; i < NTEST; i++) {
    // Initialize the random number generator
    srand(1 + i);
    
    // Pick up a random angle, right ascension and declination
    theta = 2. * rand() / RAND_MAX - 1; // [-1, 1]
    ra    = 2. * M_PI * rand() / (1. + RAND_MAX); // [0, 2*pi[
    dec   = 0.5 * M_PI * ( 2. * rand() / RAND_MAX - 1 ); // [-pi/2, pi/2]
    ndec  = rand() % 10; // [0, 9]
    iunit = rand() % 13; // [0, 12]
    
    
    // Check the conversion of right ascension from and to sexagesimal coordinates 
    ra2strn(s, ra, ndec);
    if(str2ra(&ra2, s, &send) || fabs(ra-ra2) > SECOND*pow(10,-ndec) || *send != '\0')
      fprintf(stderr, "%u) Failed to convert \"%s\" into right ascension (real: %g, found: %g, error: %g)\n", i, s, ra, ra2, fabs(ra-ra2));
    
    // Check the conversion of declination from and to sexagesimal coordinates 
    dec2strn(s, dec, ndec);
    if(str2dec(&dec2, s, &send) || fabs(dec-dec2) > ARCSEC*pow(10,-ndec) || *send != '\0')
      fprintf(stderr, "%u) Failed to convert \"%s\" into declination (real: %g, found: %g, error: %g)\n", i, s, dec, dec2, fabs(dec-dec2));
    
    // Check the conversion of string to angle
    sprintf(s, "%.*f%s", ndec, theta, units[iunit]);
    switch(iunit) {
      case 0:  theta *= M_PI/180;       break; // "degree"
      case 1:  theta *= M_PI/10800;     break; // "arcmin"
      case 2:  theta *= M_PI/648000;    break; // "arcsec"
      case 3:  theta *= M_PI/648000000; break; // "mas"
      case 4:  theta *= M_PI/12;        break; // "hour"
      case 5:  theta *= M_PI/720;       break; // "minute"
      case 6:  theta *= M_PI/43200;     break; // "second"
      case 7:  theta *= M_PI/180;       break; // "d"
      case 8:  theta *= M_PI/12;        break; // "h"
      case 9:  theta *= M_PI/720;       break; // "m"
      case 10: theta *= M_PI/43200;     break; // "s"
      case 11: theta *= M_PI/10800;     break; // "'"
      case 12: theta *= M_PI/648000;    break; // "\""
    }
    if(str2angle(&theta2, s, &send) || fabs(theta-theta2) > pow(10,-ndec) || *send != '\0')
      fprintf(stderr, "%u) Failed to convert \"%s\" into angle (real: %g, found: %g, error: %g)\n", i, s, theta, theta2, fabs(theta-theta2));
  }
  
  return 0;
}
