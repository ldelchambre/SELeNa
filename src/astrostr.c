/**
 * @file astrostr.c
 *
 * Astronomical positions and angle string manipulation, source file.
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#ifndef M_PI
// Ensure M_PI if defined (not always the case in math.h)
#define M_PI 3.14159265358979323846
#endif

/**
 * Conversion factor from degree to radian, 1 degree = DEGREE rad
 */
#define DEGREE (M_PI/180.)

/**
 * Conversion factor from arcmin to radian, 1 arcmin = ARCMIN rad
 */
#define ARCMIN (DEGREE/60.)

/**
 * Conversion factor from arcsec to radian, 1 arcsec = ARCSEC rad
 */
#define ARCSEC (ARCMIN/60.)

/**
 * Conversion factor from hour to radian, 1 hour  = HOUR rad
 */
#define HOUR   (M_PI/12.)

/**
 * Conversion factor from minute to radian, 1 minute = MINUTE rad
 */
#define MINUTE (HOUR/60.)

/**
 * Conversion factor from second to radian, 1 second = SECOND rad
 */
#define SECOND (MINUTE/60.)

/**
 * Unit of angle
 */
static const struct {
  char *unit;  // The unit name or symbol
  double conv; // The unit conversion factor
} _astroparse_units[] = {
  {"degree", DEGREE},
  {"arcmin", ARCMIN},
  {"arcsec", ARCSEC},
  {"mas", 1e-3*ARCSEC},
  {"hour", HOUR},
  {"minute", MINUTE},
  {"second", SECOND},
  {"d", DEGREE},
  {"h", HOUR},
  {"m", MINUTE},
  {"s", SECOND},
  {"'", ARCMIN},
  {"\"", ARCSEC},
  {"", 0.}
};

int str2ra(double *ra, const char *s, char **send) {
  double tmp;
  char *endptr;
  char neg = 0;

  // Convert right ascension (RA) into radian
  *ra = strtod(s,&endptr);

  // If nothing was converted, return with error
  if(s == endptr)
    return -1;

  // Depending on the value of the last character
  switch(*endptr) {

    // If s match xd or Xd, then RA is expressed in degree
    case 'd':
    case 'D':
      *ra *= DEGREE;
      endptr++;
      break;

    // If s match hh:mm:ss, then RA is expressed in sexagecimal notation
    case ':':
      // Check whether the right ascension is negative
      neg = signbit(*ra); // Avoid neg = (*ra < 0) as if *ra == -0, this test fails
      if(neg)
        *ra = fabs(*ra);

      // Convert hours into radian
      *ra *= HOUR;
      s = endptr + 1; // Skip ':'

      // Convert minutes into radian
      tmp = MINUTE*strtod(s,&endptr);
      if(signbit(tmp) || s == endptr || *endptr != ':')
        return -1;
      *ra += tmp;
      s = endptr + 1; // Skip ':'

      // Convert seconds into radian
      tmp = SECOND*strtod(s,&endptr);
      if(signbit(tmp) || s == endptr)
        return -1;
      *ra += tmp;

      // Convert take the negative value of ra if it is negative
      if(neg) *ra = copysign(*ra, -1.);
      break;
  }

  // Set send if asked to do so
  if(send)
    *send = endptr;

  // Normalize ra so as to stand within [0, 2*M_PI[
  if(*ra < 0 || 2.*M_PI <= *ra)
    *ra = (*ra > 0) ? fmod(*ra, 2.*M_PI) : 2.*M_PI - fmod(-(*ra), 2.*M_PI);

  return 0;
}

int str2dec(double *dec, const char *s, char **send) {
  double tmp;
  char *endptr;
  char neg;

  // Convert declination (DEC) into radian
  *dec = strtod(s,&endptr);

  // If nothing was converted, return with error
  if(s == endptr)
    return -1;

  // Depending on the value of the last character
  switch(*endptr) {
    // If s match xd or Xd, then DEC is expressed in degree
    case 'd':
    case 'D':
      *dec *= DEGREE;
      endptr++;
      break;
    // If s match hh:mm:ss, then DEC is expressed in sexagecimal notation
    case ':':
      // Check whether the declination is negative
      neg = signbit(*dec); // We don't use (*dec < 0) as if *dec == -0., this test fails
      if(neg)
        *dec = fabs(*dec);

      // Convert degree into radian
      *dec *= DEGREE;
      s = endptr + 1; // Skip ':'

      // Convert arcmin into radian
      tmp = ARCMIN*strtod(s,&endptr);
      if(signbit(tmp) || s == endptr || *endptr != ':')
        return -1;
      *dec += tmp;
      s = endptr + 1; // Skip ':'

      // Convert seconds into radian
      tmp = ARCSEC*strtod(s,&endptr);
      if(signbit(tmp) || s == endptr)
        return -1;
      *dec += tmp;

      // Convert take the negative value of dec if it is negative
      if(neg) *dec *= -1;
      break;
  }

  // Set send if asked to do so
  if(send)
    *send = endptr;

  // Normalize dec so as to stand within [-pi/2, pi/2]
  *dec = (*dec > 0) ? fmod(*dec, 2.*M_PI) : 2.*M_PI - fmod(-(*dec), 2.*M_PI);
  if(*dec <= M_PI)
    *dec = (0.5 * M_PI < *dec) ? M_PI - (*dec) : (*dec);
  else
    *dec = (*dec < 1.5 * M_PI) ? M_PI - (*dec) : (*dec) - 2.*M_PI;

  return 0;
}

int str2angle(double *theta, const char *s, char **send) {
  char *endptr;

  // Convert right ascension (RA) into radian
  *theta = strtod(s,&endptr);

  // If nothing was converted, return with error
  if(s == endptr)
    return -1;

  // Check if the angle was followed by a unit
  s = endptr;
  if(!isspace(*s) && *s != '\0') {
    unsigned int i;

    // Browse all angles stored in _astroparse_units
    for(i = 0; s == endptr && _astroparse_units[i].conv != 0.; i++) {
      if(strcasecmp(_astroparse_units[i].unit,endptr) == 0) {
        *theta *= _astroparse_units[i].conv;
        endptr += strlen(_astroparse_units[i].unit);
      }
    }
    s = endptr;
  }

  // Set send if asked to do so
  if(send)
    *send = endptr;

  return 0;
}

char *ra2str(char *s, const double ra) {
  return ra2strn(s,ra,3);
}

char *ra2strn(char *s, const double ra, const unsigned int ndec) {
  double ss = ra;
  const double hh = floor(ss/HOUR);
  ss = ss - hh*HOUR; if(ss < 0) ss = 0;
  const double mm = floor(ss/MINUTE);
  ss -= mm*MINUTE; if(ss < 0) ss = 0;
  ss /= SECOND;
  sprintf(s, "%02.0f:%02.0f:%0*.*f", hh, mm, (ndec > 0) ? ndec+3 : 2, ndec, ss);
  return s;
}

char *dec2str(char *s, const double dec) {
  return dec2strn(s, dec, 3);
}

char *dec2strn(char *s, const double dec, const unsigned int ndec) {
  const char neg = (dec < 0.);
  double ss = (neg) ? -dec : dec;
  const double dd = floor(ss/DEGREE);
  ss -= dd*DEGREE; if(ss < 0) ss = 0;
  const double mm = floor(ss/ARCMIN);
  ss -= mm*ARCMIN; if(ss < 0) ss = 0;
  ss /= ARCSEC;
  sprintf(s, "%+03.0f:%02.0f:%0*.*f", (neg) ? -dd : dd, mm, (ndec > 0) ? ndec+3 : 2, ndec, ss);
  return s;
}
