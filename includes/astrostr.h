/**
 * @file astrostr.h
 * 
 * Astronomical positions and angle string manipulation, header file.
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
#ifndef _ASTROSTR_H_
#define _ASTROSTR_H_

/**
 * Convert the right ascension contained in s into radian.
 * 
 *   By default, s is considered to be expressed in radian. If followed by 'd' 
 * or 'D', then s is considered to be expressed in degree (e.g. "30d"). If s
 * matches the format hh:mm:ss.ddd, then sexagesimal coordinates are assumed 
 * were hh is epxressed in hour, mm in minutes, ss in seconds and ddd in 
 * milli-seconds (e.g. "9:14:36.7845").
 * 
 *  If send is not NULL, then it will point to the character after the last 
 * character that was converted (e.g. if s = "00:59:59toto", then on return 
 * *send = "toto"). Note that we will always have that s < *send if str2ra 
 * return 0.
 * 
 * @param ra   The right ascension associated with s in radian (0 <= ra < 2*pi).
 * @param s    The string to be converted.
 * @param send If not NULL, point to the character after the last character used in the conversion.
 * 
 * @return 0 on sucess, or -1 if s cannot be converted (in which case ra and send are undefined).
 */
int str2ra(double *ra, const char *s, char **send);

/**
 * Convert the declination contained in s into radian.
 * 
 *   By default, s is considered to be expressed in radian. If followed by 'd' 
 * or 'D', then s is considered to be expressed in degree (e.g. "30d"). If s
 * matches the format dd:mm:ss.ddd, then sexagesimal coordinates are assumed 
 * were dd is epxressed in degree, mm in arcminutes, ss in arcseconds and ddd in 
 * milli-arcseconds (e.g. "-12:59:00.4").
 * 
 *  If send is not NULL, then it will point to the character after the last 
 * character that was converted (e.g. if s = "00:59:59toto", then on return 
 * *send = "toto"). Note that we will always have that s < *send if str2dec
 * return 0.
 * 
 * @param ra   The declination associated with s in radian (0 <= ra < 2*pi).
 * @param s    The string to be converted.
 * @param send If not NULL, point to the character after the last character used in the conversion.
 * 
 * @return 0 on sucess, or -1 if s cannot be converted (in which case ra and send are undefined).
 */
int str2dec(double *dec, const char *s, char **send);

/**
 * Convert angle contained in s into radian.
 * 
 * By default, s is considered to be expressed in radian. 
 * The following units can be happened after digits:
 *   degreed    : The given angle is expressed in degree (e.g. "30d" or "15degree")
 *   arcmin / ' : The given angle is expressed in arcminute (e.g. "1arcmin" or "2'")
 *   arcsec / " : The given angle is expressed in arcsecond (e.g. "3arcsec" or "4\"")
 *   mas        : The given angle is expressed in milli-arcsecond (e.g. "5mas")
 *   hour / h   : The given angle is expressed in hour (e.g. "6hour" or "7h")
 *   minute / m : The given angle is expressed in minute (e.g. "8minute" or "9m")
 *   second / s : The given angle is expressed in second (e.g. "1second" or "2s").
 * 
 *  If send is not NULL, then it will point to the character after the last 
 * character that was converted (e.g. if s = "00:59:59toto", then on return 
 * *send = "toto"). Note that we will always have that s < *send if str2dec
 * return 0.
 * 
 * @param theta The angle associated with s in radian.
 * @param s     The string to be converted.
 * @param send  If not NULL, point to the character after the last character used in the conversion.
 * 
 * @return 0 on sucess, or -1 if s cannot be converted (in which case ra and send are undefined).
 */
int str2angle(double *theta, const char *s, char **send);

/**
 * Convert a right ascension into a sexagesimal string with 3 digits of precision.
 * 
 * On return, s will contains the right ascension in the form "hh:mm:ss.sss"
 * along with a ending null character. s should hence have a minimal size of 13
 * characters.
 * 
 * @param s  The sexagesimal representation of ra (must have a size of at least 13 characters).
 * @param ra The right ascension to be converted. ra must be in the range [0, 2*pi[.
 * 
 * @return A pointer to s
 */
char *ra2str(char *s, const double ra);

/**
 * Convert a right ascension into a sexagesimal string with ndec digits of precision.
 * 
 * On return, s will contains the right ascension in the form "hh:mm:ss" where
 * ss has a precision of ndec digits (e.g. if ndec = 2 then we have 
 * "hh:mm:ss.ss", if ndec = 0 then we have "hh:mm:ss") along with a ending null 
 * character. s should hence have a minimal size of 
 *    9      , if ndec = 0
 *    10+ndec, if ndec > 0
 * characters.
 * 
 * @param s    The sexagesimal representation of ra (see remarks on sizes in the description).
 * @param ra   The right ascension to be converted. ra must be in the range [0, 2*pi[.
 * @param ndec The number of digits to use in the representation of second.
 * 
 * @return A pointer to s
 */
char *ra2strn(char *s, const double ra, const unsigned int ndec);

/**
 * Convert a declination into a sexagesimal string with 3 digits of precision.
 * 
 * On return, s will contains the declination in the form "+dd:mm:ss.sss"
 * along with a ending null character. s should hence have a minimal size of 14
 * characters. Note that the "-" or "+" sign will always be present in this 
 * representation.
 * 
 * @param s   The sexagesimal representation of dec (must have a size of at least 14 characters).
 * @param dec The declination to be converted. dec must be in the range [-pi/2, pi/2].
 * 
 * @return A pointer to s
 */
char *dec2str(char *s, const double dec);

/**
 * Convert a declination into a sexagesimal string with ndec digits of precision.
 * 
 * On return, s will contains the declination in the form "dd:mm:ss" where
 * ss has a precision of ndec digits (e.g. if ndec = 2 then we have 
 * "+dd:mm:ss.ss", if ndec = 0 then we have "dd:mm:ss") along with a ending null 
 * character. s should hence have a minimal size of 
 *    10     , if ndec = 0
 *    11+ndec, if ndec > 0
 * characters.
 * 
 * @param s    The sexagesimal representation of dec (see remarks on sizes in the description).
 * @param dec  The declination to be converted. dec must be in the range [-pi/2, pi/2].
 * @param ndec The number of digits to use in the representation of tha arcsecond.
 * 
 * @return A pointer to s
 */
char *dec2strn(char *s, const double dec, const unsigned int ndec);

#endif
