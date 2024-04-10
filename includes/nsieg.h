/**
 * @file nsieg.h
 * 
 * Non-singular isothermal ellipsoid in presence of an external shear lens model, header file.
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
#ifndef _NSIEG_H_
#define _NSIEG_H_

/**
 * Compute the deflection angle that is associated with a Non-Singular
 * Isothermal Ellipsoid lens model in presence of an external shear (NSIEg
 * lens model).
 *
 *   NSIEg lens models are characterized by a projected mas density of
 *
 *             kappa(x,y) = 0.5 * b / sqrt(s^2 + x^2 + y^2 / q^2),
 *
 * where q = 1 - e = a / b is the projected axis ratio and where a and b are
 * respectively the length of the semi-minor axis and of the semi-major axis of
 * the ellipse of constant mass density, or equivalently of constant isophot
 * contours. The mass distribution hence stands at the origin along with its
 * semi-major axis aligned with the x-axis.
 *
 * Given two dimensional cartesian coorinates taken in the image plane, theta,
 * this function retrieve the two dimensional deflection angle that is associated
 * with the provided model and that is such that the source position theta_s can
 * be computed trough
 *   theta_s = theta - alpha.
 *
 * @param alpha       The two dimensional deflection angle that is associated with theta
 * @param theta       The two dimensional position in the image plane for which alpha must be computed
 * @param b           Normalization factor. In the limit of a singular (s = 0) and spherical (e = 0) \
 *                    lens model, b stand to be the Einstein radius of the model (Keeton, 2001, p. 8, eq. 36).
 * @param e           The ellipticity of the mass distribution, or equivalently of the isophot contours.
 * @param theta_e     The inclination of the semi-major axis of the mass distribution measured East of north [rad].
 * @param gamma       The external shear strength.
 * @param theta_gamma The direction of the external shear, measured East of North [rad].
 * @param s           Core radius
 *
 * @see Keeton, C. R., "Software for Gravitational Lensing", http://www.physics.rutgers.edu/~keeton/gravlens/manual.pdf
 * @see Kormann, R. et al, "Isothermal elliptical gravitational lens models", 1994, A&A, 284(1), 285
 * @see Keeton, C. R. and Kochanek, C. S., "Gravitational lensing by spiral galaxies", 1998, ApJ, 495, 157
 * @see Keeton, C. R., "A Catalog of Mass Models for Gravitational Lensing", 2001, arXiv:astro-ph/0102341
 */
void nsieg_alpha(double alpha[2],
                 const double theta[2],
                 const double b,
                 const double e,
                 const double theta_e,
                 const double gamma,
                 const double theta_gamma,
                 const double s);

/**
 * Compute the signed amplification that is associated with a Non-Singular
 * Isothermal Ellipsoid lens model in presence of an external shear (NSIEg
 * lens model).
 * 
 *  The returned amplification is signed in the sense that it account for the 
 * swapping of the images when the source position crosses the caustics. Beware
 * that a source standing right on a caustic curve will produce an infinite 
 * amplification of the source.
 *
 * @param theta       The two dimensional position in the image plane for which the amplification must be computed
 * @param b           Normalization factor.
 * @param e           The ellipticity of the mass distribution, or equivalently of the isophot contours.
 * @param theta_e     The inclination of the semi-major axis of the mass distribution measured East of north [rad].
 * @param gamma       The external shear strength.
 * @param theta_gamma The direction of the external shear, measured East of North [rad].
 * @param s           Core radius
 *
 * @return The signed amplification in the image plane theta
 */
double nsieg_mu(const double theta[2],
                const double b,
                const double e,
                const double theta_e,
                const double gamma,
                const double theta_gamma,
                const double s);

#endif
