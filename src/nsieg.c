/**
 * @file nsieg.c
 * 
 * Non-singular isothermal ellipsoid in presence of an external shear lens model, source file.
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
#include <nsieg.h>
#include <math.h>

void nsieg_alpha(double alpha[2],
                 const double theta[2],
                 const double b,
                 const double e,
                 const double theta_e,
                 const double gamma,
                 const double theta_gamma,
                 const double s) {
  // Compute some intermediate variables
  const double se   = sin(theta_e);
  const double ce   = cos(theta_e);
  const double x    =  ce * theta[0] + se * theta[1]; // Align the ellipse such that the
  const double y    = -se * theta[0] + ce * theta[1]; // Semi-major axis stands along the x axis
  const double q    = 1. - e;             // Compute the axis ratio from ellipticity
  const double qp   = sqrt( 1. - q * q ); // Compute the eccentricity from the axis ratio
  const double psi  = sqrt( q * q * ( s * s + x * x ) + y * y );
  const double s2g  = sin(2.*theta_gamma);
  const double c2g  = cos(2.*theta_gamma);
  double alphax, alphay;

  // Compute the deflection angle with the semi-major axis of the ellipse
  // aligned to the x axis
  if(e == 0.) { // In the case of the NSISg lens model, another formulation avoid
                // the division by zero because qp = 0
    alphax = b * x / (s + psi);
    alphay = b * y / (s + psi);
  } else {
    alphax = b * q * atan( qp * x / ( psi + s )) / qp;
    alphay = b * q * atanh( qp * y / ( psi + q * q * s)) / qp;
  }
  
  // Rotate the deflection angle so as to take theta_e into account
  // and add the contribution of the external shear
  alpha[0] = ce * alphax - se * alphay + gamma * ( c2g * theta[0] + s2g * theta[1] );
  alpha[1] = se * alphax + ce * alphay + gamma * ( s2g * theta[0] - c2g * theta[1] );
}

double nsieg_mu(const double theta[2],
                const double b,
                const double e,
                const double theta_e,
                const double gamma,
                const double theta_gamma,
                const double s) {
  // Compute some intermediate variables
  const double se   = sin(theta_e);
  const double ce   = cos(theta_e);
  const double x    =  ce * theta[0] + se * theta[1]; // Align the ellipse such that the
  const double y    = -se * theta[0] + ce * theta[1]; // Semi-major axis stands along the x axis
  const double q    = 1. - e;             // Compute the axis ratio from ellipticity
  const double psi  = sqrt( q * q * ( s * s + x * x ) + y * y );
  const double s2g  = sin(2.*theta_gamma);
  const double c2g  = cos(2.*theta_gamma);
  const double q2   = q * q;
  const double qp2  = 1. - q * q;
  
  // Compute the partial derivatives of psi
  const double dpsidt0 = (q2 * ce * x - se * y) / psi;
  const double dpsidt1 = (q2 * se * x + ce * y) / psi;
  
  // Compute the partial derivatives of alphax
  const double psis       = psi + s;
  const double dalphax    = b * q / ( psis * psis + qp2 * x * x );
  const double dalphaxdt0 = dalphax * ( psis * ce - x * dpsidt0 );
  const double dalphaxdt1 = dalphax * ( psis * se - x * dpsidt1 );
  
  // Compute the partial derivative of alphay
  const double psip       = psi + q2 *s;
  const double dalphay    = b * q / ( psip * psip - qp2 * y * y );
  const double dalphaydt0 = dalphay * ( psip * (-se) - y * dpsidt0 );
  const double dalphaydt1 = dalphay * ( psip *   ce  - y * dpsidt1 );
  
  // Compute the Jacobian matrix of the lens equation theta_s = theta - alpha
  const double J00 = 1. - ( ce * dalphaxdt0 - se * dalphaydt0 + gamma * c2g );
  const double J01 = ce * dalphaxdt1 - se * dalphaydt1 + gamma * s2g;
  /* J10 = J01 = se * dalphaxdt0 + ce * dalphaydt0 + gamma * s2g; */
  const double J11 = 1. - ( se * dalphaxdt1 + ce * dalphaydt1 - gamma * c2g );
  
  // The amplification is the determinant of the inverse Jacobian matrix
  return 1./(J00*J11-J01*J01);
}
