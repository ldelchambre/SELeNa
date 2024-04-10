/**
 * @file raster.c
 * 
 * Triangle rasterization, source file.
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
#include <float.h>
#include <math.h>
#include <raster.h>

/**
 *   Compute the range of x values to explore regarding the intersection of 
 * the line of constant value y = y0 with the triangle given by (v0, v1, v2).
 * 
 * @param x   The range of x values to explore (x[0] <= x[1] but see below)
 * @param v0  The coordinates of the vertex having the lowest y
 * @param v1  The coordinates of the vertex having y inbetween v0 and v2
 * @param v2  The coordinates of the vertex having the highest y
 * @param l02 The equation of the line joining v0 to v2 (y = l02[0] * x + l02[1]) 
 * @param l01 The equation of the line joining v0 to v1 (y = l01[0] * x + l01[1]) 
 * @param l12 The equation of the line joining v1 to v2 (y = l12[0] * x + l12[1])
 * @param y0  The y value against which we have compute the intersection
 * 
 * Note that if no intersection is found, then x[0] > x[1]
 * Note that we must have v0[1] <= v1[1] <= v2[1]
 */
void _raster_xintersect(double x[2],
                        const double v0[2], const double v1[2], const double v2[2],
                        const double l02[2], const double l01[2], const double l12[2],
                        const double y0) {
  // Initialize the xrange to return
  x[0] = +DBL_MAX;
  x[1] = -DBL_MAX;
  
  // If the triangle intersect the line of constant value y = y0, compute these
  // intersections
  if(v0[1] < y0 && y0 < v2[1]) {
    // Find the x value where the line joining v0 to v2 intercept y
    x[0] = (v0[0] == v2[0]) ? v0[0] : (y0 - l02[1]) / l02[0];

    // If the line joining v0 to v1 intercept the line y = y0
    if(y0 < v1[1])
      x[1] = (v0[0] == v1[0]) ? v0[0] : (y0 - l01[1]) / l01[0];
    else // If the line joining v1 to v2 intercept the line y = y0
      x[1] = (v1[0] == v2[0]) ? v1[0] : (y0 - l12[1]) / l12[0];

    // Sort the intersections according to their x values
    if(x[1] < x[0]) {
      const double tmp = x[0];
      x[0] = x[1];
      x[1] = tmp;
    }
  }
}

unsigned int raster(const double v0[2], const double v1[2], const double v2[2],
                    const double xfrom, const double xto, const unsigned int nx,
                    const double yfrom, const double yto, const unsigned int ny,
                    int (*process)(const unsigned int, const unsigned int, void *),
                    void *arg) {
  // Sort vertices according to their y positions
  // v0s[1] <= v1s[1] <= v2s[1]
  const double *v0s, *v1s, *v2s;
  if(v0[1] < v1[1]) { // v0 < v1
    if(v1[1] < v2[1]) { // v0 < v1 < v2
      v0s = v0; v1s = v1; v2s = v2;
    } else { // v0 < v1 && v2 <= v1
      if(v0[1] < v2[1]) { // v0 < v2 <= v1
        v0s = v0; v1s = v2; v2s = v1;
      } else { // v2 <= v0 < v1
        v0s = v2; v1s = v0; v2s = v1;
      }
    }
  } else { // v1 <= v0
    if(v2[1] < v1[1]) { // v2 < v1 <= v0
      v0s = v2; v1s = v1; v2s = v0;
    } else { // v1 <= v0 && v1 <= v2
      if(v0[1] < v2[1]) { // v1 <= v0 < v2
        v0s = v1; v1s = v0; v2s = v2;
      } else { // v1 <= v2 <= v0
        v0s = v1; v1s = v2; v2s = v0;
      }
    }
  }

  // If the maximal y value of the vertices is below yfrom or the minimal y value
  // is above yto, then nothing has to be done.
  if(yto <= v0s[1] || v2s[1] < yfrom)
    return 0;

  // Compute the size of the pixels in x and y
  const double dx = ( xto - xfrom ) / nx;
  const double dy = ( yto - yfrom ) / ny;

  // Compute the vertices positions in pixel space
  const double v0p[2] = { (v0s[0] - xfrom) / dx, (v0s[1] - yfrom) / dy };
  const double v1p[2] = { (v1s[0] - xfrom) / dx, (v1s[1] - yfrom) / dy };
  const double v2p[2] = { (v2s[0] - xfrom) / dx, (v2s[1] - yfrom) / dy };
  
  // Compute the y pixel indices of all vertices
  const unsigned int iy0 = (v0p[1] < 0.) ? -1U : v0p[1];
  const unsigned int iy1 = (v1p[1] < 0.) ? -1U : v1p[1];
  const unsigned int iy2 = (v2p[1] < 0.) ? -1U : v2p[1];
  
  // Compute the y indices of the pixels to browse
  unsigned int       iy     = (iy0 == -1U) ? 0   : iy0; // max(v0p[1], 0)
  const unsigned int iyto   = (iy2 >= ny)  ? ny-1: iy2; // min(v2p[1], n-1)

  // Compute the lines equation in pixel space. Format: y = lij[0] * x + lij[1]
  // where lij is the line joining the point vip and the point vjp.
  // Beware that division by zero can happen here if vip[0] == vjp[0].
  double l02[2] = { ( v2p[1] - v0p[1] ) / ( v2p[0] - v0p[0] ), v0p[1] - l02[0] * v0p[0] };
  double l01[2] = { ( v1p[1] - v0p[1] ) / ( v1p[0] - v0p[0] ), v0p[1] - l01[0] * v0p[0] };
  double l12[2] = { ( v2p[1] - v1p[1] ) / ( v2p[0] - v1p[0] ), v1p[1] - l12[0] * v1p[0] };

  // The range of x values to explore from line y = iy (xlow) to line y = iy+1 (xhigh)
  // xrange = {min(xlow[0],xhigh[0]), max(xlow[1],xhigh[1])};
  double xlow[2], xhigh[2], xrange[2];
  
  // The x indices of the pixels to browse for line iy
  unsigned int ix, ixto;
  
  // The number of pixel we found
  unsigned int npix = 0;

  // Compute the intersection of the triangle with the axis y = iy
  _raster_xintersect(xlow, v0p, v1p, v2p, l02, l01, l12, iy);

  // Browse the y indices
  for(; iy <= iyto; iy++) {
    // Update xlow according to the vertices that may be contained between iy and iy+1
    if(iy == iy0) {
      if(v0p[0] < xlow[0]) xlow[0] = v0p[0];
      if(v0p[0] > xlow[1]) xlow[1] = v0p[0];
    }
    if(iy == iy1) {
      if(v1p[0] < xlow[0]) xlow[0] = v1p[0];
      if(v1p[0] > xlow[1]) xlow[1] = v1p[0];
    }
    if(iy == iy2) {
      if(v2p[0] < xlow[0]) xlow[0] = v2p[0];
      if(v2p[0] > xlow[1]) xlow[1] = v2p[0];
    }
    
    // Compute the intersection of the triangle with the axis y = iy + 1
    _raster_xintersect(xhigh, v0p, v1p, v2p, l02, l01, l12, 1. + iy);
    
    // Compute the x range to browse
    xrange[0] = (xlow[0] < xhigh[0]) ? xlow[0] : xhigh[0];
    xrange[1] = (xlow[1] > xhigh[1]) ? xlow[1] : xhigh[1];
    
    // If the xrange overlap [0, nx[
    if(0 <= xrange[1] && xrange[0] < nx) {
      // Compute the x indices to browse
      ix   = (xrange[0] < 0.)  ? 0.   : xrange[0];
      ixto = (xrange[1] >= nx) ? nx-1 : xrange[1];
      // Browse the ix and iy indices
      for(; ix <= ixto; ix++) {
        if(process(ix, iy, arg))
          return -1;
        npix++;
      }
    }

    // Copy the higher intersections as the upcoming lower intersections
    xlow[0] = xhigh[0];
    xlow[1] = xhigh[1];
  }
  
  return npix;
}
