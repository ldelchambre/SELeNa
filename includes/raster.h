/**
 * @file raster.h
 * 
 * Triangle rasterization, header file.
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
#ifndef _RASTER_H_
#define _RASTER_H_

/**
 * Compute the indices of the grid defined by (xfrom,xto,nx,yfrom,yto,ny) 
 * for which the triangle (v0, v1, v2) overlap the 'pixels' of the grid, this is
 * rasterization.
 * 
 * Each set of indices is returned through a callback function 
 *   int process(const unsigned int ix, const unsigned int iy, void *arg)
 * where the square defined by
 *   x0 = xfrom + ( xto - xfrom ) / nx * ix
 *   x1 = xfrom + ( xto - xfrom ) / nx * ( ix + 1 )
 *   y0 = yfrom + ( yto - yfrom ) / ny * iy
 *   y1 = yfrom + ( yto - yfrom ) / ny * ( iy + 1 )
 * overlaps the triangle (v0, v1, v2).
 * 
 *  If the callback function, return 0, then all indices will be browsed whereas
 * if the returned value differs from 0, then the browsing stop and the function
 * return -1UL.
 * 
 * @param v0    First vertex of the triangle
 * @param v1    Second vertex of the triangle
 * @param v2    Third vertex of the triangle
 * @param xfrom Starting x position of the grid
 * @param xto   Ending x position of the grid
 * @param nx    Number of 'pixels' in the x direction
 * @param yfrom Starting y position of the grid
 * @param yto   Ending y position of the grid
 * @param ny    Number of 'pixels' in the y direction
 * @param process Callback function
 * @param arg   Parameters to pass to the callback function
 * 
 * @return The number of pairs of indices we found, or -1UL if process returned something different thant 0
 */
unsigned int raster(const double v0[2], const double v1[2], const double v2[2],
                    const double xfrom, const double xto, const unsigned int nx,
                    const double yfrom, const double yto, const unsigned int ny,
                    int (*process)(const unsigned int, const unsigned int, void *),
                    void *arg);

#endif
