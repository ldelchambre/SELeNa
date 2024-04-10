/**
 * @file raster_test.c
 * 
 * Triangle rasterization, test file.
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
#include <math.h>
#include <raster.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Perform 2^20 tests
#define NTEST (1 << 20)
#define NMAX 25

#define raster_test_edge(p0,p1,x,y) ( ( y - p0[1] ) * ( p1[0] - p0[0] ) - ( x - p0[0] ) * ( p1[1] - p0[1] ) )

typedef struct raster_test_data {
  double xfrom, xto, dx, yfrom, yto, dy;
  unsigned int nx, ny;
  double v0[2], v1[2], v2[2];
  char *map;
} raster_test_data_t;

int raster_test_is_contained(const double x, const double y, const raster_test_data_t *data) {
  const double e01 = raster_test_edge(data->v0,data->v1,x,y);
  const double e12 = raster_test_edge(data->v1,data->v2,x,y);
  const double e20 = raster_test_edge(data->v2,data->v0,x,y);
  return (signbit(e01) == signbit(e12) && signbit(e12) == signbit(e20));
}

int raster_fill_map(const unsigned int ix, const unsigned int iy, void *arg) {
  raster_test_data_t *data = (raster_test_data_t *) arg;
  char (*map)[data->ny] = (char (*)[data->ny]) data->map;
  map[ix][iy] = 1;
  return 0;
}

int main(void) {
  raster_test_data_t data;
  unsigned int i, ix, iy;
  double x, y;

  // Allocate the pixels map
  data.map = (char *) calloc(NMAX*NMAX, sizeof(char));
  if(!data.map) {
    fprintf(stderr, "Unable to allocate the pixel map\n");
    return -1;
  }

  // Perform NTEST tests
  for(i = 0; i < NTEST; i++) {
    // Initialize the random number generator
    srand(i+1);

    // Generate a random grid
    data.xfrom = 2.*rand()/RAND_MAX-1.;
    data.xto   = 2.*rand()/RAND_MAX-1.;
    data.nx    = rand() % NMAX + 1;
    data.yfrom = 2.*rand()/RAND_MAX-1.;
    data.yto   = 2.*rand()/RAND_MAX-1.;
    data.ny    = rand() % NMAX + 1;

    // Ensure that xfrom <= xto and that yfrom <= yto
    if(data.xfrom > data.xto) {
      const double tmp = data.xfrom;
      data.xfrom = data.xto;
      data.xto = tmp;
    }
    if(data.yfrom > data.yto) {
      const double tmp = data.yfrom;
      data.yfrom = data.yto;
      data.yto = tmp;
    }

    // Compute the pixel sampling
    data.dx = (data.xto - data.xfrom) / data.nx;
    data.dy = (data.yto - data.yfrom) / data.ny;

    // Generate a random triangle in [-1, 1]
    data.v0[0] = 2.*rand()/RAND_MAX-1.;
    data.v0[1] = 2.*rand()/RAND_MAX-1.;
    data.v1[0] = 2.*rand()/RAND_MAX-1.;
    data.v1[1] = 2.*rand()/RAND_MAX-1.;
    data.v2[0] = 2.*rand()/RAND_MAX-1.;
    data.v2[1] = 2.*rand()/RAND_MAX-1.;

    // Reset the pixel map
    memset(data.map, 0, data.ny*data.nx);

    // Compute the rasterization of (v0,v1,v2)
    raster(data.v0, data.v1, data.v2, data.xfrom, data.xto, data.nx, data.yfrom, data.yto, data.ny, raster_fill_map, &data);

    // Browse each pixel of the map
    {
      char (*map)[data.ny] = (char (*)[data.ny]) data.map; // Re-interpret map so as to ease understanding
      for(iy = 0; iy < data.ny; iy++) {
        for(ix = 0; ix < data.nx; ix++) {
          // Pick-up a random point falling in pixel (ix,iy)
          // i.e. xfrom + ( xto - xfrom ) / nx * ix <= x < xfrom + ( xto - xfrom ) / nx * (ix+1)
          //      yfrom + ( yto - yfrom ) / ny * ix <= x < yfrom + ( yto - yfrom ) / ny * (iy+1)
          x = data.xfrom + ((double) ix + (double) rand() / (1.+RAND_MAX)) * data.dx;
          y = data.yfrom + ((double) iy + (double) rand() / (1.+RAND_MAX)) * data.dy;

          // If map[ix][iy] == 0, then (x,y) should not be contained in (v0,v1,v2)
          if(map[ix][iy] == 0 && raster_test_is_contained(x, y, &data)) {
            fprintf(stderr, "%u) Point (%g,%g) is contained in triangle (%g,%g), (%g,%g), (%g,%g) but map[%u][%u]=0 with xfrom=%g, xto=%g, nx=%u, yfrom=%g, yto=%g, ny=%u\n",
                            i, x, y, data.v0[0], data.v0[1], data.v1[0], data.v1[1], data.v2[0], data.v2[1], ix, iy, data.xfrom, data.xto, data.nx, data.yfrom, data.yto, data.ny);
            fprintf(stderr, "  %g %g %g %g %g %g %g %g\n", x, y, data.v0[0], data.v0[1], data.v1[0], data.v1[1], data.v2[0], data.v2[1]);
          }
        }
      }
    }
  }

  // Delete the rasterization map
  free(data.map);

  return 0;
}
