/**
 * @file htmdb_test.c
 * 
 * Hierarchical Triangular Mesh database, test file.
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
#include <htmdb.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NSTARS (1ULL << 20)
#define NCONSTRAINTS 128ULL
#define LEVEL 8
#define CONSTR_MAX_D M_PI // The maximal search radius to use (1 degree)

typedef struct _htmdb_test_constraint {
  double p[3];
  double cosd;
} _htmdb_test_constraint_t;

int _htmdb_test_process(const double *v, void *data, void *arg) {
  _htmdb_test_constraint_t *constr = (_htmdb_test_constraint_t *) data;
  unsigned long long *nfound = (unsigned long long *) arg;

  if(v[0]*constr->p[0]+v[1]*constr->p[1]+v[2]*constr->p[2] < constr->cosd) {
    fprintf(stderr, "Star v=(%+8.6e,%+8.6e,%+8.6e) is not contained in constraint ", v[0],v[1],v[2]);
    fprintf(stderr, "p=(%+8.6e,%+8.6e,%+8.6e), cosd=%+8.6e\n", constr->p[0], constr->p[1], constr->p[2], constr->cosd);
    return -1;
  }

  (*nfound)++;

  return 0;
}

int main(void) {
  // Create the HTM database
  htmdb_t *db = htmdb_init(LEVEL);

  // Current constraint
  _htmdb_test_constraint_t constr;

  // Number of stars belonging to this constraint
  unsigned long long nfound;

  // Star positions
  double ra, dec;
  double v[3];

  // General usage counter(s)
  unsigned long long i;

  // If the HTM database creation fails
  if(!db) {
    fprintf(stderr, "Unable to create an HTM database of level %d\n", LEVEL);
    return -1;
  }

  // Initialize the random number generator
  srand(1);

  // Add NSTARS random stars uniformly drawn on the sphere
  for(i = 0; i < NSTARS; i++) {
    // Uniformly draw a point on the sphere (in spherical coordinates)
    ra = 2.*M_PI*rand()/RAND_MAX;
    dec = asin(2.*rand()/RAND_MAX-1);

    // Convert spherical coordinates into cartesian coordinates
    htm_cart_coord(v, ra, dec);

    // Add this entry to the HTM database
    if(htmdb_add(db, v, &constr)) {
      fprintf(stderr, "%llu) Unable to add entry to the HTM database\n", i);
      htmdb_destroy(db, NULL);
      return -1;
    }
  }

  // Probe NCONSTRAINTS random constraints
  for(i = 0; i < NCONSTRAINTS; i++) {
    // Pick up a random direction for this constraint
    ra = 2.*M_PI*rand()/RAND_MAX;
    dec = asin(2.*rand()/RAND_MAX-1);
    htm_cart_coord(constr.p, ra, dec);

    // Pick up a search radius within [0, CONSTR_MAX_D]
    constr.cosd = cos(CONSTR_MAX_D*rand()/RAND_MAX);

    // Set the number of stars found to zero
    nfound = 0;

    // Process this constraint
    if(htmdb_constraint(db, constr.p, constr.cosd, _htmdb_test_process, &nfound)) {
      fprintf(stderr, "%llu) An error occured during the processing of constraint ", i);
      fprintf(stderr, "p=(%+8.6e,%+8.6e,%+8.6e), cosd=%+8.6e\n", constr.p[0], constr.p[1], constr.p[2], constr.cosd);
      htmdb_destroy(db, NULL);
      return -1;
    }

    // Process the complementary constraint
    constr.cosd *= -1; constr.p[0] *= -1; constr.p[1] *= -1; constr.p[2] *= -1;
    if(htmdb_constraint(db, constr.p, constr.cosd, _htmdb_test_process, &nfound)) {
      fprintf(stderr, "%llu) An error occured during the processing of the complementary constraint ", i);
      fprintf(stderr, "p=(%+8.6e,%+8.6e,%+8.6e), cosd=%+8.6e\n", constr.p[0], constr.p[1], constr.p[2], constr.cosd);
      htmdb_destroy(db, NULL);
      return -1;
    }

    // Check that the constraint and its complementary contain all stars
    if(nfound != NSTARS) {
      fprintf(stderr, "%llu) Not all stars where identified during the processing of the constraint and of its complementary ", i);
      fprintf(stderr, "(%llu stars identified instead of %llu)\n", nfound, NSTARS);
      htmdb_destroy(db, NULL);
      return -1;
    }
  }

  // Free the HTM database
  htmdb_destroy(db, NULL);

  return 0;
}
