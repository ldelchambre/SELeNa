/**
 * @file htm_test.c
 * 
 * Hierarchical Triangular Mesh, test file.
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
 
#include <htm.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NSTARS 65536
#define NCONSTRAINTS 1024
#define LEVEL 13

// This file is part of the official release of HTM
#include <../other/cc_aux.c>

// Compare two three-dimensional vector for equality
#define _vec_eq(u,v) (u[0] == v[0] && u[1] == v[1] && u[2] == v[2])

/**
 * Return a random number uniformly drawn between [0, 1]
 */
double randf() {
  return ((double) rand()) / RAND_MAX;
}

/**
 * Structure containing the information we have about a star
 */
typedef struct htm_test_star {
  double v[3]; ///< Cartesian coordinate of the star on the celestial sphere
  htmid id; ///< HTM identifier
  double dmax; ///< Maximal great circle distance between any two vertices of the containing HTM triangle
} htm_test_star_t;

/**
 * Structure containing the informations we have about a constraint
 */
typedef struct htm_test_constraint {
  htm_list_t *l; ///< The list of HTM identifier associated with this constraint
  char *s; ///< The array of stars that were selected using this constraint
} htm_test_constraint_t;

/**
 * Compare two star based on their HTM id
 */
int compar_htm_test_star_t(const void *a, const void *b) {
  const htmid diff = ((htm_test_star_t *) a)->id - ((htm_test_star_t *) b)->id;
  return (diff < 0) ? -1 : ((diff > 0) ? 1 : 0);
}

/**
 * Select the stars that are contained within the provided constraint.
 * s[i] == 1 if star i belong to the constraint; 0 otherwise
 */
void select_stars(char *s, const htm_test_star_t *stars, const htm_list_t *list) {
  unsigned long istar = 0;

  while(list && istar < NSTARS) {
    switch(htm_cmp(stars[istar].id,list->id)) {
      case -1:
        s[istar++] = 0;
        break;
      case 0:
        s[istar++] = 1;
        break;
      case 1:
        list = list->next;
        break;
    }
  }

  while(istar < NSTARS) s[istar++] = 0;
}

/**
 * Build a random constraint
 */
htm_list_t *build_constraint(double *p, double *d, const unsigned int level) {
  // Build the constraint
  htm_cart_coord(p, 2. * M_PI * randf(), asin( 2. * randf() - 1. ));
  *d = M_PI * randf();

  // Retrieve the list of HTM identifiers for this constraint
  return htm_constraint(p, cos(*d), level);
}

/**
 * Check that the constraint effectively match the selected stars
 */
int check_constraint(const char *s, const double *p, const double d, const htm_test_star_t *stars, const unsigned int ilevel, const unsigned int iconstr) {
  double gcirc;
  unsigned long istar;


  for(istar = 0; istar < NSTARS; istar++) {
    // Compute the great circle distance between the star and the centre of the
    // constraint
    gcirc = htm_gcirc_dist(p,stars[istar].v);

    // If the star is contained in the constraint
    if(gcirc <= d) {
      // If the star is not recognized as belonging to the constraint
      if(s[istar] == 0) {
        fprintf(stderr, "%u-%u-%lu) Star contained in constraint but not identified\n", ilevel, iconstr, istar);
        fprintf(stderr, "\tp = (%+8.6e,%+8.6e,%+8.6e), d = %+8.6e, v = (%+8.6e,%+8.6e,%+8.6e), gcirc = %+8.6e, dmax = %+8.6e\n",
                        p[0],p[1],p[2],d,stars[istar].v[0],stars[istar].v[1],stars[istar].v[2], gcirc, stars[istar].dmax);
        return -1;
      }
    // If the star is not contained in a triangle that may overlap the constraint
    } else if(gcirc > d + stars[istar].dmax) {
      // If the star is considered to be embedded in the constraint
      if(s[istar]) {
        fprintf(stderr, "%u-%u-%lu) Star not contained in constraint but identified as such\n", ilevel, iconstr, istar);
        fprintf(stderr, "\tp = (%+8.6e,%+8.6e,%+8.6e), d = %+8.6e, v = (%+8.6e,%+8.6e,%+8.6e), gcirc = %+8.6e, dmax = %+8.6e\n",
                        p[0],p[1],p[2],d,stars[istar].v[0],stars[istar].v[1],stars[istar].v[2], gcirc, stars[istar].dmax);
        return -1;
      }
    }
  }

  return 0;
}

/**
 * Check whether star a and b are both selected
 */
int constraint_and(const char a, const char b) {
  return (a && b);
}

/**
 * Check if star a or star b are selected
 */
int constraint_or(const char a, const char b) {
  return (a || b);
}

/**
 * Check that the combination of two constraints return sensible results
 *
 * cmp function is either @constraint_and or @constraint_or
 */
int check_constraint_combination(const char *s1, const char *s2, const char *scomb, int (*cmp)(const char, const char), const unsigned int ilevel, const unsigned int iconstr) {
  unsigned long istar;

  for(istar = 0; istar < NSTARS; istar++) {
    if(cmp(s1[istar],s2[istar]) != scomb[istar]) {
      fprintf(stderr, "%u-%u-%lu) Wrong combination of constraints\n", ilevel, iconstr, istar);
      fprintf(stderr, "\tcomb(%d,%d) = %d instead of %d\n", s1[istar], s2[istar], scomb[istar], cmp(s1[istar],s2[istar]));
      return -1;
    }
  }

  return 0;
}

/**
 * Build a random star
 */
void build_star(htm_test_star_t *s, const unsigned int level) {
  double ra, dec, v0[3], v1[3], v2[3];
  double d01, d12, d20;

  // Pick up random spherical coordinates uniformly distributed on the sphere
  ra = 2. * M_PI * randf();
  dec = asin( 2. * randf() - 1. );

  // Convert spherical coordinates into cartesian coordinates
  htm_cart_coord(s->v, ra, dec);

  // Compute the HTM id
  s->id = htm_get(s->v, level);

  // Get the vertices of the HTM triangle
  htm_vertices(v0, v1, v2, s->id);

  // Compute the maximal great circle distance between any two vertices
  // of the triangle
  d01 = htm_gcirc_dist(v0, v1);
  d12 = htm_gcirc_dist(v1, v2);
  d20 = htm_gcirc_dist(v2, v0);
  s->dmax = (d01 < d12) ? ((d12 < d20) ? d20 : d12) : ((d01 < d20) ? d20 : d01 );
}

/**
 * Check the HTM functions on a single star
 */
int check_star(const htm_test_star_t *s, const unsigned int level, const unsigned long star) {
  { // Check that the code produce the same result as the official code
    // Get the official HTM id
    long long id = cc_vector2ID(s->v[0],s->v[1],s->v[2],level);

    if(s->id != id) {
      fprintf(stderr, "%u-%lu) HTM id differ from the official id\n", level, star);
      fprintf(stderr, "\tv = (%+8.6e,%+8.6e,%+8.6e), id = %lld instead of %lld",
                      s->v[0], s->v[1], s->v[2], (long long) s->id, id);
      return -1;
    }
  }

  // Check that the level is correctly retrieved from the HTM identifier
  if(level != htm_level(s->id)) {
    fprintf(stderr, "%u-%lu) HTM level is false\n", level, star);
    fprintf(stderr, "\tv = (%+8.6e,%+8.6e,%+8.6e), id = %lld, level = %u\n",
                      s->v[0], s->v[1], s->v[2], (long long) s->id, htm_level(s->id));
  }

  { // Check the conversion from spherical to cartesian coordinates
    double ra, dec, v[3];

    htm_sphere_coord(&ra, &dec, s->v);
    htm_cart_coord(v, ra, dec);

    v[0] -= s->v[0];
    v[1] -= s->v[1];
    v[2] -= s->v[2];

    if(sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]) > 1e-15) {
      fprintf(stderr, "%u-%lu) Spherical coordinates conversion error\n", level, star);
      fprintf(stderr, "\tv = (%+8.6e,%+8.6e,%+8.6e), (ra,dec) = (%+8.6e pi,%+8.6e pi)\n",
                      s->v[0], s->v[1], s->v[2], ra, dec);
      return -1;
    }
  }

  // Check that the conversion between string and HTM id works
  {
    char sid[3+level];
    htmid id;

    htm_string(sid,s->id);

    id = htm_from_string(sid);

    if(s->id != id) {
      fprintf(stderr, "%u-%lu) Conversion between string and HTM id fails\n", level, star);
      fprintf(stderr, "\tv = (%+8.6e,%+8.6e,%+8.6e), id = %lld, s = \"%s\"",
                      s->v[0], s->v[1], s->v[2], (long long) s->id, sid);
      return -1;
    }
  }

  // Check that the HTM vertices are correct
  {
    double n, v[3], v0[3], v1[3], v2[3];

    // Get the vertices
    htm_vertices(v0,v1,v2,s->id);

    // Compute the mid-point
    v[0] = v0[0]+v1[0]+v2[0];
    v[1] = v0[1]+v1[1]+v2[1];
    v[2] = v0[2]+v1[2]+v2[2];
    n = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    v[0] /= n;
    v[1] /= n;
    v[2] /= n;

    // Check that the mid-point is contained in this HTM triangle
    if(s->id != htm_get(v,level)) {
      fprintf(stderr, "%u-%lu) Vertices of the HTM triangle are wrong\n", level, star);
      fprintf(stderr, "\tid = %lld, ", (long long) s->id);
      fprintf(stderr, "v0 = (%+8.6e,%+8.6e,%+8.6e), ", v0[0], v0[1], v0[2]);
      fprintf(stderr, "v1 = (%+8.6e,%+8.6e,%+8.6e), ", v1[0], v1[1], v1[2]);
      fprintf(stderr, "v2 = (%+8.6e,%+8.6e,%+8.6e)\n", v2[0], v2[1], v2[2]);
      return -1;
    }
  }

  return 0;
}

/**
 * Main test function
 */
int main(void) {
  // Allocate stars
  htm_test_star_t *stars = calloc(NSTARS, sizeof(htm_test_star_t));

  // Allocate constraints
  htm_test_constraint_t c = {NULL,NULL}, cprev = {NULL,NULL};
  c.s = calloc(NSTARS,sizeof(char));
  cprev.s = calloc(NSTARS,sizeof(char));
  char *scomb = calloc(NSTARS,sizeof(char));

  // Counters over levels, stars, constraints
  unsigned long istar, iconstr, ilevel;

  // Initialize the random number generator
  srand(1);

  // Browse each HTM level
  for(ilevel = 0; ilevel <= LEVEL; ilevel++) {
    // Build the stars catalog
    for(istar = 0; istar < NSTARS; istar++) {
      // Build the ith star
      build_star(&stars[istar], ilevel);
      // Check the ith star
      check_star(&stars[istar], ilevel, istar);
    }

    // Sort the array of stars according to their HTM id
    qsort(stars, NSTARS, sizeof(htm_test_star_t), compar_htm_test_star_t);

    // Build constraints
    for(iconstr = 0; iconstr < NCONSTRAINTS; iconstr++) {
      double p[3], d;

      // Build this constraint
      c.l = build_constraint(p, &d, ilevel);

      // Select the stars contained in this constraint
      select_stars(c.s, stars, c.l);

      // Check this constraint
      check_constraint(c.s, p, d, stars, ilevel, iconstr);

      // If a previous constraint was already built
      if(iconstr) {
        // Compute the intersection of this constraint with previous one
        // then check it
        htm_list_t *l = htm_list_and(c.l,cprev.l);
        select_stars(scomb, stars, l);
        htm_list_destroy(l);
        check_constraint_combination(c.s,cprev.s,scomb,constraint_and, ilevel, iconstr);

        // Compute the union of this constraint with previous one
        // then check it
        l = htm_list_or(c.l,cprev.l);
        select_stars(scomb, stars, l);
        htm_list_destroy(l);
        check_constraint_combination(c.s,cprev.s,scomb,constraint_or, ilevel, iconstr);

        // We do not need the previous constraint anymore
        htm_list_destroy(cprev.l);
      }

      // Store this constraint as the previous one
      {
        char *tmp = cprev.s;
        cprev.s = c.s;
        c.s = tmp;
        cprev.l = c.l;
      }
    }

    // Destroy the last computed constraint
    htm_list_destroy(cprev.l);
  }

  // Free arrays
  free(stars);
  free(scomb);
  free(c.s);
  free(cprev.s);

  return 0;
}
