/**
 * @file htm.h
 * 
 * Hierarchical Triangular Mesh, header file.
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
#ifndef _HTM_H_
#define _HTM_H_

#include <limits.h>

/**
 * Signed integer designed to contain HTM identifiers.
 *
 * @note Negative identifiers are used so as to represent errors
 */
typedef long long int htmid;

/**
 * The number of bits used in order to represent HTM identifiers
 */
#define HTMID_BITS ((int) sizeof(htmid)*CHAR_BIT)

/**
 * The maximal HTML level
 *
 * @note We set a limit of 25 because at this level, numerical errors start to appear
 */
#define HTM_LEVEL_MAX (((HTMID_BITS-6)/2 <= 25) ? (HTMID_BITS-6)/2 : 25)

/**
 * The minimal value of a valid HTM identifier
 */
#define HTMID_MIN ((htmid) 0x08LL)

/**
 * The maximal value of a valid HTM identifier
 */
#define HTMID_MAX ((htmid) (0x10LL << (2LL * HTML_LEVEL_MAX)) - 1LL)

/**
 * Value representing a bad HTM identifier
 */
#define HTMID_BAD ((htmid) -1LL)

/**
 * A list of HTM identifiers.
 *
 * @note Lists of HTM identifiers must be freed with htm_list_destroy
 */
typedef struct htm_list {
  htmid id; ///< The HTM identifier
  struct htm_list *next; ///< The next HTM identifier or NULL
} htm_list_t;

/**
 * Convert spherical coordinates into cartesian coordinates
 *
 * @param v Output vector of three dimensional cartesian coordinates
 * @param ra Right ascension [rad]
 * @param dec Declination [rad]
 */
void htm_cart_coord(double *v, const double ra, const double dec);

/**
 * Compare two HTM identifier
 *
 * @param id1 The first HTM identifier (in [HTMID_MIN, HTMID_MAX])
 * @param id2 The second HTM identifier (in [HTMID_MIN, HTMID_MAX])
 *
 * @return -1 if (id1 < id2); 1 if (id2 > id1); 0 if id1 and id2 overlap
 *
 * @note HTM identifier are compared based on the lowest level contained in id1 or id2
 */
int htm_cmp(const htmid id1, const htmid id2);

/**
 * Build a HTM constraint
 *
 * @param p The three dimensional cartesian coordinates of the direction of the constraint (must have a norm of 1)
 * @param cosd The cosine of the search radius of the constraints (in [-1, 1])
 * @param level The maximal HTM recursion level to use (in [0, HTML_LEVEL_MAX])
 *
 * @return The list of HTM identifiers overlapping this constraint ordered in increasing order
 */
htm_list_t *htm_constraint(const double *p, const double cosd, const unsigned int level);

/**
 * Convert a string representation of an HTM identifier into an integer
 *
 * @param s The string to convert
 *
 * @return The HTM identifier or HTMID_BAD if conversion fails
 */
htmid htm_from_string(const char *s);

/**
 * Copy the three dimensional vertex v into u
 *
 * @param u Three dimensional coordinates of the destination
 * @param u Three dimensional coordinates of the source
 */
#define htm_copy_vertex(u,v) do {u[0]=v[0];u[1]=v[1];u[2]=v[2];} while(0)

/**
 * Compute the great circle distance between two point on the sphere if
 * these are expressed in cartesian coordinates
 *
 * @param u Three dimensional coordinates of the first point (must have norm of 1)
 * @param v Three dimensional coordinates of the second point (must have norm of 1)
 *
 * @return The great circle distance between u and v
 */
#define htm_gcirc_dist(u,v) acos(u[0]*v[0]+u[1]*v[1]+u[2]*v[2])

/**
 * Get the HTM identifier that is associated with a point on the sphere for
 * a given level.
 *
 * @param v Three dimensional cartesian coordinates of the point (must have norm of 1)
 * @param level The HTM recursion level to use (in [0, HTM_LEVEL_MAX])
 *
 * @return The HTM identifier associated with v or HTMID_BAD if level > HTM_LEVEL_MAX
 */
htmid htm_get(const double *v, const unsigned level);

/**
 * Get the HTM level that is associated with an identifier
 *
 * @param id The HTM identifier (in [HTMID_MIN, HTMID_MAX])
 *
 * @return The HTM level that is associated with id
 *
 * @note Result is undefined if id is not in [HTMID_MIN, HTMID_MAX]
 */
unsigned int htm_level(const htmid id);

/**
 * Compute the intersection between two lists of HTM identifiers
 *
 * @param l1 First list of identifiers (must be sorted in ascending order)
 * @param l2 Second list of identifiers (must be sorted in ascending order)
 *
 * @return A newly allocated list containing the intersection between l1 and l2, sorted in ascending order
 */
htm_list_t *htm_list_and(const htm_list_t *l1, const htm_list_t *l2);

/**
 * Compute the union between two lists of HTM identifiers
 *
 * @param l1 First list of identifiers (must be sorted in ascending order)
 * @param l2 Second list of identifiers (must be sorted in ascending order)
 *
 * @return A newly allocated list containing the union of l1 and l2, sorted in ascending order
 */
htm_list_t *htm_list_or(const htm_list_t *l1, const htm_list_t *l2);

/**
 * Destroy a list of HTM identifiers
 *
 * @param list The list of HTM identifiers
 */
void htm_list_destroy(htm_list_t *list);

/**
 * Retrieve spherical coordinates from three dimensional cartesian coordinates
 *
 * @param ra Right ascension (in [0, 2 pi]) [rad]
 * @param dec Declination (in [-pi/2, pi/2]) [rad]
 * @param v Input three dimensional cartesian coordinates (must have norm of 1)
 */
void htm_sphere_coord(double *ra, double *dec, const double *v);

/**
 * Convert an HTM identifier into a readable, null-terminated string
 *
 * @param s Destination character array (must be at least of size 3 + htm_level(id))
 * @param id The HMT identifier to convert (in [HTMID_MIN, HTMID_MAX])
 *
 * @return s
 */
char *htm_string(char *s, const htmid id);

/**
 * Retrieve the vertices of the HTM triangle that are associated with the
 * provided identifier
 *
 * @param v0 Three dimensional cartesian coordinate of the first vertex of the HTM triangle
 * @param v1 Three dimensional cartesian coordinate of the second vertex of the HTM triangle
 * @param v2 Three dimensional cartesian coordinate of the third vertex of the HTM triangle
 * @param id The HTM identifier
 *
 * @return -1 if level > HTML_LEVEL_MAX; 0 otherwise
 */
int htm_vertices(double *v0, double *v1, double *v2, const htmid id);

#endif
