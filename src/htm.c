/**
 * @file htm.c
 * 
 * Hierarchical Triangular Mesh, source file.
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
#include <string.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////
// Macros
////////////////////////////////////////////////////////////////////////////////

/**
 * Numerical precision we require in order for a point to belong to
 * an HTM triangle
 */
#define _HTM_EPSILON 1e-15 // 2.0626e-4 [µas]

/**
 * Compute the cross product of the three dimensional vectors u and v
 *
 * @param dest The three dimensional destination vector (dest != u, dest != v)
 * @param u Three dimensional coordinates of the first vertex (dest != u)
 * @param v Three dimensional coordinates of the second vertex (dest != v)
 */
#define _htm_cross_product(dest,u,v) \
  do {dest[0]=u[1]*v[2]-u[2]*v[1];dest[1]=u[2]*v[0]-u[0]*v[2];dest[2]=u[0]*v[1]-u[1]*v[0];} while(0)

/**
 * Compute the dot product of the three dimensional vectors u and v
 *
 * @param u Three dimensional coordinates of the first vertex
 * @param v Three dimensional coordinates of the second vertex
 *
 * @return The dot product of u and v
 */
#define _htm_dot_product(u,v) \
  (u[0]*v[0]+u[1]*v[1]+u[2]*v[2])

/**
 * Compute the triple product of the three dimensional vectors u, v and p
 *
 * ( u X v ) . p
 *
 * @param u Three dimensional coordinates of the first vertex (dest != u)
 * @param v Three dimensional coordinates of the second vertex (dest != v)
 */
#define _htm_triple_product(u,v,p) \
  (p[0]*(u[1]*v[2]-u[2]*v[1])+p[1]*(u[2]*v[0]-u[0]*v[2])+p[2]*(u[0]*v[1]-u[1]*v[0]))

////////////////////////////////////////////////////////////////////////////////
// Type definitions
////////////////////////////////////////////////////////////////////////////////

/**
 * Structure designed to ease the building of HTM lists by keeping both
 * the start and end nodes
 */
typedef struct _htm_list {
  htm_list_t *start; ///< The start node
  htm_list_t *end; ///< The end node
} _htm_list_t;

////////////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////////////

/**
 * The initial vertices composing the top level HTM triangles
 */
static const double _htm_vtop[][3] = {
  { 0.,  0.,  1.},
  { 1.,  0.,  0.},
  { 0.,  1.,  0.},
  {-1.,  0.,  0.},
  { 0., -1.,  0.},
  { 0.,  0., -1.}
};

/**
 * The top level HTM triangles
 */
static const struct {
  char name[2]; ///< Name
  unsigned char id; ///< Identifier
  unsigned char iv0; ///< Index of the first vertex of the triangle
  unsigned char iv1; ///< Index of the second vertex of the triangle
  unsigned char iv2; ///< Index of the third vertex of the triangle
} _htm_ttop[] = {
/* 0 */ {"S2", 10, 3, 5, 4}, // x < 0 & y < 0 & z < 0
/* 1 */ {"N1", 13, 4, 0, 3}, // x < 0 & y < 0 & z > 0
/* 2 */ {"S1",  9, 2, 5, 3}, // x < 0 & y > 0 & z < 0
/* 3 */ {"N2", 14, 3, 0, 2}, // x < 0 & y > 0 & z > 0
/* 4 */ {"S3", 11, 4, 5, 1}, // x > 0 & y < 0 & z < 0
/* 5 */ {"N0", 12, 1, 0, 4}, // x > 0 & y < 0 & z > 0
/* 6 */ {"S0",  8, 1, 5, 2}, // x > 0 & y > 0 & z < 0
/* 7 */ {"N3", 15, 2, 0, 1}  // x > 0 & y > 0 & z > 0
};

/**
 * Mapping between the top level HTM identifier and indices within @_htm_ttop
 */
static const unsigned char _htm_ittop[] = {
/* 0 */ 6,
/* 1 */ 2,
/* 2 */ 0,
/* 3 */ 4,
/* 4 */ 5,
/* 5 */ 1,
/* 6 */ 3,
/* 7 */ 7
};

////////////////////////////////////////////////////////////////////////////////
// Private functions prototypes
////////////////////////////////////////////////////////////////////////////////
/**
 * Build the list of HTM identifiers overlapping the provided constraint
 *
 * @param l The resulting list of HTM identifiers
 * @param p The three dimensional cartesian coordinates of the direction of the constraint (must have a norm of 1)
 * @param cosd The cosine of the search radius of the constraints (in [-1, 1])
 * @param id The current HTM identifier associated with the triangle (v0,v1,v2)
 * @param v0 Three dimensional cartesian coordinate of the first vertex of the HTM triangle
 * @param v1 Three dimensional cartesian coordinate of the second vertex of the HTM triangle
 * @param v2 Three dimensional cartesian coordinate of the third vertex of the HTM triangle
 * @param level The current HTM recursion level
 * @param maxlevel The maximal HTM recursion level to use (in [0, HTML_LEVEL_MAX])
 *
 * @return 0 If the list was successfully updated, -1 upon error(s)
 */
int _htm_constraint(_htm_list_t *l, const double *p, const double cosd, const htmid id, const double *v0, const double *v1, const double *v2, const unsigned int level, const unsigned int maxlevel);

/**
 * Check whether the great circle joining points vi and vj intersect the
 * constraint (half-plane) defined by p and cosd.
 *
 * @param p The three dimensional cartesian coordinates of the direction of the constrain (must have a norm of 1)
 * @param cosd The cosine of the search radius to use (in [-1,1])
 * @param vi Three dimensional coordinates of the first point
 * @param vj Three dimensional coordinates of the second point
 *
 * @see Szalay, A.S. et al., 2007, "Indexing the Sphere with the Hierarchical Triangular Mesh", eprint arXiv:cs/0701164
 *
 * @return 1 If the great circle and the half plane intersect, 0 otherwise
 */
int _htm_constraint_intersect_gcirc(const double *p, const double cosd, const double *vi, const double *vj);

/**
 * Check whether a constrain (p, cosd) intersect a given HTM triangle.
 *
 * @param v0 Three dimensional cartesian coordinate of the first vertex of the HTM triangle
 * @param v1 Three dimensional cartesian coordinate of the second vertex of the HTM triangle
 * @param v2 Three dimensional cartesian coordinate of the third vertex of the HTM triangle
 * @param p Cartesian coordinates of the direction of the constraint (must have norm of 1)
 * @param cosd Cosine of the search radius
 *
 * @return 0 if the constraint do not intersect the triangle, 1 if it partially \
 * overlap or 2 if the triangle is entirey contained in it
 */
int _htm_constraint_intersect_triangle(const double *v0, const double *v1, const double *v2, const double *p, const double cosd);

/**
 * Check whether v is contained within the triangle defined by the three
 * dimensional cartesian vertices (v0, v1, v2), given in counter
 * clockwise order.
 *
 * @param v Three dimensional coordinates of the point
 * @param v0 Three dimensional coordinates of the first vertex
 * @param v1 Three dimensional coordinates of the second vertex
 * @param v2 Three dimensional coordinates of the third vertex
 *
 * @return 1 if v is contained in (v0,v1,v2); 0 otherwise
 */
int _htm_is_contained(const double *v, const double *v0, const double *v1, const double *v2);

/**
 * Add an entry to the end of the HTM list
 *
 * @param list The list to modify
 * @param id The HTM identifier to add
 *
 * @return 0 if entry was successfully added, -1 if memory allocation fails
 */
int _htm_list_add(_htm_list_t *list, const htmid id);

/**
 * Compare two nodes of HTM identifiers
 *
 * @param l1 First node
 * @param l2 Second node
 *
 * @return -1 if l1 is considered to be lower than l2; 1 if l2 is considered as lower than l1; 0 if l1 and l2 overlaps
 */
int _htm_list_cmp(const htm_list_t *l1, const htm_list_t *l2);

/**
 * Compute the mid point between three dimensional cartesian coordinates
 *
 * @param v Three dimensional coordinates of the destination
 * @param v0 Three dimensional coordinates of the first vertex
 * @param v1 Three dimensional coordinates of the second vertex
 */
void _htm_midpoint(double *v, const double *v0, const double *v1);

////////////////////////////////////////////////////////////////////////////////
// Methods implementation
////////////////////////////////////////////////////////////////////////////////

int _htm_constraint(_htm_list_t *l, const double *p, const double cosd, const htmid id, const double *v0, const double *v1, const double *v2, const unsigned int level, const unsigned int maxlevel) {
  double w0[3], w1[3], w2[3];

  // Check whether this triangle overlap the constraint
  switch(_htm_constraint_intersect_triangle(v0,v1,v2,p,cosd)) {
    // If they do not overlap, then we are done
    case 0:
      break;
    // If they partially overlap
    case 1:
      // If we reached the maximal level
      if(level == maxlevel) {
        // Add this identifier to the end of the list
        if(_htm_list_add(l,id)) return -1;
      } else {
        // Compute the intermediate vertices
        _htm_midpoint(w0,v1,v2);
        _htm_midpoint(w1,v0,v2);
        _htm_midpoint(w2,v0,v1);
        // Recursively process each sub-triangle
        if(_htm_constraint(l, p, cosd, id << 2, v0, w2, w1, level + 1, maxlevel)) return -1;
        if(_htm_constraint(l, p, cosd, (id << 2) | 0x01LL, v1, w0, w2, level + 1, maxlevel)) return -1;
        if(_htm_constraint(l, p, cosd, (id << 2) | 0x02LL, v2, w1, w0, level + 1, maxlevel)) return -1;
        if(_htm_constraint(l, p, cosd, (id << 2) | 0x03LL, w0, w1, w2, level + 1, maxlevel)) return -1;
      }
      break;
    // If the triangle is entirely contained in the constraint
    case 2:
      // Add this identifier to the end of the list
      if(_htm_list_add(l,id)) return -1;
      break;
  }

  return 0;
}

int _htm_constraint_intersect_gcirc(const double *p, const double cosd, const double *vi, const double *vj) {
    double a, b, c, d, s1, s2, costij, u2, gi, gj;

    gi = _htm_dot_product(vi,p); // gi = cos(phi)
    gj = _htm_dot_product(vj,p); // gj = cos(varphi)
    costij = _htm_dot_product(vi,vj); // costij = cos(theta)

    // u2 = tan(theta/2)^2
    u2 = ( 1. - costij ) / ( 1. + costij );

    // Quadratic equation a s^2 + b s + c = 0
    a = -u2 * ( gi + cosd );
    b = gi * ( u2 - 1. ) + gj * ( u2 + 1. );
    c = gi - cosd;

    // If this equation is linear rather than quadratic
    if(a == 0) {
      // Solve this linear equation
      s1 = -c / b;
      return (0. <= s1 && s1 <= 1.);
    }

    // Compute discriminant b^2 - 4 a c
    d = b * b - 4. * a * c;

    // If no solution exist
    if(d < 0.)
      return 0;

    { // Solve the quadratic equation (see Press et al., 2002, §5.6)
      const double q = -0.5 * ( b + ((b < 0) ? -sqrt(d) : sqrt(d)));
      s1 = q / a;
      s2 = c / q;
    }

    // If one solution stads in [0, 1] then an intersection exists
    return (( 0. <= s1 && s1 <= 1. ) || ( 0. <= s2 && s2 <= 1. ) );
}

int _htm_constraint_intersect_triangle(const double *v0, const double *v1, const double *v2, const double *p, const double cosd) {
  // If the search angle exceed pi/2, consider the complementary problem
  if(cosd < 0.) {
    const double ip[3] = {-p[0], -p[1], -p[2]};
    return 2 - _htm_constraint_intersect_triangle(v0,v1,v2,ip,-cosd);
  }

  // Check if vertices are contained within the constraint
  if(_htm_dot_product(v0,p) >= cosd) {
    return (_htm_dot_product(v1,p) >= cosd && _htm_dot_product(v2,p) >= cosd) ? 2 : 1;
  } else if(_htm_dot_product(v1,p) >= cosd || _htm_dot_product(v2,p) >= cosd)
    return 1;

  // Check whether the constaint is entirely contained within the HTM triangle
  if(_htm_is_contained(p, v0, v1, v2))
    return 1;

  // Check if constraint intersect the triangle edges
  if(_htm_constraint_intersect_gcirc(p,cosd,v0,v1)
     || _htm_constraint_intersect_gcirc(p,cosd,v1,v2)
     || _htm_constraint_intersect_gcirc(p,cosd,v0,v2))
    return 1;
  return 0;
}

int _htm_is_contained(const double *v, const double *v0, const double *v1, const double *v2) {
  if(_htm_triple_product(v0,v1,v) < -_HTM_EPSILON) return 0;
  if(_htm_triple_product(v1,v2,v) < -_HTM_EPSILON) return 0;
  if(_htm_triple_product(v2,v0,v) < -_HTM_EPSILON) return 0;
  return 1;
}

int _htm_list_add(_htm_list_t *list, const htmid id) {
  // Allocate the new node
  if(!list->start) {
    list->start = list->end = (htm_list_t *) malloc(sizeof(htm_list_t));
  } else {
    list->end->next = (htm_list_t *) malloc(sizeof(htm_list_t));
    list->end = list->end->next;
  }

  // If allocation fails, destroy the list then return -1
  if(!list->end) {
    htm_list_destroy(list->start);
    return -1;
  }

  // Assign values to node
  list->end->id = id;
  list->end->next = NULL;

  return 0;
}

int _htm_list_cmp(const htm_list_t *l1, const htm_list_t *l2) {
  return htm_cmp(l1->id, l2->id);
}

void _htm_midpoint(double *v, const double *v0, const double *v1) {
  double n;
  v[0] = v0[0] + v1[0];
  v[1] = v0[1] + v1[1];
  v[2] = v0[2] + v1[2];
  n = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if(n > 0.) {
    v[0] /= n;
    v[1] /= n;
    v[2] /= n;
  }
}

/**
 * Unachieved function designed to check whether a constaint and the
 * bounding circle of a HTM triangle overlap. Its interest is
 * questionable given the large amount of operations implied. It was
 * however present in the original algorithm.
 *
int _htm_bounding_circle_intersect(const double *p, const double cosd, const double *v0, const double *v1, const double *v2) {
    double v10[3],v21[3],vb[3];
    double cosdb, cost, n;

    // vb = ((v1-v0) X (v2-v1)) / |(v1-v0) X (v2-v1)|
    v10[0] = v1[0]-v0[0]; v10[1] = v1[1]-v0[1]; v10[2] = v1[2]-v0[2];
    v21[0] = v2[0]-v1[0]; v21[1] = v2[1]-v1[1]; v21[2] = v2[2]-v1[2];
    _htm_cross_product(vb,v10,v21);
    n = sqrt(vb[0]*vb[0]+vb[1]*vb[1]+vb[2]*vb[2]);
    vb[0] /= n; vb[1] /= n; vb[2] /= n;

    // cosdb = v0 . vb
    cosdb = _htm_dot_product(v0,vb);

    // cost = vb . p
    cost = _htm_dot_product(vb, p);

    // Check that t <= d + db with t=acos(cost); d=acos(cosd); db=acos(db)
    // Ideally this should be done without using acos
    return (acos(cost) <= acos(cosd) + acos(cosdb));
}
*/

void htm_cart_coord(double *v, const double ra, const double dec) {
  const double cosd = cos(dec);
  v[0] = cos(ra) * cosd;
  v[1] = sin(ra) * cosd;
  v[2] = sin(dec);
}

int htm_cmp(const htmid id1, const htmid id2) {
  const unsigned int level1 = htm_level(id1);
  const unsigned int level2 = htm_level(id2);
  htmid diff;

  if(level1 < level2)
    diff = id1 - ( id2 >> 2 * ( level2 - level1 ) );
  else
    diff = ( id1 >> 2 * (level1-level2) ) - id2;

  return (diff < 0) ? -1 : ((diff > 0) ? 1 : 0);
}

htm_list_t *htm_constraint(const double *p, const double cosd, const unsigned int level) {
  _htm_list_t res = {NULL, NULL};
  double v0[3], v1[3], v2[3];
  unsigned int i;

  for(i = 0; i < 8; i++) {
    // Browse each top level HTM triangle in the order of their HTM identifier
    htm_copy_vertex(v0,_htm_vtop[_htm_ttop[_htm_ittop[i]].iv0]);
    htm_copy_vertex(v1,_htm_vtop[_htm_ttop[_htm_ittop[i]].iv1]);
    htm_copy_vertex(v2,_htm_vtop[_htm_ttop[_htm_ittop[i]].iv2]);

    // Recusrsively process this triangle
    if(_htm_constraint(&res, p, cosd, _htm_ttop[_htm_ittop[i]].id, v0, v1, v2, 0, (level < HTM_LEVEL_MAX) ? level : HTM_LEVEL_MAX)) {
      htm_list_destroy(res.start);
      return NULL;
    }
  }

  // Return the start node
  return res.start;
}

void htm_list_destroy(htm_list_t *list) {
  htm_list_t *tmp;

  while(list) {
    tmp = list->next;
    free(list);
    list = tmp;
  }
}

htmid htm_from_string(const char *s) {
  htmid res = 0;
  unsigned long i;

  switch(s[0]) {
    case 'N':
    case 'n':
      res += 0x3;
      break;
    case 'S':
    case 's':
      res += 0x2;
      break;
    default:
      return HTMID_BAD;
  }

  for(i = 1; s[i] != '\0'; i++) {
    res <<= 2;
    if('0' <= s[i] && s[i] <= '3')
      res += s[i] - '0';
    else
      return HTMID_BAD;
  }

  if(i - 2 > HTM_LEVEL_MAX)
    return HTMID_BAD;

  return res;
}

htmid htm_get(const double *v, const unsigned int level) {
  // Vertices defining a triangle on the sphere
  double v0[3], v1[3], v2[3];

    // Resultling HTM identifier
  htmid res = 0;

  // Intermidate vertices defining a triangle on the sphere
  double w0[3], w1[3], w2[3];

  // Index within the _htm_ttop array
  unsigned char iv;

  // Level counter
  unsigned int ilevel;

  // If the asked level exceed the maximal allowed level, return with error
  if(level > HTM_LEVEL_MAX)
    return HTMID_BAD;

  // Compute the index within the array of initial triangle
  iv = 4 * ( v[0] > 0 ) + 2 * ( v[1] > 0 ) + ( v[2] > 0);

  // Build the level 0 HTM identifier
  res = _htm_ttop[iv].id;

  // Build the level 0 vertices
  htm_copy_vertex(v0,_htm_vtop[_htm_ttop[iv].iv0]);
  htm_copy_vertex(v1,_htm_vtop[_htm_ttop[iv].iv1]);
  htm_copy_vertex(v2,_htm_vtop[_htm_ttop[iv].iv2]);

  // Browse each level
  for(ilevel = 0; ilevel < level; ilevel++) {
    // Shift the previous id by 2 bits
    res <<= 2;

    // Compute the intermediate vertices
    _htm_midpoint(w0,v1,v2);
    _htm_midpoint(w1,v0,v2);
    _htm_midpoint(w2,v0,v1);

    // Check in which triangle the input point is contained
    if(_htm_is_contained(v,v0,w2,w1)) {
      htm_copy_vertex(v1,w2);
      htm_copy_vertex(v2,w1);
    } else if(_htm_is_contained(v,v1,w0,w2)) {
      res += 1;
      htm_copy_vertex(v0,v1);
      htm_copy_vertex(v1,w0);
      htm_copy_vertex(v2,w2);
    } else if(_htm_is_contained(v,v2,w1,w0)) {
      res += 2;
      htm_copy_vertex(v0,v2);
      htm_copy_vertex(v1,w1);
      htm_copy_vertex(v2,w0);
    } else /* if(_htm_is_contained(v,w0,w1,w2)) */ {
      res += 3;
      htm_copy_vertex(v0,w0);
      htm_copy_vertex(v1,w1);
      htm_copy_vertex(v2,w2);
    }
  }

  return res;
}

unsigned int htm_level(const htmid id) {
#ifdef __GNUC__
  return ((HTMID_BITS - 1 - __builtin_clzll((unsigned long long) id)) >> 1) - 1;
#else
#warning "Function __builtin_clzll not available, htm_level may hence experience a slow execution"
  unsigned int level;

  for(level = 0; (id >> level); level += 2)
    ;

  return (level >> 1) - 2;
#endif
}

htm_list_t *htm_list_and(const htm_list_t *l1, const htm_list_t *l2) {
  _htm_list_t res = {NULL, NULL};
  unsigned int level1, level2;

  // While both lists contain elements
  while(l1 && l2) {
    // Compare nodes from the two lists
    switch(_htm_list_cmp(l1,l2)) {
      case 0: // l1 and l2 overlaps
        // Get the level of l1 and l2
        level1 = htm_level(l1->id);
        level2 = htm_level(l2->id);
        // If l2 is contained into l1
        if(level1 < level2) {
          // Add l2 to the resulting list
          if(_htm_list_add(&res, l2->id)) {
            htm_list_destroy(res.start);
            return NULL;
          }
          // Advance in l2
          l2 = l2->next;
        // If l1 is contained into l2
        } else if(level1 > level2) {
          // Add l1 to the resulting list
          if(_htm_list_add(&res, l1->id)) {
            htm_list_destroy(res.start);
            return NULL;
          }
          // Advance in l1
          l1 = l1->next;
        } else { // If l1 == l2
          // Add l1 to the resulting list
          if(_htm_list_add(&res, l1->id)) {
            htm_list_destroy(res.start);
            return NULL;
          }
          // Advance in l1 and l2
          l1 = l1->next;
          l2 = l2->next;
        }
        break;
      case -1: // l1 is lower than l2
        l1 = l1->next; // Advance l1
        break;
      case 1: // l2 is lower than l1
        l2 = l2->next; // Advance l2
        break;
    }
  }

  return res.start;
}

htm_list_t *htm_list_or(const htm_list_t *l1, const htm_list_t *l2) {
  _htm_list_t res = {NULL, NULL};
  unsigned int level1, level2;

  // While both lists contain elements
  while(l1 && l2) {
    // Compare nodes from the two lists
    switch(_htm_list_cmp(l1,l2)) {
      case 0: // l1 and l2 overlaps
        // Get the level of l1 and l2
        level1 = htm_level(l1->id);
        level2 = htm_level(l2->id);
        // If l2 is contained into l1
        if(level1 < level2) {
          // Add l1 to the resulting list
          if(_htm_list_add(&res, l1->id)) {
            htm_list_destroy(res.start);
            return NULL;
          }
          // Advance in l2
          do {
            l2 = l2->next;
          } while(l2 && _htm_list_cmp(l1,l2) == 0);
          // Advance in l1
          l1 = l1->next;
        // If l1 is contained into l2
        } else if(level1 > level2) {
          // Add l2 to the resulting list
          if(_htm_list_add(&res, l2->id)) {
            htm_list_destroy(res.start);
            return NULL;
          }
          // Advance in l1
          do {
            l1 = l1->next;
          } while(l1 && _htm_list_cmp(l1,l2) == 0);
          // Advance in l2
          l2 = l2->next;
        } else { // If l1 == l2
          // Add l1 to the resulting list
          if(_htm_list_add(&res, l1->id)) {
            htm_list_destroy(res.start);
            return NULL;
          }
          // Advance in l1 and l2
          l1 = l1->next;
          l2 = l2->next;
        }
        break;
      case -1: // l1 is lower than l2
        // Add l1 to the resulting list
        if(_htm_list_add(&res, l1->id)) {
          htm_list_destroy(res.start);
          return NULL;
        }
        // Advance in l1
        l1 = l1->next;
        break;
      case 1: // l2 is lower than l1
        // Add l2 to the resulting list
        if(_htm_list_add(&res, l2->id)) {
          htm_list_destroy(res.start);
          return NULL;
        }
        // Advance in l2
        l2 = l2->next;
        break;
    }
  }

  // Add the rest of l1 to the list
  while(l1) {
    if(_htm_list_add(&res, l1->id)) {
      htm_list_destroy(res.start);
      return NULL;
    }
    l1 = l1->next;
  }

  // Add the rest of l2 to the list
  while(l2) {
    if(_htm_list_add(&res, l2->id)) {
      htm_list_destroy(res.start);
      return NULL;
    }
    l2 = l2->next;
  }

  return res.start;
}

void htm_sphere_coord(double *ra, double *dec, const double *v) {
  *ra = atan2(v[1],v[0]); //
  *dec = asin(v[2]);
  // Normalize ra to stand within [0,2 pi]
  if(*ra < 0.) *ra += 2.*M_PI;
}

char *htm_string(char *s, const htmid id) {
  const unsigned int level = htm_level(id);
  long long _id = id;
  unsigned int ilevel;

  if(id < 0) {
    s[0] = '?';
    s[1] = '\0';
  } else {
    s[2+level] = '\0';
    for(ilevel = 0; ilevel < level; ilevel++) {
      s[level-ilevel+1] = '0' + (_id & 0x03LL);
      _id >>= 2;
    }
    strncpy(s,_htm_ttop[_htm_ittop[_id & 0x07LL]].name,2);
  }

  return s;
}

int htm_vertices(double *v0, double *v1, double *v2, const htmid id) {
  const unsigned int level = htm_level(id);
  double w0[3], w1[3], w2[3];
  unsigned int ilevel;

  if(level > HTM_LEVEL_MAX)
    return -1;

  // Compute the initial set of vertices
  {
    const unsigned char iv = _htm_ittop[(id >> 2*level) & 0x07LL];
    htm_copy_vertex(v0,_htm_vtop[_htm_ttop[iv].iv0]);
    htm_copy_vertex(v1,_htm_vtop[_htm_ttop[iv].iv1]);
    htm_copy_vertex(v2,_htm_vtop[_htm_ttop[iv].iv2]);
  }

  // Recursively build vertices based on the HTM ID
  for(ilevel = level; ilevel; ilevel--) {
    const unsigned char it = (id >> 2*(ilevel-1)) & 0x03LL;
    switch(it) {
      case 0:
        _htm_midpoint(w1,v0,v2);
        _htm_midpoint(w2,v0,v1);
        htm_copy_vertex(v1,w2);
        htm_copy_vertex(v2,w1);
        break;
      case 1:
        _htm_midpoint(w0,v1,v2);
        _htm_midpoint(w2,v0,v1);
        htm_copy_vertex(v0,v1);
        htm_copy_vertex(v1,w0);
        htm_copy_vertex(v2,w2);
        break;
      case 2:
        _htm_midpoint(w0,v1,v2);
        _htm_midpoint(w1,v0,v2);
        htm_copy_vertex(v0,v2);
        htm_copy_vertex(v1,w1);
        htm_copy_vertex(v2,w0);
        break;
      case 3:
        _htm_midpoint(w0,v1,v2);
        _htm_midpoint(w1,v0,v2);
        _htm_midpoint(w2,v0,v1);
        htm_copy_vertex(v0,w0);
        htm_copy_vertex(v1,w1);
        htm_copy_vertex(v2,w2);
        break;
    }
  }

  return 0;
}
