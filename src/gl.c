/**
 * @file gl.c
 *
 * Gravitational lens inversion algorithm, source file.
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
#include <gl.h>
#include <math.h>
#include <nelder.h>
#include <pheap.h>
#include <raster.h>
#include <string.h>
#include <stdlib.h>

/**
 * If _GL_USE_RASTERIZATION is defined, then the tesselation scheme will use
 * a rasterization (i.e. a 'pixelization') of the source plane in order to
 * fasten the search of the triangle in the image plane which corresponds to the
 * triangle in the image plane. This rasterization induce a small overhead in
 * both time and memory consumption when building the caustic curves. It may
 * hence be de-activated by commenting the line below (e.g. we are building a
 * lot of caustics curves through gl_build_caustics but call gl_img only once
 * or twice for each of these configurations). At the exception of the
 * above-mentioned example, _GL_USE_RASTERIZATION should be defined by default.
 */
#define _GL_USE_RASTERIZATION

/**
 * If fabs(_gl_edge(p0, p1, x, y)) < _GL_EDGE_TOLERANCE, then (x,y) is
 * considered to stand on the line joining p0 to p1. This prevent numerical
 * error to stem the detection of images. A by-side effect is that a unique
 * source position can be associated with multiple adjacent triangles, leading
 * to the appearance of phantom images (see however _GL_DIST_TOLERANCE for a
 * solution to this problem).
 */
#define _GL_EDGE_TOLERANCE 1e-15

/**
 * The minimal separation that should be present between two images in order to
 * consider them as being two distinct images.
 *
 * That is, if _gl_sqdist(img1, img2) <= _GL_DIST_TOLERANCE,
 * then img1 and img2 are multiple occurrences of the same image arising from
 * numerical instabilities and one of these is a phantom image.
 */
#define _GL_DIST_TOLERANCE 1e-10

/**
 * The size of the memory chunk to use during memory allocations (in byes).
 * Allocating memory through chunks allow to lessen the number of call to
 * malloc as a way to speed up the execution time.
 *
 * @see pheap.h header file for more informations
 */
#define _GL_PHEAP_CHUNK_SIZE (1U << 24)

/**
 * If _GL_DIST_TOLERANCE is defined, then the function gl_img maintain a list of
 * refined image positions so as to check whether these are unique or not. Given
 * that we do not know in advance the number of images we will get, these should
 * be added to a dynamically allocated structure.
 *   Still, we know in advance that the number of images will rarely exceed some
 * given value, that is relatively small. The idea is then to store a fixed number
 * of images, _GL_IMG_RES_SIZE, on the stack of the gl_img function. If this number
 * is exceeded, then a new chunk of _GL_IMG_RES_SIZE image positions in allocated
 * in a dynamically allocated structure an so on.
 * This approach allow:
 *   1) To protect at best against the errors that may occur in the gl_img function
 *      as most of the time, no heap allocation will be performed (only stack will
 *      be accessed).
 *   2) To speed up processing in multi-threaded environment as each malloc/free
 *      call use a global lock.
 * Note that _GL_IMG_RES_SIZE should never be lower than 1.
 */
#define _GL_IMG_RES_SIZE 8

#ifdef __GNUC__
#define UNUSED(param) __attribute__((unused)) param
#else
#define UNUSED(param) param
#endif

////////////////////////////////////////////////////////////////////////////////
// Type definition
////////////////////////////////////////////////////////////////////////////////
/**
 * Structure containing the triangle vertices
 */
typedef struct _gl_triangle {
  double v0[2]; ///< The first vertex
  double v1[2]; ///< The second vertex
  double v2[2]; ///< The third vertex
} _gl_triangle_t;

/**
 * Pointer to a _gl_caustic_node structure. This typedef aims to ease the
 * understanding of the pointers manipulation.
 */
typedef struct _gl_caustic_node *_gl_caustic_node_ptr_t;

/**
 * Structure containing the mapping between triangles in the source plane and
 * triangles in the image plane
 */
typedef struct _gl_caustic_node {
  _gl_triangle_t src;          ///< The triangle projected in the source plane
  _gl_triangle_t img;          ///< The triangle from the image plane
  _gl_caustic_node_ptr_t next; ///< The next node
} _gl_caustic_node_t;

struct gl_caustics {
  pheap_t *heap; ///< The persistent heap that should be used in order to allocate the _gl_caustic_node structures

#ifdef _GL_USE_RASTERIZATION
  double xfrom;    ///< The starting x position in the rasterization scheme
  double xto;      ///< The ending x position in the rasterization scheme
  unsigned int nx; ///< The number of x intervals we want to use in the rasterization scheme
  double dx;       ///< The sampling in the x position [ dx = (xto - xfrom) / nx ]

  double yfrom;    ///< The starting y position in the rasterization scheme
  double yto;      ///< The ending y position in the rasterization scheme
  unsigned int ny; ///< The number of y intervals we want to use in the rasterization scheme
  double dy;       ///< The sampling in the y position [ dy = (yto - yfrom) / ny ]

  _gl_caustic_node_ptr_t *map; ///< An array of pointer of _gl_caustic_node structures used in the rasterization scheme
  _gl_caustic_node_ptr_t outliers; ///< The set of triangle not contained in the rasterization scheme
#else
  _gl_caustic_node_ptr_t nodes; ///< The set of triangles from the tesselation scheme
#endif
};

/**
 * Argument of the _gl_build_caustics function
 */
typedef struct _gl_build_caustics_arg {
  gl_caustics_t *caustics; ///< The caustics structure to be built

  void (*alpha)(double *, const double *, void *); ///< The function giving the two dimensional deflection angle
  double (*mu)(const double *, void *); ///< The amplification function or NULL
  void *arg; ///< The argument to pass to the alpha and mu functions

  double eps; ///< The separation in the image plane below which we stop the recursive splitting of the caustic tessels

} _gl_build_caustics_arg_t;

#ifdef _GL_USE_RASTERIZATION
/**
 * Argument of the _gl_raster_add function
 */
typedef struct _gl_raster_add_arg {
  gl_caustics_t *caustics; ///< The caustics structure in which we will add the mapping between the src and img triangle
  const _gl_triangle_t *src; ///< The sources positions
  const _gl_triangle_t *img; ///< The images positions
} _gl_raster_add_arg_t;
#endif

/**
 * Argument of the _gl_nelder_sdf and _gl_nelder_end functions
 */
typedef struct _gl_nelder_arg {
  void (*alpha)(double *, const double *, void *); ///< The function giving the two dimensional deflection angle
  void *alpha_arg; ///< The parameters to pass to the alpha function
  const double *src; ///< The source position
  const double sdf; ///< The objective square deviation function we want to reach
  unsigned int iter; ///< The current iteration (Nelder-Mead algorithm systematically stop after GL_NELDER_MAXITER iterations so as to prevent infinite loop)
} _gl_nelder_arg_t;

////////////////////////////////////////////////////////////////////////////////
// Macro(s)
////////////////////////////////////////////////////////////////////////////////

/**
 * The edge function return a positive value if the point (x,y) stand on the
 * left of the line joining the point p0 and p1, a negative value if it stand
 * on the right and zero if is stand on the line.
 *
 * Example:
 * If                       (x,y)
 *                            X
 *      p0 X-----------------------------X p1
 * then _gl_edge(p0,p1,x,y) > 0.
 *
 * @param p0 Two dimensional coordinates of the first point constituting the line
 * @param p1 Two dimensional coordinates of the second point constituting the line
 * @param x  X coordinates of the point
 * @param y  Y coordinates of the point
 *
 * @return A positive value if (x,y) stand on the left of the line joining p0 to p1, a negative value if it stands on the right and zero if it stands on the line
 */
#define _gl_edge(p0,p1,x,y) ( ( y - p0[1] ) * ( p1[0] - p0[0] ) - ( x - p0[0] ) * ( p1[1] - p0[1] ) )

/**
 * Compute the squared distance between two vertices v0 and v1
 *
 * @param v0 The first vertex
 * @param v1 The second vertex
 *
 * @return (v1[0]-v0[0])^2 + (v1[1]-v0[1])^2
 */
#define _gl_sqdist(v0,v1) ((v1[0]-v0[0])*(v1[0]-v0[0])+(v1[1]-v0[1])*(v1[1]-v0[1]))

////////////////////////////////////////////////////////////////////////////////
// Private function prototypes
////////////////////////////////////////////////////////////////////////////////
/**
 * Compute the mapping between a triangle in the image plane and the
 * corresponding triangle in the source plane. If the triangle projected in the
 * source plane overlap a caustic curve, it is recursively split until a
 * precision of barg->eps in the image plane is reached.
 *
 * @param barg The argument that are necessary in order to build the mapping
 * @param img  The triangle for which the mapping has to be computed
 *
 * @return 0 upon success, -1 if memory allocation fails
 */
int _gl_build_caustics(_gl_build_caustics_arg_t *barg,
    const _gl_triangle_t *img);

/**
 * Check whether we should recursively split the triangle img. The condition for
 * splitting a triangle is that the amplification function, mu, is not NULL,
 * that it returns amplifications of opposite sign for two of the image vertices
 * (i.e. the triangle overlap a caustic curve) and that the distance between
 * each vertices of the triangle stand to be higher than eps.
 *
 * @param img The triangle for which the split criterion has to be evaluated
 * @param mu  The signed amplification function, if NULL then triangle is never split
 * @param arg The argument to pass to the amplification function
 * @param eps The maximal distance between the triangle vertices below which we stop the recursive splitting
 *
 * @return 1 if the triangle img has to be recursively splitted, 0 otherwise
 */
int _gl_build_caustics_split(const _gl_triangle_t *img,
    double (*mu)(const double *, void *), void *arg, const double eps);

/**
 * Add a mapping between a triangle in the source plane and a triangle in the image plane in the caustics structure.
 *
 * @param caustics The caustics structure in which to add the mapping
 * @param src The projection of the triangle img in the source plane
 * @param img The triangle in the image plane
 *
 * @return 0 upon success, -1 if memory allocation fails
 */
int _gl_caustic_add(gl_caustics_t *caustics, const _gl_triangle_t *src,
    const _gl_triangle_t *img);

/**
 * Add a node containing the mapping between source triangle and image triangle in a list of caustic nodes
 *
 * @param heap The persistent heap used to do the memory allocation
 * @param list A pointer to the list where to add this node, node will becomes the first item of the list
 * @param src  The source triangle
 * @param img  The image triangle
 *
 * @return 0 upon success; -1 if the memory allocation fails
 */
int _gl_caustic_node_add(pheap_t *heap, _gl_caustic_node_ptr_t *list,
    const _gl_triangle_t *src, const _gl_triangle_t *img);

/**
 * Function called to check the convergence of the Nelder-Mead refinement of the
 * image positions.
 *  Convergence is reached when the squared deviation function (sdf) of the
 * image against its projection in the source plane falls below arg->sdf or
 * after GL_NELDER_MAXITER iterations of the Nelder-Mead algorithm.
 *
 * @param img The image positions computed in the last iteration of the Nelder-Mead algorithm (not used here)
 * @param sdf The square deviation function associated with each of the image positions
 * @param arg Pointer to a _gl_nelder_arg_t structure (contain the objective sdf and the actual number of iteration already performed)
 *
 * @return 1 if the Nelder-Mead algorithm should stop iterating, 0 otherwise
 */
int _gl_nelder_end(const double UNUSED(*img), const double *sdf, void *arg);

/**
 * Compute the squared deviation function (sdf) of the images position in order
 * to minimize the latter through the use of a Nelder-Mead algorithm.
 *
 * @param img The image position
 * @param arg Pointer to a _gl_nelder_arg_t structure containing the arguments that are necessary in order to compute the sdf.
 *
 * @return The squared deviation function of the image against the source position contained in arg
 *
 * @see @gl_sdf function.
 */
double _gl_nelder_sdf(const double *img, void *arg);

#ifdef _GL_USE_RASTERIZATION
/**
 * Add the mapping between the image vertices contained in arg->img and the
 * corresponding vertices projected in the source plane (arg->src) into the
 * two dimensional array contained in arg->caustics->map[ix][iy].
 *
 * @param ix  The x position index in which we will add this entry (coming from a rasterization scheme)
 * @param iy  The y position index in which we will add this entry (coming from a rasterization scheme)
 * @param arg Pointer to a _gl_raster_add_arg_t structure containing the arguments that are necessary for adding the mapping.
 *
 * @return 0 upon success, -1 if memory allocation fails
 */
int _gl_raster_add(const unsigned int ix, const unsigned int iy, void *arg);
#endif

/**
 * Get the source position based on the image position and on the deflection angle function.
 *
 * @param src   The source position
 * @param img   The image position
 * @param alpha The deflection angle function
 * @param arg   The argument to pass to the deflection angle function
 */
void _gl_src(double *src, const double *img,
    void (*alpha)(double *, const double *, void *), void *arg);

/**
 * Check whether the point p is contained within the triangle t.
 *
 * @param t The triangle coordinates
 * @param p Two dimensional coordinates of the point
 *
 * @return 1 if p is contained in t, 0 otherwise
 */
int _gl_triangle_contain(const _gl_triangle_t *t, const double p[2]);

////////////////////////////////////////////////////////////////////////////////
// Private function implementation
////////////////////////////////////////////////////////////////////////////////
int _gl_build_caustics(_gl_build_caustics_arg_t *barg,
    const _gl_triangle_t *img) {
  _gl_triangle_t t;

  // Determine whether we should recursively split this triangle
  if (_gl_build_caustics_split(img, barg->mu, barg->arg, barg->eps)) {
    // v0' = v2
    memcpy(t.v0, img->v2, 2 * sizeof(double));

    // Compute the mean position between v0 and v1
    // v2' = 0.5 * ( v0 + v1 )
    t.v2[0] = 0.5 * (img->v0[0] + img->v1[0]);
    t.v2[1] = 0.5 * (img->v0[1] + img->v1[1]);

    // Process triangle (v0', v0, v2')
    memcpy(t.v1, img->v0, 2 * sizeof(double));
    if (_gl_build_caustics(barg, &t))
      return -1;

    // Process triangle (v0', v1, v2')
    memcpy(t.v1, img->v1, 2 * sizeof(double));
    if (_gl_build_caustics(barg, &t))
      return -1;
  } else { // If we don't have to recursively split this triangle
    // Get the source positions based on the image positions
    _gl_src(t.v0, img->v0, barg->alpha, barg->arg);
    _gl_src(t.v1, img->v1, barg->alpha, barg->arg);
    _gl_src(t.v2, img->v2, barg->alpha, barg->arg);

    // Add this node to the caustic structure
    if (_gl_caustic_add(barg->caustics, &t, img))
      return -1;
  }

  return 0;
}

int _gl_build_caustics_split(const _gl_triangle_t *img,
    double (*mu)(const double *, void *), void *arg, const double eps) {
  // If the amplification function is not provided, don't split the node
  if (!mu)
    return 0;

  // Don't split node if the maximal separation between any two point of the
  // triangle in the image plane is less than eps
  {
    const double eps2 = eps * eps;
    const double d01 = _gl_sqdist(img->v0, img->v1);
    const double d02 = _gl_sqdist(img->v0, img->v2);
    const double d12 = _gl_sqdist(img->v1, img->v2);

    if (d01 < eps2 && d02 < eps2 && d12 < eps2)
      return 0;
  }

  // Don't split the node if the (signed) amplification of all vertices of the
  // triangle have the same sign (i.e. we are presumably far from a caustic curve)
  {
    const double mu0 = mu(img->v0, arg);
    const double mu1 = mu(img->v1, arg);
    const double mu2 = mu(img->v2, arg);

    if (signbit(mu0) == signbit(mu1) && signbit(mu1) == signbit(mu2))
      return 0;
  }

  // Otherwise, split the node
  return 1;
}

int _gl_caustic_add(gl_caustics_t *caustics, const _gl_triangle_t *src,
    const _gl_triangle_t *img) {
#ifdef _GL_USE_RASTERIZATION
  _gl_raster_add_arg_t arg = { caustics, src, img };

  // If src partially stands outside the grid defined by
  // (xfrom, xto, yfrom, yto) then add it to the list of outliers
  if (src->v0[0] < caustics->xfrom || src->v0[1] < caustics->yfrom
      || src->v1[0] < caustics->xfrom || src->v1[1] < caustics->yfrom
      || src->v2[0] < caustics->xfrom || src->v2[1] < caustics->yfrom
      || src->v0[0] > caustics->xto || src->v0[1] > caustics->yto
      || src->v1[0] > caustics->xto || src->v1[1] > caustics->yto
      || src->v2[0] > caustics->xto || src->v2[1] > caustics->yto) {
    if (_gl_caustic_node_add(caustics->heap, &caustics->outliers, src, img))
      return -1;
  }

  // Add this entry to caustics->map using a rasterization scheme
  if (raster(src->v0, src->v1, src->v2, caustics->xfrom, caustics->xto,
      caustics->nx, caustics->yfrom, caustics->yto, caustics->ny,
      _gl_raster_add, &arg) == -1U)
    return -1;
  return 0;
#else
  return _gl_caustic_node_add(caustics->heap, &caustics->nodes, src, img);
#endif
}

int _gl_caustic_node_add(pheap_t *heap, _gl_caustic_node_ptr_t *list,
    const _gl_triangle_t *src, const _gl_triangle_t *img) {
  // Allocate the caustic node
  _gl_caustic_node_ptr_t node = (_gl_caustic_node_ptr_t) pheap_malloc(heap,
      sizeof(_gl_caustic_node_t));

  // If the node allocation fails, return -1
  if (!node)
    return -1;

  // Copy src and img into node->src and node->img
  memcpy(&node->src, src, sizeof(_gl_triangle_t));
  memcpy(&node->img, img, sizeof(_gl_triangle_t));

  // Add this node to the start of the list
  node->next = *list;
  *list = node;

  return 0;
}

int _gl_nelder_end(const double UNUSED(*img), const double *sdf, void *arg) {
  _gl_nelder_arg_t *p = (_gl_nelder_arg_t *) arg;

  // If one of the computed square deviation function reached the desired
  // precision, then stop the Nelder-Mead algorithm
  if (sdf[0] <= p->sdf || sdf[1] <= p->sdf || sdf[2] <= p->sdf)
    return 1;

  // If we reached the maximum number of iteration that is allowed, then stop
  // the Nelder-Mead algorithm
  if (p->iter++ == GL_NELDER_MAXITER)
    return 1;

  // Otherwise, continue until convergence is reached
  return 0;
}

double _gl_nelder_sdf(const double *img, void *arg) {
  const _gl_nelder_arg_t *p = (const _gl_nelder_arg_t *) arg;

  // Wrapper around the gl_sdf function that allows it to be used within the
  // Nelder-Mead algorithm
  return gl_sdf(p->src, img, p->alpha, p->alpha_arg);
}

#ifdef _GL_USE_RASTERIZATION
int _gl_raster_add(const unsigned int ix, const unsigned int iy, void *_arg) {
  _gl_raster_add_arg_t *arg = (_gl_raster_add_arg_t *) _arg;

  // Re-interepret map so as to ease understanding
  _gl_caustic_node_ptr_t (*map)[arg->caustics->ny] =
      (_gl_caustic_node_ptr_t (*)[arg->caustics->ny]) arg->caustics->map;

  // Add this entry to map[ix][iy]
  if (_gl_caustic_node_add(arg->caustics->heap, &map[ix][iy], arg->src, arg->img))
    return -1;

  return 0;
}
#endif

void _gl_src(double *src, const double *img,
    void (*alpha)(double *, const double *, void *), void *arg) {
  double a[2];
  // Get the deflection angle
  alpha(a, img, arg);
  // Get the source position from the deflection angle and image position
  src[0] = img[0] - a[0];
  src[1] = img[1] - a[1];
}

int _gl_triangle_contain(const _gl_triangle_t *t, const double p[2]) {
#ifdef _GL_EDGE_TOLERANCE
  const double d01 = _gl_edge(t->v0, t->v1, p[0], p[1]);
  const double d12 = _gl_edge(t->v1, t->v2, p[0], p[1]);
  const double d20 = _gl_edge(t->v2, t->v0, p[0], p[1]);
  const double d01p = fabs(d01);
  const double d12p = fabs(d12);
  const double d20p = fabs(d20);
  const int s01 = (
      (d01p <= _GL_EDGE_TOLERANCE) ?
          signbit((d12p > d20p) ? d12 : d20) : signbit(d01));
  const int s12 = (
      (d12p <= _GL_EDGE_TOLERANCE) ?
          signbit((d01p > d20p) ? d01 : d20) : signbit(d12));
  const int s20 = (
      (d20p <= _GL_EDGE_TOLERANCE) ?
          signbit((d01p > d12p) ? d01 : d12) : signbit(d20));
  return (s01 == s12 && s12 == s20);
#else
  int s01, s12, s20;

  s01 = signbit((double) _gl_edge(t->v0, t->v1, p[0], p[1]));
  s12 = signbit((double) _gl_edge(t->v1, t->v2, p[0], p[1]));

  if(s01 == s12) {
    s20 = signbit((double) _gl_edge(t->v2, t->v0, p[0], p[1]));
    return (s12 == s20);
  }
  return 0;
#endif
}

////////////////////////////////////////////////////////////////////////////////
// Public function implementation
////////////////////////////////////////////////////////////////////////////////
gl_caustics_t *gl_build_caustics(const double xfrom, const double xto,
    const unsigned int nx, const double yfrom, const double yto,
    const unsigned int ny, const double eps,
    void (*alpha)(double *, const double *, void *),
    double (*mu)(const double *, void *), void *arg) {
  // The caustics structure to return
  gl_caustics_t *caustics;

  // Compute the sampling of the image plane in the x and y directions
  const double dx = (xto - xfrom) / nx;
  const double dy = (yto - yfrom) / ny;

  // The argument we will pass to the recursive building of the caustic curves
  _gl_build_caustics_arg_t barg;

  // The actual triangle we will process
  _gl_triangle_t img;

  // Iterators over the x and y positions
  unsigned long ix, iy;

  // Allocate the caustics structure
  caustics = (gl_caustics_t *) malloc(sizeof(gl_caustics_t));
  if (!caustics)
    return NULL;

  // Initialize the caustics structure
  caustics->heap = pheap_init(_GL_PHEAP_CHUNK_SIZE);
  if(!caustics->heap) {
    gl_destroy_caustics(caustics);
    return NULL;
  }
#ifdef _GL_USE_RASTERIZATION
  caustics->map = (_gl_caustic_node_ptr_t *) pheap_calloc(caustics->heap, nx * ny, sizeof(_gl_caustic_node_ptr_t));
  if (!caustics->map) {
    gl_destroy_caustics(caustics);
    return NULL;
  }
  caustics->outliers = NULL;
  caustics->xfrom = xfrom;
  caustics->xto = xto;
  caustics->dx = dx;
  caustics->nx = nx;
  caustics->yfrom = yfrom;
  caustics->yto = yto;
  caustics->dy = dy;
  caustics->ny = ny;
#else
  caustics->nodes = NULL;
#endif

  // Initialize the build arguments
  barg.caustics = caustics;
  barg.alpha = alpha;
  barg.mu = mu;
  barg.arg = arg;
  barg.eps = eps;

  // Browse the primary grid
  for (ix = 0; ix < nx; ix++) {
    img.v0[0] = xfrom + ix * dx;
    img.v1[0] = img.v0[0] + dx;
    for (iy = 0; iy < ny; iy++) {
      img.v0[1] = yfrom + iy * dy;
      img.v1[1] = img.v0[1] + dy;
      // v2 *--* v1
      //    | /
      // v0 *
      img.v2[0] = img.v0[0];
      img.v2[1] = img.v1[1];
      if (_gl_build_caustics(&barg, &img)) {
        gl_destroy_caustics(caustics);
        return NULL;
      }
      //       * v1
      //     / |
      // v0 *--* v2
      img.v2[0] = img.v1[0];
      img.v2[1] = img.v0[1];
      if (_gl_build_caustics(&barg, &img)) {
        gl_destroy_caustics(caustics);
        return NULL;
      }
    }
  }

  // Return the resulting caustic structure
  return caustics;
}

double gl_sdf(const double *src, const double *img,
    void (*alpha)(double *, const double *, void *), void *arg) {
  double esrc[2], dx, dy;

  // Compute the expected source position
  _gl_src(esrc, img, alpha, arg);

  // Compute the error in x and y
  dx = src[0] - esrc[0];
  dy = src[1] - esrc[1];

  // Return the squared deviation function
  return dx * dx + dy * dy;
}

unsigned int gl_img(const double *src, const gl_caustics_t *caustics,
    void (*alpha)(double *, const double *, void *), void *alpha_arg,
    const double sdf, void (*process)(const double *, const double *, void *),
    void *process_arg) {
  // Initialize the arguments of the Nelder-Mead algorithm, designed to refine
  // the image positions
  _gl_nelder_arg_t p = { alpha, alpha_arg, src, sdf, 0 };
#ifdef _GL_DIST_TOLERANCE
  // Linked list of refined images positions
  // Note that res stand on stack so as to avoid heap allocation
  // most of the time
  struct _gl_img_res {
    struct _gl_img_res *next;
    double img[_GL_IMG_RES_SIZE][2];
    double sdf[_GL_IMG_RES_SIZE];
  } res = {NULL}, *res_node;
  // The position within the linked list of refined image positions
  unsigned int ires = 0, ires_node;
  // The squared deviation function that is associated with the last computed
  // image position
  double _sdf;
#else
  double UNUSED(_sdf);
#endif
  // Pointer to the list of nodes to browse
  _gl_caustic_node_ptr_t node;
  // Array containing the positions of the images vertices (3 points in two dimensions)
  double img[6];
  // The number of image we found so far
  unsigned int nimg = 0;

#ifdef _GL_USE_RASTERIZATION
  // If the source position is contained within the rasterization grid
  // i.e. The rasterization grid is taken as the grid covering the image plane
  if (caustics->xfrom <= src[0] && src[0] <= caustics->xto
      && caustics->yfrom <= src[1] && src[1] <= caustics->yto) {
    // Re-interpret map so as to ease understanding
    _gl_caustic_node_ptr_t (*map)[caustics->ny] =
        (_gl_caustic_node_ptr_t (*)[caustics->ny]) caustics->map;
    // Compute the indices in the rasterization grid where the source stand
    const unsigned int ix = (src[0] - caustics->xfrom) / caustics->dx;
    const unsigned int iy = (src[1] - caustics->yfrom) / caustics->dy;
    // Get the list of nodes standing in this resterization 'pixel'
    node = map[ix][iy];
  } else
    node = caustics->outliers;
#else
  // If the rasterization is not used, then we have to browse the whole set of
  // mapping between the source triangle and the image triangle
  node = caustics->nodes;
#endif

  // Browse the node list
  while (node) {
    // If the source is contained within the projection of the image triangle
    // into the source plane, then we found an image
    if (_gl_triangle_contain(&node->src, src)) {
      // If the alpha function is provided, then refine the image positions through
      // the use of a Nelder-Mead algorithm
      if (alpha) {
        // Set the initial simplex of the Nelder-Mead algorithm (i.e. the
        // triangle from the image plane)
        img[0] = node->img.v0[0];
        img[1] = node->img.v0[1];
        img[2] = node->img.v1[0];
        img[3] = node->img.v1[1];
        img[4] = node->img.v2[0];
        img[5] = node->img.v2[1];

        // Refine the image position through a Nelder-Mead algorithm.
        // Convergence is reached when the squared deviation function of one of
        // the image falls belows sdf or when we reached GL_NELDER_MAXITER
        // iterations
        p.iter = 0; // Initialize the number of iterations performed so far
        _sdf = nelder2(img, 2, _gl_nelder_sdf, &p,
                       NELDER_DEFAULT_ALPHA, NELDER_DEFAULT_BETA, NELDER_DEFAULT_GAMMA,
                       NELDER_DEFAULT_DELTA, _gl_nelder_end, &p);

      } else { // Otherwise, take the centroid of all points as the image position
        img[0] = (node->img.v0[0] + node->img.v1[0] + node->img.v2[0]) / 3;
        img[1] = (node->img.v0[1] + node->img.v1[1] + node->img.v2[1]) / 3;
      }

#ifdef _GL_DIST_TOLERANCE
      // Check whether the image we found is a duplicate or not (this may happen
      // because of the finite resolution of the tiling we used and result in
      // the production of phantom images)
      res_node = &res;
      ires_node = 0;
      while (res_node
             && (res_node->next || ires_node < ires)
             && _GL_DIST_TOLERANCE < _gl_sqdist(res_node->img[ires_node], img)) {
        if(++ires_node == _GL_IMG_RES_SIZE) {
          res_node  = res_node->next;
          ires_node = 0;
        }
      }

      // If this is not a duplicate, add it to the list of results
      if (res_node->next == NULL && ires_node == ires) {
        // Fill the node informations
        res_node->img[ires][0] = img[0];
        res_node->img[ires][1] = img[1];
        res_node->sdf[ires]    = _sdf;

        // If we have to create a new chunk
        if(++ires == _GL_IMG_RES_SIZE) {
          // Allocate the result node
          res_node->next = (struct _gl_img_res *) malloc(sizeof(struct _gl_img_res));
          // If the node allocation fails
          if(!res_node->next) {
            // Delete the result list
            while(res.next) {
              res_node = res.next->next;
              free(res.next);
              res.next = res_node;
            }
            return -1U;
          }
          // Initialize rest_node
          res_node->next->next = NULL;
          // Reset ires to zero
          ires = 0;
        }
      } else if(_sdf < res_node->sdf[ires_node]) {
        // Update the result node with the image position having the minimal
        // squared deviation function
        res_node->img[ires_node][0] = img[0];
        res_node->img[ires_node][1] = img[1];
        res_node->sdf[ires_node]    = _sdf;
      }
#else
      // Process the image positions
      if (process)
        process(src, img, process_arg);
      nimg++;
#endif
    }
    // Explore the next node
    node = node->next;
  }

#ifdef _GL_DIST_TOLERANCE
  // Process the image positions
  res_node = &res;
  ires_node = 0;
  while (res_node->next || ires_node < ires) {
    if (process)
      process(src, res_node->img[ires_node], process_arg);
    if(++ires_node == _GL_IMG_RES_SIZE) {
      res_node = res_node->next;
      ires_node = 0;
    }
    nimg++;
  }
  // Delete the result list
  while(res.next) {
    res_node = res.next->next;
    free(res.next);
    res.next = res_node;
  }
#endif

  // Return the number of images we found
  return nimg;
}

void gl_destroy_caustics(gl_caustics_t *caustics) {
  // If caustics is NULL, then nothing has to be done
  if(!caustics)
    return;

  // Destroy the persistent heap and free all memory
  // regions that were allocated through it (i.e. caustics->map,
  // caustics->outliers, caustics->nodes)
  pheap_destroy(caustics->heap);

  // Free the caustics structure
  free(caustics);
}
