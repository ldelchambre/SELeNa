/**
 * @file gl.h
 * 
 * Gravitational lens inversion algorithm, header file.
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
#ifndef _GL_H_
#define _GL_H_

/**
 * The maximal number of iteration to perform in order to refine the images
 * positions through a Nelder-Mead algorithm.
 */
#define GL_NELDER_MAXITER 1024

/**
 * The caustic structure, this structure allows to easily and rapidly get the
 * approximate image positions corresponding to a given source position.
 * 
 * @see @gl_build_caustics for obtaining this structure
 */
typedef struct gl_caustics gl_caustics_t;

/**
 * Build a caustic structure that is associated with a given lens model 
 * specified by the deflection angle function, alpha, the signed 
 * amplification function, mu, a grid in the image plane specified by xfrom,
 * xto, nx, yfrom, yto, ny and a precision eps.
 * 
 *  The function of the deflection angle, alpha, must have a prototype of the 
 * form:
 * > void alpha(double *alpha, const double *img, void *arg);
 * where alpha is the two dimensional deflection angle we have to compute,
 * img is the image position and where arg is a user-provided argument 
 * corresponding to the arg parameter of the gl_build_caustics function.
 * 
 *  The amplification function, mu, on his hand must have a prototype of the 
 * form
 * > double mu(const double *img, void *arg);
 * where img and arg are defined in a similar way than in the alpha function.
 * Note that the amplification function must ideally return amplification of 
 * opposite sign when the source crosses a caustic curve.
 * 
 *  This function create a tesselation of the image plane according to the
 * follwoing grid:
 *                       nx     yto
 *           xfrom *---+-...+---* xto
 *                 |\  |    |\  |
 *                 | \ |    | \ |
 *                 |  \|    |  \|
 *                 +---+-...+---+ 
 *                 .            .
 *                 .            . ny
 *                 .            .
 *                 +---+-...+---+
 *                 |\  |    |\  |
 *                 | \ |    | \ |
 *                 |  \|    |  \|
 *                 *---+-...+---*
 *                              yfrom.
 * Each triangle from this grid is projected back into the source plane so as to
 * enable the approximation of the image position based on the source position. 
 * If the triangle projected into the source plane overlap a caustic curve (i.e.
 * the amplification function return amplification of different signs for some 
 * of the triangle vertices), then this triangle is recursively split until the 
 * maximal separation between any two vertices in the image plane falls below 
 * eps. If the amplification function, mu, is NULL, then no recursive splitting 
 * will be performed.
 * 
 * @param xfrom The starting x position of the grid in the image plane
 * @param xto   The ending x position of the grid in the image plane
 * @param nx    The number of intervals in the x direction
 * @param yfrom The starting y position of the grid in the image plane
 * @param yto   The ending y position of the grid in the image plane
 * @param ny    The number of intervals in the y direction
 * @param eps   The maximal separation that is allowed when resursively splitting triangles
 * @param alpha The two dimensional deflection function
 * @param mu    The amplification function (if NULL, no recursive splitting will occur)
 * @param arg   The argument to pass to the alpha and mu function
 * 
 * @return A pointer to a newly allocated gl_caustics_t structure, or NULL if some memory allocation fails. This pointer must be freed with @gl_destroy_caustics.
 */
gl_caustics_t *gl_build_caustics(
                         const double xfrom, const double xto, const unsigned int nx,
                         const double yfrom, const double yto, const unsigned int ny,
                         const double eps,
                         void (*alpha)(double *, const double *, void *),
                         double (*mu)(const double *, void *),
                         void *arg);

/**
 * Compute the squared deviation function at a given image position, img, 
 * regarding a given source position, src, using the provided deflection angle
 * function.
 * 
 *  The squared deviation function is defined as 
 * | src - (img - alpha(img)) |^2
 * where src is the two dimensional position of the source, img, the two 
 * dimensional position of the image and alpha the two dimensional deflection
 * angle associated with the image position.
 * 
 * @param src   The source position
 * @param img   The image position
 * @param alpha The deflection angle function (see @gl_build_caustics for a more complete description)
 * @param arg   Arguments to pass to the alpha function
 * 
 * @return The squared deviation function, | src - (img - alpha(img)) |^2
 */
double gl_sdf(const double *src, const double *img,
              void (*alpha)(double *, const double *, void *),
              void *arg);

/**
 * Get the images position corresponding to the source position src using the
 * caustic structure caustics.
 * 
 *  This function uses the pre-computed caustic structure, caustics, in order to
 * obtain the apporoximate image position corresponding to the source position,
 * src. If the deflection angle function, alpha, is not NULL, then these images
 * positions are further refined through the use of a Nelder-Mead optimization
 * algorithm. The images positions are then returned through the use of a 
 * callback function, process, having a prototype of the form
 * > void process(const double *src, const double *img, void *arg);
 * where src is the source position, img is the found image position and arg is 
 * a user-provided argument corresponding to process_arg. The image position are
 * refined until reaching a specified squared deviation function, sdf, as 
 * returned by the @gl_sdf function or after GL_NELDER_MAXITER iterations of the
 * Nelder-Mead algorithm.
 * 
 * @param src         The source position for which to retrieve image(s)
 * @param caustics    The caustic structure used for obtaining the approximate image positions.
 * @param alpha       The deflection angle function (see @gl_build_caustics for a more complete description). If NULL, then no refinement will be performed.
 * @param alpha_arg   The Arguments to pass to the alpha function
 * @param sdf         The objective squared deviation function we would like to obtain upon refinement.
 * @param process     The callback function for obtaining the image positions (can be NULL if not needed).
 * @param process_arg The arguments to pass to the callback function, process.
 * 
 * @return The number of images we found for the given source position or -1U if some memory allocation failed.
 */
unsigned int gl_img(const double *src, const gl_caustics_t *caustics,
                    void (*alpha)(double *, const double *, void *), void *alpha_arg, 
                    const double sdf,
                    void (*process)(const double *, const double *, void *), void *process_arg);

/**
 * Destroy a caustic structure
 * 
 * @param caustics The caustic structure to destroy
 */
void gl_destroy_caustics(gl_caustics_t *caustics);

#endif
