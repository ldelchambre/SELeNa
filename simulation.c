/**
 * @file simulation.c
 * 
 * Simulate gravitational lens configurations based on a non-singular isothermal
 * ellipsoid lens model in presence of an external shear.
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
#include <getopt.h>
#include <gl.h>
#include <math.h>
#include <nsieg.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __GNUC__
#define UNUSED(param) __attribute__((unused)) param
#else
#define UNUSED(param) param
#endif

/*************************************************************************************************
 * Program parameters
 ************************************************************************************************/
#define DEFAULT_EPS 1e-3
#define DEFAULT_SDF 1e-31
#define DEFAULT_NTHREAD 1
#define DEFAULT_PREC 6
#define DEFAULT_NMIN 1
#define DEFAULT_NMAX -1U

unsigned int nmin = DEFAULT_NMIN;
unsigned int nmax = DEFAULT_NMAX;
double eps = DEFAULT_EPS;
double sdf = DEFAULT_SDF;
unsigned int prec = DEFAULT_PREC;
unsigned int nthread = DEFAULT_NTHREAD;
char noheader = 0;

typedef struct prange {
  double from;
  double to;
  unsigned int i;
  unsigned int n;
} prange_t;

prange_t b  = {0,0,0,0};
prange_t e  = {0,0,0,0};
prange_t te = {0,0,0,0};
prange_t g  = {0,0,0,0};
prange_t tg = {0,0,0,0};
prange_t s  = {0,0,0,0};
prange_t xs = {0,0,0,0};
prange_t ys = {0,0,0,0};
prange_t x  = {0,0,0,0};
prange_t y  = {0,0,0,0};

#define PARAM_NTHREAD  65536
#define PARAM_NOHEADER 65537

static struct option options_long[] = {
/*-n*/{ "nmin", required_argument, 0, 'n'},
/*-N*/{ "nmax", required_argument, 0, 'N'},
/*-e*/{ "eps", required_argument, 0, 'e'},
/*-s*/{ "sdf", required_argument, 0, 's'},
/*-p*/{ "prec", required_argument, 0, 'p'},
/*  */{ "nthread", required_argument, 0, PARAM_NTHREAD},
/*  */{ "noheader", no_argument, 0, PARAM_NOHEADER},
/*-?*/{ "help", no_argument, 0, '?' },
/*  */{ 0, 0, 0, 0 } };

static const char options_short[] = "n:N:e:s:p:?";

/*************************************************************************************************
 * Help message
 ************************************************************************************************/
static const char
    *help_message =
        "USAGE:\n"
        "  simulation [options] b e te g tg s xs ys x y\n"
        "\n"
        "DESCRIPTION:\n"
        "  Simulate gravitational lens configurations based on a non-singular isothermal\n"
        "  ellipsoid (NSIEg) lens model in presence of an external shear (Kormann, 1994).\n"
        "\n"
        "  NSIEg lens models are characterized by a projected mass density in the lens \n"
        "  plane of\n"
        "    kappa(x,y) = 0.5 * b / sqrt(s^2 + x^2 + y^2 / (1 - e)^2).\n"
        "  The presence of an external shear (Kovner,1987) add a supplemental contribution\n"
        "  to the deflection angle, alpha(x,y), of the form\n"
        "    alpha'(x,y) = g * | cos(2*tg)  sin(2*tg) | * alpha(x,y).\n"
        "                      | sin(2*tg) -cos(2*tg) |\n"
        "\n"
        "  The present program estimates the lensed images positions of the provided NSIEg \n"
        "  model(s) through the use of a tesselation of the image plane where each tessel is\n"
        "  propagated back into the source plane so as to approximate the lensed image\n"
        "  positions. The tesselation scheme we used consist in the subdivision of the\n"
        "  rectangles from the grid specified by the x and y ranges into a set of two triangles\n"
        "  each. Once one of these triangle overlap a caustic curve, it is then recursively\n"
        "  splitted until its largest side in the image plane falls below the limit given by\n"
        "  the eps option.\n"
        "  The lensed images positions are then refined through the use of a Nelder-Mead\n"
        "  optimization algorithm (Nelder, 1965), where the square deviation function\n"
        "  (Schramm & Kaiser, 1987) is minimized so as to falls below the limit given\n"
        "  by the sdf option. A maximal number of %u iterations of the Nelder-Mead\n"
        "  algorithm is then performed.\n"
        "\n"
        "  The lens model parameters are as follow:\n"
        "    b  The normalization factor of the lense. In the limit of a singular (s = 0) and\n"
        "       spherical (e = 0) lens model, b stand to be the Einstein radius of the lens.\n"
        "    e  The ellipticity of the mass distribution, or equivalently of the isophot contours.\n"
        "    te The inclination of the semi-major axis of the mass distribution measured East of north [deg].\n"
        "    g  The strength of the external shear.\n"
        "    tg The orientation of the external shear [deg].\n"
        "    s  The radius of the core of the mass distribution.\n"
        "    xs The x position of the source.\n"
        "    ys The y position of the source.\n"
        "    x  The range of x values for creating the grid in the image plane.\n"
        "    y  The range of y values for creating the grid in the image plane.\n"
        "\n"
        "  The aforementioned parameters can be specified either as scalar values or as\n"
        "  ranges of values in the form \"[from:n:to]\" where from is the starting value\n"
        "  of the parameter, to its ending value and n is the number of interval we have\n"
        "  to create (e.g. \"[0:5:1]\" is equivalent to the set of parameters {0, 0.2, 0.4,\n"
        "  0.6, 0.8, 1}). Note that the x and y parameters are forced to be ranges and not scalar.\n"
        "\n"
        "OPTIONS:\n"
        "  -n, --nmin,     Only output configurations for which the number of lensed images is above of equal to nmin [default=%u]\n"
        "  -N, --nmax,     Only output configurations for which the number of lensed images is below nmax [default=%u]\n"
        "  -e, --eps,      The maximal separation between the vertices of a tessel taken in the\n"
        "                  image plane, if the latter overlap a caustic curve [default=%g].\n"
        "  -s, --sdf,      The objective squared deviation function to use in order to assess\n"
        "                  the convergence of the Nelder-Mead algorithm [default=%g].\n"
        "  -p, --prec,     The desired output precision [default=%u]\n"
        "      --nthread   The number of thread to use [default=%u]\n"
        "      --noheader, Don't print header line\n"
        "  -?, --help,     Print this help message\n"
        "\n"
        "REFERENCES:\n"
        "  Kormann, R. et al, \"Isothermal elliptical gravitational lens models\", 1994, A&A, 284(1), 285.\n"
        "  Kovner, I., \"The quadrupole gravitational lens\", 1987, ApJ, 312, 22.\n"
        "  Nelder, J. A. and Mead, R., \"A simplex method for Function minimization\", 1965, The Computer Journal, 7(4), 308.\n"
        "  Schramm, T. and Kayser, R., \"A simple imaging procedure for gravitational lenses\", 1987, A&A, 174, 361.\n";

void printHelp(FILE *out) {
  fprintf(out, help_message, GL_NELDER_MAXITER, DEFAULT_NMIN, DEFAULT_NMAX, DEFAULT_EPS, DEFAULT_SDF, DEFAULT_PREC, DEFAULT_NTHREAD);
}

/*******************************************************************************
 * The global mutex
 ******************************************************************************/
pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;

/*******************************************************************************
 * Functions
 ******************************************************************************/
int get_prange(prange_t *p, const char *s) {
  char *send;
  
  // If s contains a range of value in the form [from:n:to], read that range
  if(*s == '[') {
    p->from = strtod(++s, &send);
    if(s == send || *(s = send) != ':')
      return -1;
    p->i = 0;
    p->n = strtoul(++s, &send, 10);
    if(s == send || *(s = send) != ':')
      return -1;
    p->to = strtod(++s, &send);
    if(s == send || *(s = send) != ']')
      return -1;
    s++;
  } else { // If s contains a scalar value, read it
    p->from = p->to = strtod(s,&send);
    if(s == send) // If nothing was read
      return -1;
    p->i = p->n = 0;
    s = send;
  }
  
  // If it remains something to be converted, this is an error
  if(*s != '\0')
    return -1;
  
  // If to < from, this is an error
  if(p->to < p->from)
    return -1;
  
  // Range of the type [from:0:to] with from != to are incorrect
  if(p->n == 0 && p->from != p->to)
    return -1;
  
  return 0;
}

typedef struct nsieg_arg {
  double b;
  double e;
  double te;
  double g;
  double tg;
  double s;
} nsieg_arg_t;

void alpha_wrapper(double *alpha, const double *img, void *arg) {
  nsieg_arg_t *p = (nsieg_arg_t *) arg;
  nsieg_alpha(alpha, img, p->b, p->e, p->te, p->g, p->tg, p->s);
}

double mu_wrapper(const double *img, void *arg) {
  nsieg_arg_t *p = (nsieg_arg_t *) arg;
  return nsieg_mu(img, p->b, p->e, p->te, p->g, p->tg, p->s);
}

#define prange_val(p) ((p.n>0)?p.from+(double)p.i*(p.to-p.from)/p.n:p.from)

int next_gl_arg(nsieg_arg_t *gl_arg) {
  // Acquire lock
  pthread_mutex_lock(&lock);

  // If some lens model parameters are still unexplored
  if(b.i <= b.n) {
    
    // Get the current valued of all NSIEg arguments, except source position
    gl_arg->b  = prange_val(b);
    gl_arg->e  = prange_val(e);
    gl_arg->te = prange_val(te);
    gl_arg->g  = prange_val(g);
    gl_arg->tg = prange_val(tg);
    gl_arg->s  = prange_val(s);
    
    // Update the value of the NSIEg arguments
    if(s.n < ++s.i) {
      s.i = 0;
      if((g.from == 0 && g.i == 0) || tg.n < ++tg.i) {
        tg.i = 0;
        if(g.n < ++g.i) {
          g.i = 0;
          if((e.from == 0 && e.i == 0) || te.n < ++te.i) {
            te.i = 0;
            if(e.n < ++e.i) {
              e.i = 0;
              ++b.i;
            }
          }
        }
      }
    }
    // Realease the lock
    pthread_mutex_unlock(&lock);

    return 0;
  }
  
  // Realease the lock
  pthread_mutex_unlock(&lock);
  return -1;
}

typedef struct gl_images  {
  double img[2];
  double mu;
  double sdf;
  struct gl_images *next;
} gl_images_t;

void simulation_thread_destroy_caustics(void *arg) {
  gl_destroy_caustics((gl_caustics_t *) arg);
}

void simulation_thread_destroy_images(void *arg) {
  gl_images_t *node = *((gl_images_t **) arg);
  gl_images_t *tmp;
  
  while(node) {
    tmp = node->next;
    free(node);
    node = tmp;
  }
}

typedef struct simulation_thread_add_arg {
  nsieg_arg_t *gl_arg;
  gl_images_t **list;
  int *error;
} simulation_thread_add_arg;

void simulation_thread_add(const double *src, const double *img, void *arg) {
  simulation_thread_add_arg *staa = (simulation_thread_add_arg *) arg;
  gl_images_t *prev, *next;
  
  // Allocate the image structure
  gl_images_t *node = (gl_images_t *) malloc(sizeof(gl_images_t));
  if(!node) {
    fprintf(stderr, "Unable to allocate the image node\n");
    *(staa->error) = 1;
    return;
  }
  
  // Fill the node informations
  node->img[0] = img[0];
  node->img[1] = img[1];
  node->mu     = fabs(mu_wrapper(img, staa->gl_arg));
  node->sdf    = gl_sdf(src, img, alpha_wrapper, staa->gl_arg);
  
  // Add this node to the list (highest amplification first)
  prev = NULL;
  next = *(staa->list);
  while(next && next->mu > node->mu) {
    prev = next;
    next = next->next;
  }
  if(!prev) // Insert in front of list
    *(staa->list) = node;
  else // Insert into list 
    prev->next = node;
  node->next = next;
}

void *simulation_thread(void *arg) {
  const double dx = (xs.n > 0) ? (xs.to - xs.from) / xs.n : 0;
  const double dy = (ys.n > 0) ? (ys.to - ys.from) / ys.n : 0;
  volatile int *error = (int *) arg; // This avoid to bother about unused parameter warning
  
  nsieg_arg_t gl_arg;
  gl_caustics_t *caustics;
  gl_images_t *node, *list = NULL;
  simulation_thread_add_arg staa;
  unsigned int nimg;
  
  unsigned int ix, iy;
  double src[2];
  
  // While no error occurs and that some lens configurations remains to be explored
  while(*error == 0 && next_gl_arg(&gl_arg) == 0) {
    // Build the caustic structure
    caustics = gl_build_caustics(x.from, x.to, x.n, y.from, y.to, y.n, eps, alpha_wrapper, mu_wrapper, &gl_arg);
    if(!caustics) {
      fprintf(stderr, "An error occured while building the caustics curves of configuration (b=%g, e=%g, te=%g, g=%g, tg=%g, s=%g)\n",
                      gl_arg.b, gl_arg.e, gl_arg.te*180./M_PI, gl_arg.g, gl_arg.tg*180./M_PI, gl_arg.s);
      *error = 1;
      pthread_exit((int *) error);
    }
    
    // Push the procedure designed to cleanup the caustic structure
    pthread_cleanup_push(simulation_thread_destroy_caustics, caustics);
    
    // Fill the argument to pass to the simulation_thread_add function
    staa.gl_arg = &gl_arg;
    staa.list   = &list;
    staa.error  = (int *) error;
    
    // Browse each source position
    for(ix = 0; ix <= xs.n && *error == 0; ix++) {
      src[0] = xs.from + dx * ix;
      for(iy = 0; iy <= ys.n && *error == 0; iy++) {
        src[1] = ys.from + dy * iy;
        
        // Push the procedure designed to cleanup the list of images
        pthread_cleanup_push(simulation_thread_destroy_images, &list);
        
        // Get the images positions, amplification and sdf
        nimg = gl_img(src, caustics, alpha_wrapper, &gl_arg, sdf, simulation_thread_add, &staa);
        if(nimg == 0)
          fprintf(stderr, "Warning: no image found for the configuration (b=%g, e=%g, te=%g, g=%g, tg=%g, s=%g, xs=%g, ys=%g)\n",
                      gl_arg.b, gl_arg.e, gl_arg.te*180./M_PI, gl_arg.g, gl_arg.tg*180./M_PI, gl_arg.s, src[0], src[1]);
        else if(nimg == -1U) {
          fprintf(stderr, "Memory allocation failed for the configuration (b=%g, e=%g, te=%g, g=%g, tg=%g, s=%g, xs=%g, ys=%g)\n",
                                          gl_arg.b, gl_arg.e, gl_arg.te*180./M_PI, gl_arg.g, gl_arg.tg*180./M_PI, gl_arg.s, src[0], src[1]);
          *error = 1;
        } else if(nmin <= nimg && nimg < nmax) {
          pthread_mutex_lock(&lock); // Acquire the lock
          node = list;
          while(node) {
            fprintf(stdout, "%g %g %g %g %g %g %g %g %u %+*.*e %+*.*e %+*.*e %+*.*e\n",
                             gl_arg.b, 
                             gl_arg.e, 
                             gl_arg.te*180./M_PI, 
                             gl_arg.g, 
                             gl_arg.tg*180./M_PI, 
                             gl_arg.s, 
                             src[0], 
                             src[1],
                             nimg, 
                             (prec == 0) ? 6 : prec+7, prec, node->img[0], 
                             (prec == 0) ? 6 : prec+7, prec, node->img[1], 
                             (prec == 0) ? 6 : prec+7, prec, -2.5*log10(node->mu), 
                             (prec == 0) ? 6 : prec+7, prec, node->sdf);
            node = node->next;
          }
          pthread_mutex_unlock(&lock); // Release the lock
        }
        
        // Clean the list of images
        pthread_cleanup_pop(1); // Equivalent to simulation_thread_destroy_images(list)
        list = NULL;
      }
    }
    
    // Clean the caustic structure
    pthread_cleanup_pop(1); // Equivalent to simulation_thread_destroy_caustics(caustics)
  }
  
  pthread_exit((int *) error);
}

/*************************************************************************************************
 * Main program
 ************************************************************************************************/
int main(int argc, char **argv) {
  int error = 0;
  int nopt;
  pthread_t *pool;
  unsigned int i, j;
  int opt;
  
  // If we don't have at least 10 arguments after the optional arguments, then exit
  if(argc < 11) {
    fprintf(stderr, "Not enough argument on input, required argument are (b e te g tg s xs ys x y)\n");
    printHelp(stderr);
    pthread_mutex_destroy(&lock);
    exit(EXIT_FAILURE);
  }
  
  // Count the number of options that were passed to the program
  nopt = argc - 10;

  // Get the program parameters
  while ((opt = getopt_long(nopt, argv, options_short, options_long, NULL)) != -1) {
    switch (opt) {
      case 'n':
        nmin = strtoul(optarg, NULL, 10);
        break;
      case 'N':
        nmax = strtoul(optarg, NULL, 10);
        break;
      case 'e':
        eps = strtod(optarg, NULL);
        break;
      case 's':
        sdf = strtod(optarg, NULL);
        break;
      case 'p':
        prec = strtoul(optarg, NULL, 10);
        break;
      case PARAM_NTHREAD:
        nthread = strtoul(optarg, NULL, 10);
        break;
      case PARAM_NOHEADER:
        noheader = 1;
        break;
      case '?':
        printHelp(stdout);
        pthread_mutex_destroy(&lock);
        exit(EXIT_SUCCESS);
        break;
      default:
        fprintf(stderr, "Unrecognized options \"%s\"\n", argv[optind]);
        printHelp(stderr);
        pthread_mutex_destroy(&lock);
        exit(EXIT_FAILURE);
        break;
    }
  }
  
  // Get the whole set of parameters
  if(get_prange(&b, argv[nopt++])) {
    fprintf(stderr, "Error while retrieving the b parameter: unrecognized argument \"%s\"\n", argv[optind-1]);
    pthread_mutex_destroy(&lock);
    exit(EXIT_FAILURE);
  }
  if(get_prange(&e, argv[nopt++])) {
    fprintf(stderr, "Error while retrieving the e parameter: unrecognized argument \"%s\"\n", argv[optind-1]);
    pthread_mutex_destroy(&lock);exit(EXIT_FAILURE);
  }
  if(get_prange(&te, argv[nopt++])) {
    fprintf(stderr, "Error while retrieving the te parameter: unrecognized argument \"%s\"\n", argv[optind-1]);
    pthread_mutex_destroy(&lock);exit(EXIT_FAILURE);
  }
  te.from *= M_PI/180;te.to *= M_PI/180;// te is given in degree, convert it back into radian
  if(get_prange(&g, argv[nopt++])) {
    fprintf(stderr, "Error while retrieving the g parameter: unrecognized argument \"%s\"\n", argv[optind-1]);
    pthread_mutex_destroy(&lock);exit(EXIT_FAILURE);
  }
  if(get_prange(&tg, argv[nopt++])) {
    fprintf(stderr, "Error while retrieving the tg parameter: unrecognized argument \"%s\"\n", argv[optind-1]);
    pthread_mutex_destroy(&lock);exit(EXIT_FAILURE);
  }
  tg.from *= M_PI/180;tg.to *= M_PI/180;// tg is given in degree, convert it back into radian
  if(get_prange(&s, argv[nopt++])) {
    fprintf(stderr, "Error while retrieving the s parameter: unrecognized argument \"%s\"\n", argv[optind-1]);
    pthread_mutex_destroy(&lock);exit(EXIT_FAILURE);
  }
  if(get_prange(&xs, argv[nopt++])) {
    fprintf(stderr, "Error while retrieving the xs parameter: unrecognized argument \"%s\"\n", argv[optind-1]);
    pthread_mutex_destroy(&lock);exit(EXIT_FAILURE);
  }
  if(get_prange(&ys, argv[nopt++])) {
    fprintf(stderr, "Error while retrieving the ys parameter: unrecognized argument \"%s\"\n", argv[optind-1]);
    pthread_mutex_destroy(&lock);exit(EXIT_FAILURE);
  }
  if(get_prange(&x, argv[nopt++])) {
    fprintf(stderr, "Error while retrieving the x parameter: unrecognized argument \"%s\"\n", argv[optind-1]);
    pthread_mutex_destroy(&lock);exit(EXIT_FAILURE);
  }
  if(get_prange(&y, argv[nopt++])) {
    fprintf(stderr, "Error while retrieving the y parameter: unrecognized argument \"%s\"\n", argv[optind-1]);
    pthread_mutex_destroy(&lock);exit(EXIT_FAILURE);
  }
  
  // Check that the x and y arguments are ranges
  if(x.from == x.to || y.from == y.to) {
    fprintf(stderr, "Error x and y parameters must be ranges\n");
    pthread_mutex_destroy(&lock);exit(EXIT_FAILURE);
  }
  
  // Check for extra parameters
  if (nopt < argc) {
    fprintf(stderr, "WARNING: non-option ARGV-elements found:");
    while (optind < argc)
      fprintf(stderr, " %s", argv[nopt++]);
    fprintf(stderr, "\n");
  }
  
  // Print the header line if asked to do so
  if(noheader == 0)
    fprintf(stdout, "b e te g tg s xs ys nimg x y mag sdf\n");

  // Allocate the thread pool
  pool = (pthread_t *) calloc(nthread, sizeof(pthread_t));
  if(!pool) {
    fprintf(stderr, "Unable to allocate a thread pool of size %u\n", nthread);
    pthread_mutex_destroy(&lock);
    exit(EXIT_FAILURE);
  }

  // Launch all threads
  for(i = 0; i < nthread; i++) {
    // If the thread failed to launch, 
    if(pthread_create(&pool[i], NULL, simulation_thread, &error)) {
      // Cancel all previously launched threads
      for(j = 0; j < i; j++)
        pthread_cancel(pool[j]);
      // Wait for all threads to finish
      for(j = 0; j < i; j++)
        pthread_join(pool[j], NULL);
      // Destroy the mutex then exit
      pthread_mutex_destroy(&lock);
      free(pool);
      exit(EXIT_FAILURE);
    }
  }
  
  // Wait for all threads to finish
  for(i = 0; i < nthread; i++)
    pthread_join(pool[i], NULL);
    
  // If error(s) occurred
  if(error) {
    fprintf(stderr, "Some error(s) occurred during processing\n");
    pthread_mutex_destroy(&lock);
    free(pool);
    exit(EXIT_FAILURE);
  }
  
  // Destroy the mutex
  pthread_mutex_destroy(&lock);
  
  // Destroy the thread pool
  free(pool);
  
  // Exit
  exit(EXIT_SUCCESS);
}
