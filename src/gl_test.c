/**
 * @file gl_test.c
 * 
 * Gravitational lens inversion algorithm, test file.
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
#include <gl.h>
#include <nsieg.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NTEST 32
#define SDF_TOLERANCE 1e-15
#define MU_MAX 1e6

typedef struct nsieg_arg {
  double b;
  double e;
  double theta_e;
  double gamma;
  double theta_gamma;
  double s;
} nsieg_arg_t;

void alpha_wrapper(double *alpha, const double *theta, void *arg) {
  nsieg_arg_t *p = (nsieg_arg_t *) arg;
  nsieg_alpha(alpha, theta, p->b, p->e, p->theta_e, p->gamma, p->theta_gamma, p->s);
}

double mu_wrapper(const double *theta, void *arg) {
  nsieg_arg_t *p = (nsieg_arg_t *) arg;
  return nsieg_mu(theta, p->b, p->e, p->theta_e, p->gamma, p->theta_gamma, p->s);
}

void process_img(const double *src, const double *img, void *arg) {
  const nsieg_arg_t *p = (nsieg_arg_t *) arg;
  const double sdf = gl_sdf(src, img, alpha_wrapper, arg);
  const double mu = nsieg_mu(img, p->b, p->e, p->theta_e, p->gamma, p->theta_gamma, p->s);

  if(sdf > SDF_TOLERANCE && fabs(mu) < MU_MAX) {
    fprintf(stderr, "Predicted image position (%g, %g) has sdf=%g (mu=%g)\n", img[0], img[1], sdf, mu);
    fprintf(stderr, "\tb=%g, e=%g, theta_e=%g, gamma=%g, theta_gamma=%g, s=%g, theta_s=(%g, %g)\n",
                      p->b, p->e, p->theta_e, p->gamma, p->theta_gamma, p->s, src[0], src[1]);
  }
}

int main(void) {
  unsigned long itest;
  unsigned long iobs;
  gl_caustics_t *caustics;
  nsieg_arg_t p;
  double theta_s[2];
  unsigned int nimg;
  unsigned int stats[] = {0,0,0,0,0,0,0};

  // In the NSIEg lens model we use, b is a normalization factor such that we
  // can arbitrarily set it to 1 for simplicity
  p.b = 1.;

  for(itest = 0; itest < NTEST; itest++) {
    // Pick up random parameters for the NSIEg lens model
    srand(itest+5);
    p.e = 0.95 * rand() / RAND_MAX;
    p.theta_e = M_PI * rand() / RAND_MAX;
    p.gamma = 0.3 * rand() / RAND_MAX;
    p.theta_gamma = M_PI * rand() / RAND_MAX;
    p.s = 0.5 * rand() / RAND_MAX;

    // Build the caustic curve
    caustics = gl_build_caustics(-3, 3, 1000, -3, 3, 1000, 1e-2,
                                 alpha_wrapper, mu_wrapper, &p);
    if(!caustics) {
      fprintf(stderr, "Unable to create caustic curves for lens model\n");
      fprintf(stderr, "\tb=%g, e=%g, theta_e=%g, gamma=%g, theta_gamma=%g, s=%g\n",
                        p.b, p.e, p.theta_e, p.gamma, p.theta_gamma, p.s);
      exit(EXIT_FAILURE);
    }

    // Probe 1024 randomly drawn source positions
    for(iobs = 0; iobs < (1U << 10); iobs++) {
      theta_s[0] = 2.*rand()/RAND_MAX-1.;
      theta_s[1] = 2.*rand()/RAND_MAX-1.;
      nimg = gl_img(theta_s, caustics, alpha_wrapper, &p, 1e-15, process_img, &p);
      if(nimg == -1U) {
        fprintf(stderr, "Memory allocation failed while retrieving the image positions\n");
        fprintf(stderr, "\tb=%g, e=%g, theta_e=%g, gamma=%g, theta_gamma=%g, s=%g\n",
                                p.b, p.e, p.theta_e, p.gamma, p.theta_gamma, p.s);
        exit(EXIT_FAILURE);
      }
      stats[(nimg > 6) ? 6 : nimg]++;
    }

    // Destroy the caustic curves
    gl_destroy_caustics(caustics);
  }

  // Print statistics result
  for(iobs = 0; iobs < 7; iobs++)
    fprintf(stderr, "%u lens configurations have%s%lu image%s\n",
                    stats[iobs],
                    (iobs == 6) ? " more than " : " ",
                    (iobs > 5) ? 5 : iobs,
                    (iobs > 1) ? "s" : "");

  exit(EXIT_SUCCESS);
}
