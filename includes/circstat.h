#ifndef _CIRCSTAT_H_
#define _CIRCSTAT_H_

#include <stddef.h>

/**
* Compute the mean direction of a sample of circular data.
*
* @param[in] x The sample a angles [rad]
* @param[in] n The number of samples
*
* @return The mean direction of the provided samples
*/
double circstat_mean(const double *x, const size_t n);

/**
* Compute mean resultant vector length for circular data.
*
* @param[in] x The sample a angles [rad]
* @param[in] n The number of samples
*
* @return The mean resultant vector length of the provided samples
*/
double circstat_r(const double *x, const size_t n);

/**
* Compute the variance of a sample of circular data
*
* @param[in] x The sample a angles [rad]
* @param[in] n The number of samples
*
* @return The circular variances of the provided sample
*/
double circstat_var(const double *x, const size_t n);

/**
 * Compute the differences between 2 angles.
 *
 * @param[in] a First angle
 * @param[in] b Second angle
 *
 * @return The difference a-b normalized between -pi and pi.
 */
double circstat_diff(const double a, const double b);

/**
 * Compute the correlation coefficient between 2 samples of circular data.
 *
 * @param[in] x The first sample a angles [rad]
 * @param[in] x The second sample a angles [rad]
 * @param[in] n The number of samples
 *
 * @return The correlation coefficient between the 2 samples.
 */
double circstat_corr(const double *x, const double *y, const size_t n);

/**
 * Return the mean absolute error between 2 samples of curcular data.
 *
 * @param[in] x The first sample a angles [rad]
 * @param[in] x The second sample a angles [rad]
 * @param[in] n The number of samples
 *
 * @return The mean absolute error between the 2 samples.
 */
double circstat_abserr(const double *x, const double *y, const size_t n);

/**
 * Compute the confidence limit on the mean of a sample of circular data.
 *
 * @param[in] x The sample a angles [rad]
 * @param[in] n The number of samples
 * @param[in] sigma The sigma value at which the returned interval shoud stands.
 *
 * Examples:
 *   > double m = circstat_mean(x,n);
 *   > double dm = circstat_confmean(x,n,2.);
 *   will retrieve the mean value, m, and a confidence interval, dm, that is such
 *   that there is 95.45% probability (2-sigma) that the true mean stands in the
 *   interval m±dm.
 *
 * @return The confidence interval associated with the computed mean value or -1 if the confidence interval cannot be reliably computed.
 */
double circstat_confmean(const double *x, const size_t n, const double sigma);

/**
 * Estimate the kappa parameter of the Von Mises distribution based on the
 * provided samples.
 *
 * @param[in] x The array of samples
 * @param[in] n The size of the sample
 *
 * @return The kappa parameter of the underlying Von Mises distribution
 */
double circstat_kappa(const double *x, const size_t n);

#endif
