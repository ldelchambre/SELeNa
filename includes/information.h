#ifndef _INFORMATION_H_
#define _INFORMATION_H_

#include <stddef.h>

/**
 * Structure designed to contain an event count
 */
typedef struct information_event_count_t {
  int event; // The event identifier
  size_t count; // The number of times this event appeared
} information_event_count;

/**
 * Count the number of distinct events from an experiment results.
 *
 * @param[out] nevent The number of distinct event
 * @param[in] y The resperiment results
 * @param[in] n The number of samples within y
 *
 *  Note: The returned value must MUST BE FREED using the free() function.
 *
 * @return For each distinct event, the number of times the latter appears in the experiment result (size=nevent), or NULL upon error. Retuned list is sorted according to the event value.
 */
information_event_count *information_count(size_t *nevent, const int *y, const size_t n);

/**
 * Compute the mode (most frequent value) of the given experiment results.
 *
 * @param[out] mode The returned mode
 * @param[out] p If not NULL, filled with the ratio of samples resulting in the most frequent event
 * @param[in] y The experiment result
 * @param[in] n The number of samples within y
 *
 * @return a zero value upon success.
 */
int information_mode(int *mode, double *p, const int *y, const size_t n);

/**
 * Compute the natural entropy of the given experiment results.
 *
 * @param[out] entropy The natural entropy within the subset [nat]
 * @param[in] y The experiment result
 * @param[in] n The number of samples within y
 *
 * Natural entropy can be converted to another base entropy, b, by dividing
 * the natural entropy with log(b).
 *
 * Example:
 *  The Shannon entropy (base 2 entropy) is given by:
 * > information_entropy(entropy, y, n);
 * > entropy /= M_LN2;
 *
 * @return a zero value upon success.
 */
int information_entropy(double *entropy, const int *y, const size_t n);

/**
 * Compute the Gini impurity measurement of the given experiment results.
 *
 * @param[out] gini The Gini impurity measure
 * @param[in] y The experiment result
 * @param[in] n The number of samples within y
 *
 * @return a zero value upon success.
 */
int information_gini(double *gini, const int *y, const size_t n);

#endif
