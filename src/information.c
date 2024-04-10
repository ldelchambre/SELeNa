#include <information.h>

#include <stdlib.h>
#include <math.h>

/**
 * Structure designed to contains the count statistics of the event of an experiamnt
 */
struct _information_tree_count {
  struct _information_tree_count *left; // Left part of tree (left->event < event)
  struct _information_tree_count *right; // Right part of tree (event < right->event)
  int event; // The event at this node
  size_t count; // The number of times we encountered this event
};

/**
 * Create a new information tree node containing the event c
 *
 * @param c The class to add
 *
 * @return The newly created node or NULL upon failure
 */
struct _information_tree_count *_information_tree_count_new(const int event) {
  struct _information_tree_count *tc = (struct _information_tree_count *) malloc(sizeof(struct _information_tree_count));

  if(tc == NULL)
    return NULL;

  tc->left = tc->right = NULL;
  tc->event = event;
  tc->count = 1;

  return tc;
}

/**
 * Add an event to the information tree
 *
 * @param tc The information tree
 * @param event The event to add
 *
 * @return 0 is the event was already encountered within the experiment, -1 upon error, 1 otherwise
 */
int _information_tree_count_add(struct _information_tree_count *tc, const int event) {
  if(tc->event < event) {
    if(tc->left)
      return _information_tree_count_add(tc->left, event);
    if((tc->left = _information_tree_count_new(event)) == NULL)
      return -1;
    return 1; // return 1 to say that a new class has been added
  }
  if(event < tc->event) {
    if(tc->right)
      return _information_tree_count_add(tc->right, event);
    if((tc->right = _information_tree_count_new(event)) == NULL)
      return -1;
    return 1; // return 1 to say that a new class has been added
  }
  // else if(tc->c == c)
  tc->count++;
  return 0; // return 0 to say that no new class has been added
}

/**
 * Copy the information tree tc within the information count structure
 *
 * @param counts The destination information count structure
 * @param tc The source information counter tree
 *
 * @return A pointer to the next count structure where to insert any new information tree
 */
information_event_count *_information_tree_count_copy(information_event_count *counts, const struct _information_tree_count *tc) {
  if(tc == NULL) return counts;
  counts = _information_tree_count_copy(counts, tc->left);
  counts->event = tc->event;
  counts->count = tc->count;
  return _information_tree_count_copy(counts+1, tc->right);
}

/**
 * Delete the information counter tree
 *
 * @param tc the information counter tree to delete
 */
void _information_tree_count_free(struct _information_tree_count *tc) {
  if(tc->left) {
    _information_tree_count_free(tc->left);
    free(tc->left);
  }
  if(tc->right) {
    _information_tree_count_free(tc->right);
    free(tc->right);
  }
}

information_event_count *information_count(size_t *nevent, const int *y, const size_t n) {
  information_event_count *counts;
  struct _information_tree_count tc;
  int firstencountered;
  size_t i;

  // Check that we have at least one event
  if(n == 0) {
    *nevent = 0;
    return NULL;
  }

  // Initialize the information tree count with the first experiment result
  tc.left = tc.right = NULL;
  tc.event = y[0];
  tc.count = 1;
  *nevent = 1;

  // Add all classes to the tree count tree
  for(i = 1; i < n; i++) {
    // Add this class to the count tree
    firstencountered = _information_tree_count_add(&tc, y[i]);
    if(firstencountered == -1) {
      _information_tree_count_free(&tc);
      return NULL;
    }
    *nevent += firstencountered;
  }

  // Allocate the resulting count structure
  counts = (information_event_count *) calloc(*nevent,sizeof(information_event_count));
  if(counts == NULL) {
    _information_tree_count_free(&tc);
    return NULL;
  }

  // Copy the tree count within the count structre
  _information_tree_count_copy(counts, &tc);

  // Delete the tree count
  _information_tree_count_free(&tc);

  // Return the count structure
  return counts;
}

int information_mode(int *mode, double *p, const int *y, const size_t n) {
  information_event_count *counts;
  size_t nevent;
  size_t bestcount;
  size_t i;

  // Get the counts for each classes
  counts = information_count(&nevent, y, n);
  if(counts == NULL)
    return -1;

  // Find the mode of the given subset (most frequent value)
  bestcount = 0;
  for(i = 0; i < nevent; i++) {
    if(bestcount < counts[i].count) {
      bestcount = counts[i].count;
      *mode = counts[i].event;
    }
  }

  // Delete the detailed counts
  free(counts);

  // Set the portion of sampels belonging to this event
  if(p)
    *p = (double) bestcount / n;

  return 0;
}

int information_entropy(double *entropy, const int *y, const size_t n) {
  information_event_count *counts;
  size_t nevent;
  double p;
  size_t i;

  // Get the counts for each classes
  counts = information_count(&nevent, y, n);
  if(counts == NULL)
    return -1;

  // Compute the entropy of the given subset
  *entropy = 0;
  for(i = 0; i < nevent; i++) {
    p = (double) counts[i].count / n; // The probability of having this class
    *entropy -= p*log(p);
  }

  // Release the allocated memory
  free(counts);

  return 0;
}

int information_gini(double *gini, const int *y, const size_t n) {
  information_event_count *counts;
  size_t nevent;
  double p;
  size_t i;

  // Get the counts for each classes
  counts = information_count(&nevent, y, n);
  if(counts == NULL)
    return -1;

  // Compute the Gini impurity measure
  *gini = 0;
  for(i = 0; i < nevent; i++) {
    p = (double) counts[i].count / n; // The probability of having this class
    *gini += p*(1.-p);
  }

  // Release the allocated memory
  free(counts);

  return 0;
}
