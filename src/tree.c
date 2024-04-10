#include <tree.h>

#include <errno.h>
#include <pheap.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

/**
 * Comment in order to disable the debug mode
 */
// #define _TREE_DEBUG_MODE

/**
 * In debug mode, _TREE_DEBUG will print a debug message
 */
#ifdef _TREE_DEBUG_MODE
#define _TREE_DEBUG(format,...) do{fprintf(stdout, format, ##__VA_ARGS__);fflush(stdout);}while(0)
#else
#define _TREE_DEBUG(format,...)
#endif

/**
 * Common error handling
 */
#define _TREE_ERROR(format, ...)\
do {\
  int errsv = errno; \
  fprintf(stderr, "Tree - "format": ", ##__VA_ARGS__); \
  if(errsv) \
    fprintf(stderr, "%s", strerror(errsv)); \
  fprintf(stderr, "\n");\
  fflush(stderr);\
} while(0)

/**
 * The tree file magic number designed to uniquely identify a file as containing
 * a tree we created.
 */
#define _TREE_MAGIC_NUMBER (float) 0.385339

/**
 * The default chunck size to use in order to allocate heap memory from the persistent heap
 */
#define _TREE_PHEAP_CHUNCK_SIZE PHEAP_DEFAULT_CHUNCK_SIZE

/**
 * Structure containing the representation of a tree node
 */
typedef struct _tree_node {
  struct _tree_node *left;  //< The left node
  struct _tree_node *right; //< The right node
  // We use a union because the data field will never be filled up at the same time as iattr and cutval
  // i.e. if left == NULL or right == NULL then we are in a leaf node and hence, data is set
  // otherwise, we are in an inner node such that iattr and cutval are set
  union {
    struct {
      size_t iattr;  // If in a inner node, index of the attribute considered for splitting
      double cutval; // If in a inner node, value of the attribute considered for splitting
    };
    struct {
      void *data;       // If in a leaf node, value associated with this node
      size_t data_size; // If in a leaf node, size of the data in byte(s)
    };
  };
} _tree_node_t;

/**
 * Structure containing the representation of a tree
 */
struct tree {
  pheap_t *heap;                           // The persistent heap we used
  _tree_node_t *root;                      // The root node
};

/**
 * Structure containing the tree building parameters that are common to all threads
 */
struct _tree_build_common_arg {
  tree_t *tree;     // The tree that is currently being built
  const double *X;  // The input data set
  const void *y;    // The class to predict
  size_t nattr;     // The number of attribute within the input data set
  threadpool_t *pool; // The thread pool to use in order to build the tree
  pthread_mutex_t lock; // Lock for synchronizing builds (e.g. heap allocation)
  const tree_algorithm *algorithm; // The tree building algorithm
  struct {
    size_t ncurr; // The number of instances we don't further have to process
    size_t ntot;  // The total number of instances we have to process
    size_t nnode;
  } progress; // Variables maintaining progression
  int status;       // The current status of the tree building (0 if no error already occurs)
};

struct tree_build_job_t {
  struct _tree_build_common_arg *common; // Thread's common parameters
  _tree_node_t *tree;                    // The node that should be developed by the thread
  size_t *inst;                          // The index of the instances composing the subset that should be treated by the thread
  size_t ninst;                          // The number of instances within this subset
};

tree_t *_tree_init(void) {
  // Allocate the tree
  tree_t *tree = (tree_t *) malloc(sizeof(tree_t));
  if(tree == NULL) {
    _TREE_ERROR("Error while allocating the tree");
    return NULL;
  }

  // Initialize the persistent heap
  tree->heap = pheap_init(_TREE_PHEAP_CHUNCK_SIZE);
  if(tree->heap == NULL) {
    _TREE_ERROR("Error while initializing the persistent heap");
    free(tree);
    return NULL;
  }

  // Allocate the root node (through persistent heap)
  tree->root = (_tree_node_t *) pheap_malloc(tree->heap, sizeof(_tree_node_t));
  if(tree->root == NULL) {
    _TREE_ERROR("Error while allocating the root node");
    pheap_destroy(tree->heap);
    free(tree);
    return NULL;
  }

  // Initialize the root node
  tree->root->data = tree->root->left = tree->root->right = NULL;

  return tree;
}

/**
 * Develop a node (or create a leaf node) using the provided thread's parameter.
 *
 * @param[in,out] arg The threads parameters used to develop this node
 */
void _tree_split_node(void *arg) {
  tree_build_job *job = (tree_build_job *) arg; // The job that this thread have to do
  _tree_node_t *node = job->tree; // The node to split

  // If we are in an error state, do nothing
  if(job->common->status)
    return;

  // Find a split for the given set of instances
  node->iattr = job->common->algorithm->getsplit(&node->cutval, job->common->X,
                                                 job->common->y, job->inst,
                                                 job->ninst, job->common->nattr,
                                                 job->common->algorithm->arg);

  // If an error occurs during the split finding, exit
  if(node->iattr == -1UL) {
    _TREE_ERROR("Error while finding a node splitting criterion");
    job->common->status = -1;
    return;
  }

  // If we have to build a leaf node
  if(node->iattr == job->common->nattr) {
    _TREE_DEBUG("Creating leaf node with %zu instances\n", job->ninst);

    // Get the data associated with the current subset
    node->data = NULL;
    node->data_size = 0;
    {
      void *tmp = job->common->algorithm->getdata(&node->data_size,
                                                 job->common->y, job->inst, job->ninst,
                                                 job->common->algorithm->arg);
      if(tmp) {
        // Copy the data into a memory location allocated using the persistent heap
        pthread_mutex_lock(&job->common->lock);
        node->data = pheap_malloc(job->common->tree->heap, node->data_size);
        pthread_mutex_unlock(&job->common->lock);
        if(node->data)
          memcpy(node->data, tmp, node->data_size);
        free(tmp);
      }
    }
    if(node->data == NULL) {
      _TREE_ERROR("Error while allocating a leaf node prediction");
      job->common->status = -1;
      return;
    }
    if(node->data_size == 0) {
      _TREE_ERROR("Error: leaf node prediction has a size of zero");
      job->common->status = -1;
      return;
    }

    // Update progression
    pthread_mutex_lock(&job->common->lock);
    job->common->progress.ncurr += job->ninst;
    pthread_mutex_unlock(&job->common->lock);

  } else { // Otherwise, we have to recursively split the tree
    const double (*X)[job->common->nattr] = (const double (*)[job->common->nattr]) job->common->X; // The input dataset
    tree_build_job *ljob; // The left split job
    tree_build_job *rjob; // The right split job
    size_t from, to, tmp; // Indices within the instances array

    // Split the instances according to the split criterion we found
    from = 0;
    to = job->ninst;
    while(from < to) {
      while(from < to && X[job->inst[from]][node->iattr] < node->cutval)
        from++;
      while(from < to && node->cutval <= X[job->inst[to-1]][node->iattr])
        to--;
      // Swap job->inst[from] and job->inst[to-1]
      if(from < to) {
        tmp = job->inst[from];
        job->inst[from] = job->inst[to-1];
        job->inst[to-1] = tmp;
      }
    }

    // If we found an invalid split, then exit
    if(to == 0 || to == job->ninst) {
      _TREE_ERROR("The provided split lead to a node having 0 instances");
      job->common->status = -1;
      return;
    }

    // Allocate the left and right node
    pthread_mutex_lock(&job->common->lock);
    job->common->progress.nnode++;
    node->left = (_tree_node_t *) pheap_malloc(job->common->tree->heap, sizeof(_tree_node_t));
    node->right = (_tree_node_t *) pheap_malloc(job->common->tree->heap, sizeof(_tree_node_t));
    pthread_mutex_unlock(&job->common->lock);
    if(node->left == NULL || node->right == NULL) {
      _TREE_ERROR("Error while allocating child nodes");
      job->common->status = -1;
      node->left = node->right = node->data = NULL;
      return;
    }

    // Initialize the left and right node
    node->left->data = node->left->left = node->left->right = NULL;
    node->right->data = node->right->left = node->right->right = NULL;

    _TREE_DEBUG("Node created (iattr=%zu, cutval=%f, nleft=%zu, nright=%zu)\n", node->iattr, node->cutval, to, job->ninst-to);

    // Allocate the left and right thread arguments
    ljob = (tree_build_job *) malloc(sizeof(tree_build_job));
    rjob = (tree_build_job *) malloc(sizeof(tree_build_job));
    if(ljob == NULL || rjob == NULL) {
      _TREE_ERROR("Error while allocating the child thread arguments");
      job->common->status = -1;
      free(ljob);free(rjob);
      return;
    }

    // Set the left thread arguments
    ljob->common = job->common;
    ljob->tree = node->left;
    ljob->inst = job->inst;
    ljob->ninst = to;

    // Set the right thread arguments
    rjob->common = job->common;
    rjob->tree = node->right;
    rjob->inst = &job->inst[to];
    rjob->ninst = job->ninst-to;

    // Add jobs to the pool
    if(threadpool_addjob(job->common->pool, _tree_split_node, ljob, free)
       || threadpool_addjob(job->common->pool, _tree_split_node, rjob, free)) {
      _TREE_ERROR("Tree building error: Error while creating child jobs");
      job->common->status = -1;
      // We don't free ljob and rjob here because it is the responsability of the pool now
      return;
    }
  }
}

tree_build_job *tree_build_nowait(const double *X, const void *y, const size_t ninst, const size_t nattr, threadpool_t *pool, const tree_algorithm *algorithm) {
  tree_build_job *job;
  struct _tree_build_common_arg *common;
  size_t *inst;

  // Check that an algorithm was provided
  if(algorithm == NULL) {
    _TREE_ERROR("No algorithm provided");
    return NULL;
  }

  // Check that the algorithm getsplit and getdata functions are presents
  if(algorithm->getsplit == NULL || algorithm->getdata == NULL) {
    _TREE_ERROR("Algorithm should have valid getsplit and getdata functions");
    return NULL;
  }

  // Initialize the input algorithm
  if(algorithm->init && algorithm->init(X, y, ninst, nattr, algorithm->arg))
    return NULL;

  // Allocate the job; the common arguments and the array of instances
  job = (tree_build_job *) malloc(sizeof(tree_build_job));
  common = (struct _tree_build_common_arg *) malloc(sizeof(struct _tree_build_common_arg));
  inst = (size_t *) calloc(ninst, sizeof(size_t));
  if(job == NULL || inst == NULL || common == NULL) {
    _TREE_ERROR("Can not allocate the tree building job");
    free(job);free(common);free(inst);
    return NULL;
  }

  // Initialize the parameters that are common to all threads
  common->tree = _tree_init();
  if(common->tree == NULL) {
    free(job);free(common);free(inst);
    return NULL;
  }

  // Initialize the mutex
  if(pthread_mutex_init(&common->lock, NULL)) {
    _TREE_ERROR("Error while initializing the tree building mutex");
    tree_free(common->tree);free(job);free(common);free(inst);
    return NULL;
  }
  common->X = X;
  common->y = y;
  common->nattr = nattr;
  common->pool = pool;
  common->algorithm = algorithm;
  common->progress.ncurr = 0;
  common->progress.ntot  = ninst;
  common->progress.nnode = 0;
  common->status = 0;


  // Initialize the array of instances
  for(size_t i = 0; i < ninst; i++)
    inst[i] = i;

  // Initialize the parameters that are specific to the initial thread
  job->common = common;
  job->tree = common->tree->root;
  job->inst = inst;
  job->ninst = ninst;

  // Add a job to the pool in order to build the tree
  if(threadpool_addjob(pool, _tree_split_node, job, NULL)) {
    _TREE_ERROR("Can not launch the tree building process");
    tree_free(common->tree);pthread_mutex_destroy(&common->lock);
    free(job);free(common);free(inst);
    return NULL;
  }

  return job;
}

int tree_build_progess(size_t *n, size_t *ntot, size_t *nnode, tree_build_job *job) {
  if(job == NULL)
    return -1;
  pthread_mutex_lock(&job->common->lock);
  if(n != NULL)
    *n = job->common->progress.ncurr;
  if(ntot != NULL)
    *ntot = job->common->progress.ntot;
  if(nnode != NULL)
    *nnode = job->common->progress.nnode;
  pthread_mutex_unlock(&job->common->lock);
  return 0;
}

tree_t *tree_build_end(tree_build_job *job) {
  tree_t *tree;

  if(job == NULL)
    return NULL;
  tree = job->common->tree;

  // If the tree building failed
  if(job->common->status) {
    tree_free(tree);free(job->common);free(job->inst);free(job);
    return NULL;
  }

  // Free the ressources that were allocated for thr job
  pthread_mutex_destroy(&job->common->lock);
  free(job->common);free(job->inst);free(job);

  // Return the built tree
  return tree;
}

tree_t *tree_build(const double *X, const void *y, const size_t ninst, const size_t nattr, threadpool_t *pool, const tree_algorithm *algorithm) {
  // Initialize the tree building job
  tree_build_job *job = tree_build_nowait(X,y,ninst,nattr,pool,algorithm);

  // If the tree building already fails, return NULL
  if(job == NULL)
    return NULL;

  // Wait for the job to be done
  if(threadpool_wait(pool)) {
    _TREE_ERROR("Can not wait for the whole tree building jobs to end");
    tree_free(tree_build_end(job));
    return NULL;
  }

  // Get the built tree
  return tree_build_end(job);
}

void *_tree_get(size_t *data_size, const _tree_node_t *node, const double *x) {
  // If we are in a leaf node
  if(node->left == NULL /* && tree->right == NULL*/) {
    if(data_size) *data_size = node->data_size;
    return node->data; // Return the prediction
  }

  // If the instance attribute is lower than the node's attribute value
  if (x[node->iattr] < node->cutval)
    return _tree_get(data_size,node->left, x); // Browse the left part of the tree
  return _tree_get(data_size,node->right, x); // Otherwise browse the right part of the tree
}

void *tree_get(size_t *data_size, const tree_t *tree, const double *x) {
  if(tree == NULL)
    return NULL;
  return _tree_get(data_size, tree->root, x);
}

int _tree_dump(FILE *fout, const _tree_node_t *node) {
  // If we are in a leaf node
  if(node->left == NULL /* && tree->right == NULL*/) {
    // Write this node's data
    if(fwrite(&node->data_size,sizeof(size_t),1,fout) == 0
       || fwrite(node->data,node->data_size,1,fout) == 0) {
      _TREE_ERROR("Unable to write tree data into file");
      return -1;
    }
  } else {// If we are in an inner node
    const size_t n = 0;

    // Write this node
    if(fwrite(&n,sizeof(size_t),1,fout) == 0
       || fwrite(&node->iattr,sizeof(size_t),1,fout) == 0
       || fwrite(&node->cutval,sizeof(double),1,fout) == 0) {
      _TREE_ERROR("Unable to write tree inner node into file");
      return -1;
    }

    // Write the left and right trees
    if(_tree_dump(fout, node->left) || _tree_dump(fout, node->right))
      return -1;
  }

  return 0;
}

int tree_dump(FILE *fout, const tree_t *tree) {
  const float magicn = _TREE_MAGIC_NUMBER;
  const int version = TREE_SERIALIZATION_VERSION;

  if(tree == NULL)
    return -1;

  // Write the header data structure
  if(fwrite(&magicn, sizeof(float), 1, fout) == 0 // Write the magic number
     || fwrite(&version, sizeof(int), 1, fout) == 0){ // Write the Tree serialization version
    _TREE_ERROR("Unable to write binary tree file header");
    return -1;
  }

  // Write the provided tree into the file
  if(_tree_dump(fout,tree->root)) {
    _TREE_ERROR("Unable to write tree into the output file");
    return -1;
  }

  return 0;
}

int _tree_load(_tree_node_t *node, pheap_t *heap, FILE *fin) {
  // Read the node size (set to 0 if inner node)
  if(fread(&node->data_size, sizeof(size_t), 1, fin) == 0)
    return -1;

  if(node->data_size != 0) { // If we are in a leaf node
    // Say that we are in a leaf node
    node->left = node->right = NULL;

    // Allocate the data
    node->data = pheap_malloc(heap,node->data_size);
    if(node->data == NULL)
      return -1;

    // Read the data from file
    if(fread(node->data, node->data_size, 1, fin) == 0) {
      free(node->data);
      return -1;
    }
  } else { // If we are in an inner node

    // Allocate the left and right part of the tree
    node->left  = (_tree_node_t *) pheap_malloc(heap, sizeof(_tree_node_t));
    node->right = (_tree_node_t *) pheap_malloc(heap, sizeof(_tree_node_t));
    if(node->left == NULL || node->right == NULL) {
      node->left = node->right = node->data = NULL;
      return -1;
    }

    // Read the inner node's informations
    if(fread(&node->iattr,sizeof(size_t),1,fin) == 0
       || fread(&node->cutval,sizeof(double),1,fin) == 0) {
      node->left = node->right = node->data = NULL;
      return -1;
    }

    // Process the left part and right part of the tree
    if(_tree_load(node->left, heap, fin) || _tree_load(node->right, heap, fin)) {
      node->left = node->right = node->data = NULL;
      return -1;
    }
  }

  return 0;
}

tree_t *tree_load(FILE *fin) {
  tree_t *tree;
  float magicn;
  int version;

  // Read the data header
  if(fread(&magicn,sizeof(float),1,fin) == 0 ||
     fread(&version, sizeof(int), 1, fin) == 0) {
    _TREE_ERROR("Unable to read data header");
    return NULL;
  }

  // Check that the magic number match _tree_serialize_magicn
  if(magicn != _TREE_MAGIC_NUMBER) {
    _TREE_ERROR("File is not a binary tree file");
    return NULL;
  }

  // Check that this version is up to date
  if(version < TREE_SERIALIZATION_VERSION) {
    _TREE_ERROR("The provided file was written with a sofware version that is more recent that the current one");
    return NULL;
  }

  // Initialize the tree
  tree = _tree_init();

  // Load the tree from the file
  if(_tree_load(tree->root,tree->heap,fin)) {
    _TREE_ERROR("Unable to read tree");
    return NULL;
  }

  return tree;
}

void tree_free(tree_t *tree) {
  if(tree == NULL)
    return;
  pheap_destroy(tree->heap);
  free(tree);
}
