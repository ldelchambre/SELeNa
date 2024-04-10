#include <fload.h>
#include <fparse.h>
#include <pheap.h>
#include <stdlib.h>
#include <string.h>

typedef char *_fload_str_ptr;

struct fload_data {
  _fload_str_ptr **content;
  size_t nobs;
  size_t nattr;
  pheap_t *pheap;
};

typedef struct _fload_data_node {
  char **fields;
  size_t nfield;
  struct _fload_data_node *next;
} _fload_data_node_t;

typedef struct _fload_data_list {
  _fload_data_node_t *head;
  _fload_data_node_t *tail;
  size_t nobs;
  size_t nattr;
  pheap_t *pheap;
} _fload_data_list_t;

const char *_null_str = "";

int _fload_data_list_add(const unsigned int iline,
                         const char **fields,
                         const unsigned int nfield,
                         void *arg) {
  _fload_data_list_t *list = (_fload_data_list_t *) arg;
  _fload_data_node_t *node;
  unsigned int i, n;

  // Allocate node
  node = (_fload_data_node_t *) malloc(sizeof(_fload_data_node_t));
  if(node == NULL)
    return -1;

  // Allocate fields within node
  node->fields = (char **) calloc(nfield, sizeof(char *));
  if(node->fields == NULL) {
    free(node);
    return -1;
  }

  // Browse each field
  for(i = 0; i < nfield; i++) {
    // Compute the length of this field
    n = strlen(fields[i]);

    // If field[i] != ""
    if(n > 0) {
      // Allocate the field using the permanent heap
      node->fields[i] = (char *) pheap_calloc(list->pheap, n+1, sizeof(char));
      if(node->fields[i] == NULL) {
        free(node->fields);
        free(node);
        return -1;
      }

      // Copy field
      strcpy(node->fields[i], fields[i]);
    } else
      // If field[i] == "", set pointer to an empty string
      node->fields[i] = (char *) _null_str;
  }

  // Set the number of fields and pointer to the next node
  node->nfield = nfield;
  node->next = NULL;

  // Add this node to the list
  if(list->tail) {
    list->tail->next = node;
    list->tail = node;
  } else
    list->head = list->tail = node;

  // Update the nattr and nobs list attributes
  list->nobs++;
  if(list->nattr < nfield) list->nattr = nfield;

  return 0;
}

void _fload_data_list_copy(fload_data_t *res, const _fload_data_list_t *list) {
  _fload_str_ptr (*content)[res->nattr] = (_fload_str_ptr (*)[res->nattr]) res->content;
  _fload_data_node_t *node = list->head;
  size_t i, j;

  for(i = 0, node = list->head; i < res->nobs; i++, node = node->next) {
    for(j = 0; j < node->nfield; j++)
      content[i][j] = node->fields[j];
    for(; j < res->nattr; j++)
      content[i][j] = (char *) _null_str;
  }
}

void _fload_data_list_destroy(_fload_data_list_t *list) {
  _fload_data_node_t *node;

  while(list->head) {
    node = list->head;
    list->head = list->head->next;
    // Note we don't free node->fields[i]'s as it is the responsability of the permanent heap
    free(node->fields);
    free(node);
  }
}

fload_data_t *fload(FILE *in, const char sep) {
  pheap_t *pheap = NULL;
  fload_data_t *res = NULL;
  _fload_data_list_t list;

  // Initialize the permanent heap
  pheap = pheap_init(PHEAP_DEFAULT_CHUNCK_SIZE);
  if(pheap == NULL)
    return NULL;

  // Initialize the list structure
  list.pheap = pheap;
  list.head = list.tail = NULL;
  list.nobs = list.nattr = 0;


  // Parse the input file and store content in a list
  if(fparse(in, sep, _fload_data_list_add, &list)) {
    _fload_data_list_destroy(&list);
    pheap_destroy(pheap);
    return NULL;
  }

  // Allocate resulting structure
  res = pheap_malloc(pheap, sizeof(fload_data_t));
  if(res == NULL) {
    _fload_data_list_destroy(&list);
    pheap_destroy(pheap);
    return NULL;
  }

  // Initialize the resulting structure
  res->pheap = pheap;
  res->nattr = list.nattr;
  res->nobs = list.nobs;

  // Allocate content
  res->content = (_fload_str_ptr **) pheap_calloc(pheap, list.nobs*list.nattr, sizeof(_fload_str_ptr));
  if(res->content == NULL) {
    _fload_data_list_destroy(&list);
    pheap_destroy(pheap);
    return NULL;
  }

  // Copy the content of list into res
  _fload_data_list_copy(res, &list);

  // Destroy list
  _fload_data_list_destroy(&list);

  return res;
}

size_t fload_nobs(const fload_data_t *data) {
  if(data) return data->nobs;
  return -1;
}

size_t fload_nattr(const fload_data_t *data) {
  if(data) return data->nattr;
  return -1;
}

const char *fload_get(const fload_data_t *data, const size_t iobs, const size_t iattr) {
  if(data && iobs < data->nobs && iattr < data->nattr)
    return ((_fload_str_ptr (*)[data->nattr]) data->content)[iobs][iattr];
  return _null_str;
}

void fload_destroy(fload_data_t *data) {
  pheap_destroy(data->pheap);
}
