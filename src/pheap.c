/**
 * @file pheap.h
 *
 * Persistent heap, source file.
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
#include <pheap.h>

// ptrdiff_t
#include <stddef.h>

// memset
#include <string.h>

/**
 * Define _PHEAP_DEBUG in order to print debug messages
 */
// #define _PHEAP_DEBUG

/**
 * _DEBUG is a macros designed to print debug messages to standard error. If _PHEAP_DEBUG is not
 * defined, then this macros expand to nothing.
 */
#ifdef _PHEAP_DEBUG
// fprintf
#include <stdio.h>
#define _DEBUG(format,...) do{fprintf(stderr, "%s:%-6d : " format, __FILE__, __LINE__, ##__VA_ARGS__);fflush(stderr);}while(0)
#else
#define _DEBUG(format,...)
#endif

/**
 * A chunk a memory.
 *
 * Chunks of memory have the following structure:
 *   +---------------------------+ <- chunk
 *   |_pheap_chunk_t  *next      | 8 bytes
 *   +---------------------------+ <- allocated memory
 *   | ( char mem[0] )           | n bytes
 *   +---------------------------+
 */
typedef struct _pheap_chunk {
  struct _pheap_chunk *next; ///< Pointer to the next chunk of memory

#ifdef __GNUC__ /* GNU C allow for zero size array */
  char mem[0]; ///< Variable size array containing the allocated memory
#endif
} _pheap_chunk_t;

/**
 * Retrieve the pointer to the start of the allocated memory from a chunk.
 */
#ifdef __GNUC__
// In the case of GNU C, the allocated memory is simply stored at &chunk->mem
#define _pheap_chunk_mem(chunk) ((void *) &chunk->mem)
#else
// If this is not a GNU C source, then the allocated memory stands right after chunk->next
#define _pheap_chunk_mem(chunk) ((void *) (&chunk->next + 1))
#endif

struct pheap {
  _pheap_chunk_t *start; ///< Pointer to the head of the chunk list, or NULL if no chunk is already present
  _pheap_chunk_t *end;   ///< Pointer to the tail of the chunk list, or NULL if no chunk is already present
  void *imem;            ///< Pointer to the next memory area that should be returned, or &end->mem + mem_size if heap is currently full
  size_t mem_size;       ///< The chunk size minus the overhead size (i.e. the size of an empty chunk)
};

/**
 * Allocate and initialize a new chunk of memory
 *
 * @param size The size of the memory area to allocate (without overhead)
 *
 * @return A pointer to the newly allocated chunk, or NULL upon failure
 */
_pheap_chunk_t *_pheap_new_chunk(const size_t size) {
  // Allocate the new chunk
  _DEBUG("Allocating a new chunk of size %zd bytes + %zd bytes of overhead\n", size, sizeof(_pheap_chunk_t));
  _pheap_chunk_t *chunk = (_pheap_chunk_t *) malloc(sizeof(_pheap_chunk_t) + size);
  if(!chunk)
    return NULL;

  // Initialize the chunk
  chunk->next = NULL;

  return chunk;
}

/**
 * Compute the current free space from the heap (i.e. the memory space that remains in heap->end)
 *
 * @param heap The heap to use
 *
 * @return The currently available memory space in bytes
 */
size_t _pheap_free_space(const pheap_t *heap) {
  size_t ret;
  if(!heap->start)
    return 0;
  ret = heap->mem_size - (ptrdiff_t) (heap->imem - _pheap_chunk_mem(heap->end));
  return ret;
}

size_t pheap_min_chunck_size(void) {
  return sizeof(_pheap_chunk_t);
}

pheap_t *pheap_init(size_t chunck_size) {
  pheap_t *heap;

  // Do not allow a chunk size below pheap_min_chunck_size
  if(chunck_size < sizeof(_pheap_chunk_t))
    return NULL;

  // Allocate the persistent heap structure
  heap = (pheap_t *) malloc(sizeof(pheap_t));
  if(!heap)
    return NULL;

  // Initialize the persistent heap structure
  heap->start = heap->end = NULL;
  heap->imem = NULL;
  heap->mem_size = chunck_size - sizeof(_pheap_chunk_t);

  _DEBUG("Initializing heap %p (mem_size=%zd)\n", heap, heap->mem_size);

  return heap;
}

void *pheap_malloc(pheap_t *heap, size_t size) {
  void *ret;

  // If the heap was not correctly initialized or that the asked size is zero
  // then return NULL
  if(!heap || size == 0)
    return NULL;

  _DEBUG("Allocating %zd byte%s from heap %p (free space: %zd)\n", size, (size > 1) ? "s" : "", heap, _pheap_free_space(heap));

  // If not enough free space is available, allocate a new chunk of memory
  if(_pheap_free_space(heap) < size) {
    _pheap_chunk_t *chunk;

    // If the asked allocation exceed the chunk size, then allocate a chunk of
    // higher size and add it to the front of the list of chunks (i.e. this chunk
    // will be immediately full, so don't add it to the end while resetting heap->imem)
    if(heap->mem_size < size) {
      // Allocate the chunk
      chunk = _pheap_new_chunk(size);
      if(!chunk) return NULL;
      // Get the newly allocated memory location
      ret = _pheap_chunk_mem(chunk);
      // Add this chunk to the front of the list
      chunk->next = heap->start;
      heap->start = chunk;
      // If this chunk is the first that is allocated on this heap
      if(!heap->end) {
        // Also update the end and imem attributes
        heap->end = chunk;
        heap->imem = ret + heap->mem_size; // Mark this chunk as full
      }
      // Return the newly allocated memory area
      return ret;
    } else { // Otherwise, allocate a new chunk of size 'size'
      chunk = _pheap_new_chunk(heap->mem_size);
      if(!chunk) return NULL;
      // If the heap is currently empty, chunk will become the sole chunk in the heap
      if(!heap->start)
        heap->start = heap->end = chunk;
      else { // Otherwise, add it to the end of the list of chunks
             // Note that some memory might not have been allocated from heap->end
             // (i.e. if _pheap_free_space(heap) > 0). This memory will not be
             // allocated again and will be wasted.
             // /!\ However, we must not use realloc function so as to dim the
             //     available memory size at it is no guarantee that the returned
             //     pointer will stand at the same location as the initial pointer.
        heap->end->next = chunk;
        heap->end = chunk;
      }
      // Set heap->imem to point to the newly allocated memory area
      heap->imem = _pheap_chunk_mem(chunk);
    }
  }

  // The memory area to return stand at heap->imem
  ret = heap->imem;

  // Increase heap->imem because it is now allocated
  heap->imem += size;

  // Return pointer
  return ret;
}

void *pheap_calloc(pheap_t *heap, size_t nmemb, size_t size) {
  const size_t n = nmemb * size;
  void *ret = pheap_malloc(heap, n);
  if(!ret)
    return NULL;
  return memset(ret, 0, n);
}

void pheap_destroy(pheap_t *heap){
  _pheap_chunk_t *next;

  _DEBUG("Destroying heap %p\n", heap);

  if(!heap)
    return;

  while(heap->start) {
    next = heap->start->next;
    free(heap->start);
    heap->start = next;
  }

  free(heap);
}
