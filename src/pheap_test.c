/**
 * @file pheap.h
 *
 * Persistent heap, test file.
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
#include <stdio.h>
#include <stdlib.h>

#define NTEST           65536
#define NMALLOC         1024
#define CHUNK_MAX_SIZE  65536
#define MALLOC_MAX_SIZE 256

/**
 * This test program should be run in combination with the valgrind memory debugger
 */
int main(void) {
  const size_t min_size = pheap_min_chunck_size();
  pheap_t *heap;
  size_t chunk_size, n;
  char *mem;
  unsigned int i, j, k;

  for(i = 0; i < NTEST; i++) {
    // Initialize the random number generator
    srand(i+1);

    // Get the chunk size we will use
    chunk_size = min_size + rand() % ( CHUNK_MAX_SIZE + 1 - min_size );

    // Create the persistent heap
    heap = pheap_init(chunk_size);
    if(!heap) {
      fprintf(stderr, "Unable to create a persistent heap of size %zd\n", chunk_size);
      exit(EXIT_FAILURE);
    }

    // Do a lot of memory allocation/access
    for(j = 0; j < NMALLOC; j++) {
      n = 1 + rand() % MALLOC_MAX_SIZE;
      mem = (char *) pheap_malloc(heap, n);
      if(!mem) {
        fprintf(stderr, "Unable to allocate %zd bytes from heap (chunk_size=%zd, seed=%u, imalloc=%u)\n", n, chunk_size, i+1, j);
        pheap_destroy(heap);
        exit(EXIT_FAILURE);
      }
      for(k = 0; k < n; k++)
        mem[k] = k;
      // There is no pheap_free(mem), because this is persistent heap
    }

    // Destroy the heap
    pheap_destroy(heap);
  }

  exit(EXIT_SUCCESS);
}
