/**
 * @file pheap.h
 *
 * Persistent heap, header file.
 *
 *  A persistent heap is a memory structure that allow to allocate
 * memory regions if all those regions are freed at once. This
 * implementation allows to efficiently deal with those memory
 * allocations without having to rely on numerous time-consuming
 * malloc/calloc. It is particularly well adapted to multi-threaded
 * code where those malloc's constitute bottlenecks for the entire
 * program.
 *   Rather than allocating lots of small memory regions through
 * malloc/calloc, this implementation instead allocate a set of
 * typically large memory regions, that are memory chunks, where
 * each of these chunk is partitioned into smaller regions
 * corresponding to those returned by pheap_malloc.
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
#ifndef _PHEAP_H_
#define _PHEAP_H_

// size_t
#include <stdlib.h>

/**
 * The default chunk size (1048576 bytes = 1MB).
 */
#define PHEAP_DEFAULT_CHUNCK_SIZE ( 1U << 20 )

/**
 * The persistent heap structure
 */
typedef struct pheap pheap_t;

/**
 * The minimal chunk size. The lower bound on the size of the chunk is
 * imposed by the amount of overhead that is necessary in order to manage
 * chunks and heap.
 *
 * @return The minimal chunk size in bytes.
 */
size_t pheap_min_chunck_size(void);

/**
 * Initialize a persistent heap structure coming along with the specified
 * chunk size.
 *
 * @param chunk_size The chunk size to use
 *
 * @return A pointer to the newly allocated heap structure upon success, NULL otherwise
 */
pheap_t *pheap_init(size_t chunck_size);

/**
 * Allocate size bytes from the heap.
 *
 * @param heap The heap through which memory allocation is provided.
 * @param size The size of the memory to be allocated
 *
 * @return A pointer to the allocated memory, or NULL upon failure
 */
void *pheap_malloc(pheap_t *heap, size_t size);

/**
 * Allocate an array of nmemb, each of size bytes from the heap.
 *
 * Allocated memory will be initialized to zero.
 *
 * @param heap  The heap through which memory allocation is provided.
 * @param nmemb The number of elements to allocate
 * @param size  The size of each elements in bytes
 *
 * @return A pointer to the allocated memory, or NULL upon failure
 */
void *pheap_calloc(pheap_t *heap, size_t nmemb, size_t size);

/**
 * Destroy a persistent map and free all memory that was allocated through it.
 *
 * @param heap The heap to destroy
 */
void pheap_destroy(pheap_t *heap);

#endif /* _PHEAP_H_ */
