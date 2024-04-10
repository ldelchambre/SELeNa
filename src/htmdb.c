/**
 * @file htmdb.c
 * 
 * Hierarchical Triangular Mesh database, source file.
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
#include <htmdb.h>
#include <pheap.h>
#include <stdlib.h>

/**
 * The number of HTMDB entry to allocate in each chunk of the pesrsisten heap
 */
#define HTMDB_PHEAP_NENTRY (1U << 20)

/**
 * Linked list of the HTM db entries
 */
typedef struct _htmdb_entry {
  double v[3]; ///< The position of this star
  void *data; ///< The data that are associated with this star
  struct _htmdb_entry *next; ///< Pointer to the next entry
} _htmdb_entry_t;

/**
 * Structure of the HTM db
 */
struct htmdb {
  pheap_t *heap; ///< The heap where to allocate the _htmdb_entry
  _htmdb_entry_t **entries; ///< Entries of the HTM db (one per HTM id)
  unsigned int level; ///< The level of this HTM db
  htmid mask; ///< The mask to apply to the HTM id in order to get the index within @entries
};

htmdb_t *htmdb_init(const unsigned int level) {
  htmdb_t *db;
  pheap_t *heap;

  // Do not accept to overflow HTM level
  if(level > HTM_LEVEL_MAX)
    return NULL;

  // Initialize the persisten heap
  heap = pheap_init(HTMDB_PHEAP_NENTRY * sizeof(_htmdb_entry_t));
  if(!heap)
    return NULL;

  // Allocate the HTM database structure
  db = pheap_malloc(heap, sizeof(htmdb_t));
  if(!db)
    return NULL;
    
  // Set heap
  db->heap = heap;

  // Allocate the HTM entries
  db->entries = pheap_calloc(db->heap, 0x08ULL << ( 2 * level ), sizeof(_htmdb_entry_t *));
  if(!db->entries) {
    pheap_destroy(heap);
    return NULL;
  }

  // Assign level
  db->level = level;

  // Compute the mask that is necessary in order to query database
  db->mask = ( 0x08ULL << ( 2 * level ) ) - 1ULL;

  return db;
}

int htmdb_add(htmdb_t *db, const double *v, void *data) {
  // Get the HTM id for this entry
  const htmid id = htm_get(v, db->level);

  // Allocate this entry
  _htmdb_entry_t *e = pheap_malloc(db->heap, sizeof(_htmdb_entry_t));

  // If allocation fails
  if(!e)
    return -1;

  // Assign fields to this entry
  htm_copy_vertex(e->v,v);
  e->data = data;

  // Add this entry to the list
  e->next = db->entries[id & db->mask];
  db->entries[id & db->mask] = e;

  return 0;
}

int htmdb_constraint(const htmdb_t *db, const double *p, const double cosd, int (*process)(const double *, void *, void *), void *arg) {
  // Get the list of HTM identifiers for this constraint
  htm_list_t *l = htm_constraint(p,cosd,db->level);
  htm_list_t *il;

  // If the list of HTM identifiers cannot be retrieved
  if(!l)
    return -1;

  // Browse each HTM identifier of the constraint
  for(il = l; il; il = il->next) {
    // Get the minimal HTM identifier corresponding to the level of this database
    htmid id = il->id << (2 * (db->level - htm_level(il->id)));

    // Browse each HTM identifier of the database
    // corresponding to the HTM identifer of the constraint
    for(; htm_cmp(il->id, id) == 0; id++) {
      // Get the root entry
      const _htmdb_entry_t *e = db->entries[id & db->mask];
      // Browse each entry
      for(; e; e = e->next) {
        // If we have to process this entry
        if(e->v[0]*p[0]+e->v[1]*p[1]+e->v[2]*p[2] >= cosd) {
          // Process it
          if(process(e->v, e->data, arg)) {
            // If process return with a non-zero value, then we have to stop our
            // iteration
            htm_list_destroy(l);
            return 1;
          }
        }
      }
    }
  }

  // Destroy the list of HTM identifiers
  htm_list_destroy(l);
  return 0;
}

void htmdb_destroy(htmdb_t *db, void (*freedata)(void *)) {
  pheap_t *heap;
  
  // If db is NULL, then nothing has to be done
  if(!db)
    return;
    
  // Get heap
  heap = db->heap;

  // If a destructor is provided and db has some entries
  if(freedata && db->entries) {
    const unsigned long long n = db->mask + 1ULL;
    unsigned long long i;

    // Browse each entry of the database
    for(i = 0; i < n; i++) {
      // Get this entry
      _htmdb_entry_t *e = db->entries[i];

      // Browse the list of entries
      while(e) {
        // Call destructor on data
        freedata(e->data);
        // Go to next entry
        e = e->next;
      }
    }
  }

  // Delete the database
  pheap_destroy(heap);
}
