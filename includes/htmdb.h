/**
 * @file htmdb.h
 *
 * Hierarchical Triangular Mesh database, header file.
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
#ifndef _HTMDB_H_
#define _HTMDB_H_

#include <htm.h>

/**
 * The HTM database structure
 */
typedef struct htmdb htmdb_t;

/**
 * Build an HTM database of the given level
 *
 * @param level The level of the HTM database
 *
 * @return A pointer to the newly allocated HTM database (must be freed with @htmbd_destroy)
 */
htmdb_t *htmdb_init(const unsigned int level);

/**
 * Add an entry into the HTM database
 *
 * @param db The HTM database in which this entry should be added
 * @param v The three dimensional cartesian coordinate of this star (must have a norm of 1)
 * @param data The data that are associated with this star
 *
 * @return 0 upon success; -1 otherwise
 */
int htmdb_add(htmdb_t *db, const double *v, void *data);

/**
 * Perform a query on the HTM database in the form of a constraint
 *
 * @param db The HTM database in which this entry should be added
 * @param p The three dimensional cartesian coordinates of the direction of the constraint (must have a norm of 1)
 * @param cosd The cosine of the search radius of the constraints (in [-1, 1])
 * @param process Callback function designed to process stars belonging to this constraint.
 * @param arg Arguments to pass to the process function
 *
 * The implementation of the process function should be
 * > int process(const double *v, void *data, void *arg);
 * where v and data are the fields passed in the initialization of the entry (in
 * the @htmdb_add function) and arg is the last parameter from @htmdb_constraint.
 * This function should return a zero value in order for the potential rest of
 * stars belonging to this constraint to processed or return a non-zero value if
 * we want to stop iteration (e.g. we are searching for a specific object and we
 * found it).
 *
 * @return 0 upon sucess; -1 upon error; 1 if the function stops because process returned a non-zero value
 */
int htmdb_constraint(const htmdb_t *db, const double *p, const double cosd, int (*process)(const double *, void *, void *), void *arg);

/*
TODO To implement
typedef struct htmdb_query htmdb_query_t;
void htmdb_write(const htmdb_t *db, FILE *out);
htmdb_t *htmdb_read(FILE *in);
htmdb_query_t *htmdb_query_new(const double *p, const double cosd);
htmdb_query_t *htmdb_query_and(htmdb_query_t *q1, htmdb_query_t *q2);
htmdb_query_t *htmdb_query_or(htmdb_query_t *q1, htmdb_query_t *q2);
int htmdb_query_exec(const htmdb_t *db, const htmdb_query_t *query, int (*process)(const double *, void *, void *), void *arg);
void htmdb_query_destroy(htmdb_query_t *query);
*/

/**
 * Destroy an HTM database
 *
 * @param db The HTM database to destroy
 * @param freedata The function that should be called in order to free data (see htmdb_add), NULL if no function is needed
 */
void htmdb_destroy(htmdb_t *db, void (*freedata)(void *));

#endif
