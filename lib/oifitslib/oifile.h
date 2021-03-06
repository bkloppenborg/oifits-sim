/* $Id: oifile.h,v 1.4 2007-11-01 12:56:17 jsy1001 Exp $ */

/**
 * @file oifile.h
 * @ingroup oifile
 *
 * Data structure definitions and function prototypes for file-level
 * operations on OIFITS data.
 *
 * Copyright (C) 2007 John Young
 *
 *
 * This file is part of OIFITSlib.
 *
 * OIFITSlib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * OIFITSlib is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with OIFITSlib.  If not, see
 * http://www.gnu.org/licenses/
 */

/**
 * @defgroup oifile  File-level OIFITS I/O
 *
 * This module provides routines for reading and writing entire OIFITS
 * files. Use read_oi_fits() to read a file. free_oi_fits() will free
 * the storage allocated by read_oi_fits(). Use init_oi_fits() to
 * initialise an empty dataset, then after filling in the contents
 * call write_oi_fits() to write the dataset to a file.
 *
 * Functions to format and display strings summarising the file
 * contents are provided: format_oi_fits_summary() and
 * print_oi_fits_summary()

 * A set of functions oi_fits_lookup_*()
 * (e.g. oi_fits_lookup_array())) are also provided, to facilitate
 * following cross-references between OI_FITS tables.
 *
 * @{
 */

#ifndef OIFILE_H
#define OIFILE_H

#include <glib.h>
#include "exchange.h"


/*
 * Macros
 */

#define MEMDUP(dest, src, size) \
  { (dest) = malloc((size)); memcpy((dest), (src), (size)); } 


/*
 * Data structures
 */

/** Data for OIFITS file */
typedef struct {
  int numArray;              /**< Length of arrayList */
  int numWavelength;         /**< Length of wavelengthList */
  int numVis;                /**< Length of visList */
  int numVis2;               /**< Length of vis2List */
  int numT3;                 /**< Length of t3List */
  int numT4;                 /**< Length of t4List */
  oi_target targets;         /**< oi_target struct */
  GList *arrayList;          /**< Linked list of oi_array structs */
  GList *wavelengthList;     /**< Linked list of oi_wavelength structs */
  GList *visList;            /**< Linked list of oi_vis structs */
  GList *vis2List;           /**< Linked list of oi_vis2 structs */
  GList *t3List;             /**< Linked list of oi_t3 structs */
  GList *t4List;             /**< Linked list of oi_t4 structs */
  GHashTable *arrayHash;     /**< Hash table of oi_array, indexed by ARRNAME */
  GHashTable *wavelengthHash; /**< Hash table of oi_wavelength, 
				   indexed by INSNAME */
} oi_fits;


/*
 * Function prototypes, for functions from oifile.c
 */
void init_oi_fits(oi_fits *);
STATUS write_oi_fits(const char *, oi_fits, STATUS *);
STATUS read_oi_fits(const char *, oi_fits *, STATUS *);
void free_oi_fits(oi_fits *);
oi_array *oi_fits_lookup_array(const oi_fits *, const char *);
element *oi_fits_lookup_element(const oi_fits *, const char *, int);
oi_wavelength *oi_fits_lookup_wavelength(const oi_fits *, const char *);
target *oi_fits_lookup_target(const oi_fits *, int);
const char *format_oi_fits_summary(const oi_fits *);
void print_oi_fits_summary(const oi_fits *);
oi_target *dup_oi_target(const oi_target *);
oi_array *dup_oi_array(const oi_array *);
oi_wavelength *dup_oi_wavelength(const oi_wavelength *);
oi_vis *dup_oi_vis(const oi_vis *);
oi_vis2 *dup_oi_vis2(const oi_vis2 *);
oi_t3 *dup_oi_t3(const oi_t3 *);

#endif /* #ifndef OIFILE_H */

/** @} */
