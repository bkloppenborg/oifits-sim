/* $Id: demo.c,v 1.3 2008-04-24 16:13:28 jsy1001 Exp $ */

/**
 * @file demo.c
 *
 * Example program - uses oitable API to write and read a OIFITS file.
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "exchange.h"


/**
 * Read data from a user-specified text file and write out in OIFITS format.
 */
void demo_write(void)
{
  oi_array array;
  oi_target targets;
  oi_wavelength wave;
  oi_vis vis;
  oi_vis2 vis2;
  oi_t3 t3;
  char filename[FLEN_FILENAME];
  FILE *fp;
  int i, irec, iwave, itarg, status;
  fitsfile *fptr;

  /* get input filename */
  printf("Enter input TEXT filename: ");
  fgets(filename, FLEN_FILENAME, stdin);
  filename[strlen(filename)-1] = '\0'; /* zap newline */

  fp = fopen(filename, "r");

  /* Read info for OI_ARRAY table */
  fscanf(fp, "OI_ARRAY arrname %s ", array.arrname);
  fscanf(fp, "frame %s ", array.frame);
  fscanf(fp, "arrayx %lf arrayy %lf arrayz %lf ", &array.arrayx,
	 &array.arrayy, &array.arrayz);
  fscanf(fp, "nelement %d ", &array.nelement);
  array.elem = malloc(array.nelement*sizeof(element));
  for (i=0; i<array.nelement; i++) {
    
    fscanf(fp, "tel_name %s sta_name %s ", array.elem[i].tel_name,
	   array.elem[i].sta_name);
    fscanf(fp, "staxyz %lf %lf %lf diameter %f ", &array.elem[i].staxyz[0],
	   &array.elem[i].staxyz[1], &array.elem[i].staxyz[2],
	   &array.elem[i].diameter);
    array.elem[i].sta_index = i+1;
  }
  array.revision = 1;

  /* Read info for OI_TARGET table */
  fscanf(fp, "OI_TARGET ntarget %d ", &targets.ntarget);
  targets.targ = malloc(targets.ntarget*sizeof(target));
  for (itarg=0; itarg<targets.ntarget; itarg++) {
    fscanf(fp, "target_id %d ", &targets.targ[itarg].target_id);
    fscanf(fp, "target %s ", targets.targ[itarg].target);
    fscanf(fp, "raep0 %lf ", &targets.targ[itarg].raep0);
    fscanf(fp, "decep0 %lf ", &targets.targ[itarg].decep0);
    fscanf(fp, "equinox %f ", &targets.targ[itarg].equinox);
    fscanf(fp, "ra_err %lf ", &targets.targ[itarg].ra_err);
    fscanf(fp, "dec_err %lf ", &targets.targ[itarg].dec_err);
    fscanf(fp, "sysvel %lf ", &targets.targ[itarg].sysvel);
    fscanf(fp, "veltyp %s ", targets.targ[itarg].veltyp);
    fscanf(fp, "veldef %s ", targets.targ[itarg].veldef);
    fscanf(fp, "pmra %lf ", &targets.targ[itarg].pmra);
    fscanf(fp, "pmdec %lf ", &targets.targ[itarg].pmdec);
    fscanf(fp, "pmra_err %lf ", &targets.targ[itarg].pmra_err);
    fscanf(fp, "pmdec_err %lf ", &targets.targ[itarg].pmdec_err);
    fscanf(fp, "parallax %f ", &targets.targ[itarg].parallax);
    fscanf(fp, "para_err %f ", &targets.targ[itarg].para_err);
    fscanf(fp, "spectyp %s ", targets.targ[itarg].spectyp);
  }
  targets.revision = 1;

  /* Read info for OI_WAVELENGTH table */
  fscanf(fp, "OI_WAVELENGTH insname %s ", wave.insname);
  fscanf(fp, "nwave %d ", &wave.nwave);
  wave.eff_wave = malloc(wave.nwave*sizeof(float));
  wave.eff_band = malloc(wave.nwave*sizeof(float));
  fscanf(fp, "eff_wave ");
  for(i=0; i<wave.nwave; i++) {
    fscanf(fp, "%f ", &wave.eff_wave[i]);
  }
  fscanf(fp, "eff_band ");
  for(i=0; i<wave.nwave; i++) {
    fscanf(fp, "%f ", &wave.eff_band[i]);
  }
  wave.revision = 1;

  /* Read info for OI_VIS table */
  fscanf(fp, "OI_VIS date-obs %s ", vis.date_obs);
  fscanf(fp, "arrname %s insname %s ", vis.arrname, vis.insname);
  fscanf(fp, "numrec %ld ", &vis.numrec);
  vis.record = malloc(vis.numrec*sizeof(oi_vis_record));
  printf("Reading %ld vis records...\n", vis.numrec);
  /* loop over records */
  for(irec=0; irec<vis.numrec; irec++) {
    fscanf(fp, "target_id %d time %lf mjd %lf ", &vis.record[irec].target_id,
	   &vis.record[irec].time, &vis.record[irec].mjd);
    fscanf(fp, "int_time %lf visamp ", &vis.record[irec].int_time);
    vis.record[irec].visamp = malloc(wave.nwave*sizeof(DATA));
    for(iwave=0; iwave<wave.nwave; iwave++) {
      fscanf(fp, "%lf ", &vis.record[irec].visamp[iwave]);
    }
    fscanf(fp, "visamperr ");
    vis.record[irec].visamperr = malloc(wave.nwave*sizeof(DATA));
    for(iwave=0; iwave<wave.nwave; iwave++) {
      fscanf(fp, "%lf ", &vis.record[irec].visamperr[iwave]);
    }
    fscanf(fp, "visphi ");
    vis.record[irec].visphi = malloc(wave.nwave*sizeof(DATA));
    for(iwave=0; iwave<wave.nwave; iwave++) {
      fscanf(fp, "%lf ", &vis.record[irec].visphi[iwave]);
    }
    fscanf(fp, "visphierr ");
    vis.record[irec].visphierr = malloc(wave.nwave*sizeof(DATA));
    for(iwave=0; iwave<wave.nwave; iwave++) {
      fscanf(fp, "%lf ", &vis.record[irec].visphierr[iwave]);
    }
    fscanf(fp, "ucoord %lf vcoord %lf ", &vis.record[irec].ucoord,
	   &vis.record[irec].vcoord);
    fscanf(fp, "sta_index %d %d ", &vis.record[irec].sta_index[0],
	   &vis.record[irec].sta_index[1]);
    vis.record[irec].flag = malloc(wave.nwave*sizeof(char));
    for(iwave=0; iwave<wave.nwave; iwave++) {
      vis.record[irec].flag[iwave] = FALSE;
    }
  }
  vis.revision = 1;
  vis.nwave = wave.nwave;

  /* Read info for OI_VIS2 table */
  fscanf(fp, "OI_VIS2 date-obs %s ", vis2.date_obs);
  fscanf(fp, "arrname %s insname %s ", vis2.arrname, vis2.insname);
  fscanf(fp, "numrec %ld ", &vis2.numrec);
  vis2.record = malloc(vis2.numrec*sizeof(oi_vis2_record));
  printf("Reading %ld vis2 records...\n", vis2.numrec);
  /* loop over records */
  for(irec=0; irec<vis2.numrec; irec++) {
    fscanf(fp, "target_id %d time %lf mjd %lf ", &vis2.record[irec].target_id,
	   &vis2.record[irec].time, &vis2.record[irec].mjd);
    fscanf(fp, "int_time %lf vis2data ", &vis2.record[irec].int_time);
    vis2.record[irec].vis2data = malloc(wave.nwave*sizeof(DATA));
    for(iwave=0; iwave<wave.nwave; iwave++) {
      fscanf(fp, "%lf ", &vis2.record[irec].vis2data[iwave]);
    }
    fscanf(fp, "vis2err ");
    vis2.record[irec].vis2err = malloc(wave.nwave*sizeof(DATA));
    for(iwave=0; iwave<wave.nwave; iwave++) {
      fscanf(fp, "%lf ", &vis2.record[irec].vis2err[iwave]);
    }
    fscanf(fp, "ucoord %lf vcoord %lf ", &vis2.record[irec].ucoord,
	   &vis2.record[irec].vcoord);
    fscanf(fp, "sta_index %d %d ", &vis2.record[irec].sta_index[0],
	   &vis2.record[irec].sta_index[1]);
    vis2.record[irec].flag = malloc(wave.nwave*sizeof(char));
    for(iwave=0; iwave<wave.nwave; iwave++) {
      vis2.record[irec].flag[iwave] = FALSE;
    }
  }
  vis2.revision = 1;
  vis2.nwave = wave.nwave;

  /* Read info for OI_T3 table */
  fscanf(fp, "OI_T3 date-obs %s ", t3.date_obs);
  fscanf(fp, "arrname %s insname %s ", t3.arrname, t3.insname);
  fscanf(fp, "numrec %ld ", &t3.numrec);
  t3.record = malloc(t3.numrec*sizeof(oi_t3_record));
  printf("Reading %ld t3 records...\n", t3.numrec);
  /* loop over records */
  for(irec=0; irec<t3.numrec; irec++) {
    fscanf(fp, "target_id %d time %lf mjd %lf ", &t3.record[irec].target_id,
	   &t3.record[irec].time, &t3.record[irec].mjd);
    fscanf(fp, "int_time %lf t3amp ", &t3.record[irec].int_time);
    t3.record[irec].t3amp = malloc(wave.nwave*sizeof(DATA));
    for(iwave=0; iwave<wave.nwave; iwave++) {
      fscanf(fp, "%lf ", &t3.record[irec].t3amp[iwave]);
    }
    fscanf(fp, "t3amperr ");
    t3.record[irec].t3amperr = malloc(wave.nwave*sizeof(DATA));
    for(iwave=0; iwave<wave.nwave; iwave++) {
      fscanf(fp, "%lf ", &t3.record[irec].t3amperr[iwave]);
    }
    fscanf(fp, "t3phi ");
    t3.record[irec].t3phi = malloc(wave.nwave*sizeof(DATA));
    for(iwave=0; iwave<wave.nwave; iwave++) {
      fscanf(fp, "%lf ", &t3.record[irec].t3phi[iwave]);
    }
    fscanf(fp, "t3phierr ");
    t3.record[irec].t3phierr = malloc(wave.nwave*sizeof(DATA));
    for(iwave=0; iwave<wave.nwave; iwave++) {
      fscanf(fp, "%lf ", &t3.record[irec].t3phierr[iwave]);
    }
    fscanf(fp, "u1coord %lf v1coord %lf ", &t3.record[irec].u1coord,
	   &t3.record[irec].v1coord);
    fscanf(fp, "u2coord %lf v2coord %lf ", &t3.record[irec].u2coord,
	   &t3.record[irec].v2coord);
    fscanf(fp, "sta_index %d %d %d ", &t3.record[irec].sta_index[0],
	   &t3.record[irec].sta_index[1],
	   &t3.record[irec].sta_index[2]);
    t3.record[irec].flag = malloc(wave.nwave*sizeof(char));
    for(iwave=0; iwave<wave.nwave; iwave++) {
      t3.record[irec].flag[iwave] = FALSE;
    }
  }
  t3.revision = 1;
  t3.nwave = wave.nwave;

  fclose(fp);

  /* Write out FITS file */
  printf("Enter output OIFITS filename: ");
  fgets(filename, FLEN_FILENAME, stdin);
  filename[strlen(filename)-1] = '\0'; /* zap newline */
  printf("Writing FITS file %s...\n", filename);
  status = 0;
  fits_create_file(&fptr, filename, &status);
  if (status) {
    fits_report_error(stderr, status);
    exit(EXIT_FAILURE);
  }
  write_oi_target(fptr, targets, &status);
  write_oi_vis(fptr, vis, 1, &status);
  write_oi_vis2(fptr, vis2, 1, &status);
  write_oi_t3(fptr, t3, 1, &status);
  write_oi_array(fptr, array, 1, &status);
  write_oi_wavelength(fptr, wave, 1, &status);

  if (status) {
    /* Error occurred - delete partially-created file and exit */
    fits_delete_file(fptr, &status);
    exit(EXIT_FAILURE);
  } else {
    fits_close_file(fptr, &status);
  }

  /* Free storage */
  free_oi_target(&targets);
  free_oi_vis(&vis);
  free_oi_vis2(&vis2);
  free_oi_t3(&t3);
  free_oi_array(&array);
  free_oi_wavelength(&wave);
}


/**
 * Read OIFITS file.
 *
 * This code will read a general OIFITS file containing multiple
 * tables of each type (except OI_TARGET). However, its not that
 * useful since the same oi_vis/vis2/t3 object is used to receive the
 * data from all tables of the corresponding type!  See oifile.c for a
 * more useful routine (which requires GLib).
 * 
 */
void demo_read(void)
{
  oi_array array;
  oi_target targets;
  oi_wavelength wave;
  oi_vis vis;
  oi_vis2 vis2;
  oi_t3 t3;
  char filename[FLEN_FILENAME];
  int status, hdutype;
  fitsfile *fptr, *fptr2;

  printf("Enter input OIFITS filename: ");
  fgets(filename, FLEN_FILENAME, stdin);
  filename[strlen(filename)-1] = '\0'; /* zap newline */
  printf("Reading FITS file %s...\n", filename);
  status = 0;
  fits_open_file(&fptr, filename, READONLY, &status);
  if (status) {
    fits_report_error(stderr, status);
    goto except;
  }
  read_oi_target(fptr, &targets, &status);
  if (status) goto except;
  /* open 2nd connection to read referenced OI_ARRAY and OI_WAVELENGTH tables
     without losing place in file */
  fits_open_file(&fptr2, filename, READONLY, &status);

  /* Read all OI_VIS tables & corresponding OI_ARRAY/OI_WAVELENGTH tables */
  while (1==1) {
    if (read_next_oi_vis(fptr, &vis, &status)) break; /* no more OI_VIS */
    printf("Read OI_VIS  with  ARRNAME=%s INSNAME=%s\n",
	   vis.arrname, vis.insname);
    if (strlen(vis.arrname) > 0) {
      /* if ARRNAME specified, read corresponding OI_ARRAY
         Note we may have read it previously */
      read_oi_array(fptr2, vis.arrname, &array, &status);
    }
    read_oi_wavelength(fptr2, vis.insname, &wave, &status);
    if (!status) {
      /* Free storage ready to reuse structs for next table */
      free_oi_wavelength(&wave);
      if (strlen(vis.arrname) > 0) free_oi_array(&array);
      free_oi_vis(&vis);
    }
  }
  if (status != END_OF_FILE) goto except;
  status = 0;

  /* Read all OI_VIS2 tables & corresponding OI_ARRAY/OI_WAVELENGTH tables */
  fits_movabs_hdu(fptr, 1, &hdutype, &status); /* back to start */
  while (1==1) {
    if (read_next_oi_vis2(fptr, &vis2, &status)) break; /* no more OI_VIS2 */
    printf("Read OI_VIS2 with  ARRNAME=%s INSNAME=%s\n",
	   vis2.arrname, vis2.insname);
    if (strlen(vis2.arrname) > 0) {
      /* if ARRNAME specified, read corresponding OI_ARRAY */
      read_oi_array(fptr2, vis2.arrname, &array, &status);
    }
    read_oi_wavelength(fptr2, vis2.insname, &wave, &status);
    if (!status) {
      /* Free storage ready to reuse structs for next table */
      free_oi_wavelength(&wave);
      if (strlen(vis2.arrname) > 0) free_oi_array(&array);
      free_oi_vis2(&vis2);
    }
  }
  if (status != END_OF_FILE) goto except;
  status = 0;

  /* Read all OI_T3 tables & corresponding OI_ARRAY/OI_WAVELENGTH tables */
  fits_movabs_hdu(fptr, 1, &hdutype, &status); /* back to start */
  while (1==1) {
    if (read_next_oi_t3(fptr, &t3, &status)) break; /* no more OI_T3 */
    printf("Read OI_T3   with  ARRNAME=%s INSNAME=%s\n",
	   t3.arrname, t3.insname);
    if (strlen(t3.arrname) > 0) {
      /* if ARRNAME specified, read corresponding OI_ARRAY */
      read_oi_array(fptr2, t3.arrname, &array, &status);
    }
    read_oi_wavelength(fptr2, t3.insname, &wave, &status);
    if (!status) {
      /* Free storage ready to reuse structs for next table */
      free_oi_wavelength(&wave);
      if (strlen(t3.arrname) > 0) free_oi_array(&array);
      free_oi_t3(&t3);
    }
  }
  if (status != END_OF_FILE) goto except;
  status = 0;

  fits_close_file(fptr, &status);
  fits_close_file(fptr2, &status);
  free_oi_target(&targets);
  return;

 except:
  exit(EXIT_FAILURE);
}


/**
 * Main function for demonstration program.
 */
int main(int argc, char *argv[]) 
{
  demo_write();
  demo_read();
  exit(EXIT_SUCCESS);
}
