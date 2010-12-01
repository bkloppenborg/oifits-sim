/* $Id: oifits-filter.c,v 1.4 2008-06-09 16:21:04 jsy1001 Exp $ */

/**
 * @file oifits-filter.c
 *
 * Command-line OIFITS filter utility.
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

#include "oifilter.h"


/**
 * Main function for command-line filter utility.
 */
int main(int argc, char *argv[]) 
{
#ifndef HAVE_G_OPTION_GROUP
  printf("Need GLib >= 2.6\n");
#else

  GError *error;
  GOptionContext *context;
  char inFilename[FLEN_FILENAME], outFilename[FLEN_FILENAME];
  oi_fits inData, outData;
  int status;

  /* Parse command-line */
  error = NULL;
  context =
    g_option_context_new("INFILE OUTFILE - write filtered dataset to new file");
  g_option_context_set_main_group(context, get_oi_filter_option_group());
  g_option_context_parse(context, &argc, &argv, &error);
  if(error != NULL) {
    printf("Error parsing command-line options: %s\n", error->message);
    g_error_free(error);
    exit(2); /* standard unix behaviour */
  }
  if(argc != 3) {
    printf("Wrong number of command-line arguments\n"
	   "Enter '%s --help' for usage information\n", argv[0]);
    exit(2);
  }
  g_strlcpy(inFilename, argv[1], FLEN_FILENAME);
  g_strlcpy(outFilename, argv[2], FLEN_FILENAME);

  /* Read FITS file */
  status = 0;
  read_oi_fits(inFilename, &inData, &status);
  if(status) goto except;

  /* Display summary info */
  printf("=== INPUT DATA: ===\n");
  print_oi_fits_summary(&inData);
  printf("=== Applying filter: ===\n");
  print_oi_filter(get_user_oi_filter());

  /* Apply filter */
  apply_user_oi_filter(&inData, &outData);
  printf("--> OUTPUT DATA: ===\n");
  print_oi_fits_summary(&outData);

  /* Write out filtered data */
  write_oi_fits(outFilename, outData, &status);
  if(status) goto except;

  free_oi_fits(&inData);
  free_oi_fits(&outData);
  g_option_context_free(context);
  exit(EXIT_SUCCESS);

 except:
  exit(EXIT_FAILURE);
#endif /* #ifndef HAVE_G_OPTION_GROUP */
}
