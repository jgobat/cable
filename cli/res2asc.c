/*
    This file is part of WHOI Cable, a program for the static and dynamic
    analysis of oceanographic cable structures.

    Copyright (C) 1997-2016 by Woods Hole Oceanographic Institution (WHOI)
    and Jason Gobat

    WHOI Cable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    WHOI Cable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with WHOI Cable.  If not, see <http://www.gnu.org/licenses/>.
*/

/****************************************************************************
 *
 * File:        res2asc.c
 *
 * Description: The driver (which is really all there is) for a program
 *		that reads cable binary results files and spits out
 *		tabulated results in ASCII format.
 *
 * History:
 *
 ****************************************************************************/

# include <stdio.h>
# include <malloc.h>
# include <string.h>
# include <stdlib.h>
# include "compress.h"
# include "options.h"
# include "rotate.h"
# include "error.h"
# include "results.h"

void *problem = NULL;
void *environment = NULL;

static int variable_group []   = {0, DISPLACEMENT, DISPLACEMENT, DISPLACEMENT, 
                                     VELOCITY, VELOCITY, VELOCITY, 
              			     FORCE, FORCE, FORCE, 
				     MOMENT, MOMENT, MOMENT, 
				     EULER, EULER, EULER, EULER};

static char *local3D_names []  = {"s", "z", "x", "y", 
				       "u", "v", "w", 
				       "T", "Sn", "Sb",
                                       "Mt", "Mn", "Mb",
                                       "B0", "B1", "B2", "B3"};
static char *global3D_names [] = {"s", "z", "x", "y", 
				       "W", "U", "V",
				       "Fz", "Fx", "Fy",
                                       "Mz", "Mx", "My",
                                       "B0", "B1", "B2", "B3"}; 
static char *local2D_names []  = {"s", "z", "x", NULL, 
				       "u", "v", NULL, 
				       "T", "Sn", NULL, 
				       NULL, NULL, "Mb", 
				       "phi", NULL, NULL, NULL};
static char *global2D_names [] = {"s", "z", "x", NULL,
			               "W", "U", NULL,
				       "Fz", "Fx", NULL,
				       NULL, NULL, "My", 
				       "phi", NULL, NULL, NULL};

static int	default_global = 0;

static char *usage_msg = "\
usage: res2asc [-in in_file] [-variables v1 v2 ...] [output files] [options]\n\
\n\
  output files (at least one must be specified):\n\
\n\
  -static filename   output filename for static results\n\
  -time filename     output filename for nodal time series results\n\
  -snap filename     output filename for time history snapshot results\n\
\n\
  boolean options:\n\
\n\
  -global            rotate results into global (XYZ) coordinates (default=off)\n\
  -totals            present dynamic results as static+dynamic (default=off)\n\
  -lbs               convert force units from N to lbf (default=off)\n\
  -ft                convert length units from m to ft (default=off)\n\
  -header            print a header on all tables (default=on)\n\
\n\
  valid variable names are: \n\
\n\
    2D prob in local coord:  s, z, x, T, Sn, Mb, phi\n\
    2D prob in global coord: s, z, x, Fz, Fx, My, phi\n\
    3D prob in local coord:  s, z, x, y, T, Sn, Sb, Mt, Mn, Mb, B0, B1, B2, B3\n\
    3D prob in global coord: s, z, x, y, Fz, Fx, Fy, Mz, Mx, My, B0, B1, B2, B3\n\
";

static int usage(msg)
   char *msg;
{
   fprintf (stderr,"res2asc: %s\n\n", msg);
   fputs(usage_msg, stderr);
   exit (1);
}

int main (argc, argv)
   int	 argc;
   char *argv [ ];
{
   char	       **var_names;
   char		*in_name;
   char		*time_name;
   char		*static_name;
   char		*snap_name;
   char	       **variables;
   int		 num_variables;
   int		*var_pos;
   int		 n;
   int		 global;
   int		 totals;
   int		 lbs;
   int		 header;
   int		 ft;
   ResFile       in;
   FILE		*static_fp;
   FILE		*time_fp;
   FILE		*snap_fp;
   int		 i, j, k;
   double      **d_st;
   double     ***d_time;
   double     **d_snap;
   double     ***d_timeres;
   double     **d_snapres;
   Result       *res;


   n = GetSoloStringOption(argc, argv, "in", &in_name);
   if (n <= 0)
      usage ("must specify an input file (-in)");

   static_name = NULL;
   n = GetSoloStringOption(argc, argv, "static", &static_name);
   if (n < 0)
      usage ("error in static output file option");

   if (static_name && strcmp(in_name, static_name) == 0)
      usage ("input and output file names must be different");

   time_name = NULL;
   n = GetSoloStringOption(argc, argv, "time", &time_name);
   if (n < 0)
      usage ("error in time output file option");

   if (time_name && strcmp(in_name, time_name) == 0)
      usage ("input and output file names must be different");

   snap_name = NULL;
   n = GetSoloStringOption(argc, argv, "snap", &snap_name);
   if (n < 0)
      usage ("error in snap output file option");

   if (snap_name && strcmp(in_name, snap_name) == 0)
      usage ("input and output file names must be different");

   num_variables = GetStringOption (argc, argv, "variables", &variables);
   if (num_variables <= 0)
      usage ("must specify a list of variables to tabulate");

   global = default_global;
   n = GetBooleanOption (argc, argv, "global", &global);
   if (n < 0)
      usage ("error in global boolean option");

   totals = 0;
   n = GetBooleanOption (argc, argv, "totals", &totals);
   if (n < 0)
      usage ("error in totals boolean option");

   lbs = 0;
   n = GetBooleanOption (argc, argv, "lbs", &lbs);
   if (n < 0)
      usage ("error in lbs boolean option");

   header = 1;
   n = GetBooleanOption (argc, argv, "header", &header);
   if (n < 0)
      usage ("error in header boolean option");

   ft = 0;
   n = GetBooleanOption (argc, argv, "ft", &ft);
   if (n < 0)
      usage ("error in ft boolean option");

   if (!static_name && !time_name && !snap_name) {
      fprintf (stderr,"res2asc: must specify at least one kind of output file\n");
      exit (1);
   }

   in = res_open(in_name, "rb");
   if (in == NULL) {
      fprintf (stderr,"res2asc: could not open input file for reading\n");
      exit (1);
   }

   res = ReadResultsFile(in, global, totals, lbs, ft);
   if (!res)
      exit (1);


   if (static_name) {
      static_fp = fopen(static_name, "w");
      if (static_fp == NULL) 
         Fatal("res2asc: could not open static output file for reading");
   }
   else
      static_fp = NULL;

   if (time_name && res -> sample_dt) {
      time_fp = fopen(time_name, "w");
      if (time_fp == NULL) 
         Fatal("res2asc: could not open time output file for reading");
   }
   else
      time_fp = NULL;

   if (snap_name && res -> snap_dt) {
      snap_fp = fopen(snap_name, "w");
      if (snap_fp == NULL) 
         Fatal("res2asc: could not open snap output file for reading");
   }
   else
      snap_fp = NULL;

   if (global && res -> twoD)
      var_names = global2D_names;
   else if (global && !res -> twoD)
      var_names = global3D_names;
   else if (!global && res -> twoD)
      var_names = local2D_names;
   else if (!global && !res -> twoD)
      var_names = local3D_names;

   d_time = (double ***) malloc(sizeof(double **) * num_variables);
   d_snap = (double **) malloc(sizeof(double *) * num_variables);
   d_st = (double **) malloc(sizeof(double *) * num_variables);

   d_snapres = (double **) malloc(sizeof(double *) * 17);
   d_timeres = (double ***) malloc(sizeof(double **) * 17);

   d_snapres [0] = NULL;
   d_snapres [1] = res -> y;
   d_snapres [2] = res -> x;
   d_snapres [3] = res -> z;
   d_snapres [4] = res -> velocity [0];
   d_snapres [5] = res -> velocity [1];
   d_snapres [6] = res -> velocity [2];
   d_snapres [7] = res -> force [0];
   d_snapres [8] = res -> force [1];
   d_snapres [9] = res -> force [2];
   d_snapres [10] = res -> moment [0];
   d_snapres [11] = res -> moment [1];
   d_snapres [12] = res -> moment [2];
   d_snapres [13] = res -> beta [0];
   d_snapres [14] = res -> beta [1];
   d_snapres [15] = res -> beta [2];
   d_snapres [16] = res -> beta [3];
   
   d_timeres [0] = NULL; 
   d_timeres [1] = res -> y_t;
   d_timeres [2] = res -> x_t;
   d_timeres [3] = res -> z_t;
   d_timeres [4] = res -> velocity_t [0];
   d_timeres [5] = res -> velocity_t [1];
   d_timeres [6] = res -> velocity_t [2];
   d_timeres [7] = res -> force_t [0];
   d_timeres [8] = res -> force_t [1];
   d_timeres [9] = res -> force_t [2];
   d_timeres [10] = res -> moment_t [0];
   d_timeres [11] = res -> moment_t [1];
   d_timeres [12] = res -> moment_t [2];
   d_timeres [13] = res -> beta_t [0];
   d_timeres [14] = res -> beta_t [1];
   d_timeres [15] = res -> beta_t [2];
   d_timeres [16] = res -> beta_t [3];

   var_pos = (int *) malloc(sizeof(int) * num_variables);

   for (i = 0 ; i < num_variables ; i++) {
      var_pos [i] = -1;

      for (j = 0 ; j < 17 ; j++) {
         if (var_names [j] == NULL)
            continue;

         if (strcmp(var_names [j], variables [i]) == 0) {
            if (strcmp(var_names [j], "s") == 0) {
               d_time [i] = NULL;
               d_snap [i] = NULL;
               d_st [i] = res -> s;
               var_pos [i] = j;
            }
            else if (strcmp(var_names [j], "u") == 0 ||
                strcmp(var_names [j], "v") == 0 ||
                strcmp(var_names [j], "w") == 0) {

               d_time [i] = d_timeres [j];
               d_snap [i] = d_snapres [j];
               d_st [i] = NULL;
               var_pos [i] = j;

               if (!res -> output_map [variable_group [j]])
                  Fatal("variable %s not stored in output file", var_names [j]);
            }
            else {
               d_time [i] = d_timeres [j];
               d_snap [i] = d_snapres [j];
               d_st [i] = d_snapres [j];
               var_pos [i] = j;

               if (!res -> output_map [variable_group [j]])
                  Fatal("variable %s not stored in output file", var_names [j]);
           }

            break;
         }
      }

      if (var_pos [i] == -1)
         Fatal ("could not match variable name %s", variables [i]);
   }

	/*
	 * write out the static results
	 */

   if (static_fp) {
      if (header) {
         for (i = 0 ; i < num_variables ; i++) {
            if (d_st [i] == NULL)
               continue;

            if (res -> depth_ref && strcmp(var_names [var_pos [i]], "z") == 0)
               fprintf (static_fp,"    %-3s       ", "dep");
            else 
               fprintf (static_fp,"    %-3s       ", var_names [var_pos [i]]);
         }

         fprintf  (static_fp, "\n");

         for (i = 0 ; i < num_variables ; i++) 
            if (d_st [i] != NULL)
               fprintf (static_fp,"------------  ");
  
         fprintf  (static_fp, "\n");
      }

      for (i = 0 ; i < res -> npoints ; i++) {
         for (j = 0 ; j < num_variables ; j++) {
            if (d_st [j] == NULL)
               continue; 
   
            fprintf (static_fp ,"%11.5f   ", d_st [j][i]); 
         }

         fprintf (static_fp, "\n");
      } 

      fclose(static_fp);
   }

   if (!res -> dynamic && (snap_fp || time_fp)) 
      Fatal("res2asc: no dynamic results available for conversion");
  
   if (snap_fp) {
      if (header) {
         for (i = 0 ; i < num_variables ; i++) {
            if (d_snap [i] == NULL)
               continue;

            if (res -> depth_ref && strcmp(var_names [var_pos [i]], "z") == 0)
               fprintf (snap_fp,"    %-3s       ", "dep");
            else 
               fprintf (snap_fp,"    %-3s       ", var_names [var_pos [i]]);
         }

         fprintf  (snap_fp, "\n");

         for (i = 0 ; i < num_variables ; i++) {
            if (d_snap [i] != NULL)
               fprintf (snap_fp,"------------  ");
         }

         fprintf  (snap_fp, "\n");
      } 

      for (i = 0 ; i < res -> nsteps ; i++) {

         ReadResultSnapshot(res, i, totals, lbs, ft);

         for (j = 0 ; j < res -> npoints ; j++) {
            for (k = 0 ; k < num_variables ; k++) {
               if (d_snap [k] == NULL)
                  continue;

               fprintf (snap_fp,"%11.5f   ", d_snap [k][j]);
            }
            fprintf (snap_fp, "\n");
         }
         fprintf (snap_fp, "\n");
      }

      fclose(snap_fp);
   }

   if (time_fp) {
      if (header) {
         fprintf (time_fp,"    %-3s       ", "t");

         n = 0;
         for (j = 0 ; j < res -> num_output_nodes ; j++) {
            for (i = 0 ; i < num_variables ; i++) {
               if (d_time [i] == NULL) 
                  continue;

               fprintf (time_fp,"    %-3s(%3d)  ", var_names [var_pos [i]],
                        res -> output_nodes [j]);
               n++;
            }
         }

         fprintf (time_fp, "\n");

         for (i = 0 ; i < n + 1 ; i++)
            fprintf (time_fp,"------------  ");
   
         fprintf (time_fp, "\n");
      }

      for (k = 0 ; k < res -> nsamples ; k++) {
         fprintf (time_fp,"%11.5f   ", k*res -> sample_dt);

         for (j = 0 ; j < res -> num_output_nodes ; j++) {
            for (i = 0 ; i < num_variables ; i++) {
               if (d_time [i] == NULL)
                  continue;

               fprintf (time_fp,"%11.5f   ", d_time [i][j][k]);
            }
         }
         fprintf (time_fp, "\n");
      }
      fclose(time_fp);
   }

   res_close(in);


   exit (0);
}
