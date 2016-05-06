/*
    This file is part of WHOI Cable, a program for the static and dynamic
    analysis of oceanographic cable structures -> 

    Copyright (C) 1997-2016 by Woods Hole Oceanographic Institution (WHOI)

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
 * File:        res2mat.c
 *
 * Description: The driver (which is really all there is) for a program
 *              that reads cable binary results files and writes
 *		a binary Matlab (v4) file containing all of the results.
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
# include "problem.h"
# include "results.h"
# include "error.h"

typedef struct {
   int	type;
   int	mrows;
   int ncols;
   int imagf;
   int namlen;
} MATheader;

static int architecture ( )
{
   int	x = 1;

   if (*((char *) &x) == 1)
      return 0;
   else
      return 1;
}

static void MatlabMatrix (a, nr, nc, tpose, name, fp)
   double      **a;
   int		 nr;
   int		 nc;
   int		 tpose;
   char		*name;
   FILE		*fp;
{
   int		arch;
   int		mopt;
   double	x;
   unsigned	i, j;
   MATheader	h;

   arch = architecture ( );
   
   mopt = arch*1000 + 0*100 + 0*10 + 0*1;
                      /* reserved */
                              /* double precision */
                                     /* numeric full matrix */
   h.type = mopt;
   h.mrows = nr;
   h.ncols = nc;
   h.imagf = 0;
   h.namlen = strlen(name) + 1;
       
   fwrite (&h, sizeof(MATheader), 1, fp);
   fwrite (name, sizeof(char), h.namlen, fp);

   for (i = 0 ; i < nc ; i++) {
      for (j = 0 ; j < nr ; j++) {
         if (tpose)
            x = a [i][j];
         else
            x = a [j][i];

         fwrite (&x, sizeof(double), 1, fp);
      }
   }

   return;
}

static void MatlabVector (a, n, name, fp)
   double	*a;
   int		 n;
   char		*name;
   FILE		*fp;
{
   int		arch;
   int		mopt;
   double	x;
   unsigned	i;
   MATheader	h;

   arch = architecture ( );
   

   mopt = arch*1000 + 0*100 + 0*10 + 0*1;
                      /* reserved */
                              /* double precision */
                                     /* numeric full matrix */
   h.type = mopt;
   h.mrows = n;
   h.ncols = 1;
   h.imagf = 0;
   h.namlen = strlen(name) + 1;
       
   fwrite (&h, sizeof(MATheader), 1, fp);
   fwrite (name, sizeof(char), h.namlen, fp);
   
   for (i = 0 ; i < n ; i++) {
      x = a [i];
      fwrite (&x, sizeof(double), 1, fp);
   }

   return;
}

static int 
ProcessBranchInfo (Result *res, FILE *out)
{
   int	    i;
   double **info;

   info = array(res -> num_branch + 1, 2);

   info [0][0] = 1;
   for (i = 1 ; i <= res -> num_branch ; i++) {
      info [i - 1][1] = res -> branch_starts [i - 1] - 1;
      info [i][0] = res -> branch_starts [i - 1];
   }

   info [res -> num_branch][1] = res -> npoints;

   MatlabMatrix (info, res -> num_branch + 1, 2, 0, "br", out);  

   return 0; 
}

static int 
ProcessDynamic (Result *res, FILE *out, 
                int global, int lbs, int totals, int twoD)
{
   int	         i, n, j;
   double	*nodes;
   double	*tm;
   double  **force[3], **vel[3], **moment[3], **beta[4], **x, **y, **z;

   for (i = 0 ; i < 3 ; i++) {
      force[i] = vel[i] = moment[i] = beta[i] = NULL;
   }
   beta[3] = x = y = z = NULL;

   if (res -> nsamples && res -> num_output_nodes) {
      nodes = (double *) malloc(sizeof(double) *res -> num_output_nodes);
      if (nodes == NULL)
          Fatal("could not allocate memory for output nodes");

      for (i = 0 ; i < res -> num_output_nodes ; i++) 
         nodes [i] = res -> output_nodes [i];

      tm = (double *) malloc(sizeof(double) * res -> nsamples);
      if (tm == NULL)
          Fatal("could not allocate memory for time vector");

      for (i = 0 ; i < res -> nsamples ; i++)
         tm [i] = res -> t_start + i*res -> sample_dt;

      MatlabVector (nodes, res -> num_output_nodes, "nodes", out);
      MatlabVector (tm, res -> nsamples, "t", out);
      free (nodes);
      free (tm);

      MatlabVector (&(res -> sample_dt), 1, "dt", out);
   }
  
   if (res -> snap_dt && res -> nsteps) {
      MatlabVector (&(res -> snap_dt), 1, "snap_dt", out);

      if (res -> output_map[DISPLACEMENT]) {
         x = array(res -> nsteps, res -> npoints);
         y = array(res -> nsteps, res -> npoints);
         z = array(res -> nsteps, res -> npoints);
      }
      for (i = 0 ; i < 3 ; i++) {
         if (res -> output_map[VELOCITY]) {
            vel[i] = array(res -> nsteps, res -> npoints);
         }
         if (res -> output_map[FORCE]) {
            force[i] = array(res -> nsteps, res -> npoints);
         }
         if (res -> output_map[MOMENT]) {
            moment[i] = array(res -> nsteps, res -> npoints);
         }
         if (res -> output_map[EULER]) {
            beta[i] = array(res -> nsteps, res -> npoints);
         }
      }
      if (res -> output_map[EULER]) {
         beta[3] = array(res -> nsteps, res -> npoints);
      }

      n = sizeof(double) * res -> npoints;
      for (i = 0 ; i < res -> nsteps ; i++) {
          ReadResultSnapshot(res, i, totals, lbs, 0);
          if (global)
             SetSnapGlobal(res); 

          if (res -> output_map[DISPLACEMENT]) {
              memcpy(x[i], res -> x, n);
              memcpy(y[i], res -> y, n);
              memcpy(z[i], res -> z, n);
          }
          for (j = 0 ; j < 3 ; j++) {
              if (res -> output_map[VELOCITY]) {
                  memcpy(vel[j][i], res -> velocity[j], n);
              }
              if (res -> output_map[FORCE]) {
                  memcpy(force[j][i], res -> force[j], n);
              }
              if (res -> output_map[MOMENT]) {
                  memcpy(moment[j][i], res -> moment[j], n);
              }
              if (res -> output_map[EULER]) {
                  memcpy(beta[j][i], res -> beta[j], n);
              }
          }
          if (res -> output_map[EULER]) {
              memcpy(beta[3][i], res -> beta[3], n);
          }
      }
   }

   if (res -> output_map [DISPLACEMENT] && res -> num_output_nodes && res -> nsamples) {
      MatlabMatrix (res -> y_t, res -> nsamples, res -> num_output_nodes, 1, "z_t", out);
      MatlabMatrix (res -> x_t, res -> nsamples, res -> num_output_nodes, 1, "x_t", out);
      if (!twoD)
         MatlabMatrix (res -> z_t, res -> nsamples, res -> num_output_nodes, 1, "y_t", out);
   }
   if (res -> output_map [DISPLACEMENT] && res -> nsteps) {
      MatlabMatrix (y, res -> npoints, res -> nsteps, 1, "z_s", out);
      MatlabMatrix (x, res -> npoints, res -> nsteps, 1, "x_s", out);
      if (!twoD)
         MatlabMatrix (z, res -> npoints, res -> nsteps, 1, "y_s", out);
   }
   if (res -> output_map [VELOCITY] && res -> num_output_nodes && res -> nsamples) {
       MatlabMatrix (res -> velocity_t [0], res -> nsamples, res -> num_output_nodes, 1, global ? "W_t" : "u_t", out);
       MatlabMatrix (res -> velocity_t [1], res -> nsamples, res -> num_output_nodes, 1, global ? "U_t" : "w_t", out);
       if (!twoD)
           MatlabMatrix (res -> velocity_t [2], res -> nsamples, res -> num_output_nodes, 1, global ? "V_t" : "v_t", out);
   }
   if (res -> output_map [VELOCITY] && res -> nsteps && vel[0]) {
      MatlabMatrix (vel[0], res -> npoints, res -> nsteps, 1, global ? "W_s" : "u_s", out);
      MatlabMatrix (vel[1], res -> npoints, res -> nsteps, 1, global ? "U_s" : "w_s", out);
      if (!twoD)
         MatlabMatrix (vel[2], res -> npoints, res -> nsteps, 1, global ? "V_s" : "v_s", out);
   }
   if (res -> output_map [FORCE] && res -> num_output_nodes && res -> nsamples) {
      MatlabMatrix(res -> force_t [0], res -> nsamples, res -> num_output_nodes, 1, global ? "Fz_t" : "T_t", out);
      MatlabMatrix(res -> force_t [1], res -> nsamples, res -> num_output_nodes, 1, global ? "Fx_t" : "Sn_t", out);
      if (!twoD)
         MatlabMatrix(res -> force_t [2], res -> nsamples, res -> num_output_nodes, 1, global ? "Fy_t" : "Sb_t", out);
   }
   if (res -> output_map [FORCE] && res -> nsteps) {
       MatlabMatrix (force [0], res -> npoints, res -> nsteps, 1, global ? "Fz_s" : "T_s", out);
       MatlabMatrix (force [1], res -> npoints, res -> nsteps, 1, global ? "Fx_s" : "Sn_s", out);
       if (!twoD)
           MatlabMatrix (force [2], res -> npoints, res -> nsteps, 1, global ? "Fy_s" : "Sb_s", out);
   }
   if (res -> output_map [MOMENT] && res -> num_output_nodes && res -> nsamples) {
      if (!twoD) {
         MatlabMatrix (res -> moment_t [0], res -> nsamples ,res -> num_output_nodes, 1, global ? "Mz_t" : "Mt_t", out);
         MatlabMatrix (res -> moment_t [1], res -> nsamples ,res -> num_output_nodes, 1, global ? "Mx_t" : "Mn_t", out);
      }
      MatlabMatrix (res -> moment_t [2], res -> nsamples ,res -> num_output_nodes, 1, global ? "My_t" : "Mb_t", out);
   }
   if (res -> output_map [MOMENT] && res -> nsteps) {
      if (!twoD) {
         MatlabMatrix (moment [0], res -> npoints, res -> nsteps, 1, global ? "Mz_s" : "Mt_s", out);
         MatlabMatrix (moment [1], res -> npoints, res -> nsteps, 1, global ? "Mx_s" : "Mn_s", out);
      }
      MatlabMatrix (moment [2], res -> npoints, res -> nsteps, 1, global ? "My_s" : "Mb_s", out);
   }
   if (res -> output_map [EULER] &&res -> num_output_nodes && res -> nsamples) {
      if (twoD)  
         MatlabMatrix (res -> beta_t [0], res -> nsamples ,res -> num_output_nodes, 1, "phi_t", out);
      else {
         MatlabMatrix (res -> beta_t [0], res -> nsamples ,res -> num_output_nodes, 1, "B0_t", out);
         MatlabMatrix (res -> beta_t [1], res -> nsamples ,res -> num_output_nodes, 1, "B1_t", out);
         MatlabMatrix (res -> beta_t [2], res -> nsamples ,res -> num_output_nodes, 1, "B2_t", out);
         MatlabMatrix (res -> beta_t [3], res -> nsamples ,res -> num_output_nodes, 1, "B3_t", out);
      }
   }
   if (res -> ext_dt) {
      MatlabMatrix(res -> ext_t[0], res -> n_ext_steps, res -> n_ext, 0, "zthrust_t", out);
      MatlabMatrix(res -> ext_t[1], res -> n_ext_steps, res -> n_ext, 0, "xthrust_t", out);
      if (!twoD) 
         MatlabMatrix(res -> ext_t[2], res -> n_ext_steps, res -> n_ext, 0, "ythrust_t", out);
    
   }
   if (res -> output_map [EULER] && res -> nsteps) {
      if (twoD)
         MatlabMatrix (beta [0], res -> npoints, res -> nsteps, 1, "phi_s", out);
      else {
         MatlabMatrix (beta [0], res -> npoints, res -> nsteps, 1, "B0_s", out);
         MatlabMatrix (beta [1], res -> npoints, res -> nsteps, 1, "B1_s", out);
         MatlabMatrix (beta [2], res -> npoints, res -> nsteps, 1, "B2_s", out);
         MatlabMatrix (beta [3], res -> npoints, res -> nsteps, 1, "B3_s", out);
      }
   }

   if (res -> buoy_dt) {
      MatlabVector (&(res -> buoy_dt), 1, "buoy_dt", out);
    
      MatlabVector(res -> buoy [0], res -> nbuoy, "surge", out);
      MatlabVector(res -> buoy [1], res -> nbuoy, "sway", out);
      MatlabVector(res -> buoy [2], res -> nbuoy, "heave", out);
      MatlabVector(res -> buoy [3], res -> nbuoy, "roll", out);
      MatlabVector(res -> buoy [4], res -> nbuoy, "pitch", out);
      MatlabVector(res -> buoy [5], res -> nbuoy, "yaw", out);
   }

   if (res -> seg_dt) {
       MatlabMatrix(res -> seg_top_pay, res -> nseg_steps, res -> nseg, 0, "top_pay", out);
       MatlabMatrix(res -> seg_bot_pay, res -> nseg_steps, res -> nseg, 0, "bot_pay", out);
       MatlabMatrix(res -> seg_top_spooled, res -> nseg_steps, res -> nseg, 0, "top_spooled", out);
       MatlabMatrix(res -> seg_bot_spooled, res -> nseg_steps, res -> nseg, 0, "bot_spooled", out);
       MatlabMatrix(res -> seg_unstretched, res -> nseg_steps, res -> nseg, 0, "unstretched", out);
       MatlabMatrix(res -> seg_stretched, res -> nseg_steps, res -> nseg, 0, "stretched", out);
       MatlabMatrix(res -> seg_first, res -> nseg_steps, res -> nseg, 0, "first", out);
       MatlabMatrix(res -> seg_last, res -> nseg_steps, res -> nseg, 0, "last", out);
   }

   free_array(x, res -> nsteps, res -> npoints);
   free_array(y, res -> nsteps, res -> npoints);
   free_array(z, res -> nsteps, res -> npoints);
   for (i = 0 ; i < 3 ; i++) {
       free_array(force[i], res -> nsteps, res -> npoints);
       free_array(moment[i], res -> nsteps, res -> npoints);
       free_array(vel[i], res -> nsteps, res -> npoints);
       free_array(beta[i], res -> nsteps, res -> npoints);
   }
   free_array(beta[3], res -> nsteps, res -> npoints);

   return 0;
}

static int 
ProcessStatic (Result *res, FILE *out, int global, int twoD)
{
   MatlabVector (res -> s, res -> npoints, "s", out);
   if (res -> output_map [DISPLACEMENT]) {
      MatlabVector (res -> y_st, res -> npoints, "z", out);
      MatlabVector (res -> x_st, res -> npoints, "x", out);
      if (!twoD)
         MatlabVector (res -> z_st, res -> npoints, "y", out);
   }
   if (res -> output_map [FORCE]) {
      if (global) {
         MatlabVector (res -> force_st [0], res -> npoints, "Fz", out);
         MatlabVector (res -> force_st [1], res -> npoints, "Fx", out);
         if (!twoD)
            MatlabVector (res -> force_st [2], res -> npoints, "Fy", out);
      }
      else {
         MatlabVector (res -> force_st [0], res -> npoints, "T", out);
         MatlabVector (res -> force_st [1], res -> npoints, "Sn", out);
         if (!twoD)
            MatlabVector (res -> force_st [2], res -> npoints, "Sb", out);
      }
   }
   if (res -> output_map [MOMENT]) {
      if (global) {
         if (!twoD) {
            MatlabVector (res -> moment_st [0], res -> npoints, "Mz", out);
            MatlabVector (res -> moment_st [1], res -> npoints, "Mx", out);
         }
         MatlabVector (res -> moment_st [2], res -> npoints, "My", out);
      }
      else {
         if (!twoD) {
            MatlabVector (res -> moment_st [0], res -> npoints, "Mt", out);
            MatlabVector (res -> moment_st [1], res -> npoints, "Mn", out);
         }
         MatlabVector (res -> moment_st [2], res -> npoints, "Mb", out);
      }
   }
   if (res -> output_map [EULER]) {
      if (!twoD) { 
         MatlabVector (res -> beta_st [0], res -> npoints, "B0", out);
         MatlabVector (res -> beta_st [1], res -> npoints, "B1", out);
         MatlabVector (res -> beta_st [2], res -> npoints, "B2", out);
         MatlabVector (res -> beta_st [3], res -> npoints, "B3", out);
      }
      else
         MatlabVector (res -> beta_st [0], res -> npoints, "phi", out);
   }

   return 0;
}

int
SolutionToMatlab(char *resname, char *matname, int global, int lbs, int totals, int use_file_twoD, int requested_twoD)
{
    int		 twoD;
    Result   *res;
    ResFile	 in;
    FILE		*out;

    in = res_open(resname, "rb");
    if (in == NULL) {
        error("res2mat: could not open input file for reading\n");
        return 1;
    }
 
    res = ReadResultsFile(in, global, 1, lbs, 0);
    if (!res) {
        error("error reading results!");
        return 1;
    }

    out = fopen(matname, "wb");
    if (out == NULL) {
        error("res2mat: could not open output file for writing");
        return 1;
    }

    // twoD = 0;
    if (use_file_twoD)
        twoD = res -> twoD;
    else
        twoD = requested_twoD;

//    printf("use_file = %d, req_twoD = %d, file_twoD = %d\n",
//           use_file_twoD, requested_twoD, res -> twoD); 

    if (res -> depth_ref) 
        MatlabVector (&(res -> depth), 1, "depth", out);

    if (!res -> output_map [EULER] && global) {
         error("res2mat: cannot rotate into global coordinates without Euler information");
        exit (1);
    }

    if (res -> num_branch > 0)
        ProcessBranchInfo (res, out);
     
    ProcessStatic(res, out, global, twoD);

    if (!totals)
        RemoveStatic (res);

    if (res -> dynamic) 
        ProcessDynamic(res, out, global, lbs, totals, twoD);
  
    fclose(out);

    res_close(in);
 
    return 0;
}
