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
 * File:        results.c
 *
 * Description: routines to read cable binary results file.  
 *
 ****************************************************************************/

# include <math.h>
# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include "compress.h"
# include "allocate.h"
# include "error.h"
# include "results.h"
# include "rotate.h"
# include "tension.h"

double *vector(nr)
   int   nr;
{
   int	    i; 
   double  *ptr;

   ptr = (double *) malloc(sizeof(double) * nr);
   if (ptr == NULL) {
      Fatal("could not allocate memory (vector %d)", nr);
      return NULL;
   } 

   for (i = 0 ; i < nr ; i++)
      ptr [i] = 0.0;

   return ptr;
}

double **array(nr, nc)
   int   nr;
   int   nc;
{
   int	    i, j;
   double **ptr;

   ptr = (double **) malloc(sizeof(double *) * nr);
   if (ptr == NULL) {
      Fatal("could not allocate memory (array %d x %d)", nr, nc);
      return NULL;
   } 
   for (i = 0 ; i < nr ; i++) {
      ptr [i] = (double *) malloc(sizeof(double) * nc);
      if (ptr[i] == NULL) {
          Fatal("could not allocate memory row %d (array %d x %d)", i, nr, nc);
          return NULL;
      }
   }

   for (i = 0 ; i < nr ; i++)
      for (j = 0 ; j < nc ; j++)
         ptr [i][j] = 0.0;

   return ptr;
}
   
void 
free_array(double **ptr, int nr, int nc)
{
   int	    i;

   if (!ptr)
      return;

   for (i = 0 ; i < nr ; i++)
      free (ptr [i]);

   free (ptr);

   return;
}

int
ReadResultSnapshot(Result *res, int frame,
                   int totals, int lbs, int ft)
{
    int i, k, cnt;
    double d[64];
    double force, length;

    if (res -> snapTell == NULL || frame >= res -> nsteps || res -> snapTell[frame] == 0)
        return 1;

    if (lbs)
      force = 4.4482216;
    else
      force = 1.0;

    if (ft)
      length = 0.3048;
    else
      length = 1.0;

    res_seek(res -> in, res -> snapTell[frame], SEEK_SET);
    for (i = 0 ; i < res -> npoints ; i++) {
       cnt = res_read (d, sizeof(double), res -> numvars, res -> in);
# ifdef ZLIB
       cnt /= sizeof(double); /* gzread returns num bytes, not num items */
# endif
        if (cnt != res -> numvars) {
           printf("error reading snapshot data, cnt=%d, nv=%d, np=%d\n", cnt, res -> numvars, res -> npoints);
 	       return 1; 
        }
        k = 0;
        if (res -> output_map [DISPLACEMENT]) {
           res -> y [i] = d [k]/length - (1 - totals)*res -> y_st [i]; k++;
           res -> x [i] = d [k]/length - (1 - totals)*res -> x_st [i]; k++;
           res -> z [i] = d [k]/length - (1 - totals)*res -> z_st [i]; k++; 
        }
        if (res -> output_map [VELOCITY]) {
           res -> velocity [0][i] = d [k]/length; k++;
           res -> velocity [1][i] = d [k]/length; k++;
           res -> velocity [2][i] = d [k]/length; k++;
        }
        if (res -> output_map [FORCE]) {
           res -> force [0][i] = d [k]/force +
                                      totals*res -> force_st [0][i]; k++;
           res -> force [1][i] = d [k]/force +
                                      totals*res -> force_st [1][i]; k++;
           res -> force [2][i] = d [k]/force +
                                      totals*res -> force_st [2][i]; k++;
        }
        if (res -> output_map [MOMENT]) {
           res -> moment [0][i] = d [k]/length/force +
                                  totals*res -> moment_st [0][i]; k++;
           res -> moment [1][i] = d [k]/length/force +
                                  totals*res -> moment_st [1][i]; k++;
           res -> moment [2][i] = d [k]/length/force +
                                  totals*res -> moment_st [2][i]; k++;
        }
        if (res -> output_map [EULER]) {
           res -> beta [0][i] = d [k] 
	            - (1 - totals)*res -> beta_st [0][i]; k++;
           res -> beta [1][i] = d [k] 
	            - (1 - totals)*res -> beta_st [1][i]; k++;
           res -> beta [2][i] = d [k] 
	            - (1 - totals)*res -> beta_st [2][i]; k++;
           res -> beta [3][i] = d [k] 
                - (1 - totals)*res -> beta_st [3][i]; k++;
        }

    } /* end for loop over all nodes */

    return 0;
}

static int 
ProcessDynamic (Result *res, int num_nodes, 
                ResFile in, int totals, int lbs, int ft)
{
   int		 i, j;
   int		 k;
   int		 cnt;
   double	 start, end, duration;
   double	 dummy;
   char		 magic;
   int		 snap_j, samp_j, buoy_j, seg_j, ext_j;
   int		 node_i;
   int       m[2];
   double	*d;
   double	 x [6];
   double        length;
   double        force;

   if (lbs)
      force = 4.4482216;
   else
      force = 1.0;

   if (ft)
      length = 0.3048;
   else
      length = 1.0;

   res -> totals = totals;

   if (res -> dynamic & 16) 
       res_read(&start, sizeof(double), 1, in);
   else
       start = 0;

   res_read(&end, sizeof(double), 1, in);
   
   duration = end - start;
   res -> t_start = start;

   res_read(&dummy, sizeof(double), 1, in);	 // solution timestep dt

   res_read(&(res -> sample_dt), sizeof(double), 1, in);
   res_read(&(res -> snap_dt), sizeof(double), 1, in);
   
   if (res -> dynamic & 0x02) {
      res_read(&(res -> buoy_dt), sizeof(double), 1, in);
   }
   else
      res -> buoy_dt = 0.0;

   if (res -> dynamic & 0x04) {
      res_read(&(res -> seg_dt), sizeof(double), 1, in);
      res_read(&(res -> nseg), sizeof(int), 1, in);
   }
   else {
      res -> seg_dt = 0.0;
      res -> nseg = 0;
      printf("no segment info\n"); fflush(stdout);
   }

   if (res -> dynamic & 0x08) {
      res_read(&(res -> ext_dt), sizeof(double), 1, in);
      res_read(&(res -> n_ext), sizeof(int), 1, in);
      printf("next = %d, ext_dt = %g\n", res -> n_ext, res -> ext_dt);
   }
   else {
      res -> ext_dt = 0.0;
      res -> n_ext = 0;
      printf("no external force info\n"); fflush(stdout);
   }

   printf("start = %f\n", start);
   printf("end   = %f\n", end);
   printf("dt=%f, snap=%f, seg=%f, ext=%f\n", res -> sample_dt, res -> snap_dt, res -> seg_dt, res -> ext_dt);
          

   res_read(&(res -> num_output_nodes), sizeof(int), 1, in);

   if (res -> num_output_nodes) {
      res -> output_nodes = (int *) malloc(sizeof(int) * res -> num_output_nodes);
      if (res -> output_nodes == NULL) {
          Fatal("could not allocate output nodes");
      }

      res_read(res -> output_nodes, sizeof(int), res -> num_output_nodes, in);
   }
   else
      res -> output_nodes = NULL;
 
   res -> numvars = 0;
   if (res -> output_map [DISPLACEMENT])
      res -> numvars += 3; 
   if (res -> output_map [VELOCITY])
      res -> numvars += 3; 
   if (res -> output_map [FORCE])
      res -> numvars += 3; 
   if (res -> output_map [MOMENT])
      res -> numvars += 3; 
   if (res -> output_map [EULER])
      res -> numvars += 4; 

   if (res -> sample_dt > 0.0 && res -> num_output_nodes > 0) {
      res -> nsamples = (int) ((duration + res -> sample_dt/2.0) / res -> sample_dt) + 1;

      if (res -> output_map [DISPLACEMENT]) {
         res -> x_t = array(res -> num_output_nodes, res -> nsamples);
         res -> y_t = array(res -> num_output_nodes, res -> nsamples);
         res -> z_t = array(res -> num_output_nodes, res -> nsamples);
      }
      if (res -> output_map [VELOCITY]) {
         for (i = 0 ; i < 3 ; i++) {
            res -> velocity_t [i] = array(res -> num_output_nodes, res -> nsamples);
            res -> velocity [i] = vector(num_nodes);
         }

      }
      if (res -> output_map [FORCE]) {
         for (i = 0 ; i < 3 ; i++)
            res -> force_t [i] = array(res -> num_output_nodes, res -> nsamples);
      }
      if (res -> output_map [MOMENT]) {
         for (i = 0 ; i < 3 ; i++)
            res -> moment_t [i] = array(res -> num_output_nodes, res -> nsamples);
      }
      if (res -> output_map [EULER]) {
         for (i = 0 ; i < 4 ; i++)
            res -> beta_t [i] = array(res -> num_output_nodes, res -> nsamples);
      }
   } 
   else {
      res -> nsamples = 0;
      res -> sample_dt = 0.0;
   }

   if (res -> snap_dt > 0.0) {
      res -> nsteps = (int) ((duration + res -> snap_dt/2.0) / res -> snap_dt) + 1;
      res -> snapTell = (unsigned long *) malloc(sizeof(unsigned long) * res -> nsteps);

      for (i = 0 ; i < res -> nsteps ; i++)
          res -> snapTell[i] = 0;

      if (res -> snapTell == NULL) 
          Fatal("could not allocate tell points for snapshots");
   }
   else {
      res -> nsteps = 0;
      res -> snap_dt = 0.0;
      res -> snapTell = NULL;
      res -> snapTell = (unsigned long *) malloc(sizeof(unsigned long) * res -> nsteps);
   }

   if (res -> seg_dt > 0.0 && res -> nseg) {
      res -> nseg_steps = (int) ((duration + res -> seg_dt/2.0) / res -> seg_dt) + 1;
       
      res -> seg_top_pay = array(res -> nseg_steps, res -> nseg); 
      res -> seg_bot_pay = array(res -> nseg_steps, res -> nseg); 
      res -> seg_stretched = array(res -> nseg_steps, res -> nseg); 
      res -> seg_unstretched = array(res -> nseg_steps, res -> nseg); 
      res -> seg_top_spooled = array(res -> nseg_steps, res -> nseg); 
      res -> seg_bot_spooled = array(res -> nseg_steps, res -> nseg); 
      res -> seg_first = array(res -> nseg_steps, res -> nseg); 
      res -> seg_last = array(res -> nseg_steps, res -> nseg); 
   }
   if (res -> ext_dt && res -> n_ext > 0) {
      res -> n_ext_steps = (int) ((duration + res -> ext_dt/2.0) / res -> ext_dt) + 1;
       for (i = 0 ; i < 3 ; i++)  {
           res -> ext_t[i] = array(res -> n_ext_steps, res -> n_ext);
       }
   }
   if (res -> buoy_dt > 0.0) {
      res -> nbuoy = (int) ((duration + res -> buoy_dt/2.0) / res -> buoy_dt) + 1;
      for (i = 0 ; i < 6 ; i++)
         res -> buoy [i] = vector(res -> nbuoy);
   }
      
   snap_j = 0;
   samp_j = 0;
   buoy_j = 0;
   seg_j  = 0;
   ext_j  = 0;
   d = vector(res -> numvars);

   while (!res_eof(in)) {
      res_read (&magic, sizeof(char), 1, in);
      if (res_eof (in))
         break;
      
      if (magic == 't') {
         for (i = 0 ; i < res -> num_output_nodes ; i++)  {
            cnt = res_read (d, sizeof(double), res -> numvars, in);
# ifdef ZLIB
            cnt /= sizeof(double); /* gzread returns num bytes, not num items */
# endif
            if (cnt != res -> numvars) {
                error ("error reading time series data");
                return 1; 
            }


            node_i = res -> output_nodes [i] - 1;
            
            k = 0;
            if (res -> output_map [DISPLACEMENT]) {
               res -> y_t [i][samp_j] = d [k]/length 
					- (1 - totals)*res -> y_st [node_i]; k++;
               res -> x_t [i][samp_j] = d [k]/length 
					- (1 - totals)*res -> x_st [node_i]; k++;
               res -> z_t [i][samp_j] = d [k]/length 
					- (1 - totals)*res -> z_st [node_i]; k++;
            }
            if (res -> output_map [VELOCITY]) {
               res -> velocity_t [0][i][samp_j] = d [k]/length; k++;
               res -> velocity_t [1][i][samp_j] = d [k]/length; k++;
               res -> velocity_t [2][i][samp_j] = d [k]/length; k++;
            }
            if (res -> output_map [FORCE]) {
               res -> force_t [0][i][samp_j] = d [k]/force + 
                                          totals*res -> force_st [0][node_i]; k++;
               res -> force_t [1][i][samp_j] = d [k]/force +
                                          totals*res -> force_st [1][node_i]; k++;
               res -> force_t [2][i][samp_j] = d [k]/force +
                                          totals*res -> force_st [2][node_i]; k++;

            }
            if (res -> output_map [MOMENT]) {
               res -> moment_t [0][i][samp_j] = d [k]/length/force +
                                        totals*res -> moment_st [0][node_i]; k++;
               res -> moment_t [1][i][samp_j] = d [k]/length/force +
                                        totals*res -> moment_st [1][node_i]; k++;
               res -> moment_t [2][i][samp_j] = d [k]/length/force +
                                        totals*res -> moment_st [2][node_i]; k++;
            }
            if (res -> output_map [EULER]) {
               res -> beta_t [0][i][samp_j] = d [k] 
                                    - (1 - totals)*res -> beta_st [0][node_i]; k++;
               res -> beta_t [1][i][samp_j] = d [k] 
                                    - (1 - totals)*res -> beta_st [1][node_i]; k++;
               res -> beta_t [2][i][samp_j] = d [k] 
                                    - (1 - totals)*res -> beta_st [2][node_i]; k++;
               res -> beta_t [3][i][samp_j] = d [k] 
                                    - (1 - totals)*res -> beta_st [3][node_i]; k++;
            }

         } /* end for loop over output nodes */
         samp_j ++;
      }

      else if (magic == 's') {
         res -> snapTell[snap_j ++] = res_tell(in);

         for (i = 0 ; i < num_nodes ; i++) {
            cnt = res_read (d, sizeof(double), res -> numvars, in);
# ifdef ZLIB
	        cnt /= sizeof(double); /* gzread returns num bytes, not num items */
# endif
            if (cnt != res -> numvars) {
               error ("error reading snapshot data, i=%d (%d/%d), snap=%d, tell=%ld", 
                      i+1, cnt, res -> numvars, snap_j, gztell(in));
     	       return 1; 
            }
         } /* end for loop over all nodes */

      } /* end if entry was a snapshot */

      else if (magic == 'b') {
         cnt = res_read (x, sizeof(double), 6, in);
# ifdef ZLIB
         cnt /= sizeof(double); /* gzread returns num bytes, not num items */
# endif
         if (cnt != 6*sizeof(double)) {
            error ("error reading buoy motion data, cnt = %d", cnt);
            return 1; 
         }

         for (i = 0 ; i < 6 ; i++)
            res -> buoy [i][buoy_j] = x [i];

         buoy_j ++;

      } /* end if entry was a buoy motion record */
      else if (magic == 'e' && res -> nseg && res -> nseg_steps) {
         for (i = 0 ; i < res -> nseg ; i++) {
            cnt = res_read(x, sizeof(double), 6, in);
            if (cnt != 6*sizeof(double)) {
               error ("error reading segment data, cnt = %d", cnt);
               return 1; 
            }
            res -> seg_unstretched[seg_j][i] = x[0];
            res -> seg_stretched[seg_j][i]   = x[1];
            res -> seg_bot_pay[seg_j][i]     = x[2];
            res -> seg_top_pay[seg_j][i]     = x[3];
            res -> seg_bot_spooled[seg_j][i] = x[4];
            res -> seg_top_spooled[seg_j][i] = x[5];
            cnt = res_read(m, sizeof(int), 2, in);
            if (cnt != 2*sizeof(int)) {
               error ("error reading segment data, cnt = %d", cnt);
               return 1; 
            }
            res -> seg_first[seg_j][i] = m[0];
            res -> seg_last[seg_j][i] = m[1];
         }
         seg_j ++;
      }
      else if (magic == 'x' && res -> n_ext && res -> n_ext_steps) {
         for (i = 0 ; i < res -> n_ext ; i++) {
            cnt = res_read(x, sizeof(double), 3, in);
            if (cnt != 3*sizeof(double)) {
               error ("error reading external force data, cnt = %d", cnt);
               return 1; 
            }
            res -> ext_t[0][ext_j][i] = x[0];
            res -> ext_t[1][ext_j][i] = x[1];
            res -> ext_t[2][ext_j][i] = x[2];
         }
         ext_j ++;
      }
      else {
      }
   } /* end while not at end of file */

   return 0;
}

static void
AllocateRes(Result *res, int num_nodes)
{
   int i, numvars;

   res -> s = vector(num_nodes);
   numvars = 0;

   if (res -> output_map [DISPLACEMENT]) {
      res -> x = vector(num_nodes);
      res -> y = vector(num_nodes);
      res -> z = vector(num_nodes);

      res -> x_st = vector(num_nodes);
      res -> y_st = vector(num_nodes);
      res -> z_st = vector(num_nodes);

      numvars += 3;
   }
   if (res -> output_map [VELOCITY]) {
      for (i = 0 ; i < 3 ; i++) {
         res -> velocity [i] = vector(num_nodes);
      }
   }
   if (res -> output_map [FORCE]) {
      for (i = 0 ; i < 3 ; i++) {
         res -> force [i] = vector(num_nodes);
         res -> force_st [i] = vector(num_nodes);
      }

      numvars += 3;
   }
   if (res -> output_map [MOMENT]) {
      for (i = 0 ; i < 3 ; i++) {
         res -> moment [i] = vector(num_nodes);
         res -> moment_st [i] = vector(num_nodes);
      }

      numvars += 3;
   }
   if (res -> output_map [EULER]) {
      for (i = 0 ; i < 4 ; i++) {
         res -> beta [i] = vector(num_nodes); 
         res -> beta_st [i] = vector(num_nodes);
      }

      numvars += 4;
   }

    res -> numvars = numvars;

    return; 
}

static int 
ProcessStatic (Result *res, int num_nodes, ResFile in, int lbs, int ft)
{
   int		 k;
   int		 j;
   double	*d;
   double	 length;
   double	 force;

   if (lbs)
      force = 4.4482216;
   else
      force = 1.0;

   if (ft)
      length = 0.3048;
   else
      length = 1.0;

   AllocateRes(res, num_nodes);

   d = vector(res -> numvars + 1);

   for (j = 0 ; j < num_nodes ; j++) {
      res_read(d, sizeof(double), res -> numvars + 1, in);
      
      k = 0;
      res -> s [j] = d [k]/length; k++;

      if (res -> output_map [DISPLACEMENT]) {
         res -> y_st [j] = res -> y [j] = d [k]/length; k++;
         res -> x_st [j] = res -> x [j] = d [k]/length; k++;
         res -> z_st [j] = res -> z [j] = d [k]/length; k++;
      }
      if (res -> output_map [FORCE]) {
         res -> force_st [0][j] = res -> force [0][j] = d [k]/force; k++;
         res -> force_st [1][j] = res -> force [1][j] = d [k]/force; k++;
         res -> force_st [2][j] = res -> force [2][j] = d [k]/force; k++;
      }
      if (res -> output_map [MOMENT]) {
         res -> moment_st [0][j] = res -> moment [0][j] = d [k]/force/length; k++;
         res -> moment_st [1][j] = res -> moment [1][j] = d [k]/force/length; k++;
         res -> moment_st [2][j] = res -> moment [2][j] = d [k]/force/length; k++;
      }
      if (res -> output_map [EULER]) {
         res -> beta_st [0][j] = res -> beta [0][j] = d [k]; k++;
         res -> beta_st [1][j] = res -> beta [1][j] = d [k]; k++;
         res -> beta_st [2][j] = res -> beta [2][j] = d [k]; k++;
         res -> beta_st [3][j] = res -> beta [3][j] = d [k]; k++;
      }
   }

   return 0;
}


void 
SetSnapGlobal (Result *res)
{
   double      b0, b1, b2, b3;
   int         i, k;
   double      G [4];

   if (res->beta [0] == NULL)
      return;

   for (i = 0 ; i < res->npoints ; i++) {

      b0 = res->beta_st [0][i];
      b1 = res->beta_st [1][i];
      b2 = res->beta_st [2][i];
      b3 = res->beta_st [3][i];

      if (res->force_st [0] != NULL) {
         RotateToFixed(res->force_st [0][i], 
                res->force_st [1][i], 
                res->force_st [2][i], 
                b0, b1, b2, b3,
                &(G [0]), &(G [1]), &(G [2]), res->twoD);

         for (k = 0 ; k < 3 ; k++)
            res->force_st [k][i] = G [k];
      }
      if (res->moment_st [0] != NULL && !res->twoD) {
         RotateToFixed(res->moment_st [0][i], 
                res->moment_st [1][i], 
                res->moment_st [2][i], 
                b0, b1, b2, b3,
                &(G [0]), &(G [1]), &(G [2]), res->twoD);
         for (k = 0 ; k < 3 ; k++)
            res->moment_st [k][i] = G [k];
      }

      b0 = res->beta [0][i] + (1 - res->totals)*res->beta_st [0][i];
      b1 = res->beta [1][i] + (1 - res->totals)*res->beta_st [1][i];
      b2 = res->beta [2][i] + (1 - res->totals)*res->beta_st [2][i];
      b3 = res->beta [3][i] + (1 - res->totals)*res->beta_st [3][i];

         if (res->velocity [0] != NULL) {
            RotateToFixed(res->velocity [0][i], 
                   res->velocity [1][i], 
                   res->velocity [2][i], 
                   b0, b1, b2, b3,
                   &(G [0]), &(G [1]), &(G [2]), res->twoD);

            for (k = 0 ; k < 3 ; k++)
               res->velocity [k][i] = G [k];
         }
         if (res->force [0] != NULL) {
            RotateToFixed(res->force [0][i], 
                           res->force [1][i], 
                           res->force [2][i], 
                           b0, b1, b2, b3,
                           &(G [0]), &(G [1]), &(G [2]), res->twoD);

            for (k = 0 ; k < 3 ; k++)
               res->force [k][i] = G [k];
         }
         if (res->moment [0] != NULL && !res->twoD) {
            RotateToFixed(res->moment [0][i], 
                           res->moment [1][i], 
                           res->moment [2][i], 
                           b0, b1, b2, b3,
                           &(G [0]), &(G [1]), &(G [2]), res->twoD);
            for (k = 0 ; k < 3 ; k++)
               res->moment [k][i] = G [k];
         }
      }
}
 
static void 
SetSampleGlobal (Result *res)
{
   double      b0, b1, b2, b3;
   int	       nd;
   int         i, j, k;
   double      G [4];


   for (i = 0 ; i < res->nsamples ; i++) {
      for (j = 0 ; j < res->num_output_nodes ; j++)  {
         nd = res->output_nodes [j] - 1;
         b0 = res->beta_t [0][j][i] + (1 - res->totals)*res->beta_st [0][nd];
         b1 = res->beta_t [1][j][i] + (1 - res->totals)*res->beta_st [1][nd];
         b2 = res->beta_t [2][j][i] + (1 - res->totals)*res->beta_st [2][nd];
         b3 = res->beta_t [3][j][i] + (1 - res->totals)*res->beta_st [3][nd];

         if (res->velocity_t [0] != NULL) {
            RotateToFixed(res->velocity_t [0][j][i], 
                           res->velocity_t [1][j][i], 
                           res->velocity_t [2][j][i], 
                           b0, b1, b2, b3,
                           &(G [0]), &(G [1]), &(G [2]), res->twoD);

            for (k = 0 ; k < 3 ; k++)
               res->velocity_t [k][j][i] = G [k];
         }
         if (res->force_t [0] != NULL) {
            RotateToFixed(res->force_t [0][j][i], 
                           res->force_t [1][j][i], 
                           res->force_t [2][j][i], 
                           b0, b1, b2, b3,
                           &(G [0]), &(G [1]), &(G [2]), res->twoD);

            for (k = 0 ; k < 3 ; k++)
               res->force_t [k][j][i] = G [k];
         }
         if (res->moment_t [0] != NULL) {
            RotateToFixed(res->moment_t [0][j][i], 
                           res->moment_t [1][j][i], 
                           res->moment_t [2][j][i], 
                           b0, b1, b2, b3,
                           &(G [0]), &(G [1]), &(G [2]), res->twoD);

            for (k = 0 ; k < 3 ; k++)
               res->moment_t [k][j][i] = G [k];
         }
      }
   }

   return;
}

Result *
ReadResultsFile (ResFile in, int global, int totals, int lbs, int ft)
{
   int      i, j;
   int	    read_branch_info;
   int      read_input_text;
   char	    magic [6];
   char	   *title;
   char	    b;
   int	    num_nodes;
   int	    status;
   Result  *res;

   res = (Result *) malloc(sizeof(Result));
   if (res == NULL) 
       Fatal("could not allocate memory for result"); 

   res -> in = in;

	/*
	 * initialize all of the res -> lt vectors so we will 
	 * know if they were even in the input file	
	 */

   res -> s = NULL;

   res -> x = NULL;
   res -> y = NULL;
   res -> z = NULL;

   res -> x_t = NULL;
   res -> y_t = NULL;
   res -> z_t = NULL;

   res -> x_st = NULL;
   res -> y_st = NULL;
   res -> z_st = NULL;

   res -> input = NULL;

   for (i = 0 ; i < 3 ; i++) {
      res -> velocity [i] = NULL;
      res -> force [i] = NULL;
      res -> moment [i] = NULL;
      res -> beta [i] = NULL;

      res -> velocity_t [i] = NULL;
      res -> force_t [i] = NULL;
      res -> moment_t [i] = NULL;
      res -> beta_t [i] = NULL;

      res -> force_st [i] = NULL;
      res -> moment_st [i] = NULL;
      res -> beta_st [i] = NULL;
   }
   res -> beta [3] = NULL;
   res -> beta_t [3] = NULL;

   res -> ext_t[0] = res -> ext_t[1] = res -> ext_t[2] = NULL;

   res_read(magic, sizeof(char), 6, in);
   if (strncmp(magic, "cabres", 6) != 0) {
      error ("incorrect input file format, magic %6s", magic);
      return NULL;
   }

   res_read(&b, sizeof(char), 1, in);
   if (!b) {
      error ("static solution not present in file");
      return NULL; 
   }

   if (b & 2)
      res -> depth_ref = 1;
   else
      res -> depth_ref = 0;

   if (b & 4)
      res -> type = drawHorizontal;
   else if (b & 8)
      res -> type = drawTowing;
   else
      res -> type = drawStandard;
  
   if (b & 16)
      res -> twoD = 1;
   else
      res -> twoD = 0;

   if (b & 32)
      read_branch_info = 1;
   else
      read_branch_info = 0;

   if (b & 64)
      read_input_text = 1;
   else
      read_input_text = 0;

   res_read(&b, sizeof(char), 1, in);
   res -> dynamic = b;

   if (read_input_text) {
      // i = strlen() + 1 in solver/output.c for terminating null
      res_read(&i, sizeof(int), 1, in);
      res -> input = (char *) malloc(sizeof(char) * i);
      res_read(res -> input, sizeof(char), i, in);
   }

   res_read(&i, sizeof(int), 1, in);
   num_nodes = i;

   res_read(&i, sizeof(int), 1, in);
   title = (char *) malloc(sizeof(char)*i);
   res_read(title, sizeof(char), i, in);
   res -> title = title;

   if (res -> depth_ref)
      res_read(&(res -> depth), sizeof(double), 1, in);
   else
      res -> depth = 0.0;

   if (ft)
      res -> depth /= 0.3048;

   for (j = 1 ; j < MAXOUTPUT ; j++) {
      res_read(&b, sizeof(char), 1, in);
      res -> output_map [j] = b;
   }

   if (read_branch_info) {
      res_read(&i, sizeof(int), 1, in);
      res -> num_branch = i;
   
      res -> branch_starts = (int *) malloc(sizeof(int) * res -> num_branch);
      for (j = 0 ; j < res -> num_branch ; j++) {
         res_read(&i, sizeof(int), 1, in);
         res -> branch_starts [j] = i;
      }
   }
   else 
      res -> num_branch = 0; 

   res -> npoints = num_nodes;
  
   printf("nnodes = %d\n", num_nodes);

   if (res -> num_branch)
      res -> main_points = res -> branch_starts [0] - 1;
   else
      res -> main_points = num_nodes;
 
   ProcessStatic (res, num_nodes, in, lbs, ft);
   
   if (res -> dynamic)  {
      status = ProcessDynamic (res, num_nodes, in, totals, lbs, ft);
      if (status)
         return NULL;
   }
   else {
      res -> nsteps = 1;
      res -> snap_dt = res -> sample_dt = 0.0;
      res -> output_nodes = NULL;       
      res -> num_output_nodes = 0;
   }

   if (global) {
      res -> global = 1;
      SetSampleGlobal (res);
   }

   return res;
}

void
ResultsToProblem(Result *res, Problem *p, Environment *e, int twoD)
{
   int      i, j;
   Node     node, a;

   j = 0;
   for (i = 0 ; i <= p -> num_branch ; i++) {
      if (i == 0)
         a = p -> node[1];
      else 
         a = p -> branch[i] -> first;
 
      while(a) {
         if (a -> active)
            node = a;
         else if (a -> number  <= a -> segment -> first_active -> number)
            node = a -> segment -> first_active; 
         else if (a -> number >= a -> segment -> last_active -> number)
            node = a -> segment -> last_active;
         else
            node = a; // for -Wall
 
         node -> s = res -> s [j];

         node -> x = res -> y_st [j];
         node -> y = res -> x_st [j];
         node -> z = res -> z_st [j];

         if (twoD) {
             if (res -> output_map[FORCE]) {
                node -> Y[1] = Strain(res -> force_st[0][j], node -> material);
                node -> Y[2] = res -> force_st [1][j];
             }
             else {
                node -> Y[1] = 0.0;
                node -> Y[2] = 0.0;
             }
             if (res -> output_map[EULER])
                node -> Y[3] = res -> beta_st[0][j];
             else
                node -> Y[3] = 0.0;
             
             if (res -> output_map[MOMENT])
                node -> Y[4] = res -> moment_st [2][j];
             else
                node -> Y[4] = 0.0;
         }
         else {
             if (res -> output_map[FORCE]) {
                 node -> Y[1] = Strain(res -> force_st[0][j], node -> material);
                 node -> Y[2] = res -> force_st [1][j];
                 node -> Y[3] = res -> force_st [2][j];
             }
             else {
                 node -> Y[1] = node -> Y[2] = node -> Y[3] = 0.0;
             }
             if (res -> output_map[EULER]) {
                 node -> Y[4] = res -> beta_st [0][j];
                 node -> Y[5] = res -> beta_st [1][j];
                 node -> Y[6] = res -> beta_st [2][j];
                 node -> Y[7] = res -> beta_st [3][j];
             }
             else {
                 node -> Y[4] = node -> Y[5] = node -> Y[6] = node -> Y[7] = 0.0;
             }
             if (res -> output_map[MOMENT]) {
                 node -> Y[8] = res -> moment_st [0][j];
                 node -> Y[9] = res -> moment_st [1][j];
                 node -> Y[10] = res -> moment_st [2][j];
             }
             else {
                 node -> Y[8] = node -> Y[9] = node -> Y[10] = 0.0;
             }
         }

         a = a -> next;
         j ++;
      }
   }  

   return;
}



static int 
SumStaticDynamic(Result *res, int dir)
{
   int   i, j, k;
   int   nd;   

      for (i = 0 ; i < res -> npoints ; i++) {
         
            if (res -> output_map [DISPLACEMENT]) {
               res -> x [i] += dir*res -> x_st [i]; 
               res -> y [i] += dir*res -> y_st [i]; 
               res -> z [i] += dir*res -> z_st [i]; 
            }

            for (k = 0 ; k < 3 ; k++) {
               if (res -> output_map [FORCE])
                  res -> force [k][i] += dir*res -> force_st [k][i];

               if (res -> output_map [MOMENT])
                  res -> moment [k][i] += dir*res -> moment_st [k][i];

               if (res -> output_map [EULER])
                  res -> beta [k][i] += dir*res -> beta_st [k][i];
            }
            if (res -> output_map [EULER])
               res -> beta [3][i] += dir*res -> beta_st [3][i];
         }

      for (j = 0 ; j < res -> num_output_nodes ; j++)  {

         nd = res -> output_nodes [j] - 1;

         for (i = 0 ; i < res -> nsamples ; i++) {
             
            if (res -> output_map [DISPLACEMENT]) {
               res -> x_t [j][i] += dir*res -> x_st [nd]; 
               res -> y_t [j][i] += dir*res -> y_st [nd]; 
               res -> z_t [j][i] += dir*res -> z_st [nd]; 
            }

            for (k = 0 ; k < 3 ; k++) {
               if (res -> output_map [FORCE])
                  res -> force_t [k][j][i] += dir*res -> force_st [k][nd];
        
               if (res -> output_map [MOMENT]) 
                  res -> moment_t [k][j][i] += dir*res -> moment_st [k][nd];

               if (res -> output_map [EULER]) 
                  res -> beta_t [k][j][i] += dir*res -> beta_st [k][nd];
            }
            if (res -> output_map [EULER]) 
               res -> beta_t [3][j][i] += dir*res -> beta_st [3][nd];
         }
      }

   return 0;
}

int 
RemoveStatic(Result *res )
{
   if (!res -> totals)
      return 1;

   res -> totals = 0;

   return SumStaticDynamic(res, -1);
}

int 
AddStatic (Result *res )
{
   if (res -> totals)
      return 1;

   res -> totals = 1;

   return SumStaticDynamic(res, 1);
}

void
TensionToGlobal(Node node, int twoD, double *Fx, double *Fy, double *Fz)
{
    double  T;

    T = NodeTension(node);
    if (twoD)
        RotateToFixed(T, node -> Y[2], 0, node -> Y[3], 0, 0, 0, Fx, Fy, Fz, twoD);  
    else
        RotateToFixed(T, node -> Y[2], node->Y[3], 
                      node -> Y[4], node -> Y[5], node -> Y[6], node -> Y[7],
                      Fx, Fy, Fz, twoD);  
}

char
CheckMagic(char *name)
{
   ResFile in;
   char    magic[6];

   in = res_open(name, "rb");
   if (in == NULL)
       return -1;

   res_read(magic, sizeof(char), 6, in);
   res_close(in);
   if (strncmp(magic, "cabres", 6) == 0) {
      return 1;
   }

   return 0; 
}
