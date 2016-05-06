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
 * File:	output.c
 *
 * Description:	routines to write the binary results files
 *
 ***************************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <malloc.h>
# include <string.h>
# include "math.h"
# include "compress.h"
# include "problem.h"
# include "tension.h"
# include "output.h"
# include "error.h"
# include "results.h"
# include "segments.h"

static ResFile flush_fp;

extern Problem *problem;             
extern Environment *environment;

int InitializeResultsFile (fp, num_nodes, title, output_map, decimate, dynamic, twoD, input)
   ResFile	 fp;
   int		 num_nodes;
   char		*title;
   int		*output_map;
   int		 decimate;
   int		 dynamic;
   int		 twoD;
   char     *input;
{
   static char	magic [6] = "cabres";
   char	 	b;
   int		i, j;
   double	x;
   int		num_out;

   if (res_write(magic, sizeof(char), 6, fp) == 0) {
        error("write error");
        return 1;
   }

   b = 1;

   if (environment -> depth)
      b = b | 2;

   if (problem -> type == Horizontal)
      b = b | 4;
   else if (problem -> type == Towing || problem -> type == Drifter)
      b = b | 8;

   if (twoD)
      b = b | 16;

   if (problem -> num_branch > 0)
      b = b | 32;

   if (input) 
      b = b | 64;

   res_write(&b, sizeof(char), 1, fp);

   b = dynamic;
   res_write(&b, sizeof(char), 1, fp);

   if (input) {
      i = strlen(input) + 1; // add one for the terminating null
      res_write(&i, sizeof(int), 1, fp);
      res_write(input, sizeof(char), i, fp);
   }

   num_out = (int) ceil((double) num_nodes / (double) decimate);
   res_write(&num_out, sizeof(int), 1, fp);

   i = strlen(title) + 1;
   res_write(&i, sizeof(int), 1, fp);
   res_write(title, sizeof(char), i, fp);

   if (environment -> depth) {
      x = environment -> depth;
      res_write (&x, sizeof(double), 1, fp);
   }

   for (j = 1 ; j < MAXOUTPUT ; j++) {
      b = output_map [j];
      res_write(&b, sizeof(char), 1, fp);
   }

   if (problem -> num_branch > 0) {

      i = problem -> num_branch;
      res_write(&i, sizeof(int), 1, fp);

      for (j = 1 ; j <= problem -> num_branch ; j++) {
         i = problem -> branch [j] -> first -> output_number;
         
         res_write(&i, sizeof(int), 1, fp);
      }    
   }

   flush_fp = fp;
 
   return 0;
}

void
WriteDynamicExternalForces(Problem *p, ResFile out)
{
    static char     magic = 'x';
    double          d[3];
    int             i;

    if (res_write(&magic, sizeof(char), 1, out) == 0) {
        error("ext write error");
        return;
    }

    d[0] = p -> terminal[1] -> xthrust.value; 
    d[1] = p -> terminal[1] -> ythrust.value; 
    d[2] = p -> terminal[1] -> zthrust.value; 
    res_write (d, sizeof(double), 3, out);
   
    for (i = 1 ; i <= p -> num_segments ; i++) {
        if (p -> segment[i] -> connector)  {
            d[0] = p -> segment[i] -> connector_xthrust.value;
            d[1] = p -> segment[i] -> connector_ythrust.value;
            d[2] = p -> segment[i] -> connector_zthrust.value;
            res_write (d, sizeof(double), 3, out);
        }

        if (p -> segment[i] -> branch && p -> segment[i] -> branch -> terminal) {
            d[0] = p -> segment[i] -> branch -> terminal -> xthrust.value; 
            d[1] = p -> segment[i] -> branch -> terminal -> ythrust.value; 
            d[2] = p -> segment[i] -> branch -> terminal -> zthrust.value; 
        }
    }  

    d[0] = p -> terminal[2] -> xthrust.value; 
    d[1] = p -> terminal[2] -> ythrust.value;
    d[2] = p -> terminal[2] -> zthrust.value;
    res_write (d, sizeof(double), 3, out);
}
   
int WriteStaticSolution (n, num_nodes, out, output_map, decimate, twoD)
   Node		 *n;
   int		  num_nodes;
   ResFile	  out;
   int		 *output_map;
   int		  decimate;
   int		  twoD;
{
   int		 i, k;
   int		 numvars;
   double	*d;
   Node		 node, a;

   if (out == NULL)
      return 0;

   numvars = 1;
   numvars += 3*output_map [MOTION];
   numvars += 3*output_map [FORCE];
   numvars += 3*output_map [MOMENT];
   numvars += 4*output_map [EULER];

   d = (double *) malloc(sizeof(double) * numvars);

   for (i = 0 ; i <= problem -> num_branch ; i++) {
      if (i == 0)
         a = n [1];
      else
         a = problem -> branch[i] -> first;

      while (a) {
         k = 0;

         if (a -> active) 
            node = a;
         else if (a -> number  <= a -> segment -> first_active -> number)
            node = a -> segment -> first_active;
         else if (a -> number >= a -> segment -> last_active -> number)
            node = a -> segment -> last_active;
         else
            node = a; // for -Wall

         // printf("%d %d %f %f %f\n", a -> number, node -> number,
         //       node -> x, node -> y, node -> z);

         d [k] = a -> s; k++;	// should we use node here (most recent active)?

         if (output_map [MOTION]) {
            d [k] = node -> x; k++;
            d [k] = node -> y; k++;
            d [k] = node -> z; k++;
         }
         if (output_map [FORCE]) {
            d [k] = Tension(node -> Y[1], node -> material); k++;
            d [k] = node -> Y[2]; k++;
            d [k] = twoD ? 0.0 : node -> Y[3]; k++;
         }
         if (output_map [MOMENT]) {
            d [k] = twoD ? 0.0 : node -> Y[8]; k++;
            d [k] = twoD ? 0.0 : node -> Y[9]; k++;
            d [k] = (twoD ? node -> Y[4] : node -> Y[10]); k++;
         }
         if (output_map [EULER]) {
            d [k] = twoD ? node -> Y[3] : node -> Y[4]; k++;
            d [k] = twoD ? 0.0 : node -> Y[5]; k++;
            d [k] = twoD ? 0.0 : node -> Y[6]; k++;
            d [k] = twoD ? 0.0 : node -> Y[7]; k++;
         }

         res_write (d, sizeof(double), numvars, out);

         a = a -> next;
      }
   }

   free(d);

   return 0;
}

int WriteDynamicHeader (out, t_start, t_end, dt, sample_dt, snap_dt, seg_dt, buoy_dt, ext_dt, n, on, decimate, node)
   ResFile 	 out;
   double    t_start;
   double	 t_end;
   double	 dt;
   double	 sample_dt;
   double	 snap_dt;
   double    seg_dt;
   double	 buoy_dt;
   double    ext_dt;
   int		 n;
   int		*on;
   int		 decimate;
   Node		*node;
{
   double	 d;
   int		 i, j;
   int       n_ext;

   if (out == NULL)
      return 0;

   d = t_start;
   res_write(&d, sizeof(double), 1, out);
   d = t_end;
   res_write(&d, sizeof(double), 1, out);
   d = dt;
   res_write(&d, sizeof(double), 1, out);
   d = sample_dt;
   res_write(&d, sizeof(double), 1, out);
   d = snap_dt;
   res_write(&d, sizeof(double), 1, out);
   printf("wrote snap_dt = %f\n", snap_dt);

   if (buoy_dt) {
      d = buoy_dt;
      res_write(&d, sizeof(double), 1, out);
      printf("wrote buoy_dt\n");
   }

   if (seg_dt) {
      d = seg_dt;
      res_write(&d, sizeof(double), 1, out);
      i = problem -> num_segments;
      res_write(&i, sizeof(int), 1, out);
      printf("wrote seg_dt = %f, ns = %d\n", seg_dt, problem -> num_segments);
   }
   
   if (ext_dt) {
        d = ext_dt; 
        res_write(&d, sizeof(double), 1, out);

        n_ext = 2;
        for (i = 1 ; i <= problem -> num_segments ; i++) {
            if (problem -> segment[i] -> connector)
                n_ext ++;
            if (problem -> segment[i] -> branch && problem -> segment[i] -> branch -> terminal)
                n_ext ++;
        }
        i = n_ext;
        res_write(&i, sizeof(int), 1, out);
        printf("wrote ext_dt = %f, n_ext = %d\n", ext_dt, n_ext);
   }

   i = n;
   res_write(&i, sizeof(int), 1, out);
   printf("wrote n = %d\n", n);

   for (i = 0 ; i < n ; i++) {
      j = node [on [i]] -> output_number;
      j = (int) ceil((double) j / (double) decimate);

      res_write(&j, sizeof(int), 1, out);
   }

   return 0;
}

int
WriteDynamicSegmentData(Problem *p, ResFile out)
{
    double          ds, ds0;
    static Segment *seg;
    static int      ns = 0;
    static char     magic = 'e';
    double          d[6];
    int             m[2];
    int             i;
    Node            n;

    if (ns == 0) {
        seg = BuildSegmentArray(p, &ns, 0);
    }

    if (res_write(&magic, sizeof(char), 1, out) == 0) {
        error("seg write error");
        return 1;
    }

    for (i = 1 ; i <= ns ; i++) {
        n = seg[i] -> first_active;
        ds = ds0 = 0;
        while(n && n -> active && n -> segment == seg[i]) {
            ds0 += n -> ds;
            ds  += n -> ds*(1.0 + n -> Y[1]);
            n = n -> next_active;
        }
        d[0] = ds0;
        d[1] = ds;
        d[2] = seg[i] -> first_active -> pay;
        d[3] = seg[i] -> last_active -> pay;
        d[4] = seg[i] -> bottom_spooled;
        d[5] = seg[i] -> top_spooled;
        res_write(d, sizeof(double), 6, out);
        m[0] = seg[i] -> first_active -> number;
        m[1] = seg[i] -> last_active -> number;
        res_write(m, sizeof(int), 2, out);
    }

    return 0;

}

int WriteBuoyMotion(out, x)
   ResFile	 out;
   double	*x;
{
   static char	  magic = 'b';

   res_write(&magic, sizeof(char), 1, out);
   res_write (x, sizeof(double), 6, out);

   return 0;
}

int WriteDynamicSnapshot (out, output_map, n, num_nodes, decimate, twoD)
   ResFile 	 out;
   int		*output_map;
   Node		*n;
   int		 num_nodes;
   int		 decimate;
   int		 twoD;
{
   static char	  magic = 's';
   static double *d = NULL;
   static int     numvars;
   int		  i, j, k;
   Node		  node, a;
   
   if (out == NULL)
      return 1;

   if (d == NULL) {
      numvars = 0;
      numvars += 3*output_map [MOTION];
      numvars += 3*output_map [VEL];
      numvars += 3*output_map [FORCE];
      numvars += 3*output_map [MOMENT];
      numvars += 4*output_map [EULER];

      d = (double *) malloc (sizeof(double) * numvars);
   } 
   
   if (res_write(&magic, sizeof(char), 1, out) == 0) {
        error("snapshot write error");
        return 1;
   }

   j = 0;
   for (i = 0 ; i <= problem -> num_branch ; i++) {
      if (i == 0)
         a = n [1];
      else
         a = problem -> branch[i] -> first;
                                                                                
      j = 0;
      while (a) {
         k = 0;
                                                                                
         if (a -> active)
            node = a;
         else if (a -> number  <= a -> segment -> first_active -> number)
            node = a -> segment -> first_active;
         else if (a -> number >= a -> segment -> last_active -> number)
            node = a -> segment -> last_active;
         else
            node = a; // for -Wall

         if (output_map [MOTION]) {
            d [k] = node -> x; k++;
            d [k] = node -> y; k++;
            d [k] = node -> z; k++;
         }
         if (output_map [VEL]) {
            d [k] = twoD ? node -> Y[3] : node -> Y[4]; k++;
            d [k] = twoD ? node -> Y[4] : node -> Y[5]; k++;
            d [k] = twoD ? 0.0 : node -> Y[6]; k++;
         }
         if (output_map [FORCE]) {
            d [k] = Tension(node -> Y[1], node -> material) 
                    - Tension(node -> Ys[1], node -> material); k++;
            d [k] = node -> Y[2] - node -> Ys[2]; k++;
            d [k] = twoD ? 0.0 : node -> Y[3] - node -> Ys[3]; k++;
         }
         if (output_map [MOMENT]) {
            d [k] = twoD ? 0.0 : (node -> Y[11] - node -> Ys[8]); k++;
            d [k] = twoD ? 0.0 : (node -> Y[12] - node -> Ys[9]); k++;
            d [k] = 
             (twoD ? node -> Y[6] - node -> Ys[4] : node -> Y[13] - node -> Ys[10]); 
            k++;
         }
         if (output_map [EULER]) {
            d [k] = twoD ? node -> Y[5] : node -> Y[7]; k++;
            d [k] = twoD ? 0.0 : node -> Y[8]; k++;
            d [k] = twoD ? 0.0 : node -> Y[9]; k++;
            d [k] = twoD ? 0.0 : node -> Y[10]; k++;
         }

         if (res_write (d, sizeof(double), numvars, out) == 0) {
            error("snapshot write error");
            return 1;
         }

         a = a -> next;
         j ++;
      }
   }

   res_flush(out); 

   return 0;
}

int WriteDynamicResult (out, output_map, output_nodes, 
                        num_output_nodes, decimate, 
                        n, num_nodes,  twoD)
   ResFile	 out;
   int		*output_map;
   int		*output_nodes;
   int		 num_output_nodes;
   int		 decimate;
   Node		*n;
   int		 num_nodes;
   int	         twoD;
{
   static char	  magic = 't';
   static double *d = NULL;
   static int     numvars;
   int		  i, j;
   int		  k;
   Node		  node;
   
   if (out == NULL)
      return 1;

   if (d == NULL) {
      numvars = 0;
      numvars += 3*output_map [MOTION];
      numvars += 3*output_map [VEL];
      numvars += 3*output_map [FORCE];
      numvars += 3*output_map [MOMENT];
      numvars += 4*output_map [EULER];

      d = (double *) malloc (sizeof(double) * numvars);
   } 
   
   if (res_write(&magic, sizeof(char), 1, out) == 0) {
        error("sample write error");
        return 1;
   }

   for (i = 1 ; i <= num_output_nodes ; i++) {
      k = 0;

      j = output_nodes [i - 1];

	/*
	 * if this node is not active find the closest node that is
	 */

      if (n [j] -> active) 
         node = n[j];
      else if (j > n[j] -> segment -> last_active -> number)
         node = n[j] -> segment -> last_active;
      else if (j < n[j] -> segment -> first_active -> number)
         node = n[j] -> segment -> first_active;
      else
         node = n[j];

      if (output_map [MOTION]) {
         d [k] = node -> x; k++;
         d [k] = node -> y; k++;
         d [k] = node -> z; k++;
      }
      if (output_map [VEL]) {
         d [k] = twoD ? node -> Y[3] : node -> Y[4]; k++;
         d [k] = twoD ? node -> Y[4] : node -> Y[5]; k++;
         d [k] = twoD ? 0.0 : node -> Y[6]; k++;
      }
      if (output_map [FORCE]) {
         d [k] = Tension(node -> Y[1], node -> material) 
                 - Tension(node -> Ys[1], node -> material); k++;
         d [k] = node -> Y[2] - node -> Ys[2]; k++;
         d [k] = twoD ? 0.0 : node -> Y[3] - node -> Ys[3]; k++;
      }
      if (output_map [MOMENT]) {
         d [k] = twoD ? 0.0 : (node -> Y[11] - node -> Ys[8]); k++;
         d [k] = twoD ? 0.0 : (node -> Y[12] - node -> Ys[9]); k++;
         d [k] = 
           (twoD ? node -> Y[6] - node -> Ys[4] : node -> Y[13] - node -> Ys[10]); 
         k++;
      }
      if (output_map [EULER]) {
         d [k] = twoD ? node -> Y[5] : node -> Y[7]; k++;
         d [k] = twoD ? 0.0 : node -> Y[8]; k++;
         d [k] = twoD ? 0.0 : node -> Y[9]; k++;
         d [k] = twoD ? 0.0 : node -> Y[10]; k++;
      }

      res_write (d, sizeof(double), numvars, out);
   }

   return 0;
}

int FakeDynamicSnapshot (out, n, num_nodes, decimate, twoD)
   ResFile	 out;
   Node		*n;
   int		 num_nodes;
   int	         decimate;
   int		 twoD;
{
   static char	  magic = 's';
   static double *d = NULL;
   static int     numvars;
   int		  i, k;
   Node		  node, a;
   
   if (d == NULL) {
      numvars = 0;
      numvars += 3;
      numvars += 3;
      numvars += 3;
      numvars += 3;
      numvars += 4;

      d = (double *) malloc (sizeof(double) * numvars);
   } 
   
   res_write(&magic, sizeof(char), 1, out);

   for (i = 0 ; i <= problem -> num_branch ; i++) {
      if (i == 0)
         a = n [1];
      else
         a = problem -> branch[i] -> first;
                                                                                
      while (a) {
         k = 0;
                                                                                
         if (a -> active)
            node = a;
         else if (a -> number  <= a -> segment -> first_active -> number)
            node = a -> segment -> first_active;
         else if (a -> number >= a -> segment -> last_active -> number)
            node = a -> segment -> last_active;
         else
            node = a; // shut -Wall up

         if (1) {
            d [k] = node -> x; k++;
            d [k] = node -> y; k++;
            d [k] = node -> z; k++;
         }
         if (1) {
            d [k] = 0.0; k++;
            d [k] = 0.0; k++;
            d [k] = 0.0; k++;
         }
         if (1) {
            d [k] = Tension(node -> Y[1], node -> material); k++;
            d [k] = node -> Y[2]; k++;
            d [k] = twoD ? 0.0 : node -> Y[3]; k++;
         }
         if (1) {
            d [k] = twoD ? 0.0 : node -> Y [8]; k++;
            d [k] = twoD ? 0.0 : node -> Y [9]; k++;
            d [k] = (twoD ? node -> Y [4] : node -> Y [10]); k++;
         }
         if (1) {
            d [k] = twoD ? node -> Y [3] : node -> Y [4]; k++;
            d [k] = twoD ? 0.0 : node -> Y [5]; k++;
            d [k] = twoD ? 0.0 : node -> Y [6]; k++;
            d [k] = twoD ? 0.0 : node -> Y [7]; k++;
         }

         res_write (d, sizeof(double), numvars, out);

         a = a -> next;
      }
   }

   return 0;
}

void FlushOutput ( )
{
   res_flush(flush_fp);
}

void
InitializeProgressFile(Problem *problem, ResFile prog, int ndof)
{
    int    k;

    k = problem -> num_nodes;
    res_write (&k, sizeof(int), 1, prog);
    k = problem -> num_segments;
    res_write (&k, sizeof(int), 1, prog);
    k = problem -> num_nodes*(5*sizeof(int) + sizeof(double)*(29 + 6*ndof))
        + problem -> num_segments*(2*sizeof(int) + 14*sizeof(double));
    res_write (&k, sizeof(int), 1, prog);
    printf("progress record size = %d\n", k);
}

void
WriteProgress(double t, Problem *problem, ResFile prog, int ndof)
{
    int      i;
    double   d;
    int      k;

    d = t; 
    if (res_write (&d, sizeof(double), 1, prog) == 0) {
        error("progress write error");                        // 1d
        return;
    }
 
    for (i = 1 ; i <= problem -> num_nodes ; i++) {
        k = problem -> node[i] -> active_number; 
        res_write(&k, sizeof(int), 1, prog);                        // 1i

        if (problem -> node[i] -> next_active)
            k = problem -> node[i] -> next_active -> number;        
        else 
            k = 0;
        res_write(&k, sizeof(int), 1, prog);                        // 1i

        if (problem -> node[i] -> prev_active)
            k = problem -> node[i] -> prev_active -> number; 
        else
            k = 0;
        res_write(&k, sizeof(int), 1, prog);                        // 1i

        k = problem -> node[i] -> position;
        res_write(&k, sizeof(int), 1, prog);                        // 1i
        d = problem -> node[i] -> ds;
        res_write(&d, sizeof(double), 1, prog);                     // 1d
        k = problem -> node[i] -> active;
        res_write(&k, sizeof(int), 1, prog);                        // 1i

        res_write(&(problem -> node[i] -> Y[1]), sizeof(double), ndof, prog);   //6*nd*d
        res_write(&(problem -> node[i] -> Yd[1]), sizeof(double), ndof, prog); 
        res_write(&(problem -> node[i] -> Y_o[1]), sizeof(double), ndof, prog);
        res_write(&(problem -> node[i] -> Yd_o[1]), sizeof(double), ndof, prog);
        res_write(&(problem -> node[i] -> Y_f[1]), sizeof(double), ndof, prog);
        res_write(&(problem -> node[i] -> Y_o_f[1]), sizeof(double), ndof, prog);

        res_write(&(problem -> node[i] -> x), sizeof(double), 24, prog);  // 24d
        res_write(&(problem -> node[i] -> pay), sizeof(double), 4, prog); // 4d
    } 

    for (i = 1 ; i <= problem -> num_segments ; i++) {
        k = problem -> segment[i] -> first_active -> number;
        res_write(&k, sizeof(int), 1, prog);                        // 1i
        k = problem -> segment[i] -> last_active -> number;
        res_write(&k, sizeof(int), 1, prog);                        // 1i
                                                                    // 14d
        res_write(&(problem -> segment[i] -> top_wet), sizeof(double), 14, prog);
    }
}
