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
 * File:	mesh.c
 *
 * Description: contains code to "optimize" the spatial discretization
 *		of a segment ... optimal in this sense means that
 *	 	we put more nodes in areas of higher curvature to
 *		help resolve boundary layers	
 *
 * History:
 *		
 ****************************************************************************/

# include <stdio.h>
# include <math.h>
# include <malloc.h>
# include <stdlib.h>
# include "compress.h"
# include "problem.h"
# include "tension.h"
# include "error.h"
# include "output.h"
# include "solve.h"
# include "segments.h"

# define SQR(a) ((a)*(a))

extern Analysis *analysis;

static double EvaluateCurve(s, y, node, num_nodes)
   double   s;
   double  *y;
   Node    *node;
   int      num_nodes;
{
   double      x;
   int	       i1;
   int         i2;  
   static int  i_start;
   double      s1;
   double      s2;
   int	       i;

   if (s > node [num_nodes] -> s) 
      return y [num_nodes];

   if (s == node [1] -> s)
     i_start = 1;

   i = i_start;
  
   if (node [i] -> s == s)  
      return y [i];

   while (!(node [i] -> s <= s && node [i+1] -> s >= s)) {
      i ++;

      if (i >= num_nodes)
         i = 1;
   } 

   if (node [i] -> s == s) {
      i_start = i;
      return y [i];
   }

   i1 = i;
   i2 = i + 1;
   s1 = node [i1] -> s;
   s2 = node [i2] -> s;

   x = y [i1] + (y [i2] - y [i1]) / (s2 - s1)*(s - s1);

   i_start = i;

   return x;
}

static double rhs(s, gamma, Om3, node, num_nodes, delta)
   double	 s;
   double	 gamma;
   double	*Om3;   
   Node	        *node;
   int		 num_nodes;
   double	 delta;
{
   double	Om_s;
   double       result;

   Om_s = EvaluateCurve(s, Om3, node, num_nodes);
   
   result = gamma / (delta*Om_s + 1.0);

   return result;
}

static void Integrate (s, gamma, Om3, node, num_nodes, dq, delta)
   double	  *s;
   double	   gamma;
   double	  *Om3;
   Node		  *node;
   int		   num_nodes;
   double	   dq;
   double	   delta;
{
   int		   k;
   double	   st;
   double	   dsdq1, dsdq2, dsdq3, dsdq4;

   for (k = 1 ; k < num_nodes ; k++) {

      st = s [k];
      dsdq1 = rhs(st, gamma, Om3, node, num_nodes, delta);
   
      st = s [k] + dsdq1/2.0;
      dsdq2 = rhs(st, gamma, Om3, node, num_nodes, delta);

      st = s [k] + dsdq2/2.0;
      dsdq3 = rhs(st, gamma, Om3, node, num_nodes, delta);

      st = s [k] + dsdq3;
      dsdq4 = rhs(st, gamma, Om3, node, num_nodes, delta);

      s [k + 1] = s [k] + ((dsdq1 + dsdq4)/6.0 + (dsdq2 + dsdq3)/3.0)*dq;
   }

   return;
}


static double *Smooth(x, node, num_nodes, L)
   double	*x;
   Node		*node;
   int		 num_nodes;
   double	 L;
{
   int		 i, j;
   double	*y;
   double	 max;

   y = (double *) malloc(sizeof(double) * num_nodes); y --;


   y [1] = fabs(x [2] - x [1]);
   max = y [1];

   for (i = 2 ; i <= num_nodes ; i++) {
      y [i] = fabs(x [i] - x [i-1]);
      if (y [i] > max)
         max = y [i];
   }

   for (j = 1 ; j <= num_nodes ; j++) 
      y [j] = y [j] / max;

   return y;
 
/*  
   i = 1;
   max = fabs(x [i]);

   s0 = node [1] -> s;

   while ((node [i] -> s - s0) < L) {
      if (fabs(x [i]) > max)
         max = fabs(x [i]);

      i = i + 1;
   }

   i = 1;
   while ((node [i] -> s - s0) < L/2) {
      y [i] = max;
      i = i + 1;
   } 

   while (node [i] -> s < node [num_nodes] -> s - L/2) {
      j = i;
      max = fabs(x [i]);

      while (node [j] -> s < node [i] -> s + L/2 && j <= num_nodes) {
         if (fabs(x [j]) > max)
            max = fabs(x [j]);

         j ++;
      }

      j = i;
      while (node [j] -> s > node [i] -> s - L/2 && j >= 1) {
         if (fabs(x [j]) > max)
            max = fabs(x [j]);

         j --;
      }

      y [i] = max;
      i ++;
   }

   for (j = i ; j <= num_nodes ; j++)  
      y [j] = max;

   max = y [1];
   for (j = 1 ; j <= num_nodes ; j++) {
      if (y [j] > max)
         max = y [j];
   }

   for (j = 1 ; j <= num_nodes ; j++) 
      y [j] = y [j] / max;
 
   return y; 
*/
}

int 
AutoGenerateMesh (Problem *problem, Node *node)
{
   double	 *Om;
   Segment	 *seg;
   int		  nseg;
   double	 *s;
   int		  i, j;
   int		  it;
   double	  x1, x2;
   double	  fx1, fx2, fx3;
   double	  err;
   double	 *curv;
   double	  delta;
   double	  gamma;
   double	  dq;
   Node		 *n;
   double	 *Om_seg;
   int		  num_nodes;
   double	  Lerr, Lerr_frac;  

   // seg = BuildSegmentArray (problem, &nseg, 1);
   seg = problem -> segment;
   nseg = problem -> num_segments;

   num_nodes = 0;
   for (i = 1 ; i <= nseg ; i++) {
      if (seg [i] -> num_nodes > num_nodes)
         num_nodes = seg [i] -> num_nodes;
   }

   Om = (double *) malloc(sizeof(double) * num_nodes);
   Om --;

   for (i = 1 ; i <= num_nodes ; i++)
      Om [i] = node[i] -> Y[4];

   s = (double *) malloc(sizeof(double) * num_nodes);
   s --;

  
   for (i = 1 ; i <= nseg ; i++) {
      num_nodes = seg [i] -> num_nodes;
      if (num_nodes < 5)
         continue;

      Om_seg = &(Om [seg [i] -> first -> number]) - 1;
      n = &(node [seg [i] -> first -> number]) - 1;
    
      dq = seg [i] -> length / (double) (num_nodes - 1);

      curv = Smooth(Om_seg, n, num_nodes, analysis -> mesh_smooth);

      delta = analysis -> mesh_amplify;
      gamma = 1.0;

      x1 = fx1 = x2 = fx2 = fx3 = 0.0;
 
      for (it = 1 ; it <= 100 ; it++) {
   
         s [1] = n [1] -> s;

         Integrate (s, gamma, curv, n, num_nodes, dq, delta);

         if (it == 1) {
            x1 = gamma;
            fx1 = (s [num_nodes] - n [num_nodes] -> s);

            gamma = delta + 1.0;
         }
         else if (it == 2) {
            x2 = gamma;
            fx2 = (s [num_nodes] - n [num_nodes] -> s);

            gamma = (x2 + x1)/2.0;   	      /* bisection */
         }
         else {
            fx3 = (s [num_nodes] - n [num_nodes] -> s);

            if (fx3*fx1 < 0)
               x2 = gamma;
            else {
               x1 = gamma;
               fx1 = fx3;
            }

            gamma = (x1 + x2)/2.0;
         }

         if (it >= 3)
            err = fabs(fx3/dq);
         else
            err = 0;

         if (it > 3 && err <= 1e-6)
            break;
      }

      if (err > 1e-6) {
         Lerr = s [num_nodes] - n [num_nodes] -> s;
         DisplayMessage("warning: couldn't refine segment %d", i);

         Lerr_frac = Lerr / (double) (num_nodes - 1);

         for (j = 2; j <= num_nodes ; j++)
            s [j] -= Lerr_frac*(j - 1); 
      }

      for (j = 1 ; j <= num_nodes ; j++) {
         n [j] -> s = s [j];

         if (j > 1)
            n [j - 1] -> ds = n [j] -> s - n [j - 1] -> s;
      }
   }


   return 0;
}
