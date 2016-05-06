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
 * File:	load.c
 *
 * Description:	routines to read result files into the y arrays so that
 *		we can pick up where we left off (not implemented) or 
 *		re-use a static solution
 *
 ***************************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <malloc.h>
# include <math.h>
# include "compress.h"
# include "problem.h"
# include "output.h"
# include "error.h"
# include "tension.h"
# include "transforms.h"

extern Problem *problem;             
extern Environment *environment;

static void ConvertStatic2Dto3D(node, num_nodes, th)
   Node		 *node;
   int	  	  num_nodes;
   double     th;
{
   double	  y2d [5];
   int		  i, j;
   double	  phi, cphh, sphh;

   for (i = 1 ; i <= num_nodes ; i++) {
      for (j = 1 ; j <= 4 ; j++)
         y2d[j] = node[i] -> Y[j];

      node[i] -> Y [1] = y2d [1];
      node[i] -> Y [2] = y2d [2];
      node[i] -> Y [3] = 0;
 
      phi  = y2d [3];
 
      cphh = cos(0.5*phi);
      sphh = sin(0.5*phi);

      node[i] -> Y[4] = 0.0;
      node[i] -> Y[5] = cphh;
      node[i] -> Y[6] = sphh;
      node[i] -> Y[7] = 0.0;

      node[i] -> Y[8] = 0.0;
      node[i] -> Y[9] = 0.0;
      node[i] -> Y[10] = y2d [4];
   }

   return;
}

static void FillStaticY (in, node, num_nodes, twoD)
   ResFile       in;
   Node		*node;
   int		 num_nodes;
   int		 twoD;
{
   int		 k;
   int		 j;
   double	*d;
   int		 numvars;
   Node		 n;

   numvars = 1;	 /* the s array */

   numvars += 3; /* motions   */
   numvars += 3; /* forces    */
   numvars += 3; /* moments   */
   numvars += 4; /* euler     */

   d = (double *) malloc(sizeof(double)*numvars);

   for (j = 0 ; j <= problem -> num_branch ; j++) {
      if (j == 0)
         n = node[1];
      else
         n = problem -> branch[j] -> first;

      while (n) {
         res_read(d, sizeof(double), numvars, in);

         k = 0;
         n -> s = d [k]; k++;
         n -> x = d [k]; k++;
         n -> y = d [k]; k++;
         n -> z = d [k]; k++;

/*
         if (n -> prev)
            n -> prev -> ds = n -> s - n -> prev -> s;
*/
         n -> Y [1] = Strain (d [k], n -> material); k++; 
         n -> Y [2] = d [k]; k++; 
 
         if (!twoD) {
             n -> Y [3] = d [k]; k++; 
         }
         else
            k++;

         if (!twoD) {
            n -> Y [8] = d [k]; k++; 	/* Om1 */
            n -> Y [9] = d [k]; k++; 	/* Om2 */
            n -> Y [10] = d [k]; k++; 	/* Om3 */
         }
         else {
            k++;
            k++;
            n -> Y [4] = d [k]; k++; 	/* Om3 */
         }

         if (!twoD) {
            n -> Y [4] = d [k]; k++; 		/* B0 */
            n -> Y [5] = d [k]; k++; 		/* B1 */
            n -> Y [6] = d [k]; k++; 		/* B2 */
            n -> Y [7] = d [k]; k++; 		/* B3 */
         }
         else   
            n -> Y [3] = d [k]; k++;			/* phi */

         n = n -> next;
      }
   }
 
   return;
}

int LoadStaticSolution (name, node, nn, twoD, theta)
   char		*name;
   Node		*node;
   int		 nn;
   int		 twoD;
   double    theta;
{
   ResFile       in;
   char		 magic [6];
   char		*title;
   char		 b;
   int		 i;
   int		 j;
   int       read_input_text;
   int		 count;
   int		 num_nodes;
   int		 depth_ref;
   int		 branch_info;
   int		 file_twoD;
   double	 depth;
   char     *text_input;
   Node	 	 last;
   Terminal	 term;
   double	 T, Sn, Sb, B0, B1, B2, B3, phi;
   char         *var_names[] = {"", "motion", "velocity",
                                    "force", "moment", "euler"};


   if (problem -> type == Drifter) {
/*
      error ("cannot load drifter static solutions from files");
      return NULL;
*/
      problem -> terminal [2] -> yspeed.value = .8;
   }

   in = res_open(name, "rb");
   if (in == NULL) {
      error ("could not open results file %s for reading", name);
      return 1;
   }

   res_read(magic, sizeof(char), 6, in);
   if (strncmp(magic, "cabres", 6) != 0) {
      error ("incorrect results file format for %s", name);
      return 1;
   }

   res_read(&b, sizeof(char), 1, in);
   if (!b) {
      error ("static solution not present in file %s", name);
      return 1;
   }

   if (b & 2)
      depth_ref = 1;
   else
      depth_ref = 0;

   if (b & 64)
      read_input_text = 1;
   else
      read_input_text = 0;


   file_twoD = b & 16;

   if (file_twoD && !twoD)  
      DisplayMessage("loaded solution is 2D, converting to 3D");
   else if (!file_twoD && twoD) {
      error ("static solution in file is from 3D algorithm");
      return 1;
   }
   
   if (b & 32)
      branch_info = 1;
   else
      branch_info = 0;

   res_read(&b, sizeof(char), 1, in); /* dynamic flag */

   if (read_input_text) {
      // i = strlen() + 1 in solver/output.c for terminating null
      res_read(&i, sizeof(int), 1, in);
      text_input = (char *) malloc(sizeof(char) * i);
      res_read(text_input, sizeof(char), i, in);
   }

   res_read(&i, sizeof(int), 1, in);
   num_nodes = i;

   if (num_nodes != nn) {
      fprintf (stderr,"%d %d\n", num_nodes, nn);
      error ("number of nodes in %s does not match the current input", name); 
      return 1;
   }

   res_read(&i, sizeof(int), 1, in);
   title = (char *) malloc(sizeof(char)*i);
   res_read(title, sizeof(char), i, in);

   if (depth_ref) {
      res_read(&depth, sizeof(double), 1, in);
/*
      if (depth !=  environment -> depth) {
         error ("depth in %s does not match current input", name);
         return 1;
      }
*/
   }

   for (j = 1 ; j <= 5 ; j++) {
      res_read(&b, sizeof(char), 1, in);
      if (!b && j != VEL) {
         error ("solutions for %s missing from %s", var_names [j], name);
         return 1;
      }
   }

   for (j = 1 ; j <= 5 ; j++) 	/* unused output map bytes */
      res_read(&b, sizeof(char), 1, in);

   if (branch_info) {
      res_read(&i, sizeof(int), 1, in);
      count = i; 

      if (count != problem -> num_branch) {
         error ("branch info in %s does not match the current input", name); 
         return 1;
      }

      for (j = 1 ; j <= count ; j++)
         res_read(&i, sizeof(int), 1, in); 
   }
   

   FillStaticY (in, node, num_nodes, file_twoD);
   if (file_twoD && !twoD) {
      fprintf(stderr,"converting 2D to 3D, theta = %g\n", theta);
      ConvertStatic2Dto3D(node, num_nodes, theta);
   }
   last = problem -> terminal[2] -> node;

   if (problem -> type == Drifter && !environment -> depth)
      environment -> surface = last -> x;

   if (problem -> terminal [2] -> buoy && environment -> depth) 
       problem -> terminal [2] -> buoy -> draft = 
                       environment -> depth - last -> x;
   
   for (j = 0 ; j <= problem -> num_branch ; j++) {
      if (j == 0) 
         term = problem -> terminal [2];
      else
         term = problem -> branch[j] -> terminal;

      last = term -> node;

      T = Tension(last -> Y[1], last -> material);
      Sn = last -> Y[2];

      if (twoD) {
         phi = last -> Y [3];
         term -> xforce = T*cos(phi) - Sn*sin(phi); 
         term -> yforce = T*sin(phi) + Sn*cos(phi);
         term -> zforce = 0.0;
      }
      else {
         Sb = last -> Y [3];
         B0 = last -> Y [4];
         B1 = last -> Y [5];
         B2 = last -> Y [6];
         B3 = last -> Y [7];
         term -> xforce = XComponent(T, Sn, Sb, B0, B1, B2, B3);
         term -> yforce = YComponent(T, Sn, Sb, B0, B1, B2, B3);
         term -> zforce = ZComponent(T, Sn, Sb, B0, B1, B2, B3);
      }
   }

   return 0;
}

long
LoadProgressPoint(char *name, Problem *problem, double restart_t, int ndof)
{
    int     nn, ns, sz, i;
    ResFile in;    
    double  t;
    int     k;
    double  d;
    long    pos;


    in = res_open(name, "rb");
    if (in == NULL) {
        error("could not open progress file %s", name);
        return -1;
    }  

    res_read(&nn, sizeof(int), 1, in);
    res_read(&ns, sizeof(int), 1, in);
    res_read(&sz, sizeof(int), 1, in);


    if (nn != problem -> num_nodes || ns != problem -> num_segments || sz <= 0) {
        error("size mismatch in progress file cannot load\n");
        res_close(in);
        return -1;
    }

    t = -1;
    while(!res_eof(in)) {
        if (res_read(&t, sizeof(double), 1, in) == 0) 
            break;

        if (t == restart_t)
            break;

        res_seek(in, sz, SEEK_CUR);
    }

   if (t != restart_t) {
        error("restart time %f not found in progress file %s", restart_t, name);
        res_close(in);
        return -2;
    }

    printf("found restart time at %ld\n", res_tell(in));

    for (i = 1 ; i <= problem -> num_nodes ; i++) {
        res_read(&k, sizeof(int), 1, in);                       
        problem -> node[i] -> active_number = k; 

        res_read(&k, sizeof(int), 1, in);                       
        if (k == 0)
            problem -> node[i] -> next_active = NULL;
        else 
            problem -> node[i] -> next_active = problem -> node[k];

        res_read(&k, sizeof(int), 1, in);                        
        if (k == 0) 
            problem -> node[i] -> prev_active = NULL;
        else
            problem -> node[i] -> prev_active = problem -> node[k];

        res_read(&k, sizeof(int), 1, in);                       
        problem -> node[i] -> position = k;

        res_read(&d, sizeof(double), 1, in);                    
        problem -> node[i] -> ds = d;

        res_read(&k, sizeof(int), 1, in);                       
        problem -> node[i] -> active = k;

        res_read(&(problem -> node[i] -> Y[1]),     sizeof(double), ndof, in);   
        res_read(&(problem -> node[i] -> Yd[1]),    sizeof(double), ndof, in); 
        res_read(&(problem -> node[i] -> Y_o[1]),   sizeof(double), ndof, in);
        res_read(&(problem -> node[i] -> Yd_o[1]),  sizeof(double), ndof, in);
        res_read(&(problem -> node[i] -> Y_f[1]),   sizeof(double), ndof, in);
        res_read(&(problem -> node[i] -> Y_o_f[1]), sizeof(double), ndof, in);

        res_read(&(problem -> node[i] -> x), sizeof(double), 24, in);  
        res_read(&(problem -> node[i] -> pay), sizeof(double), 4, in); 
    } 

    for (i = 1 ; i <= problem -> num_segments ; i++) {

        res_read(&k, sizeof(int), 1, in);                        
        problem -> segment[i] -> first_active = problem -> node[k];
        res_read(&k, sizeof(int), 1, in);                        
        problem -> segment[i] -> last_active = problem -> node[k];
                                                                    // 14d
        res_read(&(problem -> segment[i] -> top_wet), sizeof(double), 14, in);
    }

    pos = res_tell(in);

    res_close(in);
    return pos;
}
