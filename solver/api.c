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

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <malloc.h>
# include <unistd.h>
# include <math.h>
# include "compress.h"
# include "problem.h"
# include "error.h"
# include "Tree.h"
# include "output.h"
# include "solve.h"
# include "allocate.h"
# include "tension.h"
# include "transforms.h"
# include "segments.h"

# define NE        13
# define NJ_COMPAT 3
# define NB        7
# define NB_BRANCH 4

char	      *load_name = "cabstat.res";
char          *in_name = "cab.in";
char          *out_name = NULL;
int            unlink_input = 0;
int            static_finished = 0;

static int      output_map [MAXOUTPUT];
static int      nn;
static Node    *node;
static double **s = NULL;
static double **Y;
static double **Y_o;
static double **Yd_o;
static double **Y_static;
static double  *xc, *yc, *zc;
static double  *xco, *yco, *zco;
static double  *xdot_o;
static double  *ydot_o;
static double  *zdot_o;

static double   sample_dt = 0.0;
static double   snap_dt = 0.0;

	/*
         * load_static = 1: get IC from cabstat.res, 0: generate IC
	 * nn	       = pointer to int for numnodes 
	 * xcor	       = point to double array of X coordinates (nn)
	 * ycor	       = point to double array of Y coordinates (nn)
	 * zcor	       = point to double array of Z coordinates (nn)
	 * F	       = pointer to double arrary of force at mooring top 
	 * Fcur	       = pointer to double array of current force 
         */

static void FillFortranArrays(double tm, float *xcor, float *ycor, float *zcor, 
			      float *F, float *Fcur, int stat)
{
   int		  j;
   double	  T, Sn, Sb;
   double  	  B0, B1, B2, B3;
   double	**y;
   int		  idx;
   double         uack, vack, wack;
   double         ywind, zwind;
   double         ycurr, zcurr;
   double         bdr;
   Buoy           b;


   if (stat) {
      y = Y_static;
      idx = 4;
   }
   else {
      y = Y;
      idx = 7;      
   }

   T = Tension(y [1][nn], node [nn] -> material); 
   Sn = y [2][nn]; 
   Sb = y [3][nn]; 
   B0 = y [idx][nn]; 
   B1 = y [idx+1][nn]; 
   B2 = y [idx+2][nn]; 
   B3 = y [idx+3][nn]; 

   F [0] = -YComponent(T, Sn, Sb, B0, B1, B2, B3);
   F [1] = -ZComponent(T, Sn, Sb, B0, B1, B2, B3);
   F [2] = -XComponent(T, Sn, Sb, B0, B1, B2, B3);

   for (j = 1 ; j <= nn ; j++) {
      zcor [j-1] = node [j] -> x; 
      xcor [j-1] = node [j] -> y; 
      ycor [j-1] = node [j] -> z; 
   }

   b = problem.terminal [2] -> buoy;
   b -> draft = environment.surface - node [nn] -> x;
   bdr =  0.5*environment.rho*b -> Cdn*
          ProjectedArea(b, b -> draft);

   WindDrag (tm, b, &ywind, &zwind);
   Current (tm, environment.depth, node [nn] -> y, node [nn] -> z,
            &uack, &vack, &wack);

   ycurr = bdr*vack*fabs(vack);
   zcurr = bdr*wack*fabs(wack);

   fprintf (stderr,"draft = %g, Fc = %g, Fcurr + Fwind = %g, Fcurr = %g, Fw = %g\n", b -> draft, F [0], ycurr + ywind, ycurr, ywind);

   Fcur [0] = ycurr + ywind;
   Fcur [1] = zcurr + zwind;
   Fcur [2] = 0.0;

   return;
}

void cable_init__(int *load_static, int *num_nodes, 
                 float *x, float *y, float *z, 
                 float *F, float *Fcur, int *ierr)
{
   double     **(*static_solve) ( );
   int		i, n;
   int		nan;

   DisplayMode (0, 1);

	/*
	 * read the input file and dump debugging output if requested
	 */

   n = ReadModelFile (in_name);
   if (n)  {
      *ierr = 1;
      return;
   }

   if (problem.type != Surface) {
      error("cable library can only handle surface problems");
      *ierr = 1;
      return;
   }

	/*
	 * polish off the problem definition work that isn't
	 * actually done by the parser -- build the array of
	 * nodes from the segment tree and calculcate derived
	 * material and body properties
	 */


   node = CreateNodeArray (&nn, &nan);
   if (nan != nn) {
      error ("cable library cannot handle payout problems");
      *ierr = 1;
      return;
   }
      
   *num_nodes = nn;

   FillInObjectProperties ( );

   problem -> branch  = BuildBranchArray(problem, &(problem -> num_branch));

	/*
	 * make the output variable map -- everything for API version
	 */

   for (i = 1 ; i < MAXOUTPUT ; i++)
      output_map [i] = 0;

   output_map [MOTION] = 1;
   output_map [VEL]    = 1;
   output_map [FORCE]  = 1;
   output_map [MOMENT] = 1;
   output_map [EULER]  = 1;

   n = CheckTypeParameters ( );
   n += CheckStaticParameters ( );

   n += CheckDynamicParameters ( );
   n += CheckEnvironmentParameters ( );
 
   n += CheckMaterialProperties (0);
   n += CheckBuoyProperties ( );

   n += CheckBranchTerminalProperties ( );

   if (n)
      error ("%d errors found in problem description", n);

   if (n) {
      *ierr = 1;
      return;
   }

   if (analysis.static_solution == Shooting) 
      static_solve = ShootStaticProblem3D;
   else
      static_solve = SolveStaticProblem3D;
   
	/*
	 * solve the static problem or load static solution from file
	 */

   if (*load_static) 
      Y_static = LoadStaticSolution (load_name, node, nn, 0, 0.0);
   else {
      Y_static = static_solve (NULL, node, nn, nn, NULL, NULL);
      problem.dynstat = 1; 
   }

   if (Y_static == NULL) {
      if (*load_static)
         error ("could not load static solution from file");
      else
         error ("could not get static solution");

      *ierr = 1;
      return;
   }

   FillFortranArrays(0.0, x, y, z, F, Fcur, 1);

   *ierr = 0;
    
   return;		
}

void cable_update__(float *tm, 
                    float *x, float *y, float *z,
                    float *F, float *Fcur,
                    int *ierr)
{
   static double	scalv [] = {0, 0.07, 1000.0, 1000.0, 1.0, 1.0, 1.0,
                                       0.5, 0.5, 0.5, 0.5, 0.001, 0.001, 0.001};

   int            adapt_count [11];
   static int     max_adapt = 4; 
   int            level; 
   double	  t; 
   double	  targ;	
   double	  dt;
   int		  i, j;
   int	  	  it;
   int		  dynstat_conv;
   double	  dt1g, gdt_inv;
   static double  start;


   if (s == NULL) {
      s = (double **) malloc (sizeof(double *) * NE); s--;

      for (i = 1 ; i <= NE ; i++) { 
         s [i] = (double *) malloc(sizeof(double) * (2*NE + 1));     
         s [i] --;
      }

      Y = (double **) malloc (sizeof(double *) * NE);   Y--;
      Y_o = (double **) malloc (sizeof(double *) * NE); Y_o--;
      Yd_o = (double **) malloc (sizeof(double *) * NE); Yd_o--;

      for (i = 1 ; i <= NE ; i++) { 
         Y [i] = (double *) malloc(sizeof(double) * nn); Y [i] --;
         Y_o [i] = (double *) malloc(sizeof(double) * nn); Y_o [i] --;
      }

      for (i = 1 ; i <= NE ; i++) { 
         Yd_o [i] = (double *) malloc(sizeof(double) * nn); Yd_o [i] --;
      }     
   
      xc = (double *) malloc(sizeof(double) * nn); xc --;
      yc = (double *) malloc(sizeof(double) * nn); yc --;
      zc = (double *) malloc(sizeof(double) * nn); zc --;
      xco = (double *) malloc(sizeof(double) * nn); xco --;
      yco = (double *) malloc(sizeof(double) * nn); yco --;
      zco = (double *) malloc(sizeof(double) * nn); zco --;

      xdot_o = (double *) malloc(sizeof(double) * nn); xdot_o --;
      ydot_o = (double *) malloc(sizeof(double) * nn); ydot_o --;
      zdot_o = (double *) malloc(sizeof(double) * nn); zdot_o --;
   
	/*
	 * copy the static solution 
	 */

      for (i = 1 ; i <= nn ; i++) {
         Y [1][i] = Y_o [1][i] = Y_static [1][i];          /* strain       */
         Y [2][i] = Y_o [2][i] = Y_static [2][i];          /* Sn           */
         Y [3][i] = Y_o [3][i] = Y_static [3][i];          /* Sb           */
         Y [7][i] = Y_o [7][i] = Y_static [4][i];          /* beta_0       */
         Y [8][i] = Y_o [8][i] = Y_static [5][i];          /* beta_1       */
         Y [9][i] = Y_o [9][i] = Y_static [6][i];          /* beta_2       */
         Y [10][i] = Y_o [10][i] = Y_static [7][i];         /* beta_3       */
         Y [11][i] = Y_o [11][i] = Y_static [8][i];         /* omega_1      */
         Y [12][i] = Y_o [12][i] = Y_static [9][i];         /* omega_2      */
         Y [13][i] = Y_o [13][i] = Y_static [10][i];        /* omega_3      */

         Y [4][i] = Y_o [4][i] = 0.0;
         Y [5][i] = Y_o [5][i] = 0.0;
         Y [6][i] = Y_o [6][i] = 0.0;

         xc [i] = xco [i] = node [i] -> x;
         yc [i] = yco [i] = node [i] -> y;
         zc [i] = zco [i] = node [i] -> z;

         for (j = 1 ; j <= NE ; j++)
            Yd_o [j][i] = 0.0;

         xdot_o [i] = 0.0;
         ydot_o [i] = 0.0;
         zdot_o [i] = 0.0;
      }

      start = analysis.dt;
   }

   DisplayDynamicHeader ( );

   dt = analysis.dt;
   targ = *tm;
   t = start;

   level = 0; 
   dt = analysis.dt;

   for (i = 1 ; i <= max_adapt ; i++)
      adapt_count [i] = 0;

   for (t = start ; t <= targ + dt/100 ; t += dt) {
      fprintf (stderr,"t = %g\n", t);

      if (!level)
         dt = analysis.dt;

      it = SolveDE (DynamicDifeq3D, DynamicUpdate3D,
                    &(analysis -> dynamic_it), &(analysis -> dynamic_tolerance),
                    &(analysis -> dynamic_relaxation), scalv, 
 	            NE, NB, NB_BRANCH, NJ_COMPAT, 
		    node, nn, Y, Y_o, Yd_o, s, t, dt, 0.0, 0,
                    xc, yc, zc, xco, yco, zco);
   
      if (it) {
         if (level == max_adapt) {
            error("max adaption level exceeded");
            SetError(C_MAXADAPTEXCEEDED);
            *ierr = -1;
            return;
         }

         level ++;
         adapt_count [level] = 1;
         t -= dt;
         dt /= 10.0;
  
         for (j = 1 ; j <= nn ; j++) {
            xc [j] = xco [j];
            yc [j] = yco [j];
            zc [j] = zco [j];
            for (i = 1 ; i <= NE ; i ++)
               Y [i][j] = Y_o [i][j];
         }

         DisplayMessage("adapting, dt = %g\n", dt);

         continue;
      }

      dt1g = dt*(1.0 - analysis.gamma);
      gdt_inv = 1.0/analysis.gamma/dt;

      for (i = 1 ; i <= nn ; i++) {
         for (j = 1 ; j <= NE ; j++) {
            Yd_o [j][i] = gdt_inv*(Y [j][i] - Y_o [j][i] - dt1g*Yd_o [j][i]);
            Y_o [j][i] = Y [j][i];
         }

         xdot_o [i] = gdt_inv*(xc [i] - xco [i] - dt1g*xdot_o [i]);
         ydot_o [i] = gdt_inv*(yc [i] - yco [i] - dt1g*ydot_o [i]);
         zdot_o [i] = gdt_inv*(zc [i] - zco [i] - dt1g*zdot_o [i]);
         xco [i] = xc [i]; 
         yco [i] = yc [i]; 
         zco [i] = zc [i]; 
      }

      if (level) {
         if (adapt_count [level] == 10) {
            dt *= 10.0;
            level --;

     /*
      * this is to allow for the possibility that the previous level
      * was on the last step before we decided we needed to adapt
      * down again ...
      */

            for (j = level ; j >= 1 ; j--) {
               if (adapt_count [j] >= 10) {
                  dt *= 10.0;
                  level --;
               }
            }

            DisplayMessage("adapting back, dt = %g\n", dt);

            if (level)
               adapt_count [level] ++;
         }
         else
            adapt_count [level] ++;
      }
   }

   start = t;
   FillFortranArrays(t, x, y, z, F, Fcur, 0);
   fprintf (stderr,"fcable = %g %g %g\n", F [0], F [1], F [2]);
   fprintf (stderr,"fcurnt = %g %g %g\n", Fcur [0], Fcur [1], Fcur [2]);


   if (problem.dynstat) {
      dynstat_conv = CheckDynstatConvergence(Y, Y_o, dt, scalv, nn, NE, 0);

      if (dynstat_conv) {
         *ierr = 1;
         return;
      }
   }

   *ierr = 0;

   return;
}

void cable_save_state__(int *ierr)
{
   int		i;
   ResFile	fp;

        /*
         * copy the latest dynamic solution into the static
         * solution in case the driver routines needs the
         * latest info for any reason (dynstat option for instance)
         */

   for (i = 1 ; i <= nn ; i++) {
      Y_static [1][i] = Y [1][i];
      Y_static [2][i] = Y [2][i];
      Y_static [3][i] = Y [3][i];
      Y_static [4][i] = Y [7][i];
      Y_static [5][i] = Y [8][i];
      Y_static [6][i] = Y [9][i];
      Y_static [7][i] = Y [10][i];
      Y_static [8][i] = Y [11][i];
      Y_static [9][i] = Y [12][i];
      Y_static [10][i] = Y [13][i];
   }

   fp = res_open(load_name, "wb");
   if (fp == NULL) {
      error("could not open file %s for writing", load_name);
      *ierr = 1;
      return;
   }

   InitializeResultsFile(fp, nn, problem.title, output_map, 1, 0, 0);
   WriteStaticSolution(node, nn, Y_static, fp, output_map, 1, 0);
   res_close(fp);

   *ierr = 0;

   return;
}

void cable_free_mem__(int *ierr)
{
   int	i;

   xc ++; free(xc);
   yc ++; free(yc);
   zc ++; free(zc);
   xco ++; free(xco);
   yco ++; free(yco);
   zco ++; free(zco);
   xdot_o ++; free(xdot_o);
   ydot_o ++; free(ydot_o);
   zdot_o ++; free(zdot_o);

   for (i = 1 ; i <= NE ; i++) {
      s [i] ++; free (s [i]);
   }
   for (i = 1 ; i <= NE ; i++) {
      Y [i] ++; free (Y [i]);
      Y_o [i] ++; free (Y_o [i]);
   }

   for (i = 1 ; i <= NE ; i++) {
      Yd_o [i] ++; free (Yd_o [i]);
   }

   s ++; free (s);
   Y ++; free (Y);
   Y_o ++; free (Y_o);
   Yd_o ++; free (Yd_o);

   *ierr = 0;

   return;
}
