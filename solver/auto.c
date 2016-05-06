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
 * File:        auto.c
 *
 * Description: contains the automatic static solution path algorithm
 *
 ****************************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <malloc.h>
# include <unistd.h>
# include <math.h>
# include "compress.h"
# include "problem.h"
# include "error.h"
# include "solve.h"
# include "control.h"

extern Problem     *problem;             
extern Environment *environment;
extern Analysis    *analysis;
extern Debug debug;

static int Dynstat(y, node, num_nodes, active, num_active, twoD)
   double	**y;
   Node          *node;
   int            num_nodes;
   Node		 *active;
   int            num_active;
   int		  twoD;
{
   int		  status;
 
   problem -> dynstat = 1;
   environment -> Uscale = 1.0;

   if (twoD)
      status = SolveDynamicProblem2D(node, num_nodes, active, num_active, 
                                     NULL, NULL, NULL, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 1, NULL, 0, NULL, 0);
   else
      status = SolveDynamicProblem3D(node, num_nodes, active, num_active, 
                                     NULL, NULL, NULL, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 1, NULL, 0, NULL, 0);

   problem -> dynstat = 0;

   return status;   
}

int AutoStaticSolve (node, num_nodes, active, num_active, twoD, static_file)
   Node         *node;
   int           num_nodes;
   Node         *active;
   int           num_active;
   int		 twoD;
   char         *static_file;
{
   Analysis	  an_bkup;
   int		  err;
   double	  x;
   ProblemType    pt;
   static double  factors[] = {1.0, 5.0, 10.0, 20.0, 50.0, 100.0};
   int		  num_factors;
   int		  status;

   num_factors = sizeof(factors);

   memcpy(&an_bkup, analysis, sizeof(Analysis));

   pt = problem -> type;

   analysis -> outer_tolerance      = 0.01;
   analysis -> static_tolerance     = 0.01;
   analysis -> dynamic_tolerance    = 1e-6;

   analysis -> static_it            = 500;
   analysis -> outer_it             = 50;
   analysis -> dynamic_it           = 20;
   analysis -> shooting_it          = 50;
   analysis -> dt                   = 0.05;
   analysis -> duration             = 5000.0;

   analysis -> dynamic_relaxation   = 1.0;  
   analysis -> static_relaxation     = 0.1;
   analysis -> outer_relaxation     = 0.95;

   analysis -> static_solution      = Shooting;

   RefreshDisplay(analysis);

   debug.status = 0;

   status = 1;

   if (static_file) {
       LoadStaticSolution(static_file, node, num_nodes, twoD, 0.0);
       problem -> dynstat = 1;
       Dynstat(NULL, node, num_nodes, active, num_active, twoD);
         DisplayMessage("converged, d = %g", 
 	  	         environment -> depth - node [num_nodes] -> x);

         memcpy(analysis, &an_bkup, sizeof(Analysis));
         RefreshDisplay (analysis);
         return 0;
   }

   for (x = 1.0 ; x <= 64.0 ; x *= 2.0) {

      environment -> Uscale = x;
      DisplayMessage("applying %4.1fx current", x);

      if (twoD)
         status = ShootStaticProblem2D(0, node, num_nodes, active, num_active, NULL, NULL);
      else
         status = ShootStaticProblem3D(0, node, num_nodes, active, num_active, NULL, NULL);

      if (problem -> solution -> userQuit) {
          return 1;
      }

      err = GetError();

      if (err == C_NOSHOOTING) {
         SetError(C_STATICSOLUTIONFAILED);
         DisplayMessage("auto solve requires shooting");
         memcpy(analysis, &an_bkup, sizeof(Analysis));
         RefreshDisplay(analysis);
         return 1;
      }

/*
   this should speed things up but could also cause trouble if we got into
   a siuation where the relaxation solutions were failing
*/
      if (err == C_BUOYNOTATSURFACE) {
         x = x/2*0.95;
         continue;
      }
      else if (err == C_MAXITERATIONSEXCEEDED) 
         continue;
     
      if (!status) {
         DisplayMessage("got %3.1fx shooting, d = %g", x,
 	  	         environment -> depth - node [num_nodes] -> x);


         problem -> dynstat = 1;

         if (twoD)
            status = SolveStaticProblem2D(1, node, num_nodes, active, num_active, NULL, NULL);
         else
            status = SolveStaticProblem3D(1, node, num_nodes, active, num_active, NULL, NULL);

         if (problem -> solution -> userQuit) {
            return 1;
         }

         if (status) {
            DisplayMessage("%3.1fx relaxation solution failed", x);
            continue;
         }
 
         DisplayMessage("got %3.1fx relaxation solution, d = %g", x,
 	  	         environment -> depth - node [num_nodes] -> x);

         Dynstat(NULL, node, num_nodes, active, num_active, twoD);
         if (problem -> solution -> userQuit) {
            return 1;
         }
         if (GetError() == C_MAXADAPTEXCEEDED) {
            error("max adaptation level exceeded");
            return 0;
         }

         DisplayMessage("converged, d = %g", 
 	  	         environment -> depth - node [num_nodes] -> x);

         memcpy(analysis, &an_bkup, sizeof(Analysis));
         RefreshDisplay (analysis);
         return 0;
      }
   }

   SetError(C_STATICSOLUTIONFAILED);
   DisplayMessage("auto solution failed");

   memcpy(analysis, &an_bkup, sizeof(Analysis));
   RefreshDisplay(analysis);

   return 1;
}
