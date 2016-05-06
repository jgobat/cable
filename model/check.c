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
 * File:        check.c
 *
 * Description: routines for checking the validity (and presence) of
 *		analysis parameters
 *
 * History:
 *
 ****************************************************************************/

# include <stdio.h>
# include <math.h>
# include "compress.h"
# include "problem.h"
# include "error.h"

extern Problem *problem;             
extern Environment *environment;
extern Analysis *analysis;

static int mat_count;
static int twoD;

static int CheckSegment (item)
   Item	item;
{
   Segment	s = (Segment) item;
   Material	m;
   int		i;

   m = s -> material;

   if (m -> type == Nonlinear) {
      for (i = 0 ; i < 3 ; i++) {
         if (!m -> T [i].expr) {
            error ("%d derivative of T(e) must be defined for material %s",
                   i, m -> name);
            mat_count ++;
         }
      }
   }
   if (!m -> EA) {
      error ("material %s has no axial stiffness (EA)", m -> name);
      mat_count ++;
   }
/*
   if (!m -> EI) {
      error ("material %s has no bending stiffness (EI)", m -> name);
      mat_count ++;
   }

   if (!m -> GJ && !twoD) {
      error ("material %s has no torsional stiffness (GJ)", m -> name);
      mat_count ++;
   }
*/
   return 0;
}
 
static int CheckBuoy (b)
   Buoy		b;
{
   int		count;

   count = 0;

   switch (b -> type) {

   case Cylinder:
      if (!b -> h) {
         error ("cylindrical buoys must have h defined");
         count ++;
      }
      break;

   case Sphere:   
      if (!b -> d) {
         error ("spherical buoys must have d defined");
         count ++;
      }
      break;

   case Capsule: 
      if (!b -> d || !b -> h) {
         error ("capsule buoys must have d and h defined");
         count ++;
      }
      break;
 
   case Axisymmetric:
      if (!b -> diameters || !b -> num_diameters) {
         error ("axisymmetric buoys must have diameters defined");
         count ++;
      }
      break;

   case Ship:
   case Platform:
      break;

   default:
      error ("buoy %s has undefined type", b -> name);
      count ++;
 
      break;
   }

   return count;
}

int CheckBuoyProperties ()
{
   int	count;
   int	i;

   count = 0;

   for (i = 1 ; i <= 2 ; i++)
      if (problem -> terminal [i] -> buoy && problem -> type != General)
         count += CheckBuoy (problem -> terminal [i] -> buoy);

   return count;
}

int CheckBranchTerminalProperties ()
{
   int	count;
   int  i;

   count = 0;

   for (i = 1 ; i <= problem -> num_branch ; i++) 
      if (problem -> branch [i] -> terminal -> buoy) 
         count += CheckBuoy (problem -> branch [i] -> terminal -> buoy);

   return count;
}
 
int CheckMaterialProperties (flag)
   int  flag;
{
   twoD = flag;
 
   mat_count = 0;

   TreeSetIterator (problem -> segment_tree, CheckSegment, NULL);
   TreeIterate (problem -> segment_tree);

   return mat_count;
}

int CheckTypeParameters ()
{
   int		count;

   count = 0;

   switch (problem -> type) {

   case General:
      if (!problem -> terminal [2] -> xforce || !problem -> terminal [2] -> yforce) {
         error ("you must specify x-force and y-force for general problems");
         count ++;
      }
      break;

   case Towing:
      if (!problem -> terminal [2] -> yspeed.value && !environment -> Uy.value &&
          !problem -> terminal [2] -> zspeed.value && !environment -> Uz.value && !environment -> current_file) {
          error ("ship speed must be initially non-zero in towing problems"); 
          count ++;
      }
      break;
   
   case Surface:
      if (!environment -> depth) { 
         error ("you must specify a depth value in surface problems");
         count ++;
      }
      if (!problem -> terminal [2] -> buoy) {
         error ("second terminal must be a buoy in surface problems");
         count ++;
         break;
      }
      if (!problem -> terminal [2] -> buoy -> Cdn) {
         error ("second terminal buoy must have nonzero Cdn");
         count ++;
      }
/*
      if (environment -> forcing > Force) {
         fprintf (stderr,"%d\n", environment -> forcing);
         error ("invalid forcing method for this problem type");
         count ++;
      }
*/
      break;

   case Deployment:
      if (!environment -> depth) { 
         error ("you must specify a depth value in deployment problems");
         count ++;
      }
      if (!problem -> terminal [2] -> buoy) {
         error ("second terminal must be a buoy in deployment problems");
         count ++;
         break;
      }
      if (!problem -> terminal [2] -> buoy -> Cdn) {
         error ("second terminal buoy must have nonzero Cdn");
         count ++;
      }
      if (!problem -> terminal [1] -> anchor) {
         error ("first terminal must be an anchor in deployment problems");
         count ++;
      }
      if (problem -> terminal [1] -> anchor -> wet <= 0.0) {
         error ("anchor must be negatively buoyant in deployment problems");
         count ++;
      }

      break;

   case Subsurface:
      if (!problem -> terminal [2] -> buoy) {
         error ("second terminal must be a buoy in subsurface problems");
         count ++;
         break;
      }
      if (!problem -> terminal [2] -> buoy -> Cdn) {
         error ("second terminal buoy must have nonzero Cdn");
         count ++;
      }
      break;

   case Horizontal:
   case Riser:
   case Webster:
      if ((problem -> terminal [2] -> x == problem -> terminal [1] -> x) &&
          (problem -> terminal [2] -> y == problem -> terminal [1] -> y) &&
          (problem -> terminal [2] -> z == problem -> terminal [1] -> z)) {
         error ("second terminal must be positioned away from first terminal");
         count ++;
      }

      break;

   }

   return count;
}

int CheckStaticParameters ()
{
   int	count;
   
   count = 0;

   if (analysis -> relax_up > analysis -> relax_down) {
      error ("upwards relaxation adaptive factor must be <= downwards factor");
      count ++;
   }
   if (!analysis -> static_it) {
      error ("no static iteration limit given, (static-iterations)");
      error ("you must at least specify a global maximum iteration limit (max-iterations)");
      count ++;
   }
   if (!analysis -> static_relaxation) {
      error ("no static relaxation factor given, (static-relaxation)");
      error ("you must at least specify a global relaxation factor (relaxation)");
      count ++;
   }
   if (!analysis -> static_tolerance) {
      error ("no static tolerance factor given, (static-tolerance)");
      error ("you must at least specify a global tolerance factor (tolerance)");
      count ++;
   }

   if (analysis -> mesh_amplify && !analysis -> mesh_smooth) {
      error ("you must give a mesh smoothing length to use auto meshing");
      count ++;
   }

   return count;
}

int CheckDynamicParameters ()
{
   int	count;
   double factor;

   count = 0;

   if (!analysis -> dynamic_it) {
      error ("no dynamic iteration limit given, (dynamic-iterations)");
      error ("you must at least specify a global maximum iteration limit (max-iterations)");
      count ++;
   }
   if (!analysis -> dynamic_relaxation) {
      error ("no dynamic relaxation factor given, (dynamic-relaxation)");
      error ("you must at least specify a global relaxation factor (relaxation)");
      count ++;
   }
   if (!analysis -> dynamic_tolerance) {
      error ("no dynamic tolerance factor given, (dynamic-tolerance)");
      error ("you must at least specify a global tolerance factor (tolerance)");
      count ++;
   }
   if (!analysis -> dt) {
      error ("you must provide a time step for dynamic analysis (time-step)");
      count ++;
   }
   if (!analysis -> duration) {
      error ("you must provide a duration for dynamic analysis (duration)");
      count ++;
   }
   if (analysis -> alpha_k > 0.5 || analysis -> alpha_m > 0.5) {
      error ("alpha_k and alpha_m integration parameters must be <= 0.5");
      count ++;
   }
   if (analysis -> sp_rho != -2.0 && (analysis -> sp_rho < -1.0 || analysis -> sp_rho >= 1.0)) {
      error ("spectral radius (dynamic-rho) must be >= -1 and < 1.0");
      count ++;
   }
   if (analysis -> dynamic_var_smooth.value > 1 || analysis -> dynamic_var_smooth.value <= 0) {
        error("smoothing parameter (dynamic-var-smoothing) must be > 0 and <= 1");
        count ++;
   }

   factor = analysis -> gamma + analysis -> alpha_m - analysis -> alpha_k;
   if (factor < 0.5 - 1e-6 && factor > 0.5 + 1e-6) 
      error ("warning: time integration is not O(dt^2) accurate given current parameters");
 
   return count;
}

int CheckEnvironmentParameters ()
{
   int	 i, j;
   char *dir_names[] = {NULL, "z", "x", "y"};
   int   count;

   count = 0;

   for (i = 1 ; i <= 3 ; i++) {
      for (j = 0 ; j < environment -> num_components [i] ; j++) {
         if (environment -> amplitude [i][j] && !environment -> period [i][j]) { 
            error("%s forcing component %d has non-zero amplitude but no period", dir_names [i], j + 1);
            count ++;
         }
      }
   }

   if (environment -> forcing == Velocity && environment -> wave_file) {
      error("use velocity-file (not wave-file) with velocity forcing");
      count ++;
   }

   return count;
}
