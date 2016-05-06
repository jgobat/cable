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
 * File:        setup.c
 *
 * Description: various routines for doing basic problem setup
 *
 * History:
 *
 ****************************************************************************/

# include <stdio.h>
# include <math.h>
# include "compress.h"
# include "problem.h"
# include "Tree.h"
# include "solve.h"
# include "error.h"

static int 
MaterialCalculations (Item item, void *data)
{
   Environment *e = (Environment *) data;
   Material	mat = (Material) item;
   double	g, rho;
   double	Ai, Ao;
   double	Acs;

   rho = e -> rho;
   g = e -> gravity;

   Ai = 0.25*M_PI*mat -> id*mat -> id;
   Ao = 0.25*M_PI*mat -> d*mat -> d;
   Acs = Ao - Ai;

   if (!mat -> A)
      mat -> A = Ao;

   if (!mat -> m)
      mat -> m = Acs * mat -> rho + Ai*e -> fill_density;

   if (!mat -> w)
      mat -> w = mat -> m * g;

   if (!mat -> wet)
      mat -> wet = mat -> w - (Ao * rho * g);

   mat -> rV = (mat -> m - mat -> wet/g);

   if (!mat -> amn) {
      if (mat -> Can)
         mat -> amn = rho*mat -> Can*Ao;  
      else if (mat -> Cmn)
         mat -> amn = rho*(mat -> Cmn - 1.0)*Ao;  
      else 
         mat -> amn = rho*Ao;
   }

   if (!mat -> amt) {
      if (mat -> Cat)
         mat -> amt = rho*mat -> Cat*Ao;  
      else if (mat -> Cmt)
         mat -> amt = rho*(mat -> Cmt - 1.0)*Ao;  
      else 
         mat -> amt = 0.0;
   }

   if (!mat -> EA)
      mat -> EA = mat -> E * Acs;

   if (!mat -> EI) {
      if (mat -> I)
         mat -> EI = mat -> E * mat -> I;
      else  {
          mat -> I = M_PI / 64.0 * (pow(mat -> d, 4.0) - pow(mat -> id, 4.0));
          mat -> EI = mat -> E * mat -> I;
      }
   }
   else {
      if (!mat -> I) { // EI given but no I
         if (mat -> E)
             mat -> I = mat -> EI / mat -> E;
         else if (mat -> EA)
             mat -> I = mat -> EI / mat -> EA * Acs;  
      }
   }

   if (!mat -> G) 
      mat -> G = mat -> E / 2.0 / (1.0 + mat -> nu);
   
   if (!mat -> GJ) {
      if (mat -> J) 
         mat -> GJ = mat -> G * mat -> J;
      else
         mat -> GJ = mat -> G * M_PI / 32.0 * (pow(mat -> d, 4.0) 
		  			    - pow(mat -> id, 4.0));
   }

   mat -> Cdn.value *= 0.5 * rho * mat -> d;

   mat -> Cdt.value *= 0.5 * rho * mat -> d * M_PI;

   if (mat -> type == Nonlinear && !mat -> EA)
      mat -> EA = EvalCode(mat -> T [1].expr, NULL, 0, 0.01, 
                           0, 0, 0, CURRNODEDATA)/0.01;

   mat -> mud_b = 2.0*sqrt((mat -> amn + mat -> m)*e -> bottom_stiffness)*e -> bottom_damping;

   return 0;
}

static int 
ConnectorCalculations (Item item, void *data)
{
   Environment *e = (Environment *) data;
   Connector	c = (Connector) item;
   double	g, rho;

   rho  = e -> rho;
   g = e -> gravity;

   c -> Cdn = 0.5*rho*0.25*M_PI*(c -> d * c -> d)*c -> Cdn;
   c -> Cdt = 0.5*rho*0.25*M_PI*(c -> d * c -> d)*c -> Cdt;

   if (!c -> am)
      c -> am = 0.5*M_PI*pow(c -> d, 3.0)/6.0*rho;
/*
   if (c -> length == 0.0)
      c -> length = c -> d;
*/
   return 0;
}

static int 
AnchorCalculations (Item item, void *data)
{
   Environment *e = (Environment *) data;
   Anchor	a = (Anchor) item;
   double	g, rho;

   rho  = e -> rho;
   g = e -> gravity;

   a -> Cdn = 0.5*rho*0.25*M_PI*(a -> d * a -> d)*a -> Cdn;
   a -> Cdt = 0.5*rho*0.25*M_PI*(a -> d * a -> d)*a -> Cdt;
   a -> am = 0.5*M_PI*pow(a -> d, 3.0)/6.0*rho;

   return 0;
}

static int 
BuoyCalculations (Item item, void *data)
{
   Environment *e = (Environment *) data;
   Buoy		b = (Buoy) item;
   double	d3;
   double	admass;

   if (b -> type == Ship || b -> type == Platform)
      return 0;
 
   if (!b -> w)
      b -> w = e -> gravity*b -> m;

   if (b -> type == Cylinder) 
      b -> max_draft = b -> h;
   else if (b -> type == Sphere) 
      b -> max_draft = b -> d;
   else if (b -> type == Capsule) 
      b -> max_draft = b -> d;
   else if (b -> type == Axisymmetric && b -> diameters) 
      b -> max_draft = b -> diameters [b -> num_diameters].level;
   else
      b -> max_draft = 0.0;
   
   if (!b -> S)
      b -> S = ProjectedArea (b, b -> max_draft);

   if (!b -> buoyancy) 
      b -> buoyancy = Buoyancy (b, b -> max_draft, e);

   b -> min_draft = Draft(b -> w, b, e);

   	/*
	 * derived constants for Morison forcing
	 */

   if (b -> type == Cylinder) {
      admass = M_PI*b -> d*b -> d*b -> h/4.0*e -> rho;
      b -> Mmma = b -> m + admass;
      b -> Marma = 2.0*admass;
   }
   else {
      d3 = pow(b -> d, 3.0);
      admass = 0.5*M_PI*d3/6.0*e -> rho;
   
      b -> Mmma  = b -> m + admass;
      b -> Marma = 3.0*admass;
   }

   b -> Mdr   = 0.5*e -> rho*b -> S*b -> Cdn;

	/*
	 * if no Cd is set for wind, use the Cd in water
	 */

   if (!b -> Cdw)
      b -> Cdw = b -> Cdn;

   return 0;
}

double 
dispersion(double w, double g, double H)
{
   double	c;
   double	k;

   if (!g)
      return 0.0;
     
   if (!H) {
      k = w*w/g;
      return k;
   }

   c = w*w*H / g;
   if (c > 20)
      k = c/H;
   else if (c > 2)
      k = c*(1.0 + 2.0*exp(-2.0*c) - 12.0*exp(-4.0*c))/H;
   else
      k = sqrt(c)*(1.0 + 0.169*c + 0.031*c*c)/H;

   return k;
}

void 
FillInEnvironment (Environment *e)
{
   int     i,j;
   double  w, k;

  
   for (i = 1 ; i <= 3 ; i++) {
      for (j = 0 ; j < e -> num_components [i] ; j++) {
         if (e -> period [i][j]) {
            w = 2.0*M_PI / e -> period [i][j];
            k = dispersion (w, e -> gravity, e -> depth);

            e -> omega [i][j] = w;
            e -> k [i][j]     = k;
         }
      }
   }

   e -> surface = e -> depth;

   return;
}

void 
FillInAnalysis(Problem *p, Analysis *a)
{
   double rho;

   if (!a -> dynamic_it)
      a -> dynamic_it = a -> maxit; 
   if (!a -> static_it)
      a -> static_it = a -> maxit; 
   if (!a -> outer_it) {
      if (a -> maxit)
         a -> outer_it = a -> maxit; 
      else
         a -> outer_it = a -> static_it;
   }
   if (!a -> shooting_it) 
      a -> shooting_it = a -> static_it;

   if (!a -> dynamic_tolerance)
      a -> dynamic_tolerance = a -> tolerance; 
   if (!a -> static_tolerance)
      a -> static_tolerance = a -> tolerance; 
   if (!a -> outer_tolerance) {
      if (a -> tolerance)
         a -> outer_tolerance = a -> tolerance; 
      else
         a -> outer_tolerance = a -> static_tolerance;
   }

   if (!a -> dynamic_relaxation)
       a -> dynamic_relaxation = a -> relaxation;
   if (!a -> static_relaxation)
       a -> static_relaxation = a -> relaxation;
   if (!a -> outer_relaxation) {
      if (p -> type == Surface)
         a -> outer_relaxation = 0.95;
      else if (p -> type == Horizontal)
         a -> outer_relaxation = 5.0;
      else if (p -> type == HorizontalDrifter)
         a -> outer_relaxation = 0.01;
      else if (p -> type == Webster)
         a -> outer_relaxation = 0.95;
   }     

   if ((p -> type == Horizontal || p -> type == HorizontalDrifter) && a -> static_outer_method < Fixed)
        a -> static_outer_method = Fixed; 

   if (!a -> gamma) {
      rho = a -> sp_rho;

      a -> alpha_k = rho / (rho - 1.0);
      a -> alpha_m = (3.0*rho + 1.0) / (2.0*rho - 2.0);
      a -> gamma = 0.5 - a -> alpha_m + a -> alpha_k;
   }
 
   return;
}

void 
FillInObjectProperties (Problem *p, Environment *e )
{
   TreeSetAndIterate(p -> buoy_tree, BuoyCalculations, e);
   TreeSetAndIterate(p -> material_tree, MaterialCalculations, e);
   TreeSetAndIterate(p -> connector_tree, ConnectorCalculations, e);
   TreeSetAndIterate(p -> anchor_tree, AnchorCalculations, e);

   return;
}
