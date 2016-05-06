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

/*****************************************************************************
 * 
 * File:	catenary.c
 *
 * Description:	routines to generate the initial guess for the static
 * 	 	solution algorithms.  That guess is basically an inextensible
 *		catenary solution with forcing only at the endpoints
 *		(i.e., no drag over the cable segments).
 *
 *		The 2D and 3D specific routines in static_2d.c and
 *		static_3d.c must still provide routines to actually
 *		solve the catenary problem ->   This stuff just figures out how
 *		to specify that problem -> 
 *
 *****************************************************************************/

# include <stdio.h>
# include <math.h>
# include <malloc.h>
# include "compress.h"
# include "problem.h"
# include "tension.h"
# include "error.h"
# include "solve.h"
# include "segments.h"

# define SQR(a) ((a)*(a))

extern Problem *problem;
extern Environment *environment;
extern Analysis *analysis;

void EquivalentWeight (seg, ns, w0, EA, length)
   Segment	*seg;
   int		 ns;
   double	*w0;
   double	*EA;
   double	*length;
{
   double	weight;
   double	ws;
   double	segw;
   int		i;
   int		j;
   double	k_equiv, k_total; 

   weight = 0.0;
   *length = 0.0;
   k_total = 0.0;
   
   for (i = 1 ; i <= ns ; i++) {
      segw = 0.0;

      ws = seg [i] -> length * seg [i] -> material -> wet;
      segw += ws;

      if (seg [i] -> connector != NULL)
         segw += seg [i] -> connector -> wet;


      for (j = 1 ; j <= seg [i] -> num_attach ; j++) {
         segw += seg [i] -> attach [j].object -> wet 
                 * seg [i] -> attach [j].num_nodes;         
      }

      weight += segw;
/*
      if (weight <= 0.0)
         weight = 0.0;
*/
      k_equiv = seg [i] -> material -> EA / seg [i] -> length;
      k_total += 1.0/k_equiv;

      *length += seg [i] -> length;
   }
/*
   if (weight <= 0.0) 
      weight = seg [ns] -> length * seg [ns] -> material -> w;
*/
   k_total = 1.0 / k_total;
   *EA = k_total * (*length);

   fprintf(stderr,"Equivalent weight (excluding terminals) = %g\n", weight);

   *w0 = weight / *length;

   return;
}

	/*
	 * use two-dimensional Newton-Raphson to calculate the
	 * vertical and horizontal forces associated with an
	 * inclined catenary with the top of the catenary at	
	 * vertical coordinate x and horizontal coordinate y
	 */
  
void CatenaryForces(w0, EA, L, x, y, Hr, Vr)
   double	w0;
   double	EA;
   double	L;
   double	x;
   double	y;
   double	*Hr;
   double	*Vr; 
{
   double	H, V;
   double	weight;
   int		it;
   int		converged;
   double	f, g;
   double	dfdH, dfdV;
   double	dgdH, dgdV;
   double	t1, t2;
   double	H3, V2;
   double	H2;
   double	dH;
   double	dV;
   int		flip;

   weight = w0*L;

   V = weight / 2.0;	/* our initial guesses */
   H = fabs(weight);

   converged = 0.0;

   if (y < 0.0) {
      flip = 1;
      y = fabs(y);
   }
   else
      flip = 0;

   for (it = 1 ; it <= 1000 ; it ++) {
      V2 = V*V;
      H3 = H*H*H;
      H2 = H*H;

      t1 = sqrt(1.0 + SQR(V/H));
      t2 = sqrt(1.0 + SQR((V - weight)/H));    

      f = H/w0*(asinh(V/H) - asinh((V - weight)/H)) + H*L/EA - y;
      g = H/w0*(t1 - t2) + V*L/EA - x;

      dfdH = 1.0/w0*(asinh(V/H) - asinh((V - weight)/H)) 
             + H/w0*(-V/H2/t1 + V/H2/t2) + L/EA;
      dfdV = H/w0*(1.0/H/t1 - 1.0/H/t2);

      dgdH = 1.0/w0*(t1 - t2)
             + H/w0*(-V2/H3/t1 - (-V2/H3 - 2.0*V*weight/H3)/t2);
      dgdV = H/w0*(V/H2/t1 - (V/H2 - weight/H2)/t2) + L/EA;

      dV = (-g + dgdH*f/dfdH) / (dgdV - dfdV*dgdH/dfdH);
      dH = (-f - dfdV*dV)/dfdH;

      H += dH;
      V += dV;

      if (fabs(dH/H) < analysis -> static_tolerance 
          && fabs(dV/V) < analysis -> static_tolerance) {

         converged = 1;
         break;
      }
   }

   if (!converged)
      DisplayMessage("warning: didn't converge in cat force\n");

   if (flip)
      H = -H;

   *Hr = H;
   *Vr = V;

   return;
}

	/*
	 * our very first guess at an initial solution is the
	 * catenary for an inextensible cable with no bending stiffness
	 */

int InitialGuess (solve, n, np)
   void		(*solve) ( );
   Node		*n;
   int           np;
{
   static double w0 = 0.0;
   static double length = 0.0;
   static double weight = 0.0;
   static double EA = 0.0;
   Node		 first, from;
   Buoy		 b;
   Buoy		 b1;
   Terminal	 term;
   double	 xf, yf, zf;
   double        xthrust, ythrust, zthrust;
   double	 xforst, yforst, zforst;
   double	 hforst;
   double	 th;
   Segment      *seg;
   int	 	 nsegments;
   int		 i;
   double	 uack, vack, wack;
   double	 d, dd;
   double	 xc1, xc2;
   double	 yc1, yc2;
   double	 zc1, zc2;
   double	 xdist, ydist, zdist;
   double	 L, L1;
   double	 hdist;
   double	 U, Utheta, Umax;
   double	 phi;
 
   if (w0 == 0.0) {
      // build our own version that does not include branches
      seg = BuildSegmentArray (problem, &nsegments, 0);
      EquivalentWeight(seg, nsegments, &w0, &EA, &length);

      weight = w0 * length;
   }

	/*
	 * based on problem type decide what forces to use at the
	 * two main ends
	 */

   term = problem -> terminal [2];
   b = term -> buoy;
   xf = term -> xforce;
   yf = term -> yforce;
   zf = term -> zforce;

   switch (problem -> type) {

   case General:
      xforst = xf;
      yforst = yf;
      zforst = zf;  
      break;

   case Horizontal:
   case HorizontalDrifter:
   case Riser:
      xc2 = problem -> terminal [2] -> x;
      xc1 = problem -> terminal [1] -> x;

      yc2 = problem -> terminal [2] -> y;
      yc1 = problem -> terminal [1] -> y;

      zc2 = problem -> terminal [2] -> z;
      zc1 = problem -> terminal [1] -> z;

      hdist = sqrt(SQR(yc2 - yc1) + SQR(zc2 - zc1));

      CatenaryForces(w0, EA, length, xc2 - xc1, hdist, &hforst, &xforst);
      yforst = hforst*(yc2 - yc1)/hdist;
      zforst = hforst*(zc2 - zc1)/hdist;
      break;

   case Webster:
      phi = M_PI/2.04;

      xforst = problem -> terminal [2] -> tension * cos(phi);
      yforst = problem -> terminal [2] -> tension * sin(phi);
      zforst = 0.0;

      problem -> terminal [2] -> phi = phi;
      break;

   case Deployment:
      if (environment -> depth)
         environment -> surface = environment -> depth;
      else
         environment -> surface = length;
         
      Current(0.0, environment -> depth, 0.0, 0.0, &uack, &vack, &wack);
      vack += -problem -> terminal [1] -> yspeed.value;
      wack += -problem -> terminal [1] -> zspeed.value;

      xc1 = environment -> depth;
      xc2 = environment -> depth;

      zc1 = zc2 = 0.0;

      yc1 = 0.0;
      yc2 = -length;

      problem -> terminal [1] -> x = environment -> depth;
      problem -> terminal [1] -> y = 0.0;
      problem -> terminal [1] -> z = 0.0;

      hdist = sqrt(SQR(yc2 - yc1) + SQR(zc2 - zc1));

      fprintf(stderr,"%g %g %g %g %g\n", w0, EA, length, xc2 - xc1, hdist);
      CatenaryForces(w0, EA, length, xc2 - xc1, hdist, &hforst, &xforst);

      yforst = hforst*(yc2 - yc1)/hdist;
      zforst = hforst*(zc2 - zc1)/hdist;

      xforst = Buoyancy(b, b -> draft, environment) - b -> w;
      yforst = 0.5*environment -> rho*b -> Cdn
               *vack*fabs(vack)*ProjectedArea(b, b -> draft); 
      zforst = 0.5*environment -> rho*b -> Cdn
               *wack*fabs(wack)*ProjectedArea(b, b -> draft); 
      break;

   case Towing:
      if (environment -> depth)
         environment -> surface = environment -> depth;
      else
         environment -> surface = length;

      Current(0.0, environment -> surface - length, 0.0, 0.0, &uack, &vack, &wack);
      vack += -problem -> terminal [2] -> yspeed.value;
      wack += -problem -> terminal [2] -> zspeed.value;
      fprintf(stderr,"%g %g %g\n", uack, vack, wack);
      b1 = problem -> terminal [1] -> buoy;

      if (problem -> terminal [1] -> profile.expr || 
          problem -> terminal [1] -> profile.value) {
         Thrust(0.0, problem -> terminal [1], &xthrust, &ythrust, &zthrust);
         problem -> terminal [1] -> xthrust.value = xthrust = 0.0;
      }
      else
         Thrust(0.0, problem -> terminal [1], &xthrust, &ythrust, &zthrust);
      
      problem -> terminal [1] -> xforce = b1 -> w - b1 -> buoyancy + n[1] -> segment -> bottom_wet;
      problem -> terminal [1] -> yforce = 
            0.5*environment -> rho*b1 -> S*vack*fabs(vack)*b1 -> Cdn;
      problem -> terminal [1] -> zforce = 
            0.5*environment -> rho*b1 -> S*wack*fabs(wack)*b1 -> Cdn;

      xforst = problem -> terminal [1] -> xforce + weight - xthrust;
      yforst = -problem -> terminal [1] -> yforce - ythrust;
      zforst = -problem -> terminal [1] -> zforce - zthrust;
  
      break;

   case Subsurface:
      Current(0.0, length, 0.0, 0.0, &uack, &vack, &wack);
      fprintf(stderr,"current = %g %g\n", vack, wack);
      xforst = b -> buoyancy - b -> w;
      yforst = 0.5*environment -> rho*b -> S*vack*fabs(vack)*b -> Cdn;
      zforst = 0.5*environment -> rho*b -> S*wack*fabs(wack)*b -> Cdn;
      fprintf(stderr,"forst = %g %g\n", yforst, zforst);
      break;

   case Surface:
      Current(0.0, environment -> depth, 0.0, 0.0, &uack, &vack, &wack);

      xforst = Buoyancy(b, b -> draft, environment) - b -> w;
      yforst = 0.5*environment -> rho*b -> Cdn
               *ProjectedArea(b, b -> draft)*vack*fabs(vack);
      zforst = 0.5*environment -> rho*b -> Cdn
               *ProjectedArea(b, b -> draft)*wack*fabs(wack);

      break;

   case Drifter:
      if (environment -> depth)
         environment -> surface = environment -> depth;
      else
         environment -> surface = length;

	/*
	 * set the initial drift speed based on maximum current
	 */

      Utheta = 0.0;
      Umax   = 0.0;

      if (!problem -> terminal [2] -> yspeed.value &&
          !problem -> terminal [2] -> zspeed.value) {

         dd = environment -> surface / (double) np;
         for (d = 0.0 ; d <= environment -> surface ; d += dd) { 
            Current(0.0, d, 0.0, 0.0, &uack, &vack, &wack);
            U = sqrt(vack*vack + wack*wack);
   
            if (U > Umax) {
               Umax = U;
               Utheta = atan2(wack ,vack);
            } 
         }
         Umax = 1.001*Umax; 
      }
      else {
          Umax = sqrt(problem -> terminal [2] -> yspeed.value*
                      problem -> terminal [2] -> yspeed.value +
                      problem -> terminal [2] -> zspeed.value*
                      problem -> terminal [2] -> zspeed.value);

          Utheta = atan2(problem -> terminal [2] -> zspeed.value,
                         problem -> terminal [2] -> yspeed.value);
      }
                      
      problem -> terminal [2] -> yspeed.value = Umax*cos(Utheta);
      problem -> terminal [2] -> zspeed.value = Umax*sin(Utheta);

      Current(0.0, 0.0, 0.0, 0.0, &uack, &vack, &wack);
      vack += -problem -> terminal [2] -> yspeed.value;
      wack += -problem -> terminal [2] -> zspeed.value;

      b1 = problem -> terminal [1] -> buoy;

      problem -> terminal [1] -> xforce = b1 -> w - b1 -> buoyancy;
      problem -> terminal [1] -> yforce =
            0.5*environment -> rho*b1 -> S*vack*fabs(vack)*b1 -> Cdn;
      problem -> terminal [1] -> zforce =
            0.5*environment -> rho*b1 -> S*wack*fabs(wack)*b1 -> Cdn;

      Current(0.0, environment -> depth, 0.0, 0.0, &uack, &vack, &wack);
      vack += -problem -> terminal [2] -> yspeed.value;
      wack += -problem -> terminal [2] -> zspeed.value;

      xforst = problem -> terminal [1] -> xforce + weight;

      b -> draft = Draft(xforst + b -> w, b, environment);
      if (b -> draft > b -> max_draft || b -> draft < 0) {
         DisplayMessage("not have enough draft to support system weight");
         SetError(C_CATENARYSOLUTIONFAILED);
         return 1;
      }

      b -> S = ProjectedArea(b, b -> draft);
      
      yforst = 0.5*environment -> rho*b -> Cdn*b -> S*vack*fabs(vack);
      zforst = 0.5*environment -> rho*b -> Cdn*b -> S*wack*fabs(wack);
      yforst = -problem -> terminal [1] -> yforce; // - ythrust; ????
      zforst = -problem -> terminal [1] -> zforce; // - zthrust; ????
      break;

   default:
      DisplayMessage("unrecognized problem type");
      SetError(C_CATENARYSOLUTIONFAILED);
      return 1;
   }

   if (term -> initial_xforce 
       || term -> initial_yforce 
       || term -> initial_zforce) {
      xforst = term -> initial_xforce; // 0; // 20; // 5.25;
      yforst = term -> initial_yforce; // 20; //
      zforst = term -> initial_zforce; // -20;
   }

   problem -> terminal [2] -> xforce = xforst;
   problem -> terminal [2] -> yforce = yforst;
   problem -> terminal [2] -> zforce = zforst;

   hforst = sqrt(yforst*yforst + zforst*zforst);

   if (hforst == 0.0) {
      DisplayMessage ("there must be horizontal static forcing");
      SetError(C_CATENARYSOLUTIONFAILED);
      return 1;
   }

   if (fabs(yforst) > 1e-6 && fabs(zforst) < 1e-6) 
      th = yforst < 0 ? -M_PI : 0.0;
   else if (fabs(zforst) > 1e-6 && fabs(yforst) < 1e-6)
      th = zforst < 0 ? -M_PI/2.0 : M_PI/2.0;
   else
      th = atan2(zforst, yforst);

/* 
   if (fabs(th) > 0.5*M_PI)
      hforst = -hforst;
*/

   solve(n[1], 
         problem -> terminal[1] -> x, 
         problem -> terminal[1] -> y, 
         problem -> terminal[1] -> z, 
         hforst, xforst, th, w0, length);

   for (i = 1 ; i <= problem -> num_branch ; i++) {
      if (problem -> branch [i] -> w0 == 0.0)
          EquivalentWeight (problem -> branch [i] -> segment, 
     		            problem -> branch [i] -> num_segment,
	   		    &(problem -> branch [i] -> w0),
	   		    &(problem -> branch [i] -> EA),
	   		    &(problem -> branch [i] -> length));

      first = problem -> branch [i] -> segment[1] -> first_active;
      from = problem -> branch [i] -> segment_from -> last_active;

      if (problem -> branch [i] -> terminal -> xforce ||
          problem -> branch [i] -> terminal -> yforce || 
          problem -> branch [i] -> terminal -> zforce) {

         xforst = problem -> branch [i] -> terminal -> xforce;
         yforst = problem -> branch [i] -> terminal -> yforce;
         zforst = problem -> branch [i] -> terminal -> zforce;
      }
      else if (problem -> branch [i] -> terminal -> buoy &&
               problem -> type != HorizontalDrifter) {
         b = problem -> branch [i] -> terminal -> buoy;
         Current(0.0, length, 0.0, 0.0, &uack, &vack, &wack);

         xforst = b -> buoyancy - b -> w;
         yforst = 0.5*environment -> rho*b -> S*vack*fabs(vack)*b -> Cdn;
         zforst = 0.5*environment -> rho*b -> S*wack*fabs(wack)*b -> Cdn;
      }
      else {
         xc2 = problem -> branch [i] -> segment_from -> last -> x;
         yc2 = problem -> branch [i] -> segment_from -> last -> y;
   	     zc2 = problem -> branch [i] -> segment_from -> last -> z;

         if (problem -> branch [i] -> terminal -> anchor) {
	        xc1 = problem -> branch [i] -> terminal -> x;
	        yc1 = problem -> branch [i] -> terminal -> y;
	        zc1 = problem -> branch [i] -> terminal -> z;
         }
         else {
	        xc1 = problem -> branch [i] -> terminal -> node -> x;
	        yc1 = problem -> branch [i] -> terminal -> node -> y;
	        zc1 = problem -> branch [i] -> terminal -> node -> z;
         }

         L = problem -> branch [i] -> length;
         L1 =  sqrt(SQR(xc2 - xc1) + SQR(yc2 - yc1) + SQR(zc2 - zc1));

         if (L1 > L) {
            xdist = (xc2 - xc1)*L / L1;
            ydist = (yc2 - yc1)*L / L1;
            zdist = (zc2 - zc1)*L / L1;
         }
         else {
            xdist = xc2 - xc1;
            ydist = yc2 - yc1;
            zdist = zc2 - zc1;
         }
  
         hdist = sqrt(SQR(ydist) + SQR(zdist));
         printf("branch xd=%f, yd=%f, zd=%f, hd=%f\n",
                xdist, ydist, zdist, hdist); 
         CatenaryForces(problem -> branch [i] -> w0, 
			problem -> branch [i] -> EA,
 		        problem -> branch [i] -> length, 
			-xdist, hdist, &hforst, &xforst);

         yforst = -hforst*(yc2 - yc1)/hdist;
         zforst = -hforst*(zc2 - zc1)/hdist;
      }

      hforst = sqrt(yforst*yforst + zforst*zforst);

      if (fabs(yforst) > 1e-6 && fabs(zforst) < 1e-6)      
         th = yforst < 0 ? -M_PI : 0.0;
      else if (fabs(zforst) > 1e-6 && fabs(yforst) < 1e-6)
         th = zforst < 0 ? -M_PI/2.0 : M_PI/2.0;
      else
         th = atan2(zforst, yforst);
      
      solve(first, from -> x, from -> y, from -> z, hforst, xforst, th,
	        problem -> branch [i] -> w0, problem -> branch [i] -> length);

      //xforst = 1;
      //yforst = 0.1;
      //zforst = 0.1;
      problem -> branch [i] -> terminal -> xforce = xforst;
      problem -> branch [i] -> terminal -> yforce = yforst;
      problem -> branch [i] -> terminal -> zforce = zforst;
      printf("branch xf=%f, yf=%f, zf=%f\n", xforst, yforst, zforst);
   }

   if (!environment -> depth)
      environment -> surface = n [np] -> x;

   return 0;
}
