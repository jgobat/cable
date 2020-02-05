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
 * File:	static_2d.c
 *
 * Description: contains code specific to the 2D static cable model
 *
 * History:
 *		
 ****************************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <malloc.h>
# include "compress.h"
# include "problem.h"
# include "tension.h"
# include "error.h"
# include "output.h"
# include "solve.h"
# include "segments.h"

# define SQR(a) ((a)*(a))
#define DBL_EPSILON 1e-12

extern Problem     *problem;             
extern Environment *environment;
extern Analysis    *analysis;
extern Debug debug;

static void IntegrateXY(Node start, double x0, double y0, double z0)
{
   Node		 a, a_m;
   double        sf, sf_m;
   double        ph, ph_m;
   double	 cphk, cphkm;
   double	 sphk, sphkm;

   a = start;

   a -> x = x0;
   a -> y = y0;
   a -> z = z0;

   a = a -> next_active;

   while (a && a -> active) {
   
      a_m = a -> prev_active;

      sf   = 1.0 + a -> Y[1];
      sf_m = 1.0 + a_m -> Y[1];

      ph   = a -> Y [3];
      ph_m = a_m -> Y [3];

      cphk  = cos(ph);
      cphkm = cos(ph_m);
      sphk  = sin(ph);
      sphkm = sin(ph_m);

      a -> x = a_m -> x + a_m -> ds*0.5*(sf*cphk + sf_m*cphkm);
      a -> y = a_m -> y + a_m -> ds*0.5*(sf*sphk + sf_m*sphkm);
      a -> z = 0.0;

      a = a -> next_active;
   }
}

void StaticUpdate2D (active, num_active, tm, dt) 
   Node		 *active;
   int		  num_active;
   double	  tm;			/* not used		*/
   double	  dt;			/* not used		*/
{
   Node	 	from;
   int		i, j;
   double	dx, adjustment;

   IntegrateXY(active[1], problem -> terminal [1] -> x, problem -> terminal [1] -> y, 0);
   for (i = 1 ; i <= problem -> num_branch ; i++) {
      from = problem -> branch[i] -> segment_from -> last_active;
      IntegrateXY(problem -> branch[i] -> segment[1] -> first_active, from -> x, from -> y, 0);
   }

   if (!environment -> depth)
      environment -> surface = active [num_active] -> x;

   if (problem -> terminal [1] -> buoy && environment -> depth && !problem -> terminal[2] -> x) {

      dx = active [num_active] -> x - active [1] -> x;
      adjustment = environment -> depth - dx;

      for (j = 1 ; j <= num_active ; j++) 
         active [j] -> x += adjustment;
   }

   return;
}

void StaticDifeq2D (eq_type, n, nm, ne, rhs, num_rows, s, 
		    node_unused, tm, dt, current_factor)
   EquationType	  eq_type;
   Node		  n, nm;
   int		  ne;
   int		  rhs;
   int		  num_rows;
   double	**s;
   Node		 *node_unused;
   double	  tm;			/* not used		*/
   double	  dt;			/* not used		*/
   double         current_factor;
{
   Node		   nj;	
   double	   e, Sn, phi, Om3;
   double	   e_m, Sn_m, phi_m, Om3_m;
   double	   e_j, Sn_j, phi_j;
   double	   EI, w0, drat, drap;
   double          EI_m, w0_m, drat_m, drap_m;
   double	   w_m, w;
   double	   ds;
   int		   i, j;
   double	   sign_j;
   double	   cphk, sphk, sphkm, cphkm, cphkj, sphkj;
   double          tk, tkm, tkj, tdk, tdkm, tdkj, tddk, tddkm;
   double	   sqstrk, sqstrkm;
   double	   elfak, elfakm;
   double	   ulck, vlck, ulckm, vlckm;
   double 	   uack, vack, uackm, vackm;
   double	   xforst, yforst;
   double	   xthrust, ythrust;
   double	   wind_drag;
   double	   dr;
   double	   wet;
   Connector	   c;
   Buoy		   b; 
   double	   mu;
   double	   Fb;
   double	   Fb_m;
   double	   bottom, bottom_m;
   Node	   	   nd;

   w0 = n -> material -> wet;
   EI = n -> material -> EI;

   e = n -> Y [1];
   Sn = n -> Y [2];
   phi = n -> Y [3];
   Om3 = n -> Y [4];

   tk    = Tension(e, n -> material);
   tdk   = TensionD(e, n -> material);
   tddk  = TensionDD(e, n -> material);

   cphk = cos(phi);
   sphk = sin(phi);

   if (e < -1.0)
      sqstrk = 0.0;
   else
      sqstrk  = sqrt(1.0 + e);

   elfak  = pow((1.0 + e), 3.0);


   if (eq_type == Cable || eq_type == Junction || eq_type == Connection) {
      ds = nm -> ds;
   
      w0_m = nm -> material -> wet;
      EI_m = nm -> material -> EI;

      e_m = nm -> Y [1];
      Sn_m = nm -> Y [2];
      phi_m = nm -> Y [3];
      Om3_m = nm -> Y [4];

      tkm   = Tension(e_m, nm -> material);
      tdkm  = TensionD(e_m, nm -> material);
      tddkm = TensionDD(e_m, nm -> material);

      cphkm = cos(phi_m);
      sphkm = sin(phi_m);

      if (e_m < -1.0)
         sqstrkm = 0.0;
      else
         sqstrkm  = sqrt(1.0 + e_m);

      elfakm = pow((1.0 + e_m), 3.0);
   }
   else {
      sqstrkm = elfakm = 0.0;
      cphkm = sphkm = 0.0;
      e_m = Sn_m = phi_m = Om3_m = 0.0;
      w0_m = EI_m = 0.0;

      tkm = tdkm = tddkm = 0.0;

      ds = 0.0;       
   }

   for (i = 1 ; i <= num_rows ; i++)
      for (j = 1 ; j <= rhs ; j++)
         s [i][j] = BLANK;


	/*
	 * anchor node
	 */

   if (eq_type == BottomBoundary) {
      if (problem -> type == Towing || problem -> type == Drifter) {

         s [1][4]    = 1.0;
         s [1][rhs]  = Om3;

	/*
	 * in some versions I have protected the following with
	 * if (Towing), but the FHD seemed to work better in very	
	 * strong shears with the ability to change the yforce during
	 * inner iterations
	 */
	 		
         Current(0.0, n -> x, n -> y, 0.0, &uack, &vack, NULL);

         vack += -problem -> terminal [2] -> yspeed.value;

         b = problem -> terminal [1] -> buoy;

         problem -> terminal [1] -> yforce = 
            0.5*environment -> rho*b -> S*vack*fabs(vack)*b -> Cdn;

         if (problem -> type == Towing) {
            Thrust(0.0, problem -> terminal [1], &xthrust, &ythrust, NULL);
            problem -> terminal[1] -> xthrust.value = xthrust;
            problem -> terminal[1] -> ythrust.value = ythrust;
         }
         else {
            xthrust = 0.0;
            ythrust = 0.0;
         } 

         yforst = problem -> terminal [1] -> yforce - ythrust;
         xforst = problem -> terminal [1] -> xforce - xthrust;

         s [2][1] = -tdk*cphk;
         s [2][2] = sphk;
         s [2][3] = tk*sphk + Sn*cphk;
         s [2][rhs] = xforst - (tk*cphk - Sn*sphk);

         s [3][1] = -tdk*sphk;
         s [3][2] = -cphk;
         s [3][3] = -tk*cphk + Sn*sphk;
         s [3][rhs] = -yforst - (tk*sphk + Sn*cphk);
      }
      else {
         s [1][4]    = 1.0;
         s [1][rhs]  = Om3; 
      }
   }

	/* 	
	 * buoy node
	 */

   else if (eq_type == TopBoundary) {
      if (problem -> type == Towing || problem -> type == Drifter) {
         s [1][4]    = 1.0;
         s [1][rhs]  = Om3;
      }
      else {
         s [1][4]    = 1.0;
         s [1][rhs]  = Om3;

         if (problem -> type == Surface && !problem -> dynstat) {
            b = problem -> terminal [2] -> buoy;

            xforst = problem -> terminal [2] -> xforce;
            WindDrag (0.0, b, &wind_drag, NULL);
            yforst = problem -> terminal [2] -> yforce + wind_drag;
         }
         else if (problem -> type == Subsurface && !problem -> dynstat) {
            xforst = problem -> terminal [2] -> xforce;

            b = problem -> terminal [2] -> buoy;
            Current(0.0, n -> x, n -> y, 0.0, &uack, &vack, NULL);
            yforst =  0.5*environment -> rho*b -> S*vack*fabs(vack)*b -> Cdn;
         }   
         else {
            yforst = problem -> terminal [2] -> yforce;
            xforst = problem -> terminal [2] -> xforce;
         }

         s [2][1] = -tdk*cphk;
         s [2][2] = sphk;
         s [2][3] = tk*sphk + Sn*cphk;
         s [2][rhs] = xforst - (tk*cphk - Sn*sphk);

         s [3][1] = -tdk*sphk;
         s [3][2] = -cphk;
         s [3][3] = -tk*cphk + Sn*sphk;
         s [3][rhs] = yforst - (tk*sphk + Sn*cphk);
      }
   }
   else if (eq_type == BranchStart) {
      s [1][4]    = 1.0;
      s [1][rhs]  = Om3;
   }
   else if (eq_type == BranchTerminal) {

      s [1][4]    = 1.0;
      s [1][rhs]  = Om3;

      if (n -> segment -> branch -> terminal -> loop_main_node) {
         nd = n -> segment -> branch -> terminal -> loop_main_node;

         n -> segment -> branch -> terminal -> yforce -=
             (n -> y - nd -> y)*0.01;
         n -> segment -> branch -> terminal -> xforce -=
             (n -> x - nd -> x)*0.01;

         yforst = n -> segment -> branch -> terminal -> yforce;
         xforst = n -> segment -> branch -> terminal -> xforce;
         printf("xf = %g, yf = %g, ex = %g, ey = %g\n",
                xforst, yforst, n -> y - nd -> y, n -> x - nd -> x);
      }
      else if (n -> segment -> branch -> terminal -> anchor != NULL
               || problem -> type == HorizontalDrifter) {
         yforst = n -> segment -> branch -> terminal -> yforce;
         xforst = n -> segment -> branch -> terminal -> xforce;
      }
      else {
         b = n -> segment -> branch -> terminal -> buoy;
         Current(0.0, n -> x, n -> y, 0.0, &uack, &vack, NULL);

         xforst = Buoyancy(b, environment -> surface - n -> x, environment) - b -> w;
         yforst = 0.5*environment -> rho*vack*fabs(vack)*b -> Cdn*
                  ProjectedArea(b, environment -> surface - n -> x);
      }

      s [2][1] = -tdk*cphk;
      s [2][2] = sphk;
      s [2][3] = tk*sphk + Sn*cphk;
      s [2][rhs] = xforst - (tk*cphk - Sn*sphk);

      s [3][1] = -tdk*sphk;
      s [3][2] = -cphk;
      s [3][3] = -tk*cphk + Sn*sphk;
      s [3][rhs] = yforst - (tk*sphk + Sn*cphk);
   }

	/*
	 * node between two distinct segments
	 */

   else if (eq_type == Connection || eq_type == Junction) {

      if (nm -> segment -> connector == NULL
          && nm -> segment -> connection == Spliced) {	/* no lumped mass	*/
         s [1][1]      = -tdkm;
         s [1][ne + 1] = tdk;
         s [1][rhs]    = tk - tkm;

         for (i = 2 ; i <= ne ; i++) {
            s [i][i] = -1.0;
            s [i][i + ne] = 1.0;
            s [i][rhs] = n -> Y[i] - nm -> Y[i];
         }
      }
      else {
         c = nm -> segment -> connector;
         if (c) {
             dr =  c -> Cdn;

             Current (0.0, n -> x, n -> y, 0.0, &uack, &vack, NULL);

             uack *= current_factor;
             vack *= current_factor;

             if (problem -> type == Towing || problem -> type == Drifter)
                vack -= problem -> terminal [2] -> yspeed.value;
             else if (problem -> type == Deployment)
                vack -= problem -> terminal [1] -> yspeed.value;

             wet = c -> wet + n -> segment -> bottom_wet + nm -> segment -> top_wet;
    
             if ((problem -> type == Deployment  
                  || problem -> type == Surface 
                  || problem -> type == Horizontal 
                  || problem -> type == HorizontalDrifter) && !problem -> dynstat) {

                if (n -> x - environment -> surface > DBL_EPSILON && wet < -DBL_EPSILON) {
                   wet = wet*(1.0 + tanh(50.0*(environment -> surface - n -> x)));
                   // fprintf(stderr,"wet = %f, d = %f\n",
                   //        wet, environment -> surface - n -> x);
                }

             }
        
             if (problem -> type == HorizontalDrifter) {
                 xforst = nm -> segment -> connector_xforce;
                 yforst = nm -> segment -> connector_yforce;
             }
             else {
                 xforst = yforst = 0;
             }
         }
         else {
            dr = wet = 0;
            vack = 0;
            xforst = yforst = 0;
         }

         s [1][1]      = -cphkm*tdkm;
         s [1][ne + 1] = cphk*tdk;
         s [1][2]      = sphkm;
         s [1][ne + 2] = -sphk;
         s [1][3]      = sphkm*tkm + cphkm*Sn_m;
         s [1][ne + 3] = -sphk*tk - cphk*Sn;
         s [1][rhs] = -wet - cphkm*tkm + sphkm*Sn_m + cphk*tk - sphk*Sn + xforst;

         s [2][1]      = -sphkm*tdkm;
         s [2][ne + 1] = sphk*tdk;
         s [2][2]      = -cphkm;
         s [2][ne + 2] = cphk;
         s [2][3]      = -cphkm*tkm + sphkm*Sn_m;
         s [2][ne + 3] = cphk*tk - sphk*Sn;
         s [2][rhs] = dr*vack*fabs(vack) + yforst
                         - sphkm*tkm - cphkm*Sn_m + sphk*tk + cphk*Sn;

         if (eq_type == Junction) {
            for (i = 1 ; i <= n -> segment -> junction.num_nodes ; i++) {
               nj = n -> segment -> junction.node [i];
               
               e_j = nj -> Y [1];
               phi_j = nj -> Y [3];
               Sn_j = nj -> Y [2];

               cphkj = cos(phi_j);
               sphkj = sin(phi_j);

               tkj = Tension(e_j, nj -> material);
               tdkj = TensionD(e_j, nj -> material);

               if (nj == nj -> segment -> branch -> last) 
                  sign_j = -1.0;
               else
                  sign_j = 1.0;

               s [1][(i + 1)*ne + 1] = sign_j*cphkj*tdkj;
               s [1][(i + 1)*ne + 2] = -sign_j*sphkj;
               s [1][(i + 1)*ne + 3] = sign_j*(-sphkj*tkj - cphkj*Sn_j);
               s [1][rhs] += sign_j*(cphkj*tkj - sphkj*Sn_j);
   
               s [2][(i + 1)*ne + 1] = sign_j*sphkj*tdkj;
               s [2][(i + 1)*ne + 2] = sign_j*cphkj;
               s [2][(i + 1)*ne + 3] = sign_j*(cphkj*tkj - sphkj*Sn_j);
               s [2][rhs] += sign_j*(sphkj*tkj + cphkj*Sn_j);
            }
         }

         s [3][ne + 4] = 1.0;
         s [3][rhs] = Om3;
 
         s [4][4] = 1.0;
         s [4][rhs] = Om3_m;	
      }
   }

	/*
	 * a regular old internal node, eq_type == Cable
	 */

   else {
      Current (0.0, n -> x, n -> y, 0.0, &uack, &vack, NULL);
      Current (0.0, nm -> x, nm -> y, 0.0, &uackm, &vackm, NULL);
      
      uack *= current_factor;
      vack *= current_factor;

      uackm *= current_factor;
      vackm *= current_factor;

      if (problem -> type == Towing || problem -> type == Drifter) {
         vack -= problem -> terminal [2] -> yspeed.value;
         vackm -= problem -> terminal [2] -> yspeed.value;
      }
      else if (problem -> type == Deployment) {
         vack -= problem -> terminal [1] -> yspeed.value;
         vackm -= problem -> terminal [1] -> yspeed.value;
      }

      
      ulck = cphk*uack + sphk*vack;
      vlck = -sphk*uack + cphk*vack;

      ulckm = cphkm*uackm + sphkm*vackm;
      vlckm = -sphkm*uackm + cphkm*vackm;

      DragCoeff(n, n -> material, 0.0, fabs(ulck), fabs(vlck), &drat, &drap);
      DragCoeff(nm, nm -> material, 0.0, fabs(ulckm), fabs(vlckm), &drat_m, &drap_m);

      w = n -> material -> w;
      w_m = nm -> material -> w;

      if (n -> attachment) {
         w0 += n -> attachment -> wet/ds;
	 w += n -> attachment -> m*environment -> gravity / ds;
         drat += n -> attachment -> Cdt/ds;
         drap += n -> attachment -> Cdn/ds;
      }
      if (nm -> attachment) {
         w0_m += nm -> attachment -> wet/ds;
	 w_m += nm -> attachment -> m*environment -> gravity / ds;
         drat_m += nm -> attachment -> Cdt/ds;
         drap_m += nm -> attachment -> Cdn/ds;
      }

	/*
	 * allow for the possibility of line floating on the surface
	 */

      if ((problem -> type == Surface || problem -> type == Deployment || problem -> type == Drifter) && !problem -> dynstat) {
         
         if (nm -> x - environment -> surface > DBL_EPSILON && w0_m < -DBL_EPSILON) {
            w0_m = w0_m*(1.0 + tanh(50.0*(environment -> surface - nm -> x)));
         }

         if (n -> x - environment -> surface > DBL_EPSILON && w0 < -DBL_EPSILON) {
            w0 = w0*(1.0 + tanh(50.0*(environment -> surface - n -> x)));
         }

      }

        /*
         * for riser problems there may be material coming up out of
         * the surface for which we should use the dry, not the wet,
         * weight
         */

      if (problem -> type == Riser) {
         if (nm -> x >= environment -> surface)
            w0_m = nm -> material -> w;

         if (n -> x >= environment -> surface)
            w0 = n -> material -> w;
      }

	/*
	 * for all problems there may be line on the ground
	 * and we need to account for the normal force of the
	 * bottom plus frictional forces
	 */

      bottom = Bottom(n -> y, n -> z, 0.0);
      bottom_m = Bottom(nm -> y, nm -> z, 0.0);
      mu = environment -> bottom_friction;

      if ((n -> position == TopBoundary && problem -> terminal [2] -> anchor) ||
          (n -> position == BranchTerminal && n -> segment -> branch -> terminal -> anchor)) {

         Fb_m = 0.0;
         Fb = 0.0;
      }
      else {
         if (nm -> x < bottom_m && w0_m > 0.0) {
            if (nm -> active_number == 1)
               Fb_m = w0_m;
            else
               Fb_m = (bottom_m - nm -> x)*environment -> bottom_stiffness;

            if (Fb_m > w0_m)
               Fb_m = w0_m;
         }
         else {
            Fb_m = 0.0;
         }

         if (n -> x < bottom && w0 > 0.0) {
            Fb = (bottom - n -> x)*environment -> bottom_stiffness;

            if (Fb > w0)
               Fb = w0;
         }
         else { 
            Fb = 0.0;
         }
      }
 
      s [1][1]      = tddkm*(e - e_m) - (tdk + tdkm)
                      + 0.5*ds*drat_m*ulckm*fabs(ulckm)/sqstrkm;
      s [1][ne + 1] = tddk*(e - e_m) + (tdk + tdkm)
                      + 0.5*ds*drat*ulck*fabs(ulck)/sqstrk; 
      s [1][2]      = -(phi - phi_m);
      s [1][ne + 2] = s [1][2];
      s [1][3]      = (Sn + Sn_m) + ds*((w0 - Fb_m)*sphkm 
           + 2.0*drat_m*ulckm*sign(ulckm)*sqstrkm*(-uackm*sphkm + vackm*cphkm));
      s [1][ne + 3] = -(Sn + Sn_m) + ds*((w0 - Fb)*sphk
           + 2.0*drat*sqstrk*sign(ulck)*ulck*(-uack*sphk + vack*cphk));

      s [1][rhs] = (tdk + tdkm)*(e - e_m) - (Sn + Sn_m)*(phi - phi_m)
                   - ds*((w0_m - Fb_m)*cphkm + (w0 - Fb)*cphk
			 + mu*Fb_m + mu*Fb 
 		 	 - drat_m*ulckm*fabs(ulckm)*sqstrkm 
                         - drat*ulck*fabs(ulck)*sqstrk);
 
      s [2][1]      = tdkm*(phi - phi_m) 
                      + 0.5*ds*drap_m*vlckm*fabs(vlckm)/sqstrkm;
      s [2][ne + 1] = tdk*(phi - phi_m)
                      + 0.5*ds*drap*vlck*fabs(vlck)/sqstrk;
      s [2][2]      = -2.0;
      s [2][ne + 2] = 2.0;
      s [2][3]      = -(tk + tkm) + ds*((w0_m - Fb_m)*cphkm 
           + 2.0*drap_m*sqstrkm*sign(vlckm)*vlckm*(-uackm*cphkm - vackm*sphkm));
      s [2][ne + 3] = tk + tkm + ds*((w0 - Fb)*cphk  
           + 2.0*drap*sqstrk*sign(vlck)*vlck*(-uack*cphk - vack*sphk));
     
      s [2][rhs] = 2.0*(Sn - Sn_m) + (tk + tkm)*(phi - phi_m)
                   + ds*((w0 - Fb)*sphk + (w0_m - Fb_m)*sphkm
                         + drap_m*vlckm*fabs(vlckm)*sqstrkm 
			 + drap*vlck*fabs(vlck)*sqstrk);

      s [3][3]      = -2.0;
      s [3][ne + 3] = 2.0;
      s [3][4]      = -ds;
      s [3][ne + 4] = -ds;
      s [3][rhs] = 2.0*(phi - phi_m) - ds*(Om3 + Om3_m);
     
      s [4][1]      = 3.0*ds*Sn_m*(1.0 + e_m)*(1.0 + e_m);
      s [4][ne + 1] = 3.0*ds*Sn*(1.0 + e)*(1.0 + e);
      s [4][2]      = ds*elfakm;
      s [4][ne + 2] = ds*elfak;
      s [4][4]      = -(EI + EI_m);
      s [4][ne + 4] = (EI + EI_m);
      s [4][rhs] = (EI + EI_m)*(Om3 - Om3_m)
                   + ds*(Sn*elfak + Sn_m*elfakm);
   }
       
   return;
}

void BendingDerivatives2D (Node start)
{
   Node		 np, nm;
   double	 dOm;
   double	 slopel, sloper;
   Node		 n;

   n = start;
   while (n) {
      
      np = n -> next_active;
      nm = n -> prev_active;
    
      if (n -> position == BottomBoundary ||
          n -> position == Junction ||
          n -> position == Connection ||
	  n -> position == BranchStart) 

         n -> Y [4] = (np -> Y[3] - n -> Y [3]) / n -> ds;
      else if (n -> segment -> last_active == n || n -> position == TopBoundary || np == NULL) 
         n -> Y [4] = (n -> Y [3] - nm -> Y [3]) / nm -> ds;
      else {
         sloper = (np -> Y [3] - n -> Y [3])/n -> ds;
         slopel = (n -> Y [3] - nm -> Y[3])/nm -> ds;
         n -> Y [4] = (sloper*nm -> ds + slopel*n -> ds)/
                    (nm -> ds + n -> ds);
      }

      n = n -> next_active;
   }

   n = start;

   while(n) {
      np = n -> next_active;
      nm = n -> prev_active;
    
      if (n -> position == BottomBoundary ||
          n -> position == Junction ||
          n -> position == Connection ||
	  n -> position == BranchStart) 

         dOm = (np -> Y[4] - n -> Y[4]) / n -> ds;
      
      else if (n -> segment -> last_active == n || n -> position == TopBoundary || np == NULL) {
         dOm = (n -> Y [4] - nm -> Y [4]) / nm -> ds;
      }
      else {
         sloper = (np -> Y[4] - n -> Y[4])/n -> ds;
         slopel = (n -> Y[4] - nm -> Y[4])/nm -> ds;
         dOm = (sloper*nm -> ds + slopel*n -> ds)/
                    (nm -> ds + n -> ds);
      }

      n -> Y[2] = -n -> material -> EI*dOm/pow(1.0 + n -> Y[1], 3.0);

      n = n -> next_active;
   }

   return;
}

void SolveCatenary2D(start, x0, y0, z0, hforst, xforst, th, w0, length)
   Node		 start;
   double	 x0, y0, z0;
   double	 hforst;
   double	 xforst;
   double	 th;
   double	 w0;
   double	 length;
{
   Node		n;
   double       tens;
   double       p;
   double       term0, term1, term2, term3, term4, term5;
   double       phi;
   double       sv;
   double	cth;

   if (hforst == 0.0)
      hforst = 0.1;
   if (xforst == 0.0)
      xforst = -1.0;

   if (w0 == 0.0)
      w0 = 1e-6;

   term0 = hforst / w0;
   term1 = xforst / hforst;
   term2 = term1 - length/term0;
   term4 = sqrt(1.0 + term2*term2);
   term5 = log(term2 + term4);

   cth = cos(th);

   sv = 0;
   n = start;
   while (n) {

      term3 = term2 + w0*sv/hforst;
      p = term0*( log(term3 + sqrt(term3*term3 + 1.0)) - term5 );
      tens = sqrt(hforst*hforst + SQR(xforst - w0*(length - sv)));
      if (fabs(term3) > 10.0)
         phi = 1.0 / term3;
      else {
         if (fabs(term3) < 1e-6)
            phi = atan2(1.0, term3);
         else
            phi = atan(1.0/term3);
      }

      if ((xforst - w0*(length - sv)) < 0.0) 
         phi = M_PI + phi;
/*
      if (hforst < 0.0) 
         phi  = -phi;
*/
      n -> Y[3] = phi*cth;

      n -> Y[1] = Strain (tens, n -> material);

      sv += n -> ds;

      n = n -> next_active;
   }   

   BendingDerivatives2D(start);

	// this is a problem to be dealt with ---

   IntegrateXY(start, x0, y0, z0);

   return;
}

int SolveStaticProblem2D (init_loaded, node, num_nodes, active, num_active, out, output_map)
   int		 init_loaded;
   Node		*node;
   int		 num_nodes;
   Node		*active;
   int		 num_active;
   ResFile	 out;
   int		*output_map;
{
/*
   static double	scalv [] = {0, 0.07, 100.0, M_PI_4, 0.001};
*/
   static double	scalv [] = {0, 0.07, 1000.0, 1.4, 0.0034};
   int		  profiling, maneuvering;
   int            i;
   double         current_factor;
   int	          current_steps;
   double	**s;
   int		  singular;
   int		  resolved;
   double	  resolve_err;
   int		  it;
   int		  ne, nb;
   int		  njn;
   int		  nb_branch;
   int		  nj_compat;
   double         T, phi, T1, phi1;
   Buoy           b, b1;
   double         drag, drag1;
   double         uack, vack;
   Segment       *seg;
   int            nseg;

   problem -> twoD = 1;
   problem -> dynamic = 0;

   ne        = 4;	/* number of eq's (DOF) at each node         */
   nb_branch = 1;	/* number of BC equations at BranchStart     */
   nj_compat = 0;	/* number of extra compat eq's at a junction */

   if (problem -> type == Towing || problem -> type == Drifter) 
      nb = 3;		/* number of BC equations at node 1	     */
   else
      nb = 1;

   // seg = BuildSegmentArray(problem, &nseg, 1);
   seg = problem -> segment;
   nseg = problem -> num_segments;

   njn = problem -> junction_size;

   s = (double **) malloc (sizeof(double *) * (ne + njn*nj_compat)); s--;
   for (i = 1 ; i <= ne + nj_compat*njn ; i++) {
      s [i] = (double *) malloc(sizeof(double) * ((2 + njn)*ne + 1));
      s [i] --;
   }
/*
   for (i = 1 ; i <= num_nodes ; i++) {

      node[i] -> Ys = (double *) malloc(sizeof(double) * 4); node[i] -> Ys --;
      node[i] -> Y = (double *) malloc(sizeof(double) * 6); node[i] -> Y --;
      node[i] -> Y_o = (double *) malloc(sizeof(double) * 6); node[i] -> Y_o --;
      node[i] -> Yd_o = (double *) malloc(sizeof(double) * 6); node[i] -> Yd_o --;
   }
*/
   if (!init_loaded) {
      if (analysis -> static_initial_guess == Catenary || analysis -> static_solution == Catenary) {

         if (problem -> type == Surface || problem -> type == Deployment)
            problem -> terminal [2] -> buoy -> draft = 
                            problem -> terminal [2] -> buoy -> max_draft;

         i = InitialGuess (SolveCatenary2D, active, num_active); 
         if (i)
            return 1;

         if (analysis -> static_solution == Catenary) 
         return 0;
      }
      else {
         i = ShootStaticProblem2D (0, node, num_nodes, active, num_active, out, output_map);

         if (i) {
            SetError(C_SHOOTINGSOLUTIONFAILED);
            return 1;
         }

         DisplayMessage("Shooting initial forces = %g, %g, draft=%g", problem -> terminal [2] -> xforce,
                  problem -> terminal [2] -> yforce, problem->terminal[2]->buoy ? problem->terminal[2]->buoy->draft : 0.0);
      }
   }

   if (debug.status) {
      WriteStaticSolution(node, num_nodes, out, output_map, debug.decimate, 1);

      WriteDynamicHeader (out, 0.0, analysis -> outer_it, 1.0, 
			 debug.sample_it, debug.snap_it, 0.0, 0.0, 0.0, 0, NULL, 1, node);
      FakeDynamicSnapshot(out, node, num_nodes, debug.decimate, 1);
   }

   DisplayMessage("Initial forces = %g, %g", problem -> terminal [2] -> xforce,
                  problem -> terminal [2] -> yforce);

   resolved = 0;
   singular = 0;

   maneuvering = profiling = 0;

   if (problem -> type == Towing && problem -> terminal [1] -> x && problem -> terminal [1] -> y)
      maneuvering = 1;
   if (problem -> type == Towing && (problem -> terminal [1] -> profile.expr ||
       problem -> terminal [1] -> profile.value || problem -> terminal [1] -> profile_m)) 
      profiling = 1;

   DisplayStaticHeader ( );

   current_steps = analysis -> current_steps;

   for (it = 1 ; it <= analysis -> outer_it ; it++) {

      if (it == 1)
         current_steps = analysis -> current_steps;
      else
         current_steps = 0;

      for (i = 0 ; i <= current_steps ; i++) {
         if (current_steps > 0)
            current_factor = (double) (i + 1) / (double) (current_steps + 1);
         else
            current_factor = 1.0;

         fprintf(stdout ? stdout : stderr,"F = (%.6f, %.6f), speed = (%.6f, %.6f, %.6f, %.6f), x0 = %.6f, draft = %.6f\n",
                 problem -> terminal[2] -> xforce,
                 problem -> terminal[2] -> yforce,
                 problem -> terminal[1] -> xspeed.value,
                 problem -> terminal[1] -> yspeed.value,
                 problem -> terminal[2] -> xspeed.value,
                 problem -> terminal[2] -> yspeed.value, 
                 problem -> terminal[1] -> x,
                 problem -> terminal[2] -> buoy ? 
                    problem -> terminal[2] -> buoy -> draft : 0.0);

         singular = SolveDE (StaticDifeq2D, StaticUpdate2D,
                             &(analysis -> static_it), &(analysis -> static_tolerance),
                             &(analysis -> static_relaxation), scalv, 
                             ne, nb, nb_branch, nj_compat,
                             node, num_nodes, active, num_active, s, 
                             0.0, 0.0, current_factor, it);

        
         if (problem -> solution -> userQuit) {
            return 1;
         }

         if (singular) 
            break;
      }


      if (singular && problem -> type != Surface && problem -> type != Drifter)
        break;

      if (problem -> dynstat) {
         resolved = 1;
         break;
      }
   
      if (debug.status)
         FakeDynamicSnapshot(out, node, num_nodes, debug.decimate, 1);

      if (problem -> type == Surface || problem -> type == Deployment) {
         resolve_err = ResolveBuoy (active, num_active, it, 1);  
     
         if (resolve_err < 0.0)
            break; 
         else if (resolve_err < analysis -> outer_tolerance)
            resolved = 1;

#if defined (GUI) || defined (WINGUI)
         ControlPlotSnaps(problem, environment);
#endif 
/*
         else if ((problem -> type == Surface || problem -> type == Deployment) && analysis -> static_initial_guess == Catenary) {
            InitialGuess (SolveCatenary2D, active, num_active);
         }
*/ 
         DisplayMessage("(%d) draft = %6.4f, F = (%7.2f, %7.2f)", it,
                  problem -> terminal [2] -> buoy -> draft,
                  problem -> terminal [2] -> xforce, 
                  problem -> terminal [2] -> yforce);

         IntegrateXY(active[1], problem -> terminal[1] -> x, 0, 0);
      }

      else if (problem -> type == Webster) {
         resolved = ResolveWebster (it, active, num_active, 1);

         DisplayMessage("z [n] = %7.3f, x [n] = %7.3f", 
                  active [num_active] -> x, 
		  active [num_active] -> y);
      }
      else if (problem -> type == Horizontal || problem -> type == HorizontalDrifter || problem -> type == Riser) {
         resolved = ResolveAnchor (active, num_active, seg, nseg, 1);

         DisplayMessage("z [n] = %7.3f, x [n] = %7.3f", 
                  active [num_active] -> x, 
		  active [num_active] -> y);
      }
      else if (problem -> type == Drifter) {
         resolved = ResolveSpeed (active, num_active, it);
/*
         if (resolve_err < 0.0)
            break;
         else if (resolve_err < analysis -> outer_tolerance)
            resolved = 1;
         else if (!resolved && analysis -> static_initial_guess == Catenary)
            InitialGuess (SolveCatenary2D, node, num_active);

*/
         DisplayMessage("speed = %6.4f, F = (%7.2f, %7.2f)", 
                  problem -> terminal [2] -> yspeed.value,
                  problem -> terminal [2] -> xforce, 
                  problem -> terminal [2] -> yforce);
      }
      else if (profiling) {
         resolved = ResolveStartDepth(active, num_active, it);
         DisplayMessage("depth = %.1f", 
                        (environment -> depth ? 
                         environment -> depth - active [1] -> x :
                         active [num_active] -> x));
         if (resolved == -1) { /* my isn't this an ugly little break-out ... */
            resolved = 0;
            break;
         }
      }
      else
         resolved = 1;

      if (resolved)
         break;

   }

   if (problem -> dynstat && singular) 
      return 1;

   if (singular && singular != analysis -> static_it + 1) {
      error ("could not get static solution");
      SetError(C_STATICSOLUTIONFAILED);
      return 1;
   }

   if (!resolved) {
      SetError(C_MAXITERATIONSEXCEEDED);
      error ("never converged on outer iterations");
   }

   for (i = 1 ; i <= ne + nj_compat*njn ; i++) {
      s [i] ++; free (s [i]);
   }
   s ++; free (s);
   

   if (problem -> type == Horizontal || problem -> type == Riser 
       || problem -> type == HorizontalDrifter
       || problem -> type == Webster) {

      T = Tension(active[1] -> Y[1], active [1] -> material);
      phi = active[1] -> Y [3];

      problem -> terminal[1] -> xforce = - T*cos(phi) + active[1] -> Y[2]*sin(phi);
      problem -> terminal[1] -> yforce = - T*sin(phi) - active[1] -> Y[2]*cos(phi);

      DisplayMessage("Done. Fh1= %.0f, Fv1= %.0f, Fh2= %.0f, Fv2= %.0f",
              problem -> terminal[1] -> yforce,
              problem -> terminal[1] -> xforce,
              problem -> terminal [2] -> yforce,
              problem -> terminal [2] -> xforce);
   }
   else if (problem -> type == Towing) {
      T = Tension(active[num_active] -> Y[1], active [num_active] -> material);
      phi = active[num_active] -> Y [3];
      b = problem -> terminal [2] -> buoy;
      b -> draft = Draft(b -> w + T*cos(phi) - active[num_active] -> Y[2]*sin(phi), b, environment);
      Current(0.0, environment -> depth, 0, 0, &uack, &vack, NULL);
      drag = 0.5*environment -> rho*ProjectedArea(b, b -> draft)*b -> Cdn*(vack - problem -> terminal[2] -> yspeed.value)*fabs(vack - problem -> terminal[2] -> yspeed.value);
      DisplayMessage("Done. Tow T= %.1f, Tow H= %.1f, Ship T= %.1f, Ship Fz= %.2f, Ship Fx= %.2f, Ship drag= %.2f, Tow X=%.1f, Ship X=%.1f",
               Tension(active[1] -> Y[1], active [1] -> material),
               (environment -> depth ? environment -> depth - active [1] -> x :
                active [num_active] -> x),
               T, T*cos(phi) - active[num_active] -> Y[2]*sin(phi),
               T*sin(phi) + active[num_active] -> Y[2]*cos(phi), drag,
               active[1] -> y, active[num_active] -> y);

      printf("%g %g\n", environment -> depth ? environment -> depth - active [1] -> x : active[num_active] -> x, active[num_active] -> y);
   }
   else if (problem -> type == Surface || problem -> type == Deployment) {
      DisplayMessage("Done. Buoy draft= %.3f, Swet = %.3f",
               problem -> terminal [2] -> buoy -> draft,
               ProjectedArea(problem -> terminal [2] -> buoy,
                             environment -> surface - active [num_active] -> x));
      T = Tension(active[1] -> Y[1], active [1] -> material);
      phi = active[1] -> Y [3];
      DisplayMessage(" Anchor Fx=%.1f, Fz=%.1f, Wmin(mu=%g, SF=%g)=%.0f lbs (%.0f lbs steel in air)",
                     T*sin(phi) + active[1] -> Y[2]*cos(phi),
                     T*cos(phi) - active[1] -> Y[2]*sin(phi),
                     problem -> terminal[1]->friction, problem -> terminal[1]->safety,
                     ((T*sin(phi) + active[1] -> Y[2]*cos(phi))*problem -> terminal[1]->safety/problem -> terminal[1]->friction + 
                     T*cos(phi) - active[1] -> Y[2]*sin(phi))/4.4482216,
                     ((T*sin(phi) + active[1] -> Y[2]*cos(phi))*problem -> terminal[1]->safety/problem -> terminal[1]->friction + 
                     T*cos(phi) - active[1] -> Y[2]*sin(phi))/4.4482216/0.87);
   }
   else if (problem -> type == Subsurface) {
      T = Tension(active[1] -> Y[1], active [1] -> material);
      phi = active[1] -> Y [3];
      DisplayMessage("Buoy depth=%.1f, Anchor Fx=%.1f, Fz=%.1f, Wmin(mu=%g, SF=%g)=%.0f lbs (%.0f lbs steel in air)",
                     environment -> depth - active[num_active] -> x,
                     T*sin(phi) + active[1] -> Y[2]*cos(phi),
                     T*cos(phi) - active[1] -> Y[2]*sin(phi),
                     problem -> terminal[1]->friction, problem -> terminal[1]->safety,
                     ((T*sin(phi) + active[1] -> Y[2]*cos(phi))*problem -> terminal[1]->safety/problem -> terminal[1]->friction + 
                     T*cos(phi) - active[1] -> Y[2]*sin(phi))/4.4482216,
                     ((T*sin(phi) + active[1] -> Y[2]*cos(phi))*problem -> terminal[1]->safety/problem -> terminal[1]->friction + 
                     T*cos(phi) - active[1] -> Y[2]*sin(phi))/4.4482216/0.87);
      for (i = 0 ; i < problem -> solution -> num_output_nodes ; i++) {
         DisplayMessage("node %d depth=%.1f", problem -> solution -> output_nodes[i], environment -> depth - active[problem -> solution -> output_nodes[i]] -> x);
      }
   }
   else if (problem -> type == Drifter) {
      T = Tension(active[num_active] -> Y[1], active [num_active] -> material);
      T1 = Tension(active[1] -> Y[1], active [1] -> material);
      phi = active[num_active] -> Y[3];
      phi1 = active[1] -> Y[3];

      b = problem -> terminal [2] -> buoy;
      b1 = problem -> terminal [1] -> buoy;

      Current(0.0, environment -> depth, 0, 0, &uack, &vack, NULL);
      drag = 0.5*environment -> rho*ProjectedArea(b, b -> draft)*b -> Cdn*(vack - problem -> terminal[2] -> yspeed.value)*fabs(vack - problem -> terminal[2] -> yspeed.value); 
      // fprintf(stderr,"draft = %g, vack = %g, drag = %g\n", b -> draft, vack, drag);

      Current(0.0, active[1] -> x, 0, 0, &uack, &vack, NULL);
      drag1 = 0.5*environment -> rho*b1 -> S*b1 -> Cdn*(vack - problem -> terminal[2] -> yspeed.value)*fabs(vack - problem -> terminal[2] -> yspeed.value); 
      DisplayMessage("Done. Speed= %.2f, Sink H= %.0f, Sink Fx= %.2f, Sink drag= %.2f, Buoy Fx = %.2f, Buoy drag = %.2f", 
               problem -> terminal [2] -> yspeed.value,
               environment -> surface - active [1] -> x,
               T1*sin(phi1) + active[1] -> Y[2]*cos(phi1), drag1,
               T*sin(phi) + active[num_active] -> Y[2]*cos(phi), drag);
      fprintf(stdout ? stdout : stderr, "--> Done. Speed= %.2f, Sink H= %.0f, Sink Fx= %.2f, Sink drag= %.2f, Buoy Fx = %.2f, Buoy drag = %.2f\n", 
               problem -> terminal [2] -> yspeed.value,
               environment -> surface - active [1] -> x,
               T1*sin(phi1) + active[1] -> Y[2]*cos(phi1), drag1,
               T*sin(phi) + active[num_active] -> Y[2]*cos(phi), drag);
   }

   return 0;
}
