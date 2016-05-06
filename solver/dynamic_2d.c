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
 * File:	dynamic_2d.c
 *
 * Description: contains code specific to the solution of the 2D dynamic model
 *
 * History:
 *		
 ****************************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <malloc.h>
# include <stdlib.h>
# include <string.h>
# include "compress.h"
# include "problem.h"
# include "tension.h"
# include "output.h"
# include "error.h"
# include "solve.h"
# include "control.h"
# include "segments.h"

# define SQR(a) ((a)*(a))

# define NE	6
# define NB	3
# define NJ_COMPAT 2
# define NB_BRANCH 1

extern Problem     *problem;             
extern Environment *environment;
extern Analysis    *analysis;
extern Debug debug;

# if (defined HAVELAMP && !defined API)
extern int UpdateLAMP(Node *, double **, int, double, ResFile, double, int);
# endif

static double	ak, omak;
static double	am, omam;
static double	omam2, amam2, am2, omg_g;

static void 
GetXYDot(Node a, double *xdot, double *ydot)
{
   double	u, v, u_o, v_o;
   double	ph_o, ph;
   double	xdot_o, ydot_o;
   double	cphk, cphok, sphk, sphok;

   u    = a -> Y [3];
   u_o  = a -> Y_o [3];
   v    = a -> Y [4];
   v_o  = a -> Y_o [4];
   ph   = a -> Y [5];
   ph_o = a -> Y_o [5];
 
   xdot_o = a -> xdot_o;
   ydot_o = a -> ydot_o;

   cphk  = cos(ph);
   sphk  = sin(ph);
   cphok = cos(ph_o);
   sphok = sin(ph_o);

   *xdot = (omak*(u*cphk - v*sphk) + ak*(u_o*cphok - v_o*sphok)
              - am*xdot_o)/omam;
   *ydot = (omak*(u*sphk + v*cphk) + ak*(u_o*sphok + v_o*cphok)
               - am*ydot_o)/omam;
   return;
}

static void 
IntegrateSpatialXY(Node start)
{
   Node		  a, a_m;
   double	  sf, sf_m, sf_o, sf_om;
   double	  cphk, cphkm, cphok, cphokm;
   double	  sphk, sphkm, sphok, sphokm;
   double	  ph;
   double	  ph_o;
   double	  ph_om;
   double	  ph_m;
   double	  ds;
   double     xdot, ydot;

   a = start;

   while (a && a -> active) {	
       a_m = a -> prev_active;
       ds = a_m -> ds / 2.0;

       sf    = (1.0 + a -> Y [1])*ds;
       sf_o  = (1.0 + a -> Y_o [1])*ds;
       sf_m  = (1.0 + a_m -> Y [1])*ds;
       sf_om = (1.0 + a_m -> Y_o [1])*ds;
     
       ph    = a -> Y [5];
       ph_o  = a -> Y_o [5];
       ph_m  = a_m -> Y [5];
       ph_om = a_m -> Y_o [5];

       cphk   = cos(ph);
       cphkm  = cos(ph_m);
       cphok  = cos(ph_o);
       cphokm = cos(ph_om);

       sphk   = sin(ph);
       sphkm  = sin(ph_m);
       sphok  = sin(ph_o);
       sphokm = sin(ph_om);
       
       a -> x = a_m -> x 
                + (omak*(sf*cphk + sf_m*cphkm) 
	          + ak*(sf_o*cphok + sf_om*cphokm)
                  - ak*(a -> x_o - a_m -> x_o))/omak;
       a -> y = a_m -> y 
                + (omak*(sf*sphk + sf_m*sphkm) 
 	        + ak*(sf_o*sphok + sf_om*sphokm)
                - ak*(a -> y_o - a_m -> y_o))/omak;
       a -> z = 0.0;

       GetXYDot(a, &xdot, &ydot);
       a -> xdot = xdot;
       a -> ydot = ydot;

       a = a -> next_active;
    }

}

static void 
IntegrateTemporalXY(Node start, double dt, double xdot, double ydot)
{
   Node		  a, a_m;
   double	  cphk, cphkm, cphok, cphokm;
   double	  sphk, sphkm, sphok, sphokm;
   double	  u, v, ph;
   double	  u_o, v_o, ph_o;
   double	  u_om, v_om, ph_om;
   double	  u_m, v_m, ph_m;
   double	  xdot_m, ydot_m;
   double	  xdot_o, ydot_o;
   double	  xdot_om, ydot_om;
   double         gamma, gdt, dt1g;
   double	  ds, ds0;

   gamma = analysis -> gamma;
   gdt = gamma*dt;
   dt1g = dt*(1.0 - gamma);

   a = start;
   ds = ds0 = 0; 
   while (a && a -> active) {	
       a_m = a -> prev_active;

       ds0 += a_m -> ds;
       ds += a_m -> ds*(1.0 + 0.5*(a -> Y[1] + a_m -> Y[1]));

       xdot_m = xdot;
       ydot_m = ydot;
       
       u    = a -> Y [3];
       u_o  = a -> Y_o [3];
       u_m  = a_m -> Y [3];
       u_om = a_m -> Y_o [3];

       v    = a -> Y [4];
       v_o  = a -> Y_o [4];
       v_m  = a_m -> Y [4];
       v_om = a_m -> Y_o [4];

       ph    = a -> Y [5];
       ph_o  = a -> Y_o [5];
       ph_m  = a_m -> Y [5];
       ph_om = a_m -> Y_o [5];

       xdot_o = a -> xdot_o;
       ydot_o = a -> ydot_o;

       xdot_om = a_m -> xdot_o;
       ydot_om = a_m -> ydot_o;

       cphk   = cos(ph);
       cphkm  = cos(ph_m);
       cphok  = cos(ph_o);
       cphokm = cos(ph_om);

       sphk   = sin(ph);
       sphkm  = sin(ph_m);
       sphok  = sin(ph_o);
       sphokm = sin(ph_om);

       xdot = (omak*(u*cphk - v*sphk + u_m*cphkm - v_m*sphkm)
                  + ak*(u_o*cphok - v_o*sphok + u_om*cphokm - v_om*sphokm)
                  - am*(xdot_o + xdot_om))/omam - xdot_m;
       ydot = (omak*(u*sphk + v*cphk + u_m*sphkm + v_m*cphkm)
                  + ak*(u_o*sphok + v_o*cphok + u_om*sphokm + v_om*cphokm)
                  - am*(ydot_o + ydot_om))/omam - ydot_m;

       a -> x = a -> x_o + dt1g*xdot_o + gdt*xdot;
       a -> y = a -> y_o + dt1g*ydot_o + gdt*ydot;
       a -> z = 0.0;

       a -> xdot = xdot;
       a -> ydot = ydot;

       a = a -> next_active;
   }

   return;
}

void 
DynamicUpdate2D (
   Node		 *active,
   int		  num_active,
   double	  tm,
   double	  dt)
{
   int		  i;
   Node		  a;
   double	  gamma;
   double	  dt1g;
   double	  gdt;
   double	  xdot, ydot;
   Node		  from;

   gamma = analysis -> gamma;
   gdt = gamma*dt;
   dt1g = dt*(1.0 - gamma);

   a = active[1];

   if (problem -> type == Deployment || 
       problem -> type == Towing ||
       problem -> type == HorizontalDrifter ||
       problem -> type == Drifter ||
       (problem -> terminal [1] -> anchor && 
        problem -> terminal [1] -> anchor -> mu) ||
       (problem -> terminal [1] -> release > 0.0 && 
           tm >= problem -> terminal [1] -> release && !problem -> dynstat)) {

      GetXYDot(a, &xdot, &ydot);

      a -> x = a -> x_o + dt1g*a -> xdot_o + gdt*xdot;
      a -> y = a -> y_o + dt1g*a -> ydot_o + gdt*ydot;
      a -> z = 0.0; 
      a -> xdot = xdot;
      a -> ydot = ydot;
   }
   else {
      a -> xdot = xdot = 0.0;
      a -> ydot = ydot = 0.0;

      a -> x = problem -> terminal [1] -> x;
      a -> y = problem -> terminal [1] -> y;
      a -> z = 0.0;
   }


   if (analysis -> integration == Temporal) 
      IntegrateTemporalXY(a -> next_active, dt, xdot, ydot);
   else
      IntegrateSpatialXY(a -> next_active);

   for (i = 1 ; i <= problem -> num_branch ; i++) {
      a = problem -> branch [i] -> segment[1] -> first_active;
      from = problem -> branch [i] -> segment_from -> last_active;

      a -> x = from -> x; 
      a -> y = from -> y;
      a -> z = 0.0;
      GetXYDot(from, &xdot, &ydot);
      a -> xdot = xdot;
      a -> ydot = ydot;
      a -> zdot = 0;

      if (analysis -> integration == Temporal) 
         IntegrateTemporalXY(a -> next_active, dt, xdot, ydot);
      else 
         IntegrateSpatialXY(a -> next_active);
   }

   return;
}

void 
DynamicDifeq2D (
   EquationType	  eq_type,
   Node		  n,
   Node		  nm,
   int		  ne,
   int		  rhs, 
   int		  num_rows,
   double	**s,
   Node		 *node,
   double	  tm, 
   double	  dt,
   double	  current_factor)		/* not used		*/
{
   static double   ulcpre, vlcpre;
   static double   uacpre, vacpre;
   static double   ulcopre, vlcopre;
   static double   uacopre, vacopre;
   static double   uldpre, vldpre;
   static double   uldopre, vldopre;
   Node		   nj;	
   Buoy		   buoy;
   Anchor	   a;
   Connector	   c;
   double	   e, Sn, phi, Om3;
   double	   e_m, Sn_m, phi_m, Om3_m; 
   double          e_o, Sn_o, phi_o, Om3_o;
   double          e_om, Sn_om, phi_om, Om3_om;
   double	   e_j, Sn_j, phi_j;
   double          u, v;
   double          u_m, v_m;
   double	   cphk, sphk;
   double	   cphkm, sphkm;
   double	   cphkj, sphkj;
   double          u_j, v_j;
   double	   u_o, v_o;
   double	   u_om, v_om;
   double	   u_oj, v_oj, phi_oj;
   double	   cphok, sphok, cphokm, sphokm, cphokj, sphokj;
   double	   e_d,  phi_d, u_d, v_d;
   double	   e_d_m, phi_d_m, u_d_m, v_d_m;
   double	   e_d_o, phi_d_o, u_d_o, v_d_o;
   double	   e_d_om, phi_d_om, u_d_om, v_d_om;
   double          U, V, U_o, V_o;
   double	   Ur, Vr;
   double	   EI, w0, drat, drap, drat_o, drap_o;
   double	   w0_m, drat_m, drap_m, drat_om, drap_om;
   double	   damp_t;
   double	   damp_n;
   double	   ds;
   int		   i, j;
   int		   njn;
   double	   sign_j;
   double          tk, tkm, tkj, tdk, tdkm, tdkj, tddk, tddkm;
   double          tok, tokm, tdok, tdokm;
   double	   sqstrk, sqstrkm, sqstrok, sqstrokm;
   double	   elfak, elfakm, elfaok, elfaokm;
   double	   cam, camma_n, camma_t, caarma_n, caarma_t;
   double	   cam_m, camma_n_m, camma_t_m, caarma_n_m, caarma_t_m;
   double	   urk, vrk, urkm, vrkm;
   double	   ulck, vlck, ulckm, vlckm;
   double	   uldk, vldk, uldkm, vldkm;
   double	   ulcok, ulcokm, vlcok, vlcokm;
   double	   uldok, uldokm, vldok, vldokm;
   double 	   uack, vack, uackm, vackm;
   double	   uacok, vacok, uacokm, vacokm;
   double	   urok, urokm, vrok, vrokm;
   double	   urka, urkma, uroka, urokma;
   double	   vrka, vrkma, vroka, vrokma;
   double	   t_ds;
   double	   Uspd, Vspd;
   double	   pay, pay_m, pay_o, pay_om;
   double	   Uw, Vw, Uw_o, Vw_o, Uw_om, Vw_om, Uw_m, Vw_m;
   double	   Udw, Vdw, Udw_m, Vdw_m;
   double	   Udw_o, Vdw_o, Udw_om, Vdw_om;
   double          bmma, barma, bdr, bdr_n, bdr_t;
   double	   wet;
   double	   mma, dr_n, dr_t;
   double	   carma;
   double	   mud_b, mud_bm;
   double	   mud_bo, mud_bom;
   double	   mu;
   double	   bottom, bottom_m;
   double	   bottom_o, bottom_om;
   double	   Fb, Fb_m;
   double	   Fb_o, Fb_om;
   double	   B;
   double	   xforst, yforst;
   double	   xthrust, ythrust;
   double	   wave_b;
   double	   Fex [4];
   double	   nav;
   double	   wind_drag;
   double	   gdt;
   double	   omak2_ds, akak2_ds, ak2_ds;
   double          ftr;
   ForcingMethod   forcing;
   double	   Fv, Fh;
   double	   friction;

   forcing = environment -> forcing;

   EI     = n -> material -> EI;
   damp_n = n -> material -> bn;
   damp_t = n -> material -> bt;

   w0     = n -> material -> wet;
   cam    = n -> material -> m;
   camma_n  = cam + n -> material -> amn;
   camma_t  = cam + n -> material -> amt;
   caarma_n   = n -> material -> amn + n -> material -> rV;
   if (n -> material -> amt)
      caarma_t = n -> material -> amt + n -> material -> rV;
   else
      caarma_t = 0.0;

   e  = n -> Y [1];
   Sn = n -> Y [2];
   u  = n -> Y [3];
   v  = n -> Y [4];
   phi = n -> Y [5];
   Om3 = n -> Y [6];

   e_o  = n -> Y_o [1];
   Sn_o = n -> Y_o [2];
   u_o  = n -> Y_o [3];
   v_o  = n -> Y_o [4];
   phi_o = n -> Y_o [5];
   Om3_o = n -> Y_o [6];

   tk     = Tension(e, n -> material);
   tdk    = TensionD(e, n -> material);
   tddk   = TensionDD(e, n -> material);
   tok     = Tension(e_o, n -> material);
   tdok    = TensionD(e_o, n -> material);

   sqstrk  = sqrt(1.0 + e);
   elfak   = pow((1.0 + e), 3.0);
   sqstrok = sqrt(1.0 + e_o);
   elfaok  = pow((1.0 + e_o), 3.0);

   cphk = cos(phi);
   sphk = sin(phi);
   cphok = cos(phi_o);
   sphok = sin(phi_o);

   ds = nm -> ds;

	/* 
	 * for the straight material these will be the same as
	 * the next node, but we may need to modify them for 
	 * attachments and bottom interaction
	 */
   
   w0_m     = w0;
   cam_m    = cam;
   camma_t_m = camma_t;
   camma_n_m = camma_n;
   caarma_n_m  = caarma_n;
   caarma_t_m  = caarma_t;

   e_m = nm -> Y [1];
   Sn_m = nm -> Y [2];
   u_m  = nm -> Y [3];
   v_m  = nm -> Y [4];
   phi_m = nm -> Y [5];
   Om3_m = nm -> Y [6];

   e_om = nm -> Y_o [1];
   Sn_om = nm -> Y_o [2];
   u_om  = nm -> Y_o [3];
   v_om  = nm -> Y_o [4];
   phi_om = nm -> Y_o [5];
   Om3_om = nm -> Y_o [6];

   tkm    = Tension(e_m, nm -> material);
   tdkm   = TensionD(e_m, nm -> material);
   tddkm  = TensionDD(e_m, nm -> material);
   tokm    = Tension(e_om, nm -> material);
   tdokm   = TensionD(e_om, nm -> material);


   sqstrkm  = sqrt(1.0 + e_m);
   elfakm   = pow((1.0 + e_m), 3.0);
   sqstrokm = sqrt(1.0 + e_om);
   elfaokm  = pow((1.0 + e_om), 3.0);

   cphkm = cos(phi_m);
   sphkm = sin(phi_m);
   cphokm = cos(phi_om);
   sphokm = sin(phi_om);

   for (i = 1 ; i <= num_rows ; i++)
      for (j = 1 ; j <= rhs ; j++)
         s [i][j] = BLANK;

	/*
	 * anchor node
	 */

   if (eq_type == BottomBoundary) {
      s [1][6] = 1.0;
      s [1][rhs] = Om3;

	/*
	 * zero force case
	 */

      // there are possible edge cases in which we would switch 
      // back and forth between the BC cases in the following chain
      // of ifs. Just in case this happens, we zero out all the entries
      // so that the nspiv map of active entries does not change

      for (i = 1 ; i <= 5 ; i++) {
         s [2][i] = 0.0;
         s [3][i] = 0.0;
      }

      if (problem -> terminal [1] -> release > 0.0 && 
          tm >= problem -> terminal [1] -> release && !problem -> dynstat) {


         s [2][1] = 1.0;
         s [2][rhs] = e;

         s [3][2] = 1.0;
         s [3][rhs] = Sn;
      }

	/*
	 * anchored problems where the anchor might drag
	 */

      else if (problem -> terminal [1] -> anchor &&
               problem -> terminal [1] -> anchor -> mu) {

          Fv = tk*cphk - Sn*sphk;
          Fh = tk*sphk + Sn*cphk;
      
          a    = problem -> terminal [1] -> anchor;
          bmma = a -> m + a -> am;
          mu   = a -> mu;
          wet  = a -> wet;
 
          friction = -sign(Fh)*mu*(wet - Fv);

          if (fabs(friction) > fabs(Fh)) {		/* no motion */
             s [2][3] = 1.0;	
             s [2][rhs] = u;

             s [3][4] = 1.0;
             s [3][rhs] = v;
          }
          else {
             V   = u*sphk + v*cphk;
             V_o = u_o*sphok + v_o*cphok;
             U   = u*cphk - v*sphk;
/*
             fprintf(stderr,"Fh = %7.1f, Fv = %7.1f, friction = %7.1f, V = %4.2f, V_o = %4.2f, U = %4.2f\n", Fh, Fv, friction, V, V_o, U);
*/
             s [2][1]   = -tdk*sphk;
             s [2][2]   = -cphk;
             s [2][3]   = bmma/dt*sphk;
             s [2][4]   = bmma/dt*cphk;
             s [2][5]   = bmma/dt*(u*cphk - v*sphk) - tk*cphk + Sn*sphk;
             s [2][rhs] = bmma/dt*(V - V_o) - Fh - friction;

             s [3][3]  = cphk;
             s [3][4]  = -sphk;
             s [3][5]  = -u*sphk - v*cphk;
             s [3][rhs] = u*cphk - v*sphk;
          }
      }

	/*
	 * zero velocity cases
	 */

      else if ((problem -> dynstat && problem -> terminal [1] -> anchor) ||
               (problem -> terminal [1] -> anchor && problem -> type != Deployment) ||
               (problem -> type == Deployment && n -> x <= Bottom(n -> y, n -> z, tm)) ||
               (problem -> type == Towing && environment -> depth 
                && n -> x <= Bottom(n -> y, n -> z, tm))) {


         s [2][3] = 1.0;	
         s [2][rhs] = u;

         s [3][4] = 1.0;
         s [3][rhs] = v;
      }


	/*
	 * specified velocity (ROV maneuvering most likely)
	 */
	
      else if (problem -> type == Towing && 
               ((problem -> terminal [1] -> xspeed.value || 
                 problem -> terminal [1] -> xspeed.expr) ||
                (problem -> terminal [1] -> yspeed.value ||
                 problem -> terminal [1] -> yspeed.expr))) {


         Speed (tm, problem -> terminal [1], &Uspd, &Vspd, NULL);

         s [2][3]   = 1.0;
         s [2][5]   = Uspd*sphk - Vspd*cphk;
         s [2][rhs] = u - Uspd*cphk - Vspd*sphk;

         s [3][4]   = 1.0;
         s [3][5]   = Uspd*cphk + Vspd*sphk;
         s [3][rhs] = v + Uspd*sphk - Vspd*cphk;
      }
      else {
         if (problem -> type == Deployment) {
            a = problem -> terminal [1] -> anchor;
            bmma = a -> m + a -> am;
            barma = 0.0;
            bdr_n = a -> Cdn;
            bdr_t = a -> Cdt;
            wet = a -> wet;

            Uw = Vw = 0.0;
            Udw = Vdw = 0.0;

            xthrust = 0.0;
            ythrust = 0.0;
         }
         else {
            buoy = problem -> terminal [1] -> buoy;

            bmma = buoy -> Mmma;
            barma = buoy -> Marma;
            bdr_n = buoy -> Mdr;
            bdr_t = buoy -> Mdr*buoy -> Cdt/buoy -> Cdn;

            wet = buoy -> w - buoy -> buoyancy + n -> segment -> bottom_wet;

            if (forcing != Velocity && forcing != Force) {

               WaveParticleVelocity (tm, n -> x, n -> y, 0.0, &Uw, &Vw, NULL);
               WaveParticleAcceleration (tm, n -> x, n -> y, 0.0, 
                                         &Udw, &Vdw, NULL);
            }
            else {
               Uw = Vw = 0.0;
               Udw = Vdw = 0.0;
            }

            if (problem -> type == Towing)  {
               Thrust(tm, problem -> terminal [1], &xthrust, &ythrust, NULL);
               problem -> terminal[1] -> xthrust.value = xthrust;
               problem -> terminal[1] -> ythrust.value = ythrust;
            }
            else {
               xthrust = 0.0;
               ythrust = 0.0;
            } 
         }

         U   = u*cphk - v*sphk;
         U_o = u_o*cphok - v_o*sphok;
         V   = u*sphk + v*cphk;
         V_o = u_o*sphok + v_o*cphok;

         Current (tm, n -> x, n -> y, n -> z, &uack, &vack, NULL);

         Ur   = U - Uw;
         Vr   = V - vack - Vw;

         if (problem -> type == HorizontalDrifter) { 
            s [2][3]   = cphk;
            s [2][4]   = -sphk;
            s [2][5]   = -u*sphk - v*cphk;
            s [2][rhs] = u*cphk - v*sphk;
         }
         else {
            ftr = bmma/dt + 2.0*bdr_t*sign(Ur)*Ur;
            s [2][1] = -tdk*cphk;
            s [2][2] = sphk;
            s [2][3] = cphk*ftr;
            s [2][4] = -sphk*ftr;
            s [2][5] = (-u*sphk - v*cphk)*ftr + Sn*cphk + tk*sphk;
            s [2][rhs] = bmma/dt*(U - U_o)
                - barma*Udw - tk*cphk + Sn*sphk 
                + wet + bdr_t*Ur*fabs(Ur) - xthrust;
         }
         ftr = bmma/dt + 2.0*bdr_n*sign(Vr)*Vr;
         s [3][1] = -tdk*sphk;
         s [3][2] = -cphk;
         s [3][3] = sphk*ftr;
         s [3][4] = cphk*ftr;
         s [3][5] = (u*cphk - v*sphk)*ftr + Sn*sphk - tk*cphk;
         s [3][rhs] = bmma/dt*(V - V_o)
                - barma*Vdw
                - tk*sphk - Sn*cphk
                + bdr_n*Vr*fabs(Vr) - ythrust;
      }
   }

	/* 	
	 * buoy node (or second anchor node)
	 */

   else if (eq_type == TopBoundary) {
      for (i = 1 ; i <= 5 ; i++) {
         s [1][i] = 0.0;
	     s [2][i] = 0.0;
      }
  
      if (problem -> terminal [2] -> release > 0.0 && 
          tm >= problem -> terminal [2] -> release && !problem -> dynstat) {

         s [1][1] = 1.0;
         s [1][rhs] = e;
        
         s [2][2] = 1.0;
         s [2][rhs] = Sn;
      }
      else if (problem -> terminal [2] -> anchor) {

         s [1][3] = 1.0;
         s [1][rhs] = u;

         s [2][4] = 1.0;
         s [2][rhs] = v;
      }
      else if (problem -> type == Towing) {

         Speed (tm, problem -> terminal [2], &Uspd, &Vspd, NULL);
         InputVelocity (tm, n -> x, n -> y, 0.0, &Uw, &Vw, NULL);

         s [1][3]   = 1.0;
         s [1][5]   = (Uspd + Uw)*sphk - (Vspd + Vw)*cphk;
         s [1][rhs] = u - (Uspd + Uw)*cphk - (Vspd + Vw)*sphk;

         s [2][4]   = 1.0;
         s [2][5]   = (Uspd + Uw)*cphk + (Vspd + Vw)*sphk;
         s [2][rhs] = v + (Uspd + Uw)*sphk - (Vspd + Vw)*cphk;
      }
      else if (forcing == Force) { 
         InputVelocity (tm, n -> x, n -> y, 0.0, &Uw, &Vw, NULL);

         yforst = problem -> terminal [2] -> yforce + Vw;
         xforst = problem -> terminal [2] -> xforce + Uw;

         s [1][1] = -tdk*cphk;
         s [1][2] = sphk;
         s [1][5] = tk*sphk + Sn*cphk;
         s [1][rhs] = xforst - (tk*cphk - Sn*sphk);

         s [2][1] = -tdk*sphk;
         s [2][2] = -cphk;
         s [2][5] = -tk*cphk + Sn*sphk;
         s [2][rhs] = yforst - (tk*sphk + Sn*cphk);
      }
      else if (forcing == Morison || problem -> type == Deployment
               || problem -> type == HorizontalDrifter
               || ((problem -> type == Surface || problem -> type == Subsurface) 
                   && problem -> dynstat)) {

         buoy = problem -> terminal [2] -> buoy;
         bmma = buoy -> Mmma;

         if (problem -> type == Deployment || problem -> dynstat) {
            buoy -> draft = environment -> surface - n -> x;
            barma = 0.0;
            bdr     =  0.5*environment -> rho*buoy -> Cdn*
                       ProjectedArea(buoy, buoy -> draft);

            Uw = Vw = 0.0;
            Udw = Vdw = 0.0;
         }
         else {
            buoy -> draft = buoy -> max_draft; // added for drifter work
            bdr = buoy -> Mdr;
            barma = buoy -> Marma;

            WaveParticleVelocity (tm, n -> x, n -> y, 0.0, &Uw, &Vw, NULL);
            WaveParticleAcceleration (tm, n -> x, n -> y, 0.0, 
				      &Udw, &Vdw, NULL);
         }


         if (problem -> type == Surface || problem -> type == Drifter) {
	        WindDrag (tm, buoy, &wind_drag, NULL);
            Current (tm, environment->surface, n -> y, n -> z, &uack, &vack, NULL); 
         }
         else {
            wind_drag = 0.0;
            Current (tm, n -> x, n->y, n->z, &uack, &vack, NULL);
         }

         U   = u*cphk - v*sphk;
         U_o = u_o*cphok - v_o*sphok;
         V   = u*sphk + v*cphk;
         V_o = u_o*sphok + v_o*cphok;

         Ur   = U - Uw;
         Vr   = V - vack - Vw;

         // B = Buoyancy(buoy, environment -> depth - n -> x);
         B = Buoyancy(buoy, buoy -> draft, environment); // changed for drifter work

         if (problem -> type == HorizontalDrifter) {
            s [1][3]   = cphk;
            s [1][4]   = -sphk;
            s [1][5]   = -u*sphk - v*cphk;
            s [1][rhs] = u*cphk - v*sphk;
         }
         else {
            ftr = bmma/dt + 2.0*bdr*sign(Ur)*Ur;

            s [1][1] = tdk*cphk;
            s [1][2] = -sphk;
            s [1][3] = cphk*ftr;
            s [1][4] = -sphk*ftr;
            s [1][5] = (-u*sphk - v*cphk)*ftr - tk*sphk - Sn*cphk;
            s [1][rhs] = bmma*(U - U_o)/dt
                      - barma*Udw
                      + bdr*Ur*fabs(Ur)
                      + tk*cphk - Sn*sphk
                      - (B - buoy -> w);
         } 
         ftr = bmma/dt + 2.0*bdr*sign(Vr)*Vr;

         s [2][1] = tdk*sphk;
         s [2][2] = cphk;
         s [2][3] = sphk*ftr;
         s [2][4] = cphk*ftr;
         s [2][5] = (u*cphk - v*sphk)*ftr - Sn*sphk + tk*cphk;
         s [2][rhs] = bmma/dt*(V - V_o)
                - barma*Vdw
                + tk*sphk + Sn*cphk
 		        - wind_drag
                + bdr*Vr*fabs(Vr);

/*
          fprintf (stderr,"draft = %g, Fc = %g, Fcurr + Fw = %g", buoy -> draft,
		   tk*sphk + Sn*cphk, wind_drag - bdr*Vr*fabs(Vr));
*/
      }
      else if (forcing == WaveFollower || forcing == Velocity || forcing == LAMP) {


         if (forcing == WaveFollower)
            WaveSurfaceVelocity (tm, n -> x, n -> y, 0.0, &Uw, &Vw, NULL);
         else 
            InputVelocity (tm, n -> x, n -> y, 0.0, &Uw, &Vw, NULL);
         
	     if (problem -> type == Towing) {	// don't use speed for drifters
            Speed (tm, problem -> terminal [2], &Uspd, &Vspd, NULL);
            Uw += Uspd;
            Vw += Vspd;
   	     }

         s [1][3] = -cphk;
         s [1][4] = sphk;
         s [1][5] = u*sphk + v*cphk;
         s [1][rhs] = Uw - (u*cphk - v*sphk);

/* force BC
         s [1][1] = -tdk*cphk;
         s [1][2] = sphk;
         s [1][3] = 0.0;
         s [1][4] = 0.0;
         s [1][5] = tk*sphk + Sn*cphk;
         s [1][rhs] = problem -> terminal[2] -> xforce - (tk*cphk - Sn*sphk);
*/

         if (forcing == WaveFollower && 
             (problem -> type == Surface || problem -> type == Drifter)) {

	        buoy = problem -> terminal [2] -> buoy;

            // buoy -> draft = environment -> surface - n -> x;
            buoy -> draft = Draft(tk*cphk - Sn*sphk, buoy, environment);

	        bdr = 0.5*environment -> rho*buoy -> Cdn*ProjectedArea(buoy, buoy -> draft);
	        bmma = buoy -> Mmma;
 
	        WindDrag (tm, buoy, &wind_drag, NULL);

            Current (tm, environment->surface - buoy->draft, n -> y, n -> z, 
                     &uack, &vack, NULL);

            V   = u*sphk + v*cphk;
            V_o = u_o*sphok + v_o*cphok;

            Vr   = V - vack;

            ftr = bmma/dt + 2.0*bdr*sign(Vr)*Vr;
            s [2][1] = tdk*sphk;
            s [2][2] = cphk;
            s [2][3] = sphk*ftr;
            s [2][4] = cphk*ftr;
            s [2][5] = (u*cphk - v*sphk)*ftr - Sn*sphk + tk*cphk;
            s [2][rhs] = bmma/dt*(V - V_o)
		                 + tk*sphk + Sn*cphk
                         + bdr*Vr*fabs(Vr)
		                 - wind_drag;
         }
         else {
            s [2][3] = -sphk;
            s [2][4] = -cphk;
            s [2][5] = -u*cphk + v*sphk;
            s [2][rhs] = Vw - (u*sphk + v*cphk);
         }
      }
      else if (forcing == FroudeKrylov) {
         buoy = problem -> terminal [2] -> buoy;

         bmma = buoy -> m + buoy -> am;
         bdr  =  0.5*environment -> rho*buoy -> Cdn*
                        ProjectedArea(buoy, environment -> surface - n -> x);

         FroudeKrylovCoefficients (tm, buoy, n -> x, n -> y, 0.0,
                                   Fex, &wave_b, &B);
         
         WaveParticleVelocity (tm, n -> x, n -> y, 0.0, &Uw, &Vw, NULL);

         Current (tm, n -> x, n -> y, n -> z, &uack, &vack, NULL);

         U   = u*cphk - v*sphk;
         U_o = u_o*cphok - v_o*sphok;
         V   = u*sphk + v*cphk;
         V_o = u_o*sphok + v_o*cphok;

         Ur   = U - Uw;
         Vr   = V - vack - Vw;

         ftr = bmma/dt + 2.0*bdr*sign(Ur)*Ur + wave_b;
        
         s [1][1] = tdk*cphk;
         s [1][2] = -sphk;
         s [1][3] = cphk*ftr;
         s [1][4] = -sphk*ftr;
         s [1][5] = -tk*sphk - Sn*cphk + (-u*sphk - v*cphk)*ftr;
         s [1][rhs] = bmma*(U - U_o)/dt
                      - Fex [1] + wave_b*U + bdr*Ur*fabs(Ur)
                      + (tk*cphk - Sn*sphk) - B;

         s [2][3]    = sphk;
         s [2][4]    = cphk;
         s [2][5]    = u*cphk - v*sphk;
         s [2][rhs]  = u*sphk + v*cphk;
      }

      s [3][6]   = 1.0;
      s [3][rhs] = Om3;
   }
   else if (eq_type == BranchStart) {
      s [1][6]   = 1.0;
      s [1][rhs] = Om3;
   }
   else if (eq_type == BranchTerminal) {

      s [1][6]   = 1.0;
      s [1][rhs] = Om3;

      if (n -> segment -> branch -> terminal -> buoy) {
 
         buoy = n -> segment -> branch -> terminal -> buoy;

         bmma    = buoy -> Mmma; 
         bdr     = buoy -> Mdr;

         if (forcing != Velocity && forcing != Force) {
            barma   = buoy -> Marma;
            WaveParticleVelocity (tm, n -> x, n -> y, 0.0, &Uw, &Vw, NULL);
            WaveParticleAcceleration (tm, n -> x, n -> y, 0.0, &Udw, &Vdw,NULL);
         }
         else {
            barma = 0.0;
            Uw = Vw = 0.0;
            Udw = Vdw = 0.0;
         }

         Current (tm, n -> x, n -> y, n -> z, &uack, &vack, NULL);

         U   = u*cphk - v*sphk;
         U_o = u_o*cphok - v_o*sphok;
         V   = u*sphk + v*cphk;
         V_o = u_o*sphok + v_o*cphok;

         Ur   = U - Uw;
         Vr   = V - vack - Vw;

         B = Buoyancy(buoy, environment -> depth - n -> x, environment);

         ftr = bmma/dt + 2.0*bdr*sign(Ur)*Ur;

         if (problem -> type == HorizontalDrifter) {
            s [2][3]   = cphk;
            s [2][4]   = -sphk; 
            s [2][5]   = -u*sphk - v*cphk;
            s [2][rhs] = u*cphk - v*sphk;
         }
         else { 
            s [2][1] = tdk*cphk;
            s [2][2] = -sphk;
            s [2][3] = cphk*ftr;
            s [2][4] = -sphk*ftr;
            s [2][5] = (-u*sphk - v*cphk)*ftr - tk*sphk - Sn*cphk;
            s [2][rhs] = bmma*(U - U_o)/dt
                      - barma*Udw
                      + bdr*Ur*fabs(Ur)
                      + (tk*cphk - Sn*sphk)
                      - (B - buoy -> w);
         }

         ftr = bmma/dt + 2.0*bdr*sign(Vr)*Vr;

         s [3][1] = tdk*sphk;
         s [3][2] = cphk;
         s [3][3] = sphk*ftr;
         s [3][4] = cphk*ftr;
         s [3][5] = (u*cphk - v*sphk)*ftr - Sn*sphk + tk*cphk;
         s [3][rhs] = bmma/dt*(V - V_o)
                - barma*Vdw
                + tk*sphk + Sn*cphk
                + bdr*Vr*fabs(Vr);
      }
      else if (n -> segment -> branch -> terminal -> anchor) {	
         s [2][3] = 1.0;
         s [2][rhs] = u;

         s [3][4] = 1.0;
         s [3][rhs] = v;
      }
      else { // if (n -> segment -> branch -> terminal -> loop_main_node) {
	     /* empty -- compatibility already expressed at Junction */
      }
   }

	/*
	 * node between two distinct segments
	 */

   else if (eq_type == Connection || eq_type == Junction) {
      
      if (nm -> segment -> connector == NULL 
          && nm -> segment -> connection == Spliced) {	/* no lumped mass 	  */
         s [1][1]      = -tdkm;
         s [1][ne + 1] = tdk;
         s [1][rhs]    = tk - tkm;

         for (i = 2 ; i <= ne ; i++) {
            s [i][i] = -1.0;
            s [i][i + ne] = 1.0;
            s [i][rhs] = n -> Y[i] - nm -> Y[i];
         }
      }
      else {				/* pinned connection lumped mass may be present */
         c = nm -> segment -> connector;
         if (c) {
             mma = c -> m + c -> am;
             dr_n = c -> Cdn;
             dr_t = c -> Cdt;

             wet = c -> wet + nm -> segment -> top_wet + n -> segment -> bottom_wet;
             carma   = 1.5*pow(c -> d, 3.0)*M_PI*environment -> rho/6.0;
         }
         else {
             mma = dr_n = dr_t = wet = carma = 0.0;
         }

         if (problem -> type == Deployment  || problem -> type == Surface)
             if (n -> x > environment -> surface && wet < 0.0)
                wet = wet*(1.0 + tanh(50.0*(environment -> surface - n -> x)));

         njn = 0;

         s [1][6]   = 1.0;
         s [1][rhs] = Om3_m;

         s [2][ne + 6] = 1.0;
         s [2][rhs]    = Om3;

         nav = 2.0;

         u_om   = nm -> Y_o [3];
         v_om   = nm -> Y_o [4];
         phi_om = nm -> Y_o [5];

         cphokm = cos(phi_om);
         sphokm = sin(phi_om);

         U   = (u_m*cphkm - v_m*sphkm + u*cphk - v*sphk)/nav;
         U_o = (u_om*cphokm - v_om*sphokm + u_o*cphok - v_o*sphok)/nav;
         V   = (u_m*sphkm + v_m*cphkm + u*sphk + v*cphk)/nav;
         V_o = (u_om*sphokm + v_om*cphokm  + u_o*sphok + v_o*cphok)/nav;

         if (eq_type == Junction) {
            njn = n -> segment -> junction.num_nodes;

            U   = nav*U;
            U_o = nav*U_o;
            V   = nav*V;
            V_o = nav*V_o;

            for (i = 1 ; i <= njn ; i++) {
               nj = n -> segment -> junction.node [i];

               u_j    = nj -> Y [3];
               v_j    = nj -> Y [4];
               u_oj   = nj -> Y_o [3];
               v_oj   = nj -> Y_o [4];

               phi_j  = nj -> Y [5];
               phi_oj = nj -> Y_o [5];

               cphkj = cos(phi_j);
               sphkj = sin(phi_j);
               cphokj = cos(phi_oj);
               sphokj = sin(phi_oj);

               U   += u_j*cphkj - v_j*sphkj;
               U_o += u_oj*cphokj - v_oj*sphokj;
               V   += u_j*sphkj + v_j*cphkj;
               V_o += u_oj*sphokj + v_oj*cphokj;
            }

            nav += njn;

            U   = U / nav;
            U_o = U_o / nav;
            V   = V / nav;
            V_o = V_o / nav;
         }

         Current (tm, n -> x, n -> y, n -> z, &uack, &vack, NULL);
         if (problem -> type == Deployment && problem -> dynstat) 
            vack -= problem -> terminal [1] -> yspeed.value;
        
         ConnectorThrust(tm, nm -> segment, nm, &xthrust, &ythrust, NULL); 
         nm -> segment -> connector_xthrust.value = xthrust;
         nm -> segment -> connector_ythrust.value = ythrust;

         if (forcing != Velocity && forcing != Force) {
            WaveParticleVelocity (tm, n -> x, n -> y, 0.0, &Uw, &Vw, NULL);

            WaveParticleAcceleration (tm, n -> x, n -> y, 0.0, 
                                      &Udw, &Vdw, NULL);

         }
         else {
            carma = 0.0;
            Uw = Vw = 0.0;
            Udw = Vdw = 0.0;
         }

         Ur   = U - Uw;
         Vr   = V - vack - Vw;

	/*
 	 * sum of forces
	 */

         if (problem -> type == HorizontalDrifter) {
            s [3][3]      = cphkm;
            s [3][4]      = -sphkm;
            s [3][5]      = -u_m*sphkm - v_m*cphkm;
            s [3][rhs]    = u_m*cphkm - v_m*sphkm;
         }
         else {
            ftr = (mma/dt + 2.0*dr_t*sign(Ur)*Ur)/nav;

            s [3][1]        = tdkm*cphkm;
            s [3][ne + 1]   = -tdk*cphk;
            s [3][2]        = -sphkm;
            s [3][ne + 2]   = sphk;
            s [3][3]        = cphkm*ftr;
            s [3][ne + 3]   = cphk*ftr;
            s [3][4]        = -sphkm*ftr;
            s [3][ne + 4]   = -sphk/nav*ftr;
            s [3][5]        = ftr*(-u_m*sphkm - v_m*cphkm) 
                           - Sn_m*cphkm - tkm*sphkm;
            s [3][ne + 5]   = ftr*(-u*sphk - v*cphk) + Sn*cphk + tk*sphk;
            s [3][rhs] = mma*(U - U_o)/dt
               - carma*Udw
               + (-Sn_m*sphkm + tkm*cphkm + Sn*sphk - tk*cphk)
               + dr_t*Ur*fabs(Ur) 
               + wet - xthrust;
         }
        
         ftr = (mma/dt + 2.0*dr_n*sign(Vr)*Vr)/nav;

         s [4][1]        = tdkm*sphkm;
         s [4][ne + 1]   = -tdk*sphk;
         s [4][2]        = cphkm;
         s [4][ne + 2]   = -cphk;
         s [4][3]        = sphkm*ftr;
         s [4][ne + 3]   = sphk*ftr;
         s [4][4]        = cphkm*ftr;
         s [4][ne + 4]   = cphk*ftr;
         s [4][5]        = ftr*(u_m*cphkm - v_m*sphkm) + tkm*cphkm - Sn_m*sphkm;
         s [4][ne + 5]   = ftr*(u*cphk - v*sphk) - tk*cphk + Sn*sphk;
         s [4][rhs] = mma/dt*(V - V_o)
                - carma*Vdw - ythrust
                + (Sn_m*cphkm + tkm*sphkm - Sn*cphk - tk*sphk)
                + dr_n*Vr*fabs(Vr);

	/*
 	 * compatibility of nodes on the main line
	 */

         s [5][3]      = cphkm;
         s [5][ne + 3] = -cphk;
         s [5][4]      = -sphkm;
         s [5][ne + 4] = sphk;
         s [5][5]      = -u_m*sphkm - v_m*cphkm;
         s [5][ne + 5] = u*sphk + v*cphk;
         s [5][rhs]    = u_m*cphkm - v_m*sphkm 
                         - u*cphk + v*sphk;

         s [6][3]      = sphkm;
         s [6][ne + 3] = -sphk;
         s [6][4]      = cphkm;
         s [6][ne + 4] = -cphk;
         s [6][5]      = u_m*cphkm - v_m*sphkm;
         s [6][ne + 5] = -u*cphk + v*sphk;
         s [6][rhs]    = u_m*sphkm + v_m*cphkm 
                         - u*sphk - v*cphk;

         j = 7;

         for (i = 1 ; i <= njn ; i++) {
            nj = n -> segment -> junction.node [i];
        
            e_j    = nj -> Y [1];
            Sn_j   = nj -> Y [2];

            u_j    = nj -> Y [3];
            v_j    = nj -> Y [4];

            tkj   = Tension(e_j, nj -> material);
            tdkj  = TensionD(e_j, nj -> material);
 
            phi_j  = nj -> Y [5];

            cphkj = cos(phi_j);
            sphkj = sin(phi_j);

            if (nj == nj -> segment -> branch -> last) 
               sign_j = -1.0;
            else
               sign_j = 1.0;

	/*
 	 * add any force terms for this junction node
	 */

            ftr = (mma/dt + 2.0*dr_t*sign(Ur)*Ur)/nav;
            s [3][(i + 1)*ne + 1] = -sign_j*cphkj*tdkj;
            s [3][(i + 1)*ne + 2] = sign_j*sphkj;
            s [3][(i + 1)*ne + 3] = cphkj*ftr;
            s [3][(i + 1)*ne + 4] = -sphkj*ftr;
            s [3][(i + 1)*ne + 5] = ftr*(-u_j*sphkj - v_j*cphkj)
                              + sign_j*(Sn_j*cphkj + tkj*sphkj);
            s [3][rhs] += sign_j*(Sn_j*sphkj - tkj*cphkj);

            ftr = (mma/dt + 2.0*dr_n*sign(Vr)*Vr)/nav;
            s [4][(i + 1)*ne + 1] = -sign_j*tdkj*sphkj;
            s [i][(i + 1)*ne + 2] = -sign_j*cphkj;
            s [4][(i + 1)*ne + 3] = sphkj*ftr;
            s [4][(i + 1)*ne + 4] = cphkj*ftr;
            s [4][(i + 1)*ne + 5] = ftr*(u_j*cphkj - v_j*sphkj)
                               + sign_j*(Sn_j*sphkj - tkj*cphkj);
            s [4][rhs] += sign_j*(-Sn_j*cphkj - tkj*sphkj);

	/*
 	 * compatibility of any junction nodes
	 */
 
            s [j][ne + 3]         = cphk;
            s [j][(i + 1)*ne + 3] = -cphkj;
            s [j][ne + 4]         = -sphk;
            s [j][(i + 1)*ne + 4] = sphkj;
            s [j][ne + 5]         = -u*sphk - v*cphk;
            s [j][(i + 1)*ne + 5] = u_j*sphkj + v_j*cphkj;
            s [j][rhs] = u*cphk - v*sphk - u_j*cphkj + v_j*sphkj;
            j ++;

            s [j][ne + 3]         = sphk;
            s [j][(i + 1)*ne + 3] = -sphkj;
            s [j][ne + 4]         = cphk;
            s [j][(i + 1)*ne + 4] = -cphkj;
            s [j][ne + 5]         = u*cphk - v*sphk;
            s [j][(i + 1)*ne + 5] = -u_j*cphkj + v_j*sphkj;
            s [j][rhs] = u*sphk + v*cphk - u_j*sphkj - v_j*cphkj;
            j ++;
         }
      }
   }

	/*
	 * a regular old internal node, eq_type == Cable, SegmentTop
	 */

   else {
      gdt = analysis -> gamma*dt;
      t_ds = 2.0/ds;


      omak2_ds = omak*omak/ds;
      akak2_ds = ak*omak/ds;
      ak2_ds   = ak*ak / ds;

      e_d_o   = n -> Yd_o[1];
      // Sn_d_o  = n -> Yd_o[2];
      u_d_o   = n -> Yd_o[3];
      v_d_o   = n -> Yd_o[4];
      phi_d_o = n -> Yd_o[5];
 
      e_d_om   = nm -> Yd_o[1];
      /// Sn_d_om  = nm -> Yd_o[2];
      u_d_om   = nm -> Yd_o[3];
      v_d_om   = nm -> Yd_o[4];
      phi_d_om = nm -> Yd_o[5];

      e_d   = n -> Yd[1] = ((e - e_o)/gdt - omg_g*e_d_o);
      // Sn_d  = n -> Yd[2] = ((Sn - Sn_o)/gdt - omg_g*Sn_d_o);
      u_d   = n -> Yd[3] = ((u - u_o)/gdt - omg_g*u_d_o);
      v_d   = n -> Yd[4] = ((v - v_o)/gdt - omg_g*v_d_o);
      phi_d = n -> Yd[5] = ((phi - phi_o)/gdt - omg_g*phi_d_o);

      e_d_m   = nm -> Yd[1] = ((e_m - e_om)/gdt - omg_g*e_d_om);
      // Sn_d_m  = nm -> Yd[2] = ((Sn_m - Sn_om)/gdt - omg_g*Sn_d_om);
      u_d_m   = nm -> Yd[3] = ((u_m - u_om)/gdt - omg_g*u_d_om);
      v_d_m   = nm -> Yd[4] = ((v_m - v_om)/gdt - omg_g*v_d_om);
      phi_d_m = nm -> Yd[5] = ((phi_m - phi_om)/gdt - omg_g*phi_d_om);

      Current (tm, n -> x, n -> y, n -> z, &uack, &vack, NULL);
      Current (tm - dt, n -> x_o, n -> y_o, n -> z_o, &uacok, &vacok, NULL);
      if (problem -> type == Deployment && problem -> dynstat) {
         vack -= problem -> terminal [1] -> yspeed.value;
         vacok -= problem -> terminal [1] -> yspeed.value;
      }

      if (forcing == WaveFollower || forcing == Morison || forcing == LAMP) {
         WaveParticleVelocity (tm, n -> x, n -> y, 0.0, &Uw, &Vw, NULL);
         WaveParticleAcceleration (tm, n -> x, n -> y, 0.0, &Udw, &Vdw, NULL);

         WaveParticleVelocity (tm - dt, n -> x_o, n -> y_o, 0.0, &Uw_o, &Vw_o, NULL);
         WaveParticleAcceleration (tm - dt, n -> x_o, n -> y_o, 0.0, &Udw_o, &Vdw_o, NULL);

         uldk = cphk*Udw + sphk*Vdw;
         vldk = -sphk*Udw + cphk*Vdw;

         uldok = cphok*Udw_o + sphok*Vdw_o;
         vldok = -sphok*Udw_o + cphok*Vdw_o;
      }
      else {
         Uw = Vw = Uw_o = Vw_o = 0.0;
         Udw = Vdw = Udw_o = Vdw_o = 0.0;

         uldk = vldk = uldok = vldok = 0.0;
      }

      ulck = cphk*(uack + Uw) + sphk*(vack + Vw);
      vlck = -sphk*(uack + Uw) + cphk*(vack + Vw);

      ulcok = cphok*(uacok + Uw_o) + sphok*(vacok + Vw_o);
      vlcok = -sphok*(uacok + Uw_o) + cphk*(vacok + Vw_o);

      if (n -> active_number == 2 || nm -> active_number != n -> active_number - 1) {
         Current (tm, nm -> x, nm -> y, nm -> z, &uackm, &vackm, NULL);
         Current (tm - dt, nm -> x_o, nm -> y_o, nm -> z_o, &uacokm, &vacokm, NULL);

         if (problem -> type == Deployment && problem -> dynstat) {
            vackm -= problem -> terminal [1] -> yspeed.value;
            vacokm -= problem -> terminal [1] -> yspeed.value;
         }

         if (forcing == WaveFollower || forcing == Morison || forcing == LAMP) {
            WaveParticleVelocity (tm, nm -> x, nm -> y, 0.0, &Uw_m, &Vw_m, NULL);
            WaveParticleAcceleration (tm, nm -> x, nm -> y, 0.0, &Udw_m, &Vdw_m, NULL);

            WaveParticleVelocity (tm - dt, nm -> x_o, nm -> y_o, 0.0, &Uw_om, &Vw_om, NULL);
            WaveParticleAcceleration (tm - dt, nm -> x_o, nm -> y_o, 0.0,
			              &Udw_om, &Vdw_om, NULL);

            uldkm = cphkm*Udw_m + sphkm*Vdw_m;
            vldkm = -sphkm*Udw_m + cphkm*Vdw_m;

            uldokm = cphokm*Udw_om + sphokm*Vdw_om;
            vldokm = -sphokm*Udw_om + cphokm*Vdw_om;
         }
         else {
            Uw_m = Vw_m = Uw_om = Vw_om = 0.0;
            Udw_m = Vdw_m = Udw_om = Vdw_om = 0.0;
 
            uldkm = vldkm = uldokm = vldokm = 0.0;
         }

         ulckm = cphkm*(uackm + Uw_m) + sphkm*(vackm + Vw_m);
         vlckm = -sphkm*(uackm + Uw_m) + cphkm*(vackm + Vw_m);

         ulcokm = cphokm*(uacokm + Uw_om) + sphokm*(vacokm + Vw_om);
         vlcokm = -sphokm*(uacokm + Vw_om) + cphokm*(vacokm + Vw_om);
      }
      else {
         ulckm  = ulcpre;
         vlckm  = vlcpre;
         uackm  = uacpre;
         vackm  = vacpre;

         ulcokm  = ulcopre;
         vlcokm  = vlcopre;
         uacokm  = uacopre;
         vacokm  = vacopre;

         uldkm  = uldpre;
         vldkm  = vldpre;
           
         uldokm = uldopre;
         vldokm = vldopre;
      }

      ulcpre = ulck;
      vlcpre = vlck;
      uacpre = uack;
      vacpre = vack;

      uldpre = uldk;
      vldpre = vldk;

      ulcopre = ulcok;
      vlcopre = vlcok;
      uacopre = uacok;
      vacopre = vacok;

      uldopre = uldok;
      vldopre = vldok;

      urk  = u - ulck;
      vrk  = v - vlck;
      urkm = u_m - ulckm;
      vrkm = v_m - vlckm;

      urka = fabs(urk);
      vrka = fabs(vrk);
      urkma = fabs(urkm);
      vrkma = fabs(vrkm);

      urok = u_o - ulcok;
      vrok = v_o - vlcok;
      urokm = u_om - ulcokm;
      vrokm = v_om - vlcokm;

      uroka = fabs(urok);
      vroka = fabs(vrok);
      urokma = fabs(urokm);
      vrokma = fabs(vrokm);

      DragCoeff(n, n -> material, tm, urka, vrka, &drat, &drap);
      DragCoeff(nm, nm -> material, tm, urkma, vrkma, &drat_m, &drap_m);

      DragCoeff(n, n -> material, tm-dt, uroka, vroka, &drat_o, &drap_o);
      DragCoeff(nm, nm -> material, tm-dt, urokma, vrokma, &drat_om, &drap_om);

      if (n -> attachment) {
         w0    += n -> attachment -> wet/ds;
         drat  += n -> attachment -> Cdt/ds;
         drap  += n -> attachment -> Cdn/ds;
         drat_o  += n -> attachment -> Cdt/ds;
         drap_o  += n -> attachment -> Cdn/ds;
         cam   += n -> attachment -> m / ds;
         camma_t += (n -> attachment -> m + n -> attachment -> am)/ds;
         camma_n += (n -> attachment -> m + n -> attachment -> am)/ds;
         caarma_n += n -> attachment -> am/ds;
      }
      if (nm -> attachment) {
         w0_m    += nm -> attachment -> wet/ds;
         drat_m  += nm -> attachment -> Cdt/ds;
         drap_m  += nm -> attachment -> Cdn/ds;
         drat_om  += nm -> attachment -> Cdt/ds;
         drap_om  += nm -> attachment -> Cdn/ds;
         cam_m   += nm -> attachment -> m / ds;
         camma_t_m += (nm -> attachment -> m + nm -> attachment -> am)/ds;
         camma_n_m += (nm -> attachment -> m + nm -> attachment -> am)/ds;
         caarma_n_m += nm -> attachment -> am/ds;
      }

      if (problem -> type == Surface || problem -> type == Deployment) {
         if (nm -> x >= environment -> surface && w0_m < 0.0)
            w0_m = w0_m*(1.0 + tanh(50.0*(environment -> surface - nm -> x)));

         if (n -> x >= environment -> surface && w0 < 0.0)
            w0 = w0*(1.0 + tanh(50.0*(environment -> surface - n -> x)));
      }

      if (problem -> type == Riser) {
         if (nm -> x >= environment -> surface)
            w0_m = nm -> material -> w;

         if (n -> x >= environment -> surface)
            w0 = n -> material -> w;
      }

      if ((n -> position == TopBoundary && problem -> terminal [2] -> anchor) ||
         (n -> position == BranchTerminal && 
	  n -> segment -> branch -> terminal -> anchor) || !environment -> depth) {

         Fb_m = Fb = Fb_om = Fb_o = 0.0;
         mud_b = mud_bm = mud_bo = mud_bom = 0.0;
         mu = 0.0;
      }
      else {
         bottom    = Bottom(n -> y, n -> z, tm);
         bottom_m  = Bottom(nm -> y, nm -> z, tm);
         bottom_o  = Bottom(n -> y_o, n -> z_o, tm);
         bottom_om = Bottom(nm -> y_o, nm -> z_o, tm);

         mu = environment -> bottom_friction;

         if (nm -> x < bottom_m && w0_m > 0.0) {

            if (nm -> active_number == 1) 
               Fb_m = w0_m;
            else 
               Fb_m = (bottom_m - nm -> x)*environment -> bottom_stiffness;
/*
            if (Fb_m > w0_m && problem -> type != Towing)
               Fb_m = w0_m;
*/
            mud_bm = nm -> material -> mud_b;
         }
         else {
            Fb_m = mud_bm = 0.;
         }
    
         if (n -> x < bottom && w0 > 0.0) {
            Fb = (bottom - n -> x)*environment -> bottom_stiffness;
            mud_b = n -> material -> mud_b;
/*
            if (Fb > w0 && problem -> type != Towing)
               Fb = w0;
*/
         }
         else {
            Fb = mud_b =  0.0;
         }

         if (nm -> x_o < bottom_om && w0_m > 0.0) {

            if (nm -> active_number == 1)
               Fb_om = w0_m;
            else
               Fb_om = (bottom_om - nm -> x_o)*environment -> bottom_stiffness;
/*
            if (Fb_om > w0_m && problem -> type != Towing)
               Fb_om = w0_m;
*/
            mud_bom = nm -> material -> mud_b;
         }
         else {
            Fb_om = mud_bom = 0.0;
         }

         if (n -> x_o < bottom_o && w0 > 0.0) {
            Fb_o = (bottom_o - n -> x_o)*environment -> bottom_stiffness;
            mud_bo = n -> material -> mud_b;
/*
            if (Fb_o > w0 && problem -> type != Towing)
               Fb_o = w0;
*/
         }
         else {
            Fb_o = mud_bo = 0.0;
         }
      }

      pay = pay_o = pay_om = pay_m = 0;

      if (nm == n -> segment -> first_active) {
          pay_m = cam_m*n -> segment -> bottom_pay_speed;
          pay_om = cam_m*n -> segment -> bottom_pay_speed_o;
      }
      else if (n == n -> segment -> last_active) {
          pay = cam*n -> segment -> top_pay_speed;
          pay_o = cam*n -> segment -> top_pay_speed_o;
      }


      s [1][1]      =   omak2_ds*(tddkm*(e - e_m) - (tdk + tdkm))
                      + akak2_ds*(tddkm*(e_o - e_om) - (tdok + tdokm)) 
                      - omak*0.5*drat_m*urkm*urkma/sqstrkm;
      s [1][ne + 1] =   omak2_ds*(tddk*(e - e_m) + (tdk + tdkm))
    		      + akak2_ds*(tddk*(e_o - e_om) + (tdok + tdokm))
                      - omak*0.5*drat*urk*urka/sqstrk;
      s [1][2]      = -omak2_ds*(phi - phi_m) - akak2_ds*(phi_o - phi_om);
      s [1][ne + 2] = s [1][2];
      s [1][3]      = -camma_t_m*(omam2 + amam2)/gdt
                + omak2_ds*(pay_m + pay) + akak2_ds*(pay_o + pay_om)
                - omak*(damp_t + mud_bm + 2.0*drat_m*urkm*sign(urkm)*sqstrkm);
      s [1][ne + 3] = -camma_t*(omam2 + amam2)/gdt
                - omak2_ds*(pay + pay_m) - akak2_ds*(pay_o + pay_om)
                - omak*(damp_t + mud_b + 2.0*drat*urk*sign(urk)*sqstrk);
      s [1][4]      = cam_m*(omam2*phi_d_m + amam2*phi_d_om)
                      + omak2_ds*pay_m*(phi - phi_m) + akak2_ds*pay_m*(phi_o - phi_om);
      s [1][ne + 4] = cam*(omam2*phi_d + amam2*phi_d_o)
                      + omak2_ds*pay*(phi - phi_m) + akak2_ds*pay*(phi_o - phi_om);
      s [1][5]      = omak2_ds*(Sn + Sn_m) + akak2_ds*(Sn_o + Sn_om)
                      - omak2_ds*(pay_m*v_m + pay*v) - akak2_ds*(pay_om*v_om + pay_o*v_o)
                      + omam2*(cam_m*v_m/gdt + caarma_t_m*(vlckm/gdt
                                           - ulckm*phi_d_m))
                      + amam2*(cam_m*v_om/gdt + caarma_t_m*(vlcokm/gdt
                                           - ulckm*phi_d_om))
                      + omak*((w0_m - Fb_m)*sphkm
              - 2.0*drat_m*sign(urkm)*urkm*(sphkm*uackm - cphkm*vackm)*sqstrkm);
      s [1][ne + 5] = -omak2_ds*(Sn + Sn_m) - akak2_ds*(Sn_o + Sn_om)
                      + omak2_ds*(pay_m*v_m + pay*v) + akak2_ds*(pay_om*v_om + pay_o*v_o)
                      + omam2*(cam*v/gdt + caarma_t*(vlck/gdt
                                           - ulck*phi_d))
                      + amam2*(cam*v_o/gdt + caarma_t*(vlcok/gdt
                                           - ulck*phi_d_o))
                      + omak*((w0 - Fb)*sphk
                      - 2.0*drat*sign(urk)*urk*sqstrk*(sphk*uack - cphk*vack));
      s [1][rhs] =  
       omak2_ds*(  (tdk + tdkm)*(e - e_m) - (Sn + Sn_m)*(phi - phi_m)
                 - (pay_m + pay)*(u - u_m) + (pay_m*v_m + pay*v)*(phi - phi_m))
     + akak2_ds*(  (tdok + tdokm)*(e - e_m) - (Sn_o + Sn_om)*(phi - phi_m)
                 + (tdk + tdkm)*(e_o - e_om) - (Sn + Sn_m)*(phi_o - phi_om)
                 - (pay_om + pay_o)*(u - u_m) + (pay_om*v_om + pay_o*v_o)*(phi - phi_m)
                 - (pay_m + pay)*(u_o - u_om) + (pay_m*v_m + pay*v)*(phi_o - phi_om))
     + ak2_ds*((tdok + tdokm)*(e_o - e_om) - (Sn_o + Sn_om)*(phi_o - phi_om)
               - (pay_om + pay_o)*(u_o - u_om) + (pay_om*v_om + pay_o*v_o)*(phi_o - phi_om))
     + omam2*(-camma_t_m*u_d_m - camma_t*u_d 
	      + (cam_m*v_m + caarma_t_m*vlckm)*phi_d_m 
              + (cam*v + caarma_t*vlck)*phi_d)
     + amam2*(-camma_t_m*u_d_m - camma_t*u_d 
              + (cam_m*v_om + caarma_t_m*vlcokm)*phi_d_m 
              + (cam*v_o + caarma_t*vlcok)*phi_d
              -camma_t_m*u_d_om - camma_t*u_d_o 
              + (cam_m*v_m + caarma_t_m*vlckm)*phi_d_om 
              + (cam*v + caarma_t*vlck)*phi_d_o)
     + am2*(-camma_t_m*u_d_om - camma_t*u_d_o 
            + (cam_m*v_om + caarma_t_m*vlcokm)*phi_d_om 
            + (cam*v_o + caarma_t*vlcok)*phi_d_o)
     - omak*((w0 - Fb)*cphk + (w0_m - Fb_m)*cphkm 
             + mu*sign(urk)*Fb + mu*sign(urkm)*Fb_m
             + (damp_t + mud_b)*u + (damp_t + mud_bm)*u_m
             + drat*urk*urka*sqstrk + drat_m*urkm*urkma*sqstrkm
             - caarma_t*uldk - caarma_t_m*uldkm)
     - ak*((w0 - Fb_o)*cphok + (w0_m - Fb_om)*cphokm 
           + mu*sign(urok)*Fb_o + mu*sign(urokm)*Fb_om
           + (damp_t + mud_bo)*u_o + (damp_t + mud_bom)*u_om
           + drat_o*urok*uroka*sqstrok + drat_om*urokm*urokma*sqstrokm
           - caarma_t*uldok - caarma_t*uldokm);


      s [2][1]      = tdkm*(omak2_ds*(phi - phi_m) + akak2_ds*(phi_o - phi_om))
                      - omak*0.5*drap_m*vrkm*vrkma/sqstrkm;
      s [2][ne + 1] = tdk*(omak2_ds*(phi - phi_m) + akak2_ds*(phi_o - phi_om))
                      - omak*0.5*drap*vrk*vrka/sqstrk;
      s [2][2]      = -2.0*(omak2_ds + akak2_ds);
      s [2][ne + 2] = -s [2][2];
      s [2][3]      = -cam_m*(omam2*phi_d_m + amam2*phi_d_om)
                       - omak2_ds*pay_m*(phi - phi_m) - akak2_ds*pay_m*(phi_o - phi_om);
      s [2][ne + 3] = -cam*(omam2*phi_d + amam2*phi_d_o)
                       - omak2_ds*pay*(phi - phi_m) - akak2_ds*pay*(phi_o - phi_om);
      s [2][4]      = -camma_n_m*(omam2 + amam2)/gdt 
                 + omak2_ds*(pay_m + pay) + akak2_ds*(pay_o + pay_om)
                 - omak*(2.0*sign(vrkm)*vrkm*drap_m*sqstrkm + damp_n + mud_bm);
      s [2][ne + 4] = -camma_n*(omam2 + amam2)/gdt
                      - omak2_ds*(pay + pay_m) - akak2_ds*(pay_o + pay_om)
                      - omak*(2.0*sign(vrk)*vrk*drap*sqstrk + damp_n + mud_b);
      s [2][5]      = -omak2_ds*(tk + tkm) - akak2_ds*(tok + tokm)
                      + omak2_ds*(pay_m*u_m + pay*u) + akak2_ds*(pay_om*u_om + pay_o*u_o)
     - omam2*(cam_m*u_m/gdt + caarma_n_m*(ulckm/gdt 
                                    + vlckm*phi_d_m))
     - amam2*(cam_m*u_om/gdt + caarma_n_m*(ulcokm/gdt 
                                     + vlckm*phi_d_om))
     + omak*((w0_m - Fb_m)*cphkm
             - 2.0*drap_m*sign(vrkm)*vrkm*sqstrkm*(cphkm*uackm + sphkm*vackm));
      s [2][ne + 5] = omak2_ds*(tk + tkm) + akak2_ds*(tok + tokm)
                      - omak2_ds*(pay_m*u_m + pay*u) - akak2_ds*(pay_om*u_om + pay_o*u_o)
        - omam2*(cam*u/gdt 
	         + caarma_n*(ulck/gdt + vlck*phi_d))
        - amam2*(cam*u_o/gdt 
                 + caarma_n*(ulcok/gdt + vlck*phi_d_o))
        + omak*((w0 - Fb)*cphk 
                - 2.0*drap*sign(vrk)*vrk*sqstrk*(cphk*uack + sphk*vack));
      s [2][rhs] =  
       omak2_ds*(2.0*(Sn - Sn_m) + (tk + tkm)*(phi - phi_m)
                 - (pay_m + pay)*(v - v_m) - (pay_m*u_m + pay*u)*(phi - phi_m))
     + akak2_ds*(  2.0*(Sn - Sn_m) + (tok + tokm)*(phi - phi_m)
               + 2.0*(Sn_o - Sn_om) + (tk + tkm)*(phi_o - phi_om)
               - (pay_om + pay_o)*(v - v_m) - (pay_om*u_om + pay_o*u_o)*(phi - phi_m)
               - (pay_m + pay)*(v_o - v_om) - (pay_m*u_m + pay*u)*(phi_o - phi_om))
     + ak2_ds*(2.0*(Sn_o - Sn_om) + (tok + tokm)*(phi_o - phi_om)
               - (pay_om + pay_o)*(v_o - v_om) - (pay_om*u_om + pay_o*u_o)*(phi_o - phi_om))
     - omam2*(camma_n*v_d + camma_n_m*v_d_m 
              + (cam*u + caarma_n*ulck)*phi_d 
              + (cam_m*u_m + caarma_n_m*ulckm)*phi_d_m)
     - amam2*(camma_n*v_d + camma_n_m*v_d_m
              + (cam*u_o + caarma_n*ulcok)*phi_d
              + (cam_m*u_om + caarma_n_m*ulcokm)*phi_d_m
              + camma_n*v_d_o + camma_n_m*v_d_om
              + (cam*u + caarma_n*ulck)*phi_d_o
              + (cam_m*u_m + caarma_n_m*ulckm)*phi_d_om)
     - am2*(camma_n*v_d_o + camma_n_m*v_d_om
            + (cam*u_o + caarma_n*ulcok)*phi_d_o
            + (cam_m*u_om + caarma_n_m*ulcokm)*phi_d_om)
     + omak*((w0 - Fb)*sphk + (w0_m - Fb_m)*sphkm
             - (damp_n + mud_b)*v - (damp_n + mud_bm)*v_m
             - drap*vrk*vrka*sqstrk - drap_m*vrkm*vrkma*sqstrkm
             + caarma_n*vldk + caarma_n_m*vldkm)
     + ak*((w0 - Fb_o)*sphok + (w0_m - Fb_om)*sphokm
           - (damp_n + mud_bo)*v_o - (damp_n + mud_bom)*v_om
           - drap_o*vrok*vroka*sqstrok - drap_om*vrokm*vrokma*sqstrokm
           + caarma_n*vldok + caarma_n_m*vldokm);

      if (nm == n -> segment -> first_active) {
          pay_m = n -> segment -> bottom_pay_speed;
          pay_om = n -> segment -> bottom_pay_speed_o;
      }
      else if (n == n -> segment -> last_active) {
          pay = nm -> segment -> top_pay_speed;
          pay_o = nm -> segment -> top_pay_speed_o;
      }

      s [3][1]      = -omam/gdt;
      s [3][ne + 1] = -omam/gdt;
      s [3][3]      = -2.0*(omak2_ds + akak2_ds);
      s [3][ne + 3] = -s [3][3];
      s [3][4]      = -omak2_ds*(phi - phi_m) - akak2_ds*(phi_o - phi_om);
      s [3][ne + 4] = s [3][4];
      s [3][5]      = omak2_ds*(v + v_m) + akak2_ds*(v_o + v_om);
      s [3][ne + 5] = -s [3][5];
      s [3][rhs] = 
        omak2_ds*(2.0*(u - u_m) +2.0*(pay - pay_m) - (v + v_m)*(phi - phi_m))
      + akak2_ds*(  2.0*(u - u_m) +2.0*(pay - pay_m)- (v_o + v_om)*(phi - phi_m)
                  + 2.0*(u_o - u_om) +2.0*(pay_o - pay_om) - (v + v_m)*(phi_o - phi_om))
      + ak2_ds*(2.0*(u_o - u_om) +2.0*(pay_o - pay_om) - (v_o + v_om)*(phi_o - phi_om))
      - omam*(e_d + e_d_m) - am*(e_d_o + e_d_om); 

      s [4][1]      = -omam2*phi_d_m - amam2*phi_d_om;
      s [4][ne + 1] = -omam2*phi_d - amam2*phi_d_o;
      s [4][3]      = omak2_ds*(phi - phi_m) + akak2_ds*(phi_o - phi_om);
      s [4][ne + 3] = s [4][3];
      s [4][4]      = -2.0*omak2_ds - 2.0*akak2_ds;
      s [4][ne + 4] = -s [4][4];
      s [4][5]      = -omak2_ds*(u + u_m + pay + pay_m) 
                      - akak2_ds*(u_o + u_om + pay_o + pay_om)
                      - (omam2*(1.0 + e_m) + amam2*(1.0 + e_om))/gdt;
      s [4][ne + 5] = omak2_ds*(u + u_m + pay + pay_m) 
                      + akak2_ds*(u_o + u_om + pay_o + pay_om)
                      - (omam2*(1.0 + e) + amam2*(1.0 + e_o))/gdt;
      s [4][rhs] = 
        omak2_ds*(2.0*(v - v_m) + (u + u_m + pay_m + pay)*(phi - phi_m))
      + akak2_ds*(  2.0*(v - v_m) + (u_o + u_om + pay_o + pay_om)*(phi - phi_m)
                  + 2.0*(v_o - v_om) + (u + u_m + pay + pay_o)*(phi_o - phi_om))
      + ak2_ds*(2.0*(v_o - v_om) + (u_o + u_om + pay_o + pay_om)*(phi_o - phi_om))
      - omam2*((1.0 + e)*phi_d + (1.0 + e_m)*phi_d_m)
      - amam2*(  (1.0 + e_o)*phi_d + (1.0 + e_om)*phi_d_m
               + (1.0 + e)*phi_d_o + (1.0 + e_m)*phi_d_om)
      - am2*((1.0 + e_o)*phi_d_o + (1.0 + e_om)*phi_d_om);

      s [5][5]      = -omak*t_ds;
      s [5][ne + 5] = omak*t_ds;
      s [5][6]      = -omak;
      s [5][ne + 6] = -omak;
      s [5][rhs] = 
        omak*t_ds*(phi - phi_m) + ak*t_ds*(phi_o - phi_om)
      - omak*(Om3 + Om3_m) - ak*(Om3_o + Om3_om);

      s [6][1]      = 3.0*omak*Sn_m*(1.0 + e_m)*(1.0 + e_m);
      s [6][ne + 1] = 3.0*omak*Sn*(1.0 + e)*(1.0 + e);
      s [6][2]      = omak*elfakm;
      s [6][ne + 2] = omak*elfak;
      s [6][6]      = -omak*EI*t_ds;
      s [6][ne + 6] = omak*EI*t_ds;
      s [6][rhs] = 
        omak*t_ds*EI*(Om3 - Om3_m) + ak*EI*t_ds*(Om3_o - Om3_om)
      + omak*(Sn*elfak + Sn_m*elfakm) + ak*(Sn_o*elfaok + Sn_om*elfaokm);
                   
   }
       
   return;
}

static void InitialVelocity (phi, u, v, U, V)
   double	phi;
   double      *u;
   double      *v;
   double      *U, *V;
{
   if (problem -> type != Towing && problem -> type != Drifter && problem -> type != Deployment) {
      *u = 0.0;
      *v = 0.0;
      *U = 0.0;
      *V = 0.0;
    
      return;
   }

   if (problem -> type == Deployment)
      *V = problem -> terminal [1] -> yspeed.value;
   else
      *V = problem -> terminal [2] -> yspeed.value;

   *u = (*V)*sin(phi);
   *v = (*V)*cos(phi);
 
   *U = 0.0;
 
   return;
}

int SolveDynamicProblem2D (node, num_nodes, active, num_active, out, 
                           output_map, output_nodes, num_output_nodes, 
			               sample_dt, snap_dt, seg_dt, buoy_dt, ext_dt, 
                           decimate, prog_name, prog_dt, restart_name, restart_t)
   Node		*node;
   int	   	 num_nodes;
   Node		*active;
   int		 num_active;
   ResFile	 out;
   int		*output_map;
   int		*output_nodes;
   int		 num_output_nodes;
   double	 sample_dt;
   double	 snap_dt;
   double    seg_dt;
   double	 buoy_dt;
   double    ext_dt;
   int		 decimate;
   char     *prog_name;
   double    prog_dt;
   char     *restart_name;
   double    restart_t;
{
/*
   static double	scalv [] = {0, 0.07, 1000.0, 1.0, 1.0, 0.1, 0.1};
*/
   static double	scalv [] = {0, 0.07, 1000.0, 1.0, 1.0, 1.4, 0.0034};
   int           *adapt_count;
   int		  max_adapt;
   int		  level;
   double	**s;
   double	  t; 
   double	  dt;
   double	  u, v, U, V;
   int		  i, j;
   int	  	  it;
   char		  buffer [256];
   int		  displ_msg;
   int		  njn;
   double	  dt1g, gdt_inv;
   int		  dynstat_conv;
   Node		  close;
   Segment   *seg;
   int        nseg;
   ResFile    prog;
   double     xthrust, ythrust, zthrust;

   problem -> twoD = 1;
   problem -> dynamic = 1;

   ak    = analysis -> alpha_k;
   omak  = (1.0 - ak);
   am    = analysis -> alpha_m;
   omam  = (1.0 - am);
   omam2 = omam*omam;
   amam2 = am*omam;
   am2   = am*am;
   omg_g = (1.0 - analysis -> gamma) / analysis -> gamma;

   // seg = BuildSegmentArray(problem, &nseg, 1);
   seg = problem -> segment;
   nseg = problem -> num_segments;

   njn = problem -> junction_size;

   s = (double **) malloc (sizeof(double *) * (NE + NJ_COMPAT*njn)); s--;

   for (i = 1 ; i <= NE + NJ_COMPAT*njn ; i++) {
      s [i] = (double *) malloc(sizeof(double) * ((2 + njn)*NE + 1));     
      s [i] --;
   }
   
   environment -> surface = environment -> depth; // active[num_active] -> x;
//   fprintf(stderr,"%g %g %g %g\n", environment -> depth, environment -> surface,
//              active[1] -> x, active[num_active] -> x);

	/*
	 * copy the static solution 
	 */
   
   for (i = 1 ; i <= num_active ; i++) 
      for (j = 1 ; j <= 4 ; j++)
         active[i] -> Ys[j] = active[i] -> Y[j];

	// make sure we have sensible numbers for Ys in all nodes as they
	// get used in the output routines

   for (i = 1 ; i <= num_nodes ; i++) {
      if (!node[i] -> active) {
         if (node[i] -> number < node[i] -> segment -> first_active -> number) 
            close = node[i] -> segment -> first_active;
         else
            close = node[i] -> segment -> last_active;
      
         for (j = 1 ; j <= NE ; j++)
            node[i] -> Ys[j] = close -> Ys[j];
      }
   }

   for (i = 1 ; i <= num_active ; i++) {
      active[i] -> Y[1] = active[i] -> Y_o[1] 
                        = active[i] -> Y_f[1]
                        = active[i] -> Y_o_f[1] = active[i] -> Ys[1]; // e
      active[i] -> Y[2] = active[i] -> Y_o[2] 
                        = active[i] -> Y_f[2]
                        = active[i] -> Y_o_f[2] = active[i] -> Ys[2]; // Sn
      active[i] -> Y[5] = active[i] -> Y_o[5] 
                        = active[i] -> Y_f[5]
                        = active[i] -> Y_o_f[5] = active[i] -> Ys[3]; // phi
      active[i] -> Y[6] = active[i] -> Y_o[6] 
                        = active[i] -> Y_f[6]
                        = active[i] -> Y_o_f[6] = active[i] -> Ys[4]; // Om3

      InitialVelocity (active[i] -> Ys[3], &u, &v, &U, &V); /* handle towing problems */

      active[i] -> Y[3] = active[i] -> Y_o[3] = u;
      active[i] -> Y[4] = active[i] -> Y_o[4] = v;

      active[i] -> Y_f[3] = active[i] -> Y_o_f[3] = u;
      active[i] -> Y_f[4] = active[i] -> Y_o_f[4] = v;

      active[i] -> x_o = active[i] -> x_f = active[i] -> x_o_f = active [i] -> x;
      active[i] -> y_o = active[i] -> y_f = active[i] -> y_o_f = active [i] -> y;

      for (j = 1 ; j <= NE ; j++) {
         active[i] -> Yd_o[j] = 0.0;
      }

      active[i] -> xdot = active[i] -> xdot_f = active[i] -> xdot_o = active[i] -> xdot_o_f = U;
      active[i] -> ydot = active[i] -> ydot_f = active[i] -> ydot_o = active[i] -> ydot_o_f = V;
   }

   Thrust(restart_t, problem -> terminal [1], &xthrust, &ythrust, &zthrust);
   problem -> terminal[1] -> xthrust.value = xthrust;
   problem -> terminal[1] -> ythrust.value = ythrust;
   problem -> terminal[1] -> zthrust.value = zthrust;
   Thrust(restart_t, problem -> terminal [2], &xthrust, &ythrust, &zthrust);
   problem -> terminal[2] -> xthrust.value = xthrust;
   problem -> terminal[2] -> ythrust.value = ythrust;
   problem -> terminal[2] -> zthrust.value = zthrust;


   if (problem -> type == Surface || problem -> type == Drifter || problem -> type == Deployment)
      problem -> terminal [2] -> buoy -> draft = environment -> depth - active[num_active] -> x;

   dt = analysis -> dt;

   if (sample_dt == 0.0)
      sample_dt = dt;

   WriteDynamicHeader (out, 0.0, analysis -> duration, dt, sample_dt, snap_dt, 
                       seg_dt, buoy_dt, ext_dt,
                       num_output_nodes, output_nodes, decimate, node);

   if (prog && prog_dt)
      InitializeProgressFile(problem, prog, NE);

   if (snap_dt > 0.0)  {
      WriteDynamicSnapshot (out, output_map, 
                            node, num_nodes, decimate, 1);
   }
   if (seg_dt > 0.0) {
      WriteDynamicSegmentData(problem, out);
   }
   if (ext_dt > 0.0)
      WriteDynamicExternalForces(problem, out);

   WriteDynamicResult (out, output_map, output_nodes, 
                       num_output_nodes, decimate, 
		               node, num_nodes, 1);   

   if (analysis -> adapt_factor > 1) {
      max_adapt = analysis -> adapt_levels;	
      adapt_count = (int *) malloc(sizeof(int) * max_adapt);
      adapt_count --;

			/* this is the maximum number of times we will */
			/* try to reduce dt to get around problems     */
			/* before bailing out and giving up            */
   
      for (i = 1 ; i <= max_adapt ; i++)
         adapt_count [i] = 0;
   }
   else {
      max_adapt = 0;
      adapt_count = NULL;
   }

   level = 0;

	/*	
	 * the result at t = 0.0 is simply the static solution
	 * and we have already written it out so start at t = dt
	 */

   DisplayDynamicHeader ( );
 
   dynstat_conv = 0;

#if (defined HAVELAMP && !defined API)
   if (environment -> forcing == LAMP && !problem -> dynstat) 
      UpdateLAMP(node, Y, num_active, 0.0, out, buoy_dt, 1);
#endif

   for (t = dt ; t <= analysis -> duration + dt/2.0 ; t += dt) {
 
#if (defined HAVELAMP && !defined API)
      if (environment -> forcing == LAMP && !problem -> dynstat) 
         UpdateLAMP(node, Y, num_active, t, out, buoy_dt, 1);
#endif

      if (!level)
         dt = analysis -> dt;

      if (!problem -> dynstat) {
         num_active = ProcessSpools(problem, node, num_nodes, active, num_active, t, dt, NE, 1);
         if (num_active < 2) {
            DisplayMessage("not enough active nodes - all nodes spooled");
            SetError(C_NOMORENODES);
            return 1;
         }
      }

      it = SolveDE (DynamicDifeq2D, DynamicUpdate2D, 
                    &(analysis -> dynamic_it), &(analysis -> dynamic_tolerance),
                    &(analysis -> dynamic_relaxation), scalv, 
		            NE, NB, NB_BRANCH, NJ_COMPAT,
                    node, num_nodes, active, num_active, 
		            s, t, dt, 0.0, 0);

      if (problem -> solution -> userQuit) {
          return 1;
      }

      if (it && max_adapt) {
         if (level == max_adapt) {
            DisplayMessage("max adaptation level exceeded");
            SetError(C_MAXADAPTEXCEEDED);
            return 1;
         }

         level ++;
         adapt_count [level] = 1;
         t -= dt;
         dt /= (double) analysis -> adapt_factor;
      
         for (j = 1 ; j <= nseg ; j++) {
            seg[j] -> top_pay_speed = seg[j] -> top_pay_speed_o;
            seg[j] -> bottom_pay_speed = seg[j] -> bottom_pay_speed_o;
         }
         for (j = 1 ; j <= num_active ; j++) {
            active[j] -> x = active[j] -> x_o;
            active[j] -> y = active[j] -> y_o;

            active[j] -> x_f = active[j] -> x_o_f;
            active[j] -> y_f = active[j] -> y_o_f;

            for (i = 1 ; i <= NE ; i ++)  {
               active[j] -> Y[i] = active[j] -> Y_o[i];
               active[j] -> Y_f[i] = active[j] -> Y_o_f[i];
            }
         }

         DisplayMessage("adapting, dt = %g", dt);

         continue;
      }
      else if (it) {
         DisplayMessage("not adapting - count not proceed");
         return 1;
      }


#ifdef BUILD_TOWDEPTH_CONTROL
      if (problem -> type == Towing && (problem -> terminal [1] -> profile.expr || 
          problem -> terminal [1] -> profile.value || problem -> terminal [1] -> profile_m))
         TowDepthControl (t, dt, environment -> surface - active[1] -> x,
                          environment -> surface - active[1] -> x_o, 
                          active[1] -> y, active[1] -> y_o, active [num_active] -> s);
      else if (problem -> type == Towing && problem -> terminal [1] -> flap_file)
         FlapForces(t,  dt, environment -> surface - active[1] -> x,
                    environment -> surface - active[1] -> x_o, 
                    active[1] -> y, active[1] -> y_o);
#endif

      sprintf (buffer, "t = %g,", t);
      displ_msg = 0;
      if (num_output_nodes &&  check(t, dt, sample_dt) <= analysis -> dt/2000) {
         WriteDynamicResult (out, output_map, output_nodes, 
                             num_output_nodes, decimate, 
                             node, num_nodes, 1); 
         strcat (buffer, " result");
         displ_msg = 1;
#if (defined GUI || defined WINGUI)
         if (problem -> solution -> plotProgress) {
            ControlPlotTime(t, problem, environment);
         } 
#endif
      }

      if (snap_dt && check(t, dt, snap_dt) <= analysis -> dt/2000) {
         WriteDynamicSnapshot (out, output_map, 
                               node, num_nodes, decimate, 1);
         strcat (buffer, " snapshot");
         displ_msg = 1;
#if (defined GUI || defined WINGUI)
         if (problem -> solution -> plotProgress) {
            ControlPlotSnaps(problem, environment);
         } 
#endif
      }

      if (seg_dt && check(t, dt, seg_dt) <= analysis -> dt/2000) {
          WriteDynamicSegmentData(problem, out);
          strcat(buffer, " seg");
          displ_msg = 1;
      }
      if (ext_dt && check(t, dt, ext_dt) <= analysis -> dt/2000) {
          WriteDynamicExternalForces(problem, out);
          strcat(buffer, " ext");
          displ_msg = 1;
      }
      if (prog && prog_dt && check(t, dt, prog_dt) <= analysis -> dt/2000) {
          WriteProgress(t, problem, prog, NE);
          strcat(buffer, " prog");
          displ_msg = 1;
      }


      
      if (displ_msg)
         DisplayMessage (buffer);

      if (problem -> dynstat) {
         dynstat_conv = CheckDynstatConvergence(active, num_active, dt, scalv, NE, 1);
         if (dynstat_conv)
            break;
      }

      dt1g = dt*(1.0 - analysis -> gamma);
      gdt_inv = 1.0/analysis -> gamma/dt;
      for (i = 1 ; i <= nseg ; i++) {
         seg[i] -> top_pay_speed_o = seg[i] -> top_pay_speed;
         seg[i] -> bottom_pay_speed_o = seg[i] -> bottom_pay_speed;
      }

      for (i = 1 ; i <= num_active ; i++) {
         for (j = 1 ; j <= NE ; j++) {
            active[i] -> Yd_o[j] = gdt_inv*(active[i] -> Y [j] - active[i] -> Y_o [j] - dt1g*active[i] -> Yd_o [j]);
            active[i] -> Y_o[j] = active[i] -> Y [j];
            active[i] -> Y_o_f[j] = active[i] -> Y_f [j];
            
         }

         active[i] -> xdot_o =  gdt_inv*(active[i] -> x - active[i] -> x_o - dt1g*active[i] -> xdot_o);
         active[i] -> ydot_o =  gdt_inv*(active[i] -> y - active[i] -> y_o - dt1g*active[i] -> ydot_o);
        
         active[i] -> xdot_o_f = active[i] -> xdot_f;
         active[i] -> ydot_o_f = active[i] -> ydot_f;

         active[i] -> x_o = active[i] -> x;
         active[i] -> y_o = active[i] -> y;

         active[i] -> x_o_f = active[i] -> x_f; 
         active[i] -> y_o_f = active[i] -> y_f; 

         active[i] -> pay_o = active[i] -> pay;

         active[i] -> pay_o_f = active[i] -> pay_f;
      }

      SmoothNodeData(t, active, num_active, 1);

      if (level) {
         if (adapt_count [level] == analysis -> adapt_factor) {
            dt *= analysis -> adapt_factor;
            level --;

        /*
         * this is to allow for the possibility that a previous level
         * or a succession of previous levels was on the last step
         * before we decided we needed to adapt down again ... if any
         * level was not on the last step then we return to that level
         */

            for (j = level ; j >= 1 ; j--) {
               if (adapt_count [j] >= analysis -> adapt_factor) {
                  dt *= analysis -> adapt_factor;
                  level --;
               }
               else
                  break;
            }

            DisplayMessage("adapting back, dt = %g", dt);

            if (level)
               adapt_count [level] ++;
         }
         else
            adapt_count [level] ++;
      }
   }

	/*
	 * copy the latest dynamic solution into the static
	 * solution in case the driver routines needs the
	 * latest info for any reason (dynstat option for instance)
	 */

   for (i = 1 ; i <= num_nodes ; i++) {
      node[i] -> Ys[1] = node[i] -> Y[1];
      node[i] -> Ys[2] = node[i] -> Y[2];
      node[i] -> Ys[3] = node[i] -> Y[3] = node[i] -> Y[5];
      node[i] -> Ys[4] = node[i] -> Y[4] = node[i] -> Y[6];
   }


   for (i = 1 ; i <= NE + njn*NJ_COMPAT ; i++) {
      s [i] ++; free (s [i]);
   }

   s ++; free (s);

   return 0;
}
