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
  * File:	dynamic.c
  *
  * Description: contains the main solver for the cable model
  *
  * History:
  *		
  ****************************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <malloc.h>
# include <math.h>
# include <string.h>
# include "compress.h"
# include "problem.h"
# include "tension.h"
# include "output.h"
# include "error.h"
# include "solve.h"
# include "transforms.h"
# include "segments.h"

# define SQR(a) ((a)*(a))

# define NE	   13
# define NJ_COMPAT 3
# define NB	   7
# define NB_BRANCH 4

# define TOLERANCE 1e-20

extern Problem     *problem;             
extern Environment *environment;
extern Analysis    *analysis;

# if (defined HAVELAMP && !defined API)
extern int UpdateLAMP(Node *, double **, int, double, ResFile, double, int);
# endif

static   double        ak, omak;
static   double        am, omam;
static   double	       omam2, amam2, am2, omg_g, omak_2;

static void 
GetXYZDot(Node a, double *xdot, double *ydot, double *zdot)
{
   double	u, v, w;
   double	u_o, v_o, w_o;
   double	B0, B1, B2, B3;
   double	B0_o, B1_o, B2_o, B3_o;
   double	xdot_o, ydot_o, zdot_o;

   u    = a -> Y [4];
   u_o  = a -> Y_o [4];
   v    = a -> Y [5];
   v_o  = a -> Y_o [5];
   w    = a -> Y [6];
   w_o  = a -> Y_o [6];
   B0   = a -> Y [7];
   B0_o = a -> Y_o [7];
   B1   = a -> Y [8];
   B1_o = a -> Y_o [8];
   B2   = a -> Y [9];
   B2_o = a -> Y_o [9];
   B3   = a -> Y [10];
   B3_o = a -> Y_o [10];

   xdot_o = a -> xdot_o;
   ydot_o = a -> ydot_o;
   zdot_o = a -> zdot_o;

   *xdot = (omak*XComponent(u,v,w,B0,B1,B2,B3)
                + ak*XComponent(u_o,v_o,w_o,B0_o,B1_o,B2_o,B3_o)
                - am*xdot_o) / omam;
   *ydot = (omak*YComponent(u,v,w,B0,B1,B2,B3)
                + ak*YComponent(u_o,v_o,w_o,B0_o,B1_o,B2_o,B3_o)
                - am*ydot_o) / omam;
   *zdot = (omak*ZComponent(u,v,w,B0,B1,B2,B3)
                + ak*ZComponent(u_o,v_o,w_o,B0_o,B1_o,B2_o,B3_o)
                - am*zdot_o) / omam;

   return;
}

static void 
IntegrateSpatialXYZ(Node start)
{
   double	B0, B1, B2, B3;
   double	B0_o, B1_o, B2_o, B3_o;
   double	B0_m, B1_m, B2_m, B3_m;
   double	B0_om, B1_om, B2_om, B3_om;
   double	sf, sf_m, sf_om, sf_o;
   Node		a, a_m;
   double	ds;
   double   xdot, ydot, zdot;

   a = start;
   while (a && a -> active) {
       a_m = a -> prev_active;

       ds = 0.5*a_m -> ds;

       B0 = a -> Y [7];
       B1 = a -> Y [8];
       B2 = a -> Y [9];
       B3 = a -> Y [10];
       B0_o = a -> Y_o [7];
       B1_o = a -> Y_o [8];
       B2_o = a -> Y_o [9];
       B3_o = a -> Y_o [10];
       B0_m = a_m -> Y [7];
       B1_m = a_m -> Y [8];
       B2_m = a_m -> Y [9];
       B3_m = a_m -> Y [10];
       B0_om = a_m -> Y_o [7];
       B1_om = a_m -> Y_o [8];
       B2_om = a_m -> Y_o [9];
       B3_om = a_m -> Y_o [10];

       sf    = (1.0 + a -> Y [1])*ds;
       sf_o  = (1.0 + a_m -> Y_o [1])*ds;
       sf_m  = (1.0 + a_m -> Y [1])*ds;
       sf_om = (1.0 + a_m -> Y_o [1])*ds;

       a -> x = a_m -> x  
                       + (omak*(XtComponent(sf,B0,B1,B2,B3) 
                                + XtComponent(sf_m,B0_m,B1_m,B2_m,B3_m))
                          + ak*(XtComponent(sf_o,B0_o,B1_o,B2_o,B3_o)
                                + XtComponent(sf_om,B0_om,B1_om,B2_om,B3_om))
                          - ak*(a -> x_o - a_m -> x_o)) / omak;

       a -> y = a_m -> y  
                       + (omak*(YtComponent(sf,B0,B1,B2,B3)
                                + YtComponent(sf_m,B0_m,B1_m,B2_m,B3_m))
                          + ak*(YtComponent(sf_o,B0_o,B1_o,B2_o,B3_o)
                                + YtComponent(sf_om,B0_om,B1_om,B2_om,B3_om))
                          - ak*(a -> y_o - a_m -> y_o)) / omak;

       a -> z = a_m -> z  
 		       + (omak*(ZtComponent(sf,B0,B1,B2,B3)
    			        + ZtComponent(sf_m,B0_m,B1_m,B2_m,B3_m))
                          + ak*(ZtComponent(sf_o,B0_o,B1_o,B2_o,B3_o)
                                + ZtComponent(sf_om,B0_om,B1_om,B2_om,B3_om))
                          - ak*(a -> z_o - a_m -> z_o)) / omak;

    
      GetXYZDot(a, &xdot, &ydot, &zdot);
      a -> xdot = xdot;
      a -> ydot = ydot;
      a -> zdot = zdot;

      a = a -> next_active;
   }

   return;
}

static void 
IntegrateTemporalXYZ(Node start, double dt, double xdot, double ydot, double zdot)
{
   double	B0, B1, B2, B3;
   double	B0_o, B1_o, B2_o, B3_o;
   double	B0_m, B1_m, B2_m, B3_m;
   double	B0_om, B1_om, B2_om, B3_om;
   double	xdot_m, ydot_m, zdot_m;
   double	xdot_om, ydot_om, zdot_om;
   double	xdot_o, ydot_o, zdot_o;
   double	u_o, v_o, w_o;
   double	u, v, w;
   double	u_m, v_m, w_m;
   double	u_om, v_om, w_om;
   Node		a, a_m;
   double	gamma, gdt, dt1g;
   double   ds, ds0;

   gamma = analysis -> gamma;
   gdt = gamma*dt;
   dt1g = dt*(1.0 - gamma);

   a = start;
   while (a && a -> active) {
       a_m = a -> prev_active;


       ds0 += a_m -> ds;
       ds += a_m -> ds*(1.0 + 0.5*(a -> Y[1] + a_m -> Y[1]));

       xdot_m = xdot;
       ydot_m = ydot;
       zdot_m = zdot;

       B0 = a -> Y [7];
       B1 = a -> Y [8];
       B2 = a -> Y [9];
       B3 = a -> Y [10];
       B0_o = a -> Y_o [7];
       B1_o = a -> Y_o [8];
       B2_o = a -> Y_o [9];
       B3_o = a -> Y_o [10];
       B0_m = a_m -> Y [7];
       B1_m = a_m -> Y [8];
       B2_m = a_m -> Y [9];
       B3_m = a_m -> Y [10];
       B0_om = a_m -> Y_o [7];
       B1_om = a_m -> Y_o [8];
       B2_om = a_m -> Y_o [9];
       B3_om = a_m -> Y_o [10];

       u = a -> Y [4];
       v = a -> Y [5];
       w = a -> Y [6];
       u_m = a_m -> Y [4];
       v_m = a_m -> Y [5];
       w_m = a_m -> Y [6];
       u_o = a -> Y_o [4];
       v_o = a -> Y_o [5];
       w_o = a -> Y_o [6];
       u_om = a_m -> Y_o [4];
       v_om = a_m -> Y_o [5];
       w_om = a_m -> Y_o [6];

       xdot_om = a_m -> xdot_o;
       ydot_om = a_m -> ydot_o;
       zdot_om = a_m -> zdot_o;
       xdot_o = a -> xdot_o;
       ydot_o = a -> ydot_o;
       zdot_o = a -> zdot_o;
    
       xdot = (omak*(XComponent(u,v,w,B0,B1,B2,B3)
                        + XComponent(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m))
                  + ak*(XComponent(u_o,v_o,w_o,B0_o,B1_o,B2_o,B3_o)
                        + XComponent(u_om,v_om,w_om,B0_om,B1_om,B2_om,B3_om))
                  - am*(xdot_o + xdot_om)) / omam - xdot_m;
       ydot = (omak*(YComponent(u,v,w,B0,B1,B2,B3)
                        + YComponent(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m))
                  + ak*(YComponent(u_o,v_o,w_o,B0_o,B1_o,B2_o,B3_o)
                        + YComponent(u_om,v_om,w_om,B0_om,B1_om,B2_om,B3_om))
                  - am*(ydot_o + ydot_om)) / omam - ydot_m;
       zdot = (omak*(ZComponent(u,v,w,B0,B1,B2,B3)
                        + ZComponent(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m))
                  + ak*(ZComponent(u_o,v_o,w_o,B0_o,B1_o,B2_o,B3_o)
                        + ZComponent(u_om,v_om,w_om,B0_om,B1_om,B2_om,B3_om))
                  - am*(zdot_o + zdot_om)) / omam - zdot_m;

       a -> x = a -> x_o + dt1g*xdot_o + gdt*xdot;
       a -> y = a -> y_o + dt1g*ydot_o + gdt*ydot;
       a -> z = a -> z_o + dt1g*zdot_o + gdt*zdot;

       a -> xdot = xdot;
       a -> ydot = ydot;
       a -> zdot = zdot;

       a = a -> next_active;
   }

   // fprintf(stderr,"total active lengths: unstretched = %g, stretched = %g\n", ds0, ds);

   return;
}

void 
DynamicUpdate3D (
   Node		 *active,
   int	      num_active,
   double	  tm,
   double	  dt)
{
   Node	  	  first, from;
   int		  i;
   double	  xdot, ydot, zdot;
   double         gamma;
   double	  gdt;
   double         dt1g;
   Node		  a;
   

   gamma = analysis -> gamma;
   gdt = gamma*dt;
   dt1g = dt*(1.0 - gamma);


   a = active[1];

   if (problem -> type == Deployment || 
       problem -> type == Towing || 
       problem -> type == HorizontalDrifter ||
       problem -> type == Drifter ||
       (problem -> terminal [1] -> release > 0.0 &&
           tm >= problem -> terminal [1] -> release && !problem -> dynstat)) {

       GetXYZDot(a, &xdot, &ydot, &zdot);

       a -> x = a -> x_o + dt1g*a -> xdot_o + gdt*xdot;
       a -> y = a -> y_o + dt1g*a -> ydot_o + gdt*ydot;
       a -> z = a -> z_o + dt1g*a -> zdot_o + gdt*zdot;
       a -> xdot = xdot;
       a -> ydot = ydot;
       a -> zdot = zdot;
   }
   else {
      a -> xdot = xdot = 0.0;
      a -> ydot = ydot = 0.0;
      a -> zdot = zdot = 0.0;

      a -> x = problem -> terminal [1] -> x;
      a -> y = problem -> terminal [1] -> y;
      a -> z = problem -> terminal [1] -> z;
   }

   if (analysis -> integration == Spatial) 
      IntegrateSpatialXYZ(a -> next_active);
   else 
      IntegrateTemporalXYZ(a -> next_active, dt, xdot, ydot, zdot);

   for (i = 1 ; i <= problem -> num_branch ; i++) {
      first = problem -> branch [i] -> segment[1] -> first_active;
      from = problem -> branch[i] -> segment_from -> last_active;

      first -> x = from -> x;
      first -> y = from -> y;
      first -> z = from -> z;

      GetXYZDot(from, &xdot, &ydot, &zdot);
      first -> xdot = xdot;
      first -> ydot = ydot;
      first -> zdot = zdot;

      if (analysis -> integration == Spatial) 
         IntegrateSpatialXYZ(first -> next_active);
      else { 
         IntegrateTemporalXYZ(first -> next_active, dt, xdot, ydot, zdot);
      }
   }
/*
   printf("seg 1: %f %f %f\n", tm, active[1] -> ds*(1.0 + active[1] -> Y[1]),
          sqrt(pow(active[1] -> x - active[2] -> x,  2.0) 
               + pow(active[1] -> y - active[2] -> y, 2.0)
               + pow(active[1] -> z - active[2] -> z, 2.0)));
   printf("seg n: %f %f %f\n", tm, active[num_active-1] -> ds*(1.0 + active[num_active-1] -> Y[1]),
          sqrt(pow(active[num_active-1] -> x - active[num_active] -> x,  2.0) 
               + pow(active[num_active-1] -> y - active[num_active] -> y, 2.0)
               + pow(active[num_active-1] -> z - active[num_active] -> z, 2.0)));
*/ 
   return;
}

void DynamicDifeq3D (
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
   double         current_factor)		/* not used		*/
{
   static double   ulcpre, vlcpre, wlcpre;
   static double   uacpre, vacpre, wacpre;
   static double   ulcopre, vlcopre, wlcopre;
   static double   uacopre, vacopre, wacopre;
   // static double   uldpre, vldpre, wldpre;
   // static double   uldopre, vldopre, wldopre;
   Node		   nj;
   Connector	   c;
   Buoy		   buoy;
   Anchor	   a;
   double	   e, Sn, Sb, B0, B1, B2, B3, Om1, Om2, Om3;
   double	   e_m, Sn_m, Sb_m; 
   double          B0_m, B1_m, B2_m, B3_m, Om1_m, Om2_m, Om3_m;
   double	   e_j, Sn_j, Sb_j; 
   double          B0_j, B1_j, B2_j, B3_j;
   double	   e_o, Sn_o, Sb_o; 
   double          B0_o, B1_o, B2_o, B3_o, Om1_o, Om2_o, Om3_o;
   double	   e_om, Sn_om, Sb_om; 
   double          B0_om, B1_om, B2_om, B3_om, Om1_om, Om2_om, Om3_om;
   // double	   e_oj, Sn_oj, Sb_oj; 
   double          B0_oj, B1_oj, B2_oj, B3_oj;
   double          u, v, w, u_o, v_o, w_o;
   double          u_m, v_m, w_m, u_om, v_om, w_om;
   double	   u_j, v_j, w_j, u_oj, v_oj, w_oj;
   double	   e_d, u_d, v_d, w_d; // Sn_d, Sb_d
   double	   B0_d, B1_d, B2_d, B3_d;
   double	   e_d_m, u_d_m, v_d_m, w_d_m; // Sn_d_m, Sb_d_m
   double	   B0_d_m, B1_d_m, B2_d_m, B3_d_m;
   double	   e_d_o, u_d_o, v_d_o, w_d_o; // Sn_d_o, Sb_d_o
   double	   B0_d_o, B1_d_o, B2_d_o, B3_d_o;
   double	   e_d_om, u_d_om, v_d_om, w_d_om; // Sn_d_om, Sb_d_om
   double  	   B0_d_om, B1_d_om, B2_d_om, B3_d_om;
   double          U, V, W, U_o, V_o, W_o;
   double          U_m, V_m, W_m, U_om, V_om, W_om;
   double          U_j, V_j, W_j, U_oj, V_oj, W_oj;
   double	   Ua, Ua_o, Va, Va_o, Wa, Wa_o;
   double	   Ur, Wr, Vr;
   double	   EI, GJ, w0, drat, drap, drat_o, drap_o;
   double          w0_m, drat_m, drap_m, drat_om, drap_om;
   double	   damp_t;
   double	   damp_n;
   double	   ds;
   double	   xforst, yforst, zforst;
   double          xthrust, ythrust, zthrust;
   int		   i, j;
   int		   njn;
   double	   sign_j;
   double	   ep, ep_m, ep_o, ep_om;
   double          tk, tkm, tkj, tdk, tdkm, tdkj, tddk, tddkm;
   double          tok, tokm, tdok, tdokm; // tokj
   double	   sqstrk, sqstrkm, sqstrok, sqstrokm;
   double	   elfak, elfakm, elfaok, elfaokm;
   double	   urk, vrk, wrk, urkm, vrkm, wrkm;
   double	   urok, vrok, wrok, urokm, vrokm, wrokm;
   double	   ulck, vlck, wlck, ulckm, wlckm, vlckm;
   // double          uldk, vldk, vldkm, wldk, wldkm; uldkm
   double 	   uack, vack, wack, uackm, vackm, wackm;
   double	   ulcok, vlcok, wlcok, ulcokm, wlcokm, vlcokm;
   // double          uldok, vldok, vldokm, wldok, wldokm; // uldokm
   double 	   uacok, vacok, wacok, uacokm, vacokm, wacokm;
   double	   pervk, pervkm;
   double	   pervok, pervokm;
   double	   ftr, ftr_m;
   double	   ftrb, ftrn, ftrb_m, ftrn_m;
   double	   urka, urkma, urokma, uroka;
   double          cam, camma_t, camma_n, caarma_n, caarma_t;
   double          cam_m, camma_n_m, camma_t_m,caarma_n_m, caarma_t_m;
   double	   t_ds;
   double	   wet;
   double	   Fex [4];
   double	   wave_b, B;
   double	   Uw, Vw, Ww, Uw_m, Vw_m, Ww_m;
   double	   Udw, Vdw, Wdw, Udw_m, Vdw_m, Wdw_m;
   double	   Uw_o, Vw_o, Ww_o, Uw_om, Vw_om, Ww_om;
   double	   Udw_o, Vdw_o, Wdw_o, Udw_om, Vdw_om, Wdw_om;
   double	   nav;
   double          barma, bmma, bdr;
   double	   Uspd, Vspd, Wspd;
   double      pay, pay_m, pay_o, pay_om;
   double	   Fy_wind, Fz_wind;
   double	   mma, dr_t, dr_n;
   double          carma;
   double          mud_b, mud_bo;
   double          mud_bm, mud_bom;
   double	   Fb, Fb_om, Fb_m, Fb_o;
   double	   bottom_m, bottom, bottom_o, bottom_om;
   double	   mu;
   double          gdt;
   double          omak2_ds, akak2_ds, ak2_ds;
   ForcingMethod   forcing;

   forcing = environment -> forcing; 

   ds = nm -> ds;
   
   w0 = n -> material -> wet;
   EI = n -> material -> EI;
   GJ = n -> material -> GJ;
   damp_t = n -> material -> bt;
   damp_n = n -> material -> bn;

   e  = n -> Y [1];
   Sn = n -> Y [2];
   Sb = n -> Y [3];
   u  = n -> Y [4];
   v  = n -> Y [5];
   w  = n -> Y [6];
   B0 = n -> Y [7];
   B1 = n -> Y [8];
   B2 = n -> Y [9];
   B3 = n -> Y [10];
   Om1 = n -> Y [11];
   Om2 = n -> Y [12];
   Om3 = n -> Y [13];

   e_o  = n -> Y_o [1];
   Sn_o = n -> Y_o [2];
   Sb_o = n -> Y_o [3];
   u_o  = n -> Y_o [4];
   v_o  = n -> Y_o [5];
   w_o  = n -> Y_o [6];
   B0_o = n -> Y_o [7];
   B1_o = n -> Y_o [8];
   B2_o = n -> Y_o [9];
   B3_o = n -> Y_o [10];
   Om1_o = n -> Y_o [11];
   Om2_o = n -> Y_o [12];
   Om3_o = n -> Y_o [13];

   tk     = Tension(e, n -> material);
   tok    = Tension(e_o, n -> material);
   tdk    = TensionD(e, n -> material);
   tdok   = TensionD(e_o, n -> material);
   tddk   = TensionDD(e, n -> material);

   ep = 1.0 + e;
   ep_o = 1.0 + e_o;

   sqstrk  = sqrt(ep);
   sqstrok  = sqrt(ep_o);
   elfak  = pow(ep, 3.0);
   elfaok  = pow(ep_o, 3.0);

   cam = n -> material -> m;
   camma_n = n -> material -> m + n -> material -> amn;
   camma_t = n -> material -> m + n -> material -> amt;
   caarma_n = n -> material -> amn + n -> material -> rV;
   if (n -> material -> amt)
      caarma_t = n -> material -> amt + n -> material -> rV;
   else 
      caarma_t = 0.0;

   w0_m = nm -> material -> wet;

   e_m = nm -> Y [1];
   Sn_m = nm -> Y [2];
   Sb_m = nm -> Y [3];
   u_m  = nm -> Y [4];
   v_m  = nm -> Y [5];
   w_m  = nm -> Y [6];
   B0_m = nm -> Y [7];
   B1_m = nm -> Y [8];
   B2_m = nm -> Y [9];
   B3_m = nm -> Y [10];
   Om1_m = nm -> Y [11];
   Om2_m = nm -> Y [12];
   Om3_m = nm -> Y [13];

   e_om  = nm -> Y_o [1];
   Sn_om = nm -> Y_o [2];
   Sb_om = nm -> Y_o [3];
   u_om  = nm -> Y_o [4];
   v_om  = nm -> Y_o [5];
   w_om  = nm -> Y_o [6];
   B0_om = nm -> Y_o [7];
   B1_om = nm -> Y_o [8];
   B2_om = nm -> Y_o [9];
   B3_om = nm -> Y_o [10];
   Om1_om = nm -> Y_o [11];
   Om2_om = nm -> Y_o [12];
   Om3_om = nm -> Y_o [13];

   tkm    = Tension(e_m, nm -> material);
   tokm   = Tension(e_om, nm -> material);
   tdkm   = TensionD(e_m, nm -> material);
   tdokm  = TensionD(e_om, nm -> material);
   tddkm  = TensionDD(e_m, nm -> material);

   cam_m    = nm -> material -> m;
   camma_n_m  = nm -> material -> m + nm -> material -> amn;
   camma_t_m  = nm -> material -> m + nm -> material -> amt;
   caarma_n_m = nm -> material -> amn + nm -> material -> rV;
   if (nm -> material -> amt)
      caarma_t_m = nm -> material -> amt + nm -> material -> rV;
   else 
      caarma_t_m = 0.0;

   ep_m = 1.0 + e_m;
   ep_om = 1.0 + e_om;

   sqstrkm = sqrt(ep_m);
   elfakm = pow(ep_m, 3.0);
   sqstrokm = sqrt(ep_om);
   elfaokm = pow(ep_om, 3.0);

   for (i = 1 ; i <= num_rows ; i++)
      for (j = 1 ; j <= rhs ; j++)
         s [i][j] = BLANK;


	/*
	 * anchor node
	 */

   if (eq_type == BottomBoundary) {
      s [1][7] = 1.0;
      s [1][rhs] = B0;

      for (i = 1 ; i <= 10 ; i++) {
         s [2][i] = 0.0;
         s [3][i] = 0.0;
         s [4][i] = 0.0;
      }

      if (problem -> terminal [1] -> release > 0.0 &&
          tm >= problem -> terminal [1] -> release && !problem -> dynstat) {


         s [2][1] = 1.0;
         s [2][rhs] = e;
    
         s [3][2] = 1.0;
         s [3][rhs] = Sn;
    
         s [4][3] = 1.0;
         s [4][rhs] = Sb;
     }
     else if ((problem -> dynstat && problem -> terminal [1] -> anchor) ||
	      (problem -> terminal [1] -> anchor && problem -> type != Deployment) ||
              (problem -> type == Deployment && n -> x <= Bottom(n -> y, n -> z, tm)) ||
              (problem -> type == Towing && environment -> depth
               && n -> x <= Bottom(n -> y, n -> z, tm))) {


         s [2][4] = 1.0;
         s [2][rhs] = u;
    
         s [3][5] = 1.0;
         s [3][rhs] = v;
    
         s [4][6] = 1.0;
         s [4][rhs] = w;
      }
      else {
         if (problem -> type == Deployment) {
            a = problem -> terminal [1] -> anchor;
            bmma = a -> m + a -> am;
            barma = 0.0;
            bdr = a -> Cdn;
            wet = a -> wet;

            Uw = Vw = Ww = 0.0;
            Udw = Vdw = Wdw = 0.0;

            xthrust = ythrust = zthrust = 0.0;
         }
         else {
            buoy = problem -> terminal [1] -> buoy;

            bmma = buoy -> Mmma;
            barma = buoy -> Marma;
            bdr = buoy -> Mdr;

            wet = buoy -> w - buoy -> buoyancy;

            if (forcing != Velocity && forcing != Force) {
               WaveParticleVelocity (tm, n -> x, n -> y, n -> z, &Uw, &Vw, &Ww);
               WaveParticleAcceleration (tm, n -> x, n -> y, n -> z,
                                         &Udw, &Vdw, &Wdw);
            }
            else {
               Uw = Vw = Ww = 0.0;
               Udw = Vdw = Wdw = 0.0;
            }

            if (problem -> type == Towing || problem -> type == HorizontalDrifter) {
               Thrust(tm, problem -> terminal [1], &xthrust, &ythrust, &zthrust);
               problem -> terminal[1] -> xthrust.value = xthrust;
               problem -> terminal[1] -> ythrust.value = ythrust;
               problem -> terminal[1] -> zthrust.value = zthrust;
            }
            else { 
               xthrust = 0.0;
               ythrust = 0.0;
               zthrust = 0.0;
            }
         }

         Current (tm, n -> x, n -> y, n -> z, &uack, &vack, &wack);

         U   = XComponent(u,v,w,B0,B1,B2,B3);
         U_o = XComponent(u_o,v_o,w_o,B0_o,B1_o,B2_o,B3_o);

         V   = YComponent(u,v,w,B0,B1,B2,B3);
         V_o = YComponent(u_o,v_o,w_o,B0_o,B1_o,B2_o,B3_o);

         W   = ZComponent(u,v,w,B0,B1,B2,B3);
         W_o = ZComponent(u_o,v_o,w_o,B0_o,B1_o,B2_o,B3_o);

         Ur   = U - Uw;
         Vr   = V - vack - Vw;
         Wr   = W - wack - Ww;

         if (problem -> type == HorizontalDrifter) {
            s [2][4]  = dXdt(B0,B1,B2,B3);
            s [2][5]  = dXdn(B0,B1,B2,B3);
            s [2][6]  = dXdb(B0,B1,B2,B3);
            s [2][7]  = dXdB0(u,v,w,B0,B1,B2,B3);
            s [2][8]  = dXdB1(u,v,w,B0,B1,B2,B3);
            s [2][9]  = dXdB2(u,v,w,B0,B1,B2,B3);
            s [2][10]  = dXdB3(u,v,w,B0,B1,B2,B3);
            s [2][rhs] = XComponent(u,v,w,B0,B1,B2,B3);
         }
         else {
             ftr = bmma/dt + 2.0*bdr*sign(Ur)*Ur;

             s [2][1] = -tdk*dXdt(B0,B1,B2,B3);
             s [2][2] = -dXdn(B0,B1,B2,B3);
             s [2][3] = -dXdb(B0,B1,B2,B3);
             s [2][4] = ftr*dXdt(B0,B1,B2,B3);
             s [2][5] = ftr*dXdn(B0,B1,B2,B3);
             s [2][6] = ftr*dXdb(B0,B1,B2,B3);
             s [2][7] = ftr*dXdB0(u,v,w,B0,B1,B2,B3) - dXdB0(tk,Sn,Sb,B0,B1,B2,B3);
             s [2][8] = ftr*dXdB1(u,v,w,B0,B1,B2,B3) - dXdB1(tk,Sn,Sb,B0,B1,B2,B3);
             s [2][9] = ftr*dXdB2(u,v,w,B0,B1,B2,B3) - dXdB2(tk,Sn,Sb,B0,B1,B2,B3);
             s [2][10] = ftr*dXdB3(u,v,w,B0,B1,B2,B3) - dXdB3(tk,Sn,Sb,B0,B1,B2,B3);
             s [2][rhs] = bmma/dt*(U - U_o)
              - barma*Udw
              - XComponent(tk,Sn,Sb,B0,B1,B2,B3)
              + wet + bdr*fabs(Ur)*Ur - xthrust;
         }

         ftr = bmma/dt + 2.0*bdr*sign(Vr)*Vr;

         s [3][1] = -tdk*dYdt(B0,B1,B2,B3);
         s [3][2] = -dYdn(B0,B1,B2,B3);
         s [3][3] = -dYdb(B0,B1,B2,B3);
         s [3][4] = ftr*dYdt(B0,B1,B2,B3);
         s [3][5] = ftr*dYdn(B0,B1,B2,B3);
         s [3][6] = ftr*dYdb(B0,B1,B2,B3);
         s [3][7] = ftr*dYdB0(u,v,w,B0,B1,B2,B3) - dYdB0(tk,Sn,Sb,B0,B1,B2,B3);
         s [3][8] = ftr*dYdB1(u,v,w,B0,B1,B2,B3) - dYdB1(tk,Sn,Sb,B0,B1,B2,B3);
         s [3][9] = ftr*dYdB2(u,v,w,B0,B1,B2,B3) - dYdB2(tk,Sn,Sb,B0,B1,B2,B3);
         s [3][10] = ftr*dYdB3(u,v,w,B0,B1,B2,B3) - dYdB3(tk,Sn,Sb,B0,B1,B2,B3);
         s [3][rhs] = bmma/dt*(V - V_o) - ythrust
              - barma*Vdw - YComponent(tk,Sn,Sb,B0,B1,B2,B3) + bdr*fabs(Vr)*Vr;

         ftr = bmma/dt + 2.0*bdr*sign(Wr)*Wr;

         s [4][1] = -tdk*dZdt(B0,B1,B2,B3);
         s [4][2] = -dZdn(B0,B1,B2,B3);
         s [4][3] = -dZdb(B0,B1,B2,B3);
         s [4][4] = ftr*dZdt(B0,B1,B2,B3);
         s [4][5] = ftr*dZdn(B0,B1,B2,B3);
         s [4][6] = ftr*dZdb(B0,B1,B2,B3);
         s [4][7] = ftr*dZdB0(u,v,w,B0,B1,B2,B3) - dZdB0(tk,Sn,Sb,B0,B1,B2,B3);
         s [4][8] = ftr*dZdB1(u,v,w,B0,B1,B2,B3) - dZdB1(tk,Sn,Sb,B0,B1,B2,B3);
         s [4][9] = ftr*dZdB2(u,v,w,B0,B1,B2,B3) - dZdB2(tk,Sn,Sb,B0,B1,B2,B3);
         s [4][10] = ftr*dZdB3(u,v,w,B0,B1,B2,B3) - dZdB3(tk,Sn,Sb,B0,B1,B2,B3);
         s [4][rhs] = bmma/dt*(W - W_o) - zthrust
              - barma*Wdw - ZComponent(tk,Sn,Sb,B0,B1,B2,B3) + bdr*fabs(Wr)*Wr;
      }

      s [5][11] = 1.0;
      s [5][rhs] = Om1;

      s [6][12] = 1.0;
      s [6][rhs] = Om2;

      s [7][13] = 1.0;
      s [7][rhs] = Om3;
   }

	/* 	
	 * buoy node
	 */

   else if (eq_type == TopBoundary) {

      s [1][7]  = 2.0*B0;
      s [1][8]  = 2.0*B1;
      s [1][9]  = 2.0*B2;
      s [1][10] = 2.0*B3;
      s [1][rhs] = B0*B0 + B1*B1 + B2*B2 + B3*B3 - 1.0;
/*
      s [1][11] = 1;
      s [1][rhs] = Om1;
*/
      for (i = 1 ; i <= 10 ; i++) {
         s [2][i] = 0.0;
         s [3][i] = 0.0;
         s [4][i] = 0.0;
      }

      if (problem -> terminal [2] -> release > 0.0 &&
          tm >= problem -> terminal [2] -> release && !problem -> dynstat) {


         s [2][1] = 1.0;
         s [2][rhs] = e;

         s [3][2] = 1.0;
         s [3][rhs] = Sn;

         s [4][3] = 1.0;
         s [4][rhs] = Sb;
      }
      else if (problem -> terminal [2] -> anchor) {

         s [2][4]    = 1.0;
         s [2][rhs]  = u;

         s [3][5]    = 1.0;
         s [3][rhs]  = v;

         s [4][6]    = 1.0;
         s [4][rhs]  = w;
      }
      else if (problem -> type == Towing) {

         Speed (tm, problem -> terminal [2], &Uspd, &Vspd, &Wspd);
         InputVelocity (tm, n -> x, n -> y, n -> z, &Uw, &Vw, &Ww);
         

         U = Uspd + Uw;
         V = Vspd + Vw;
         W = Wspd + Ww;

         s [2][4]   = 1.0;
         s [2][7]   = -dtdB0(U,V,W,B0,B1,B2,B3);
         s [2][8]   = -dtdB1(U,V,W,B0,B1,B2,B3);
         s [2][9]   = -dtdB2(U,V,W,B0,B1,B2,B3);
         s [2][10]   = -dtdB3(U,V,W,B0,B1,B2,B3);
         s [2][rhs] = u - tComponent(U,V,W,B0,B1,B2,B3);

         s [3][5]   = 1.0;
         s [3][7]   = -dndB0(U,V,W,B0,B1,B2,B3);
         s [3][8]   = -dndB1(U,V,W,B0,B1,B2,B3);
         s [3][9]   = -dndB2(U,V,W,B0,B1,B2,B3);
         s [3][10]   = -dndB3(U,V,W,B0,B1,B2,B3);
         s [3][rhs] = v - nComponent(U,V,W,B0,B1,B2,B3);

         s [4][6]   = 1.0;
         s [4][7]   = -dbdB0(U,V,W,B0,B1,B2,B3);
         s [4][8]   = -dbdB1(U,V,W,B0,B1,B2,B3);
         s [4][9]   = -dbdB2(U,V,W,B0,B1,B2,B3);
         s [4][10]   = -dbdB3(U,V,W,B0,B1,B2,B3);
         s [4][rhs] = w - bComponent(U,V,W,B0,B1,B2,B3);
      }
      else if (forcing == Force) {

         InputVelocity (tm, n -> x, n -> y, n -> z, &Uw, &Vw, &Ww);

         xforst = problem -> terminal [2] -> xforce + Uw;
         yforst = problem -> terminal [2] -> yforce + Vw;
         zforst = problem -> terminal [2] -> xforce + Ww;

         s [2][1]  = -tdk*dXdt(B0,B1,B2,B3);
         s [2][2]  = -dXdn(B0,B1,B2,B3);
         s [2][3]  = -dXdb(B0,B1,B2,B3);            
         s [2][7]  = -dXdB0(tk,Sn,Sb,B0,B1,B2,B3);
         s [2][8]  = -dXdB1(tk,Sn,Sb,B0,B1,B2,B3);
         s [2][9]  = -dXdB2(tk,Sn,Sb,B0,B1,B2,B3);
         s [2][10]  = -dXdB3(tk,Sn,Sb,B0,B1,B2,B3);
         s [2][rhs] = xforst - XComponent(tk,Sn,Sb,B0,B1,B2,B3);

         s [3][1]  = -tdk*dYdt(B0,B1,B2,B3);
         s [3][2]  = -dYdn(B0,B1,B2,B3);
         s [3][3]  = -dYdb(B0,B1,B2,B3);            
         s [3][7]  = -dYdB0(tk,Sn,Sb,B0,B1,B2,B3);
         s [3][8]  = -dYdB1(tk,Sn,Sb,B0,B1,B2,B3);
         s [3][9]  = -dYdB2(tk,Sn,Sb,B0,B1,B2,B3);
         s [3][10]  = -dYdB3(tk,Sn,Sb,B0,B1,B2,B3);
         s [3][rhs] = yforst - YComponent(tk,Sn,Sb,B0,B1,B2,B3);

         s [4][1]  = -tdk*dZdt(B0,B1,B2,B3);
         s [4][2]  = -dZdn(B0,B1,B2,B3);
         s [4][3]  = -dZdb(B0,B1,B2,B3);            
         s [4][7]  = -dZdB0(tk,Sn,Sb,B0,B1,B2,B3);
         s [4][8]  = -dZdB1(tk,Sn,Sb,B0,B1,B2,B3);
         s [4][9]  = -dZdB2(tk,Sn,Sb,B0,B1,B2,B3);
         s [4][10]  = -dZdB3(tk,Sn,Sb,B0,B1,B2,B3);
         s [4][rhs] = zforst - ZComponent(tk,Sn,Sb,B0,B1,B2,B3);
      }
      else if (forcing == Morison || problem -> type == Deployment
               || problem -> type == HorizontalDrifter
               || ((problem -> type == Surface || problem -> type == Subsurface)
                   && problem -> dynstat)) {

         buoy = problem -> terminal [2] -> buoy;
         bmma  = buoy -> Mmma; 

         if (problem -> type == Deployment || problem -> dynstat) {
            buoy -> draft = environment -> surface - n -> x;

            barma = 0.0;
            bdr     =  0.5*environment -> rho*buoy -> Cdn*
                       ProjectedArea(buoy, buoy -> draft);

	        Uw = Vw = Ww = 0.0;
	        Udw = Vdw = Wdw = 0.0;
         }
         else {
            bdr   = buoy -> Mdr;
            barma = buoy -> Marma;

            WaveParticleVelocity (tm, n -> x, n -> y, n -> z, &Uw, &Vw, &Ww);
            WaveParticleAcceleration (tm, n -> x, n -> y, n -> z, &Udw, &Vdw, &Wdw);
         }

         if (problem -> type == HorizontalDrifter) {
            Thrust(tm, problem -> terminal [2], &xthrust, &ythrust, &zthrust);
            problem -> terminal[2] -> xthrust.value = xthrust;
            problem -> terminal[2] -> ythrust.value = ythrust;
            problem -> terminal[2] -> zthrust.value = zthrust;
         }
         else {
            xthrust = 0.0;
            ythrust = 0.0;
            zthrust = 0.0;
         }

         if (problem -> type == Surface)
            WindDrag (tm, buoy, &Fy_wind, &Fz_wind);
         else {
            Fy_wind = 0.0;
            Fz_wind = 0.0;
         }

         Current (tm, n -> x, n -> y, n -> z, &uack, &vack, &wack);

         B = Buoyancy(buoy, environment -> depth - n -> x, environment);

         U   = XComponent(u,v,w,B0,B1,B2,B3);
         U_o = XComponent(u_o,v_o,w_o,B0_o,B1_o,B2_o,B3_o);

         V   = YComponent(u,v,w,B0,B1,B2,B3);
         V_o = YComponent(u_o,v_o,w_o,B0_o,B1_o,B2_o,B3_o);

         W   = ZComponent(u,v,w,B0,B1,B2,B3);
         W_o = ZComponent(u_o,v_o,w_o,B0_o,B1_o,B2_o,B3_o);


         Ur   = U - Uw;
         Vr   = V - vack - Vw;
         Wr   = W - wack - Ww;

         ftr = bmma/dt + 2.0*bdr*sign(Ur)*Ur;
         if (problem -> type == HorizontalDrifter) {
            s [2][4]  = dXdt(B0,B1,B2,B3);
            s [2][5]  = dXdn(B0,B1,B2,B3);
            s [2][6]  = dXdb(B0,B1,B2,B3);
            s [2][7]  = dXdB0(u,v,w,B0,B1,B2,B3);
            s [2][8]  = dXdB1(u,v,w,B0,B1,B2,B3);
            s [2][9]  = dXdB2(u,v,w,B0,B1,B2,B3);
            s [2][10]  = dXdB3(u,v,w,B0,B1,B2,B3);
            s [2][rhs] = XComponent(u,v,w,B0,B1,B2,B3);
         }
         else {
            s [2][1] = tdk*dXdt(B0,B1,B2,B3);
            s [2][2] = dXdn(B0,B1,B2,B3);
            s [2][3] = dXdb(B0,B1,B2,B3);
            s [2][4] = ftr*dXdt(B0,B1,B2,B3);
            s [2][5] = ftr*dXdn(B0,B1,B2,B3);
            s [2][6] = ftr*dXdb(B0,B1,B2,B3);
            s [2][7] = ftr*dXdB0(u,v,w,B0,B1,B2,B3) + dXdB0(tk,Sn,Sb,B0,B1,B2,B3);
            s [2][8] = ftr*dXdB1(u,v,w,B0,B1,B2,B3) + dXdB1(tk,Sn,Sb,B0,B1,B2,B3);
            s [2][9] = ftr*dXdB2(u,v,w,B0,B1,B2,B3) + dXdB2(tk,Sn,Sb,B0,B1,B2,B3);
            s [2][10] = ftr*dXdB3(u,v,w,B0,B1,B2,B3) + dXdB3(tk,Sn,Sb,B0,B1,B2,B3);
            s [2][rhs] = bmma/dt*(U - U_o)
              - barma*Udw + XComponent(tk,Sn,Sb,B0,B1,B2,B3)
              - (B - buoy -> w) + bdr*fabs(Ur)*Ur - xthrust;
         }
         ftr = bmma/dt + 2.0*bdr*sign(Vr)*Vr;

         s [3][1] = tdk*dYdt(B0,B1,B2,B3);
         s [3][2] = dYdn(B0,B1,B2,B3);
         s [3][3] = dYdb(B0,B1,B2,B3);
         s [3][4] = ftr*dYdt(B0,B1,B2,B3);
         s [3][5] = ftr*dYdn(B0,B1,B2,B3);
         s [3][6] = ftr*dYdb(B0,B1,B2,B3);
         s [3][7] = ftr*dYdB0(u,v,w,B0,B1,B2,B3) + dYdB0(tk,Sn,Sb,B0,B1,B2,B3);
         s [3][8] = ftr*dYdB1(u,v,w,B0,B1,B2,B3) + dYdB1(tk,Sn,Sb,B0,B1,B2,B3);
         s [3][9] = ftr*dYdB2(u,v,w,B0,B1,B2,B3) + dYdB2(tk,Sn,Sb,B0,B1,B2,B3);
         s [3][10] = ftr*dYdB3(u,v,w,B0,B1,B2,B3) + dYdB3(tk,Sn,Sb,B0,B1,B2,B3);
         s [3][rhs] = bmma/dt*(V - V_o)
              - barma*Vdw + YComponent(tk,Sn,Sb,B0,B1,B2,B3) + bdr*fabs(Vr)*Vr
              - Fy_wind - ythrust;

         ftr = bmma/dt + 2.0*bdr*sign(Wr)*Wr;

         s [4][1] = tdk*dZdt(B0,B1,B2,B3);
         s [4][2] = dZdn(B0,B1,B2,B3);
         s [4][3] = dZdb(B0,B1,B2,B3);
         s [4][4] = ftr*dZdt(B0,B1,B2,B3);
         s [4][5] = ftr*dZdn(B0,B1,B2,B3);
         s [4][6] = ftr*dZdb(B0,B1,B2,B3);
         s [4][7] = ftr*dZdB0(u,v,w,B0,B1,B2,B3) + dZdB0(tk,Sn,Sb,B0,B1,B2,B3);
         s [4][8] = ftr*dZdB1(u,v,w,B0,B1,B2,B3) + dZdB1(tk,Sn,Sb,B0,B1,B2,B3);
         s [4][9] = ftr*dZdB2(u,v,w,B0,B1,B2,B3) + dZdB2(tk,Sn,Sb,B0,B1,B2,B3);
         s [4][10] = ftr*dZdB3(u,v,w,B0,B1,B2,B3) + dZdB3(tk,Sn,Sb,B0,B1,B2,B3);
         s [4][rhs] = bmma/dt*(W - W_o)
              - barma*Wdw + ZComponent(tk,Sn,Sb,B0,B1,B2,B3) + bdr*fabs(Wr)*Wr
              - Fz_wind - zthrust;
      }
      else if (forcing == WaveFollower || forcing == Velocity || forcing == LAMP) {

         if (forcing == WaveFollower) 
            WaveSurfaceVelocity (tm, n -> x, n -> y, n -> z, &Uw, &Vw, &Ww);
         else 
            InputVelocity (tm, n -> x, n -> y, n -> z, &Uw, &Vw, &Ww);

	/*
 	 * Just in case someone is trying something funny with
	 * a speed specification
	 */

         Speed (tm, problem -> terminal [2], &Uspd, &Vspd, &Wspd);

         Uw += Uspd;
         Vw += Vspd;
         Ww += Wspd;
         fprintf(stderr,"Uw = %g, Vw = %g, Ww =%g\n", Uw, Vw, Ww);

         s [2][4]  = -dXdt(B0,B1,B2,B3);
         s [2][5]  = -dXdn(B0,B1,B2,B3);
         s [2][6]  = -dXdb(B0,B1,B2,B3);
         s [2][7]  = -dXdB0(u,v,w,B0,B1,B2,B3);
         s [2][8]  = -dXdB1(u,v,w,B0,B1,B2,B3);
         s [2][9]  = -dXdB2(u,v,w,B0,B1,B2,B3);
         s [2][10]  = -dXdB3(u,v,w,B0,B1,B2,B3);
         s [2][rhs] = Uw - XComponent(u,v,w,B0,B1,B2,B3);

         if (forcing == WaveFollower && 
             (problem -> type == Surface || problem -> type == Drifter)) {

            buoy = problem -> terminal [2] -> buoy;
            buoy -> draft = environment -> surface - n -> x;

            bmma  = buoy -> Mmma; 
            bdr   = 0.5*environment -> rho*buoy -> Cdn*ProjectedArea(buoy, buoy -> draft);

            Current (tm, n -> x, n -> y, n -> z, &uack, &vack, &wack);
            Current (tm - dt, n -> x_o, n -> y_o, n -> z_o, &uacok, &vacok, &wacok);
            WindDrag(tm, buoy, &Fy_wind, &Fz_wind);

            V   = YComponent(u,v,w,B0,B1,B2,B3);
            V_o = YComponent(u_o,v_o,w_o,B0_o,B1_o,B2_o,B3_o);

            W   = ZComponent(u,v,w,B0,B1,B2,B3);
            W_o = ZComponent(u_o,v_o,w_o,B0_o,B1_o,B2_o,B3_o);


            Vr   = V - vack;
            Wr   = W - wack;

            ftr = bmma/dt + 2.0*bdr*sign(Vr)*Vr;

            s [3][1] = tdk*dYdt(B0,B1,B2,B3);
            s [3][2] = dYdn(B0,B1,B2,B3);
            s [3][3] = dYdb(B0,B1,B2,B3);
            s [3][4] = ftr*dYdt(B0,B1,B2,B3);
            s [3][5] = ftr*dYdn(B0,B1,B2,B3);
            s [3][6] = ftr*dYdb(B0,B1,B2,B3);
            s [3][7] = ftr*dYdB0(u,v,w,B0,B1,B2,B3) 
			+ dYdB0(tk,Sn,Sb,B0,B1,B2,B3);
            s [3][8] = ftr*dYdB1(u,v,w,B0,B1,B2,B3) 
			+ dYdB1(tk,Sn,Sb,B0,B1,B2,B3);
            s [3][9] = ftr*dYdB2(u,v,w,B0,B1,B2,B3) 
			+ dYdB2(tk,Sn,Sb,B0,B1,B2,B3);
            s [3][10] = ftr*dYdB3(u,v,w,B0,B1,B2,B3) 
			+ dYdB3(tk,Sn,Sb,B0,B1,B2,B3);
            s [3][rhs] = bmma/dt*(V - V_o)
                 + YComponent(tk,Sn,Sb,B0,B1,B2,B3)
                 + bdr*fabs(Vr)*Vr - Fy_wind;

            ftr = bmma/dt + 2.0*bdr*sign(Wr)*Wr;

            s [4][1] = tdk*dZdt(B0,B1,B2,B3);
            s [4][2] = dZdn(B0,B1,B2,B3);
            s [4][3] = dZdb(B0,B1,B2,B3);
            s [4][4] = ftr*dZdt(B0,B1,B2,B3);
            s [4][5] = ftr*dZdn(B0,B1,B2,B3);
            s [4][6] = ftr*dZdb(B0,B1,B2,B3);
            s [4][7] = ftr*dZdB0(u,v,w,B0,B1,B2,B3)
                       + dZdB0(tk,Sn,Sb,B0,B1,B2,B3);
            s [4][8] = ftr*dZdB1(u,v,w,B0,B1,B2,B3)
                       + dZdB1(tk,Sn,Sb,B0,B1,B2,B3);
            s [4][9] = ftr*dZdB2(u,v,w,B0,B1,B2,B3)
                       + dZdB2(tk,Sn,Sb,B0,B1,B2,B3);
            s [4][10] = ftr*dZdB3(u,v,w,B0,B1,B2,B3)
                       + dZdB3(tk,Sn,Sb,B0,B1,B2,B3);
            s [4][rhs] = bmma/dt*(W - W_o)
                + ZComponent(tk,Sn,Sb,B0,B1,B2,B3)
                + bdr*fabs(Wr)*Wr - Fz_wind;
         }
         else {
            s [3][4]  = -dYdt(B0,B1,B2,B3);
            s [3][5]  = -dYdn(B0,B1,B2,B3);
            s [3][6]  = -dYdb(B0,B1,B2,B3);
            s [3][7]  = -dYdB0(u,v,w,B0,B1,B2,B3);
            s [3][8]  = -dYdB1(u,v,w,B0,B1,B2,B3);
            s [3][9]  = -dYdB2(u,v,w,B0,B1,B2,B3);
            s [3][10]  = -dYdB3(u,v,w,B0,B1,B2,B3);
            s [3][rhs] = Vw - YComponent(u,v,w,B0,B1,B2,B3);

            s [4][4]  = -dZdt(B0,B1,B2,B3);
            s [4][5]  = -dZdn(B0,B1,B2,B3);
            s [4][6]  = -dZdb(B0,B1,B2,B3);
            s [4][7]  = -dZdB0(u,v,w,B0,B1,B2,B3);
            s [4][8]  = -dZdB1(u,v,w,B0,B1,B2,B3);
            s [4][9]  = -dZdB2(u,v,w,B0,B1,B2,B3);
            s [4][10]  = -dZdB3(u,v,w,B0,B1,B2,B3);
            s [4][rhs] = Ww - ZComponent(u,v,w,B0,B1,B2,B3);
         }
      }
      else if (forcing == FroudeKrylov) {
         buoy = problem -> terminal [2] -> buoy;

         bmma = buoy -> m + buoy -> am;
         bdr  = buoy -> Mdr; 

         FroudeKrylovCoefficients (tm, buoy, n -> x, n -> y, n -> z,
                                   Fex, &wave_b, &wet);

         WaveParticleVelocity (tm, n -> x, n -> y, n -> z, &Uw, &Vw, &Ww);

         Current (tm, n -> x, n -> y, n -> z, &uack, &vack, &wack);

         U   = XComponent(u,v,w,B0,B1,B2,B3);
         U_o = XComponent(u_o,v_o,w_o,B0_o,B1_o,B2_o,B3_o);

         V   = YComponent(u,v,w,B0,B1,B2,B3);
         V_o = YComponent(u_o,v_o,w_o,B0_o,B1_o,B2_o,B3_o);

         W   = ZComponent(u,v,w,B0,B1,B2,B3);
         W_o = ZComponent(u_o,v_o,w_o,B0_o,B1_o,B2_o,B3_o);

         Ur   = U - Uw;
         Vr   = V - vack - Vw;
         Wr   = W - wack - Ww;

         ftr = bmma/dt + 2.0*bdr*sign(Ur)*Ur + wave_b;

         s [2][1] = tdk*dXdt(B0,B1,B2,B3);
         s [2][2] = dXdn(B0,B1,B2,B3);
         s [2][3] = dXdb(B0,B1,B2,B3);
         s [2][4] = ftr*dXdt(B0,B1,B2,B3);
         s [2][5] = ftr*dXdn(B0,B1,B2,B3);
         s [2][6] = ftr*dXdb(B0,B1,B2,B3);
         s [2][7] = ftr*dXdB0(u,v,w,B0,B1,B2,B3) + dXdB0(tk,Sn,Sb,B0,B1,B2,B3);
         s [2][8] = ftr*dXdB1(u,v,w,B0,B1,B2,B3) + dXdB1(tk,Sn,Sb,B0,B1,B2,B3);
         s [2][9] = ftr*dXdB2(u,v,w,B0,B1,B2,B3) + dXdB2(tk,Sn,Sb,B0,B1,B2,B3);
         s [2][10] = ftr*dXdB3(u,v,w,B0,B1,B2,B3) + dXdB3(tk,Sn,Sb,B0,B1,B2,B3);
         s [2][rhs] = bmma/dt*(U - U_o)
              - Fex [1] + wave_b*U 
              + XComponent(tk,Sn,Sb,B0,B1,B2,B3)
              - wet + bdr*fabs(Ur)*Ur;

         s [3][4]   = dYdt(B0,B1,B2,B3);
         s [3][5]   = dYdn(B0,B1,B2,B3);
         s [3][6]   = dYdb(B0,B1,B2,B3);
         s [3][7]   = dYdB0(u,v,w,B0,B1,B2,B3);
         s [3][8]   = dYdB1(u,v,w,B0,B1,B2,B3);
         s [3][9]   = dYdB2(u,v,w,B0,B1,B2,B3);
         s [3][10]  = dYdB3(u,v,w,B0,B1,B2,B3);
         s [3][rhs] = YComponent(u,v,w,B0,B1,B2,B3);

         s [4][4]   = dZdt(B0,B1,B2,B3);
         s [4][5]   = dZdn(B0,B1,B2,B3);
         s [4][6]   = dZdb(B0,B1,B2,B3);
         s [4][7]   = dZdB0(u,v,w,B0,B1,B2,B3);
         s [4][8]   = dZdB1(u,v,w,B0,B1,B2,B3);
         s [4][9]   = dZdB2(u,v,w,B0,B1,B2,B3);
         s [4][10]  = dZdB3(u,v,w,B0,B1,B2,B3);
         s [4][rhs] = ZComponent(u,v,w,B0,B1,B2,B3);
      }

      s [5][12]   = 1.0;
      s [5][rhs]  = Om2;

      s [6][13]   = 1.0;
      s [6][rhs]  = Om3;
   } 
   else if (eq_type == BranchStart) {
      s [1][7] = 1.0;
      s [1][rhs] = B0; // - n -> Ys[4]; // B0;

      s [2][11] = 1.0;
      s [2][rhs] = Om1;

      s [3][12] = 1.0;
      s [3][rhs] = Om2;

      s [4][13] = 1.0;
      s [4][rhs] = Om3;
   }
   else if (eq_type == BranchTerminal) {
      s [1][7]  = 2.0*B0;
      s [1][8]  = 2.0*B1;
      s [1][9]  = 2.0*B2;
      s [1][10] = 2.0*B3;
      s [1][rhs] = B0*B0 + B1*B1 + B2*B2 + B3*B3 - 1.0;

      s [2][12]   = 1.0;
      s [2][rhs]  = Om2;

      s [3][13]   = 1.0;
      s [3][rhs]  = Om3;

      if (n -> segment -> branch -> terminal -> buoy) {
         buoy = n -> segment -> branch -> terminal -> buoy;

         bmma  = buoy -> Mmma;
         bdr   = buoy -> Mdr;

         Current (tm, n -> x, n -> y, n -> z, &uack, &vack, &wack);

         if (forcing != Velocity && forcing != Force) {
            barma = buoy -> Marma;
            WaveParticleVelocity (tm, n -> x, n -> y, n -> z, &Uw, &Vw, &Ww);
            WaveParticleAcceleration (tm, n -> x, n -> y, n -> z,
                                      &Udw, &Vdw, &Wdw);
         }
         else {
	        barma = 0.0;
            Uw = Vw = Ww = 0.0;
            Udw = Vdw = Wdw = 0.0;
         }

         wet = buoy -> buoyancy - buoy -> w;

         U   = XComponent(u,v,w,B0,B1,B2,B3);
         U_o = XComponent(u_o,v_o,w_o,B0_o,B1_o,B2_o,B3_o);
         V   = YComponent(u,v,w,B0,B1,B2,B3);
         V_o = YComponent(u_o,v_o,w_o,B0_o,B1_o,B2_o,B3_o);
         W   = ZComponent(u,v,w,B0,B1,B2,B3);
         W_o = ZComponent(u_o,v_o,w_o,B0_o,B1_o,B2_o,B3_o);

         Ur   = U - Uw;
         Vr   = V - vack - Vw;
         Wr   = W - wack - Ww;

         ftr = bmma/dt + 2.0*sign(Ur)*Ur;

         s [4][1] = tdk*dXdt(B0,B1,B2,B3);
         s [4][2] = dXdn(B0,B1,B2,B3);
         s [4][3] = dXdb(B0,B1,B2,B3);
         s [4][4] = ftr*dXdt(B0,B1,B2,B3);
         s [4][5] = ftr*dXdn(B0,B1,B2,B3);
         s [4][6] = ftr*dXdb(B0,B1,B2,B3);
         s [4][7] = ftr*dXdB0(u,v,w,B0,B1,B2,B3) 
                    + ftr*dXdB0(tk,Sn,Sb,B0,B1,B2,B3);
         s [4][8] = ftr*dXdB1(u,v,w,B0,B1,B2,B3) 
                    + ftr*dXdB1(tk,Sn,Sb,B0,B1,B2,B3);
         s [4][9] = ftr*dXdB2(u,v,w,B0,B1,B2,B3) 
                    + ftr*dXdB2(tk,Sn,Sb,B0,B1,B2,B3);
         s [4][10] = ftr*dXdB3(u,v,w,B0,B1,B2,B3) 
                    + ftr*dXdB3(tk,Sn,Sb,B0,B1,B2,B3);
         s [4][rhs] = bmma/dt*(U - U_o)
              - barma*Udw
              + XComponent(tk,Sn,Sb,B0,B1,B2,B3)
              - wet + bdr*fabs(Ur)*Ur;

         ftr = bmma/dt + 2.0*sign(Vr)*Vr;

         s [5][1] = tdk*dYdt(B0,B1,B2,B3);
         s [5][2] = dYdn(B0,B1,B2,B3);
         s [5][3] = dYdb(B0,B1,B2,B3);
         s [5][4] = ftr*dYdt(B0,B1,B2,B3);
         s [5][5] = ftr*dYdn(B0,B1,B2,B3);
         s [5][6] = ftr*dYdb(B0,B1,B2,B3);
         s [5][7] = ftr*dYdB0(u,v,w,B0,B1,B2,B3)
                    + dYdB0(tk,Sn,Sb,B0,B1,B2,B3);
         s [5][8] = ftr*dYdB1(u,v,w,B0,B1,B2,B3)
                    + dYdB1(tk,Sn,Sb,B0,B1,B2,B3);
         s [5][9] = ftr*dYdB2(u,v,w,B0,B1,B2,B3)
                    + dYdB2(tk,Sn,Sb,B0,B1,B2,B3);
         s [5][10] = ftr*dYdB3(u,v,w,B0,B1,B2,B3)
                    + dYdB3(tk,Sn,Sb,B0,B1,B2,B3);
         s [5][rhs] = bmma/dt*(V - V_o)
              - barma*Vdw
              + YComponent(tk,Sn,Sb,B0,B1,B2,B3)
              + bdr*fabs(Vr)*Vr;

         ftr = bmma/dt + 2.0*sign(Wr)*Wr;

         s [6][1] = tdk*dZdt(B0,B1,B2,B3);
         s [6][2] = dZdn(B0,B1,B2,B3);
         s [6][3] = dZdb(B0,B1,B2,B3);
         s [6][4] = ftr*dZdt(B0,B1,B2,B3);
         s [6][5] = ftr*dZdn(B0,B1,B2,B3);
         s [6][6] = ftr*dZdb(B0,B1,B2,B3);
         s [6][7] = ftr*dZdB0(u,v,w,B0,B1,B2,B3)
                    + dZdB0(tk,Sn,Sb,B0,B1,B2,B3);
         s [6][8] = ftr*dZdB1(u,v,w,B0,B1,B2,B3)
                    + dZdB1(tk,Sn,Sb,B0,B1,B2,B3);
         s [6][9] = ftr*dZdB2(u,v,w,B0,B1,B2,B3)
                    + dZdB2(tk,Sn,Sb,B0,B1,B2,B3);
         s [6][10] = ftr*dZdB3(u,v,w,B0,B1,B2,B3)
                    + dZdB3(tk,Sn,Sb,B0,B1,B2,B3);
         s [6][rhs] = bmma/dt*(W - W_o)
              - barma*Wdw
              + YComponent(tk,Sn,Sb,B0,B1,B2,B3)
              + bdr*fabs(Wr)*Wr;
      }
      else if (n -> segment -> branch -> terminal -> anchor) {
         s [4][4]    = 1.0;
         s [4][rhs]  = u;

         s [5][5]    = 1.0;
         s [5][rhs]  = v;

         s [6][6]    = 1.0;
         s [6][rhs]  = w;
      }
      else { // if (n -> segment -> branch -> terminal -> node) {
	/* empty -- compatibility already expressed at Junction */
      }

   }


	/*
	 * node between two distinct segments
	 */

   else if (eq_type == Connection || eq_type == Junction) {
      if (nm -> segment -> connector == NULL
          && nm -> segment -> connection == Spliced) {  /* no lumped mass     */
         s [1][1]      = -tdkm;
         s [1][ne + 1] = tdk;
         s [1][rhs]    = tk - tkm;

         for (i = 2 ; i <= ne ; i++) {
            s [i][i] = -1.0;
            s [i][i + ne] = 1.0;
            s [i][rhs] = n -> Y[i] - nm -> Y[i];
         }
      }
      else {				/* lumped mass is present */
         s [1][ne + 7] = 1.0;
         s [1][rhs] = B0; // - n -> Ys[4]; // 0; // B0; 

         s [2][ne + 11] = 1.0;		/* release the moments 	*/	
         s [2][rhs] = Om1;		/* above the mass	*/

         s [3][ne + 12] = 1.0;
         s [3][rhs] = Om2;

         s [4][ne + 13] = 1.0;
         s [4][rhs] = Om3;

         s [5][7]  = 2.0*B0_m;	/* enforce Euler condition */	
         s [5][8]  = 2.0*B1_m;	/* below the mass	   */
         s [5][9]  = 2.0*B2_m;
         s [5][10] = 2.0*B3_m;
         s [5][rhs] = B0_m*B0_m + B1_m*B1_m + B2_m*B2_m + B3_m*B3_m - 1.0;

         s [6][12] = 1.0;		        /* release the moments 	   */
         s [6][rhs] = Om2_m;		/* below the mass	   */

         s [7][13] = 1.0;
         s [7][rhs] = Om3_m;

         c = nm -> segment -> connector;
         if (c) {
            mma  = c -> m + c -> am;
            dr_n   =  c -> Cdn;
            dr_t   =  c -> Cdt;
            wet  = c -> wet;
            if (problem -> type == Deployment  || problem -> type == Surface)
               if (n -> x > environment -> surface && wet < 0.0)
                  wet = wet*(1.0 + tanh(50.0*(environment -> surface - n -> x)));
         }
         else {
            dr_n = dr_t = wet = mma = 0;
         }


	/*
	 * HACK HACK HACK -- DEOS spar buoy hack -- HACK HACK HACK
	 */
/*
         if (problem -> type == Horizontal) {
            wet = -(5059.3*(5000.0 - n -> x) - 33366.0);
            dr = 0.5*environment -> rho*1.0*0.8*(5000.0 - n -> x);
         }
*/
	 nav = 2.0;

         U    = XComponent(u,   v,   w,   B0,   B1,   B2,   B3);
         U_m  = XComponent(u_m, v_m, w_m, B0_m, B1_m, B2_m, B3_m);
         U_o  = XComponent(u_o, v_o, w_o, B0_o, B1_o, B2_o, B3_o);
         U_om = XComponent(u_om,v_om,w_om,B0_om,B1_om,B2_om,B3_om);

         Ua   = (U + U_m)/nav;
         Ua_o = (U_o + U_om)/nav; 

         V    = YComponent(u,   v,   w,   B0,   B1,   B2,   B3);
         V_m  = YComponent(u_m, v_m, w_m, B0_m, B1_m, B2_m, B3_m);
         V_o  = YComponent(u_o, v_o, w_o, B0_o, B1_o, B2_o, B3_o);
         V_om = YComponent(u_om,v_om,w_om,B0_om,B1_om,B2_om,B3_om);
         
         Va   = (V + V_m)/nav;
         Va_o = (V_o + V_om)/nav;

         W    = ZComponent(u,   v,   w,   B0,   B1,   B2,   B3);
         W_m  = ZComponent(u_m, v_m, w_m, B0_m, B1_m ,B2_m, B3_m);
         W_o  = ZComponent(u_o, v_o, w_o, B0_o, B1_o, B2_o, B3_o);
         W_om = ZComponent(u_om,v_om,w_om,B0_om,B1_om,B2_om,B3_om);
                     
         Wa   = (W + W_m)/nav;
         Wa_o = (W_o + W_om)/nav;

	     njn = n -> segment -> junction.num_nodes;
         if (eq_type == Junction) {
	        Wa = Wa*nav;
	        Wa_o = Wa_o*nav;

	        Va = Va*nav;
	        Va_o = Va_o*nav;

	        Ua = Ua*nav;
	        Ua_o = Ua_o*nav;

            for (i = 1 ; i <= njn ; i++) {
               nj = n -> segment -> junction.node [i];
               
               u_j = nj -> Y [4];
               v_j = nj -> Y [5];
               w_j = nj -> Y [6];
               u_oj = nj -> Y_o [4];
               v_oj = nj -> Y_o [5];
               w_oj = nj -> Y_o [6];

               B0_j = nj -> Y [7];
               B1_j = nj -> Y [8];
               B2_j = nj -> Y [9];
               B3_j = nj -> Y [10];
               B0_oj = nj -> Y_o [7];
               B1_oj = nj -> Y_o [8];
               B2_oj = nj -> Y_o [9];
               B3_oj = nj -> Y_o [10];

               U_j  = XComponent(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j);
               U_oj = XComponent(u_oj,v_oj,w_oj,B0_oj,B1_oj,B2_oj,B3_oj);

               V_j  = YComponent(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j);
               V_oj = YComponent(u_oj,v_oj,w_oj,B0_oj,B1_oj,B2_oj,B3_oj);

               W_j  = ZComponent(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j);
               W_oj = ZComponent(u_oj,v_oj,w_oj,B0_oj,B1_oj,B2_oj,B3_oj);


	           Wa += W_j;
	           Wa_o += W_oj;

	           Va += V_j;
	           Va_o += V_oj;

	           Ua += U_j;
	           Ua_o += U_oj;
            }

            nav += njn;
 
            Wa = Wa / nav;
            Wa_o = Wa_o / nav;

            Va = Va / nav;
            Va_o = Va_o / nav;

            Ua = Ua / nav;
            Ua_o = Ua_o / nav;
         }

         Current (tm, n -> x, n -> y, n -> z, &uack, &vack, &wack);
         if (problem -> type == Deployment && problem -> dynstat) {
            vack -= problem -> terminal [1] -> yspeed.value;
            wack -= problem -> terminal [1] -> zspeed.value;
         }
        
         ConnectorThrust(tm, nm -> segment, nm, &xthrust, &ythrust, &zthrust); 
         
         nm -> segment -> connector_xthrust.value = xthrust;
         nm -> segment -> connector_ythrust.value = ythrust;
         nm -> segment -> connector_zthrust.value = zthrust;
         
         if (forcing != Velocity && forcing != Force) {
            WaveParticleAcceleration (tm, n -> x, n -> y, n -> z, 
                                      &Udw, &Vdw, &Wdw);
            WaveParticleVelocity (tm, n -> x, n -> y, n -> z, 
                                      &Uw, &Vw, &Ww);
 
            carma   = 1.5*pow(c -> d, 3.0)*M_PI*environment -> rho/6.0;
         }
         else {
            carma = 0.0;
            Uw = Vw = Ww = 0.0;
            Udw = Vdw = Wdw = 0.0;
         }

         Ur   = U - Uw;
         Vr   = V - vack - Vw;
         Wr   = W - wack - Ww;

	    // sum of forces in the global X direction

         if (problem -> type == HorizontalDrifter) {
             s [8][4]      = dXdt(B0_m,B1_m,B2_m,B3_m);
             s [8][5]      = dXdn(B0_m,B1_m,B2_m,B3_m);
             s [8][6]      = dXdb(B0_m,B1_m,B2_m,B3_m);
             s [8][7]      = dXdB0(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m);
             s [8][8]      = dXdB1(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m);
             s [8][9]      = dXdB2(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m);
             s [8][10]      = dXdB3(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m);
             s [8][rhs] = U_m;
         }
         else {
            ftr = (mma/dt + 2.0*dr_t*sign(Ur)*Ur)/nav;
            s [8][1]      = tdkm*dXdt(B0_m,B1_m,B2_m,B3_m);
            s [8][ne + 1] = -tdk*dXdt(B0,B1,B2,B3);
            s [8][2]      = dXdn(B0_m,B1_m,B2_m,B3_m);
            s [8][ne + 2] = -dXdn(B0,B1,B2,B3);
            s [8][3]      = dXdb(B0_m,B1_m,B2_m,B3_m);
            s [8][ne + 3] = -dXdb(B0,B1,B2,B3);
            s [8][4]      = ftr*dXdt(B0_m,B1_m,B2_m,B3_m);
            s [8][ne + 4] = ftr*dXdt(B0,B1,B2,B3);
            s [8][5]      = ftr*dXdn(B0_m,B1_m,B2_m,B3_m);
            s [8][ne + 5] = ftr*dXdn(B0,B1,B2,B3);
            s [8][6]      = ftr*dXdb(B0_m,B1_m,B2_m,B3_m);
            s [8][ne + 6] = ftr*dXdb(B0,B1,B2,B3);
            s [8][7]      = ftr*dXdB0(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m)
                 + dXdB0(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
            s [8][ne + 7] = ftr*dXdB0(u,v,w,B0,B1,B2,B3) 
                - dXdB0(tk,Sn,Sb,B0,B1,B2,B3);
            s [8][8]      = ftr*dXdB1(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m)
                + dXdB1(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
            s [8][ne + 8] = ftr*dXdB1(u,v,w,B0,B1,B2,B3) 
                - dXdB1(tk,Sn,Sb,B0,B1,B2,B3);
            s [8][9]      = ftr*dXdB2(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m)
                + dXdB2(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
            s [8][ne + 9] = ftr*dXdB2(u,v,w,B0,B1,B2,B3) 
                - dXdB2(tk,Sn,Sb,B0,B1,B2,B3);
            s [8][10]      = ftr*dXdB3(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m)
                + dXdB3(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
            s [8][ne + 10] = ftr*dXdB3(u,v,w,B0,B1,B2,B3) 
                - dXdB3(tk,Sn,Sb,B0,B1,B2,B3);
            s [8][rhs] = mma*(Ua - Ua_o)/dt + wet - xthrust
             - carma*Udw 
             + dr_t*(Ur)*fabs(Ur) 
             + XComponent(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m)
             - XComponent(tk,Sn,Sb,B0,B1,B2,B3);
         }
  
	/*
	 * sum of forces in the global y direction
	 */

         ftr = (mma/dt + 2.0*dr_n*sign(Vr)*Vr)/nav;
         
         s [9][1]      = tdkm*dYdt(B0_m,B1_m,B2_m,B3_m);
         s [9][ne + 1] = -tdk*dYdt(B0,B1,B2,B3);
         s [9][2]      = dYdn(B0_m,B1_m,B2_m,B3_m);
         s [9][ne + 2] = -dYdn(B0,B1,B2,B3);
         s [9][3]      = dYdb(B0_m,B1_m,B2_m,B3_m);
         s [9][ne + 3] = -dYdb(B0,B1,B2,B3);
         s [9][4]      = ftr*dYdt(B0_m,B1_m,B2_m,B3_m);
         s [9][ne + 4] = ftr*dYdt(B0,B1,B2,B3);
         s [9][5]      = ftr*dYdn(B0_m,B1_m,B2_m,B3_m);
         s [9][ne + 5] = ftr*dYdn(B0,B1,B2,B3);
         s [9][6]      = ftr*dYdb(B0_m,B1_m,B2_m,B3_m);
         s [9][ne + 6] = ftr*dYdb(B0,B1,B2,B3);
         s [9][7]      = ftr*dYdB0(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m)
                + dYdB0(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
         s [9][ne + 7] = ftr*dYdB0(u,v,w,B0,B1,B2,B3) 
                - dYdB0(tk,Sn,Sb,B0,B1,B2,B3);
         s [9][8]      = ftr*dYdB1(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m)
                + dYdB1(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
         s [9][ne + 8] = ftr*dYdB1(u,v,w,B0,B1,B2,B3) 
                - dYdB1(tk,Sn,Sb,B0,B1,B2,B3);
         s [9][9]      = ftr*dYdB2(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m)
                + dYdB2(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
         s [9][ne + 9] = ftr*dYdB2(u,v,w,B0,B1,B2,B3) 
                - dYdB2(tk,Sn,Sb,B0,B1,B2,B3);
         s [9][10]      = ftr*dYdB3(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m)
                + dYdB3(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
         s [9][ne + 10] = ftr*dYdB3(u,v,w,B0,B1,B2,B3) 
                - dYdB3(tk,Sn,Sb,B0,B1,B2,B3);
         s [9][rhs] = mma*(Va - Va_o)/dt 
          - carma*Vdw - ythrust
          + dr_n*(Vr)*fabs(Vr) 
          + YComponent(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m)
          - YComponent(tk,Sn,Sb,B0,B1,B2,B3);

	/*
	 * sum of forces in the global Z direction	
	 */

         ftr = (mma/dt + 2.0*dr_n*sign(Wr)*Wr)/nav;
         
         s [10][1]      = tdkm*dZdt(B0_m,B1_m,B2_m,B3_m);
         s [10][ne + 1] = -tdk*dZdt(B0,B1,B2,B3);
         s [10][2]      = dZdn(B0_m,B1_m,B2_m,B3_m);
         s [10][ne + 2] = -dZdn(B0,B1,B2,B3);
         s [10][3]      = dZdb(B0_m,B1_m,B2_m,B3_m);
         s [10][ne + 3] = -dZdb(B0,B1,B2,B3);
         s [10][4]      = ftr*dZdt(B0_m,B1_m,B2_m,B3_m);
         s [10][ne + 4] = ftr*dZdt(B0,B1,B2,B3);
         s [10][5]      = ftr*dZdn(B0_m,B1_m,B2_m,B3_m);
         s [10][ne + 5] = ftr*dZdn(B0,B1,B2,B3);
         s [10][6]      = ftr*dZdb(B0_m,B1_m,B2_m,B3_m);
         s [10][ne + 6] = ftr*dZdb(B0,B1,B2,B3);
         s [10][7]      = ftr*dZdB0(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m)
                + dZdB0(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
         s [10][ne + 7] = ftr*dZdB0(u,v,w,B0,B1,B2,B3) 
                - dZdB0(tk,Sn,Sb,B0,B1,B2,B3);
         s [10][8]      = ftr*dZdB1(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m)
                + dZdB1(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
         s [10][ne + 8] = ftr*dZdB1(u,v,w,B0,B1,B2,B3) 
                - dZdB1(tk,Sn,Sb,B0,B1,B2,B3);
         s [10][9]      = ftr*dZdB2(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m)
                + dZdB2(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
         s [10][ne + 9] = ftr*dZdB2(u,v,w,B0,B1,B2,B3) 
                - dZdB2(tk,Sn,Sb,B0,B1,B2,B3);
         s [10][10]      = ftr*dZdB3(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m)
                + dZdB3(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
         s [10][ne + 10] = ftr*dZdB3(u,v,w,B0,B1,B2,B3) 
                - dZdB3(tk,Sn,Sb,B0,B1,B2,B3);
         s [10][rhs] = mma*(Wa - Wa_o)/dt 
          - carma*Wdw - zthrust
          + dr_n*(Wr)*fabs(Wr) 
          + ZComponent(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m)
          - ZComponent(tk,Sn,Sb,B0,B1,B2,B3);

	/*
	 * enforce compatibility of the X, Y, Z velocity of the
	 * nodes above and below the mass
	 */

         s [11][4]      = -dXdt(B0_m,B1_m,B2_m,B3_m);
         s [11][ne + 4] = dXdt(B0,B1,B2,B3);
         s [11][5]      = -dXdn(B0_m,B1_m,B2_m,B3_m);
         s [11][ne + 5] = dXdn(B0,B1,B2,B3);
         s [11][6]      = -dXdb(B0_m,B1_m,B2_m,B3_m);
         s [11][ne + 6] = dXdb(B0,B1,B2,B3);
         s [11][7]      = -dXdB0(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m);
         s [11][ne + 7] = dXdB0(u,v,w,B0,B1,B2,B3);
         s [11][8]      = -dXdB1(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m);
         s [11][ne + 8] = dXdB1(u,v,w,B0,B1,B2,B3);
         s [11][9]      = -dXdB2(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m);
         s [11][ne + 9] = dXdB2(u,v,w,B0,B1,B2,B3);
         s [11][10]      = -dXdB3(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m);
         s [11][ne + 10] = dXdB3(u,v,w,B0,B1,B2,B3);
         s [11][rhs] = U - U_m;

         s [12][4]      = -dYdt(B0_m,B1_m,B2_m,B3_m);
         s [12][ne + 4] = dYdt(B0,B1,B2,B3);
         s [12][5]      = -dYdn(B0_m,B1_m,B2_m,B3_m);
         s [12][ne + 5] = dYdn(B0,B1,B2,B3);
         s [12][6]      = -dYdb(B0_m,B1_m,B2_m,B3_m);
         s [12][ne + 6] = dYdb(B0,B1,B2,B3);
         s [12][7]      = -dYdB0(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m);
         s [12][ne + 7] = dYdB0(u,v,w,B0,B1,B2,B3);
         s [12][8]      = -dYdB1(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m);
         s [12][ne + 8] = dYdB1(u,v,w,B0,B1,B2,B3);
         s [12][9]      = -dYdB2(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m);
         s [12][ne + 9] = dYdB2(u,v,w,B0,B1,B2,B3);
         s [12][10]      = -dYdB3(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m);
         s [12][ne + 10] = dYdB3(u,v,w,B0,B1,B2,B3);
         s [12][rhs] = V - V_m;

         s [13][4]      = -dZdt(B0_m,B1_m,B2_m,B3_m);
         s [13][ne + 4] = dZdt(B0,B1,B2,B3);
         s [13][5]      = -dZdn(B0_m,B1_m,B2_m,B3_m);
         s [13][ne + 5] = dZdn(B0,B1,B2,B3);
         s [13][6]      = -dZdb(B0_m,B1_m,B2_m,B3_m);
         s [13][ne + 6] = dZdb(B0,B1,B2,B3);
         s [13][7]      = -dZdB0(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m);
         s [13][ne + 7] = dZdB0(u,v,w,B0,B1,B2,B3);
         s [13][8]      = -dZdB1(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m);
         s [13][ne + 8] = dZdB1(u,v,w,B0,B1,B2,B3);
         s [13][9]      = -dZdB2(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m);
         s [13][ne + 9] = dZdB2(u,v,w,B0,B1,B2,B3);
         s [13][10]      = -dZdB3(u_m,v_m,w_m,B0_m,B1_m,B2_m,B3_m);
         s [13][ne + 10] = dZdB3(u,v,w,B0,B1,B2,B3);
         s [13][rhs] = W - W_m;

         j = 14;

         for (i = 1 ; i <= njn ; i++) {
            nj = n -> segment -> junction.node [i];
            
            e_j    = nj -> Y [1];
            Sn_j   = nj -> Y [2];
            Sb_j   = nj -> Y [3];
            u_j    = nj -> Y [4];
            v_j    = nj -> Y [5];
            w_j    = nj -> Y [6];
            B0_j   = nj -> Y [7];
            B1_j   = nj -> Y [8];
            B2_j   = nj -> Y [9];
            B3_j   = nj -> Y [10];

            //e_oj    = nj -> Y_o [1];
            //Sn_oj   = nj -> Y_o [2];
            //Sb_oj   = nj -> Y_o [3];
            u_oj    = nj -> Y_o [4];
            v_oj    = nj -> Y_o [5];
            w_oj    = nj -> Y_o [6];
            B0_oj   = nj -> Y_o [7];
            B1_oj   = nj -> Y_o [8];
            B2_oj   = nj -> Y_o [9];
            B3_oj   = nj -> Y_o [10];

            tkj   = Tension(e_j, nj -> material);
            // tokj  = Tension(e_oj, nj -> material);
            tdkj  = TensionD(e_j, nj -> material);

	    if (nj == nj -> segment -> last)
               sign_j = -1.0;
	    else
	       sign_j = 1.0;

            ftr = (mma/dt + 2.0*dr_t*sign(Ur)*Ur)/nav;

            s [8][(i + 1)*ne + 1] = -sign_j*tdkj*dXdt(B0_j,B1_j,B2_j,B3_j);
            s [8][(i + 1)*ne + 2] = -sign_j*dXdn(B0_j,B1_j,B2_j,B3_j);
            s [8][(i + 1)*ne + 3] = -sign_j*dXdb(B0_j,B1_j,B2_j,B3_j);
            s [8][(i + 1)*ne + 4] = ftr*dXdt(B0_j,B1_j,B2_j,B3_j);
            s [8][(i + 1)*ne + 5] = ftr*dXdn(B0_j,B1_j,B2_j,B3_j);
            s [8][(i + 1)*ne + 6] = ftr*dXdb(B0_j,B1_j,B2_j,B3_j);
            s [8][(i + 1)*ne + 7] = ftr*dXdB0(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j)
                - sign_j*dXdB0(tkj,Sn_j,Sb_j,B0_j,B1_j,B2_j,B3_j);
            s [8][(i + 1)*ne + 8] = ftr*dXdB1(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j)
                - sign_j*dXdB1(tkj,Sn_j,Sb_j,B0_j,B1_j,B2_j,B3_j);
            s [8][(i + 1)*ne + 9] = ftr*dXdB2(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j)
                - sign_j*dXdB2(tkj,Sn_j,Sb_j,B0_j,B1_j,B2_j,B3_j);
            s [8][(i + 1)*ne + 10] = ftr*dXdB3(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j)
                - sign_j*dXdB3(tkj,Sn_j,Sb_j,B0_j,B1_j,B2_j,B3_j);

            s [8][rhs] += -sign_j*XComponent(tkj,Sn_j,Sb_j,B0_j,B1_j,B2_j,B3_j);

            ftr = (mma/dt + 2.0*dr_n*sign(Vr)*Vr)/nav;

            s [9][(i + 1)*ne + 1] = -sign_j*tdkj*dYdt(B0_j,B1_j,B2_j,B3_j);
            s [9][(i + 1)*ne + 2] = -sign_j*dYdn(B0_j,B1_j,B2_j,B3_j);
            s [9][(i + 1)*ne + 3] = -sign_j*dYdb(B0_j,B1_j,B2_j,B3_j);
            s [9][(i + 1)*ne + 4] = ftr*dYdt(B0_j,B1_j,B2_j,B3_j);
            s [9][(i + 1)*ne + 5] = ftr*dYdn(B0_j,B1_j,B2_j,B3_j);
            s [9][(i + 1)*ne + 6] = ftr*dYdb(B0_j,B1_j,B2_j,B3_j);
            s [9][(i + 1)*ne + 7] = ftr*dYdB0(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j)
                - sign_j*dYdB0(tkj,Sn_j,Sb_j,B0_j,B1_j,B2_j,B3_j);
            s [9][(i + 1)*ne + 8] = ftr*dYdB1(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j)
                - sign_j*dYdB1(tkj,Sn_j,Sb_j,B0_j,B1_j,B2_j,B3_j);
            s [9][(i + 1)*ne + 9] = ftr*dYdB2(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j)
                - sign_j*dYdB2(tkj,Sn_j,Sb_j,B0_j,B1_j,B2_j,B3_j);
            s [9][(i + 1)*ne + 10] = ftr*dYdB3(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j)
                - sign_j*dYdB3(tkj,Sn_j,Sb_j,B0_j,B1_j,B2_j,B3_j);

            s [9][rhs] += -sign_j*YComponent(tkj,Sn_j,Sb_j,B0_j,B1_j,B2_j,B3_j);

            ftr = (mma/dt + 2.0*dr_n*sign(Wr)*Wr)/nav;

            s [10][(i + 1)*ne + 1] = -sign_j*tdkj*dZdt(B0_j,B1_j,B2_j,B3_j);
            s [10][(i + 1)*ne + 2] = -sign_j*dZdn(B0_j,B1_j,B2_j,B3_j);
            s [10][(i + 1)*ne + 3] = -sign_j*dZdb(B0_j,B1_j,B2_j,B3_j);
            s [10][(i + 1)*ne + 4] = ftr*dZdt(B0_j,B1_j,B2_j,B3_j);
            s [10][(i + 1)*ne + 5] = ftr*dZdn(B0_j,B1_j,B2_j,B3_j);
            s [10][(i + 1)*ne + 6] = ftr*dZdb(B0_j,B1_j,B2_j,B3_j);
            s [10][(i + 1)*ne + 7] = ftr*dZdB0(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j)
                - sign_j*dZdB0(tkj,Sn_j,Sb_j,B0_j,B1_j,B2_j,B3_j);
            s [10][(i + 1)*ne + 8] = ftr*dZdB1(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j)
                - sign_j*dZdB1(tkj,Sn_j,Sb_j,B0_j,B1_j,B2_j,B3_j);
            s [10][(i + 1)*ne + 9] = ftr*dZdB2(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j)
                - sign_j*dZdB2(tkj,Sn_j,Sb_j,B0_j,B1_j,B2_j,B3_j);
            s [10][(i + 1)*ne + 10] = ftr*dZdB3(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j)
                - sign_j*dZdB3(tkj,Sn_j,Sb_j,B0_j,B1_j,B2_j,B3_j);

            s [10][rhs] += -sign_j*ZComponent(tkj,Sn_j,Sb_j,B0_j,B1_j,B2_j,B3_j);

            U_j = XComponent(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j);
            V_j = YComponent(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j);
            W_j = ZComponent(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j);

            s [j][ne + 4]         = dXdt(B0,B1,B2,B3);
            s [j][(i + 1)*ne + 4] = -dXdt(B0_j,B1_j,B2_j,B3_j);
            s [j][ne + 5]         = dXdn(B0,B1,B2,B3);
            s [j][(i + 1)*ne + 5] = -dXdn(B0_j,B1_j,B2_j,B3_j);
            s [j][ne + 6]         = dXdb(B0,B1,B2,B3);
            s [j][(i + 1)*ne + 6] = -dXdb(B0_j,B1_j,B2_j,B3_j);
            s [j][ne + 7]         = dXdB0(u,v,w,B0,B1,B2,B3);
            s [j][(i + 1)*ne + 7] = -dXdB0(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j);
            s [j][ne + 8]         = dXdB1(u,v,w,B0,B1,B2,B3);
            s [j][(i + 1)*ne + 8] = -dXdB1(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j);
            s [j][ne + 9]         = dXdB2(u,v,w,B0,B1,B2,B3);
            s [j][(i + 1)*ne + 9] = -dXdB2(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j);
            s [j][ne + 10]         = dXdB3(u,v,w,B0,B1,B2,B3);
            s [j][(i + 1)*ne + 10] = -dXdB3(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j);
            s [j][rhs] = U - U_j;
      	    j ++;

            s [j][ne + 4]         = dYdt(B0,B1,B2,B3);
            s [j][(i + 1)*ne + 4] = -dYdt(B0_j,B1_j,B2_j,B3_j);
            s [j][ne + 5]         = dYdn(B0,B1,B2,B3);
            s [j][(i + 1)*ne + 5] = -dYdn(B0_j,B1_j,B2_j,B3_j);
            s [j][ne + 6]         = dYdb(B0,B1,B2,B3);
            s [j][(i + 1)*ne + 6] = -dYdb(B0_j,B1_j,B2_j,B3_j);
            s [j][ne + 7]         = dYdB0(u,v,w,B0,B1,B2,B3);
            s [j][(i + 1)*ne + 7] = -dYdB0(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j);
            s [j][ne + 8]         = dYdB1(u,v,w,B0,B1,B2,B3);
            s [j][(i + 1)*ne + 8] = -dYdB1(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j);
            s [j][ne + 9]         = dYdB2(u,v,w,B0,B1,B2,B3);
            s [j][(i + 1)*ne + 9] = -dYdB2(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j);
            s [j][ne + 10]         = dYdB3(u,v,w,B0,B1,B2,B3);
            s [j][(i + 1)*ne + 10] = -dYdB3(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j);
            s [j][rhs] = V - V_j;
      	    j ++;

            s [j][ne + 4]         = dZdt(B0,B1,B2,B3);
            s [j][(i + 1)*ne + 4] = -dZdt(B0_j,B1_j,B2_j,B3_j);
            s [j][ne + 5]         = dZdn(B0,B1,B2,B3);
            s [j][(i + 1)*ne + 5] = -dZdn(B0_j,B1_j,B2_j,B3_j);
            s [j][ne + 6]         = dZdb(B0,B1,B2,B3);
            s [j][(i + 1)*ne + 6] = -dZdb(B0_j,B1_j,B2_j,B3_j);
            s [j][ne + 7]         = dZdB0(u,v,w,B0,B1,B2,B3);
            s [j][(i + 1)*ne + 7] = -dZdB0(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j);
            s [j][ne + 8]         = dZdB1(u,v,w,B0,B1,B2,B3);
            s [j][(i + 1)*ne + 8] = -dZdB1(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j);
            s [j][ne + 9]         = dZdB2(u,v,w,B0,B1,B2,B3);
            s [j][(i + 1)*ne + 9] = -dZdB2(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j);
            s [j][ne + 10]         = dZdB3(u,v,w,B0,B1,B2,B3);
            s [j][(i + 1)*ne + 10] = -dZdB3(u_j,v_j,w_j,B0_j,B1_j,B2_j,B3_j);
            s [j][rhs] = W - W_j;
      	    j ++;
         }
      }
   }

	/*
	 * a regular old internal node
	 */

   else {
      gdt = analysis -> gamma*dt;
      t_ds = 2.0/ds;


      omak2_ds = omak*omak/ds;
      akak2_ds = ak*omak/ds;
      ak2_ds   = ak*ak / ds;


      e_d_o  = n -> Yd_o [1];
      //Sn_d_o = n -> Yd_o [2];
      //Sb_d_o = n -> Yd_o [3];
      u_d_o  = n -> Yd_o [4];
      v_d_o  = n -> Yd_o [5];
      w_d_o  = n -> Yd_o [6];
      B0_d_o = n -> Yd_o [7];
      B1_d_o = n -> Yd_o [8];
      B2_d_o = n -> Yd_o [9];
      B3_d_o = n -> Yd_o [10];

      e_d_om  = nm -> Yd_o [1];
      // Sn_d_om = nm -> Yd_o [2];
      // Sb_d_om = nm -> Yd_o [3];
      u_d_om  = nm -> Yd_o [4];
      v_d_om  = nm -> Yd_o [5];
      w_d_om  = nm -> Yd_o [6];
      B0_d_om = nm -> Yd_o [7];
      B1_d_om = nm -> Yd_o [8];
      B2_d_om = nm -> Yd_o [9];
      B3_d_om = nm -> Yd_o [10];

      e_d  = n -> Yd[1] = ((e - e_o)/gdt - omg_g*e_d_o);
      //Sn_d = n -> Yd[2] = ((Sn - Sn_o)/gdt - omg_g*Sn_d_o);
      //Sb_d = n -> Yd[3] = ((Sb - Sb_o)/gdt - omg_g*Sb_d_o);
      u_d  = n -> Yd[4] = ((u - u_o)/gdt - omg_g*u_d_o);
      v_d  = n -> Yd[5] = ((v - v_o)/gdt - omg_g*v_d_o);
      w_d  = n -> Yd[6] = ((w - w_o)/gdt - omg_g*w_d_o);
      B0_d = n -> Yd[7] = ((B0 - B0_o)/gdt - omg_g*B0_d_o);
      B1_d = n -> Yd[8] = ((B1 - B1_o)/gdt - omg_g*B1_d_o);
      B2_d = n -> Yd[9] = ((B2 - B2_o)/gdt - omg_g*B2_d_o);
      B3_d = n -> Yd[10] = ((B3 - B3_o)/gdt - omg_g*B3_d_o);

      e_d_m  = nm -> Yd[1] = ((e_m - e_om)/gdt - omg_g*e_d_om);
      // Sn_d_m = nm -> Yd[2] = ((Sn_m - Sn_om)/gdt - omg_g*Sn_d_om);
      // Sb_d_m = nm -> Yd[3] = ((Sb_m - Sb_om)/gdt - omg_g*Sb_d_om);
      u_d_m  = nm -> Yd[4] = ((u_m - u_om)/gdt - omg_g*u_d_om);
      v_d_m  = nm -> Yd[5] = ((v_m - v_om)/gdt - omg_g*v_d_om);
      w_d_m  = nm -> Yd[6] = ((w_m - w_om)/gdt - omg_g*w_d_om);
      B0_d_m = nm -> Yd[7] = ((B0_m - B0_om)/gdt - omg_g*B0_d_om);
      B1_d_m = nm -> Yd[8] = ((B1_m - B1_om)/gdt - omg_g*B1_d_om);
      B2_d_m = nm -> Yd[9] = ((B2_m - B2_om)/gdt - omg_g*B2_d_om);
      B3_d_m = nm -> Yd[10] = ((B3_m - B3_om)/gdt - omg_g*B3_d_om);

      Current (tm, n -> x, n -> y, n -> z, &uack, &vack, &wack);
      Current (tm - dt, n -> x_o, n -> y_o, n -> z_o, &uacok, &vacok, &wacok);
      if (problem -> type == Deployment && problem -> dynstat) {
         vack -= problem -> terminal [1] -> yspeed.value;
         wack -= problem -> terminal [1] -> zspeed.value;
         vacok -= problem -> terminal [1] -> yspeed.value;
         wacok -= problem -> terminal [1] -> zspeed.value;
      }


      if (forcing == WaveFollower || forcing == Morison || forcing == LAMP) {
         WaveParticleVelocity (tm, n -> x, n -> y, n -> z, &Uw, &Vw, &Ww);
         WaveParticleAcceleration (tm, n -> x, n -> y, n -> z,
                                   &Udw, &Vdw, &Wdw);

         WaveParticleVelocity (tm - dt, n -> x_o, n -> y_o, n -> z_o,
                                                       &Uw_o, &Vw_o, &Ww_o);
         WaveParticleAcceleration (tm - dt, n -> x_o, n -> y_o, n -> z_o,
                                                       &Udw_o, &Vdw_o, &Wdw_o);

         //uldk = tComponent(Udw,Vdw,Wdw,B0,B1,B2,B3);
         //vldk = nComponent(Udw,Vdw,Wdw,B0,B1,B2,B3);
         //wldk = bComponent(Udw,Vdw,Wdw,B0,B1,B2,B3);

         //uldok = tComponent(Udw_o,Vdw_o,Wdw_o,B0_o,B1_o,B2_o,B3_o);
         //vldok = nComponent(Udw_o,Vdw_o,Wdw_o,B0_o,B1_o,B2_o,B3_o);
         //wldok = bComponent(Udw_o,Vdw_o,Wdw_o,B0_o,B1_o,B2_o,B3_o);
      }
      else {
         Uw = Vw = Ww = Uw_o = Vw_o = Ww_o = 0.0;
         Udw = Vdw = Wdw = Udw_o = Vdw_o = Wdw_o = 0.0;

         // uldk = vldk = wldk = uldok = vldok = wldok = 0.0;
      }

      uack += Uw;	uacok += Uw_o;
      vack += Vw;	vacok += Vw_o;
      wack += Ww;	wacok += Ww_o;

      ulck = tComponent(uack,vack,wack,B0,B1,B2,B3);
      vlck = nComponent(uack,vack,wack,B0,B1,B2,B3);
      wlck = bComponent(uack,vack,wack,B0,B1,B2,B3);

      ulcok = tComponent(uacok,vacok,wacok,B0_o,B1_o,B2_o,B3_o);
      vlcok = nComponent(uacok,vacok,wacok,B0_o,B1_o,B2_o,B3_o);
      wlcok = bComponent(uacok,vacok,wacok,B0_o,B1_o,B2_o,B3_o);

      if (n -> active_number == 2 || nm -> active_number != n -> active_number - 1) {
         Current (tm,      nm -> x,   nm -> y,   nm -> z,   &uackm,  &vackm,  &wackm);
         Current (tm - dt, nm -> x_o, nm -> y_o, nm -> z_o, &uacokm, &vacokm, &wacokm);

         if (problem -> type == Deployment && problem -> dynstat) {
            vackm -= problem -> terminal [1] -> yspeed.value;
            wackm -= problem -> terminal [1] -> zspeed.value;
            vacokm -= problem -> terminal [1] -> yspeed.value;
            wacokm -= problem -> terminal [1] -> zspeed.value;
         }

         if (forcing == WaveFollower || forcing == Morison || forcing == LAMP) {
            WaveParticleVelocity (tm, nm -> x, nm -> y, nm -> z,
                                  &Uw_m, &Vw_m, &Ww_m);
            WaveParticleAcceleration (tm, nm -> x, nm -> y, nm -> z,
                                      &Udw_m, &Vdw_m, &Wdw_m);

            WaveParticleVelocity (tm - dt, nm -> x_o, nm -> y_o, nm -> z_o,
	   		          &Uw_om, &Vw_om, &Ww_om);
            WaveParticleAcceleration (tm - dt, nm -> x_o, nm -> y_o, nm -> z_o,
			              &Udw_om, &Vdw_om, &Wdw_om);

            // uldkm = tComponent(Udw_m,Vdw_m,Wdw_m,B0_m,B1_m,B2_m,B3_m);
            // vldkm = nComponent(Udw_m,Vdw_m,Wdw_m,B0_m,B1_m,B2_m,B3_m);
            // wldkm = bComponent(Udw_m,Vdw_m,Wdw_m,B0_m,B1_m,B2_m,B3_m);

            // uldokm = tComponent(Udw_om,Vdw_om,Wdw_om,B0_om,B1_om,B2_om,B3_om);
            // vldokm = nComponent(Udw_om,Vdw_om,Wdw_om,B0_om,B1_om,B2_om,B3_om);
            // wldokm = bComponent(Udw_om,Vdw_om,Wdw_om,B0_om,B1_om,B2_om,B3_om);
         }
         else {
            Uw_m = Vw_m = Ww_m = Uw_om = Vw_om = Ww_om = 0.0;
            Udw_m = Vdw_m = Wdw_m = Udw_om = Vdw_om = Wdw_om = 0.0;
 
            // vldkm = wldkm = vldokm = wldokm = 0.0;
         }

         uackm += Uw_m;		uacokm += Uw_om;
         vackm += Vw_m;		vacokm += Vw_om;
         wackm += Ww_m;		wacokm += Ww_om;

         ulckm = tComponent(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m);
         vlckm = nComponent(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m);
         wlckm = bComponent(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m);

         ulcokm = tComponent(uacokm,vacokm,wacokm,B0_om,B1_om,B2_om,B3_om);
         vlcokm = nComponent(uacokm,vacokm,wacokm,B0_om,B1_om,B2_om,B3_om);
         wlcokm = bComponent(uacokm,vacokm,wacokm,B0_om,B1_om,B2_om,B3_om);
      } 
      else {
         ulckm  = ulcpre;
         vlckm  = vlcpre;
         wlckm  = wlcpre;
         uackm  = uacpre;
         vackm  = vacpre;
         wackm  = wacpre;

         ulcokm  = ulcopre;
         vlcokm  = vlcopre;
         wlcokm  = wlcopre;
         uacokm  = uacopre;
         vacokm  = vacopre;
         wacokm  = wacopre;

         // uldkm = uldpre;
         // vldkm = vldpre;
         // wldkm = wldpre;
 
         // uldokm = uldopre;
         // vldokm = vldopre;
         // wldokm = wldopre;
      }

      ulcpre = ulck;
      vlcpre = vlck;
      wlcpre = wlck;
      uacpre = uack;
      vacpre = vack;
      wacpre = wack;

      // uldpre = uldk;
      // vldpre = vldk;
      // wldpre = wldk;

      ulcopre = ulcok;
      vlcopre = vlcok;
      wlcopre = wlcok;
      uacopre = uacok;
      vacopre = vacok;
      wacopre = wacok;

      // uldopre = uldok;
      // vldopre = vldok;
      // wldopre = wldok;

      urk  = u - ulck;
      vrk  = v - vlck;
      wrk  = w - wlck;
      urkm = u_m - ulckm;
      vrkm = v_m - vlckm;
      wrkm = w_m - wlckm;

      urok = u_o - ulcok;
      vrok = v_o - vlcok;
      wrok = w_o - wlcok;
      urokm = u_om - ulcokm;
      vrokm = v_om - vlcokm;
      wrokm = w_om - wlcokm;

      urka = fabs(urk);
      urkma = fabs(urkm);

      uroka = fabs(urok);
      urokma = fabs(urokm);

      pervk  = sqrt(vrk*vrk + wrk*wrk);
      pervkm = sqrt(vrkm*vrkm + wrkm*wrkm);
      pervok  = sqrt(vrok*vrok + wrok*wrok);
      pervokm = sqrt(vrokm*vrokm + wrokm*wrokm);

      DragCoeff(n, n -> material, tm, urka, pervk, &drat, &drap);
      DragCoeff(nm, nm -> material, tm, urkma, pervkm, &drat_m, &drap_m);
      DragCoeff(n, n -> material, tm-dt, uroka, pervok, &drat_o, &drap_o);
      DragCoeff(nm, nm -> material, tm-dt, urokma, pervokm, &drat_om, &drap_om);

      //printf("a %f %f %f %f %f %f %f %f : %f %f %f %f %f %f %f %f\n",
      //       drat, drat_m, drat_o, drat_om,
      //       drap, drap_m, drap_o, drap_om,
      //       urk, urkm, urok, urokm, pervk, pervkm, pervok, pervokm);
 
      if (n -> attachment) {
         w0     += n -> attachment -> wet/ds;
         drat   += n -> attachment -> Cdt/ds;
         drap   += n -> attachment -> Cdn/ds;
         drat_o += n -> attachment -> Cdt/ds;
         drap_o += n -> attachment -> Cdn/ds;
         cam    += n -> attachment -> m/ds;
         camma_n  += (n -> attachment -> m + n -> attachment -> am)/ds;
         camma_t  += (n -> attachment -> m + n -> attachment -> am)/ds;
         caarma_n += n -> attachment -> am/ds;
      }
      if (nm -> attachment) {
         w0_m     += nm -> attachment -> wet/ds;
         drat_m   += nm -> attachment -> Cdt/ds;
         drap_m   += nm -> attachment -> Cdn/ds;
         drat_om  += nm -> attachment -> Cdt/ds;
         drap_om  += nm -> attachment -> Cdn/ds;
         cam_m    += nm -> attachment -> m/ds;
         camma_t_m  += (nm -> attachment -> m + nm -> attachment -> am)/ds;
         camma_n_m  += (nm -> attachment -> m + nm -> attachment -> am)/ds;
         caarma_n_m += nm -> attachment -> am/ds;
      }

      if (pervk < TOLERANCE) {
         drap = 0.0;
         pervk = 1.0;
      }

      if (pervok < TOLERANCE) {
         drap_o = 0.0;
         pervok = 1.0;
      }
  
      if (pervokm < TOLERANCE) {
         drap_om = 0.0;
         pervokm = 1.0;
      }

      if (pervkm < TOLERANCE) {
         drap_m = 0.0;
         pervkm = 1.0;
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


      //printf("b %f %f %f %f %f %f %f %f : %f %f %f %f %f %f %f %f\n",
      //       drat, drat_m, drat_o, drat_om,
      //       drap, drap_m, drap_o, drap_om,
      //       urk, urkm, urok, urokm, pervk, pervkm, pervok, pervokm);

      ftr   = drat*urk*sign(urk)*sqstrk;
      ftr_m = drat_m*urkm*sign(urkm)*sqstrkm;
      
      s [1][1]       = omak2_ds*(tddkm*(e - e_m) - (tdk + tdkm))
                       + akak2_ds*(tddkm*(e_o - e_om) - (tdok + tdokm))
                       - 0.5*omak*drat_m*(urkm)*fabs(urkm)/sqstrkm;
      s [1][ne + 1]  = omak2_ds*(tddk*(e - e_m) + (tdk + tdkm))
                       + akak2_ds*(tddk*(e_o - e_om) + (tdok + tdokm))
                       - 0.5*omak*drat*(urk)*fabs(urk)/sqstrk;
      s [1][2]       = -omak*Om3_m;
      s [1][ne + 2]  = -omak*Om3;
      s [1][3]       = omak*Om2_m;
      s [1][ne + 3]  = omak*Om2; 
      s [1][4]       = -camma_t_m/gdt*(omam2 + amam2) 
                 - omak*(mud_bm + damp_t + 2.0*ftr_m)
                 + omak2_ds*(pay + pay_m)
                 + akak2_ds*(pay_om + pay_o);
      s [1][ne + 4]  = -camma_t/gdt*(omam2 + amam2)
                 - omak*(mud_b + damp_t + 2.0*ftr)
                 - omak2_ds*(pay + pay_m)
                 - akak2_ds*(pay_o + pay_om);
      s [1][5]       = omak*pay_m*Om3_m 
                       -2.0*cam_m*(omam2*(B3_m*B0_d_m - B2_m*B1_d_m 
                                          + B1_m*B2_d_m - B0_m*B3_d_m)
                                   + amam2*(B3_m*B0_d_om - B2_m*B1_d_om
                                            + B1_m*B2_d_om - B0_m*B3_d_om));
      s [1][ne + 5]  = omak*pay*Om3
                       -2.0*cam*(omam2*(B3*B0_d - B2*B1_d + B1*B2_d - B0*B3_d)
                                 + amam2*(B3*B0_d_o - B2*B1_d_o
                                          + B1*B2_d_o - B0*B3_d_o));
      s [1][6]       = -omak*pay_m*Om2_m 
                       -2.0*cam_m*(omam2*(-B2_m*B0_d_m - B3_m*B1_d_m
                                          + B0_m*B2_d_m + B1_m*B3_d_m)
                                   + amam2*(-B2_m*B0_d_om - B3_m*B1_d_om
                                            + B0_m*B2_d_om + B1_m*B3_d_om));
      s [1][ne + 6]  = -omak*pay*Om2 
                       -2.0*cam*(omam2*(-B2*B0_d - B3*B1_d
                                        + B0*B2_d + B1*B3_d)
                                 + amam2*(-B2*B0_d_o - B3*B1_d_o
                                          + B0*B2_d_o + B1*B3_d_o));
      s [1][7]       = -2.0*omam2*(cam_m*((v_m*B3_m - w_m*B2_m)/gdt
                                           + w_m*B2_d_m - v_m*B3_d_m)
                 - caarma_t_m*((B0_m*uackm + B3_m*vackm - B2_m*wackm)/gdt 
                               + uackm*B0_d_m - wackm*B2_d_m + vackm*B3_d_m))
                 - 2.0*amam2*(cam_m*((v_om*B3_om - w_om*B2_om)/gdt
                                      + w_m*B2_d_om - v_m*B3_d_om)
                 - caarma_t_m*((B0_om*uacokm + B3_om*vacokm - B2_om*wacokm)/gdt
                               + uackm*B0_d_om - wackm*B2_d_om + vackm*B3_d_om))
                 + 2.0*omak*(ftr_m*dtdB0(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m)
                             - (w0_m - Fb_m)*B0_m);
      s [1][ne + 7]  = -2.0*omam2*(cam*((v*B3 - w*B2)/gdt + w*B2_d - v*B3_d)
                 - caarma_t*((B0*uack + B3*vack - B2*wack)/gdt    
                              + uack*B0_d - wack*B2_d + vack*B3_d))
                 - 2.0*amam2*(cam*((v_o*B3_o - w_o*B2_o)/gdt 
                                   + w*B2_d_o - v*B3_d_o)
                 - caarma_t*((B0_o*uacok + B3_o*vacok - B2_o*wacok)/gdt
                             + uack*B0_d_o - wack*B2_d_o + vack*B3_d_o))
                  + 2.0*omak*(ftr*dtdB0(uack,vack,wack,B0,B1,B2,B3)
                              - (w0 - Fb)*B0_m);
      s [1][8]       = 2.0*omam2*(cam_m*((w_m*B3_m + v_m*B2_m)/gdt
                                         - v_m*B2_d_m - w_m*B3_d_m)
                + caarma_t_m*((B1_m*uackm + B2_m*vackm + B3_m*wackm)/gdt
			      + uackm*B1_d_m + vackm*B2_d_m + wackm*B3_d_m))
                + 2.0*amam2*(cam_m*((w_om*B3_om + v_om*B2_om)/gdt
                                    - v_m*B2_d_om - w_m*B3_d_om)
                + caarma_t_m*((B1_om*uacokm + B2_om*vacokm + B3_om*wacokm)/gdt
                              + uackm*B1_d_om + vackm*B2_d_om + wackm*B3_d_om))
                + 2.0*omak*(ftr_m*dtdB1(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m)
                            - (w0_m - Fb_m)*B1_m);
      s [1][ne + 8]  = 2.0*omam2*(cam*((w*B3 + v*B2)/gdt - v*B2_d - w*B3_d)
                + caarma_t*((B1*uack + B2*vack + B3*wack)/gdt  
                            + uack*B1_d + vack*B2_d + wack*B3_d))
                + 2.0*amam2*(cam*((w_o*B3_o + v_o*B2_o)/gdt
                                  - v*B2_d_o - w*B3_d_o)
                + caarma_t*((B1_o*uacok + B2_o*vacok + B3_o*wacok)/gdt
                              + uack*B1_d_o + vack*B2_d_o + wack*B3_d_o))
                + 2.0*omak*(ftr*dtdB1(uack,vack,wack,B0,B1,B2,B3)
                            - (w0 - Fb)*B1);
      s [1][9]       = -2.0*omam2*(cam_m*(-w_m*B0_d_m - v_m*B1_d_m
                                          + (w_m*B0_m + v_m*B1_m)/gdt)
             + caarma_t_m*(wackm*B0_d_m - vackm*B1_d_m + uackm*B2_d_m
                           + (B2_m*uackm - B1_m*vackm + B0_m*wackm)/gdt))
             - 2.0*amam2*(cam_m*(-w_m*B0_d_om - v_m*B1_d_om
                                 + (w_om*B0_om + v_om*B1_om)/gdt)
             + caarma_t_m*(wackm*B0_d_om - vackm*B1_d_om + uackm*B2_d_om
                           + (B2_om*uacokm - B1_om*vacokm + B0_om*wacokm)/gdt))
             + 2.0*omak*(ftr_m*dtdB2(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m)
                         + (w0_m - Fb_m)*B2_m);
      s [1][ne + 9]  = -2.0*omam2*(cam*(-w*B0_d - v*B1_d + (w*B0 + v*B1)/gdt)
             + caarma_t*(wack*B0_d - vack*B1_d + uack*B2_d
                           + (B2*uack - B1*vack + B0*wack)/gdt))
             - 2.0*amam2*(cam*(-w*B0_d_o - v*B1_d_o
                                 + (w_o*B0_o + v_o*B1_o)/gdt)
             + caarma_t*(wack*B0_d_o - vack*B1_d_o + uack*B2_d_o
                           + (B2_o*uacok - B1_o*vacok + B0_o*wacok)/gdt))
             + 2.0*omak*(ftr*dtdB2(uack,vack,wack,B0,B1,B2,B3)
                         + (w0 - Fb)*B2);
      s [1][10]      = -2.0*omam2*(cam_m*(v_m*B0_d_m - w_m*B1_d_m
                                          + (w_m*B1_m - v_m*B0_m)/gdt)
             - caarma_t_m*(vackm*B0_d_m + wackm*B1_d_m - uackm*B3_d_m
                           + (-B3_m*uackm + B0_m*vackm + B1_m*wackm)/gdt))
             - 2.0*amam2*(cam_m*(v_m*B0_d_om - w_m*B1_d_om
                                 + (w_om*B1_om - v_om*B0_om)/gdt)
             - caarma_t_m*(vackm*B0_d_om + wackm*B1_d_om - uackm*B3_d_om
                           + (-B3_om*uacokm + B0_om*vacokm + B1_om*wacokm)/gdt))
             + 2.0*omak*(ftr_m*dtdB3(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m)
                         + (w0_m - Fb_m)*B3_m);
      s [1][ne + 10] = -2.0*omam2*(cam*(v*B0_d - w*B1_d + (w*B1 - v*B0)/gdt)
             - caarma_t*(vack*B0_d + wack*B1_d - uack*B3_d
                         + (-B3*uack + B0*vack + B1*wack)/gdt))
             - 2.0*amam2*(cam*(v*B0_d_o - w*B1_d_o
                               + (w_o*B1_o - v_o*B0_o)/gdt)
             - caarma_t*(vack*B0_d_o + wack*B1_d_o - uack*B3_d_o
                         + (-B3_o*uacok + B0_o*vacok + B1_o*wacok)/gdt))
             + 2.0*omak*(ftr*dtdB3(uack,vack,wack,B0,B1,B2,B3) + (w0 - Fb)*B3);
      s [1][12]      = omak*(Sb_m - pay_m*w_m);
      s [1][ne + 12] = omak*(Sb - pay*w);
      s [1][13]      = omak*(pay_m*v_m - Sn_m);
      s [1][ne + 13] = omak*(pay*v - Sn); 
                                                         
      s [1][rhs] = 
       omak2_ds*((tdk + tdkm)*(e - e_m) - (pay_m + pay)*(u - u_m))
     + akak2_ds*((tdok + tdokm)*(e - e_m)
                 + (tdk + tdkm)*(e_o - e_om)
                 - (pay_om + pay_o)*(u - u_m)
                 - (pay_m + pay)*(u_o - u_om))
     + ak2_ds*((tdok + tdokm)*(e_o - e_om) - (pay_om + pay_o)*(u_o - u_om))
     - omam2*(camma_t*u_d + camma_t_m*u_d_m)
     - amam2*(camma_t*u_d + camma_t_m*u_d_m + camma_t*u_d_o + camma_t_m*u_d_om)
     - am2*(camma_t*u_d_o + camma_t_m*u_d_om)
     - 2.0*(cam*(v*B3 - w*B2) 
           - caarma_t*(B0*uack + B3*vack - B2*wack))*(omam2*B0_d + amam2*B0_d_o)
     + 2.0*(cam*(w*B3 + v*B2) 
           + caarma_t*(B1*uack + B2*vack + B3*wack))*(omam2*B1_d + amam2*B1_d_o)
     - 2.0*(cam*(w*B0 + v*B1)
          - caarma_t*(-B2*uack + B1*vack - B0*wack))*(omam2*B2_d + amam2*B2_d_o)
     - 2.0*(cam*(w*B1 - v*B0)
          - caarma_t*(-B3*uack + B0*vack + B1*wack))*(omam2*B3_d + amam2*B3_d_o)
     - 2.0*(cam_m*(v_m*B3_m - w_m*B2_m) 
            - caarma_t_m*(B0_m*uackm + B3_m*vackm 
                               - B2_m*wackm))*(omam2*B0_d_m + amam2*B0_d_om)
     + 2.0*(cam_m*(w_m*B3_m + v_m*B2_m) 
            + caarma_t_m*(B1_m*uackm + B2_m*vackm 
                               + B3_m*wackm))*(omam2*B1_d_m + amam2*B1_d_om)
     - 2.0*(cam_m*(w_m*B0_m + v_m*B1_m)
            - caarma_t_m*(-B2_m*uackm + B1_m*vackm 
                               - B0_m*wackm))*(omam2*B2_d_m + amam2*B2_d_om)
     - 2.0*(cam_m*(w_m*B1_m - v_m*B0_m)
            - caarma_t_m*(-B3_m*uackm + B0_m*vackm 
                               + B1_m*wackm))*(omam2*B3_d_m + amam2*B3_d_om)
     - 2.0*(cam*(v_o*B3_o - w_o*B2_o) 
            - caarma_t*(B0_o*uacok + B3_o*vacok 
                               - B2_o*wacok))*(amam2*B0_d + am2*B0_d_o)
     + 2.0*(cam*(w_o*B3_o + v_o*B2_o) 
            + caarma_t*(B1_o*uacok + B2_o*vacok 
                               + B3_o*wacok))*(amam2*B1_d + am2*B1_d_o)
     - 2.0*(cam*(w_o*B0_o + v_o*B1_o)
            - caarma_t*(-B2_o*uacok + B1_o*vacok 
                               - B0_o*wacok))*(amam2*B2_d + am2*B2_d_o)
     - 2.0*(cam*(w_o*B1_o - v_o*B0_o)
            - caarma_t*(-B3_o*uacok + B0_o*vacok 
                               + B1_o*wacok))*(amam2*B3_d + am2*B3_d_o)
     - 2.0*(cam_m*(v_om*B3_om - w_om*B2_om) 
            - caarma_t_m*(B0_om*uacokm + B3_om*vacokm 
                               - B2_om*wacokm))*(amam2*B0_d_m + am2*B0_d_om)
     + 2.0*(cam_m*(w_om*B3_om + v_om*B2_om) 
            + caarma_t_m*(B1_om*uacokm + B2_om*vacokm 
                               + B3_om*wacokm))*(amam2*B1_d_m + am2*B1_d_om)
     - 2.0*(cam_m*(w_om*B0_om + v_om*B1_om)
            - caarma_t_m*(-B2_om*uacokm + B1_om*vacokm 
                               - B0_om*wacokm))*(amam2*B2_d_m + am2*B2_d_om)
     - 2.0*(cam_m*(w_om*B1_om - v_om*B0_om)
            - caarma_t_m*(-B3_om*uacokm + B0_om*vacokm 
                               + B1_om*wacokm))*(amam2*B3_d_m + am2*B3_d_om)
     + omak*(Om2*Sb - Om3*Sn + Om2_m*Sb_m - Om3_m*Sn_m
             - (w0 - Fb)*(B0*B0 + B1*B1 - B2*B2 - B3*B3)
             - (w0_m - Fb_m)*(B0_m*B0_m + B1_m*B1_m - B2_m*B2_m - B3_m*B3_m)
             - (mud_b + damp_t)*u - (mud_bm + damp_t)*u_m
             - drat*(urk)*fabs(urk)*sqstrk 
             - drat_m*(urkm)*fabs(urkm)*sqstrkm
             - pay*(w*Om2 - v*Om3) - pay_m*(w_m*Om2_m - v_m*Om3_m))
     + ak*(Om2_o*Sb_o - Om3_o*Sn_o + Om2_om*Sb_om - Om3_om*Sn_om
           - (w0 - Fb_o)*(B0_o*B0_o + B1_o*B1_o - B2_o*B2_o - B3_o*B3_o)
           - (w0_m - Fb_om)*(B0_om*B0_om + B1_om*B1_om 
			     - B2_om*B2_om - B3_om*B3_om)
           - (mud_bo + damp_t)*u_o - (mud_bom + damp_t)*u_om
           - drat*(urok)*fabs(urok)*sqstrok 
           - drat_m*(urokm)*fabs(urokm)*sqstrokm
           - pay_o*(w_o*Om2_o - v_o*Om3_o) - pay_om*(w_om*Om2_om - v_om*Om3_om));
   
      ftrn_m = drap_m*sqstrkm*(pervkm + vrkm*vrkm/pervkm);
      ftrn   = drap*sqstrk*(pervk + vrk*vrk/pervk);
      ftrb_m = drap_m*sqstrkm*vrkm*wrkm/pervkm;
      ftrb   = drap*sqstrk*vrk*wrk/pervk;
 
      s [2][1]       = omak*(tdkm*Om3_m - 0.5*drap_m*vrkm*pervkm/sqstrkm);
      s [2][ne + 1]  = omak*(tdk*Om3 - 0.5*drap*vrk*pervk/sqstrk);
      s [2][2]       = -omak*t_ds;
      s [2][ne + 2]  = omak*t_ds;
      s [2][3]       = -omak*Om1_m;
      s [2][ne + 3]  = -omak*Om1;
      s [2][4]       = -omak*pay_m*Om3_m
                       -(omam2 + amam2)*2.0*cam_m*(-B3_m + B2_m - B1_m + B0_m);
      s [2][ne + 4]  = -omak*pay*Om3
                       -(omam2 + amam2)*2.0*cam*(-B3 + B2 - B1 + B0);
      s [2][5]       = -camma_n_m/gdt*(omam2 + amam2)
                       - omak*(mud_bm + damp_n + ftrn_m)
                       + omak2_ds*(pay + pay_m)
                       + akak2_ds*(pay_om + pay_o);
      s [2][ne + 5]   = -camma_n/gdt*(omam2 + amam2)
                        - omak*(mud_b + damp_n + ftrn)
                       - omak2_ds*(pay + pay_m)
                       - akak2_ds*(pay_om + pay_o);
      s [2][6]       = omak*pay_m*Om1_m
                       -(omam2 + amam2)*2.0*cam_m*(B1_m - B0_m - B3_m + B2_m)
                       - omak*ftrb_m;
      s [2][ne + 6]  = omak*pay*Om1
                       -(omam2 + amam2)*2.0*cam*(B1 - B0 - B3 + B2)
                       - omak*ftrb;
      s [2][7]       = -omam2*2.0*(cam_m*((-w_m*B1_d_m + u_m*B3_d_m) 
                                          + (w_m*B1_m - u_m*B3_m)/gdt)
                                   + caarma_n_m*((B3_m*uackm - B0_m*vackm 
                                                      - B1_m*wackm)/gdt 
                                               - vackm*B0_d_m - wackm*B1_d_m
                                               + uackm*B3_d_m))
                       - amam2*2.0*(cam_m*((-w_m*B1_d_om + u_m*B3_d_om)
                                           + (w_om*B1_om - u_om*B3_om)/gdt)
                                    + caarma_n_m*((B3_om*uacokm - B1_om*vacokm 
                                                         - B0_om*wacokm)/gdt
                                                - vackm*B0_d_om - wackm*B1_d_om
                                                + uackm*B3_d_om))
                       + 2.0*omak*(w0_m - Fb_m)*B3_m
                       + omak*(dndB0(uackm,vackm,wackm,
                                     B0_m,B1_m,B2_m,B3_m)*ftrn_m
                               + dbdB0(uackm,vackm,wackm,
                                       B0_m,B1_m,B2_m,B3_m)*ftrb_m);
      s [2][ne + 7]  = -omam2*2.0*(cam*((-w*B1_d + u*B3_d) + (w*B1 - u*B3)/gdt)
                                   + caarma_n*((B3*uack - B0*vack - B1*wack)/gdt  
                                             - vack*B0_d - wack*B1_d
                                             + uack*B3_d))
                       - amam2*2.0*(cam*((-w*B1_d_o + u*B3_d_o)
                                         + (w_o*B1_o - u_o*B3_o)/gdt)
                                    + caarma_n*((B3_o*uacok - B1_o*vacok
                                                             - B0_o*wacok)/gdt
                                              - vack*B0_d_o - wack*B1_d_o
                                              + uack*B3_d_o))
                       + 2.0*omak*(w0 - Fb)*B3 
                       + omak*(dndB0(uack,vack,wack,B0,B1,B2,B3)*ftrn
                               + dbdB0(uack,vack,wack, B0,B1,B2,B3)*ftrb);
      s [2][8]       = -omam2*2.0*(cam_m*((w_m*B0_d_m - u_m*B2_d_m)
                                          + (u_m*B2_m - w_m*B0_m)/gdt)
                                   + caarma_n_m*(-(B2_m*uackm - B1_m*vackm 
                                                              + B0_m*wackm)/gdt
                                               - wackm*B0_d_m + vackm*B1_d_m 
                                               - uackm*B2_d_m))
                       - amam2*2.0*(cam_m*((w_m*B0_d_om - u_m*B2_d_om)
                                           + (u_om*B2_om - w_om*B0_om)/gdt)
                                    + caarma_n_m*(-(B2_om*uacokm - B1_om*vacokm
                                                           + B0_om*wacokm)/gdt
                                                - wackm*B0_d_om + vackm*B1_d_om
                                                - uackm*B2_d_om))
                       - omak*2.0*(w0_m - Fb_m)*B2_m
                       + omak*(dndB1(uackm,vackm,wackm,
                                     B0_m,B1_m,B2_m,B3_m)*ftrn_m
                               + dbdB1(uackm,vackm,wackm,
                                       B0_m,B1_m,B2_m,B3_m)*ftrb_m);
      s [2][ne + 8]  = -omam2*2.0*(cam*((w*B0_d - u*B2_d) + (u*B2 - w*B0)/gdt)
                                   + caarma_n*(-(B2*uack - B1*vack + B0*wack)/gdt
                                             - wack*B0_d + vack*B1_d 
                                             - uack*B2_d))
                       - amam2*2.0*(cam*((w*B0_d_o - u*B2_d_o)
                                         + (u_o*B2_o - w_o*B0_o)/gdt)
                                    + caarma_n*(-(B2_o*uacok - B1_o*vacok
                                                           + B0_o*wacok)/gdt
                                              - wack*B0_d_o + vack*B1_d_o
                                              - uack*B2_d_o))
                       - omak*2.0*(w0 - Fb)*B2
                       + omak*(dndB1(uack,vack,wack,B0,B1,B2,B3)*ftrn
                               + dbdB1(uack,vack,wack, B0,B1,B2,B3)*ftrb);
      s [2][9]       = -omam2*2.0*(cam_m*((u_m*B1_d_m + w_m*B3_d_m)
                                          - (u_m*B1_m + w_m*B3_m)/gdt)
                                   + caarma_n_m*(-(B1_m*uackm + B2_m*vackm 
                                                              + B3_m*wackm)/gdt
                                               - uackm*B1_d_m - vackm*B2_d_m
                                               - wackm*B3_d_m)) 
                       - amam2*2.0*(cam_m*((u_m*B1_d_om + w_m*B3_d_om)
                                           - (u_om*B1_om + w_om*B3_om)/gdt)
                                    + caarma_n_m*(-(B1_om*uacokm + B2_om*vacokm
                                                             + B3_om*wacokm)/gdt
                                                - uackm*B1_d_om - vackm*B2_d_om
                                                - wackm*B3_d_om))
                       - omak*2.0*(w0_m - Fb_m)*B1_m
                       + omak*(dndB2(uackm,vackm,wackm,
                                     B0_m,B1_m,B2_m,B3_m)*ftrn_m
                               + dbdB2(uackm,vackm,wackm,
                                       B0_m,B1_m,B2_m,B3_m)*ftrb_m);
      s [2][ne + 9]  = -omam2*2.0*(cam*((u*B1_d + w*B3_d) - (u*B1 + w*B3)/gdt)
                                   + caarma_n*(-(B1*uack + B2*vack + B3*wack)/gdt
                                               - uack*B1_d - vack*B2_d
                                               - wack*B3_d))
                       - amam2*2.0*(cam*((u*B1_d_o + w*B3_d_o)
                                           - (u_o*B1_o + w_o*B3_o)/gdt)
                                    + caarma_n*(-(B1_o*uacok + B2_o*vacok
                                                             + B3_o*wacok)/gdt
                                              - uack*B1_d_o - vack*B2_d_o
                                              - wack*B3_d_o))
                       - omak*2.0*(w0 - Fb)*B1
                       + omak*(dndB2(uack,vack,wack,B0,B1,B2,B3)*ftrn
                               + dbdB2(uack,vack,wack,B0,B1,B2,B3)*ftrb);
      s [2][10]      = -omam2*2.0*(cam_m*((-u_m*B0_d_m - w_m*B2_d_m)
                                          + (u_m*B0_m + w_m*B2_m)/gdt)
                                   + caarma_n_m*((B0_m*uackm + B3_m*vackm 
                                                            - B2_m*wackm)/gdt
                                               + uackm*B0_d_m - wackm*B2_d_m
                                               + vackm*B3_d_m))
                       - amam2*2.0*(cam_m*((-u_m*B0_d_om - w_m*B2_d_om)
                                           + (u_om*B0_om + w_om*B2_om)/gdt)
                                    + caarma_n_m*((B0_om*uacokm + B3_om*vacokm
                                                             - B2_om*wacokm)/gdt
                                              + uackm*B0_d_om - wackm*B2_d_om
                                              + vackm*B3_d_om))
                       + omak*2.0*(w0_m - Fb_m)*B0_m
                       + omak*(dndB3(uackm,vackm,wackm,
                                     B0_m,B1_m,B2_m,B3_m)*ftrn_m
                               + dbdB3(uackm,vackm,wackm,
                                       B0_m,B1_m,B2_m,B3_m)*ftrb_m);
      s [2][ne + 10] = -omam2*2.0*(cam*((-u*B0_d - w*B2_d)
                                        + (u*B0 + w*B2_m)/gdt)
                                   + caarma_n*((B0*uack + B3*vack - B2*wack)/gdt
                                               + uack*B0_d - wack*B2_d
                                               + vack*B3_d))
                       - amam2*2.0*(cam*((-u*B0_d_o - w*B2_d_o)
                                         + (u_o*B0_o + w_o*B2_o)/gdt)
                                    + caarma_n*(-(B0_o*uacok + B3_o*vacok
                                                           - B2_o*wacok)/gdt
                                              + uack*B0_d_o - wack*B2_d_o
                                              + vack*B3_d_o))
                       + omak*2.0*(w0 - Fb)*B0
                       + omak*(dndB3(uack,vack,wack,B0,B1,B2,B3)*ftrn
                               + dbdB3(uack,vack,wack, B0,B1,B2,B3)*ftrb);
      s [2][11]      = omak*(w_m*pay_m - Sb_m);
      s [2][ne + 11] = omak*(w*pay - Sb);
      s [2][13]      = omak*(tkm - u_m*pay_m);
      s [2][ne + 13] = omak*(tk - u*pay);
      s [2][rhs] =
       omak*t_ds*(Sn - Sn_m) + ak*t_ds*(Sn_o - Sn_om)
     - omak2_ds*(pay_m + pay)*(v - v_m)
     - akak2_ds*((pay_om + pay_o)*(v - v_m)
                 + (pay_m + pay)*(v_o - v_om))
     - ak2_ds*(pay_om + pay_o)*(v_o - v_om)
     - omam2*(camma_n*v_d + camma_n_m*v_d_m)
     - amam2*(camma_n*v_d + camma_n_m*v_d_m + camma_n*v_d_o + camma_n_m*v_d_om)
     - am2*(camma_n*v_d_o + camma_n_m*v_d_om)
     - 2.0*(cam*(w*B1 - u*B3)
           + caarma_n*(B3*uack - B0*vack - B1*wack))*(omam2*B0_d + amam2*B0_d_o)
     - 2.0*(cam_m*(w_m*B1_m - u_m*B3_m)
            + caarma_n_m*(B3_m*uackm - B0_m*vackm 
                               - B1_m*wackm))*(omam2*B0_d_m + amam2*B0_d_om)
     - 2.0*(cam*(u*B2 - w*B0)
           - caarma_n*(B2*uack - B1*vack + B0*wack))*(omam2*B1_d + amam2*B1_d_o)
     - 2.0*(cam_m*(u_m*B2_m - w_m*B0_m)
            - caarma_n_m*(B2_m*uackm - B1_m*vackm 
                               + B0_m*wackm))*(omam2*B1_d_m + amam2*B1_d_om)
     + 2.0*(cam*(u*B1 + w*B3)
           + caarma_n*(B1*uack + B2*vack + B3*wack))*(omam2*B2_d + amam2*B2_d_o)
     + 2.0*(cam_m*(u_m*B1_m + w_m*B3_m)
            + caarma_n_m*(B1_m*uackm + B2_m*vackm 
                               + B3_m*wackm))*(omam2*B2_d_m + amam2*B2_d_om)
     - 2.0*(cam*(u*B0 + w*B2)
           + caarma_n*(B0*uack + B3*vack - B2*wack))*(omam2*B3_d + amam2*B3_d_o)
     - 2.0*(cam_m*(u_m*B0_m + w_m*B2_m)
            + caarma_n_m*(B0_m*uackm + B3_m*vackm 
                               - B2_m*wackm))*(omam2*B3_d_m + amam2*B3_d_om)
     - 2.0*(cam*(w_o*B1_o - u_o*B3_o)
            + caarma_n*(B3_o*uacok - B0_o*vacok 
                               - B1_o*wacok))*(amam2*B0_d + am2*B0_d_o)
     - 2.0*(cam_m*(w_om*B1_om - u_om*B3_om)
            + caarma_n_m*(B3_om*uacokm - B0_om*vacokm 
                               - B1_om*wacokm))*(amam2*B0_d_m + am2*B0_d_om)
     - 2.0*(cam*(u_o*B2_o - w_o*B0_o)
            - caarma_n*(B2_o*uacok - B1_o*vacok 
                               + B0_o*wacok))*(amam2*B1_d + am2*B1_d_o)
     - 2.0*(cam_m*(u_om*B2_om - w_om*B0_om)
               - caarma_n_m*(B2_om*uacokm - B1_om*vacokm 
                               + B0_om*wacokm))*(amam2*B1_d_m + am2*B1_d_om)
     + 2.0*(cam*(u_o*B1_o + w_o*B3_o)
               + caarma_n*(B1_o*uacok + B2_o*vacok 
                               + B3_o*wacok))*(amam2*B2_d + am2*B2_d_o)
     + 2.0*(cam_m*(u_om*B1_om + w_om*B3_om)
               + caarma_n_m*(B1_om*uacokm + B2_om*vacokm 
                               + B3_om*wacokm))*(amam2*B2_d_m + am2*B2_d_om)
     - 2.0*(cam*(u_o*B0_o + w_o*B2_o)
               + caarma_n*(B0_o*uacok + B3_o*vacok 
                               - B2_o*wacok))*(amam2*B3_d + am2*B3_d_o)
     - 2.0*(cam_m*(u_om*B0_om + w_om*B2_om)
               + caarma_n_m*(B0_om*uacokm + B3_om*vacokm 
                               - B2_om*wacokm))*(amam2*B3_d_m + am2*B3_d_om)
      + omak*(tk*Om3 + tkm*Om3_m - Sb*Om1 - Sb_m*Om1_m
              - 2.0*(w0 - Fb)*(B1*B2 - B0*B3) 
	      - 2.0*(w0_m - Fb_m)*(B1_m*B2_m - B0_m*B3_m)
              - (mud_b + damp_n)*v - (mud_bm + damp_n)*v_m
              - drap*vrk*pervk*sqstrk - drap_m*vrkm*pervkm*sqstrkm
              - pay*(u*Om3 - w*Om1) - pay_m*(u_m*Om3_m - w_m*Om1_m))
      + ak*(tok*Om3_o + tokm*Om3_om - Sb_o*Om1_o - Sb_om*Om1_om
            - 2.0*(w0 - Fb_o)*(B1_o*B2_o - B0_o*B3_o) 
	    - 2.0*(w0_m - Fb_om)*(B1_om*B2_om - B0_om*B3_om)
            - (mud_bo + damp_n)*v_o - (mud_bom + damp_n)*v_om
            - drap_o*vrok*pervok*sqstrok - drap_om*vrokm*pervokm*sqstrokm
            - pay_o*(u_o*Om3_o - w_o*Om1_o) - pay_om*(u_om*Om3_om - w_om*Om1_om));

      ftrb_m = drap_m*sqstrkm*(pervkm + wrkm*wrkm/pervkm);
      ftrb   = drap*sqstrk*(pervk + wrk*wrk/pervk);
      ftrn_m = drap_m*sqstrkm*vrkm*wrkm/pervkm;
      ftrn   = drap*sqstrk*vrk*wrk/pervk;

      s [3][1]       = -omak*(tdkm*Om2_m + 0.5*drap_m*wrkm*pervkm/sqstrkm);
      s [3][ne + 1]  = -omak*(tdk*Om2 + 0.5*drap*wrk*pervk/sqstrk);
      s [3][2]       = omak*Om1_m;
      s [3][ne + 2]  = omak*Om1;
      s [3][3]       = -omak*t_ds;
      s [3][ne + 3]  = omak*t_ds;
      s [3][4]       = omak*pay_m*Om2_m
                       -(omam2 + amam2)*2.0*cam_m*(B2_m + B3_m - B0_m - B1_m);
      s [3][ne + 4]  = omak*pay*Om2
                       -(omam2 + amam2)*2.0*cam*(B2 + B3 - B0 - B1);
      s [3][5]       = -omak*pay_m*Om1_m
                       -(omam2 + amam2)*2.0*cam_m*(-B1_m + B0_m + B3_m - B2_m)
                       - omak*ftrn_m;
      s [3][ne + 5]  = -omak*pay*Om1
                       -(omam2 + amam2)*2.0*cam*(-B1 + B0 + B3 - B2)
                       - omak*ftrn;
      s [3][6]       = -(omam2 + amam2)*camma_n_m/gdt
                       - omak*(mud_bm + damp_n + ftrb_m)
                       + omak2_ds*(pay + pay_m)
                       + akak2_ds*(pay_om + pay_o);
      s [3][ne + 6]  = -(omam2 + amam2)*camma_n/gdt
                       - omak*(mud_b + damp_n + ftrb)
                       - omak2_ds*(pay + pay_m)
                       - akak2_ds*(pay_om + pay_o);
      s [3][7]       = -omam2*2.0*(cam_m*(-(v_m*B1_m - u_m*B2_m)/gdt
                                          + v_m*B1_d_m - u_m*B2_d_m)
                                   + caarma_n_m*(-(B2_m*uackm - B1_m*vackm 
                                                            + B0_m*wackm)/gdt
                                               - wackm*B0_d_m + vackm*B1_d_m
                                               - uackm*B2_d_m))
                       - amam2*2.0*(cam_m*(-(v_om*B1_om - u_om*B2_om)/gdt
                                           + v_m*B1_d_om - u_m*B2_d_om)
                                    + caarma_n_m*(-(B2_om*uacokm - B1_om*vacokm
                                                            + B0_om*wacokm)/gdt
                                                - wackm*B0_d_om + vackm*B1_d_om
                                                - uackm*B2_d_om))
                       - omak*2.0*(w0_m - Fb_m)*B2_m
                       + omak*(dbdB0(uackm,vackm,wackm,
                                     B0_m,B1_m,B2_m,B3_m)*ftrb_m
                               + dndB0(uackm,vackm,wackm,
                                       B0_m,B1_m,B2_m,B3_m)*ftrn_m);
     s [3][ne + 7]   = -omam2*2.0*(cam*(-(v*B1 - u*B2)/gdt + v*B1_d - u*B2_d)
                                   + caarma_n*(-(B2*uack - B1*vack + B0*wack)/gdt
                                               - wack*B0_d + vackm*B1_d
                                               - uack*B2_d))
                       - amam2*2.0*(cam*(-(v_o*B1_o - u_o*B2_o)/gdt
                                         + v*B1_d_o - u*B2_d_o)
                                    + caarma_n*(-(B2_o*uacok - B1_o*vacok
                                                           + B0_o*wacok)/gdt
                                                - wack*B0_d_o + vack*B1_d_o
                                                - uack*B2_d_o))
                       - omak*2.0*(w0 - Fb)*B2
                       + omak*(dbdB0(uack,vack,wack,B0,B1,B2,B3)*ftrb
                               + dndB0(uack,vack,wack, B0,B1,B2,B3)*ftrn);
      s [3][8]       = -omam2*2.0*(cam_m*((-v_m*B0_d_m - u_m*B3_d_m)
                                          + (v_m*B0_m + u_m*B3_m)/gdt)
                                   + caarma_n_m*(-(B3_m*uackm - B0_m*vackm 
                                                            - B1_m*wackm)/gdt
                                               + vackm*B0_d_m + wackm*B1_d_m
                                               - uackm*B3_d_m))
                       - amam2*2.0*(cam_m*((-v_m*B0_d_om - u_m*B3_d_om)
                                           + (v_om*B0_om + u_om*B3_om)/gdt)
                                    + caarma_n_m*(-(B3_om*uacokm - B0_om*vacokm
                                                           - B1_om*wacokm)/gdt  
                                                + vackm*B0_d_om + wackm*B1_d_om
                                                - uackm*B3_d_om))
                       - omak*2.0*(w0_m - Fb_m)*B3_m
                       + omak*(dbdB1(uackm,vackm,wackm,
                                     B0_m,B1_m,B2_m,B3_m)*ftrb_m
                               + dndB1(uackm,vackm,wackm,
                                       B0_m,B1_m,B2_m,B3_m)*ftrn_m);
      s [3][ne + 8]  = -omam2*2.0*(cam*((-v*B0_d - u*B3_d) + (v*B0 + u*B3)/gdt)
                                   + caarma_n*(-(B3*uack - B0*vack - B1*wack)/gdt
                                             + vack*B0_d + wack*B1_d
                                             - uack*B3_d))
                       - amam2*2.0*(cam*((-v*B0_d_o - u*B3_d_o)
                                         + (v_o*B0_o + u_o*B3_o)/gdt)
                                    + caarma_n*(-(B3_o*uacok - B0_o*vacok
                                                           - B1_o*wacok)/gdt
                                              + vack*B0_d_o + wack*B1_d_o
                                              - uack*B3_d_o))
                       - omak*2.0*(w0 - Fb)*B3
                       + omak*(dbdB1(uack,vack,wack,B0,B1,B2,B3)*ftrb
                               + dndB1(uack,vack,wack, B0,B1,B2,B3)*ftrn);
      s [3][9]       = -omam2*2.0*(cam_m*((u_m*B0_d_m - v_m*B3_d_m)
                                          - (u_m*B0_m - v_m*B3_m)/gdt)
                                   + caarma_n_m*(-(B0_m*uackm + B3_m*vackm 
                                                            - B2_m*wackm)/gdt
                                               - uackm*B0_d_m + wackm*B2_d_m
                                               - vackm*B3_d_m))
                       - amam2*2.0*(cam_m*((u_m*B0_d_om - v_m*B3_d_om)
                                           - (u_om*B0_om - v_om*B3_om)/gdt)
                                    + caarma_n_m*(-(B0_om*uacokm + B3_om*vacokm
                                                          - B2_om*wacokm)/gdt
                                              - uackm*B0_d_om + wackm*B2_d_om
                                              - vackm*B3_d_om))
                       - omak*2.0*(w0_m - Fb_m)*B0_m
                       + omak*(dbdB2(uackm,vackm,wackm,
                                     B0_m,B1_m,B2_m,B3_m)*ftrb_m
                               + dndB2(uackm,vackm,wackm,
                                       B0_m,B1_m,B2_m,B3_m)*ftrn_m);
      s [3][ne + 9]  = -omam2*2.0*(cam*((u*B0_d - v*B3_d) - (u*B0 - v*B3)/gdt)
                                   + caarma_n*(-(B0*uack + B3*vack - B2*wack)/gdt
                                             - uack*B0_d + wack*B2_d
                                             - vack*B3_d))
                       - amam2*2.0*(cam*((u*B0_d_o - v*B3_d_o)
                                         - (u_o*B0_o - v_o*B3_o)/gdt)
                                    + caarma_n*(-(B0_o*uacok + B3_o*vacok
                                                           - B2_o*wacok)/gdt
                                              - uack*B0_d_o + wack*B2_d_o
                                              - vack*B3_d_o))
                       - omak*2.0*(w0 - Fb)*B0
                       + omak*(dbdB2(uack,vack,wack,B0,B1,B2,B3)*ftrb
                               + dndB2(uack,vack,wack, B0,B1,B2,B3)*ftrn);
      s [3][10]      = -omam2*2.0*(cam_m*((u_m*B1_d_m + v_m*B2_d_m)
                                          - (v_m*B2_m + u_m*B1_m)/gdt)
                                   + caarma_n_m*(-(B1_m*uackm + B2_m*vackm 
                                                            + B3_m*wackm)/gdt 
                                               - uackm*B1_d_m - vackm*B2_d_m
                                               - wackm*B3_d_m))
                       + amam2*2.0*(cam_m*((u_m*B1_d_om + v_m*B2_d_om)
                                           - (v_om*B2_om + u_om*B1_om)/gdt)
                                    + caarma_n_m*(-(B1_om*uacokm + B2_om*vacokm
                                                          + B3_om*wacokm)/gdt
                                              - uackm*B1_d_om - vackm*B2_d_om
                                              - wackm*B3_d_om))
                       - omak*2.0*(w0_m - Fb_m)*B1_m
                       + omak*(dbdB3(uackm,vackm,wackm,
                                     B0_m,B1_m,B2_m,B3_m)*ftrb_m
                               + dndB3(uackm,vackm,wackm,
                                       B0_m,B1_m,B2_m,B3_m)*ftrn_m);
      s [3][ne + 10] = -omam2*2.0*(cam*((u*B1_d + v*B2_d) - (v*B2 + u*B1)/gdt)
                                   + caarma_n*(-(B1*uack + B2*vack + B3*wack)/gdt 
                                             - uack*B1_d - vack*B2_d
                                             - wack*B3_d))
                       + amam2*2.0*(cam*((u*B1_d_o + v*B2_d_o)
                                           - (v_o*B2_o + u_o*B1_o)/gdt)
                                    + caarma_n*(-(B1_o*uacok + B2_o*vacok
                                                           + B3_o*wacok)/gdt
                                              - uack*B1_d_o - vack*B2_d_o
                                              - wack*B3_d_o))
                       - omak*2.0*(w0 - Fb)*B1
                       + omak*(dbdB3(uack,vack,wack,B0,B1,B2,B3)*ftrb
                               + dndB3(uack,vack,wack, B0,B1,B2,B3)*ftrn);
      s [3][11]      = omak*(Sn_m - pay_m*v_m);
      s [3][ne + 11] = omak*(Sn - pay*v);
      s [3][12]      = omak*(u_m*pay_m - tkm);
      s [3][ne + 12] = omak*(u*pay - tk);
      s [3][rhs] = 
       omak*t_ds*(Sb - Sb_m) + ak*t_ds*(Sb_o - Sb_om)
     - omak2_ds*(pay_m + pay)*(w - w_m)
     - akak2_ds*((pay_om + pay_o)*(w - w_m)
                 + (pay_m + pay)*(w_o - w_om))
     - ak2_ds*(pay_om + pay_o)*(w_o - w_om)
     - omam2*(camma_n*w_d + camma_n_m*w_d_m)
     - amam2*(camma_n*w_d + camma_n_m*w_d_m + camma_n*w_d_o + camma_n_m*w_d_om)
     - am2*(camma_n*w_d_o + camma_n_m*w_d_om)
     + 2.0*(cam*(v*B1 - u*B2)
           + caarma_n*(B2*uack - B1*vack + B0*wack))*(omam2*B0_d + amam2*B0_d_o)
     + 2.0*(cam_m*(v_m*B1_m - u_m*B2_m)
            + caarma_n_m*(B2_m*uackm - B1_m*vackm 
                               + B0_m*wackm))*(omam2*B0_d_m + amam2*B0_d_om)
     - 2.0*(cam*(v*B0 + u*B3)
           - caarma_n*(B3*uack - B0*vack - B1*wack))*(omam2*B1_d + amam2*B1_d_o)
     - 2.0*(cam_m*(v_m*B0_m + u_m*B3_m)
            - caarma_n_m*(B3_m*uackm - B0_m*vackm 
                               - B1_m*wackm))*(omam2*B1_d_m + amam2*B1_d_om)
     + 2.0*(cam*(u*B0 - v*B3)
           + caarma_n*(B0*uack + B3*vack - B2*wack))*(omam2*B2_d + amam2*B2_d_o)
     + 2.0*(cam_m*(u_m*B0_m - v_m*B3_m)
            + caarma_n_m*(B0_m*uackm + B3_m*vackm 
                               - B2_m*wackm))*(omam2*B2_d_m + amam2*B2_d_om)
     + 2.0*(cam*(v*B2 + u*B1)
           + caarma_n*(B1*uack + B2*vack + B3*wack))*(omam2*B3_d + amam2*B3_d_o)
     + 2.0*(cam_m*(v_m*B2_m + u_m*B1_m)
            + caarma_n_m*(B1_m*uackm + B2_m*vackm 
                               + B3_m*wackm))*(omam2*B3_d_m + amam2*B3_d_om)    
     + 2.0*(cam*(v_o*B1_o - u_o*B2_o)
            + caarma_n*(B2_o*uacok - B1_o*vacok 
                               + B0_o*wacok))*(amam2*B0_d + am2*B0_d_o)
     + 2.0*(cam_m*(v_om*B1_om - u_om*B2_om)
            + caarma_n_m*(B2_om*uacokm - B1_om*vacokm 
                               + B0_om*wacokm))*(amam2*B0_d_m + am2*B0_d_om)
     - 2.0*(cam*(v_o*B0_o + u_o*B3_o)
            - caarma_n*(B3_o*uacok - B0_o*vacok 
                               - B1_o*wacok))*(amam2*B1_d + am2*B1_d_o)
     - 2.0*(cam_m*(v_om*B0_om + u_om*B3_om)
            - caarma_n_m*(B3_om*uacokm - B0_om*vacokm 
                               - B1_om*wacokm))*(amam2*B1_d_m + am2*B1_d_om)
     + 2.0*(cam*(u_o*B0_o - v_o*B3_o)
            + caarma_n*(B0_o*uacok + B3_o*vacok 
                               - B2_o*wacok))*(amam2*B2_d + am2*B2_d_o)
     + 2.0*(cam_m*(u_om*B0_om - v_om*B3_om)
            + caarma_n_m*(B0_om*uacokm + B3_om*vacokm 
                               - B2_om*wacokm))*(amam2*B2_d_m + am2*B2_d_om)
     + 2.0*(cam*(v_o*B2_o + u_o*B1_o)
            + caarma_n*(B1_o*uacok + B2_o*vacok 
                               + B3_o*wacok))*(amam2*B3_d + am2*B3_d_o) 
     + 2.0*(cam_m*(v_om*B2_om + u_om*B1_om)
            + caarma_n_m*(B1_om*uacokm + B2_om*vacokm 
                               + B3_om*wacokm))*(amam2*B3_d_m + am2*B3_d_om)
     + omak*(Sn*Om1 + Sn_m*Om1_m - tk*Om2 - tkm*Om2_m
             - 2.0*(w0 - Fb)*(B1*B3 + B0*B2)
             - 2.0*(w0_m - Fb_m)*(B1_m*B3_m + B0_m*B2_m)
             - (damp_n + mud_b)*w - (damp_n + mud_bm)*w_m
             - drap*pervk*wrk*sqstrk - drap_m*pervkm*wrkm*sqstrkm
             - pay*(v*Om1 - u*Om2) - pay_m*(v_m*Om1_m - u_m*Om2_m))
     + ak*(Sn_o*Om1_o + Sn_om*Om1_om - tok*Om2_o - tokm*Om2_om
           - 2.0*(w0 - Fb_o)*(B1_o*B3_o + B0_o*B2_o)
           - 2.0*(w0_m - Fb_om)*(B1_om*B3_om + B0_om*B2_om)
           - (damp_n + mud_bo)*w_o - (damp_n + mud_bom)*w_om
           - drap_o*pervok*wrok*sqstrok - drap_om*pervokm*wrokm*sqstrokm
           - pay_o*(v_o*Om1_o - u_o*Om2_o) - pay_om*(v_om*Om1_om - u_om*Om2_om));

      if (nm == n -> segment -> first_active) {
          pay_m = n -> segment -> bottom_pay_speed;
          pay_om = n -> segment -> bottom_pay_speed_o;
      }
      else if (n == n -> segment -> last_active) {
          pay = nm -> segment -> top_pay_speed;
          pay_o = nm -> segment -> top_pay_speed_o;
      }

      s [4][1]       = -omam/gdt;
      s [4][ne + 1]  = -omam/gdt;
      s [4][4]       = -omak*t_ds;
      s [4][ne + 4]  = omak*t_ds;
      s [4][5]       = -omak*Om3_m;
      s [4][ne + 5]  = -omak*Om3;
      s [4][6]       = omak*Om2_m;
      s [4][ne + 6]  = omak*Om2;
      s [4][12]      = omak*w_m;
      s [4][ne + 12] = omak*w;
      s [4][13]      = -omak*v_m;
      s [4][ne + 13] = -omak*v;
      s [4][rhs] = 
       omak*t_ds*((u - u_m) + (pay - pay_m)) 
     + ak*t_ds*((u_o - u_om) + (pay_o - pay_om))
     - omam*(e_d + e_d_m) - am*(e_d_o + e_d_om)
     + omak*(Om2*w + Om2_m*w_m - Om3*v - Om3_m*v_m)
     + ak*(Om2_o*w_o + Om2_om*w_om - Om3_o*v_o - Om3_om*v_om);

      s [5][1]       = 2.0*omam2*(B3_m*B0_d_m - B2_m*B1_d_m
                                  + B1_m*B2_d_m - B0_m*B3_d_m)
                       + 2.0*amam2*(B3_m*B0_d_om - B2_m*B1_d_om
                                    + B1_m*B2_d_om - B0_m*B3_d_om);
      s [5][ne + 1]  = 2.0*omam2*(B3*B0_d - B2*B1_d  
                                  + B1*B2_d - B0*B3_d)
                       + 2.0*amam2*(B3*B0_d_o - B2*B1_d_o 
                                    + B1*B2_d_o - B0*B3_d_o);
      s [5][4]       = omak*Om3_m;
      s [5][ne + 4]  = omak*Om3;
      s [5][5]       = -omak*t_ds;
      s [5][ne + 5]  = omak*t_ds;
      s [5][6]       = -omak*Om1_m;
      s [5][ne + 6]  = -omak*Om1;
      s [5][7]       = 2.0*omam2*ep_m*(B3_m/gdt - B3_d_m)
                       + 2.0*amam2*(ep_om*B3_om/gdt - ep_m*B3_d_om);
      s [5][ne + 7]  = 2.0*omam2*ep*(B3/gdt - B3_d)
                       + 2.0*amam2*(ep_o*B3_o/gdt - ep*B3_d_o);
      s [5][8]       = 2.0*omam2*ep_m*(B2_d_m - B2_m/gdt)
                       + 2.0*amam2*(ep_m*B2_d_om - ep_om*B2_om/gdt);
      s [5][ne + 8]  = 2.0*omam2*ep*(B2_d - B2/gdt)
                       + 2.0*amam2*(ep*B2_d_o - ep_o*B2_o/gdt);
      s [5][9]       = 2.0*omam2*ep_m*(B1_m/gdt - B1_d_m)
                       + 2.0*amam2*(ep_om*B1_om/gdt - ep_m*B1_d_om);
      s [5][ne + 9]  = 2.0*omam2*ep*(B1/gdt - B1_d)
                       + 2.0*amam2*(ep_o*B1_o/gdt - ep*B1_d_o);
      s [5][10]      = 2.0*omam2*ep_m*(B0_d_m - B0_m/gdt)
                       + 2.0*amam2*(ep_m*B0_d_om - ep_om*B0_om/gdt);
      s [5][ne + 10] = 2.0*omam2*ep*(B0_d - B0/gdt)
                       + 2.0*amam2*(ep*B0_d_o - ep_o*B0_o/gdt);
      s [5][11]      = -omak*w_m;
      s [5][ne + 11] = -omak*w;
      s [5][13]      = omak*(u_m + pay_m);
      s [5][ne + 13] = omak*(u + pay);
      s [5][rhs] = 
       omak*t_ds*(v - v_m) + ak*t_ds*(v_o - v_om)
     + 2.0*omam2*(ep*(B3*B0_d - B2*B1_d 
                      + B1*B2_d - B0*B3_d)
                  + ep_m*(B3_m*B0_d_m - B2_m*B1_d_m 
                          + B1_m*B2_d_m - B0_m*B3_d_m))
     + 2.0*amam2*(ep_o*(B3_o*B0_d - B2_o*B1_d
                        + B1_o*B2_d - B0_o*B3_d)
                  + ep_om*(B3_om*B0_d_m - B2_om*B1_d_m
                           + B1_om*B2_d_m - B0_om*B3_d_m)
                  + ep*(B3*B0_d_o - B2*B1_d_o 
                        + B1*B2_d_o - B0*B3_d_o)
                  + ep_m*(B3_m*B0_d_om - B2_m*B1_d_om
                          + B1_m*B2_d_om - B0_m*B3_d_om))
     + 2.0*am2*(ep_o*(B3_o*B0_d_o - B2_o*B1_d_o 
                      + B1_o*B2_d_o - B0_o*B3_d_o)
                + ep_om*(B3_om*B0_d_om - B2_om*B1_d_om
                         + B1_om*B2_d_om - B0_om*B3_d_om))
     + omak*(Om3*(u + pay) + Om3_m*(u_m + pay_m) - Om1*w - Om1_m*w_m)
     + ak*(Om3_o*(u_o + pay_o) + Om3_om*(u_om + pay_om) - Om1_o*w_o - Om1_om*w_om);

      s [6][1]       = -2.0*omam2*(B2_m*B0_d_m + B3_m*B1_d_m
                                   - B0_m*B2_d_m - B1_m*B3_d_m)
 		       - 2.0*amam2*(B2_m*B0_d_om + B3_m*B1_d_om
                                    - B0_m*B2_d_om - B1_m*B3_d_om);
      s [6][ne + 1]  = -2.0*omam2*(B2*B0_d + B3*B1_d 
			           - B0*B2_d - B1*B3_d)
		       - 2.0*amam2*(B2*B0_d_o + B3*B1_o
                                    - B0*B2_d_o - B1*B3_d_o);
      s [6][4]       = -omak*Om2_m;
      s [6][ne + 4]  = -omak*Om2;
      s [6][5]       = omak*Om1_m;
      s [6][ne + 5]  = omak*Om1;
      s [6][6]       = -omak*t_ds;
      s [6][ne + 6]  = omak*t_ds;
      s [6][7]       = -2.0*omam2*ep_m*(B2_m/gdt - B2_d_m)
                       - 2.0*amam2*(ep_om*B2_om/gdt - ep_m*B2_d_om);
      s [6][ne + 7]  = -2.0*omam2*ep*(B2/gdt - B2_d)
                       - 2.0*amam2*(ep_o*B2_o/gdt - ep*B2_d_o);
      s [6][8]       = -2.0*omam2*ep_m*(B3_m/gdt - B3_d_m)
                       - 2.0*amam2*(ep_om*B3_om/gdt - ep_m*B3_d_om);
      s [6][ne + 8]  = -2.0*omam2*ep*(B3/gdt - B3_d)
                       - 2.0*amam2*(ep_o*B3_o/gdt - ep*B3_d_o);
      s [6][9]       = -2.0*omam2*ep_m*(B0_d_m - B0_m/gdt)
                       - 2.0*amam2*(ep_m*B0_d_om - ep_om*B0_om/gdt);
      s [6][ne + 9]  = -2.0*omam2*ep*(B0_d - B0/gdt)
                       - 2.0*amam2*(ep*B0_d_o - ep_o*B0_o/gdt);
      s [6][10]       = -2.0*omam2*ep_m*(B1_d_m - B1_m/gdt)
                       - 2.0*amam2*(ep_m*B1_d_om - ep_om*B1_om/gdt);
      s [6][ne + 10]  = -2.0*omam2*ep*(B1_d - B1/gdt)
                       - 2.0*amam2*(ep*B1_d_o - ep_o*B1_o/gdt);
      s [6][11]      = omak*v_m;
      s [6][ne + 11] = omak*v;
      s [6][12]      = -omak*(u_m + pay_m);
      s [6][ne + 12] = -omak*(u + pay);
   
      s [6][rhs]     = 
       omak*t_ds*(w - w_m) + ak*t_ds*(w_o - w_om)
     - 2.0*omam2*(ep*(B2*B0_d + B3*B1_d 
			     - B0*B2_d - B1*B3_d)
                  + ep_m*(B2_m*B0_d_m + B3_m*B1_d_m
                                 - B0_m*B2_d_m - B1_m*B3_d_m))
     - 2.0*amam2*(ep_o*(B2_o*B0_d + B3_o*B1_d
                               - B0_o*B2_d - B1_o*B3_d)
                  + ep_om*(B2_om*B0_d_m + B3_om*B1_d_m
                                  - B0_om*B2_d_m - B1_om*B3_d_m)
                  + ep*(B2*B0_d_o + B3*B1_d_o 
                               - B0*B2_d_o - B1*B3_d_o)
                  + ep_m*(B2_m*B0_d_om + B3_m*B1_d_om
                                 - B0_m*B2_d_om - B1_m*B3_d_om))
     - 2.0*am2*(ep_o*(B2_o*B0_d_o + B3_o*B1_d_o
                      - B0_o*B2_d_o - B1_o*B3_d_o)
                + ep_om*(B2_om*B0_d_om + B3_om*B1_d_om
                         - B0_om*B2_d_om - B1_om*B3_d_om))
     + omak*(Om1*v + Om1_m*v_m - Om2*(u + pay) - Om2_m*(u_m + pay_m) )
     + ak*(Om1_o*v_o + Om1_om*v_om - Om2_o*(u_o + pay_o) - Om2_om*(u_om + pay_om) );

      s [7][7]      = -omak*t_ds;
      s [7][ne + 7] = omak*t_ds;
      s [7][8]      = omak_2*Om1_m;
      s [7][ne + 8] = omak_2*Om1;
      s [7][9]      = omak_2*Om2_m;
      s [7][ne + 9] = omak_2*Om2;
      s [7][10]      = omak_2*Om3_m;
      s [7][ne + 10] = omak_2*Om3;
      s [7][11]      = omak_2*B1_m;
      s [7][ne + 11] = omak_2*B1;
      s [7][12]      = omak_2*B2_m;
      s [7][ne + 12] = omak_2*B2;
      s [7][13]      = omak_2*B3_m;
      s [7][ne + 13] = omak_2*B3;

      s [7][rhs]    = 
       omak*t_ds*(B0 - B0_m) + ak*t_ds*(B0_o - B0_om)
     + omak_2*(B1*Om1 + B2*Om2 + B3*Om3
               + B1_m*Om1_m + B2_m*Om2_m + B3_m*Om3_m)
     + 0.5*ak*(B1_o*Om1_o + B2_o*Om2_o + B3_o*Om3_o
               + B1_om*Om1_om + B2_om*Om2_om + B3_om*Om3_om);

      s [8][7]      = -omak_2*Om1_m;
      s [8][ne + 7] = -omak_2*Om1;
      s [8][8]      = -omak*t_ds;
      s [8][ne + 8] = omak*t_ds;
      s [8][9]      = -omak_2*Om3_m;
      s [8][ne + 9] = -omak_2*Om3;
      s [8][10]      = omak_2*Om2_m;
      s [8][ne + 10] = omak_2*Om2;
      s [8][11]      = -omak_2*B0_m;
      s [8][ne + 11] = -omak_2*B0;
      s [8][12]      = omak_2*B3_m;
      s [8][ne + 12] = omak_2*B3;
      s [8][13]      = -omak_2*B2_m;
      s [8][ne + 13] = -omak_2*B2;

      s [8][rhs] = 
       omak*t_ds*(B1 - B1_m) + ak*t_ds*(B1_o - B1_om)
     - omak_2*(B0*Om1 - B3*Om2 + B2*Om3
               + B0_m*Om1_m - B3_m*Om2_m + B2_m*Om3_m)
     - 0.5*ak*(B0_o*Om1_o - B3_o*Om2_o + B2_o*Om3_o
               + B0_om*Om1_om - B3_om*Om2_om + B2_om*Om3_om);

      s [9][7]      = -omak_2*Om2_m;
      s [9][ne + 7] = -omak_2*Om2;
      s [9][8]      = omak_2*Om3_m;
      s [9][ne + 8] = omak_2*Om3;
      s [9][9]      = -omak*t_ds;
      s [9][ne + 9] = omak*t_ds;
      s [9][10]      = -omak_2*Om1_m;
      s [9][ne + 10] = -omak_2*Om1;
      s [9][11]      = -omak_2*B3_m;
      s [9][ne + 11] = -omak_2*B3;
      s [9][12]      = -omak_2*B0_m;
      s [9][ne + 12] = -omak_2*B0;
      s [9][13]      = omak_2*B1_m;
      s [9][ne + 13] = omak_2*B1;

      s [9][rhs] = 
       omak*t_ds*(B2 - B2_m) + ak*t_ds*(B2_o - B2_om)
     - omak_2*(B3*Om1 + B0*Om2 - B1*Om3
               + B3_m*Om1_m + B0_m*Om2_m - B1_m*Om3_m)
     - 0.5*ak*(B3_o*Om1_o + B0_o*Om2_o - B1_o*Om3_o
               + B3_om*Om1_om + B0_om*Om2_om - B1_om*Om3_om);

      s [10][7]      = -omak_2*Om3_m;
      s [10][ne + 7] = -omak_2*Om3;
      s [10][8]      = -omak_2*Om2_m;
      s [10][ne + 8] = -omak_2*Om2;
      s [10][9]      = omak_2*Om1_m;
      s [10][ne + 9] = omak_2*Om1;
      s [10][10]      = -omak*t_ds;
      s [10][ne + 10] = omak*t_ds;
      s [10][11]      = omak_2*B2_m; 
      s [10][ne + 11] = omak_2*B2;
      s [10][12]      = -omak_2*B1_m;
      s [10][ne + 12] = -omak_2*B1;
      s [10][13]      = -omak_2*B0_m;
      s [10][ne + 13] = -omak_2*B0;

      s [10][rhs] = 
       omak*t_ds*(B3 - B3_m) + ak*t_ds*(B3_o - B3_om)
     + omak_2*(B2*Om1 - B1*Om2 - B0*Om3
               + B2_m*Om1_m - B1_m*Om2_m - B0_m*Om3_m)
     + 0.5*ak*(B2_o*Om1_o - B1_o*Om2_o - B0_o*Om3_o
               + B2_om*Om1_om - B1_om*Om2_om - B0_om*Om3_om);

      s [11][11]      = -omak;
      s [11][ne + 11] = omak;
      s [11][rhs]     = omak*(Om1 - Om1_m) + ak*(Om1_o - Om1_om);

      s [12][1]       = -3.0*omak*Sb_m*ep_m*ep_m;
      s [12][ne + 1]  = -3.0*omak*Sb*ep*ep;
      s [12][3]       = -omak*elfakm;
      s [12][ne + 3]  = -omak*elfak;
      s [12][11]      = omak*(GJ - EI)*Om3_m;
      s [12][ne + 11] = omak*(GJ - EI)*Om3;
      s [12][12]      = -omak*EI*t_ds;
      s [12][ne + 12] = omak*EI*t_ds;
      s [12][13]      = omak*(GJ - EI)*Om1_m;
      s [12][ne + 13] = omak*(GJ - EI)*Om1;
      s [12][rhs] = 
       omak*EI*t_ds*(Om2 - Om2_m) + ak*EI*t_ds*(Om2_o - Om2_om)
     + omak*((GJ - EI)*(Om1*Om3 + Om1_m*Om3_m)
             - Sb*elfak - Sb_m*elfakm)
     + ak*((GJ - EI)*(Om1_o*Om3_o + Om1_om*Om3_om)
           - Sb_o*elfaok - Sb_om*elfaokm);

      s [13][1]       = 3.0*omak*Sn_m*ep_m*ep_m;
      s [13][ne + 1]  = 3.0*omak*Sn*ep*ep;
      s [13][2]       = omak*elfakm;
      s [13][ne + 2]  = omak*elfak;
      s [13][11]      = omak*(EI - GJ)*Om2_m;
      s [13][ne + 11] = omak*(EI - GJ)*Om2;
      s [13][12]      = omak*(EI - GJ)*Om1_m;
      s [13][ne + 12] = omak*(EI - GJ)*Om1;
      s [13][13]      = -omak*EI*t_ds;
      s [13][ne + 13] = omak*EI*t_ds;
      s [13][rhs] = 
       omak*EI*t_ds*(Om3 - Om3_m) + ak*EI*t_ds*(Om3_o - Om3_om)
     + omak*((EI - GJ)*(Om1*Om2 + Om1_m*Om2_m)
             + Sn*elfak + Sn_m*elfakm)
     + ak*((EI - GJ)*(Om1_o*Om2_o + Om1_om*Om2_om)
           + Sn_o*elfaok + Sn_om*elfaokm);
   }
       
   return;
}

static void 
InitialVelocity (B0, B1, B2, B3, u, v, w, U, V, W)
   double	B0;
   double	B1;
   double	B2;
   double	B3;
   double      *u;
   double      *v;
   double      *w;
   double      *U;
   double      *V;
   double      *W;
{

   if (problem -> type != Towing && problem -> type != Drifter && problem -> type != Deployment) {
      *u = *U = 0.0;
      *v = *V = 0.0;
      *w = *W = 0.0;
   
      return;
   }

   if (problem -> type == Deployment) {
      *V = problem -> terminal [1] -> yspeed.value;
      *W = problem -> terminal [1] -> zspeed.value;
   }
   else {
      *V = problem -> terminal [2] -> yspeed.value;
      *W = problem -> terminal [2] -> zspeed.value;
   }

   *U = 0.0;

   *u = tComponent(0.0,(*V),(*W),B0,B1,B2,B3);
   *v = nComponent(0.0,(*V),(*W),B0,B1,B2,B3);
   *w = bComponent(0.0,(*V),(*W),B0,B1,B2,B3);

   return;
}

int SolveDynamicProblem3D (node, num_nodes, active, num_active, 
			   out, output_map, output_nodes, num_output_nodes, 
			   sample_dt, snap_dt, seg_dt, buoy_dt, ext_dt, decimate,
               prog_name, prog_dt, restart_name, restart_t)
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
   double	 buoy_dt;
   double	 ext_dt;
   double    seg_dt;
   int		 decimate;
   char     *prog_name;
   double    prog_dt;
   char     *restart_name;
   double    restart_t;
{
   static double	scalv [] = {0, 0.07, 1000.0, 1000.0, 1.0, 1.0, 1.0,
                                       0.5, 0.5, 0.5, 0.5, 0.001, 0.001, 0.001};

   int            adapt_count [11];
   int            max_adapt; 
   int            level; 
   double	**s;
   double	  t; 
   double	  dt;
   double	  u, v, w, U, V, W;
   int		  i, j;
   int	  	  it;
   char		  buffer [256];
   int		  displ_msg;
   int		  njn;
   int		  dynstat_conv;
   double	  dt1g, gdt_inv;
   Node		  close;
   Segment   *seg;
   int        nseg;
   double     xthrust, ythrust, zthrust;
   ResFile    prog;
   long       prog_pos;

   problem -> dynamic = 1;
   problem -> twoD = 0;

   ak = analysis -> alpha_k;
   am = analysis -> alpha_m;
   omam = 1.0 - am;
   omak = 1.0 - ak;

   omam2 = omam*omam;
   amam2 = am*omam;
   am2   = am*am;
   omg_g = (1.0 - analysis -> gamma) / analysis -> gamma;
   omak_2 = 0.5*omak;

   // seg = BuildSegmentArray(problem, &nseg, 1);
   seg = problem -> segment;
   nseg = problem -> num_segments;

   njn = problem -> junction_size;

   s = (double **) malloc (sizeof(double *) * (NE + njn*NJ_COMPAT)); s--;

   for (i = 1 ; i <= NE + NJ_COMPAT*njn ; i++) { 
      s [i] = (double *) malloc(sizeof(double) * ((2 + njn)*NE + 1));     
      s [i] --;
   }

	 // copy the static solution 

   for (i = 1 ; i <= num_active ; i++)
      for (j = 1 ; j <= 10 ; j++)
         active[i] -> Ys[j] = active[i] -> Y[j];

        // make sure we have sensible numbers for Ys in all nodes as they
        // get used in the output routines

   for (i = 1 ; i <= num_nodes ; i++) {
      if (!node[i] -> active) {
         if (node[i] -> number < node[i] -> segment -> first_active -> number)
            close = node[i] -> segment -> first_active;
         else
            close = node[i] -> segment -> last_active;

         for (j = 1 ; j <= 10 ; j++)
            node[i] -> Ys[j] = close -> Ys[j];

         node[i] -> x = close -> x;
         node[i] -> y = close -> y;
         node[i] -> z = close -> z;
      }
   }

   for (i = 1 ; i <= num_active ; i++) {
      active[i] -> Y[1] = active[i] -> Y_o[1] = active[i] -> Y_f[1] = active[i] -> Y_o_f[1] = active[i] -> Ys[1]; // e
      active[i] -> Y[2] = active[i] -> Y_o[2] = active[i] -> Y_f[2] = active[i] -> Y_o_f[2] = active[i] -> Ys[2]; // Sn
      active[i] -> Y[3] = active[i] -> Y_o[3] = active[i] -> Y_f[3] = active[i] -> Y_o_f[3] = active[i] -> Ys[3]; // Sb
      active[i] -> Y[7] = active[i] -> Y_o[7] = active[i] -> Y_f[7] = active[i] -> Y_o_f[7] = active[i] -> Ys[4]; // B0
      active[i] -> Y[8] = active[i] -> Y_o[8] = active[i] -> Y_f[8] = active[i] -> Y_o_f[8] = active[i] -> Ys[5]; // B1
      active[i] -> Y[9] = active[i] -> Y_o[9] = active[i] -> Y_f[9] = active[i] -> Y_o_f[9] = active[i] -> Ys[6]; // B2
      active[i] -> Y[10] = active[i] -> Y_o[10] = active[i] -> Y_f[10] = active[i] -> Y_o_f[10] = active[i] -> Ys[7]; // B3
      active[i] -> Y[11] = active[i] -> Y_o[11] = active[i] -> Y_f[11] = active[i] -> Y_o_f[11] = active[i] -> Ys[8]; // Om1
      active[i] -> Y[12] = active[i] -> Y_o[12] = active[i] -> Y_f[12] = active[i] -> Y_o_f[12] = active[i] -> Ys[9]; // Om2
      active[i] -> Y[13] = active[i] -> Y_o[13] = active[i] -> Y_f[13] = active[i] -> Y_o_f[13] = active[i] -> Ys[10]; // Om3

      InitialVelocity (active[i] -> Y [7], active[i] -> Y [8], active[i] -> Y [9], active[i] -> Y [10], &u, &v, &w, &U, &V, &W);

      active[i] -> Y[4] = active[i] -> Y_o[4] = active[i] -> Y_f[4] = active[i] -> Y_o_f[4] = u;
      active[i] -> Y[5] = active[i] -> Y_o[5] = active[i] -> Y_f[5] = active[i] -> Y_o_f[5] = v;
      active[i] -> Y[6] = active[i] -> Y_o[6] = active[i] -> Y_f[6] = active[i] -> Y_o_f[6] = w;

      active[i] -> x_o = active[i] -> x_f = active[i] -> x_o_f = active [i] -> x;
      active[i] -> y_o = active[i] -> y_f = active[i] -> y_o_f = active [i] -> y;
      active[i] -> z_o = active[i] -> z_f = active[i] -> z_o_f = active [i] -> z;

      for (j = 1 ; j <= NE ; j++)
         active[i] -> Yd_o [j] = 0.0;

      active[i] -> xdot = active[i] -> xdot_f = active[i] -> xdot_o_f = active[i] -> xdot_o = U;
      active[i] -> ydot = active[i] -> ydot_f = active[i] -> ydot_o_f = active[i] -> ydot_o = V;
      active[i] -> zdot = active[i] -> zdot_f = active[i] -> zdot_o_f = active[i] -> zdot_o = W;
   }

   if (restart_name && restart_t) {
       prog_pos = LoadProgressPoint(restart_name, problem, restart_t, NE);
       if (prog_pos <= 0) {
            error("could not load progress point from %s", prog_name);
            return 1;
        }
        num_active = MakeNextPrevActive(problem, node, num_nodes, active);
   }

   if (prog_name && prog_dt) {
      prog = res_open(prog_name, "wb");
      if (prog == NULL) {
         error("could not open progress file %s for writing", prog_name);
         return 1;    
      }
   }
   else
      prog = NULL;


   for (i = 1 ; i <= num_active ; i++) {
        // evaluate this now that we have values filled in for any
        // variables we might need to evaluate in the expressions so
        // that the timestep t = start values are correct
      if (active[i] -> next_active && active[i] -> next_active -> position == Connection && active[i] -> segment -> connector) {
         ConnectorThrust(restart_t, active[i] -> segment, active[i], &xthrust, &ythrust, &zthrust); 
         active[i] -> segment -> connector_xthrust.value = xthrust;
         active[i] -> segment -> connector_ythrust.value = ythrust;
         active[i] -> segment -> connector_zthrust.value = zthrust;
         // printf("%d init thrust = %f %f %f\n", i, xthrust, ythrust, zthrust); 
      }
      else if (active[i] -> position == BranchTerminal) {
         Thrust(restart_t, active[i] -> segment -> branch -> terminal, &xthrust, &ythrust, &zthrust);
         active[i] -> segment -> branch -> terminal -> xthrust.value = xthrust;
         active[i] -> segment -> branch -> terminal -> ythrust.value = ythrust;
         active[i] -> segment -> branch -> terminal -> zthrust.value = zthrust;
      }
   }

   Thrust(restart_t, problem -> terminal [1], &xthrust, &ythrust, &zthrust);
   problem -> terminal[1] -> xthrust.value = xthrust;
   problem -> terminal[1] -> ythrust.value = ythrust;
   problem -> terminal[1] -> zthrust.value = zthrust;
   // printf("init thrust = %f %f %f\n", xthrust, ythrust, zthrust); 

   Thrust(restart_t, problem -> terminal [2], &xthrust, &ythrust, &zthrust);
   problem -> terminal[2] -> xthrust.value = xthrust;
   problem -> terminal[2] -> ythrust.value = ythrust;
   problem -> terminal[2] -> zthrust.value = zthrust;
   // printf("init thrust = %f %f %f\n", xthrust, ythrust, zthrust); 

   

   if (problem -> type == Surface || problem -> type == Drifter || problem -> type == Deployment)
      problem -> terminal [2] -> buoy -> draft = environment -> depth - active[num_active] -> x;

   dt = analysis -> dt;

   if (sample_dt == 0.0)
      sample_dt = dt;

   WriteDynamicHeader (out, restart_t, analysis -> duration, dt, sample_dt, 
                       snap_dt, seg_dt, buoy_dt, ext_dt,
                       num_output_nodes, output_nodes, decimate, node);

   if (prog_dt && prog)
      InitializeProgressFile(problem, prog, NE);

   if (snap_dt > 0.0) 
      WriteDynamicSnapshot (out, output_map, node, num_nodes, decimate, 0);
   if (seg_dt > 0.0)
      WriteDynamicSegmentData(problem, out);
   if (ext_dt > 0.0)
      WriteDynamicExternalForces(problem, out);
   
 
   WriteDynamicResult (out, output_map, output_nodes, 
                       num_output_nodes, decimate, 
		               node, num_nodes, 0);  

   max_adapt = analysis->adapt_levels;       /* this is the maximum number of times we will */
                        /* try to reduce dt to get around problems     */
                        /* before bailing out and giving up            */
   
   for (i = 1 ; i <= max_adapt ; i++)
      adapt_count [i] = 0;
         
   level = 0; 
         
	/*	
	 * the result at t = 0.0 is simply the static solution
	 * and we have already written it out so start at t = dt
	 */

   DisplayDynamicHeader ( );


   dynstat_conv = 0;

#if (defined HAVELAMP && !defined API)
   if (environment -> forcing == LAMP)
      UpdateLAMP(node, Y, num_active, 0.0, out, buoy_dt, 0);
#endif

   for (t = restart_t + dt ; t <= analysis -> duration + dt/2.0 ; t += dt) {

#if (defined HAVELAMP && !defined API)
      if (environment -> forcing == LAMP)
         UpdateLAMP(node, Y, num_active, t, out, buoy_dt, 0);
#endif

      if (!level)
         dt = analysis -> dt;

      num_active = ProcessSpools(problem, node, num_nodes, active, num_active, t, dt, NE, 0);
      if (num_active < 2) {
         DisplayMessage("not enough active nodes - all nodes spooled");
         SetError(C_NOMORENODES);
         return 1;
      }

      it = SolveDE (DynamicDifeq3D, DynamicUpdate3D,
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
            DisplayMessage("max adaption level exceeded");
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
            active[j] -> z = active[j] -> z_o;

            active[j] -> x_f = active[j] -> x_o_f;
            active[j] -> y_f = active[j] -> y_o_f;
            active[j] -> z_f = active[j] -> z_o_f;

            for (i = 1 ; i <= NE ; i ++) {
               active[j] -> Y[i] = active[j] -> Y_o[i];
               active[j] -> Y_f[i] = active[j] -> Y_o_f[i];
            }
         }


         DisplayMessage("adapting, dt = %g", dt);

         continue;
      }

      sprintf (buffer,"t = %g,", t);
      displ_msg = 0;

      if (num_output_nodes &&  check(t, dt, sample_dt) <= analysis -> dt/2000) {
         WriteDynamicResult (out, output_map, output_nodes, 
                             num_output_nodes, decimate, 
			     node, num_nodes, 0);
         
         strcat (buffer, " result");
         displ_msg = 1;
      }

      if (snap_dt && check(t, dt, snap_dt) <= analysis -> dt/2000) {
         WriteDynamicSnapshot (out, output_map, 
                               node, num_nodes, decimate, 0);
         strcat (buffer, " snapshot");
         displ_msg = 1;
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
    

      
      if (displ_msg)
         DisplayMessage (buffer);

      if (problem -> dynstat) {
         dynstat_conv = CheckDynstatConvergence(active, num_active, dt, scalv, NE, 0);
         if (dynstat_conv)
            break;
      }

      dt1g = dt*(1.0 - analysis -> gamma);
      gdt_inv = 1.0/analysis -> gamma/dt;

      for (i = 1 ; i <= num_active ; i++) {
         for (j = 1 ; j <= NE ; j++) {
            active[i] -> Yd_o[j] = gdt_inv*(active[i] -> Y [j] - active[i] -> Y_o [j] - dt1g*active[i] -> Yd_o [j]);
            active[i] -> Y_o[j] = active[i] -> Y [j];
            active[i] -> Y_o_f[j] = active[i] -> Y_f [j];
         }

         active[i] -> xdot_o =  gdt_inv*(active[i] -> x - active[i] -> x_o - dt1g*active[i] -> xdot_o);
         active[i] -> ydot_o =  gdt_inv*(active[i] -> y - active[i] -> y_o - dt1g*active[i] -> ydot_o);
         active[i] -> zdot_o =  gdt_inv*(active[i] -> z - active[i] -> z_o - dt1g*active[i] -> zdot_o);

         active[i] -> xdot_o_f = active[i] -> xdot_f;
         active[i] -> ydot_o_f = active[i] -> ydot_f;
         active[i] -> zdot_o_f = active[i] -> zdot_f;

         active[i] -> x_o = active[i] -> x;
         active[i] -> y_o = active[i] -> y;
         active[i] -> z_o = active[i] -> z;

         active[i] -> x_o_f = active[i] -> x_f;
         active[i] -> y_o_f = active[i] -> y_f;
         active[i] -> z_o_f = active[i] -> z_f;

         active[i] -> pay_o = active[i] -> pay;
         active[i] -> pay_o_f = active[i] -> pay_f;
      }

      SmoothNodeData(t, active, num_active, 1);

      if (prog && prog_dt && check(t, dt, prog_dt) <= analysis -> dt/2000) {
          WriteProgress(t, problem, prog, NE);
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

            DisplayMessage("adapting back, dt = %g", dt);

            if (level)
               adapt_count [level] ++;
         }
         else
            adapt_count [level] ++;
      }
   }

   if (prog)
       res_close(prog);

   for (i = 1 ; i <= NE + njn*NJ_COMPAT ; i++) {
      s [i] ++; free (s [i]);
   }

        /*
         * copy the latest dynamic solution into the static
         * solution in case the driver routines needs the
         * latest info for any reason (dynstat option for instance)
         */

   for (i = 1 ; i <= num_nodes ; i++) {
      node[i] -> Ys[1] = node[i] -> Y[1];
      node[i] -> Ys[2] = node[i] -> Y[2];
      node[i] -> Ys[3] = node[i] -> Y[3];
      node[i] -> Ys[4] = node[i] -> Y[4] = node[i] -> Y[7];
      node[i] -> Ys[5] = node[i] -> Y[5] = node[i] -> Y[8];
      node[i] -> Ys[6] = node[i] -> Y[6] = node[i] -> Y[9];
      node[i] -> Ys[7] = node[i] -> Y[7] = node[i] -> Y[10];
      node[i] -> Ys[8] = node[i] -> Y[8] = node[i] -> Y[11];
      node[i] -> Ys[9] = node[i] -> Y[9] = node[i] -> Y[12];
      node[i] -> Ys[10] = node[i] -> Y[10] = node[i] -> Y[13];
   }

   s ++; free (s);

   return 0;
}
