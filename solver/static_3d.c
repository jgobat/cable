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
 * File:	static_3d.c
 *
 * Description: contains the 3d static solution specific code 
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
# include "transforms.h"
# include "segments.h"

# define SQR(a) ((a)*(a))

# define TOLERANCE  1e-10

extern Problem     *problem;             
extern Environment *environment;
extern Analysis    *analysis;
extern Debug debug;

static void IntegrateXY(Node start, double x0, double y0, double z0)
{
   Node		 a, a_m;
   double	 stk;
   double	 stkm;
   double	 ds;
   double	 B0, B1, B2, B3;
   double	 B0_m, B1_m, B2_m, B3_m;
 
   a = start;
 
   a -> x = x0;
   a -> y = y0;
   a -> z = z0;
   
   a = a -> next_active;

   while (a && a -> active) {
       a_m = a -> prev_active;

       ds = 0.5*a_m -> ds;

       stkm = ds*(1.0 + a_m -> Y [1]);
       stk  = ds*(1.0 + a -> Y [1]);

       B0   = a -> Y [4];
       B1   = a -> Y [5];
       B2   = a -> Y [6];
       B3   = a -> Y [7];
       B0_m = a_m -> Y [4];
       B1_m = a_m -> Y [5];
       B2_m = a_m -> Y [6];
       B3_m = a_m -> Y [7];

       x0 = a -> x = a_m -> x + XtComponent(stkm,B0_m,B1_m,B2_m,B3_m)
                             + XtComponent(stk,B0,B1,B2,B3);                 
       y0 = a -> y = a_m -> y + YtComponent(stkm,B0_m,B1_m,B2_m,B3_m)
                             + YtComponent(stk,B0,B1,B2,B3);
       z0 = a -> z = a_m -> z + ZtComponent(stkm,B0_m,B1_m,B2_m,B3_m)
                             + ZtComponent(stk,B0,B1,B2,B3);
      a = a -> next_active;
   } 
   
    
   return;
}

void StaticUpdate3D (active, num_active, tm, dt) 
   Node		 *active; 
   int		  num_active; 
   double	  tm;			/* not used		*/
   double	  dt;			/* not used		*/
{
   double	 dx;
   double	 adjustment;
   int		 j;
   Node	 	 from;

   IntegrateXY(active[1], problem -> terminal[1] -> x, problem -> terminal[1] -> y, problem -> terminal[1] -> z);

   for (j = 1 ; j <= problem -> num_branch ; j++) {
      from = problem -> branch[j] -> segment_from -> last_active; 
      IntegrateXY(problem -> branch[j] -> segment[1] -> first_active, from -> x, from -> y, from -> z);
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


void StaticDifeq3D (eq_type, n, nm, ne, rhs, num_rows, s, 
		            node, tm, dt, current_factor)
   EquationType	  eq_type;
   Node		      n, nm;
   int            num_rows;
   int            rhs;
   int            ne;
   double       **s;
   Node          *node;
   double         tm;                   /* not used             */
   double         dt;                   /* not used             */
   double         current_factor;
{
   static double   ulcpre, vlcpre, wlcpre;
   static double   uacpre, vacpre, wacpre;
   Node		   nj;
   double	   e, Sn, Sb, B0, B1, B2, B3, Om1, Om2, Om3;
   double	   e_m, Sn_m, Sb_m, B0_m, B1_m, B2_m, B3_m, Om1_m, Om2_m, Om3_m;
   double	   e_j, Sn_j, Sb_j, B0_j, B1_j, B2_j, B3_j;
   double	   EI, GJ, w0, drat, drap;
   double      w0_m, drat_m, drap_m;
   double	   w, w_m;
   double	   sign_j;
   double	   ds;
   double	   t_ds;
   int		   i, j;
   double      tk, tkm, tdk, tdkm, tddk, tddkm, tkj, tdkj;
   double	   sqstrk, sqstrkm;
   double	   elfak, elfakm;
   double	   ulck, vlck, wlck, ulckm, wlckm, vlckm;
   double 	   uack, vack, wack, uackm, vackm, wackm;
   double	   ftr, ftr_m, ftrn, ftrn_m, ftrb, ftrb_m;
   double	   pervk, pervkm;
   double	   xforst, yforst, zforst;
   double	   xthrust, ythrust, zthrust;
   double	   wet;
   double	   y_wind_drag, z_wind_drag;
   double	   dr;
   Connector   c;
   Buoy		   b;
   double  	   bottom, bottom_m;
   double	   mu;
   double	   Fb, Fb_m;
   Node		   nd;


   w0 = n -> material -> wet;
   EI = n -> material -> EI;
   GJ = n -> material -> GJ;

   e = n -> Y [1];
   Sn = n -> Y [2];
   Sb = n -> Y [3];
   B0 = n -> Y [4];
   B1 = n -> Y [5];
   B2 = n -> Y [6];
   B3 = n -> Y [7];
   Om1 = n -> Y [8];
   Om2 = n -> Y [9];
   Om3 = n -> Y [10];

   tk    = Tension(e, n -> material);
   tdk   = TensionD(e, n -> material);
   tddk  = TensionDD(e, n -> material);

   sqstrk  = sqrt(1.0 + n -> Y[1]);
   elfak  = pow((1.0 + n -> Y[1]), 3.0);
  

   if (eq_type == Cable || eq_type == Junction || eq_type == Connection) {
    
      ds = nm -> ds;
   
      w0_m = nm -> material -> wet;

      e_m = nm -> Y [1];
      Sn_m = nm -> Y [2];
      Sb_m = nm -> Y [3];
      B0_m = nm -> Y [4];
      B1_m = nm -> Y [5];
      B2_m = nm -> Y [6];
      B3_m = nm -> Y [7];
      Om1_m = nm -> Y [8];
      Om2_m = nm -> Y [9];
      Om3_m = nm -> Y [10];

      tkm   = Tension(e_m, nm -> material);
      tdkm  = TensionD(e_m, nm -> material);
      tddkm = TensionDD(e_m, nm -> material);


      if (ds)
         t_ds = 2.0/ds;
      else
         t_ds = 0.0;

      sqstrkm = sqrt(1.0 + nm -> Y [1]);
      elfakm = pow((1.0 + nm -> Y[1]), 3.0);
   }
   else {
      nm = NULL;
      ds = t_ds = 0.0;
      tkm = tdkm = tddkm = 0.0; 
      sqstrkm = elfakm = 0.0;
  
      Om1_m = Om2_m = Om3_m = 0.0;
      B0_m = B1_m = B2_m = B3_m = 0.0;
      e_m = Sn_m = Sb_m = 0.0;

      w0_m = 0.0;
   }

   for (i = 1 ; i <= num_rows ; i++)
      for (j = 1 ; j <= rhs ; j++)
         s [i][j] = BLANK;


	/*
	 * anchor node
	 */

   if (eq_type == BottomBoundary) {
      s [1][4]    = 1.0;
      s [1][rhs]  = B0;

      s [2][8]    = 1.0;
      s [2][rhs]  = Om1;

      s [3][9]    = 1.0;
      s [3][rhs]  = Om2;

      s [4][10]   = 1.0;
      s [4][rhs]  = Om3;

/*
      s[1][4] = 1.0;
      s[1][rhs] = B0 - 0.5;
      s[2][5] = 1.0;
      s[2][rhs] = B1 - 0.5;
      s[3][6] = 1.0;
      s[3][rhs] = B2 - 0.5;
      s[4][7] = 1.0;
      s[4][rhs] = B2 - 0.5;
*/

      if (problem -> type == Towing || problem -> type == Drifter) {

         if (problem -> type == Towing) {
            Current(0.0, node[1] -> x, node[1] -> y, node[1] -> z, &uack, &vack, &wack);

            vack += -problem -> terminal [2] -> yspeed.value;
            wack += -problem -> terminal [2] -> zspeed.value;

            b = problem -> terminal [1] -> buoy;

            problem -> terminal [1] -> yforce =
               -0.5*environment -> rho*b -> S*vack*fabs(vack)*b -> Cdn;
            problem -> terminal [1] -> zforce =
               -0.5*environment -> rho*b -> S*wack*fabs(wack)*b -> Cdn;

            Thrust(0.0, problem -> terminal [1], &xthrust, &ythrust, &zthrust);
            problem -> terminal[1] -> xthrust.value = xthrust;
            problem -> terminal[1] -> ythrust.value = ythrust;
            problem -> terminal[1] -> zthrust.value = zthrust;
         }
         else {
            xthrust = 0.0;
            ythrust = 0.0;
            zthrust = 0.0;
         }
 
         xforst = problem -> terminal [1] -> xforce - xthrust;
         yforst = problem -> terminal [1] -> yforce - ythrust;
         zforst = problem -> terminal [1] -> zforce - zthrust;

         s [5][1] = -tdk*dXdt(B0,B1,B2,B3);
         s [5][2] = -dXdn(B0,B1,B2,B3);
         s [5][3] = -dXdb(B0,B1,B2,B3);
         s [5][4] = -dXdB0(tk,Sn,Sb,B0,B1,B2,B3);
         s [5][5] = -dXdB1(tk,Sn,Sb,B0,B1,B2,B3);
         s [5][6] = -dXdB2(tk,Sn,Sb,B0,B1,B2,B3);
         s [5][7] = -dXdB3(tk,Sn,Sb,B0,B1,B2,B3);
         s [5][rhs] = xforst - XComponent(tk,Sn,Sb,B0,B1,B2,B3);

         s [6][1] = -tdk*dYdt(B0,B1,B2,B3);
         s [6][2] = -dYdn(B0,B1,B2,B3);
         s [6][3] = -dYdb(B0,B1,B2,B3);
         s [6][4] = -dYdB0(tk,Sn,Sb,B0,B1,B2,B3);
         s [6][5] = -dYdB1(tk,Sn,Sb,B0,B1,B2,B3);
         s [6][6] = -dYdB2(tk,Sn,Sb,B0,B1,B2,B3);
         s [6][7] = -dYdB3(tk,Sn,Sb,B0,B1,B2,B3);
         s [6][rhs] = yforst - YComponent(tk,Sn,Sb,B0,B1,B2,B3);

         s [7][1] = -tdk*dZdt(B0,B1,B2,B3);
         s [7][2] = -dZdn(B0,B1,B2,B3);
         s [7][3] = -dZdb(B0,B1,B2,B3);
         s [7][4] = -dZdB0(tk,Sn,Sb,B0,B1,B2,B3);
         s [7][5] = -dZdB1(tk,Sn,Sb,B0,B1,B2,B3);
         s [7][6] = -dZdB2(tk,Sn,Sb,B0,B1,B2,B3);
         s [7][7] = -dZdB3(tk,Sn,Sb,B0,B1,B2,B3);
         s [7][rhs] = zforst - ZComponent(tk,Sn,Sb,B0,B1,B2,B3);
      }
   }


	/* 	
	 * buoy node
	 */

   else if (eq_type == TopBoundary) {
      if (problem -> type == Towing || problem -> type == Drifter) {
         s [1][9]   = 1.0;
         s [1][rhs] = Om2;

         s [2][10]  = 1.0;
         s [2][rhs] = Om3;
/*
         s [3][8]  = 1.0;
         s [3][rhs] = Om1;
*/
         s [3][4]   = -2.0*B0;
         s [3][5]   = -2.0*B1;
         s [3][6]   = -2.0*B2;
         s [3][7]   = -2.0*B3;
         s [3][rhs] = 1.0 - B0*B0 - B1*B1 - B2*B2 - B3*B3;

      }
      else {     
         if (problem -> type == Subsurface && !problem -> dynstat) {
            xforst = problem -> terminal [2] -> xforce;

            b = problem -> terminal [2] -> buoy;
            Current(0.0, n -> x, n -> y, n -> z, &uack, &vack, &wack);

            yforst =  0.5*environment -> rho*b -> S*vack*fabs(vack)*b -> Cdn; 
            zforst =  0.5*environment -> rho*b -> S*wack*fabs(wack)*b -> Cdn; 
         }
         else if (problem -> type == Surface && !problem -> dynstat) {
            xforst = problem -> terminal [2] -> xforce;

            b = problem -> terminal [2] -> buoy;
            WindDrag (0.0, b, &y_wind_drag, &z_wind_drag);

            yforst = problem -> terminal [2] -> yforce + y_wind_drag;
            zforst = problem -> terminal [2] -> zforce + z_wind_drag;
         }
         else {
            xforst = problem -> terminal [2] -> xforce;
            yforst = problem -> terminal [2] -> yforce;
            zforst = problem -> terminal [2] -> zforce;
         }

         s [1][9]   = 1.0;
         s [1][rhs] = Om2;

         s [2][10]  = 1.0;
         s [2][rhs] = Om3;

         s [3][4]   = -2.0*B0;
         s [3][5]   = -2.0*B1;
         s [3][6]   = -2.0*B2;
         s [3][7]   = -2.0*B3;
         s [3][rhs] = 1.0 - B0*B0 - B1*B1 - B2*B2 - B3*B3;

         s [4][1] = -tdk*dXdt(B0,B1,B2,B3);
         s [4][2] = -dXdn(B0,B1,B2,B3);
         s [4][3] = -dXdb(B0,B1,B2,B3);
         s [4][4] = -dXdB0(tk,Sn,Sb,B0,B1,B2,B3);
         s [4][5] = -dXdB1(tk,Sn,Sb,B0,B1,B2,B3);
         s [4][6] = -dXdB2(tk,Sn,Sb,B0,B1,B2,B3);
         s [4][7] = -dXdB3(tk,Sn,Sb,B0,B1,B2,B3);
         s [4][rhs] = xforst - XComponent(tk,Sn,Sb,B0,B1,B2,B3);

         s [5][1] = -tdk*dYdt(B0,B1,B2,B3);
         s [5][2] = -dYdn(B0,B1,B2,B3);
         s [5][3] = -dYdb(B0,B1,B2,B3);
         s [5][4] = -dYdB0(tk,Sn,Sb,B0,B1,B2,B3);
         s [5][5] = -dYdB1(tk,Sn,Sb,B0,B1,B2,B3);
         s [5][6] = -dYdB2(tk,Sn,Sb,B0,B1,B2,B3);
         s [5][7] = -dYdB3(tk,Sn,Sb,B0,B1,B2,B3);
         s [5][rhs] = yforst - YComponent(tk,Sn,Sb,B0,B1,B2,B3);

         s [6][1] = -tdk*dZdt(B0,B1,B2,B3);
         s [6][2] = -dZdn(B0,B1,B2,B3);
         s [6][3] = -dZdb(B0,B1,B2,B3);
         s [6][4] = -dZdB0(tk,Sn,Sb,B0,B1,B2,B3);
         s [6][5] = -dZdB1(tk,Sn,Sb,B0,B1,B2,B3);
         s [6][6] = -dZdB2(tk,Sn,Sb,B0,B1,B2,B3);
         s [6][7] = -dZdB3(tk,Sn,Sb,B0,B1,B2,B3);
         s [6][rhs] = zforst - ZComponent(tk,Sn,Sb,B0,B1,B2,B3);
      }
   }
   else if (eq_type == BranchStart) {
      s [1][4]   = 1.0;        
      s [1][rhs] = B0;

      s [2][8]   = 1.0;         
      s [2][rhs] = Om1;        

      s [3][9]   = 1.0;
      s [3][rhs] = Om2;

      s [4][10]  = 1.0;
      s [4][rhs] = Om3;
   }
   else if (eq_type == BranchTerminal) {
      s [1][9]   = 1.0;
      s [1][rhs] = Om2;

      s [2][10]  = 1.0;
      s [2][rhs] = Om3;

      s [3][4]   = -2.0*B0;
      s [3][5]   = -2.0*B1;
      s [3][6]   = -2.0*B2;
      s [3][7]   = -2.0*B3;
      s [3][rhs] = 1.0 - B0*B0 - B1*B1 - B2*B2 - B3*B3;


      if (n -> segment -> branch -> terminal -> loop_main_node) {
         nd = n -> segment -> branch -> terminal -> loop_main_node;

         n -> segment -> branch -> terminal -> yforce -=
             (n -> y - nd -> y);
         n -> segment -> branch -> terminal -> xforce -=
             (n -> x - nd -> x);
         n -> segment -> branch -> terminal -> zforce -=
             (n -> z - nd -> z);

         xforst = n -> segment -> branch -> terminal -> xforce;
         yforst = n -> segment -> branch -> terminal -> yforce;
         zforst = n -> segment -> branch -> terminal -> zforce;
      }
      else if (n -> segment -> branch -> terminal -> anchor != NULL
               || problem -> type == HorizontalDrifter) {

         xforst = n -> segment -> branch -> terminal -> xforce;
         yforst = n -> segment -> branch -> terminal -> yforce;
         zforst = n -> segment -> branch -> terminal -> zforce;
      }
      else {
         b = n -> segment -> branch -> terminal -> buoy;
         Current(0.0, n -> x, n -> y, n -> z, &uack, &vack, &wack);

         xforst = Buoyancy(b, environment -> surface - n -> x, environment) - b -> w;
         yforst = 0.5*environment -> rho*vack*fabs(vack)*b -> Cdn*
                  ProjectedArea(b, environment -> surface - n -> x);
         zforst = 0.5*environment -> rho*wack*fabs(wack)*b -> Cdn*
                  ProjectedArea(b, environment -> surface - n -> x);
      }

      s [4][1] = -tdk*dXdt(B0,B1,B2,B3);
      s [4][2] = -dXdn(B0,B1,B2,B3);
      s [4][3] = -dXdb(B0,B1,B2,B3);
      s [4][4] = -dXdB0(tk,Sn,Sb,B0,B1,B2,B3);
      s [4][5] = -dXdB1(tk,Sn,Sb,B0,B1,B2,B3);
      s [4][6] = -dXdB2(tk,Sn,Sb,B0,B1,B2,B3);
      s [4][7] = -dXdB3(tk,Sn,Sb,B0,B1,B2,B3);
      s [4][rhs] = xforst - XComponent(tk,Sn,Sb,B0,B1,B2,B3);

      s [5][1] = -tdk*dYdt(B0,B1,B2,B3);
      s [5][2] = -dYdn(B0,B1,B2,B3);
      s [5][3] = -dYdb(B0,B1,B2,B3);
      s [5][4] = -dYdB0(tk,Sn,Sb,B0,B1,B2,B3);
      s [5][5] = -dYdB1(tk,Sn,Sb,B0,B1,B2,B3);
      s [5][6] = -dYdB2(tk,Sn,Sb,B0,B1,B2,B3);
      s [5][7] = -dYdB3(tk,Sn,Sb,B0,B1,B2,B3);
      s [5][rhs] = yforst - YComponent(tk,Sn,Sb,B0,B1,B2,B3);

      s [6][1] = -tdk*dZdt(B0,B1,B2,B3);
      s [6][2] = -dZdn(B0,B1,B2,B3);
      s [6][3] = -dZdb(B0,B1,B2,B3);
      s [6][4] = -dZdB0(tk,Sn,Sb,B0,B1,B2,B3);
      s [6][5] = -dZdB1(tk,Sn,Sb,B0,B1,B2,B3);
      s [6][6] = -dZdB2(tk,Sn,Sb,B0,B1,B2,B3);
      s [6][7] = -dZdB3(tk,Sn,Sb,B0,B1,B2,B3);
      s [6][rhs] = zforst - ZComponent(tk,Sn,Sb,B0,B1,B2,B3);
   }


	/*
	 * node between two distinct segments
	 */


   else if (eq_type == Connection || eq_type == Junction) {
      if (nm -> segment -> connector == NULL &&
          nm -> segment -> connection == Spliced) {	/* no lumped mass	*/
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
            dr = c -> Cdn;

            Current (0.0, n -> x, n -> x, n -> y, &uack, &vack, &wack);

            uack *= current_factor;
            vack *= current_factor;
            wack *= current_factor;

            if (problem -> type == Towing || problem -> type == Drifter) {
               vack -= problem -> terminal [2] -> yspeed.value;
               wack -= problem -> terminal [2] -> zspeed.value;
            }
            else if (problem -> type == Deployment) {
               vack -= problem -> terminal [1] -> yspeed.value;
               wack -= problem -> terminal [1] -> zspeed.value;
            }

            wet = c -> wet;
            if ((problem -> type == Deployment  || problem -> type == Surface 
                 || problem -> type == Horizontal || problem -> type == HorizontalDrifter) && !problem -> dynstat) {

               if (n -> x > environment -> surface && wet < 0.0)
                  wet = wet*(1.0 + tanh(50.0*(environment -> surface - n -> x)));
            }

            if (problem -> type == HorizontalDrifter) {
                xforst = nm -> segment -> connector_xforce;
                yforst = nm -> segment -> connector_yforce;
                zforst = nm -> segment -> connector_zforce;
            }
            else {
                xforst = yforst = zforst = 0;
            }
         }
         else {
            wet = dr = 0;
            uack = vack = wack = 0;
            xforst = yforst = zforst = 0;
         }

         s [1][1]      = -tdkm*dXdt(B0_m,B1_m,B2_m,B3_m);
         s [1][ne + 1] = tdk*dXdt(B0,B1,B2,B3);
         s [1][2]      = -dXdn(B0_m,B1_m,B2_m,B3_m);
         s [1][ne + 2] = dXdn(B0,B1,B2,B3);
         s [1][3]      = -dXdb(B0_m,B1_m,B2_m,B3_m);
         s [1][ne + 3] = dXdb(B0,B1,B2,B3);
         s [1][4]      = -dXdB0(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
         s [1][ne + 4] = dXdB0(tk,Sn,Sb,B0,B1,B2,B3);
         s [1][5]      = -dXdB1(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
         s [1][ne + 5] = dXdB1(tk,Sn,Sb,B0,B1,B2,B3);
         s [1][6]      = -dXdB2(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
         s [1][ne + 6] = dXdB2(tk,Sn,Sb,B0,B1,B2,B3);
         s [1][7]      = -dXdB3(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
         s [1][ne + 7] = dXdB3(tk,Sn,Sb,B0,B1,B2,B3);
         s [1][rhs] = XComponent(tk,Sn,Sb,B0,B1,B2,B3)
                      - XComponent(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m)
                      - wet + xforst; 

         s [2][1]      = -tdkm*dYdt(B0_m,B1_m,B2_m,B3_m);
         s [2][ne + 1] = tdk*dYdt(B0,B1,B2,B3);
         s [2][2]      = -dYdn(B0_m,B1_m,B2_m,B3_m);
         s [2][ne + 2] = dYdn(B0,B1,B2,B3);
         s [2][3]      = -dYdb(B0_m,B1_m,B2_m,B3_m);
         s [2][ne + 3] = dYdb(B0,B1,B2,B3);
         s [2][4]      = -dYdB0(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
         s [2][ne + 4] = dYdB0(tk,Sn,Sb,B0,B1,B2,B3);
         s [2][5]      = -dYdB1(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
         s [2][ne + 5] = dYdB1(tk,Sn,Sb,B0,B1,B2,B3);
         s [2][6]      = -dYdB2(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
         s [2][ne + 6] = dYdB2(tk,Sn,Sb,B0,B1,B2,B3);
         s [2][7]      = -dYdB3(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
         s [2][ne + 7] = dYdB3(tk,Sn,Sb,B0,B1,B2,B3);
         s [2][rhs] = dr*vack*fabs(vack) + yforst
                      - YComponent(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m)
                      + YComponent(tk,Sn,Sb,B0,B1,B2,B3);

         s [3][1]      = -tdkm*dZdt(B0_m,B1_m,B2_m,B3_m);
         s [3][ne + 1] = tdk*dZdt(B0,B1,B2,B3);
         s [3][2]      = -dZdn(B0_m,B1_m,B2_m,B3_m);
         s [3][ne + 2] = dZdn(B0,B1,B2,B3);
         s [3][3]      = -dZdb(B0_m,B1_m,B2_m,B3_m);
         s [3][ne + 3] = dZdb(B0,B1,B2,B3);
         s [3][4]      = -dZdB0(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
         s [3][ne + 4] = dZdB0(tk,Sn,Sb,B0,B1,B2,B3);
         s [3][5]      = -dZdB1(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
         s [3][ne + 5] = dZdB1(tk,Sn,Sb,B0,B1,B2,B3);
         s [3][6]      = -dZdB2(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
         s [3][ne + 6] = dZdB2(tk,Sn,Sb,B0,B1,B2,B3);
         s [3][7]      = -dZdB3(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m);
         s [3][ne + 7] = dZdB3(tk,Sn,Sb,B0,B1,B2,B3);
         s [3][rhs] = dr*wack*fabs(wack) + zforst
                      - ZComponent(tkm,Sn_m,Sb_m,B0_m,B1_m,B2_m,B3_m)
                      + ZComponent(tk,Sn,Sb,B0,B1,B2,B3);

         
         for (i = 1 ; i <= n -> segment -> junction.num_nodes ; i++) {
            
            nj = n -> segment -> junction.node [i];
            
            e_j = nj -> Y [1];
            Sn_j = nj -> Y [2];
            Sb_j = nj -> Y [3];
            B0_j = nj -> Y [4];
            B1_j = nj -> Y [5];
            B2_j = nj -> Y [6];
            B3_j = nj -> Y [7];

            tkj    = Tension(e_j, nj -> material);
            tdkj   = TensionD(e_j, nj -> material);

            if (nj == nj -> segment -> branch -> last) 
               sign_j = -1.0;
            else
               sign_j = 1.0;
 
            s [1][(i + 1)*ne + 1] = sign_j*tdkj*dXdt(B0_j,B1_j,B2_j,B3_j);
            s [1][(i + 1)*ne + 2] = sign_j*dXdn(B0_j,B1_j,B2_j,B3_j);
            s [1][(i + 1)*ne + 3] = sign_j*dXdb(B0_j,B1_j,B2_j,B3_j);
            s [1][(i + 1)*ne + 4] = sign_j*dXdB0(tkj,Sn_j,Sb_j,
                                                 B0_j,B1_j,B2_j,B3_j);
            s [1][(i + 1)*ne + 5] = sign_j*dXdB1(tkj,Sn_j,Sb_j,
                                                 B0_j,B1_j,B2_j,B3_j);
            s [1][(i + 1)*ne + 6] = sign_j*dXdB2(tkj,Sn_j,Sb_j,
                                                 B0_j,B1_j,B2_j,B3_j);
            s [1][(i + 1)*ne + 7] = sign_j*dXdB3(tkj,Sn_j,Sb_j,
                                                 B0_j,B1_j,B2_j,B3_j);
            s [1][rhs] += sign_j*XComponent(tkj,Sn_j,Sb_j,B0_j,B1_j,B2_j,B3_j);

            s [2][(i + 1)*ne + 1] = sign_j*tdkj*dYdt(B0_j,B1_j,B2_j,B3_j);
            s [2][(i + 1)*ne + 2] = sign_j*dYdn(B0_j,B1_j,B2_j,B3_j);
            s [2][(i + 1)*ne + 3] = sign_j*dYdb(B0_j,B1_j,B2_j,B3_j);
            s [2][(i + 1)*ne + 4] = sign_j*dYdB0(tkj,Sn_j,Sb_j,
                                                 B0_j,B1_j,B2_j,B3_j);
            s [2][(i + 1)*ne + 5] = sign_j*dYdB1(tkj,Sn_j,Sb_j,
                                                 B0_j,B1_j,B2_j,B3_j);
            s [2][(i + 1)*ne + 6] = sign_j*dYdB2(tkj,Sn_j,Sb_j,
                                                 B0_j,B1_j,B2_j,B3_j);
            s [2][(i + 1)*ne + 7] = sign_j*dYdB3(tkj,Sn_j,Sb_j,
                                                 B0_j,B1_j,B2_j,B3_j);
            s [2][rhs] += sign_j*YComponent(tkj,Sn_j,Sb_j,B0_j,B1_j,B2_j,B3_j);

            s [3][(i + 1)*ne + 1] = sign_j*tdkj*dZdt(B0_j,B1_j,B2_j,B3_j);
            s [3][(i + 1)*ne + 2] = sign_j*dZdn(B0_j,B1_j,B2_j,B3_j);
            s [3][(i + 1)*ne + 3] = sign_j*dZdb(B0_j,B1_j,B2_j,B3_j);
            s [3][(i + 1)*ne + 4] = sign_j*dZdB0(tkj,Sn_j,Sb_j,
                                                 B0_j,B1_j,B2_j,B3_j);
            s [3][(i + 1)*ne + 5] = sign_j*dZdB1(tkj,Sn_j,Sb_j,
                                                 B0_j,B1_j,B2_j,B3_j);
            s [3][(i + 1)*ne + 6] = sign_j*dZdB2(tkj,Sn_j,Sb_j,
                                                 B0_j,B1_j,B2_j,B3_j);
            s [3][(i + 1)*ne + 7] = sign_j*dZdB3(tkj,Sn_j,Sb_j,
                                                 B0_j,B1_j,B2_j,B3_j);
            s [3][rhs] += sign_j*ZComponent(tkj,Sn_j,Sb_j,B0_j,B1_j,B2_j,B3_j);
         }

         s [4][ne + 4] = 1.0;        /* reset B0 above the mass */
         s [4][rhs] = B0;

         s [5][ne + 8] = 1.0;        /* release the moments  */
         s [5][rhs] = Om1;           /* above the mass       */

         s [6][ne + 9] = 1.0;
         s [6][rhs] = Om2;
 
         s [7][ne + 10] = 1.0;
         s [7][rhs] = Om3;

         s [8][4] = 2.0*B0_m;        /* enforce Euler condition */
         s [8][5] = 2.0*B1_m;        /* below the mass          */
         s [8][6] = 2.0*B2_m;
         s [8][7] = 2.0*B3_m;
         s [8][rhs] = B0_m*B0_m + B1_m*B1_m + B2_m*B2_m + B3_m*B3_m - 1.0;

         s [9][9] = 1.0;                 /* release the moments     */
         s [9][rhs] = Om2_m;             /* below the mass          */
 
         s [10][10] = 1.0;
         s [10][rhs] = Om3_m;
      }
   }

	/*
	 * a regular old internal node
	 */

   else {
      Current (0.0, n -> x, n -> y, n -> z, &uack, &vack, &wack);
      
      uack *= current_factor;
      vack *= current_factor;
      wack *= current_factor;

      if (problem -> type == Towing || problem -> type == Drifter) {
         vack -= problem -> terminal [2] -> yspeed.value;
         wack -= problem -> terminal [2] -> zspeed.value;
      }
      else if (problem -> type == Deployment) {
         vack -= problem -> terminal [1] -> yspeed.value;
         wack -= problem -> terminal [1] -> zspeed.value;
      }

      ulck = tComponent(uack,vack,wack,B0,B1,B2,B3);
      vlck = nComponent(uack,vack,wack,B0,B1,B2,B3);
      wlck = bComponent(uack,vack,wack,B0,B1,B2,B3);

      if (n -> active_number == 2 || nm -> active_number != n -> active_number - 1) {
         Current (0.0, nm -> x, nm -> y, nm -> z, &uackm, &vackm, &wackm);

         uackm *= current_factor;
         vackm *= current_factor;
         wackm *= current_factor;

         if (problem -> type == Towing || problem -> type == Drifter) {
            vackm -= problem -> terminal [2] -> yspeed.value;
            wackm -= problem -> terminal [2] -> zspeed.value;
         }
         else if (problem -> type == Deployment) {
            vackm -= problem -> terminal [1] -> yspeed.value;
            wackm -= problem -> terminal [1] -> zspeed.value;
         }

         ulckm = tComponent(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m);
         vlckm = nComponent(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m);
         wlckm = bComponent(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m);
      }
      else {
         ulckm=ulcpre;
         vlckm=vlcpre;
         wlckm=wlcpre;
         uackm=uacpre;
         vackm=vacpre;
         wackm=wacpre;
      }

	// this doesn't feel like the right treatment for drag with
	// a dependence on normal velocity

      DragCoeff(n, n -> material, 0.0, fabs(ulck), sqrt(vlck*vlck + wlck*wlck), &drat, &drap);
      DragCoeff(nm, nm -> material, 0.0, fabs(ulckm), sqrt(vlckm*vlckm + wlckm*wlckm), &drat_m, &drap_m);

      ulcpre=ulck;
      vlcpre=vlck;
      wlcpre=wlck;
      uacpre=uack;
      vacpre=vack;
      wacpre=wack;

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

      pervk = sqrt(wlck*wlck + vlck*vlck);
      if (pervk < TOLERANCE) {
         drap = 0.0;
         pervk = 1.0;
      }

      pervkm = sqrt(wlckm*wlckm + vlckm*vlckm);
      if (pervkm < TOLERANCE) {
         drap_m = 0.0;  
         pervkm = 1.0;
      }

        /*
         * allow for the possibility of line floating on the surface
         */

      if ((problem -> type == Surface || problem -> type == Deployment) && !problem -> dynstat) {
         if (nm -> x >= environment -> surface && w0_m < 0.0)
            w0_m = w0_m*(1.0 + tanh(50.0*(environment -> surface - nm -> x)));

         if (n -> x >= environment -> surface && w0 < 0.0)
            w0 = w0*(1.0 + tanh(50.0*(environment -> surface - n -> x)));
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

      bottom = Bottom(n -> y, n -> z, 0);
      bottom_m = Bottom(nm -> y, nm -> z, 0);
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

      ftr   = 2.0*drat*sqstrk*sign(ulck)*ulck;
      ftr_m = 2.0*drat_m*sqstrkm*sign(ulckm)*ulckm;

      s [1][1]       = (tddkm*(e - e_m) - (tdk + tdkm))/ds
                       + 0.5*drat_m*sqstrkm*ulckm*fabs(ulckm)/sqstrkm;
      s [1][ne + 1]  = (tddk*(e - e_m) + (tdk + tdkm))/ds
                       + 0.5*drat*sqstrk*ulck*fabs(ulck)/sqstrk;
      s [1][2]       = -Om3_m;
      s [1][ne + 2]  = -Om3;
      s [1][3]       = Om2_m;
      s [1][ne + 3]  = Om2;
      s [1][4]       = -2.0*(w0_m - Fb_m)*B0_m
                       + ftr_m*dtdB0(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m);
      s [1][ne + 4]  = -2.0*(w0 - Fb)*B0
                       + ftr*dtdB0(uack,vack,wack,B0,B1,B2,B3);
      s [1][5]       = -2.0*(w0_m - Fb_m)*B1_m
                       + ftr_m*dtdB1(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m);
      s [1][ne + 5]  = -2.0*(w0 - Fb)*B1
                       + ftr*dtdB1(uack,vack,wack,B0,B1,B2,B3);
      s [1][6]       = 2.0*(w0_m - Fb_m)*B2_m
                       + ftr_m*dtdB2(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m);
      s [1][ne + 6]  = 2.0*(w0 - Fb)*B2
                       + ftr*dtdB2(uack,vack,wack,B0,B1,B2,B3);
      s [1][7]       = 2.0*(w0_m - Fb_m)*B3_m
                       + ftr_m*dtdB3(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m);
      s [1][ne + 7]  = 2.0*(w0 - Fb)*B3
                       + ftr*dtdB3(uack,vack,wack,B0,B1,B2,B3);
      s [1][9]       = Sb_m;
      s [1][ne + 9]  = Sb;
      s [1][10]      = -Sn_m;
      s [1][ne + 10] = -Sn;

      s [1][rhs] = (tdk + tdkm)/ds*(e - e_m) 
         + Sb*Om2 - Sn*Om3 + Sb_m*Om2_m - Sn_m*Om3_m
         - (w0 - Fb)*(B0*B0 + B1*B1 - B2*B2 - B3*B3)
         - (w0_m - Fb_m)*(B0_m*B0_m + B1_m*B1_m - B2_m*B2_m - B3_m*B3_m)
	 - mu*(Fb + Fb_m)
         + drat_m*ulckm*fabs(ulckm)*sqstrkm
         + drat*ulck*fabs(ulck)*sqstrk; 

      ftrn_m = drap_m*sqstrkm*(pervkm + vlckm*vlckm/pervkm);
      ftrn   = drap*sqstrk*(pervk + vlck*vlck/pervk);
      ftrb_m = drap_m*sqstrkm*vlckm*wlckm/pervkm;
      ftrb   = drap*sqstrk*vlck*wlck/pervk;
 
      s [2][1]      = tdkm*Om3_m + drap_m*vlckm*pervkm*0.5/sqstrkm;
      s [2][ne + 1] = tdk*Om3 + drap*vlck*pervk*0.5/sqstrk;
      s [2][2]      = -t_ds;
      s [2][ne + 2] = t_ds;
      s [2][3]      = -Om1_m;
      s [2][ne + 3] = -Om1;
      s [2][4]      = 2.0*(w0_m - Fb_m)*B3_m
                      + ftrn_m*dndB0(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m)
                      + ftrb_m*dbdB0(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m);
      s [2][ne + 4] = 2.0*(w0 - Fb)*B3
                      + ftrn*dndB0(uack,vack,wack,B0,B1,B2,B3)
                      + ftrb*dbdB0(uack,vack,wack,B0,B1,B2,B3);
      s [2][5]      = -2.0*(w0_m - Fb_m)*B2_m
                      + ftrn_m*dndB1(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m)
                      + ftrb_m*dbdB1(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m);
      s [2][ne + 5] = -2.0*(w0 - Fb)*B2
                      + ftrn*dndB1(uack,vack,wack,B0,B1,B2,B3)
                      + ftrb*dbdB1(uack,vack,wack,B0,B1,B2,B3);
      s [2][6]      = -2.0*(w0_m - Fb_m)*B1_m
                      + ftrn_m*dndB2(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m)
                      + ftrb_m*dbdB2(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m);
      s [2][ne + 6] = -2.0*(w0 - Fb)*B1
                      + ftrn*dndB2(uack,vack,wack,B0,B1,B2,B3)
                      + ftrb*dbdB2(uack,vack,wack,B0,B1,B2,B3);
      s [2][7]      = 2.0*(w0_m - Fb_m)*B0_m
                      + ftrn_m*dndB3(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m)
                      + ftrb_m*dbdB3(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m);
      s [2][ne + 7] = 2.0*(w0 - Fb)*B0
                      + ftrn*dndB3(uack,vack,wack,B0,B1,B2,B3)
                      + ftrb*dbdB3(uack,vack,wack,B0,B1,B2,B3);
      s [2][8]       = -Sb_m;
      s [2][ne + 8]  = -Sb;
      s [2][10]      = tkm;
      s [2][ne + 10] = tk;
      s [2][rhs] = t_ds*(Sn - Sn_m) 
		  + tk*Om3 - Sb*Om1 + tkm*Om3_m - Sb_m*Om1_m
                  - 2.0*(w0 - Fb)*(B1*B2 - B0*B3) 
                  - 2.0*(w0_m - Fb_m)*(B1_m*B2_m - B0_m*B3_m) 
                  + drap_m*vlckm*pervkm*sqstrkm
                  + drap*vlck*pervk*sqstrk;
 
      ftrb_m = drap_m*sqstrkm*(pervkm + wlckm*wlckm/pervkm);
      ftrb   = drap*sqstrk*(pervk + wlck*wlck/pervk);
      ftrn_m = drap_m*sqstrkm*vlckm*wlckm/pervkm;
      ftrn   = drap*sqstrk*vlck*wlck/pervk;

      s [3][1]      = -tdkm*Om2_m + drap_m*wlckm*pervkm*0.5/sqstrkm;
      s [3][ne + 1] = -tdk*Om2 + drap*wlck*pervk*0.5/sqstrk;
      s [3][2]      = Om1_m;
      s [3][ne + 2] = Om1;
      s [3][3]      = -t_ds;
      s [3][ne + 3] = t_ds;
      s [3][4]      = -2.0*(w0_m - Fb_m)*B2_m
                      + ftrb_m*dbdB0(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m)
                      + ftrn_m*dndB0(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m);
      s [3][ne + 4] = -2.0*(w0 - Fb)*B2
                      + ftrb*dbdB0(uack,vack,wack,B0,B1,B2,B3)
                      + ftrn*dndB0(uack,vack,wack,B0,B1,B2,B3);
      s [3][5]      = -2.0*(w0_m - Fb_m)*B3_m
                      + ftrb_m*dbdB1(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m)
                      + ftrn_m*dndB1(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m);
      s [3][ne + 5] = -2.0*(w0 - Fb)*B3
                      + ftrb*dbdB1(uack,vack,wack,B0,B1,B2,B3)
                      + ftrn*dndB1(uack,vack,wack,B0,B1,B2,B3);
      s [3][6]      = -2.0*(w0_m - Fb_m)*B0_m
                      + ftrb_m*dbdB2(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m)
                      + ftrn_m*dndB2(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m);
      s [3][ne + 6] = -2.0*(w0 - Fb)*B0
                      + ftrb*dbdB2(uack,vack,wack,B0,B1,B2,B3)
                      + ftrn*dndB2(uack,vack,wack,B0,B1,B2,B3);
      s [3][7]      = -2.0*(w0_m - Fb_m)*B1_m
                      + ftrb_m*dbdB3(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m)
                      + ftrn_m*dndB3(uackm,vackm,wackm,B0_m,B1_m,B2_m,B3_m);
      s [3][ne + 7] = -2.0*(w0 - Fb)*B1
                      + ftrb*dbdB3(uack,vack,wack,B0,B1,B2,B3)
                      + ftrn*dndB3(uack,vack,wack,B0,B1,B2,B3);
      s [3][8]      = Sn_m;
      s [3][ne + 8] = Sn;
      s [3][9]      = -tkm;
      s [3][ne + 9] = -tk;
      s [3][rhs] = t_ds*(Sb - Sb_m) 
                  + Sn*Om1 - tk*Om2 + Sn_m*Om1_m - tkm*Om2_m
                  - 2.0*(w0 - Fb)*(B1*B3 + B0*B2)
                  - 2.0*(w0_m - Fb_m)*(B1_m*B3_m + B0_m*B2_m)
                  + drap_m*wlckm*pervkm*sqstrkm
                  + drap*wlck*pervk*sqstrk;

      s [4][4]       = -t_ds;
      s [4][ne + 4]  = t_ds;
      s [4][5]       = 0.5*Om1_m;
      s [4][ne + 5]  = 0.5*Om1;
      s [4][6]       = 0.5*Om2_m;
      s [4][ne + 6]  = 0.5*Om2;
      s [4][7]       = 0.5*Om3_m;
      s [4][ne + 7]  = 0.5*Om3;
      s [4][8]       = 0.5*B1_m;
      s [4][ne + 8]  = 0.5*B1;
      s [4][9]       = 0.5*B2_m;
      s [4][ne + 9]  = 0.5*B2;
      s [4][10]      = 0.5*B3_m;
      s [4][ne + 10] = 0.5*B3;
      s [4][rhs] = t_ds*(B0 - B0_m) 
                   + 0.5*(B1*Om1 + B2*Om2 + B3*Om3 
                          + B1_m*Om1_m + B2_m*Om2_m + B3_m*Om3_m);
 
      s [5][4]       = -0.5*Om1_m;
      s [5][ne + 4]  = -0.5*Om1;
      s [5][5]       = -t_ds;
      s [5][ne + 5]  = t_ds;
      s [5][6]       = -0.5*Om3_m;
      s [5][ne + 6]  = -0.5*Om3;
      s [5][7]       = 0.5*Om2_m;
      s [5][ne + 7]  = 0.5*Om2;
      s [5][8]       = -0.5*B0_m;
      s [5][ne + 8]  = -0.5*B0;
      s [5][9]       = 0.5*B3_m;
      s [5][ne + 9]  = 0.5*B3;
      s [5][10]      = -0.5*B2_m;
      s [5][ne + 10] = -0.5*B2;
      s [5][rhs] = t_ds*(B1 - B1_m)
                  - 0.5*(B0*Om1 - B3*Om2 + B2*Om3
                         + B0_m*Om1_m - B3_m*Om2_m + B2_m*Om3_m);
 
      s [6][4]      = -0.5*Om2_m;
      s [6][ne + 4] = -0.5*Om2;
      s [6][5]      = 0.5*Om3_m;
      s [6][ne + 5] = 0.5*Om3;
      s [6][6]      = -t_ds;
      s [6][ne + 6] = t_ds;
      s [6][7]      = -0.5*Om1_m;
      s [6][ne + 7] = -0.5*Om1;
      s [6][8]      = -0.5*B3_m;
      s [6][ne + 8] = -0.5*B3;
      s [6][9]      = -0.5*B0_m;
      s [6][ne + 9] = -0.5*B0;
      s [6][10]      = 0.5*B1_m;
      s [6][ne + 10] = 0.5*B1;
      s [6][rhs] = t_ds*(B2 - B2_m)
                  - 0.5*(B3*Om1 + B0*Om2 - B1*Om3
                         + B3_m*Om1_m + B0_m*Om2_m - B1_m*Om3_m);
 
      s [7][4]      = -0.5*Om3_m;
      s [7][ne + 4] = -0.5*Om3;
      s [7][5]      = -0.5*Om2_m;
      s [7][ne + 5] = -0.5*Om2;
      s [7][6]      = 0.5*Om1_m;
      s [7][ne + 6] = 0.5*Om1;
      s [7][7]      = -t_ds;
      s [7][ne + 7] = t_ds;
      s [7][8]      = 0.5*B2_m;
      s [7][ne + 8] = 0.5*B2;
      s [7][9]      = -0.5*B1_m;
      s [7][ne + 9] = -0.5*B1;
      s [7][10]      = -0.5*B0_m;
      s [7][ne + 10] = -0.5*B0;
      s [7][rhs] = t_ds*(B3 - B3_m)
                  + 0.5*(B2*Om1 - B1*Om2 - B0*Om3
                         + B2_m*Om1_m - B1_m*Om2_m - B0_m*Om3_m);
 
      s [8][8]      = -1.0;
      s [8][ne + 8] = 1.0;
      s [8][rhs]    = Om1 - Om1_m;

      s [9][1]       = -Sb_m*3.0*(1.0 + e_m)*(1.0 + e_m);
      s [9][ne + 1]  = -Sb*3.0*(1.0 + e)*(1.0 + e);
      s [9][3]       = -elfakm;
      s [9][ne + 3]  = -elfak;
      s [9][8]       = (GJ - EI)*Om3_m;
      s [9][ne + 8]  = (GJ - EI)*Om3;
      s [9][9]       = -t_ds*EI;
      s [9][ne + 9]  = t_ds*EI;
      s [9][10]      = (GJ - EI)*Om1_m;
      s [9][ne + 10] = (GJ - EI)*Om1;
      s [9][rhs] = t_ds*EI*(Om2 - Om2_m) 
                   + (GJ - EI)*(Om1*Om3 + Om1_m*Om3_m) 
                   - (Sb*elfak + Sb_m*elfakm);
 
      s [10][1]      = Sn_m*3.0*(1.0 + e_m)*(1.0 + e_m);
      s [10][ne + 1] = Sn*3.0*(1.0 + e)*(1.0 + e);
      s [10][2]      = elfakm;
      s [10][ne + 2] = elfak;
      s [10][8]      = (EI - GJ)*Om2_m;
      s [10][ne + 8] = (EI - GJ)*Om2;
      s [10][9]      = (EI - GJ)*Om1_m;
      s [10][ne + 9] = (EI - GJ)*Om1;
      s [10][10]      = -t_ds*EI;
      s [10][ne + 10] = t_ds*EI;
      s [10][rhs] = t_ds*EI*(Om3 - Om3_m)
                   + (EI - GJ)*(Om1*Om2 + Om1_m*Om2_m) 
                   + (Sn*elfak + Sn_m*elfakm);
   }

   return;
}

void BendingDerivatives3D (Node start)
{
   Node	    n, np, nm;
   double   db0ds, db1ds, db2ds, db3ds;
   double   slopel, sloper;
   double   dOm2, dOm3;
   double   EI;
   double   elfak;

   n = start;

   while (n) {
      np = n -> next_active;
      nm = n -> prev_active;

      if (n -> position == BottomBoundary ||
          n -> position == Junction ||
          n -> position == Connection ||
          n -> position == BranchStart) {

         db0ds = (np -> Y [4] - n -> Y [4])/n -> ds;
         db1ds = (np -> Y [5] - n -> Y [5])/n -> ds;
         db2ds = (np -> Y [6] - n -> Y [6])/n -> ds;
         db3ds = (np -> Y [7] - n -> Y [7])/n -> ds;
      }
      else if (n -> segment -> last_active == n ||
               n -> position == TopBoundary || np == NULL) {
         db0ds = (n -> Y [4] - nm -> Y [4])/nm -> ds;
         db1ds = (n -> Y [5] - nm -> Y [5])/nm -> ds;
         db2ds = (n -> Y [6] - nm -> Y [6])/nm -> ds;
         db3ds = (n -> Y [7] - nm -> Y [7])/nm -> ds;
      }
      else {
         db0ds = (np -> Y [4] - nm -> Y [4])/(nm -> ds + n -> ds);
         db1ds = (np -> Y [5] - nm -> Y [5])/(nm -> ds + n -> ds);
         db2ds = (np -> Y [6] - nm -> Y [6])/(nm -> ds + n -> ds);
         db3ds = (np -> Y [7] - nm -> Y [7])/(nm -> ds + n -> ds);
      }

      n -> Y [8] = 2.0*(-n -> Y [5]*db0ds + n -> Y [4]*db1ds
                        + n -> Y [7]*db2ds - n -> Y [6]*db3ds);
      n -> Y [9] = 2.0*(-n -> Y [6]*db0ds - n -> Y [7]*db1ds 
                        + n -> Y [4]*db2ds + n -> Y [5]*db3ds);
      n -> Y [10] = 2.0*(-n -> Y [7]*db0ds + n -> Y [6]*db1ds 
                         - n -> Y [5]*db2ds + n -> Y [4]*db3ds);

      n = n -> next_active;
   }

   n = start;
   while (n) {

      np = n -> next_active;
      nm = n -> prev_active;

      if (n -> position == BottomBoundary ||
          n -> position == Junction ||
          n -> position == Connection ||
          n -> position == BranchStart) {

         dOm2 = (np -> Y [9] - n -> Y [9]) / n -> ds;
         dOm3 = (np -> Y [10] - n -> Y [10]) / n -> ds;
      }
      else if (n -> segment -> last_active == n || 
               n -> position == TopBoundary || np == NULL) {

         dOm2 = (n -> Y [9] - nm -> Y [9]) / nm -> ds;
         dOm3 = (n -> Y [10] - nm -> Y [10]) / nm -> ds;
      }
      else {
         sloper = (np -> Y [9] - n -> Y [9])/n -> ds;
         slopel = (n -> Y [9] - nm -> Y [9])/nm -> ds;
         dOm2 = (sloper*nm -> ds + slopel*n -> ds)/
                    (nm -> ds + n -> ds);
         sloper = (np -> Y [10] - n -> Y [10])/n -> ds;
         slopel = (n -> Y [10] - nm -> Y [10])/nm -> ds;
         dOm3 = (sloper*nm -> ds + slopel*n -> ds)/
                    (nm -> ds + n -> ds);
      }

      EI = n -> material -> EI;
      elfak = pow(1.0 + n -> Y[1], 3.0);

      n -> Y [2] = -EI*dOm2/elfak;
      n -> Y [3] = EI*dOm3/elfak;

      n = n -> next_active;
   }

   return;
}

void SolveCatenary3D(start, x0, y0, z0, hforst, xforst, th, w0, length)
   Node          start;
   double	 x0, y0, z0;
   double        hforst;
   double        xforst;
   double        th;
   double        w0;
   double        length;
{
   Node		n;
   double	tens;
   double	p;
   double	term0, term1, term2, term3, term4, term5;
   double	phi, sthh, cthh, cphh, sphh;
   double	sv;

   if (hforst == 0.0)
      hforst = 0.1;
   if (xforst == 0.0)
      xforst = -1.0;

   if (w0 == 0.0)
      w0 = 1e-6;
 
      
   sthh = sin(0.5*th);
   cthh = cos(0.5*th);

   term0 = hforst / w0;
   term1 = xforst / hforst;
   term2 = term1 - length/term0;
   term4 = sqrt(1.0 + term2*term2);
   term5 = log(term2 + term4);

   sv = 0;
   n = start;
   while(n) {


      term3 = term2 + w0*sv/hforst;
      p = term0*( log(term3 + sqrt(term3*term3 + 1.0)) - term5 );
      tens = sqrt(hforst*hforst + SQR(xforst - w0*(length - sv)));

      phi = atan2(1.0, term3);

/*
      if ((xforst - w0*(length - sv)) < 0.0 && !mainline) 
         phi = M_PI + phi;
*/
      cphh = cos(0.5*phi);
      sphh = sin(0.5*phi);
      n -> Y [4] = -cphh*sthh;
      n -> Y [5] = cphh*cthh;
      n -> Y [6] = sphh*cthh;
      n -> Y [7] = sphh*sthh;

      n -> Y [1] = Strain (tens, n -> material);
      n -> Y [2] = 0.0;
      n -> Y [3] = 0.0;

      sv += n -> ds;

      n = n -> next;
   }   

   BendingDerivatives3D(start);
   IntegrateXY(start, x0, y0, z0);
      
   return;
}

int SolveStaticProblem3D (init_loaded, node, num_nodes, active, num_active, out, output_map)
   int		 init_loaded;
   Node		*node;
   int		 num_nodes;
   Node		*active;
   int		 num_active;
   ResFile	 out;
   int		*output_map;
{
   static double	scalv [] = {0, 0.07, 1000.0, 1000.0, 
                                       0.5, 0.5, 0.5, 0.5, 0.001, 0.001, 0.001};

   int            i;
   double         current_factor;
   double	**s;
   int		  singular;
   int		  resolved;
   double	  resolve_err;
   int		  it;
   int		  ne, nb;
   int		  nj_compat, nb_branch;
   int		  current_steps;
   int		  njn;
   Segment   *seg;
   int        nseg;

   problem -> twoD = 0;
   problem -> dynamic = 0;

   ne = 10;
   nj_compat = 0;
   nb_branch = 4;

   njn = problem -> junction_size;

   // seg = BuildSegmentArray(problem, &nseg, 1);
   seg = problem -> segment;
   nseg = problem -> num_segments;

   if (problem -> type == Towing || problem -> type == Drifter)
      nb = 7;
   else
      nb = 4;

   s = (double **) malloc (sizeof(double *) * (ne + njn*nj_compat)); s--;

   for (i = 1 ; i <= ne + njn*nj_compat ; i++) {
      s [i] = (double *) malloc(sizeof(double) * ((2 + njn)*ne + 1));
      s [i] --;
   }

   if (!init_loaded) {
      if (analysis -> static_initial_guess == Catenary || analysis -> static_solution == Catenary) {

         if (problem -> type == Surface || problem -> type == Deployment)
            problem -> terminal [2] -> buoy -> draft = 
                            problem -> terminal [2] -> buoy -> max_draft;

         i = InitialGuess (SolveCatenary3D, active, num_active); 
         if (i)
           return 1;

         if (analysis -> static_solution == Catenary) 
            return 0;
      }
      else {
         i = ShootStaticProblem3D (0, node, num_nodes, active, num_active, out, output_map);

         if (i) {
            SetError(C_INITIALSOLUTIONFAILED);
            return 1;
         }
      }
   }

   if (debug.status) {
      WriteStaticSolution(node, num_nodes, out, output_map, debug.decimate, 0);
      WriteDynamicHeader (out, 0.0, analysis -> static_it, 1.0,
                         debug.sample_it, debug.snap_it, 0.0, 0.0, 0.0, 0, NULL, 1, node);
      FakeDynamicSnapshot(out, node, num_nodes, debug.decimate, 0);
   }

   DisplayMessage("Initial forces = %g, %g, %g", 
            problem -> terminal [2] -> xforce,
            problem -> terminal [2] -> yforce,
            problem -> terminal [2] -> zforce);

   resolved = 0;
   singular = 0;

   DisplayStaticHeader ( );


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

         singular = SolveDE (StaticDifeq3D, StaticUpdate3D,  
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

      if (singular && problem -> type != Surface)
        break;

      if (problem -> dynstat) {
         resolved = 1;
         break;
      } 

      if (problem -> type == Surface || problem -> type == Deployment) {
         resolve_err = ResolveBuoy (active, num_active, it, 0);

         if (resolve_err < 0.0)
            break;
         else if (resolve_err < analysis -> outer_tolerance)
            resolved = 1;
/*
         else if (problem -> type == Surface 
                  && analysis -> static_initial_guess == Catenary)
            InitialGuess (SolveCatenary3D, active, num_active);
*/
         DisplayMessage("draft = %6.4f, F = (%7.2f, %7.2f)",
                  problem -> terminal [2] -> buoy -> draft,
                  problem -> terminal [2] -> xforce,
                  problem -> terminal [2] -> yforce);
      }
      else if (problem -> type == Horizontal 
               || problem -> type == HorizontalDrifter
               || problem -> type == Riser) {
         resolved = ResolveAnchor (active, num_active, seg, nseg, 0);

         DisplayMessage("z [n] = %6.2f, x [n] = %6.2f, y [n] = %6.2f",
                  active [num_active] -> x, 
                  active [num_active] -> y, 
                  active [num_active] -> z);
      }
      else if (problem -> type == Drifter) {
         resolved = ResolveSpeed3D (node, NULL, num_active, it);

         DisplayMessage("Vsp = %6.4f, Wsp = %6.4f, F = (%7.2f, %7.2f, %7.2f)",
                  problem -> terminal [2] -> yspeed.value,
                  problem -> terminal [2] -> zspeed.value,
                  problem -> terminal [2] -> xforce,
                  problem -> terminal [2] -> yforce,
		  problem -> terminal [2] -> zforce);
      }
      else 
         resolved = 1;
      
      if (resolved)
         break;
   }

   if (singular && singular != analysis -> static_it + 1) {
      DisplayMessage("singularity in static solution");
      SetError(C_STATICSOLUTIONFAILED);
      return 1;
   }

   if (!resolved) {
      SetError(C_MAXITERATIONSEXCEEDED);
      DisplayMessage("never converged on outer iterations");
   }

   for (i = 1 ; i <= ne + njn*nj_compat ; i++) {
      s [i] ++; free (s [i]);
   }
   s ++; free (s);

   if (problem -> type == Horizontal || problem -> type == HorizontalDrifter || problem -> type == Riser || problem -> type == General)  {
      DisplayMessage("Done. Fh2= %.0f, Fv2= %.0f, xyz=%g,%g,%g",
              problem -> terminal [2] -> yforce,
              problem -> terminal [2] -> xforce,
              problem -> terminal [2] -> node -> y,
              problem -> terminal [2] -> node -> z,
              problem -> terminal [2] -> node -> x);
   }

   if (problem -> type == Towing) 
      DisplayMessage("Done. Tow T= %.0f, Tow H= %.2f, Ship T= %.0f",
               Tension(node[1] -> Y[1], node [1] -> material),
               (environment -> depth ? environment -> depth - node [1] -> x :
                node [num_active] -> x),
               Tension(node[num_active] -> Y[1], node [num_active] -> material));
   

   if (problem -> type == Surface || problem -> type == Deployment) 
      DisplayMessage("Done. Buoy draft= %.3f", 
                     problem -> terminal [2] -> buoy -> draft);
   

   if (problem -> type == Drifter) 
      DisplayMessage("Done. |U|= %.2f, Ts=%.1f, Hs= %.0f, Tb= %.0f",
               sqrt(SQR(problem -> terminal [2] -> yspeed.value) +
                    SQR(problem -> terminal [2] -> zspeed.value)),
               Tension(node[1] -> Y[1], node[1] -> material),
               node [num_active] -> x,
               Tension(node [num_active] -> Y [1], node [num_active] -> material));

   return 0;
}
