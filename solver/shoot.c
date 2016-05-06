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
 * File:	shoot.c
 *
 * Description: contains code specific to shooting method
 *		solutions for both 2D and 3D static cable models
 *
 *		We abuse the node -> Y array here for ease of use with
 *		the shooting integrations.  In shooting problems, the
 *		x,y,z coordinates of the nodes are part of the main
 *		integration loop, but shear, bending, and euler parameters 
 *		are not (though euler angles are) so we have some extra
 *	  	space.  
 *
 *		Every problem is inherently 3D herein.  We simply make sure 
 *		that out-of-plane forcing is zero when 2D solutions are
 *		requested. The contents in the node -> Y array are:
 *
 *		e, phi, theta, x, y, z
 *
 * History:
 *		
 ****************************************************************************/

# include <stdio.h>
# include <math.h>
# include <malloc.h>
# include "compress.h"
# include "problem.h"
# include "tension.h"
# include "error.h"
# include "output.h"
# include "solve.h"
# include "control.h"
# include "segments.h"

# define SQR(a) ((a)*(a))

extern Problem     *problem;             
extern Environment *environment;
extern Analysis    *analysis;
extern Debug debug;

static int td_node;

static void CheckDebug(Node *node, int num_nodes, ResFile out, int *output_map)
{
   static int	init = 0;

   if (debug.status) {
      if (!init) {
         WriteStaticSolution(node, num_nodes, out, output_map, debug.decimate, 0);

         WriteDynamicHeader (out, 0.0, analysis -> shooting_it, 1.0,
                             debug.sample_it, debug.snap_it, 0.0, 0.0, 0.0, 0, NULL, 1, node);
         init = 1;
      }

      FakeDynamicSnapshot(out, node, num_nodes, debug.decimate, 0);
   }
}

static void rhs(y, n, w0, drat, drap, k, twoD)
   double	y [7];
   Node		n;
   double	w0;
   double	drat;
   double	drap;
   double	k [7];
   int		twoD;
{
   double	e, phi, th, xd, yd, zd;
   double	cp, sp, ct, st;
   double	U, V, W;
   double	T, Tp;
   double	u, v, w;
   double	ep1;
   double	sqep1;
   double	perv;

   e   = y [1];
   phi = y [2];
   th  = y [3];
   xd  = y [4];
   yd  = y [5];
   zd  = y [6];

   cp = cos(phi);
   sp = sin(phi);
   ct = cos(th);
   st = sin(th);

   ep1 = 1.0 + e;
   sqep1 = sqrt(ep1);
 
   Current(0.0, xd, yd, zd, &U, &V, &W);

   if (problem -> type == Towing || problem -> type == Drifter) {
      V -= problem -> terminal [2] -> yspeed.value;
      W -= problem -> terminal [2] -> zspeed.value;
   }
   else if (problem -> type == Deployment) {
      V -= problem -> terminal [1] -> yspeed.value;
      W -= problem -> terminal [1] -> zspeed.value;
   }

   if (twoD)
      W = 0.0;

   T = Tension(e, n -> material);
   Tp = TensionD(e, n -> material);

   u = U*cp*ct + V*sp*ct - W*st;
   v = -sp*U + cp*V;
   w = U*cp*st + V*sp*st + W*ct;

   perv = sqrt(v*v + w*w);

   k [1] = (w0*cp*ct - drat*u*fabs(u)*sqep1) / Tp;
   k [2] = (-w0*sp - drap*v*perv*sqep1) / T / ct;
   k [3] = (-w0*cp*st + drap*w*perv*sqep1) / T;  
   k [4] = ep1*cp*ct;
   k [5] = ep1*sp*ct;
   k [6] = -ep1*st;
 
   return;
}

void 
StaticIntegrate3D (Node *node, int num_nodes, int dir, int twoD)
{
   Node		   n, nm;	
   double	   w0, drat, drap;
   double	   ds;
   int		   i, k;
   Connector	   c;
   Segment         seg;
   double	   w;
   double	   y1 [7];
   double 	   y2 [7];
   double	   y3 [7];
   double	   y4 [7];
   double	   yt [7];
   double	   dragY, dragZ, wet;
   double	   T;
   double	   Fx, Fy, Fz;
   double	   U, V, W;
   double	   ct;
   int		   start, stop;

   td_node = 0;

   if (dir == -1) {
      start = num_nodes;
      stop  = 1;
   }
   else {
      start = 1;
      stop = num_nodes;
   }

   k = start;
   while (k != stop) {
      n = node [k];

      if (dir == -1) {
         nm = n -> prev_active; 
         ds = nm -> ds;
      }
      else {
         ds = n -> ds;
         nm = n -> next_active;
      }
 
      if ((dir == -1 && n -> position == Connection) ||
          (dir == 1 && nm -> position == Connection)) {

         if (dir == -1) {
            c = nm -> segment -> connector;
            seg = nm -> segment;
         }
         else {
            c = n -> segment -> connector;
            seg = n -> segment;
         }

         if (c == NULL) {
            nm -> Y[1] = Strain(Tension(n -> Y[1], n -> material), nm -> material);
            for (i = 2 ; i <= 6 ; i++) 
               nm -> Y[i] = n -> Y[i];
         }
         else {
            Current(0.0, n -> Y[4], n -> Y[5], n -> Y[6], &U, &V, &W);

            if (problem -> type == Towing || problem -> type == Drifter) {
               V -= problem -> terminal [2] -> yspeed.value;
               W -= problem -> terminal [2] -> zspeed.value;
            }
            else if (problem -> type == Deployment) {
               V -= problem -> terminal [1] -> yspeed.value;
               W -= problem -> terminal [1] -> zspeed.value;
            }

            if (twoD)
               W = 0.0;

            dragY = c -> Cdn*V*fabs(V);
            dragZ = c -> Cdn*W*fabs(W);

            wet = c -> wet;

            if (problem -> type == Deployment  || problem -> type == Surface) {
               if (n -> Y[4] > environment -> depth && wet < 0.0)
                  wet = wet*(1.0 + tanh(50.0*(environment -> depth - n -> Y[4]))); 
            }
            else if (n -> Y[4] > environment -> depth)
               wet = c -> m*environment -> gravity;

            ct = cos(n -> Y[3]);

            T = Tension(n -> Y[1], n -> material);

            Fx = -wet - dir*T*cos(n -> Y[2])*ct;
            Fy = dragY - dir*T*sin(n -> Y[2])*ct;
            Fz = dragZ + dir*T*sin(n -> Y[3]);

            if (problem -> type == HorizontalDrifter) {
                Fx += seg -> connector_xforce;
                Fy += seg -> connector_yforce;
                Fz += seg -> connector_zforce;
            }
            nm -> Y[1] = Strain(sqrt(Fx*Fx + Fy*Fy + Fz*Fz), nm -> material);
            nm -> Y[2] = atan2(-dir*Fy, -dir*Fx);
            nm -> Y[3] = asin(dir*Fz/sqrt(Fx*Fx + Fy*Fy + Fz*Fz));
            nm -> Y[4] = n -> Y[4];
            nm -> Y[5] = n -> Y[5];
            nm -> Y[6] = n -> Y[6];
         }
      }
      else if ((n -> Y[4] < 0.0 && environment -> bottom_stiffness && problem -> type != Deployment) || td_node) { 
         if (k > td_node) 
            td_node = k;

         nm -> Y[1] = n -> Y[1];
         nm -> Y[2] = sign(node[td_node] -> Y[2])*M_PI/2.0;
         nm -> Y[3] = node[td_node] -> Y[3];
         nm -> Y[4] = problem -> terminal [1] -> x;
         nm -> Y[5] = n -> Y[5]
                      + dir*sin(n -> Y[2])*cos(n -> Y[3])*ds*(1.0 + n -> Y[1]);
         nm -> Y[6] = n -> Y[6] - dir*sin(n -> Y[3])*ds*(1.0 + n -> Y[1]);
      }
      else {
         w0 = (n -> material -> wet + nm -> material -> wet)/2.0;
         w = (n -> material -> w + nm -> material -> w)/2.0;

         drat   = (n -> material -> Cdt.value + nm -> material -> Cdt.value)/2.0;
         drap   = (n -> material -> Cdn.value + nm -> material -> Cdn.value)/2.0;

         if (n -> attachment) {
            w0 += n -> attachment -> wet/ds;
	        w += n -> attachment -> m*environment -> gravity / ds;
            drat += n -> attachment -> Cdt/ds;
            drap += n -> attachment -> Cdn/ds;
         }

         if (problem -> type == Deployment  || problem -> type == Surface) 
            if (n -> Y[4] > environment -> depth && w0 < 0.0)
               w0 = w0*(1.0 + tanh(50.0*(environment -> depth - n -> Y[4]))); 

         for (i = 1 ; i <= 6 ; i++) yt [i] = n -> Y[i];
         rhs(yt, n, w0, drat, drap, y1, twoD);
         
         for (i = 1 ; i <= 6 ; i++) yt [i] = n -> Y[i] + y1 [i]/2.0;
         rhs(yt, n, w0, drat, drap, y2, twoD);

         for (i = 1 ; i <= 6 ; i++) yt [i] = n -> Y[i] + y2 [i]/2.0;
         rhs(yt, n, w0, drat, drap, y3, twoD);

         for (i = 1 ; i <= 6 ; i++) yt [i] = n -> Y[i] + y3 [i];
         rhs(yt, n, w0, drat, drap, y4, twoD);

         for (i = 1 ; i <= 6 ; i++)
            nm -> Y[i] = n -> Y[i]
                + dir*ds*((y1 [i]  + y4 [i]) / 6.0 + (y2  [i] + y3 [i]) / 3.0);
      }

      k += dir;
   }

   if (!td_node)
      td_node = 1;

   for (i = 1 ; i <= num_nodes ; i++) {
      node [i] -> x = node [i] -> Y [4];
      node [i] -> y = node [i] -> Y [5];
      node [i] -> z = node [i] -> Y [6];
   }

   return;
}

static int 
IterateDraft(Node *node, int num_nodes, Buoy b, int step, int twoD, double Vw, double Ww)
{
   double	  err;
   double     xf, yf, zf;
   double	  x1, x2;
   double	  fx1, fx2, fx3;
   int		  it;
   double	  T, e, phi, theta;
   double	  anchor_x;
   double     Fz_wind, Fy_wind;

   anchor_x = problem -> terminal [1] -> x;


   x1 = fx1 = x2 = fx2 = fx3 = 0.0;

   b -> draft = b -> max_draft;

   for (it = 1 ; it <= analysis -> shooting_it ; it ++) { 
      xf = Buoyancy(b, b -> draft, environment) - b -> w;
      fprintf(stdout,"buoyancy = %g, w=%g, F=%g, d=%.11f, dscale=%.11f\n", 
              Buoyancy(b, b -> draft, environment), b -> w, xf, 
              b -> draft, 
              (b -> draft - b-> min_draft) / (b->max_draft - b -> min_draft));

      WindDrag (0.0, b, &Fy_wind, &Fz_wind);
      if (twoD)
         Fz_wind = 0.0;

      yf = Fy_wind
           + 0.5*environment -> rho*b -> Cdn*ProjectedArea(b, b -> draft)*Vw*fabs(Vw);
      if (twoD)
         zf = 0.0;
      else
         zf = Fz_wind
              + 0.5*environment -> rho*b -> Cdn*ProjectedArea(b, b -> draft)*Ww*fabs(Ww);
       T     = sqrt(xf*xf + yf*yf + zf*zf);
       phi   = atan2(yf, xf);
       theta = asin(-zf/T);

       e = Strain(T, node [num_nodes] -> material);

      // fprintf(stderr,"x_start = %g\n", x_start);
      node[num_nodes] -> Y [1] = e;
      node[num_nodes] -> Y [2] = phi;
      node[num_nodes] -> Y [3] = theta;
      node[num_nodes] -> Y [4] = environment -> depth - b -> draft;
      node[num_nodes] -> Y [5] = 0.0;
      node[num_nodes] -> Y [6] = 0.0;

      StaticIntegrate3D (node, num_nodes, -1, twoD);


      if (it >= 3) 
         err = (node[td_node] -> Y[4] - anchor_x); // node [1] -> ds;
      else
         err = 0;

      DisplayInfo (step, 0, it, err, 0.0);

# if defined (GUI) || defined (WINGUI)
      ControlProcessEvents (analysis -> solution);
      if (analysis -> solution -> userQuit) {
          return 1;
      }
# endif
      if (it > 3 && fabs(err) < analysis -> static_tolerance)
         return 0;

      if (it == 1) {
         x1 = b -> draft;
         fx1 = node[td_node] -> Y[4] - anchor_x;

         b -> draft = b -> min_draft*1.0005;
      }
      else if (it == 2) {
         x2 = b -> draft;
         fx2 = node[td_node] -> Y[4] - anchor_x;

         b -> draft = (x2 + x1)/2.0;   	      /* bisection */
      }
      else {
         fx3 = node[td_node] -> Y[4] - anchor_x; 
   
         if (fx3*fx1 < 0)
            x2 = b -> draft;
         else {
            x1 = b -> draft; 
            fx1 = fx3;
         }

         b -> draft = (x1 + x2)/2.0;
      }

   }

   SetError(C_MAXITERATIONSEXCEEDED);
   return 1;
}


static int 
IterateVertical(Node *node, int num_nodes, double xf, double yf, double zf, int step, int twoD, double xc0, double xc1, double yc0, double zc0)
{
   double	  err;
   double	  x_start;
   double	  x1, x2;
   double	  fx1, fx2, fx3;
   int		  it;
   double	  T, e, phi, theta;
   double	  anchor_x;

   T     = sqrt(xf*xf + yf*yf + zf*zf);
   phi   = atan2(yf, xf);
   theta = asin(-zf/T);

   e = Strain(T, node [num_nodes] -> material);

   anchor_x = problem -> terminal [1] -> x;

   x_start = xc0;

   x1 = fx1 = x2 = fx2 = fx3 = 0.0;

   for (it = 1 ; it <= analysis -> shooting_it ; it ++) { 

      // fprintf(stderr,"x_start = %g\n", x_start);
      node[num_nodes] -> Y [1] = e;
      node[num_nodes] -> Y [2] = phi;
      node[num_nodes] -> Y [3] = theta;
      node[num_nodes] -> Y [4] = x_start;
      node[num_nodes] -> Y [5] = yc0;
      node[num_nodes] -> Y [6] = zc0;

      StaticIntegrate3D (node, num_nodes, -1, twoD);


      if (it >= 3) 
         err = fabs((node[td_node] -> Y[4] - anchor_x)/node [1] -> ds);
      else
         err = 0;

      DisplayInfo (step, 0, it, err, 0.0);

# if defined (GUI) || defined (WINGUI)
      ControlProcessEvents (analysis -> solution);
      if (analysis -> solution -> userQuit) {
          return 1;
      }
# endif
      if (it > 3 && err < analysis -> static_tolerance)
         return 0;

      if (it == 1) {
         x1 = x_start;
         fx1 = node[td_node] -> Y[4] - anchor_x;

         x_start = xc1;

         // fprintf(stdout,"td_node = %d (%g), x1 = %g, new x_start will be %g, fx1 = %g\n", td_node, node[td_node] -> Y[4], x1, x_start, fx1);
      }
      else if (it == 2) {
         x2 = x_start;
         fx2 = node[td_node] -> Y[4] - anchor_x;

         x_start = (x2 + x1)/2.0;   	      /* bisection */
         // fprintf(stdout,"td_node = %d (%g), x2 = %g, new x_start will be %g, fx2 = %g\n", td_node, node[td_node] -> Y[4], x2, x_start, fx2);
      }
      else {
         fx3 = node[td_node] -> Y[4] - anchor_x; 
   
         if (fx3*fx1 < 0)
            x2 = x_start;
         else {
            x1 = x_start; 
            fx1 = fx3;
         }

         x_start = (x1 + x2)/2.0;
         // fprintf(stdout,"x1 = %g, x2 = %g, new x_start will be %g, fx3 = %g\n", x1, x2, x_start, fx3);
      }

   }

   SetError(C_MAXITERATIONSEXCEEDED);
   return 1;
}

static int 
IterateTowing(Node *node, int num_nodes, double xf, double yf, double zf, int step, int twoD)
{
   double	  err;
   double	  x_start;
   double	  x1, x2;
   double	  fx1, fx2, fx3;
   int		  it;
   double	  T, e, phi, theta;
   double	  ship_x;

   T     = sqrt(xf*xf + yf*yf + zf*zf);
   phi   = atan2(yf, xf);
   theta = asin(-zf/T);

   e = Strain(T, node [1] -> material);

   ship_x = environment -> depth;

   x_start = 0.8*(environment -> depth - node [num_nodes] -> s);
   
   x1 = fx1 = x2 = fx2 = fx3 = 0.0;

   for (it = 1 ; it <= analysis -> shooting_it ; it ++) { 

      node[1] -> Y[1] = e;
      node[1] -> Y[2] = phi;
      node[1] -> Y[3] = theta;
      node[1] -> Y[4] = x_start;
      node[1] -> Y[5] = 0.0;
      node[1] -> Y[6] = 0.0;

      StaticIntegrate3D (node, num_nodes, 1, twoD);

      if (it >= 3) 
         err = fabs((node[num_nodes] -> Y[4] - ship_x)/node [num_nodes - 1] -> ds);
      else
         err = 0;

      DisplayInfo (step, 0, it, err, 0.0);

# if defined (GUI) || defined (WINGUI)
      ControlProcessEvents (analysis -> solution);
      if (analysis -> solution -> userQuit) {
          return 1;
      }
# endif
      if (it > 3 && err < analysis -> static_tolerance)
         return 0;

      if (it == 1) {
         x1 = x_start;
         fx1 = node[num_nodes] -> Y[4] - ship_x;

         x_start = environment -> depth;
      }
      else if (it == 2) {
         x2 = x_start;
         fx2 = node[num_nodes] -> Y[4] - ship_x;

         x_start = (x2 + x1)/2.0;   	      /* bisection */
      }
      else {
         fx3 = node[num_nodes] -> Y[4] - ship_x; 
   
         if (fx3*fx1 < 0)
            x2 = x_start;
         else {
            x1 = x_start; 
            fx1 = fx3;
         }

         x_start = (x1 + x2)/2.0;
      }
   }

   SetError(C_MAXITERATIONSEXCEEDED);
   return 1;
}

static void 
ConvertResults2D(Node *node, int num_nodes)
{
   double	  off_y;
   double	  off_x;
   int		  i;

   off_y = node [1] -> y - problem -> terminal [1] -> y;

   if (!environment -> depth)
      off_x = node [1] -> x;
   else
      off_x = 0.0;

   for (i = 1 ; i <= num_nodes ; i++) {
      node[i] -> Y[3] = node[i] -> Y[2];	// move phi 

      node [i] -> x = node [i] -> x - off_x;
      node [i] -> y = node [i] -> y - off_y;
      node [i] -> z = 0.0;
   }

   BendingDerivatives2D(node[1]); /* fills in 2,4 */

   return;
}

static void 
ConvertResults3D(Node *node, int num_nodes)
{
   double	  off_y;
   double	  off_z;
   double	  off_x;
   int		  i;
   double	  th, phi, sthh, cthh, cphh, sphh;

   off_y = node [1] -> y - problem -> terminal [1] -> y;
   off_z = node [1] -> z - problem -> terminal [1] -> z;

   if (!environment -> depth)
      off_x = node [1] -> x;
   else
      off_x = 0.0;

   for (i = 1 ; i <= num_nodes ; i++) {

#if 0
    // this leads to better dynamic results in some out-of-plane cases
    // but singularities in others
     double     a1, a2, a3;
     double     c1, c2, c3;
     double     s1, s2, s3;

      a1 = 0;
      a2 = node[i] -> Y[2];
      a3 = node[i] -> Y[3];

      c1 = cos(a1/2.0);
      c2 = cos(a2/2.0);
      c3 = cos(a3/2.0);
      s1 = sin(a1/2.0);
      s2 = sin(a2/2.0);
      s3 = sin(a3/2.0);

      node[i] -> Y[4] = s1*s2*s3 + c1*c2*c3;
      node[i] -> Y[5] = s1*c2*c3 - s2*s3*c1;
      node[i] -> Y[6] = -s1*s2*c3 + s3*c1*c2;
      node[i] -> Y[7] = s1*s3*c2 + s2*c1*c3;


      double     sp, st, ct, cp;
      double     B0, B1, B2, B3;

      cp = cos(phi);
      ct = cos(th);
      sp = sin(phi);
      st = sin(th);

      B0 = -0.5*sqrt(cp*ct + cp + ct + 1); // -sqrt(cp*ct + cp + ct + 1);
      if (i == 1) {
          B0_1 = B0;
          printf("phi1 = %g, th1 = %g\n", phi, th);
      }
      B0 -= B0_1;

      B1 = -sp*st/4.0/B0;
      B2 = (cp*st + st)/4.0/B0;
      B3 = (sp*ct + sp)/4.0/B0;

        // was - + + + 
      node[i] -> Y[4] = node[i] -> Y[2]; // B0; // -cphh*sthh;
      node[i] -> Y[5] = node[i] -> Y[3]; //  cphh*cthh;
      node[i] -> Y[6] =  sphh*cthh;
      node[i] -> Y[7] =  sphh*sthh;
#endif

        // the way we have always done the conversion. Can lead to 
        // funny pay-out problems with out of 2D plane problems

      phi = node[i] -> Y[2];
      th  = node[i] -> Y[3];
  
      th = -th; 
      sthh = sin(0.5*th);
      cthh = cos(0.5*th);
      cphh = cos(0.5*phi);
      sphh = sin(0.5*phi);

      node[i] -> Y[4] = -cphh*sthh;
      node[i] -> Y[5] =  cphh*cthh;
      node[i] -> Y[6] =  sphh*cthh;
      node[i] -> Y[7] =  sphh*sthh;

      node [i] -> y = node [i] -> y - off_y;
      node [i] -> z = node [i] -> z - off_z;
      node [i] -> x = node [i] -> x - off_x;
   }

   BendingDerivatives3D(node[1]); /* fills in 2,3,8,9,10 */

   return;
}

static void 
ConvertResults(Node *node, int num_nodes, int twoD)
{
   if (twoD)
      ConvertResults2D(node, num_nodes);
   else
      ConvertResults3D(node, num_nodes);

   return;
}

static int 
ShootHorizontal3Dold (node, num_nodes, active, num_active, out, output_map, twoD)
   Node		*node;
   int		 num_nodes;
   Node		*active;
   int		 num_active;
   ResFile	 out;
   int		*output_map;
   int		 twoD;
{
   double	  hf, xf, yf, zf;
   int		  it;
   double	  x_err;
   double	  y_err;
   double	  z_err;
   double	  err;
   double	  rlx;
   double	  ds;
   double	  EA, w0;
   double	  length;
   Segment	 *seg;
   int		  nsegments;
   double	  yc1, yc2;
   double	  zc1, zc2;
   double	  xdist, hdist;
   double     anchor_x;
   double     xc0, xc1, yc0, zc0; // integration start coordinates

   // build local version of seg array without branches
   seg = BuildSegmentArray(problem, &nsegments, 0);
   EquivalentWeight(seg, nsegments, &w0, &EA, &length);

   yc2 = problem -> terminal [2] -> y;
   yc1 = problem -> terminal [1] -> y;

   zc2 = problem -> terminal [2] -> z;
   zc1 = problem -> terminal [1] -> z;

   if (twoD)
      hdist = yc2 - yc1;
   else
      hdist = sqrt(SQR(yc2 - yc1) + SQR(zc2 - zc1));

   xdist = problem -> terminal [2] -> x - problem -> terminal [1] -> x;

   CatenaryForces (w0, EA, length, xdist, hdist, &hf, &xf);
   yf = hf*(yc2 - yc1)/hdist;
   if (twoD)
      zf = 0.0;
   else
      zf = hf*(zc2 - zc1)/hdist;

   
   ds = active [1] -> ds;

   DisplayStaticHeader ( );

   anchor_x = problem -> terminal[1] -> x;
   xc0 = anchor_x + 2.0*active [num_active] -> s;
   xc1 = anchor_x - 2*node[num_nodes]->s; 
   yc0 = problem -> terminal [2] -> y;
   zc0 = problem -> terminal [2] -> z;
 
   for (it = 1 ; it <= analysis -> outer_it ; it++) {
      printf("xf=%g, yf=%g, zf=%g\n", xf, yf, zf);

      if (IterateVertical(active, num_active, xf, yf, zf, it, twoD, 
                          xc0, xc1, yc0, zc0))
         return 1;

      x_err = active[num_active] -> Y[4] - problem -> terminal [2] -> x;
      y_err = (active[num_active] -> Y[5] - active[1] -> Y[5]) - (yc2 - yc1);
      if (twoD)
         z_err = 0.0;
      else
         z_err = (active[num_active] -> Y[6] - active[1] -> Y [6]) - (zc2 - zc1);

      err = sqrt(x_err*x_err + y_err*y_err + z_err*z_err)/ds;
      printf("x=%g, y=%g, z=%g, err=%g\n", active[num_active]->Y[4],
                                   active[num_active]->Y[5] - active[1]->Y[5],
                                   active[num_active]->Y[6] - active[1]->Y[6],
                                   err);

      if (err < analysis -> outer_tolerance)
         break;

      CheckDebug(node, num_nodes, out, output_map);


      rlx = analysis -> outer_relaxation;

      xf -= rlx*x_err;
      yf -= rlx*y_err;
      zf -= rlx*z_err;

      if (problem -> type == HorizontalDrifter) 
          ResolvePositionedConnectors(seg, nsegments, rlx, twoD);
   }

   seg ++;
   free(seg);

   problem -> terminal [2] -> xforce = xf;
   problem -> terminal [2] -> yforce = yf;
   problem -> terminal [2] -> zforce = zf;

   ConvertResults(active, num_active, twoD);

   return 0;
}

static int 
ShootHorizontal3D(node, num_nodes, active, num_active, out, output_map, twoD)
   Node		*node;
   int		 num_nodes;
   Node		*active;
   int		 num_active;
   ResFile	 out;
   int		*output_map;
   int		 twoD;
{
   double	  err;
   int		  it;
   double	  T, e;
   double     phi, theta;
   double	  hf, xf, yf, zf;
   double	  x_err;
   double	  y_err;
   double	  z_err;
   double	  rlx;
   double	  ds;
   double	  EA, w0;
   double	  length;
   Segment	 *seg;
   int		  nsegments;
   double	  yc1, yc2;
   double	  zc1, zc2;
   double	  xdist, hdist;
   double     xc0, yc0, zc0; // integration start coordinates

   // build local version of seg array without branches
   seg = problem -> segment;
   nsegments = problem -> num_segments;

   EquivalentWeight(seg, nsegments, &w0, &EA, &length);

   yc2 = problem -> terminal [2] -> y;
   yc1 = problem -> terminal [1] -> y;

   zc2 = problem -> terminal [2] -> z;
   zc1 = problem -> terminal [1] -> z;

   if (twoD)
      hdist = yc2 - yc1;
   else
      hdist = sqrt(SQR(yc2 - yc1) + SQR(zc2 - zc1));

   xdist = problem -> terminal [2] -> x - problem -> terminal [1] -> x;

   printf("w0 = %g, xd = %g, hd = %g\n", w0, xdist, hdist);
   CatenaryForces (w0, EA, length, xdist, hdist, &hf, &xf);
   printf("cat hf=%g, xf=%g\n", hf, xf);

   yf = hf*(yc2 - yc1)/hdist;
   if (twoD)
      zf = 0.0;
   else
      zf = hf*(zc2 - zc1)/hdist;

   ds = active [1] -> ds;

   DisplayStaticHeader ( );

   xc0 = problem -> terminal[2] -> x;
   yc0 = problem -> terminal[2] -> y;
   zc0 = problem -> terminal[2] -> z;
 
   for (it = 1 ; it <= analysis -> shooting_it ; it++) {
     
      T     = sqrt(xf*xf + yf*yf + zf*zf);
      phi   = atan2(yf, xf);
      theta = asin(-zf/T);

      e = Strain(T, active [num_active] -> material);

      active[num_active] -> Y [1] = e;

      active[num_active] -> Y [2] = phi;
      active[num_active] -> Y [3] = theta;

      active[num_active] -> Y [4] = xc0;
      active[num_active] -> Y [5] = yc0;
      active[num_active] -> Y [6] = zc0;

      StaticIntegrate3D (active, num_active, -1, twoD);

      x_err = active[1] -> Y[4] - problem -> terminal [1] -> x;     
      y_err = active[1] -> Y[5] - problem -> terminal [1] -> y;
      if (twoD)
         z_err = 0.0;
      else
         z_err = active[1] -> Y[6] - problem -> terminal[1] -> z;

      err = sqrt(x_err*x_err + y_err*y_err + z_err*z_err)/ds;
      printf("pos=%g,%g,%g: x=%g, y=%g, z=%g, err=%g\n", 
             problem -> terminal[1] -> x,
             problem -> terminal[1] -> y,
             problem -> terminal[1] -> z,
             active[1]->Y[4],
             active[1]->Y[5],
             active[1]->Y[6],
             err);

      if (err < analysis -> outer_tolerance)
         break;

      CheckDebug(node, num_nodes, out, output_map);

      rlx = analysis -> outer_relaxation;

      xf += rlx*x_err;
      yf += rlx*y_err;
      zf += rlx*z_err;

      if (problem -> type == HorizontalDrifter) 
          ResolvePositionedConnectors(seg, nsegments, rlx, twoD);

      DisplayInfo (1, 0, it, err, 0.0);
   }

   problem -> terminal [2] -> xforce = xf;
   problem -> terminal [2] -> yforce = yf;
   problem -> terminal [2] -> zforce = zf;

   printf("F = %g,%g,%g\n", xf, yf, zf);

   printf("tnode(shoot) = %d %d\n",
          problem -> terminal[2] -> node -> number,
          active[num_active] -> number);
   printf("xyz(shoot) = %f %f %f\n",
           active[num_active] -> x, 
           active[num_active] -> y, 
           active[num_active] -> z);
   ConvertResults(active, num_active, twoD);

   return 0;
}

static int 
ShootWebster3D (node, num_nodes, active, num_active, out, output_map, twoD)
   Node		*node;
   int		 num_nodes;
   Node		*active;
   int		 num_active;
   ResFile	 out;
   int		*output_map;
   int		 twoD;
{
   double	  xf, yf, zf;
   double	  phi;
   int		  it;
   double	  x1, x2;
   double	  fx1, fx2, fx3;
   double	  err;
   double     xc0, xc1;

   phi = M_PI / 2.1;

   x1 = fx1 = x2 = fx2 = fx3 = 0.0;
   
   DisplayStaticHeader ( );

   xc0 = 2.0*active [num_active] -> s;
   xc1 = 1.0*active [num_active] -> s;
   for (it = 1 ; it <= analysis -> outer_it ; it++) {

      xf = problem -> terminal [2] -> xforce 
         = problem -> terminal [2] -> tension * cos(phi);
      yf = problem -> terminal [2] -> yforce 
         = problem -> terminal [2] -> tension * sin(phi);
      if (twoD) {
         zf = problem -> terminal [2] -> zforce = 0.0;
      }
      else {
         zf = problem -> terminal [2] -> zforce 
            = problem -> terminal [2] -> tension * sin(phi);
      }
      if (IterateVertical(active, num_active, xf, yf, zf, it, twoD,
                          xc0, xc1, 0, 0))
         return 1;
 
      err = fabs((active[num_active] -> Y[4] - problem -> terminal [2] -> x)/active [1] -> ds);
      if (err < analysis -> static_tolerance)
         break;

      CheckDebug(node, num_nodes, out, output_map);

    
      if (it == 1) {
         x1 = phi;
         fx1 = (active[num_active] -> Y[4] - problem -> terminal [2] -> x);

         phi = 0.5;
      }
      else if (it == 2) {
         x2 = phi;
         fx2 = (active[num_active] -> Y[4] - problem -> terminal [2] -> x);

         phi = (x2 + x1)/2.0;   	      /* bisection */
      }
      else {
         fx3 = (active[num_active] -> Y[4] - problem -> terminal [2] -> x);

         if (fx3*fx1 < 0)
            x2 = phi;
         else {
            x1 = phi;
            fx1 = fx3;
         }

         phi = (x1 + x2)/2.0;
      }

   }

   problem -> terminal [2] -> phi = phi;

   ConvertResults(active, num_active, twoD); 

   return 0;
}

static int 
ShootGeneral3D (node, num_nodes, active, num_active, out, output_map, twoD)
   Node		*node;
   int		 num_nodes;
   Node		*active;
   int		 num_active;
   ResFile	 out;
   int		*output_map;
   int		 twoD;
{
   double	  xf, yf, zf;
   double     xc0, xc1;

   xf = problem -> terminal [2] -> xforce;
   yf = problem -> terminal [2] -> yforce;
   zf = problem -> terminal [2] -> zforce;

   xc0 = 2*active[num_active] -> s;
   xc1 = active[num_active] -> s;
   if (IterateVertical(active, num_active, xf, yf, twoD ? zf : 0.0, 0, twoD, 
                       xc0, xc1, 0, 0))
      return 1;

   ConvertResults(active, num_active, twoD);

   return 0;
}

static int 
ShootSurfaceTest (node, num_nodes, active, num_active, out, output_map, twoD)
   Node		*node;
   int		 num_nodes;
   Node		*active;
   int		 num_active;
   ResFile	 out;
   int		*output_map;
   int		 twoD;
{
   Buoy	 	  b;
   double	  Uw, Vw, Ww;
   double	  xf, yf, zf;
   double	  Fy_wind, Fz_wind;
   int		  it;

   b = problem -> terminal [2] -> buoy;
   if (problem -> type == Deployment)
      problem -> terminal [1] -> x = environment -> depth; 

   b -> draft = b -> max_draft;
   Current(0.0, environment -> depth, 0.0, 0.0, &Uw, &Vw, &Ww);
   if (problem -> type == Deployment) {
      Vw -= problem -> terminal [1] -> yspeed.value;
      Ww -= problem -> terminal [1] -> zspeed.value;
   }

   if (twoD)
      Ww = 0.0;

   DisplayStaticHeader ( );

   it = IterateDraft(active, num_active, b, 1, twoD, Vw, Ww);

   ConvertResults(active, num_active, twoD);

   xf = Buoyancy(b, b -> draft, environment) - b -> w;

   WindDrag (0.0, b, &Fy_wind, &Fz_wind);
   yf = 0.5*environment -> rho*b -> Cdn*ProjectedArea(b, b -> draft)*Vw*fabs(Vw);
   zf = 0.5*environment -> rho*b -> Cdn*ProjectedArea(b, b -> draft)*Ww*fabs(Ww);

   problem -> terminal [2] -> xforce = xf;
   problem -> terminal [2] -> yforce = yf + Fy_wind;
   problem -> terminal [2] -> zforce = zf + Fz_wind;

   environment -> surface = environment -> depth;

   if (it == 0)
       DisplayMessage("converged, draft = %g", b -> draft);

   if (it == 0)
       DisplayMessage("converged, draft = %g, %g", b -> draft, b -> draft*1.05);

   return it;
}

static int 
ShootSurface3D (node, num_nodes, active, num_active, out, output_map, twoD)
   Node		*node;
   int		 num_nodes;
   Node		*active;
   int		 num_active;
   ResFile	 out;
   int		*output_map;
   int		 twoD;
{
   Buoy	 	  b;
   double	  Uw, Vw, Ww;
   double	  xf, yf, zf;
   double	  Fy_wind, Fz_wind;
   double	  calculated_draft;
   int		  it;
   double	  z1, z2;
   double	  fz1, fz2, fz3;
   double	  err;
   double     xc0, xc1;

   b = problem -> terminal [2] -> buoy;
   if (problem -> type == Deployment)
      problem -> terminal [1] -> x = environment -> depth; 

   b -> draft = b -> max_draft;
   Current(0.0, environment -> depth, 0.0, 0.0, &Uw, &Vw, &Ww);
   if (problem -> type == Deployment) {
      Vw -= problem -> terminal [1] -> yspeed.value;
      Ww -= problem -> terminal [1] -> zspeed.value;
   }

   if (twoD)
      Ww = 0.0;

   z1 = fz1 = z2 = fz2 = fz3 = 0.0;

   DisplayStaticHeader ( );

   xc0 = 2*active[num_active] -> s;
   xc1 = active[num_active] -> s;
   for (it = 1 ; it <= analysis -> outer_it ; it++) {
      xf = Buoyancy(b, b -> draft, environment) - b -> w;
      fprintf(stdout,"buoyancy = %g\n", Buoyancy(b, b -> draft, environment));
      fprintf(stdout,"weight = %g\n", b -> w);

      WindDrag (0.0, b, &Fy_wind, &Fz_wind);
      if (twoD)
         Fz_wind = 0.0;

      yf = Fy_wind
           + 0.5*environment -> rho*b -> Cdn*ProjectedArea(b, b -> draft)*Vw*fabs(Vw);
      if (twoD)
         zf = 0.0;
      else
         zf = Fz_wind
              + 0.5*environment -> rho*b -> Cdn*ProjectedArea(b, b -> draft)*Ww*fabs(Ww);
 
      fprintf(stdout,"%g %g %g\n", xf, yf, zf);
 
      if (IterateVertical (active, num_active, xf, yf, zf, it, twoD,
                           xc0, xc1, 0, 0))
         return 1;
 
      CheckDebug(active, num_active, out, output_map);

      calculated_draft = environment -> depth - active[num_active] -> Y[4];

      err = fabs((calculated_draft - b -> draft) / b -> max_draft);

      DisplayMessage("d' = %g, d = %g, e = %g",
                     b -> draft, calculated_draft, err);

      if (err < analysis -> outer_tolerance)
         break;

      if (it == 1) {
         z1 = b -> draft;
         fz1 = calculated_draft - b -> draft;
         if (fz1 > 0) {
            DisplayMessage("buoy submerged at max draft");            
	    b -> draft = calculated_draft;            
	    SetError(C_BUOYNOTATSURFACE);
            break;         
 	 }

         b -> draft = 1.2*b -> min_draft;
      }
      else if (it == 2) {
        z2 = b -> draft;
        fz2 = calculated_draft - b -> draft;

         b -> draft = (z2 + z1)/2.0;
      }
      else {
         fz3 = calculated_draft - b -> draft;
   
         if (fz3*fz1 < 0) 
            z2 = b -> draft;
         else {
            z1 = b -> draft;
            fz1 = fz3;
         }

         b -> draft = (z2 + z1)/2.0;
      }
   }

   ConvertResults(active, num_active, twoD);

   xf = Buoyancy(b, b -> draft, environment) - b -> w;

   WindDrag (0.0, b, &Fy_wind, &Fz_wind);
   yf = 0.5*environment -> rho*b -> Cdn*ProjectedArea(b, b -> draft)*Vw*fabs(Vw);
   zf = 0.5*environment -> rho*b -> Cdn*ProjectedArea(b, b -> draft)*Ww*fabs(Ww);

   problem -> terminal [2] -> xforce = xf;
   problem -> terminal [2] -> yforce = yf + Fy_wind;
   problem -> terminal [2] -> zforce = zf + Fz_wind;

   environment -> surface = environment -> depth;

   DisplayMessage("converged, draft = %g, %d", b -> draft, b -> draft);

   return 0;
}

static int  ShootSubsurface3D (node, num_nodes, active, num_active, out, output_map, twoD)
   Node		*node;
   int		 num_nodes;
   Node		*active;
   int		 num_active;
   ResFile	 out;
   int		*output_map;
   int		 twoD;
{
   double	  xf, yf, zf;
   int		  it;
   double	  x1, x2;
   double	  fx1, fx2, fx3;
   double	  err;
   double	  buoy_depth;
   Buoy		  b;
   double	  U, V, W;
   double     xc0, xc1;

   b = problem -> terminal [2] -> buoy;

   xf = b -> buoyancy - b -> w;
   yf = zf = 0.0;

   buoy_depth = environment -> depth;
  
   x1 = fx1 = x2 = fx2 = fx3 = 0.0;

   DisplayStaticHeader ( );

   xc0 = 2*active[num_active] -> s;
   xc1 = active[num_active] -> s;
   for (it = 1 ; it <= analysis -> outer_it ; it++) {

      Current(0.0, buoy_depth, 0.0, 0.0, &U, &V, &W);
      yf = 0.5*environment -> rho*b -> S*V*fabs(V)*b -> Cdn;
      if (twoD) {
         W = 0.0;
         zf = 0.0;
      }
      else 
         zf = 0.5*environment -> rho*b -> S*W*fabs(W)*b -> Cdn;

      if (IterateVertical(active, num_active, xf, yf, zf, it, twoD, 
                          xc0, xc1, 0, 0))
         return 1;
  
      err = fabs((active[num_active] -> Y[4] - buoy_depth) / active [num_active - 1] -> ds);
      if (err < analysis -> static_tolerance)
         break; 

      CheckDebug(node, num_nodes, out, output_map);

      if (it == 1) {
         x1 = buoy_depth;
         fx1 = (active[num_active] -> Y[4] - buoy_depth);

         buoy_depth = problem -> terminal [1] -> x;
      }
      else if (it == 2) {
         x2 = buoy_depth;
         fx2 = (active[num_active] -> Y[4] - buoy_depth);

         buoy_depth = (x2 + x1)/2.0;   	      /* bisection */
      }
      else {
         fx3 = (active[num_active] -> Y[4] - buoy_depth);

         if (fx3*fx1 < 0)
            x2 = buoy_depth;
         else {
            x1 = buoy_depth;
	    fx1 = fx3;
 	 }

         buoy_depth = (x1 + x2)/2.0;
      }
   }

   problem -> terminal [2] -> xforce = xf;
   problem -> terminal [2] -> yforce = yf;
   problem -> terminal [2] -> zforce = zf;

   ConvertResults(active, num_active, twoD);

   return 0;
}

static int ShootTowing3D (node, num_nodes, active, num_active, out, output_map, twoD)
   Node		*node;
   int		 num_nodes;
   Node		*active;
   int		 num_active;
   ResFile	 out;
   int		*output_map;
   int		 twoD;
{
   double	  xf, yf, zf;
   int		  it;
   double	  x1, x2;
   double	  fx1, fx2, fx3;
   double	  err;
   double	  tow_start;
   Buoy		  b;
   double	  U, V, W;
   double	  xthrust, ythrust, zthrust;

   b = problem -> terminal [1] -> buoy;

   tow_start = 0.8*(environment -> depth - node [num_nodes] -> s);

   Thrust(0.0, problem -> terminal [1], &xthrust, &ythrust, &zthrust);
   xf = b -> w - b -> buoyancy - xthrust;
  
   x1 = fx1 = x2 = fx2 = fx3 = 0.0;

   DisplayStaticHeader ( );

   yf = zf = 0;

   for (it = 1 ; it <= analysis -> outer_it ; it++) {

      Current(0.0, tow_start, 0.0, 0.0, &U, &V, &W);

      V -= problem -> terminal [2] -> yspeed.value;
      W -= problem -> terminal [2] -> zspeed.value;

      yf = -0.5*environment -> rho*b -> S*V*fabs(V)*b -> Cdn - ythrust;

      if (twoD) {
         W = 0.0;
         zf = 0.0;
      }
      else {
         zf = -0.5*environment -> rho*b -> S*W*fabs(W)*b -> Cdn - zthrust;
      }

      if (IterateTowing(active, num_active, xf, yf, zf, it, twoD))
         return 1;

      err = fabs((active[1] -> Y[4] - tow_start) / active [1] -> ds);
      if (err < analysis -> static_tolerance)
         break;

      CheckDebug(node, num_nodes, out, output_map);

      if (it == 1) {
         x1 = tow_start;
         fx1 = (active[1] -> Y[4] - tow_start);

         tow_start = environment -> depth;
      }
      else if (it == 2) {
         x2 = tow_start;
         fx2 = (active[1] -> Y[4] - tow_start);

         tow_start = (x2 + x1)/2.0;   	      /* bisection */
      }
      else {
         fx3 = (active[1] -> Y[4] - tow_start);

         if (fx3*fx1 < 0)
            x2 = tow_start;
         else {
            x1 = tow_start;
	    fx1 = fx3;
 	 }

         tow_start = (x1 + x2)/2.0;
      }
   }

   problem -> terminal [1] -> xforce = xf;
   problem -> terminal [1] -> yforce = yf;
   problem -> terminal [1] -> zforce = zf;

   ConvertResults(active, num_active, twoD);

   return 0;
}

static int ShootStaticProblem(node, num_nodes, active, num_active, out, output_map, twoD)
   Node		*node;
   int		 num_nodes;
   Node		*active;
   int		 num_active;
   ResFile	 out;
   int		*output_map;
   int		 twoD;
{
   int  status;

   if (problem -> num_branch) {
      error ("cannot get a shooting solution for problems with branches");
      SetError(C_NOSHOOTING);
      return 1;
   }

   switch (problem -> type) {

   case Surface:
   case Deployment:
      status = ShootSurfaceTest(node, num_nodes, active , num_active, out, output_map, twoD);
      return status;
//      if (status == 0) 
//          return 0;
//      else
//        return ShootSurface3D(node, num_nodes, active , num_active, out, output_map, twoD);
      break;

   case General:
      return ShootGeneral3D(node, num_nodes, active , num_active, out, output_map, twoD);
      break;

   case Subsurface:
      return ShootSubsurface3D(node, num_nodes, active , num_active, out, output_map, twoD);
      break;

   case Webster:
      return ShootWebster3D(node, num_nodes, active , num_active, out, output_map, twoD);
      break;

   case Towing:
      return ShootTowing3D(node, num_nodes, active , num_active, out, output_map, twoD);
      break;

   case Horizontal:
   case HorizontalDrifter:
   case Riser:
      return ShootHorizontal3D(node, num_nodes, active , num_active, out, output_map, twoD);
      break;
 
   default:
      error ("no shooting available for this problem type (%d)", problem -> type); 
      SetError(C_NOSHOOTING);
      return 1;
      break;
   }

   return 0;
}

int ShootStaticProblem3D(init, node, num_nodes, active, num_active, out, output_map)
   int		 init;		// not used - here for compat w/SolveStatic
   Node		*node;
   int		 num_nodes;
   Node		*active;
   int		 num_active;
   ResFile	 out;
   int		*output_map;
{
   problem -> twoD = 0;
   problem -> dynamic = 0;

   return ShootStaticProblem(node, num_nodes, active, num_active, out, output_map, 0);
}

int ShootStaticProblem2D(init, node, num_nodes, active, num_active, out, output_map)
   int		 init;		// not used - here for compat w/SolveStatic
   Node		*node;
   int		 num_nodes;
   Node		*active;
   int		 num_active;
   ResFile	 out;
   int		*output_map;
{
   problem -> twoD = 1;
   problem -> dynamic = 0;

   return ShootStaticProblem(node, num_nodes, active, num_active, out, output_map, 1);
}
