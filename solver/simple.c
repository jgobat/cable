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
 * File:        simple.c
 *
 * Description: contains code to evaluate the dynamic tension prediction
 *		of the simple model described in Gobat's PhD thesis.
 *
 *		Given a static solution from WHOI Cable, this code
 *		evaluates equivalent Cdn, Cdt to calculate model Cd,
 *		nominal mass M_T0 and varphi to calculate model M,
 *		nominal tension T0 and static tension Tst to calculate
 *		tau and dtau, and surface velocity input to calculate
 *		sigA and sigVV.
 *
 * History:
 *
 ****************************************************************************/

# include <stdio.h>
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
# include "segments.h"

extern Problem *problem;             
extern Environment *environment;
extern Analysis *analysis;

int SimpleDynamics (node, num_nodes)
   Node		*node;
   int	   	 num_nodes;
{
   int		 i;
   double	 T0;
   double	 Tst;
   Segment	*seg;
   int		 nseg;
   double	 L;
   double	 M0, M;
   double	 Cdn, Cdt, Cd;
   double	 tau, dtau;
   double	 Lseg;
   Material	 mat;
   double	 H, d;
   double	 varphi;
   double	 t;
   double	 U, V, W;
   double	 prev_U;
   double	 sum_a, sum_a2;
   double	 a, vv;
   double	 sum_vv, sum_vv2;
   Node	         topnode;
   int		 n;
   double	 sigA, sigVV;
   double	 sigT;
   double	 rho;

   // local version of seg array with no branches
   seg = BuildSegmentArray(problem, &nseg, 0);

   H = environment -> depth;
   rho = environment -> rho;

   topnode = node [num_nodes];

   T0 = 0.0;
   Cdn = Cdt = 0.0;
   M0 = 0.0;

   L = environment -> depth - topnode -> x;

   mat = NULL;
   d = 0.0;

   for (i = nseg ; i >= 1 ; i--) {
      Lseg = seg [i] -> length;
      mat  = seg [i] -> material;
      d = mat -> d; 

      if (seg [i] -> connector)
         T0 += seg [i] -> connector -> wet; 

      if (L + Lseg > H)
         break;

      Cdn += Lseg*mat -> d*mat -> Cdn.value/(0.5*rho*d);
      Cdt += Lseg*mat -> d*mat -> Cdt.value/(0.5*rho*d*M_PI);
      M0  += Lseg*(mat -> m + mat -> amt);
      L   += Lseg;
      T0  += Lseg * mat -> wet;
   }

   Lseg = H - L;
   Cdn += Lseg*mat -> d*mat -> Cdn.value/(0.5*rho*d);
   Cdt += Lseg*mat -> d*mat -> Cdt.value/(0.5*rho*d*M_PI);
   M0  += Lseg*(mat -> m + mat -> amt);
   T0  += Lseg * mat -> wet;

   Cdn = Cdn/(d * H);
   Cdt = Cdt/(d * H);

   Cd = 3.79*M_PI*Cdt + 0.46*Cdn;

   Tst = Tension(node[num_nodes] -> Y[1], topnode -> material);
   
   tau = Tst/T0;
   dtau = tau - 1.0;

   L = 0.0;
   for (i = num_nodes ; i >= 1 ; i--) {
      if (node [i] -> x < 0.0)
         break;

      L += node [i] -> ds;
   }

   varphi = H*H/L;

   M = M0 + varphi*(-0.156*(mat -> m + mat -> amt) 
		    + 0.102*(mat -> m + mat -> amn));

   DisplayMessage("tau = %5.3f, Cd = %5.3f, M = %g", tau, Cd, M);

   sum_a = sum_a2 = sum_vv = sum_vv2 = 0.0;
   prev_U = 0.0;

   n = 0;

   for (t = 0.0 ; t <= analysis -> duration ; t += analysis -> dt) {

      if (environment -> forcing== Velocity)
         InputVelocity(t, topnode -> x, topnode -> y, topnode -> z, &U, &V, &W);
      else if (environment -> forcing== WaveFollower)
         WaveSurfaceVelocity(t, topnode -> x, topnode -> y, topnode -> z, &U, &V, &W); 
     
      vv = U*fabs(U);
      sum_vv += vv;
      sum_vv2 += vv*vv;
      n ++;
 
      if (t > 0.0) {
         a = (U - prev_U)/analysis -> dt;

         sum_a += a;
         sum_a2 += a*a; 
      }

      prev_U = U;
   }       

   sum_a /= (double) (n - 1.0);
   sigA = sqrt(sum_a2/(double) (n - 1.0) - sum_a*sum_a); 
 
   sum_vv /= (double) n;
   sigVV = sqrt(sum_vv2/(double) n - sum_vv*sum_vv);

   DisplayMessage("sigVV = %5.3f, sigA = %5.3f", sigVV, sigA);

   sigT = M*tau*sigA + 0.5*environment -> rho*d*H*Cd*dtau*sigVV;

   DisplayMessage("sigT = %g", sigT);

   return 0;
}
