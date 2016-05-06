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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
# include "compress.h"
# include "problem.h"
# include "tension.h"
# include "output.h"
# include "error.h"
# include "solve.h"
# include "transforms.h"

#ifdef HAVELAMP

/* prototypes for LAMP routines */
void initializelamp_(float*,int*,float*,float*,float*,float*,float*,float*,int*);
void advancelamp_(int*,float*,float*,float*,float*,float*,float*,int*);

static float rattch [3];
static double xc0, yc0, zc0;

void TranslateLAMP(x, y, z, lampos)
   double    x, y, z;
   float     lampos [3];
{
   lampos [0] = y - rattch [0] - yc0;		/* horizontal in-plane     */ 
   lampos [1] = z - rattch [1] - zc0;		/* horizontal out-of-plane */
   lampos [2] = x - rattch [2] - xc0;		/* vertical		   */

   return;
}

#ifndef API

int UpdateLAMP(node, Y, nn, tm, out, buoy_dt, twoD)
   Node	    *node;
   double  **Y;
   int       nn;
   double    tm;
   ResFile   out;
   double    buoy_dt;
   int	     twoD;
{
   static int    init = 0;
   static float  lamp_t;
   static int	 istep;
   static int    nstep;
   static float	 tinit, dth;
   Buoy	  	 b;
   double	 B0, B1, B2, B3;
   double	 cph, sph;
   double	 T, Sn, Sb;
   float         fcurnt [3], fcable [3];
   float	 rattch0 [3], vattch [3];
   double	 uack, vack, wack;
   double	 ywind, zwind;
   double	 ycurr, zcurr;
   double	 Fx, Fy, Fz;
   double	 bdr;
   int	 	 ierr;
   float	 motion [6];
   double	 x [6];
   int		 i;
   float	 mass, disp;

   if (!init) {
      rattch0 [0] = 0.0;
      rattch0 [1] = 0.0;
      rattch0 [2] = -1.68;

      initializelamp_(rattch0, &nstep, &tinit, &dth, rattch, vattch, 
		      &mass, &disp, &ierr);
      if (ierr == 1)
         return 1;

      init = 1;
      lamp_t = 0.0;
      istep = 1;

      xc0 = node [nn] -> x;
      yc0 = node [nn] -> y;
      zc0 = node [nn] -> z;
   }

   if (tm - 1e-6 <= (istep - 2)*dth && istep > 1) 
      return 0;

   b = problem.terminal [2] -> buoy;
   b -> draft = environment.surface - node [nn] -> x;
   bdr =  0.5*environment.rho*b -> Cdn*
          ProjectedArea(b, b -> draft);

   WindDrag (tm, b, &ywind, &zwind);
   Current (tm, node [nn] -> x, node [nn] -> y, node [nn] -> z, 
	    &uack, &vack, &wack);

   ycurr = bdr*vack*fabs(vack);
   zcurr = bdr*wack*fabs(wack);

   fcurnt [0] = ycurr + ywind;
   fcurnt [1] = zcurr + zwind;
   fcurnt [2] = 0.0; 

   T = Tension(Y [1][nn], node [nn] -> material);
   Sn = Y [2][nn];

   if (twoD) {
      cph = cos(Y [5][nn]);   
      sph = sin(Y [5][nn]);   

      Fx = T*cph - Sn*sph;
      Fy = T*sph + Sn*cph;
      Fz = 0.0;
   }
   else {
      B0 = Y [7][nn];
      B1 = Y [8][nn];
      B2 = Y [9][nn];
      B3 = Y [10][nn];

      Sb = Y [3][nn];

      Fx = XComponent(T, Sn, Sb, B0, B1, B2, B3); 
      Fy = YComponent(T, Sn, Sb, B0, B1, B2, B3); 
      Fz = ZComponent(T, Sn, Sb, B0, B1, B2, B3); 
   }

   fcable [0] = -Fy;
   fcable [1] = -Fz;
   fcable [2] = -Fx;


   advancelamp_(&istep, fcurnt, fcable, 
               &lamp_t, motion, rattch, vattch, &ierr);
   if (ierr)
      return ierr;

   if (buoy_dt && check(lamp_t, dth, buoy_dt) <= analysis.dt/2000) {
      for (i= 0 ; i < 6 ; i++)
         x [i] = motion [i];

      WriteBuoyMotion(out, x);
   }

   istep ++;

   return 0;
}
#endif /* ifndef API */

#endif /* ifdef HAVELAMP */
