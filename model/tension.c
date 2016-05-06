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
 * File:        tension.c
 *
 * Description: defines the constitutive relationships
 *
 * History:
 *
 ****************************************************************************/

# include <stdio.h>
# include <math.h>
# include "compress.h"
# include "problem.h"
# include "tension.h"
# include "code.h"

	/*
	 * T = f(e)
	 */


double Tension(e, m)
   double	e;
   Material	m;
{
   if (m -> type == Linear) 
      return e*m -> EA;
   else if (m -> type == Nonlinear) 
      return EvalCode(m -> T [0].expr, NULL, 0.0, e, 0, 0, 0, CURRNODEDATA);

   return 0.0;
}

double 
NodeTension(Node n)
{
    return Tension(n -> Y[1], n -> material);
}

	/*
	 * dT/de
	 */

double TensionD(e, m)
   double	e;
   Material	m;
{
   if (m -> type == Linear) 
      return m -> EA;
   else if (m -> type == Nonlinear) 
      return EvalCode(m -> T [1].expr, NULL, 0.0, e, 0, 0, 0, CURRNODEDATA);

   return 0.0;
}

	/*
	 * d2T/de2
	 */

double TensionDD(e, m)
   double	e;
   Material	m;
{
   if (m -> type == Linear)
      return 0.0;
   else if (m -> type == Nonlinear) 
      return EvalCode(m -> T [2].expr, NULL, 0.0, e, 0, 0, 0, CURRNODEDATA);

   return 0.0;
}

double Strain (T, m)
   double	T;
   Material	m;
{
   double	e0;
   double	T0;
   double	i;

   if (m -> type == Linear)
      return T / m -> EA;
   else if (m -> type == Nonlinear) {
      e0 = 0.1;
      T0 = Tension(e0, m);
      for (i = 1 ; i <= 50 ; i++) {
         e0 = (T - T0)/TensionD(e0, m) + e0;
         T0 = Tension(e0, m);
         if (fabs(T0 - T) <= T*0.001)
            break;
      }
   
      return e0;
   }

   return 0.0;
}
