/*
    This file is part of WHOI Cable, a program for the static and dynamic
    analysis of oceanographic cable structures.

    Copyright (C) 1997-2016 by Woods Hole Oceanographic Institution (WHOI)

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
 * File:        rotate.c
 *
 * Description: routines to rotate vector results between local and global
 *		coordinate systems for two- and three-dimensionals problems
 *
 ****************************************************************************/

# include <math.h>

void RotateToFixed (g1, g2, g3, b0, b1, b2, b3, gx, gy, gz, twoD)
   double	g1, g2, g3;
   double	b0, b1, b2, b3;
   double	*gx, *gy, *gz;
   int		twoD;
{
   if (twoD) {
      *gx = g1*cos(b0) - g2*sin(b0);
      *gy = g1*sin(b0) + g2*cos(b0);
      *gz = 0.0;
   }
   else {
      *gx = (b0*b0 + b1*b1 - b2*b2 - b3*b3)*g1 + 
            2.0*(b1*b2 - b0*b3)*g2 +
            2.0*(b1*b3 + b0*b2)*g3;

      *gy = 2.0*(b1*b2 + b0*b3)*g1 +
            (b0*b0 - b1*b1 + b2*b2 - b3*b3)*g2 +
            2.0*(b2*b3 - b0*b1)*g3;

      *gz = 2.0*(b1*b3 - b0*b2)*g1 +
            2.0*(b2*b3 + b0*b1)*g2 +
            (b0*b0 - b1*b1 - b2*b2+ b3*b3)*g3;
   }
 
   return;
} 

void RotateToLocal (gx, gy, gz, b0, b1, b2, b3, g1, g2, g3, twoD)
   double	gx, gy, gz;
   double	b0, b1, b2, b3;
   double	*g1, *g2, *g3;
   int		twoD;
{
   if (twoD) {
      *g1 = gx*cos(b0) + gy*sin(b0);
      *g2 = -gx*sin(b0) + gy*cos(b0);
      *g3 = 0.0;   
   }
   else {
      *g1 = (b0*b0 + b1*b1 - b2*b2 - b3*b3)*gx + 
            2.0*(b1*b2 + b0*b3)*gy +
            2.0*(b1*b3 - b0*b2)*gz;

      *g2 = 2.0*(b1*b2 - b0*b3)*gx +
            (b0*b0 - b1*b1 + b2*b2 - b3*b3)*gy +
            2.0*(b2*b3 + b0*b1)*gz;

      *g3 = 2.0*(b1*b3 + b0*b2)*gx +
            2.0*(b2*b3 - b0*b1)*gy +
            (b0*b0 - b1*b1 - b2*b2+ b3*b3)*gz ;
   }
 
   return;
} 
