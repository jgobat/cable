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
 * File:        buoyancy.c
 *
 * Description: routines for evaluating the buoyancy of a buoy object
 *
 * History:
 *
 ****************************************************************************/

# include <stdio.h>
# include <math.h>
# include "compress.h"
# include "problem.h"

# define SQR(a) ((a)*(a))


double Buoyancy (b, draft, e)
   Buoy		b;		/* the buoy in question			   */
   double	draft;
   struct environment *e;
{
   double	rho_g;
   double	r;
   double	theta;
   double	A;
   double	r1, r2, h;
   double	V;
   int		i;
   int		waterline;
   double	m;

   rho_g = e -> rho*e -> gravity;

   if (draft < 0.0) 
      return 0.0;

   if (draft > b -> max_draft)
      return b -> buoyancy;

   switch (b -> type) {
      case Cylinder:
         return rho_g*draft*0.25*M_PI*SQR(b -> d);
 
         break;

      case Sphere:
         return rho_g*(1.5*b -> d - draft)*SQR(draft)*M_PI/3.0;
        
         break;

      case Capsule:	
          r = 0.5*b -> d;
          theta = 2.0*acos((r  - draft)/r);
          A =  0.5*r*r*(theta - sin(theta));
          return  rho_g*(1.5*b -> d - draft)*SQR(draft)*M_PI/3.0
                  + rho_g*A*(b -> h - b -> d);

      case Axisymmetric:
          V = 0.0;
          waterline = 0;

          for (i = 1 ; i < b -> num_diameters ; i++) {
             r1 = 0.5*b -> diameters [i].d;
             if (b -> diameters [i + 1].level > draft) {
                h = draft - b -> diameters [i].level;
                m = (b -> diameters [i + 1].d - b -> diameters [i].d) /
                    (b -> diameters [i + 1].level - b -> diameters [i].level); 
                
                r2 = 0.5*(b -> diameters [i].d + m*h);

                waterline = 1;
             }
             else {
                r2 = 0.5*b -> diameters [i + 1].d;
                h = b -> diameters [i + 1].level - b -> diameters [i].level;
             }
       
             V += M_PI/3.0*h*(r1*r1 + r1*r2 + r2*r2);
             if (waterline)
                break;
          }

          return rho_g*V;

      default:
         return 0.0;
   }
         
   return 0.0;
}

double Draft (buoyancy, b, e)
   double	buoyancy;	/* desired buoyancy			   */
   Buoy		b;		/* the buoy in question			   */
   struct environment *e;
{
   double	h1, h2, h3;
   double	B1, B3;
   double	rho_g;

   if (buoyancy > b -> buoyancy)
      return -1.0;

   if (buoyancy < 0.0)
      return 0.0;

   rho_g = e -> rho*e -> gravity;

   switch (b -> type) {
      case Cylinder:
         return buoyancy / b -> buoyancy*b -> h;
 
         break;

      case Sphere:
      case Capsule:
      case Axisymmetric:
         h1 = 0.0;
         h2 = b -> max_draft;
         while (1) {
            h3 = 0.5*(h1 + h2); 
            B3 = Buoyancy(b, h3, e) - buoyancy;
            B1 = Buoyancy(b, h1, e) - buoyancy;
            if (B3*B1 < 0.0)
               h2 = h3;
            else
               h1 = h3;

            if (fabs((h1 - h2) / b -> max_draft) < 1e-8 || B3 == 0.0) 
               return 0.5*(h1 + h2);  
         }
         
         break;

      default:
         return 0.0;
   }
         
   return 0.0;
}

double ProjectedArea (b, draft)
   Buoy		b;		/* the buoy in question			   */
   double	draft;
{
   double	r;
   double	theta;
   double	d1, d2;
   int		waterline;
   double	h;
   double	A;
   double	m;
   int		i;
 
   if (draft < 0.0) 
      return 0.0;

   if (draft > b -> max_draft)
      return b -> S;

   switch (b -> type) {
      case Cylinder:
         return draft * b -> d;
 
         break;

      case Sphere:
      case Capsule:
         r = 0.5*b -> d;
         theta = 2.0*acos((r  - draft)/r);
         return 0.5*r*r*(theta - sin(theta));

      case Axisymmetric:
          A = 0.0;
          waterline = 0;

          for (i = 1 ; i < b -> num_diameters ; i++) {
             d1 = b -> diameters [i].d;
             if (b -> diameters [i + 1].level > draft) {
                h = draft - b -> diameters [i].level;
                m = (b -> diameters [i + 1].d - b -> diameters [i].d) /
                    (b -> diameters [i + 1].level - b -> diameters [i].level); 
                d2 = b -> diameters [i].d + m*h;

                waterline = 1;
             }
             else {
                d2 = b -> diameters [i + 1].d;
                h = b -> diameters [i + 1].level - b -> diameters [i].level;
             }
        
             A += 0.5*h*(d1 + d2);
             if (waterline)
                break;
         }
         
         return A;
 
         break;

      default:
         return 0.0;
   }
         
   return 0.0;
}

double WaterplaneArea (b, draft)
   Buoy	   b;
   double  draft;
{
   return 0.25*M_PI*(b -> d)*(b -> d);
}

