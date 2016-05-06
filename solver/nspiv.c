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
 * File:        nspiv.c
 *
 * Description: contains code that implements Sherman's NSPIV algorithm
 *	        for Gaussian elimination and solution of sparse linear 
 *		systems
 *
 ****************************************************************************/

# include <math.h>

int nspiv(n, ia, ja, a, b, max, r, c, ic, x, y, p, iu, ju, u, info)
   int		n;
   int		*ia, *ja;
   double	*a, *b, *x;
   int		max;
   int		*r, *c, *ic;
   double	*y, *u;
   int		*p, *iu, *ju;
   int		*info;
{
   double	dk, lki, xpv, xpvmax = 0.0, yk;
   int		ck, pk, ppk, pv, v, vi, vj, vk;
   int		i, j, k;
   int		juptr;
   int		jmin, jmax;
   int		maxc, maxcl;
   int		nzcnt;
   

   for (i = 1 ; i <= n ; i++)
      x [i] = 0.0;

   iu [1] = 1;
   juptr = 0;

   for (k = 1 ; k <= n ; k++) {
      p [n+1] = n + 1;
      vk = r [k]; 

      jmin = ia [vk];
      jmax = ia [vk + 1] - 1;

      if (jmin > jmax) {
         *info = k;
         return 2;
      }

      j = jmax;
      do {
         vj = ic [ja [j]];
         x [vj] = a [j];
         ppk = n + 1;

         do {
            pk = ppk;
            ppk = p [pk];
         } while (ppk - vj < 0);

         if (ppk - vj == 0) {
            *info = k;
            return 3;
         }

         p [vj] = ppk;
         p [pk] = vj;
 
         j --;
      } while (j >= jmin);

      vi = n + 1;
      yk = b [vk];

      while (1) {
         vi = p [vi];
         if (vi >= k)
            break;

         lki = -x [vi];
         x [vi] = 0.0;

         yk = yk + lki*y [vi];
         ppk = vi;
         jmin = iu [vi];
         jmax = iu [vi + 1] - 1;

         if (jmin > jmax)
            continue;

         for (j = jmin ; j <= jmax ; j++) {
            vj = ic [ju [j]];

            if (x [vj] != 0.0 || vj == ppk) {
               x [vj] += lki*u [j];
               continue;
            }

            if (vj < ppk)
               ppk = vi;

            do {
               pk = ppk; 
               ppk = p [pk];
            } while (ppk - vj < 0);

            if (ppk > vj) {
               p [vj] = ppk;
               p [pk] = vj;
               ppk = vj;
            }

            x [vj] += lki*u [j];
         }
      } 

      if (vi > n) {
         *info = k;
         return 4;
      }

      xpvmax = fabs(x [vi]);
      maxc = vi;
      nzcnt = 0;

      maxcl = vi;

      pv = vi;

      while (1) {
         v = pv;
         pv = p [pv];

         if (pv > n)
            break;

         nzcnt ++;
         xpv = fabs(x [pv]);

         if (xpv <= xpvmax)
            continue;

         xpvmax = xpv;
         maxc = pv;
         maxcl = v;
      } 

      if (xpvmax == 0.0) {
         *info = k;
         return 4;
      }

      if (vi == k || vi == maxc) 
         vi = p [vi];
      else
         p [maxcl] = p [maxc];

      dk = 1.0/x [maxc];

      x [maxc] = x [k];

      i = c [k];
      c [k] = c [maxc];
      c [maxc] = i;

      ck = c [k];
      ic [ck] = k;
      ic [i] = maxc;
      x [k] = 0.0;

      y [k] = yk * dk;

      iu [k+1] = iu [k] + nzcnt;

      if (iu [k+1] > max + 1) {
         *info = k;
         return 5;
      }

      if (vi <= n) {
         j = vi;
         do {
            juptr ++;
            ju [juptr] = c [j];
            u [juptr] = x [j] * dk;
            x [j] = 0.0;
            j = p [j];
         } while (j <= n);
      }
   }
   k = n;

   for (i = 1 ; i <= n ; i++) {
      yk = y [k];
      jmin = iu [k];
      jmax = iu [k+1] - 1;

      if (jmin <= jmax) 
         for (j = jmin ; j <= jmax ; j++) 
            yk = yk - u [j]*y [ic [ju [j]]];

      y [k] = yk;
      ck = c [k];
      x [ck] = yk;
      k --;
   }
         
   *info = 0;
   return 0;
}
