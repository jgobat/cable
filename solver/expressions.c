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
 * File:        expressions.c
 *
 * Description: routines for evaluating variable expressions
 *
 * History:
 *
 ****************************************************************************/

# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <string.h>
# include "compress.h"
# include "problem.h"
# include "code.h"
# include "error.h"

extern Problem *problem;             
extern Environment *environment;

void Speed (t, term, u, v, w)
   double	 t;
   Terminal	 term;
   double	*u;
   double	*v;
   double	*w;
{
   Node      n;

   n = term -> node;

   if (term -> xspeed.expr)
      *u = EvalCode(term -> xspeed.expr, n, t, t, 0, 0, 0, CURRNODEDATA);
   else
      *u = term -> xspeed.value;

   if (term -> yspeed.expr)
      *v = EvalCode(term -> yspeed.expr, n, t, t, 0, 0, 0, CURRNODEDATA);
   else
      *v = term -> yspeed.value;

   if (w) {
      if (term -> zspeed.expr)
         *w = EvalCode(term -> zspeed.expr, n, t, t, 0, 0, 0, CURRNODEDATA);
 
      else
         *w = term -> zspeed.value;
   }
}

void ConnectorThrust (t, seg, node, Fx, Fy, Fz)
   double	t;
   Segment  seg;
   Node     node;
   double      *Fx;
   double      *Fy;
   double      *Fz;
{
   if (seg -> connector_xthrust.expr)
      *Fx = EvalCode(seg -> connector_xthrust.expr, 
                     node, t, t, 0, 0, 0, CURRNODEDATA);

   else
      *Fx = seg -> connector_xthrust.value;

   if (seg -> connector_ythrust.expr)
      *Fy = EvalCode(seg -> connector_ythrust.expr, 
                     node, t, t, 0, 0, 0, CURRNODEDATA);

   else
      *Fy = seg -> connector_ythrust.value;

   if (Fz) {
      if (seg -> connector_zthrust.expr)
         *Fz = EvalCode(seg -> connector_zthrust.expr, 
                        node, t, t, 0, 0, 0, CURRNODEDATA);

      else
         *Fz = seg -> connector_zthrust.value;
   }
}

void Thrust (t, term, Fx, Fy, Fz)
   double	t;
   Terminal     term;
   double      *Fx;
   double      *Fy;
   double      *Fz;
{
   Node node;

   node = term -> node;

   if (term -> xthrust.expr)
      *Fx = EvalCode(term -> xthrust.expr, node, t, t, 0, 0, 0, CURRNODEDATA);

   else
      *Fx = term -> xthrust.value;

   if (term -> ythrust.expr)
      *Fy = EvalCode(term -> ythrust.expr, node, t, t, 0, 0, 0, CURRNODEDATA);

   else
      *Fy = term -> ythrust.value;

   if (Fz) {
      if (term -> zthrust.expr)
         *Fz = EvalCode(term -> zthrust.expr, node, t, t, 0, 0, 0, CURRNODEDATA);

      else
         *Fz = term -> zthrust.value;
   }
}

static double
triterp(double *** U, 
        int ix,
        int iy,
        int iz,
        int ix1,
        int iy1,
        int iz1,
        double xf, 
        double yf, 
        double zf)
{
   double c00, c01, c10, c11, c0, c1, c;

   c00 = U[iz][iy][ix]*(1.0 - zf) + U[iz1][iy][ix]*zf;
   c10 = U[iz][iy1][ix]*(1.0 - zf) + U[iz1][iy1][ix]*zf;
   c01 = U[iz][iy][ix1]*(1.0 - zf) + U[iz1][iy][ix1]*zf;
   c11 = U[iz][iy1][ix1]*(1.0 - zf) + U[iz1][iy1][ix1]*zf;

   c0 = c00*(1.0 - yf) + c10*yf;
   c1 = c01*(1.0 - yf) + c11*yf;

   c = c0*(1.0 - xf) + c1*xf;

   return c;
}

static void CurrentFromFile(double tm, 
                            double x /* depth */, 
                            double y, 
                            double z, 
                            double *u, 
                            double *v, 
                            double *w)
{
   FILE		      *fp;
   static double   mint, minx, miny, minz;
   static double   maxt, maxx, maxy, maxz;
   static double   dz = 0, dt = 0, dx = 0, dy = 0;
   static int	   nz = 0, nt = 0, nx = 0, ny = 0;
   int             n;
   int             it, ix, iy, iz;
   int             it1, ix1, iy1, iz1;
   double	      tf, xf, yf, zf;
   double         ca, cb;
   static double ****U = NULL;
   static double ****V = NULL;
   static double ****W = NULL;
 
   *u = 0;
   *v = 0;
   *w = 0;
   
   if (!nz) {

      fp = fopen(environment -> current_file, "r");
      if (fp == NULL) {
         error("could not open current-file \"%s\"", environment -> current_file);
         return;
      }

      n = fscanf(fp, "%lf %lf %d %lf %lf %d %lf %lf %d %lf %lf %d", 
                 &dt, &mint, &nt, 
                 &dx, &minx, &nx, 
                 &dy, &miny, &ny, 
                 &dz, &minz, &nz);
      if (n != 12 
          || nt == 0 || nx == 0 || ny == 0 || nz == 0 
          || dt == 0 || dx == 0 || dy == 0 || dz == 0) {
         error("error parsing header line - should be 'dt mint nt dz minz nz dx minx nx dy miny ny'");
         error("dt, dx, dy, dz and nt, nx, ny, nz must all be non-zero");
         fclose(fp); 
         return;
      }

      maxt = mint + dt*(nt - 1);
      maxx = minx + dx*(nx - 1);
      maxy = miny + dy*(ny - 1);
      maxz = minz + dz*(nz - 1);

      U = (double ****) malloc(sizeof(double ***) * nt); 
      V = (double ****) malloc(sizeof(double ***) * nt); 
      W = (double ****) malloc(sizeof(double ***) * nt); 
 
      for (it = 0 ; it < nt ; it ++) {
         U[it] = (double ***) malloc(sizeof(double**) * nz); 
         V[it] = (double ***) malloc(sizeof(double**) * nz); 
         W[it] = (double ***) malloc(sizeof(double**) * nz); 

         for (iz = 0 ; iz < nz ; iz ++) {
            U[it][iz] = (double **) malloc(sizeof(double*) * ny); 
            V[it][iz] = (double **) malloc(sizeof(double*) * ny); 
            W[it][iz] = (double **) malloc(sizeof(double*) * ny); 

            for (iy = 0 ; iy < ny ; iy ++) {
               U[it][iz][iy] = (double *) malloc(sizeof(double) * nx); 
               V[it][iz][iy] = (double *) malloc(sizeof(double) * nx); 
               W[it][iz][iy] = (double *) malloc(sizeof(double) * nx); 
 
               for (ix = 0 ; ix < nx ; ix++) {
                   n = fscanf(fp, "%lf %lf %lf", &(V[it][iz][iy][ix]), &(W[it][iz][iy][ix]), &(U[it][iz][iy][ix]));
                   if (n != 3) {
                      error("error parsing current table, %d %d %d %d", it, iz, iy, ix);
                      fclose(fp);
                      return;
                   }
               }
            }
         }
      }
   }

   if (tm > maxt)
      tm = maxt;
   else if (tm < mint)
      tm = mint;

   if (x > maxx)
      x = maxx;
   else if (x < minx)
      x = minx;

   if (y > maxy)
      y = maxy;
   else if (y < miny)
      y = miny;

   if (z > maxz)
      z = maxz;
   else if (z < minz)
      z = minz;

   it = (int) ((tm - mint) / dt);
   ix = (int) ((x - minx) / dx);
   iy = (int) ((y - miny) / dy);
   iz = (int) ((z - minz) / dz);

   if (it == nt - 1)
      it1 = it;
   else
      it1 = it + 1;

   if (ix == nx - 1)
      ix1 = ix;
   else
      ix1 = ix + 1;

   if (iy == ny - 1)
      iy1 = iy;
   else
      iy1 = iy + 1;

   if (iz == nz - 1)
      iz1 = iz;
   else
      iz1 = iz + 1;

   tf = ((tm - mint) - it*dt) / dt;
   xf = ((x - minx) - ix*dx) / dx;
   yf = ((y - miny) - iy*dy) / dy;
   zf = ((z - minz) - iz*dz) / dz;

   ca = triterp(U[it], ix, iy, iz, ix1, iy1, iz1, xf, yf, zf);
   cb = triterp(U[it1], ix, iy, iz, ix1, iy1, iz1, xf, yf, zf);
   *u = ca + (cb - ca)*tf;

   ca = triterp(V[it], ix, iy, iz, ix1, iy1, iz1, xf, yf, zf);
   cb = triterp(V[it1], ix, iy, iz, ix1, iy1, iz1, xf, yf, zf);
   *v = ca + (cb - ca)*tf;

   if (w) {
      ca = triterp(W[it], ix, iy, iz, ix1, iy1, iz1, xf, yf, zf);
      cb = triterp(W[it1], ix, iy, iz, ix1, iy1, iz1, xf, yf, zf);
      *w = ca + (cb - ca)*tf;
   }
   

   return;
}

void 
Current (t, x, y, z, u, v, w)
   double	t;
   double   x; 
   double   y; 
   double   z; 
   double  *u;
   double  *v;
   double  *w;
{
   double   vrot, wrot;
   double   th;
   double	depth;
   double	modulation;
   double	xscale, yscale, zscale;

   depth = environment -> surface - x;

   if (environment -> depth && x < 0.0) {
      *v = 0.0;
      if (w) *w = 0.0;
      *u = 0.0;

      return;
   }

   if (depth < 0.0)
      depth = 0.0;

   if (environment -> Uxscale != 1)
      xscale = environment -> Uxscale;
   else
      xscale = environment -> Uscale;

   if (environment -> Uyscale != 1)
      yscale = environment -> Uyscale;
   else
      yscale = environment -> Uscale;

   if (environment -> Uzscale != 1)
      zscale = environment -> Uzscale;
   else
      zscale = environment -> Uscale;

   if (environment -> current_file)  {
      CurrentFromFile(t, depth, y, z, u, v, w);
      *u *= xscale;
      *v *= yscale;
      if (w) *w *= zscale;

      return;
   }

   if (w) {
      if (environment -> Uz_mod.expr)
          modulation = zscale*EvalCode(environment -> Uz_mod.expr, 
                                       NULL, t, t, x, y, z, CURRNODEDATA);

      else
          modulation = zscale*environment -> Uz_mod.value;

      if (environment -> Uz.expr)
         *w = modulation*EvalCode(environment -> Uz.expr, 
                                  NULL, t, depth, x, y, z, CURRNODEDATA);

      else
         *w = modulation*environment -> Uz.value;
   }

   if (environment -> Uy_mod.expr)
       modulation = yscale*EvalCode(environment -> Uy_mod.expr, 
                                    NULL, t, t, x, y, z, CURRNODEDATA);
   else
       modulation = yscale*environment -> Uy_mod.value;

   // modulation = 0.6;

   if (environment -> Uy.expr)
      *v = modulation*EvalCode(environment -> Uy.expr, 
                               NULL, t, depth, x, y, z, CURRNODEDATA);
   else
      *v = modulation*environment -> Uy.value;

   if (environment -> Ux_mod.expr)
       modulation = xscale*EvalCode(environment -> Ux_mod.expr, 
                                    NULL, t, t, x, y, z, CURRNODEDATA);
   else
       modulation = xscale*environment -> Ux_mod.value;

   if (environment -> Ux.expr)
      *u = modulation*EvalCode(environment -> Ux.expr, 
                               NULL, t, depth, x, y, z, CURRNODEDATA);
   else
      *u = modulation*environment -> Ux.value;

   // rotate the horizontal components if requested
   if ((th = environment -> current_rotation)) {
      vrot = (*v)*cos(th) - (*w)*sin(th);
      wrot = (*v)*sin(th) + (*w)*cos(th);
      *v = vrot;
      *w = wrot;
      // printf("rot: 1=%f, 2=%f\n", *v, *w);
   }
   //else 
   //   printf("no rot: 1=%f, 2=%f\n", *v, *w);


   return;
}

void 
DragCoeff(Node node, Material mat, double t, double u, double v, double *Cdt, double *Cdn)
{
   static double	t_fac, n_fac = 0;

   if (!n_fac) {
      n_fac = 0.5*environment -> rho;
      t_fac = n_fac*M_PI;
   }

   if (!mat -> Cdt.expr) 
      *Cdt = mat -> Cdt.value; 
   else {
      *Cdt = EvalCode(mat -> Cdt.expr, NULL, 
                      t, u, u, v, 0, CURRNODEDATA)*t_fac*mat -> d;
   }
 
   if (!node -> drag) {
      if (!mat -> Cdn.expr)
         *Cdn = mat -> Cdn.value;
      else
         *Cdn = n_fac*mat -> d*EvalCode(mat -> Cdn.expr, NULL, 
                                        t, v, u, v, 0, CURRNODEDATA);
   }
   else {
      *Cdn = node -> drag*n_fac*mat -> d;
   }
 
   return;       
}

void WindDrag (t, b, Fy, Fz)
   double  t;
   Buoy	   b;
   double *Fy;
   double *Fz;
{
   Node          node;
   static double air_rho = 0.0;
   double		 Cd_factor;
   double		 A;
   double	     V;
   double		 W;

   node = problem -> terminal[2] -> node;

   if (air_rho == 0.0)  
      air_rho = environment -> rho * 1.26/1025.9;

   if (b -> Sw)
      A = b -> Sw;
   else 
      A = ProjectedArea(b, b -> max_draft) - ProjectedArea(b, b -> draft);

   if (A <= 0.0) {
      *Fy = 0.0;
      if (Fz)
         *Fz = 0.0;

      return;
   }

   Cd_factor = 0.5*air_rho*A*b -> Cdw;

   if (environment -> y_wind.expr)
      V = EvalCode(environment -> y_wind.expr, node, t, t, 0, 0, 0, CURRNODEDATA);
   else
      V = environment -> y_wind.value;

   *Fy = Cd_factor*V*fabs(V);

   if (Fz) {
      if (environment -> z_wind.expr)
         W = EvalCode(environment -> z_wind.expr, node, t, t, 0, 0, 0, CURRNODEDATA);
      else
         W = environment -> z_wind.value;

      *Fz = Cd_factor*W*fabs(W);
   }
  
   return; 
}

double 
Bottom(double y, double z, double t)
{
   if (environment -> bottom_elevation.expr)
      return EvalCode(environment -> bottom_elevation.expr, 
                      NULL, t, y, 0, y, z, CURRNODEDATA);
   else
      return environment -> bottom_elevation.value;	/* this wouldn't  */
							/* make sense for */								/* any value <> 0 */
}

static double 
LinearParabolicProfile(double y)
{
   static int		init = 0;
   static double	L;
   static double	L1;
   static double	A;
   static double	z1;
   static double	H;
   static double	m;
   double		z;
   double		xl;

   if (!init) {
      z1 = problem -> terminal [1] -> profile_turn;
      m = problem -> terminal [1] -> profile_m;
      H = -problem -> terminal [1] -> profile_H;
      A = m*m/4.0/z1;
      L1 = 4.0*z1/m;

      L = L1 + (-H - z1*2)/m;

      init = 1;
   }

   xl = fmod(y, 2.0*L);
   
   if (xl < L) {

      if (xl < L1/2)
         z = -A*xl*xl;
      else if (xl > L - L1/2)
         z = H + A*(xl - L)*(xl - L);
      else
         z = -z1 - m*(xl - L1/2);
   }
   else if (xl >= L) {

      if (xl < L + L1/2) 
         z = H + A*(xl - L)*(xl - L);
      else if (xl > 2.0*L - L1/2)
         z = -A*(xl - 2*L)*(xl - 2*L);
      else
         z = (H + z1) + m*(xl - L - L1/2);
   }
   else
      z = 0;

   return -z;
} 

double 
Profile(double t)
{
   Node node;

   node = problem -> terminal[1] -> node;

   if (problem -> terminal [1] -> profile.expr)
      return EvalCode(problem -> terminal [1] -> profile.expr, 
                      node, t, node -> y, 0, 0, 0, CURRNODEDATA);
   else if (problem -> terminal [1] -> profile_m && problem -> terminal [1] -> profile_H && problem -> terminal [1] -> profile_turn)
      return LinearParabolicProfile(node -> y);
   else
      return problem -> terminal [1] -> profile.value; 
}
