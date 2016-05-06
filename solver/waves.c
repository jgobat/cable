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
 * File:        waves.c
 *
 * Description: routines to calculate the dynamic inputs for the various
 *		algorithms - surface velocity, particle velocity, Froude-Krylov
 *		forces, etc. 
 *
 * History:
 *
 ****************************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <malloc.h>
# include "compress.h"
# include "allocate.h"
# include "problem.h"
# include "solve.h"
# include "error.h"

extern Problem     *problem;             
extern Environment *environment;
extern Analysis    *analysis;

# ifdef HAVELAMP
extern void pvmoor_(float*,float*,float*);
extern void vfmoor_(float*,float*,float*);
extern void afmoor_(float*,float*,float*);

extern void TranslateLAMP(double, double, double, float *);
# endif

# define NFREQ 200
# define MAXFREQ 1000

static double   d_om = 0.05;	/* was 0.2 for NFREQ = 50 */
static double   om_0 = 0.2;

/*
static double   d_om = 0.1;
static double   om_0 = 0.01;
*/

static double   rand_a [MAXFREQ][4];
static double   rand_k [MAXFREQ];
static double   rand_phi [MAXFREQ];

static int 	nfreq;

static int	initialized_file = 0;
static int	initialized_random = 0;

static double   phi0 [] = { 1.69316, -2.8771, 0.983096, -1.21617, 
                            -1.97345, -0.159904, -1.08959, -1.32703, 
                            -0.027736, 0.506554, -1.88744, 0.899007, 
                            1.16416, 1.24696, -0.893483, 0.166244, 
                            0.692317, -2.02269, -1.2134, 1.91341, 
                            -1.2149, 2.25347, -0.76603, 1.01154, 
                            -3.13839, -0.615954, -1.22978, 0.469424, 
                            1.82322, -0.919958, -2.11945, -1.3447, 
                            -2.64821, -0.844951, 1.74941, 2.64124, 
                            -0.817751, -1.59875, 2.70081, 0.147633, 
                            2.9987, -2.81528, -1.23859, 1.17337, 
                            1.96069, 0.666586, -1.50923, -2.09592, 
                            1.56957, 2.80025, 1.77042, 1.30903, 
                            0.322802, 1.12641, 2.74825, 0.933896, 
                            2.23197, 2.13716, -3.03776, -2.45981, 
                            -1.4215, 0.21948, -2.47758, 0.288545, 
                            -0.638365, -3.04535, 2.0939, 1.53917, 
                            1.73716, -1.69458, -0.54061, -1.85774, 
                            2.35641, 2.46192, 2.39476, 0.986263, 
                            1.22726, -0.0429505, 0.647653, -1.8483, 
                            1.18727, -1.41551, 1.87553, 1.36981, 
                            0.0964465, 2.60829, 2.14514, 2.66289, 
                            1.42943, -2.19229, -0.0606957, 2.73811, 
                            2.03007, -0.26052, 1.83975, -2.91452, 
                            2.42538, -2.71877, 0.220767, -2.61116, 
                            -2.50973, -1.42307, -0.951817, -2.59167, 
                            2.45965, -0.60258, 0.61414, -1.8815, 
                            -1.79629, -0.532984, -2.791, -0.959463, 
                            0.569178, -0.36951, 2.94765, -0.173411, 
                            -0.187082, -2.56604, 0.874993, -0.464312, 
                            0.532696, 0.484642, -3.08891, 1.009, 
                            -1.33665, 2.44978, 2.45378, -1.49976, 
                            -2.52419, 2.9132, 1.79793, 1.84775, 
                            -0.66208, 0.931712, 3.13308, -0.227562, 
                            -2.83883, 0.710332, -0.966687, 1.41606, 
                            2.42773, 2.91031, -1.06106, -0.597485, 
                            -0.937891, -1.90212, 1.93653, 1.64567, 
                            -0.812116, 1.78993, 3.06123, -1.49177, 
                            2.11328, -1.60414, 0.801184, -1.18862, 
                            0.181638, -0.0815412, 2.85964, -1.96816, 
                            -0.188591, 2.7189, -1.40612, -2.11166, 
                            -1.68712, 1.06854, -1.88436, 0.745122, 
                            0.245877, -2.49089, -0.750716, 1.38101, 
                            -0.673653, -1.54699, 0.468333, 0.730207, 
                            -0.398784, 2.00976, 3.10273, 1.42297, 
                            -0.7391, -2.37902, -3.06496, -1.24698, 
                            0.609394, 2.36221, 0.286012, -0.0888212, 
                            2.22244, 1.96596, 0.293549, 0.276843, 
                            1.9336, 2.37203, -0.647477, 1.16581, 
                            1.62915, -2.07146, 2.63449, 1.01329 };

/* used these when NFREQ = 50
static double   phi0 [] = {2.30698, 0.733712, -1.96964, -2.13934,
                            0.161887, -0.0342021, 2.84504, -1.66304,
                            3.01803, 0.52101, -1.82689, -0.677369,
                           -0.782422, 0.450329, 1.6233, -0.152927,
                           -1.93811, -0.334018, -2.76687, 2.0241,
                           -1.61191, -2.32674, 1.61207, -1.00234,
                           -0.378339, -0.463647, -1.22827, 0.88083,
                           -2.35563, -2.65927, 0.0806407, -2.23175,
                           2.89801, -0.0250725, -1.11447, -0.58421,
                           0.433688, -2.68729, -0.57556, 0.179637,
                           -2.49879, -1.28793, 0.953925, 0.303592,
                           1.20579, -2.16877, -1.08928, -0.591202,
                           -1.66428, -2.26863};
*/


static double *Uf;
static double *Vf;
static double *Wf;
static double *tf;
static double  dt_f;

static void InitializeVelocityFile ( )
{
   FILE		*fp;
   int		 i, n;
   int		 max;

   max = 5000;
   Uf = (double *) malloc(sizeof(double) * max); Uf --; 
   Vf = (double *) malloc(sizeof(double) * max); Vf --; 
   Wf = (double *) malloc(sizeof(double) * max); Wf --; 
   tf = (double *) malloc(sizeof(double) * max); tf --; 

   fp = fopen(environment -> velocity_file, "r");
   if (fp == NULL)
      ExitErr ("could not open file %s for reading", environment -> velocity_file);

   i = 0;
   while(!feof(fp)) {
      i ++;
      n = fscanf(fp, "%lf %lf %lf %lf", &(tf [i]), &(Vf [i]), &(Wf [i]), &(Uf [i]));
      if (n != 4)
         break;      

      if (i >= max) {
         max += 5000;

         Uf ++; Uf = (double *) realloc((void *) Uf, sizeof(double) * max); Uf --;
         Vf ++; Vf = (double *) realloc((void *) Vf, sizeof(double) * max); Vf --;
         Wf ++; Wf = (double *) realloc((void *) Wf, sizeof(double) * max); Wf --;
         tf ++; tf = (double *) realloc((void *) tf, sizeof(double) * max); tf --;
      }
   }

   i --;

   fclose (fp);

   if (tf [i] < analysis -> duration) 
      ExitErr ("not enough wave input to cover full simulation (%g)", tf [i]);

   dt_f = tf [2] - tf [1];

   initialized_file = 1;

   return;
}

static void InitializeWaveFile ( )
{
   FILE	 *fp;
   int	  i;
   double om;
   double S;
   double phi;
   int    idx;


   fp = fopen(environment -> wave_file, "r");
   if (!fp) 
      ExitErr ("could not open file %s for reading", environment -> wave_file);

   if (environment -> forcing == Velocity)
      idx = 1;
   else 
      idx = 2;

   environment -> amplitude [1][0] = 0;
   environment -> amplitude [2][0] = 0;
   environment -> amplitude [3][0] = 0;

   environment -> amplitude [idx][0] = 1;

   i = 0;
   while (!feof(fp)) {
      fscanf(fp, "%lf %lf %lf", &om, &S, &phi);

      if (i == 0) 
         om_0 = om;
      else if (i == 1)
         d_om = om - om_0;

      rand_a [i][idx] = sqrt(2.0*S*d_om);
      rand_k [i] = dispersion(om, environment -> gravity, environment -> depth);
      rand_phi [i] = phi;

      i ++;
   }


   fclose(fp);

   nfreq = i - 1;

   return;
}

static void InitializeRandom ( )
{
   double  om;
   double  om_m, H2, S;
   int     i, j;

   if (environment -> wave_file) {
      InitializeWaveFile ( );
      initialized_random = 1;
      return;
   }
    
   nfreq = NFREQ;
 
   om = om_0;
   for (i = 0 ; i < NFREQ ; i++)  {
      rand_k [i] = dispersion(om, environment -> gravity, environment -> depth);
      rand_phi [i] = phi0 [i];
      om += d_om;
   }

   for (j = 1 ; j <= 3 ; j++) {
      if (environment -> amplitude [j][0]) {
         om_m = environment -> omega [j][0];
         H2 = 4.0*environment -> amplitude [j][0]*environment -> amplitude [j][0];

         om = om_0;
         for (i = 0 ; i < NFREQ ; i++) {
            S = 0.3125*pow(om_m, 4.0)/
                       pow(om, 5.0)*H2*exp(-1.25*pow(om_m/om, 4.0));
            rand_a [i][j] = sqrt(2.0*S*d_om);

            om += d_om;
         }
      }
   }

   initialized_random = 1;
}

void WaveParticleMotion(t, x, y, z, u, v, w, ud, vd, wd)
   double	 t;
   double	 x; 
   double	 y; 
   double	 z;
   double	*u;
   double	*v;
   double	*w;
   double	*ud;
   double	*vd;
   double	*wd;
{
#ifdef HAVELAMP
   float	lamp_t;
   float	rpos [3];
   float	vlamp [3];
   float	alamp [3];

   if (environment -> forcing== LAMP) {
      
      lamp_t = t;

      TranslateLAMP(x, y, z, rpos);
    
      vfmoor_(&lamp_t, rpos, vlamp);
      *u = vlamp [2];
      *v = vlamp [0];
      if (w)
         *w = vlamp [1];

      afmoor_(&lamp_t, rpos, alamp); 
      *ud = alamp [2];
      *vd = alamp [0];
      if (wd)
         *wd = alamp [1];

      return;
   }
#endif

   return;
}

void 
WaveStokesVelocity (x, v, w)
   double	 x; 
   double	*v;
   double	*w;
{
   double  ampl;
   double  attn;
   int	   i;
   double  om;
   double  k;
   double  skH;
   double  ckx;

   *v = 0.0;

   if (w)
      *w = 0.0;

   if (x > environment -> depth || x < 0) {
      *v = 0.0;
      if (w) *w = 0.0;
      return;
   }

   if (environment -> input_type == Regular) {
    
      for (i = 0 ; i < environment -> num_components [2] ; i++) {
         k = environment -> k [2][i];
         skH = sinh(k*environment -> depth);
          
         if (!isinf(skH)) {
            om = environment -> omega [2][i];
            ampl = environment -> amplitude [2][i]*environment -> amplitude[2][i]*k*om;

            ckx = cosh(2.0*k*x);

            attn = ckx/2.0/skH/skH;
            *v += ampl*attn;
            printf("ampl=%f, attn=%f, om=%f, k=%f\n", ampl, attn, om, k);
         }
      }

      if (w) {
         for (i = 0 ; i < environment -> num_components [3] ; i++) {
            k = environment -> k [3][i];
            skH = sinh(k*environment -> depth);
 
            if (!isinf(skH)) {
               om = environment -> omega [3][i];
               ampl = environment -> amplitude [3][i]*environment -> amplitude[3][i]*k*om;

               ckx = cosh(2.0*k*x);
     
               attn = ckx / 2.0 / skH / skH;
               *w += ampl*attn;
            }
         }
      }
   }
   else if (environment -> input_type == Random) {
      if (!initialized_random) 
         InitializeRandom ( );

      if (environment -> amplitude [2][0]) {
         om = om_0;

         for (i = 0 ; i < nfreq ; i++) {
            k = rand_k [i];
            if (environment -> depth - x > M_PI/k)
                break;

            if (environment -> depth < M_PI/k) {
               skH = sinh(k*environment -> depth);
 
               if (isinf(skH)) 
                  continue;
 
               ckx = cosh(2.0*k*x);
               attn = ckx / 2.0 / skH / skH;
            }
            else {
                attn = exp(2.0*k*(x - environment -> depth));    
            }
            ampl = rand_a [i][2]*rand_a[i][2]*k*om;  

            *v += ampl*attn;

            printf("ampl=%f, attn=%f, om=%f, k=%f\n", ampl, attn, om, k);

            om += d_om;
         }
      }
      if (environment -> amplitude [3][0] && w) {
         om = om_0;

         for (i = 0 ; i < nfreq ; i++) {
            k = rand_k [i];
            if (environment -> depth - x > M_PI/k)
                break;

            if (environment -> depth < M_PI/k) {
               skH = sinh(k*environment -> depth);
 
               if (isinf(skH)) 
                  continue;
 
               ckx = cosh(2.0*k*x);

               attn = ckx / 2.0 / skH / skH;
            }
            else {
                attn = exp(2.0*k*(x - environment -> depth));    
            }
            ampl = rand_a [i][3]*rand_a[i][3]*k*om;  
            *w += ampl*attn;

            om += d_om;
         }
      }
   }

   return;
}

void WaveParticleVelocity (t, x, y, z, u, v, w)
   double	 t;
   double	 x; 
   double	 y; 
   double	 z;
   double	*u;
   double	*v;
   double	*w;
{
   double  ampl;
   double  phas;
   double  attn;
   int	   i;
   double  factor;
   double  om;
   double  k;
   double  g;
   double  ckH;
   double  skx, ckx;

#ifdef HAVELAMP
   float   lamp_t;
   float   rpos [3];
   float   vlamp [3];

   if (environment -> forcing == LAMP) {
      *u = *v = 0.0;
      if (w)
         *w = 0.0;

      return;

      lamp_t = t;

      TranslateLAMP(x, y, z, rpos);
    
      vfmoor_(&lamp_t, rpos, vlamp);
      *u = vlamp [2];
      *v = vlamp [0];
      if (w)
         *w = vlamp [1];

      return;
   }
#endif
   
   if (analysis -> ramp_time && t < analysis -> ramp_time) 
      factor = t / analysis -> ramp_time;
   else 
      factor = 1.0;

   *u = 0.0;
   *v = 0.0;

   if (w)
      *w = 0.0;

   if (problem -> dynstat)
      return;

   if (x > environment -> depth || x < 0) {
      *u = *v = 0.0;
      if (w) *w = 0.0;
      return;
   }

   g = environment -> gravity;

   if (environment -> input_type == Regular) {
    
      for (i = 0 ; i < environment -> num_components [2] ; i++) {
         k = environment -> k [2][i];
         ckH = cosh(k*environment -> depth);
          
         if (!isinf(ckH)) {
            om = environment -> omega [2][i];
            ampl = factor*g*environment -> amplitude [2][i]*k/om;

	    skx = sinh(k*x);
	    ckx = cosh(k*x);

            phas = k*y - om*t - environment -> phase [2][i];
     
            attn = skx / ckH;
            *u += ampl*attn*sin(phas);

            attn = ckx / ckH;
            *v += ampl*attn*cos(phas);
         }
      }

      if (w) {
         for (i = 0 ; i < environment -> num_components [3] ; i++) {
            k = environment -> k [3][i];
            ckH = cosh(k*environment -> depth);
 
            if (!isinf(ckH)) {
               om = environment -> omega [3][i];
               ampl = factor*g*environment -> amplitude [3][i]*k/om;

               skx = sinh(k*x);
               ckx = cosh(k*x);

               phas = k*z - om*t - environment -> phase [3][i];
     
               attn = skx / ckH;
               *u += ampl*attn*sin(phas);

               attn = ckx / ckH;
               *w += ampl*attn*cos(phas);
            }
         }
      }
   }
   else if (environment -> input_type == Random) {
      if (!initialized_random) 
         InitializeRandom ( );

      if (environment -> amplitude [2][0]) {
         om = om_0;

         for (i = 0 ; i < nfreq ; i++) {
            k = rand_k [i];
            if (environment -> depth - x > M_PI/k)
                break;

            ckH = cosh(k*environment -> depth);
 
            if (isinf(ckH)) 
               continue;
 
            phas = k*y - om*t - rand_phi [i];
            ampl = factor*g*rand_a [i][2]*k/om;  

            skx = sinh(k*x);
            ckx = cosh(k*x);

            attn = skx / ckH;
            *u += ampl*attn*sin(phas);

            attn = ckx / ckH;
            *v += ampl*attn*cos(phas);

            om += d_om;
         }
      }
      if (environment -> amplitude [3][0] && w) {
         om = om_0;

         for (i = 0 ; i < nfreq ; i++) {
            k = rand_k [i];
            if (environment -> depth - x > M_PI/k)
                break;

            ckH = cosh(k*environment -> depth);
 
            if (isinf(ckH))
    	       continue;

            phas = k*z - om*t - rand_phi [i];
            ampl = factor*g*rand_a [i][3]*k/om;

	    ckx = cosh(k*x);
	    skx = sinh(k*x);

            attn = skx / ckH;
            *u += ampl*attn*sin(phas);

            attn = ckx / ckH;
            *w += ampl*attn*cos(phas);

            om += d_om;
         }
      }
   }

   return;
}

void WaveParticleAcceleration (t, x, y, z, u, v, w)
   double	 t;
   double	 x; 
   double	 y; 
   double	 z;
   double	*u;
   double	*v;
   double	*w;
{
   double  ampl;
   double  phas;
   double  attn;
   int	   i;
   double  factor;
   double  om;
   double  k;
   double  g;
   double  ckH;
   double  skx, ckx;

#ifdef HAVELAMP
   float   lamp_t;
   float   rpos [3];
   float   alamp [3];

   if (environment -> forcing == LAMP) {
      *u = *v = 0.0;
      if (w)
         *w = 0.0;
  
      return;

      lamp_t = t;

      TranslateLAMP(x, y, z, rpos);
    
      afmoor_(&lamp_t, rpos, alamp);
      *u = alamp [2];
      *v = alamp [0];
      if (w)
         *w = alamp [1];

      return;
   }
#endif

   if (analysis -> ramp_time && t < analysis -> ramp_time) 
      factor = t / analysis -> ramp_time;
   else 
      factor = 1.0;

   *u = 0.0;
   *v = 0.0;

   if (w)
      *w = 0.0;

   if (problem -> dynstat)
      return;

   if (x > environment -> depth || x < 0) {
      *u = *v = 0.0;
      if (w) *w = 0.0;
      return;
   }

   g = environment -> gravity;

   if (environment -> input_type == Regular) {

      for (i = 0 ; i < environment -> num_components [2] ; i++) {
	 k = environment -> k [2][i];

 	 ckH = cosh(k*environment -> depth);
	  
   	 if (!isinf(ckH)) {
            om = environment -> omega [2][i];

            ampl = factor*g*environment -> amplitude [2][i]*k;
            phas = k*y - om*t - environment -> phase [2][i];

	    skx = sinh(k*x);
	    ckx = cosh(k*x);

            attn = skx / ckH;
            *u += -ampl*attn*cos(phas);

            attn = ckx / ckH;
            *v += ampl*attn*sin(phas);
         }
      }

      if (w) {
         for (i = 0 ; i <= environment -> num_components [3] ; i++) {
            k = environment -> k [3][i];
	    ckH = cosh(k*environment -> depth);

      	    if (!isinf(ckH)) {
               om = environment -> omega [3][i];

               ampl = factor*g*environment -> amplitude [3][i]*k;
               phas = k*z - om*t - environment -> phase [3][i];

	       skx = sinh(k*x);
	       ckx = cosh(k*x);

               attn = skx / ckH;
               *u += -ampl*attn*cos(phas);

               attn = ckx / ckH;
               *w += ampl*attn*sin(phas);
            }
         }
      }
   }
   else if (environment -> input_type == Random) {
      if (!initialized_random)
         InitializeRandom ( );

      if (environment -> amplitude [2][0]) {
         om = om_0;

         for (i = 0 ; i < nfreq ; i++) {
            k = rand_k [i];
            if (environment -> depth - x > M_PI/k)
                break;

            ckH = cosh(k*environment -> depth);

            if (isinf(ckH))
               continue;

            phas = k*y - om*t - rand_phi [i];
            ampl = factor*rand_a [i][2]*g*k;

	    skx = sinh(k*x);
            ckx = cosh(k*x);
    
            attn = skx / ckH; 
            *u += -ampl*attn*cos(phas);

            attn = ckx / ckH;
            *v += ampl*attn*sin(phas);

            om += d_om;
         }
      }
      if (environment -> amplitude [3][0] && w) {
         om = om_0;

         for (i = 0 ; i < nfreq ; i++) {
            k = rand_k [i];
            if (environment -> depth - x > M_PI/k)
                break;

	    ckH = cosh(k * environment -> depth);

            if (isinf(ckH))
	       continue;

            phas = k*z - om*t - rand_phi [i];
            ampl = factor*rand_a [i][3]*g*k;

	    ckx = cosh(k*x);
	    skx = sinh(k*x);

            attn = skx / ckH;
            *u += -ampl*attn*cos(phas);

            attn = ckx / ckH;
            *w += ampl*attn*sin(phas);

            om += d_om;
         }
      }
   }

   return;
}

void WaveSurfaceVelocity (t, x, y, z, u, v, w)
   double	 t;
   double	 x; 
   double	 y; 
   double	 z;
   double	*u;
   double	*v;
   double	*w;
{
   double  ampl;
   double  phas;
   int     i; 
   double  factor;
   double  om;  
   double  k;
         
   if (analysis -> ramp_time && t < analysis -> ramp_time)
      factor = t / analysis -> ramp_time;
   else  
      factor = 1.0;
         
   *u = 0.0;

   *v = 0.0;
   if (w)       
      *w = 0.0;

   if (problem -> dynstat)
      return;

   if (environment -> velocity_file != NULL) {
      if (!initialized_file)
         InitializeVelocityFile ( );

      i = (int) (t/dt_f) + 1;
      *u = factor*((Uf [i+1] - Uf [i])*(t/dt_f - (i - 1)) + Uf [i]);
      *v = factor*((Vf [i+1] - Vf [i])*(t/dt_f - (i - 1)) + Vf [i]);
      if (w)
         *w = factor*((Wf [i+1] - Wf [i])*(t/dt_f - (i - 1)) + Wf [i]);
   }
   else if (environment -> input_type == Regular) {
      for (i = 0 ; i < environment -> num_components [2] ; i++) {
	 k = environment -> k [2][i];

         ampl = factor*environment -> amplitude [2][i]*environment -> omega [2][i];
         phas = k*y - environment -> omega [2][i]*t - environment -> phase [2][i];
     
         *u += ampl*sin(phas);
      }

      for (i = 0 ; i < environment -> num_components [3] ; i++) {
	 k = environment -> k [3][i];

         ampl = factor*environment -> amplitude [3][i]*environment -> omega [3][i];
         phas = k*z - environment -> omega [3][i]*t - environment -> phase [3][i];
     
         *u += ampl*sin(phas);
      }
   }
   else if (environment -> input_type == Random) {
      if (!initialized_random)
         InitializeRandom ( );

      if (environment -> amplitude [2][0]) {
         om = om_0;

         for (i = 0 ; i < nfreq ; i++) {
            k = rand_k [i];

            *u += factor*rand_a [i][2]*om*sin(k*y - om*t - rand_phi [i]);
            om += d_om;
         }
      }
      if (w && environment -> amplitude [3][0]) {
         om = om_0;

         for (i = 0 ; i < nfreq ; i++) {
            k = rand_k [i];

            *u += factor*rand_a [i][3]*om*sin(k*z - om*t - rand_phi [i]);
            om += d_om;
         }
      }
   }  

   return;
}

void InputVelocity (t, x, y, z, u, v, w)
   double	 t;
   double	 x; 
   double	 y; 
   double	 z;
   double	*u;
   double	*v;
   double	*w;
{
   double  ampl;
   double  phas;
   int     i; 
   double  factor;
   double  om;  

#ifdef HAVELAMP 
   float   lamp_t;
   float   Rlamp [3], Vlamp [3];
    
   if (environment -> forcing == LAMP) {
      lamp_t = t;
      pvmoor_(&lamp_t, Rlamp, Vlamp);

      *u = Vlamp [2];
      *v = Vlamp [0];
      if (w)
         *w = Vlamp [1];

      return;
   }
#endif
    
   if (analysis -> ramp_time && t < analysis -> ramp_time) 
      factor = t / analysis -> ramp_time;
   else  
      factor = 1.0;
         
   *u = 0.0;
   *v = 0.0;
         
   if (w)       
      *w = 0.0;

   if (problem -> dynstat)
      return;

   if (environment -> velocity_file != NULL) {
      if (!initialized_file)         
         InitializeVelocityFile ( );

      i = (int) (t/dt_f) + 1;      
      *u = factor*((Uf [i+1] - Uf [i])*(t/dt_f - (i - 1)) + Uf [i]);      
      *v = factor*((Vf [i+1] - Vf [i])*(t/dt_f - (i - 1)) + Vf [i]);
      if (w)
         *w = factor*((Wf [i+1] - Wf [i])*(t/dt_f - (i - 1)) + Wf [i]);
   }
   else if (environment -> input_type == Regular) {
      for (i = 0 ; i < environment -> num_components [1] ; i++) {
         ampl = factor*environment -> amplitude [1][i]*environment -> omega [1][i];
         phas = environment -> omega [1][i]*t - environment -> phase [1][i];
     
         *u += ampl*cos(phas);

         if (factor < 1)
            *u += environment -> amplitude [1][i]/analysis -> ramp_time*sin(phas);

      }

      for (i = 0 ; i < environment -> num_components [2] ; i++) {
         ampl = factor*environment -> amplitude [2][i]*environment -> omega [2][i];
         phas = environment -> omega [2][i]*t - environment -> phase [2][i];
     
         *v += ampl*cos(phas);

         if (factor < 1)
            *v += environment -> amplitude [2][i]/analysis -> ramp_time*sin(phas);

      }

      if (w) {
         for (i = 0 ; i < environment -> num_components [3] ; i++) {
            ampl = factor*environment -> amplitude [3][i]*environment -> omega [3][i];
            phas = environment -> omega [3][i]*t - environment -> phase [3][i];
     
            *w += ampl*cos(phas);

            if (factor < 1)
               *w += environment -> amplitude [3][i]/analysis -> ramp_time*sin(phas);

         }
      }
   }
   else if (environment -> input_type == Random) {
      if (!initialized_random)
         InitializeRandom ( );

      if (environment -> amplitude [1][0]) {
         om = om_0;

         for (i = 0 ; i < nfreq ; i++) {
            *u += factor*rand_a [i][1]*om*cos(om*t + rand_phi [i]);
            om += d_om;
         }
      }
      if (environment -> amplitude [2][0]) {
         om = om_0;

         for (i = 0 ; i < nfreq ; i++) {
            *v += factor*rand_a [i][2]*om*cos(om*t + rand_phi [i]);
            om += d_om;
         }
      }
      if (w && environment -> amplitude [3][0]) {
         om = om_0;

         for (i = 0 ; i < nfreq ; i++) {
            *w += factor*rand_a [i][3]*om*cos(om*t + rand_phi [i]);
            om += d_om;
         }
      }
   }  

   return;
}

void FroudeKrylovCoefficients (t, b, x, y, z, F, B, netw)
   double	 t;
   Buoy		 b;
   double	 x, y, z;
   double	 F [4];
   double	*B;
   double	*netw;
{
   double	   factor;
   int		   i, j;
   double	   z1, z2;
   double	   S1, S2;
   double	   ekz1, ekz2;
   double	   term_z2;
   double	   term_z1;
   double	   m;
   double	   dz;
   double	   S;
   double	   X;
   double	   d0;
   double          integral;
   double          Vg;
   static double  *level;
   static double   rho_g;
   int		   waterline;
   double	   eta;
   double	   om;
   double	   k;
   double	   draft;

   if (b -> type != Axisymmetric) 
      ExitErr ("Froude-Krylov model only implemented for Axisymmetric buoys");

   if (!initialized_random) {
      if (environment -> input_type == Random)
         InitializeRandom ( );

      rho_g = environment -> rho * environment -> gravity;

      level = (double *) malloc(sizeof(double) * b -> num_diameters);
      level --;
   }

   if (analysis -> ramp_time && t < analysis -> ramp_time)
      factor = t / analysis -> ramp_time;
   else     
      factor = 1.0;


   if ((draft = environment -> depth - x) < 0)
      ExitErr ("error in Froude Krylov forcing model - buoy out of water");

   *netw = Buoyancy(b, b -> draft, environment) - b -> w;

   F [1] = F [2] = F [3] = 0.0;
   *B = 0.0;

   d0 = 0.0;
   for (j = 1 ; j <= b -> num_diameters ; j++) {
      level [j] = b -> diameters [j].level - draft;

      if (level [j] > 0) {
          m = (b -> diameters [j].d - b -> diameters [j-1].d) /
              (level [j] - level [j-1]);

         d0 = b -> diameters [j-1].d + m*(-level [j-1]);
      }
   }

   if (environment -> input_type == Random) {

      om = om_0;
      for (i = 0 ; i < NFREQ ; i++) {
         Vg  = 0.5*om / rand_k [i]*(1.0 + 2.0*rand_k [i]*environment -> depth / 
                                        sinh(2.0*rand_k [i]*environment -> depth));

         S = 0.25*M_PI*pow(d0, 2.0);
         integral = S;

         S = 0.25*M_PI*pow(b -> diameters [1].d, 2.0);
         integral -= S*exp(rand_k [i]*level [1]);

         waterline = 0;
         for (j = 1 ; j < b -> num_diameters ; j++) {
            S1 = 0.25*M_PI*pow(b -> diameters [j].d, 2.0);
            z1 = level [j];
            ekz1 = exp(rand_k [i]*z1);

            if (level [j + 1] > 0) {
               z2 = 0;
               S2 = 0.25*M_PI*d0*d0;
               waterline = 1;
            }
            else {
               S2 = 0.25*M_PI*pow(b -> diameters [j + 1].d, 2.0);
               z2 = level [j + 1];
            }

            ekz2 = exp(rand_k [i]*z2);

            dz = z2 - z1;

            m = (S2 - S1) / dz;
   
            term_z2 = ekz2*((S1 - z1*m) + m*(z2 - 1.0/rand_k [i]));
            term_z1 = ekz1*((S1 - z1*m) + m*(z1 - 1.0/rand_k [i]));
            integral  -= (term_z2 - term_z1); 
   
            if (waterline)
               break;
         }
       
         X = rho_g * (integral + 0.25*M_PI*d0*d0);
         eta = factor*rand_a [i][2]*sin(om*t + rand_phi [i]);
         F [1] += X*eta;
         *B += rand_k [i]/(4.0*Vg*rho_g)*X*X;

         om += d_om;
      }
   }
   else {
      om = environment -> omega [2][0];
      k = environment -> k [2][0];
      Vg  = 0.5*om / k*(1.0 + 2.0*k*environment -> depth /
                              sinh(2.0*k*environment -> depth));

      S = 0.25*M_PI*pow(d0, 2.0);
      integral = S;

      S = 0.25*M_PI*pow(b -> diameters [1].d, 2.0);
      integral -= S*exp(k*level [1]);

      waterline = 0;
      for (j = 1 ; j < b -> num_diameters ; j++) {
         S1 = 0.25*M_PI*pow(b -> diameters [j].d, 2.0);
         z1 = level [j];
         ekz1 = exp(k*z1);

         if (level [j + 1] > 0) {
            z2 = 0;
            S2 = 0.25*M_PI*d0*d0;
            waterline = 1;
         }
         else {
            S2 = 0.25*M_PI*pow(b -> diameters [j + 1].d, 2.0);
            z2 = level [j + 1];
         }

         ekz2 = exp(k*z2);

         dz = z2 - z1;

         m = (S2 - S1) / dz;

         term_z2 = ekz2*((S1 - z1*m) + m*(z2 - 1.0/k));
         term_z1 = ekz1*((S1 - z1*m) + m*(z1 - 1.0/k));
   
         integral  -= (term_z2 - term_z1); 

         if (waterline)
            break;
      }
       
      X = rho_g * (integral + 0.25*M_PI*d0*d0);
      eta = factor*environment -> amplitude [2][0]
            *sin(om*t - environment -> phase [2][0]);
      F [1] = X*eta;
      *B = k/(4.0*Vg*rho_g)*X*X;
   }



   return;
}
