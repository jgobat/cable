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

/************************************************************************
 *									*
 * File:	solve.c							*
 * 									*
 * Description:	contains functions which are generic to all solution	*
 *		types (2D vs. 3D and static vs. dynamic)		*
 *									*
 ************************************************************************/
 
# include <stdio.h>
# include <stdlib.h>
# include <malloc.h>
# include <math.h>
# include "compress.h"
# include "problem.h"
# include "solve.h"
# include "error.h"
# include "tension.h"
# include "output.h"
# include "control.h"

# define BLANK (-99999999)

# define SQR(a) ((a)*(a))

# define MAX(a,b) ((a) > (b) ? (a) : (b))

extern Problem *problem;
extern Environment *environment;
extern Analysis *analysis;
extern Debug debug;

# ifdef DEBUG

static void
dumpRHS(char *fname, double *b, int n)
{
    int   i;
    FILE *fp;
    fp = fopen(fname, "w");
    if (fp == NULL)
        return; 

    for (i = 1 ; i <= n ; i++)
        fprintf(fp, "%g\n", b[i]);

    fclose(fp);
}

static void
dumpCSR(char *fname, double *a, int *rowptr, int *colind, int n, int nnz)
{
    int   i,j;
    FILE *fp;
    fp = fopen(fname, "w");
    if (fp == NULL)
        return; 

    fprintf(fp, "%d %d\n", n, nnz);
    for (i = 1 ; i <= n+1 ; i++)
        fprintf(fp, "%d ", rowptr[i]);
    fprintf(fp, "\n");
    for (i = 1 ; i <= nnz ; i++)
        fprintf(fp, "%d %g\n", colind[i], a[i]);

    fclose(fp); 
}


static void dump (fname, a, rowptr, colind, n)
   char		*fname;
   double	*a;
   int		*rowptr;
   int		*colind;
   int		n;
{
   int		i, j;
   FILE		*fp;

   fp = fopen(fname, "w");
   if (fp == NULL)
      return;
 
   for (i = 1 ; i <= n ; i++) {
      for (j = rowptr [i] ; j <= rowptr [i + 1] - 1 ; j++)
         fprintf (fp, "%d %d %g\n", i, colind [j], a [j]);
   }

   fclose (fp);

   return;
}

# endif

extern int nspiv ( );

int 
SolveDE (difeq, update, itmax, conv, base_relax, scalv, 
             ne, nb, nb_branch, nj_compat, all_node, all_m, node, m, s, 
             tm , dt, current_factor, step)
   void		(*difeq) ( );
   void		(*update) ( );
   int	       *itmax;
   double      *conv;
   double      *base_relax;
   double    	scalv [];
   int	     	ne, nb;
   int		nb_branch, nj_compat;
   Node	       *all_node;
   int		all_m;
   Node	       *node;
   int		m;
   double     **s;
   double	tm;
   double	dt;
   double       current_factor;
   int	        step;
{
   int               viva_passes; 
   static int    	 nvars;
   static int	 	 nnz;
   int           	 it;
   double        	 err, errj;
   int		 	 sing;
   int		 	 stall_count;
   double        	 relax;
   double	 	 prev_err;
   double	 	 prev_relax;
   double		 start_relax;
   double		 adapt_up, adapt_down;
   double		 data;
   EquationType	 	 eq_type;
   int			 need_top_bc_eq;
   int			 num_rows, num_cols;
   int			 row_count;
   int			 nz_count;
   int		 	 i, j, k, n, jj, mm;
   Node			 nm;
   int			 km, kj;
   int			 njn;
   int		   	 col;
   int			 rhs_col;
   int			 num_blocks;
   static int		*block_col = NULL;
   static int		 jsf;
   static double	*a = NULL;
   static double	*rhs;
   static double	*sol;
   static int		*colind;
   static int		*rowptr;
   static int		 max;
   static int		*cols;
   static int		*rows;  
   static int		*icols;
   static double        *yu, *u;
   static int		*p, *iu, *ju;
   int			 info;
   int			 mode;
   static int	 	 prev_mode = 0;
   static int		 prev_nvars;
   int			 initialize;
 
   if (block_col == NULL) {
      block_col = (int *) malloc(sizeof(int) * (2 + problem -> junction_size));
      block_col --;
   }

   viva_passes = 0;
 
   nvars = ne*m;
   initialize = 0;

   mode = (step ? 1 : 2);
   if (mode != prev_mode) {
      prev_mode = mode;
      initialize = 1;
   }

   if (nvars != prev_nvars) {
      prev_nvars = nvars;
      initialize = 1;
   }   
 
 
   if (initialize) {

      jsf = 2*ne + 1;

      if (a) {				/* free all the work arrays because */
         a++; free(a);			/* we have switched from a static   */
         colind++; free(colind);	/* to a dynamic problem		    */
         rowptr++; free(rowptr);
         rhs++; free(rhs);
         sol++; free(sol);

         yu++; free(yu);
         iu++; free(iu);
         ju++; free(ju);
         u++; free(u);
         p++; free(p);

         rows ++; free(rows); 
         cols ++; free(cols);
         icols ++; free(icols);
      }


      rowptr = (int *) malloc(sizeof(int) * (nvars + 1)); rowptr --;


      n = 1;
      rowptr [n] = 1;
      n ++;

	/*
	 * to get the true map and count of non-zeroes, we make a pass
	 * through all the nodes, getting the Jacobian at each
	 * and scanning it for valid entries
	 */

      nnz = 0;

      need_top_bc_eq = 0;
      for (k = 1 ; k <= m ; k++) {

         if (node [k] -> position == Junction) {
            nm = node [k] -> prev_active;
            km = nm -> active_number;
   
            njn = node [k] -> segment -> junction.num_nodes;

            num_cols = (2 + njn)*ne;
            rhs_col = num_cols + 1;

            num_rows = ne + nj_compat*njn;
            
            eq_type = Junction;
         }
         else if (node [k] -> position == BottomBoundary) {
            rhs_col = jsf;
            num_rows = nb;
            num_cols = ne;

    	    km = 1;
            nm = node[1];
            eq_type = BottomBoundary;
         }
         else if (node [k] -> position == BranchStart) {
            rhs_col = jsf;
            num_rows = nb_branch;
            num_cols = ne;

    	    km = k;
            nm = node[k];
            eq_type = BranchStart;
         }
         else if (need_top_bc_eq && (node [k] -> position == BranchTerminal ||
				     node [k] -> position == TopBoundary)) {
            rhs_col = jsf;

            if (node [k] -> position == BranchTerminal)
               if (node [k] -> segment -> branch -> terminal -> loop_main_node)
                  num_rows = ne - nb_branch - 2*nj_compat;
               else
                  num_rows = ne - nb_branch - nj_compat;
            else
               num_rows = ne - nb;

            num_cols = ne;
            km = k;
            nm = node[k];

            eq_type = node [k] -> position;
            need_top_bc_eq = 0;
         }
         else {	/* Connection, Cable, TopBoundary, BranchTerminal */
            rhs_col = jsf;
            num_rows = ne;

            num_cols = 2*ne;
           
            nm = node[k] -> prev_active;
            km = nm -> active_number;

            if (node [k] -> position == TopBoundary ||
                node [k] -> position == BranchTerminal) {

	           need_top_bc_eq = 1;	/* still need BC eqs at this node */
               eq_type = Cable;
            }
   	        else
               eq_type = node [k] -> position;
         }
      
         difeq(eq_type, node[k], nm, ne, rhs_col, num_rows, s, node, 
               tm, dt, current_factor);

         for (i = 1 ; i <= num_rows ; i++) {
            nz_count = 0;
            for (j = 1 ; j <= num_cols ; j++)  {
               if (s [i][j] != BLANK) {
                  nz_count ++;
               } 
            }
            rowptr [n] = rowptr [n-1] + nz_count;
            nnz += nz_count;

            n ++;
         }

         if (need_top_bc_eq)
	        k --;
      }

      a      = (double *) malloc(sizeof(double) * nnz);   a--;       
      colind = (int *) malloc(sizeof(int) * nnz);         colind --;        
      rhs    = (double *) malloc(sizeof(double) * nvars); rhs --;   
      sol    = (double *) malloc(sizeof(double) * nvars); sol --;    

      max = 10*nnz;

      yu = (double *) malloc(sizeof(double) * nvars); yu --;
      u  = (double *) malloc(sizeof(double) * max);   u --;
     
      iu = (int *) malloc(sizeof(int) * (nvars + 1)); iu --;
      ju = (int *) malloc(sizeof(int) * max);         ju --;
      p  = (int *) malloc(sizeof(int) * (nvars + 1)); p --;

      rows = (int *) malloc(sizeof(int) * nvars); rows --;
      cols = (int *) malloc(sizeof(int) * nvars); cols --;
      icols = (int *) malloc(sizeof(int) * nvars); icols --;

      for (j = 1 ; j <= nvars ; j++) {
         rows [j] = j;
         cols [j] = j;
         icols [j] = j;
      }
   }

   relax          = *base_relax;
   prev_relax     = relax;
   start_relax    = relax;

   prev_err = 1e9;
   stall_count = 0;

   adapt_up = analysis -> relax_up;
   adapt_down = analysis -> relax_down;
   
   for  (it = 1; it <= *itmax ; it++) {
      n = 1;

      difeq(BottomBoundary, node[1], node[1], ne, jsf, nb, s, node, 
            tm, dt, current_factor);
      
      for (i = 1 ; i <= nb ; i++) {
         for (j = 1 ; j <= ne ; j++) {
            if (s [i][j] != BLANK) {
               a [n] = s [i][j];
               colind [n] = j;
               n ++;
            }
         }
         rhs [i] = s [i][jsf]; 
      }
      
      row_count = nb;
      need_top_bc_eq = 0;

      for (k = 2 ; k <= m ; k++) {

         if (node [k] -> position == Junction) {
            nm = node[k] -> prev_active;
            km = nm -> active_number;

            njn = node [k] -> segment -> junction.num_nodes;

            rhs_col = (2 + njn)*ne + 1;
            num_blocks = (2 + njn);
          
            num_rows = ne + nj_compat*njn;
            
            block_col [1] = (km - 1)*ne;
            block_col [2] = (k - 1)*ne;

            for (i = 1 ; i <= njn ; i++) {
               kj = node [k] -> segment -> junction.node [i] -> active_number;
               block_col [i + 2] = (kj - 1)*ne;
            }
            eq_type = Junction;
         }
         else if (node [k] -> position == BranchStart) {
            rhs_col = jsf;
            num_rows = nb_branch;

            num_blocks = 1;

            block_col [1] = (k - 1)*ne;
            nm = node[k];
            km = k;
            eq_type = BranchStart;
         }
         else if (need_top_bc_eq && (node [k] -> position == BranchTerminal ||
				     node [k] -> position == TopBoundary)) {
            rhs_col = jsf;
            if (node [k] -> position == BranchTerminal)
               if (node [k] -> segment -> branch -> terminal -> loop_main_node)
                  num_rows = ne - nb_branch - 2*nj_compat;
               else
                  num_rows = ne - nb_branch - nj_compat;
            else
               num_rows = ne - nb;

            num_blocks = 1;

            block_col [1] = (k - 1)*ne;
            nm = node[k];
            km = k;
            eq_type = node [k] -> position;
            need_top_bc_eq = 0;
         }
         else {	/* Connection, Cable, TopBoundary, BranchTerminal */
            rhs_col = jsf;
            num_rows = ne;

            num_blocks = 2;
           
            nm = node[k] -> prev_active;
            km = nm -> active_number; 

            block_col [1] = (km - 1)*ne;
            block_col [2] = (k - 1)*ne;
            

            if (node [k] -> position == TopBoundary ||
                node [k] -> position == BranchTerminal) {

               need_top_bc_eq = 1;	/* still need BC eqs at this node */
               eq_type = Cable;
            }
   	        else
               eq_type = node [k] -> position;
         }
         
         difeq(eq_type, node[k], nm, ne, rhs_col, num_rows, s, node, 
               tm, dt, current_factor);

         for (i = 1 ; i <= num_rows ; i++) {
            for (j = 1 ; j <= num_blocks ; j++) {

               for (jj = 1 ; jj <= ne ; jj++) {
                  data = s [i][(j - 1)*ne + jj];
                  if (data != BLANK) {
                     col = block_col [j] + jj;
                     a [n] = data;
                     colind [n] = col;
                     n ++; 
                  }
               }
            }

            rhs [row_count + i] = s [i][rhs_col];
         }

         row_count += num_rows;

         if (need_top_bc_eq)
            k --;
      }


/*
      for (i = 1 ; i <= row_count ; i++) 
         fprintf (stderr, "%d, %d : %g\n", step, i, rhs [i]);

      if (step == 0 && it > 1) {
         dump("dump.dat", a, rowptr, colind, nvars);
         for (i = 1 ; i <= row_count ; i++) 
            printf("%g %g\n", rhs [i], sol [i]);

         exit (0);
      }
*/

/*
      dumpCSR("csr.dat", a, rowptr, colind, nvars, nnz);
      dumpRHS("rhs.dat", rhs, nvars);
*/
      sing = nspiv(nvars, rowptr, colind, a, rhs, max, rows, cols, icols, sol,
                   yu, p, iu, ju, u, &info);

      if (sing)  {
         DisplayMessage("singularity (%d @ %d)", sing, info);
/*
         for (i = rowptr [info] ; i <= rowptr [info] + 2*ne ; i++) 
            fprintf (stderr,"%g ", a [i]);

         fprintf (stderr,"\n");

         for (i = rowptr [info] ; i <= rowptr [info] + 2*ne ; i++) 
            fprintf (stderr,"%d ", colind [i]);

         fprintf (stderr,"\n");
*/

         return sing;
      }

      err = 0.0;

      for (j = 1 ; j <= ne ; j++) {
         errj = 0.0;

         for (k = 1 ; k <= m ; k++) {
            data = sol [(k - 1)*ne + j];

            if (!finite(data)) {
               DisplayMessage("NaN/Inf error, node = %d, DOF = %d", k, j);
               SetError(C_NANINF);
               return 1;
            }

            errj += fabs(data);
         }

         err += errj / scalv [j];
      }

      err /= nvars;

      if (start_relax != *base_relax) {
         relax = *base_relax;
         start_relax = *base_relax;
         stall_count = 0;
      }
      else if (err > prev_err) {
         relax = prev_relax / adapt_down;
         if (relax < *base_relax/1000.0) {
            relax = *base_relax/1000.0;
            stall_count ++;
         }
      }
      else {
         relax = relax * adapt_up;
         stall_count = 0;
         if (relax > *base_relax)
            relax = *base_relax;
      }

      prev_relax = relax;
      prev_err = err;

      for (k = 1 ; k <= m ; k++) {
         for (j = 1 ; j <= ne ; j++)  {
            data = sol [(k - 1)*ne + j];

            node [k] -> Y[j] -= relax*data;
         }

         if (node[k] -> Y[1] < -1.0) {
            DisplayMessage("compression error (%d)", k);
            SetError(C_COMPRESSION);
            return 1;
         }
      }

      update (node, m, tm, dt);

      DisplayInfo (step, tm, it, err, relax);
      // DisplayInfo (step, tm, it, 0.0, relax);
  
      if (!finite(err)) {
         SetError(C_NANINF);
         return 1;
      }
      

# if defined (GUI) || defined (WINGUI)
     ControlProcessEvents (problem -> solution);
     if (problem -> solution -> userQuit) {
         return 1;
     }
# endif

        /*
         * do the convergence check only after we flush the X queue 
         * so that we actually see the error go below the tolerance
         */

      if (debug.status && 0)
         FakeDynamicSnapshot(debug.out, 
			     all_node, all_m, debug.decimate, 1);

      if (err < *conv) {
#ifdef HAVE_VIVA
         if (viva_passes ++ < analysis -> viva_iterations) {
            DisplayMessage("VIVA - updating drag coefficients");
            vivaWrapper(all_node, all_m, node, m, problem, analysis, environment);
         }

         if (viva_passes >= analysis -> viva_iterations)
             return 0;
#else
        return 0;
#endif
      }
      
      if (stall_count > analysis -> stall_limit) {
         prev_relax = relax = *base_relax;
         DisplayMessage ("iterations stalled");
         stall_count = 0;
      }
   }

   DisplayMessage ("iterations exceeded");
   SetError(C_MAXITERATIONSEXCEEDED);

   return *itmax + 1;
}

double check (t, dt, sample)        
   double       t;
   double       dt;
   double       sample; 
{  
   double       n;
   
   n = floor((t + dt/2.0) / sample);
   
   return t - n*sample;
}

 
double sign(x)       
   double       x;
{
   return (x < 0 ? -1.0 : 1.0);
}

static void ResolveAnchoredBranches (node, rlx, twoD)
   Node   *node;
   double  rlx;
   int     twoD;
{
   int    i;
   Branch br;
   double x;
   double y;
   double z;

   for (i = 1 ; i <= problem -> num_branch ; i++) {
      br = problem -> branch [i];

      if (br -> terminal -> anchor) {
         x = br -> terminal -> x;
         y = br -> terminal -> y;
         if (!twoD)
            z = br -> terminal -> z;
         else
            z = 0.0;

         br -> terminal -> xforce -= rlx*(br -> last -> x - x);
         br -> terminal -> yforce -= rlx*(br -> last -> y - y);
         if (!twoD)
            br -> terminal -> zforce -= rlx*(br -> last -> z - z);

         DisplayMessage("branch %d anchor: z = %g, x = %g, y = %g", i,
                  br -> last -> x, br -> last -> y,
		  br -> last -> z);
      }
   }

   return;
}

int 
ResolvePositionedConnectors(Segment *seg, int nseg, double rlx, int twoD)
{
    int     j;
    double  x, y, z;
    double  e1,e2,e3;
    double  tol;
    int     converged;

    tol = analysis -> outer_tolerance;
    converged = 1;

    for (j = 1 ; j <= nseg ; j++) {
        if (seg[j] -> connector) {
            x = seg[j] -> connector_x;
            y = seg[j] -> connector_y;
            if (!twoD)
                z = seg[j] -> connector_z;
            else
                z = 0.0;
                
            e1 = fabs(seg[j] -> last_active -> x - x) / seg[j] -> last_active -> prev -> prev -> ds; 
            e2 = fabs(seg[j] -> last_active -> y - y) / seg[j] -> last_active -> prev -> prev -> ds; 
            e3 = fabs(seg[j] -> last_active -> z - z) / seg[j] -> last_active -> prev -> prev -> ds; 
            seg[j] -> connector_xforce -= rlx*(seg[j] -> last_active -> x - x);
            seg[j] -> connector_yforce -= rlx*(seg[j] -> last_active -> y - y);
            if (!twoD)
                seg[j] -> connector_zforce -= rlx*(seg[j] -> last_active -> z - z);
            if (e1 > tol || e2 > tol || (!twoD && e3 > tol))
                converged = 0;

            DisplayMessage("connector %d: pos=%g,%g,%g: x=%g, y=%g, z=%g", j,
                  x, y, z, seg[j] -> last_active -> x, seg[j] -> last_active -> y,
		          seg[j] -> last_active -> z);
            DisplayMessage("             xf=%g, yf=%g, zf=%g", 
                           seg[j] -> connector_xforce, 
                           seg[j] -> connector_yforce, 
                           seg[j] -> connector_zforce);
        }
    }

    return converged;
}

int 
ResolveAnchor(Node *node, int nn, Segment *seg, int nseg, int twoD)
{
   double rlx_x, rlx_y, rlx_z;
   static double prev_x, prev_y, prev_z;
   static double dfx, dfy, dfz;
   double dx, dy, dz;
   static int init = 0;
   double	x, y, z;
   double	rlx;
   double	tol;
   double	e1, e2, e3;
   int      converged;

   rlx = analysis -> outer_relaxation;
   tol = analysis -> outer_tolerance;

   x = problem -> terminal [2] -> x;
   y = problem -> terminal [2] -> y;
   if (!twoD)
      z = problem -> terminal [2] -> z;
   else
      z = 0.0;
  
   fprintf(stderr,"nn = %d, node_x=%f, anchor_x=%f, ds=%f, node_y=%f, anchor_y=%f, node_z=%g, anchor_z=%g\n", 
           nn,
           node[nn] -> x, x, 
           node [nn-1] -> ds,
           node[nn] -> y, y,
           node[nn] -> z, z);
   e1 = fabs(node [nn] -> x - x) / node [nn - 1] -> ds; 
   e2 = fabs(node [nn] -> y - y) / node [nn - 1] -> ds; 
   if (!twoD)
      e3 = fabs(node [nn] -> z - z) / node [nn - 1] -> ds; 
   else
      e3 = 0.0;

   fprintf(stderr,"e1 = %f, e2 = %f, e3 = %f\n", e1, e2, e3);
   DisplayAuxError(e1 > e2 ? (e1 > e3 ? e1 : e3) : (e2 > e3 ? e2 : e3));

   if (init && analysis -> static_outer_method == Variable) {
       dx = prev_x - node[nn] -> x;
       dy = prev_y - node[nn] -> y;
       dz = prev_z - node[nn] -> z;

       rlx_x = 0.05*dfx / dx;
       rlx_y = 0.05*dfy / dy;
       rlx_z = 0.05*dfz / dz;
   }
   else {
       init = 1;
       rlx_x = rlx;
       rlx_y = rlx;
       rlx_z = rlx;
   }

   fprintf(stderr,"rlx_x = %f, rlx_y = %f, rlx_z = %f\n", rlx_x, rlx_y, rlx_z);
   prev_x = node[nn] -> x;
   prev_y = node[nn] -> y;
   prev_z = node[nn] -> z;
   


   fprintf(stderr,"xf=%f, yf=%f, zf=%f\n", 
           problem -> terminal[2] -> xforce,
           problem -> terminal[2] -> yforce,
           problem -> terminal[2] -> zforce);
   
   dfx = rlx_x*(node [nn] -> x - x);
   dfy = rlx_y*(node [nn] -> y - y);
   problem -> terminal [2] -> xforce -= dfx;
   problem -> terminal [2] -> yforce -= dfy;
   if (!twoD) {
      dfz = rlx_z*(node [nn] -> z - z);
      problem -> terminal [2] -> zforce -= dfz;
   }

   fprintf(stderr,"xf'=%f, yf'=%f, zf'=%f\n", 
           problem -> terminal[2] -> xforce,
           problem -> terminal[2] -> yforce,
           problem -> terminal[2] -> zforce);

   ResolveAnchoredBranches (node, rlx, twoD);
    
   if (problem -> type == HorizontalDrifter)
       converged = ResolvePositionedConnectors(seg, nseg, rlx, twoD);
   else
       converged = 1;

   if (twoD) {
      if (e1 <= tol && e2 <= tol && converged)
        return 1;
   }
   else {
      if (e1 <= tol && e2 <= tol && e3 <= tol && converged)
         return 1;
   }

   return 0;
}

#define DBL_EPSILON 1e-12

double ResolveBuoy (node, nn, step, twoD)
   Node		*node;
   int           nn;
   int		 step;
   int		 twoD;
{
   static double x1, x2, x3;
   static double fx1, fx2, fx3;
   static int    found_second_limit;
   Buoy		 b;
   double	 calculated_draft;
   double	 Uw, Vw, Ww;
   double    Vs, Ws;
   double	 err;

   b = problem -> terminal [2] -> buoy;
   calculated_draft = environment -> depth - node [nn] -> x;

   err = (calculated_draft - b -> draft);
   
   DisplayMessage("d' = %g, d = %g, e = %g", b -> draft, calculated_draft, err);

   err = fabs(err / b -> draft);
   DisplayAuxError(err);
   if (err < analysis -> outer_tolerance)
      return err;

   if (step == 1) {
      x1 = b -> draft;
      fx1 = calculated_draft - x1;

	/*
	 * if the initial guess came from the catenary solution we
	 * must have used max_draft - if it came from a shooting 
	 * solution then the guessed draft was probably pretty
	 * darn close to the real deal
	 */

      if (analysis -> static_initial_guess == Catenary) {
         if (fx1 > 0) { 
            DisplayMessage("buoy submerged at max draft");
            SetError(C_BUOYNOTATSURFACE);
            return 0;
         }

         b -> draft = analysis -> outer_relaxation*b -> draft;
      }
      else if (fx1 < -DBL_EPSILON) 
         b -> draft *= analysis -> outer_relaxation;
      else
         b -> draft /= analysis -> outer_relaxation;

      found_second_limit = 0;
   }
   else if (!found_second_limit) {
      x2 = b -> draft;
      fx2 = calculated_draft - x2;

      if (b -> draft == b -> min_draft && fx2 < -DBL_EPSILON) {
         DisplayMessage("buoy not in water at min draft");
         SetError(C_BUOYNOTINWATER);
         return 0;
      }
      else if (fx2*fx1 > DBL_EPSILON) {
         if (fx2 < -DBL_EPSILON) {
            x1 = b -> draft;
            b -> draft *= analysis -> outer_relaxation;
         }
         else
            b -> draft /= analysis -> outer_relaxation;

         if (b -> draft < b -> min_draft) 
            b -> draft = b -> min_draft;
      }
      else {
         found_second_limit = 1;
         if (analysis -> static_outer_method == Secant)
            b -> draft = x2 - fx2*(x2 - x1)/(fx2 - fx1);
         else if (analysis -> static_outer_method == Bisection)
            b -> draft = (x2 + x1)/2.0;

         DisplayMessage("2nd limit: x1 = %g, x2 = %g, draft = %g",
                        x1, x2, b -> draft);
      }
   }
   else {
      fx3 = calculated_draft - b -> draft;
      x3 = b -> draft;
      fprintf(stdout,"Before: fx3 = %g, fx1 = %g, x1 = %g, x2 = %g, b -> draft = %g\n", fx3, fx1, x1, x2, b -> draft);

      if (analysis -> static_outer_method == Secant) {
         b -> draft = x3 - fx3*(x3 - x2)/(fx3 - fx2);

         x2 = x3;
         fx2 = fx3; 
      }
      else if (analysis -> static_outer_method == Bisection) {
          if (fx3*fx1 < -DBL_EPSILON) 
             x2 = b -> draft;
          else {
             x1 = b -> draft;
             fx1 = fx3;
          }

          b -> draft = (x2 + x1)/2.0;
      }

      fprintf(stdout,"After: fx3 = %g, fx1 = %g, x1 = %g, x2 = %g, b -> draft = %g\n", fx3, fx1, x1, x2, b -> draft);
   }

   if (node [nn] -> x - environment -> depth > DBL_EPSILON) {
 //     WaveStokesVelocity(environment -> depth - b->max_draft , &Vs, &Ws);
      Current(0.0, environment -> depth, node [nn] -> y, node [nn] -> z, 
              &Uw, &Vw, &Ww);
   }
   else {
//      WaveStokesVelocity(node[nn] -> x, &Vs, &Ws);
      Current(0.0, node [nn] -> x, node [nn] -> y, node [nn] -> z, 
              &Uw, &Vw, &Ww);
   }
   Vs = Ws = 0;
   fprintf(stdout, "current velocity %f, %f\n", Vw, Ww); 
   fprintf(stdout, "stokes velocity %f, %f\n", Vs, Ws); 
   Vw += Vs;
   Ww += Ws;

   if (problem -> type == Deployment) {
      Vw -= problem -> terminal [1] -> yspeed.value;
      Ww -= problem -> terminal [1] -> zspeed.value;
   }

   problem -> terminal [2] -> xforce = Buoyancy(b, b -> draft, environment) - b -> w;
   problem -> terminal [2] -> yforce = 
         0.5*environment -> rho*b -> Cdn*ProjectedArea(b, b -> draft)*Vw*fabs(Vw);
   if (!twoD)
      problem -> terminal [2] -> zforce = 
         0.5*environment -> rho*b -> Cdn*ProjectedArea(b, b -> draft)*Ww*fabs(Ww);

   ResolveAnchoredBranches (node, 20.0, twoD);

   return err;
}

int ResolveSpeed (node, nn, step)
   Node		*node;
   int       nn;
   int		 step;
{
   static double x1, x2;
   static double fx1, fx2, fx3;
   static double prev_fx;
   Buoy		 b;
   Buoy	 	 b1;
   double	 uack, vack, wack;
   double	 err;
   double	 Fy, T;
   double	 Fy_wind;
   double	 drag;
   double	 speed;
   double	 phi;
   double        v_rel;
   static int	 found_second_limit;
   

   b = problem -> terminal [2] -> buoy;
   b1 = problem -> terminal [1] -> buoy;

   speed = problem -> terminal [2] -> yspeed.value;

   T = Tension(node[nn] -> Y[1], node [nn] -> material);
   phi = node[nn] -> Y [3];
   if ((b -> draft = Draft(b -> w + T*cos(phi) - node[nn] -> Y[2]*sin(phi), 
                           b, environment)) < 0)
      b -> draft = b -> max_draft;
   Current(0.0, environment -> depth, 0.0, 0.0, &uack, &vack, &wack);

   WindDrag(0.0, problem -> terminal [2] -> buoy, &Fy_wind, NULL);

	/* 
	 * drag force on the buoy - correct sign
	 */

   drag = 0.5*environment -> rho*ProjectedArea(b, b -> draft)*b -> Cdn*(vack - problem -> terminal[2] -> yspeed.value)*fabs(vack - problem -> terminal[2] -> yspeed.value);

	/*
	 * horizontal force applied to the buoy from the line - correct sign
	 */

   Fy = -(T*sin(phi) + node[nn] -> Y[2]*cos(phi)); // -Fy_wind ???

	/*
	 * sum of horizontal forces should be zero so that's our metric
	 */

   err = fabs((Fy + drag) / b1 -> w);
   fprintf(stderr,"Fy = %g, drag = %g, V = %g, vack = %g, draft = %g, Fx = %g\n", Fy, drag, speed, vack, b -> draft, T*cos(phi) - node[nn]->Y[2]*sin(phi));
   DisplayMessage("err = %g, V' = %g", err, speed);

   DisplayAuxError(err);
   if (err < analysis -> outer_tolerance)
      return 1;

   if (step == 1) {
      x1  = speed;
      fx1 = Fy + drag;

      found_second_limit = 0;
      speed = analysis -> outer_relaxation*speed;
   }
   else if (!found_second_limit) {
      x2 = speed;
      fx2 = Fy + drag;
      if (fx2*fx1 > 0) {
         x1 = speed;
         speed = speed*analysis -> outer_relaxation;
      }
      else {
         fprintf(stderr,"step %d, found second limit: V = %g\n", step, speed);
         found_second_limit = 1;
         speed = x2 - fx2*(x2 - x1)/(fx2 - fx1); 
         speed = (x1 + x2)/2.0;
         prev_fx = fx1;
      }
   }
   else {
      fx3 = Fy + drag;

      if (fx3*fx1 < 0.0) {
         x2 = speed;
         fx2 = fx3;
/*
         if (fx3*prev_fx > 0.0)
            fx1 = fx1/2.0;
*/
      }
      else {
         x1 = speed;
         fx1 = fx3;
/*
         if (fx3*prev_fx > 0.0)
            fx2 = fx2/2.0;
*/
      }

      prev_fx = fx3;

      speed = x2 - fx2*(x2 - x1)/(fx2 - fx1);
      speed = (x1 + x2)/2.0;
   }

   Current(0.0, node[1] -> x, 0.0, 0.0, &uack, &vack, &wack);

   v_rel = vack - speed;
   problem -> terminal [1] -> yforce =
       0.5*environment -> rho*b1 -> S*v_rel*fabs(v_rel)*b1 -> Cdn;

   problem -> terminal [2] -> yspeed.value = speed;
 
   return 0;
}

int ResolveSpeed3D (node, Y, nn, step)
   Node		*node;
   double      **Y;
   int           nn;
   int		 step;
{
   static double Vx1, Vx2;
   static double Vfx1, Vfx2, Vfx3;
   static double Vprev_fx;
   static double Wx1, Wx2;
   static double Wfx1, Wfx2, Wfx3;
   static double Wprev_fx;
   Buoy		 b;
   Buoy	 	 b1;
   double	 uack, vack, wack;
   double	 err;
   double	 Fy, Fz;
   double	 Fy_wind, Fz_wind;
   double        v_rel;
   double	 w_rel;
   double	 Vcalc;
   double	 Wcalc;
   double	 Vspeed, Wspeed;
   double	 T, Sn, Sb;
   double	 B0, B1, B2, B3;

   b = problem -> terminal [2] -> buoy;
   b1 = problem -> terminal [1] -> buoy;


   Current(0.0, environment -> surface, 0.0, 0.0, &uack, &vack, &wack);

   WindDrag(0.0, problem -> terminal [2] -> buoy, &Fy_wind, &Fz_wind);

   B0 = Y [4][nn];
   B1 = Y [5][nn];
   B2 = Y [6][nn];
   B3 = Y [7][nn];

   T = Tension(Y [1][nn], node [nn] -> material);
   Sn = Y [2][nn];
   Sb = Y [3][nn];

   Fy =   2.0*T*(B1*B2 + B0*B3) 
        + Sn*(B0*B0 - B1*B1 + B2*B2 - B3*B3)
        + 2.0*Sb*(B2*B3 - B0*B1);
   Fz =   2.0*T*(B1*B3 - B0*B2) 
        + 2.0*Sn*(B2*B3 + B0*B1) 
        + Sb*(B0*B0 - B1*B1 - B2*B2 + B3*B3);

   Fy = Fy - Fy_wind;
   Fz = Fz - Fz_wind;

   Vspeed = problem -> terminal [2] -> yspeed.value;
   Wspeed = problem -> terminal [2] -> zspeed.value;

   Vcalc = vack - sign(Fy)*sqrt(2.0*fabs(Fy)/environment -> rho/b -> Cdn/b -> S);
   Wcalc = wack - sign(Fz)*sqrt(2.0*fabs(Fz)/environment -> rho/b -> Cdn/b -> S);

   if (Wspeed != 0.0 && Vspeed != 0.0)
      err = fabs((Vcalc - Vspeed) / Vspeed) + fabs((Wcalc - Wspeed) / Wspeed);
   else if (Wspeed != 0.0)
      err = fabs((Wcalc - Wspeed) / Wspeed);
   else if (Vspeed != 0.0)
      err = fabs((Vcalc - Vspeed) / Vspeed);
   else
      return 1;

   DisplayMessage("e=%g, Vg=%g, Vs=%g, Wg=%g, Ws=%g", 
                  err, Vspeed, Vcalc, Wspeed, Wcalc);          

   DisplayAuxError(err);
   if (err < analysis -> outer_tolerance)
      return 1;


   if (step == 1) {
      Vx1  = Vspeed;
      Vfx1 = Vcalc - Vx1;

      Vspeed = 0.05*Vspeed;

      Wx1  = Wspeed;
      Wfx1 = Wcalc - Wx1;

      Wspeed = 0.05*Wspeed;
   }
   else if (step == 2) {
      Vx2 = Vspeed;
      Vfx2 = Vcalc - Vx2;

      Vspeed = Vx2 - Vfx2*(Vx2 - Vx1)/(Vfx2 - Vfx1); 

      Vprev_fx = Vfx1;

      Wx2 = Wspeed;
      Wfx2 = Wcalc - Wx2;

      Wspeed = Wx2 - Wfx2*(Wx2 - Wx1)/(Wfx2 - Wfx1); 

      Wprev_fx = Wfx1;
   }
   else {
      Vfx3 = Vcalc - Vspeed;

      if (Vfx3*Vfx1 < 0.0) {
         Vx2 = Vspeed;
         Vfx2 = Vfx3;
         if (Vfx3*Vprev_fx > 0.0)
            Vfx1 = Vfx1/2.0;
      }
      else {
         Vx1 = Vspeed;
         Vfx1 = Vfx3;
         if (Vfx3*Vprev_fx > 0.0)
            Vfx2 = Vfx2/2.0;
      }

      Vprev_fx = Vfx3;

      if (Vfx2 == Vfx1)
         Vspeed = 0.0;
      else
         Vspeed = Vx2 - Vfx2*(Vx2 - Vx1)/(Vfx2 - Vfx1);

      Wfx3 = Wcalc - Wspeed;

      if (Wfx3*Wfx1 < 0.0) {
         Wx2 = Wspeed;
         Wfx2 = Wfx3;
         if (Wfx3*Wprev_fx > 0.0)
            Wfx1 = Wfx1/2.0;
      }
      else {
         Wx1 = Wspeed;
         Wfx1 = Wfx3;
         if (Wfx3*Wprev_fx > 0.0)
            Wfx2 = Wfx2/2.0;
      }

      Wprev_fx = Wfx3;

      if (Wfx2 == Wfx1)
         Wspeed = 0.0;
      else
         Wspeed = Wx2 - Wfx2*(Wx2 - Wx1)/(Wfx2 - Wfx1);
   }

   Current(0.0, 0.0, 0.0, 0.0, &uack, &vack, &wack);

   problem -> terminal [2] -> yspeed.value = Vspeed;
   problem -> terminal [2] -> zspeed.value = Wspeed;

   v_rel = vack - Vspeed;
   w_rel = wack - Wspeed;

   problem -> terminal [1] -> yforce =
          0.5*environment -> rho*b1 -> S*v_rel*fabs(v_rel)*b1 -> Cdn;
   problem -> terminal [1] -> zforce =
          0.5*environment -> rho*b1 -> S*w_rel*fabs(w_rel)*b1 -> Cdn;

   return 0;
}

int ResolveWebster (it, node, nn, twoD)
   int		 it;
   Node         *node;
   int           nn;
   int		 twoD;
{
   double	x;
   double	tol;
   double	phi;
   static double fx1, fx2, fx3;
   static double x1, x2;
   double	 phi0; 
   double	 e;

   tol = analysis -> outer_tolerance;

   x = problem -> terminal [2] -> x;
  
   e = fabs(node [nn] -> x - x); 
   DisplayAuxError(e);

   if (e <= tol)
     return 1;


   phi0 = problem -> terminal [2] -> phi;

   if (it == 1) {
      x1 = problem -> terminal [2] -> phi;
      fx1 = (node [nn] -> x - problem -> terminal [2] -> x);

      if (analysis -> static_initial_guess == Catenary)
         problem -> terminal [2] -> phi = 0.5;
      else if (fx1 < 0)
         problem -> terminal [2] -> phi *= analysis -> outer_relaxation;
      else
         problem -> terminal [2] -> phi /= analysis -> outer_relaxation;
   }
   else if (it == 2) {
      x2 = problem -> terminal [2] -> phi;
      fx2 = (node [nn] -> x - problem -> terminal [2] -> x);

      problem -> terminal [2] -> phi = x2 - fx2*((x2 - x1)/(fx2 - fx1));
   }
   else {
      fx3 = (node [nn] -> x - problem -> terminal [2] -> x);

      if (fx3*fx1 < 0) {
         x2 = problem -> terminal [2] -> phi;
         fx2 = fx3;
      }
      else {
         x1 = problem -> terminal [2] -> phi;
         fx1 = fx3;
      }

      problem -> terminal [2] -> phi = x2 - fx2*((x2 - x1)/(fx2 - fx1));
   }

   phi = problem -> terminal [2] -> phi;

   DisplayMessage("x = %g, phi_o = %g, phi = %g", node [nn] -> x, phi0, phi);

   problem -> terminal [2] -> xforce = cos(phi)*problem -> terminal [2] -> tension;   
   problem -> terminal [2] -> yforce = sin(phi)*problem -> terminal [2] -> tension;   

   return 0;
}

int CheckDynstatConvergence(node, num_nodes, dt, scalv, ndof, twoD)
   Node		 *node;
   int		  num_nodes;
   double	  dt;
   double	 *scalv;
   int		  ndof;
{
   int		i;
   int		j;
   double	err1, err2;

   err1 = err2 = 0.0;

   if (problem -> type == Surface || problem -> type == Deployment) {
      if (twoD)
         err1 = (fabs(node[num_nodes] -> Y [3]) + fabs(node[num_nodes] -> Y[4]))/scalv [3];
      else 
         err1 = (fabs(node[num_nodes] -> Y [4]) 
                + fabs(node[num_nodes] -> Y [5])
                + fabs(node[num_nodes] -> Y [6]))/scalv [4];
   }

   for (i = 1 ; i <= num_nodes ; i++) 
      for (j = 1 ; j <= ndof ; j++) 
         err2 += fabs(node[i] -> Y[j] - node[i] -> Y_o[j])/scalv [j]/dt;
      
   err2 /= (num_nodes * ndof); 
   
   DisplayAuxError (err1 > err2 ? err1 : err2);

   if (err1 <= analysis -> static_tolerance && err2 <= analysis -> static_tolerance)
      return 1;
   else
      return 0;
}

int ResolveStartDepth(node, num_nodes, it)
   Node	    *node;
   int	     num_nodes;
   int	     it;
{
   static double x1, x2;
   static double fx1, fx2, fx3;
   static double prev_fx;
   double	target;
   double	x;
   double	err;
   double	thr;

   x = environment -> surface - node [1] -> x;
   /* target = problem -> terminal [1] -> profile.value; */
   target = Profile(0.0);
   thr = problem -> terminal [1] -> xthrust.value;
   
   err = fabs(x - target) / node [num_nodes] -> s;
   DisplayAuxError(err);
   if (err < analysis -> outer_tolerance)
      return 1;

   if (it == 1) {
      x1 = thr;
      fx1 = x - target;

      if (fx1 < 0)
         thr = -problem -> terminal [1] -> buoy -> Cl*problem -> terminal [2] -> yspeed.value*problem -> terminal [2] -> yspeed.value;
      else
         thr = problem -> terminal [1] -> buoy -> Cl*problem -> terminal [2] -> yspeed.value*problem -> terminal [2] -> yspeed.value;
   }
   else if (it == 2) {
      x2 = thr;
      fx2 = x - target;

      if (fx2 < 0 && fx1 < 0) {
         DisplayMessage("cannot reach target depth, max depth is %g", x);
         SetError(C_INSUFFICIENTLIFT);
         return -1;
      }
      else if (fx2 > 0 && fx1 > 0) {
	 DisplayMessage("cannot reach target depth, min depth is %g", x);
         SetError(C_INSUFFICIENTLIFT);
         return -1;
      }
         
      thr = x2 - fx2*(x2 - x1)/(fx2 - fx1); 

      prev_fx = fx1;
   }
   else {
      fx3 = x - target;

      if (fx3*fx1 < 0.0) {
         x2 = thr;
         fx2 = fx3;
         if (fx3*prev_fx > 0.0)
            fx1 = fx1/2.0;
      }
      else {
         x1 = thr;
         fx1 = fx3;
         if (fx3*prev_fx > 0.0)
            fx2 = fx2/2.0;
      }

      prev_fx = fx3;

      thr = x2 - fx2*(x2 - x1)/(fx2 - fx1);
   }

   problem -> terminal [1] -> xthrust.value = thr;
   return 0;
}

void
SmoothNodeData(double t, Node *active, int num_active, int twoD)
{
    int     i, k;
    double  a, oma;
    int     ne;

    ne = twoD ? 6 : 13;

    if (analysis -> dynamic_var_smooth.expr)
        a = EvalCode(analysis -> dynamic_var_smooth.expr, NULL, t, t,
                     0, 0, 0, CURRNODEDATA);
    else
        a = analysis -> dynamic_var_smooth.value;

    oma = 1.0 - a;

    for (i = 1 ; i <= num_active ; i++) {
        for (k = 1 ; k <= ne ; k++) 
            active[i] -> Y_f[k] = a*active[i] -> Y[k] + oma*active[i] -> Y_f[k];

        active[i] -> x_f = a*active[i] -> x + oma*active[i] -> x_f;
        active[i] -> y_f = a*active[i] -> y + oma*active[i] -> y_f;
        if (twoD)
            active[i] -> z_f = a*active[i] -> z + oma*active[i] -> z_f;

        active[i] -> xdot_f = a*active[i] -> xdot + oma*active[i] -> xdot_f;
        active[i] -> ydot_f = a*active[i] -> ydot + oma*active[i] -> ydot_f;
        if (twoD)
            active[i] -> zdot_f = a*active[i] -> zdot + oma*active[i] -> zdot_f;

        active[i] -> pay_f = a*active[i] -> pay + oma*active[i] -> pay_f;
    }            

    return;
}


