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

# ifndef _RESULTS_H
# define _RESULTS_H

# include "compress.h"

# define DISPLACEMENT   1
# define VELOCITY	2
# define FORCE		3
# define MOMENT		4
# define EULER		5

#ifndef MAXOUTPUT
#define MAXOUTPUT 11
#endif

typedef enum {
   drawStandard   = 1,
   drawHorizontal = 2,
   drawTowing     = 3
} DrawType;

typedef struct {
   DrawType	type;
   int	        global;
   int		totals;
   int		dynamic;
   int		twoD;
   int		nsteps;
   int		nsamples;
   int		nbuoy;
   int		npoints;
   int		main_points;
   int		num_branch;
   int      nseg;
   int      nseg_steps;
   int      n_ext;
   int      n_ext_steps;
   int	   *branch_starts;
   char	   *title;
   double	sample_dt;
   double	snap_dt;
   double	buoy_dt;
   double   seg_dt;
   double   ext_dt;
   int		num_output_nodes;
   int	   *output_nodes;
   double	depth;
   int		depth_ref;
   int		output_map [MAXOUTPUT];
   int      numvars;
   double     t_start;
   double      *s;
   double     *x;
   double     *y;
   double     *z;
   double     *velocity [3];
   double     *force [3] ;
   double     *moment [3];
   double     *beta [4];
   double      *x_st;
   double      *y_st;
   double      *z_st;
   double      *force_st [3];
   double      *moment_st [3];
   double      *beta_st [4];
   double     **x_t;
   double     **y_t;
   double     **z_t;
   double     **velocity_t [3];
   double     **force_t [3] ;
   double     **moment_t [3];
   double     **beta_t [4];
   double     **ext_t [3];
   double      *buoy [6];
   char        *input;
   unsigned long *snapTell;
   double     **seg_top_pay;
   double     **seg_bot_pay;
   double     **seg_stretched;
   double     **seg_unstretched;
   double     **seg_top_spooled;
   double     **seg_bot_spooled;
   double     **seg_first;
   double     **seg_last;
   ResFile     in;
} Result;


extern Result *ReadResultsFile (
   ResFile,
   int,			/* flag indicating local -> global transformation */
   int,         /* flag indicating to process dynamic + static    */
   int,         /* flag indicating N to lbs force conversion      */
   int          /* flag indicating m to ft length conversion      */
);

extern int ReadResultSnapshot(Result *, int, int, int, int);

extern int RemoveStatic (Result *);
extern int AddStatic (Result *);
extern double **array (int, int);
extern void free_array(double **, int, int);
extern char CheckMagic(char *); // -1 could not determine; 0 no; 1=results file

extern void ResultsToProblem(Result *res, Problem *p, Environment *e, int twoD);
extern void SetSnapGlobal(Result *);


# endif /* _RESULTS_H */
