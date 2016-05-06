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
 * File:	problem.h						*
 *									*
 * Description:	This file contains the public type and function		*
 *		declarations for the problem instance.			*
 ************************************************************************/

# ifndef _PROBLEM_H
# define _PROBLEM_H
# include "Tree.h"
# include "code.h"
# include "model.h"
# include "compress.h"

typedef struct {
    char    *input;
    char	*title;			/* problem title	*/
    char	*filename;		/* file name		*/
    Terminal	 terminal [3];		/* anchors and buoys    */
    Tree	 segment_tree;		/* segment tree		*/
    Tree	 branch_tree;		/* 			*/
    Tree	 branch_segment_tree;
    Tree	 buoy_tree;		/* buoy tree		*/
    Tree	 material_tree;		/* material tree	*/
    Tree	 anchor_tree;		/* anchor tree		*/
    Tree	 connector_tree;	/* connector tree	*/
    // Tree	 layout_tree;		/* layout tree		*/
    unsigned	 num_errors;		/* number of errors	*/
    unsigned	 line;			/* current line number	*/
    int		     type;			/* problem type		*/
    int          dynstat;		/* dynstat flag		*/
    Branch	    *branch;
    int		     num_branch;
    int		     junction_size;		/* max # of junct nodes */
    Node        *node;
    int          num_nodes;
    Node     *active;
    int      num_active;
    Segment  *segment;
    int      num_segments;
    int      twoD;
    int      dynamic;
    struct solution *solution;
} Problem;

#define OUTPUT_FIRST 1
#define OUTPUT_LAST 2
#define OUTPUT_TERMINALS 4
#define OUTPUT_CONNECTORS 8

typedef struct analysis {
    double	tolerance;
    double	outer_tolerance;
    double	static_tolerance;
    double	dynamic_tolerance;
    double	duration;
    double	dt;
    double	relaxation;
    unsigned	maxit;
    double	ramp_time;
    int		adapt_factor;
    int		adapt_levels;
    int		stall_limit;
    int		static_it;
    int		dynamic_it;
    int		outer_it;
    int         shooting_it;
    double	alpha_k;
    double	alpha_m;
    double	gamma;
    double	sp_rho;
    double	mesh_smooth;
    double	mesh_amplify;
    double	relax_up;
    double	relax_down;
    double	dynamic_relaxation;
    double	static_relaxation;
    double	outer_relaxation;
    int		current_steps;
    Integration integration;
    SolutionMethod static_initial_guess;
    SolutionMethod static_solution; 
    OuterMethod static_outer_method;
    struct solution *solution;
    double  viva_dt;
    int     viva_decimate;
    int     viva_iterations;
    VarExpr dynamic_var_smooth;
} Analysis;


typedef struct environment {
    InputType	   input_type;
    ForcingMethod  forcing;
    int		   num_components [4];
    double	   amplitude [4][10];
    double     period [4][10];
    double     phase [4][10];
    double	   omega [4][10];
    double	   k [4][10];
    double	   bottom_stiffness;
    double	   bottom_damping;
    double	   bottom_friction;
    char	  *velocity_file;
    char      *wave_file;
    char	  *current_file;
    double	   rho;
    double	   depth;
    double	   surface;
    double	   gravity;
    VarExpr	   bottom_elevation;
    VarExpr	   z_wind;
    VarExpr	   y_wind;
    VarExpr	   Ux;
    VarExpr	   Uz;
    VarExpr	   Uy;
    VarExpr	   Uz_mod;
    VarExpr	   Uy_mod;
    VarExpr	   Ux_mod;
    double     Uscale;
    double     Uxscale;
    double     Uyscale;
    double     Uzscale;
    double     current_rotation;
    double	   fill_density;
    struct solution *solution;
} Environment;

typedef struct solution {
    Problem     *problem;
    Environment *environment;
    Analysis    *analysis;
    int      static_only;
    int      twoD;
    char    *in_name;
    char    *out_name;
    char    *static_file;
    char    *initial_file;
    char    *dynstat_name;
    char    *progress_file;
    double   progress_dt;
    char    *restart_file;
    double   restart_t;
    double   output_dt;
    double   snapshot_dt;
    double   buoy_dt;
    double   segment_dt;
    double   ext_dt;
    int      output_map[11];
    int      decimate;
    int     *output_nodes;
    int      num_output_nodes;
    int      output_special;
    char    *tmpdir;
    int      debug_input;
    int      bill_of_materials; 
    int      refine;
    int      simple;
    int      auto_static;
    int      quiet;
    int      quit;
    int      unlink_input;
    int      X;
    int      saved;
    int      userQuit;
    int      solutionComplete;
    int      plotProgress;
    void    *controls;
    void    *snapPlot[11];
    void    *tsPlot[11];
    void    *table;
    void    *results;
    void    *playback;
    char    *table_name;
    double   rotation;
} Solution;

typedef struct {
    Problem     *problem;
    Analysis    *analysis;
    Environment *environment;
} ParserControl;

typedef struct debug {
   int	    status;
   double   snap_it;
   double   sample_it;
   ResFile  out;
   int     *output_map;
   int      num_output_nodes;
   int     *output_nodes;
   int	    decimate;
} Debug;


extern char      *copy_input	  (int);
extern void	  CableInitLexer	  (FILE *);
extern int	  Cable_parse	  (void);
extern int    ObjectCompare(Item, Item);
extern int	  ParseCppOptions (int *, char **);
extern int    SetupCpp(char *cpp, char **include, char **define);
extern int    ReadDatabaseFile(char *, Tree mat, Tree buo, Tree con, Tree anc);
extern int	  ReadModelFile	  (char *, Problem *, Analysis *, Environment *, int);
extern int	  WriteModelFile  (char *, Problem *, Analysis *, Environment *);
extern int	  DumpModelFile	  (char *, Problem *, Analysis *, Environment *);
extern int	  fWriteModelFile (FILE *, Problem *, Analysis *, Environment *);
extern int	  fDumpModelFile  (FILE *, Problem *, Analysis *, Environment *);
extern char  *ProblemTypeName(int);

extern void FillInObjectProperties (Problem *, Environment *);
extern void FillInEnvironment(Environment *);
extern void FillInAnalysis(Problem *, Analysis *);

extern void FreeEnvironment(Environment *);
extern void FreeProblem(Problem *);
extern void FreeAnalysis(Analysis *);
extern void FreeSolution(Solution *);

extern double Draft (
   double,
   Buoy,
   Environment *
);

extern double Buoyancy (
   Buoy,
   double,
   Environment *
);

extern void CreateNodeArray (
   Problem *p,
   Node **,
   int  *,      /* number of nodes     */
   Node **,
   int *        /* number of active nodes  */
);

extern int MakeNextPrevActive (
   Problem *,
   Node *,
   int,
   Node *
);


# endif /* _PROBLEM_H */
