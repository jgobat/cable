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
 * File:	solve.h							*
 *									*
 * Description:	This file contains the primary function			*
 *		declarations for the solution routines for the		*
 *		2- and 3-D static and dynamic mooring models		*
 ************************************************************************/

# ifndef _SOLVEL_H
# define _SOLVE_H
# include "model.h"

# define BLANK (-99999999)

typedef void (*func) ();

extern int SolveStaticProblem3D (
   int,			/* initial guess is loaded */
   Node *,		/* array of nodes	  */
   int,			/* number of nodes	  */
   Node *,
   int,			/* number of active nodes */
   ResFile,		/* output file pointer	  */
   int *   		/* output map	  	  */
);

extern int SolveStaticProblem2D (
   int,
   Node *,
   int,			/* number of nodes	  */
   Node *,		// array of active nodes
   int,			/* number of active nodes */
   ResFile,		/* output file pointer	  */
   int *  		/* output map	  	  */
);

extern int ShootStaticProblem2D (
   int,
   Node *,		/* array of nodes	  */
   int,			/* number of nodes	  */
   Node *,		// array of active nodes
   int,			/* number of active nodes */
   ResFile,		/* output file pointer	  */
   int *  		/* output map	  	  */
);

extern int SolveDynamicProblem3D (
   Node *,		/* array of nodes	       */
   int,			/* number of nodes	       */
   Node *,		// array of active nodes
   int,			/* number of active nodes */
   ResFile,		/* output file pointer         */
   int *,		/* output map		       */
   int *,		/* array of output node number */
   int,			/* number of output nodes      */
   double,		/* sample rate		       */
   double,		/* snapshot increment  	       */
   double,      /* segment increment */
   double,		/* buoy motion increment       */
   double,      /* external force increment */
   int,			/* decimation factor	       */
   char *,      // progress save file
   double,      // progress time step
   char *,      // restart read file
   double       // progress restart time 
);

extern int ShootStaticProblem3D (
   int,
   Node *,		/* array of nodes	  */
   int,			/* number of nodes	  */
   Node *,		// array of active nodes
   int,			/* number of active nodes */
   ResFile,		/* output file pointer	  */
   int *  		/* output map	  	  */
);

extern int SolveDynamicProblem2D (
   Node *,		/* array of nodes	       */
   int,			/* number of nodes	       */
   Node *,		// array of active nodes
   int,			/* number of active nodes */
   ResFile,		/* output file pointer         */
   int *,		/* output map		       */
   int *,		/* array of output node number */
   int,			/* number of output nodes      */
   double,		/* sample rate		       */
   double,		/* snapshot increment  	       */
   double,      /* segment increment */
   double,		/* buoy motion increment       */
   double,      /* external force increment */
   int,		    /* decimation factor	       */
   char *,      // progress save file
   double,      // progress time step
   char *,      // restart read file
   double       // progress restart time 
);

extern int AutoStaticSolve (
   Node *,
   int,
   Node *,		// array of active nodes
   int,
   int,
   char *
);

extern int SolveDE (
   void (*) ( ),	/* difeq function (static and dynami)		*/
   void (*) ( ),	/* coordinate update func (static and dynamic)	*/
   int *,      		/* itmax (static and dynamic)			*/
   double *,   		/* conv (static and dynamic)			*/
   double *,   		/* base_relax (static and dynamic)		*/
   double *,   		/* scalv (static and dynamic)			*/
   int,       		/* ne (static and dynamic)			*/
   int,       		/* nb (static and dynamic)			*/
   int,			/* nb_branch (static and dynamic)		*/
   int,			/* nj_compat (static and dynamic)		*/
   Node *,
   int,
   Node *,     		/* *node (static and dynamic) (active only)	*/
   int,       		/* m (static and dynamic) (active only)		*/
   double **,  		/* s (static and dynamic)			*/
   double,   		/* tm (dynamic)					*/
   double,    		/* dt (dynamic)					*/
   double,		/* current_factor (static)			*/
   int 			/* outer iteration step (static)		*/
);

extern void EquivalentWeight (
   Segment *,		/* segment array	    */
   int,			/* number of segments       */
   double *,		/* address of w0 result     */
   double *,		/* address of EA result     */
   double *		/* address of length result */
);

extern void CatenaryForces (
   double,		/* weight per length		       */
   double,		/* stiffness EA			       */
   double,		/* length			       */
   double,		/* vertical offset of top	       */
   double,		/* horizontal offset of top	       */
   double *,		/* address for horizontal force return */
   double *		/* address for vertical force return   */
);

extern int InitialGuess (
   void (*) ( ),	/* catenary solution function	*/
   Node	*,		/* array of nodes		*/
   int 			/* number of nodes		*/
);

extern void SolveCatenary2D (
   Node,	// start node
   double, 	// x0
   double,   	// y0
   double,	// z0
   double,     	// hforst 
   double,     	// xforst
   double,    	// th
   double,	// w0
   double    	// length
);
   
extern void SolveCatenary3D (
   Node,
   double,           
   double,          
   double,
   double,      
   double,     
   double,    
   double,
   double    
);
   
extern void StaticDifeq2D (
   EquationType,/*    eq_type,			*/
   Node,	/*    n				*/
   Node,	/*    nm			*/
   int,		/*    ne,			*/
   int,        	/*    rhs,			*/
   int,        	/*    num_rows,			*/
   double **,  	/*  **s,			*/
   Node  *,  	/*   *node,			*/
   double,     	/*    tm,			*/
   double,    	/*    dt,			*/
   double 	/*    current factor		*/
);

extern void DynamicDifeq2D (
   EquationType,/*    eq_type,			*/
   Node,       	/*    n,			*/
   Node,	/*    nm,			*/
   int,		/*    ne,			*/
   int,        	/*    rhs,			*/
   int,        	/*    num_rows,			*/
   double **,  	/*  **s,			*/
   Node  *,  	/*   *node,			*/
   double,     	/*    tm,			*/
   double,    	/*    dt,			*/
   double 	/*    current factor		*/
);


extern void StaticDifeq3D (
   EquationType,/*    eq_type,			*/
   Node,
   Node,
   int,		/*    ne,			*/
   int,        	/*    rhs,			*/
   int,        	/*    num_rows,			*/
   double **,  	/*  **s,			*/
   Node  *,  	/*   *node,			*/
   double,     	/*    tm,			*/
   double,    	/*    dt,			*/
   double 	/*    current factor		*/
);

extern void DynamicDifeq3D (
   EquationType,/*    eq_type,			*/
   Node,
   Node,
   int,		/*    ne,			*/
   int,        	/*    rhs,			*/
   int,        	/*    ne,			*/
   double **,  	/*  **s,			*/
   Node  *,  	/*   *node,			*/
   double,     	/*    tm,			*/
   double,    	/*    dt,			*/
   double 	/*    current factor		*/
);

extern void StaticUpdate3D (
   Node *,
   int,
   double,
   double
);

extern void DynamicUpdate3D (
   Node *,
   int,
   double,
   double
);

extern void StaticUpdate2D (
   Node *,
   int, 
   double,
   double
);

extern void DynamicUpdate2D (
   Node *,
   int,
   double,
   double
);

void BendingDerivatives2D (
   Node
);

void BendingDerivatives3D (
   Node
);

extern double check (
   double,	/* time		*/
   double,	/* time step	*/
   double	/* sample rate	*/
);

extern int ResolveWebster (
   int,		/* iteration number	*/
   Node *,
   int,		/* number of nodes	*/
   int		/* twoD flag		*/
);

extern int ResolvePositionedConnectors(Segment *, int, double, int);

extern int ResolveAnchor (
   Node *,	/* array of nodes	*/
   int,		/* number of nodes	*/
   Segment *,
   int,
   int		/* twoD flag		*/
);

extern int ResolveSpeed (
   Node *,	/* array of nodes	*/
   int,		/* number of nodes	*/
   int 		/* step number	   	*/
);

extern int ResolveSpeed3D (
   Node *,	/* array of nodes	*/
   double **,	/* y information array  */
   int,		/* number of nodes	*/
   int 		/* step number	   	*/
);

extern double ResolveBuoy (
   Node *,	/* array of nodes	 */
   int,		/* number of nodes	 */
   int,		/* outer iteration count */
   int		/* twoD flag		 */
);

extern int ResolveStartDepth (
   Node *,	/* array of nodes	 */
   int,		/* number of nodes	 */
   int 		/* outer iteration count */
);

extern int ResolvePayRate (
   double,	/* current time		  */
   double,	/* time step		  */
   Node *,	/* array of nodes	  */
   int, 	/* number of active nodes */
   int 		/* number of eqs in y     */
);

extern double sign (
   double
);

extern double dispersion (
   double, 
   double, // gravity
   double  // depth
);

extern int CheckDynstatConvergence (
   Node *,
   int,
   double,
   double *,
   int,
   int
);

extern int SimpleDynamics (
   Node *,
   int
);

extern int UpdateLAMP (
   Node *,			/* array of nodes		*/
   double **,			/* y information array		*/
   int,				/* number of nodes 		*/
   double,			/* current WHOI Cable time	*/
   ResFile,			/* output file			*/
   double,			/* buoy motion increment 	*/
   int				/* twoD solution flag		*/
);

extern int ProcessSpools(Problem *p, Node *node, int num_nodes, 
                         Node *active, int num_active,
                         double t, double dt, int ne, int twoD);


extern void SmoothNodeData(double t, Node *active, int num_active, int twoD);

# endif /* _SOLVE_H */
