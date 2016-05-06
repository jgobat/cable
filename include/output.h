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
 * File:	output.h    
 *	
 * Description:	
 ************************************************************************/

# ifndef _OUTPUT_H
# define _OUTPUT_H

# define MOTION	1
# define VEL	2
# define FORCE	3
# define MOMENT	4
# define EULER	5
# define EXTERNAL 6

extern int InitializeResultsFile (
   ResFile,		/* output file pointer			*/
   int,			/* number of nodes in problem		*/
   char *,		/* problem title			*/
   int *,		/* output map				*/
   int,			/* decimation factor			*/
   int,			/* dynamic solution present flag 	*/
   int,			/* 2D solution flag			*/
   char *       /* string containing input file text */
);

extern int WriteStaticSolution (
   Node *,		/* node array		*/
   int,			/* number of nodes	*/
   ResFile,		/* output file pointer	*/
   int *,		/* output map		*/
   int,			/* decimation factor 	*/
   int			/* 2D y array flag	*/
);

extern int WriteDynamicSegmentData (
   Problem *,
   ResFile
);

extern void WriteDynamicExternalForces (
   Problem *,
   ResFile
);

extern int WriteDynamicHeader (
   ResFile,		/* output file pointer	        */
   double,      /* t start */
   double,		/* problem duration (t end)		*/
   double,		/* problem dt		        */
   double,		/* output sample dt	        */
   double,		/* snapshot dt		        */
   double,      /* segment dt                   */
   double,		/* buoy motion dt		*/
   double,      /* external force dt */
   int,			/* number of output nodes       */
   int *,		/* array of output node numbers */
   int,			/* decimation factor		*/
   Node *		/* array of nodes		*/
);

extern int WriteBuoyMotion (
   ResFile,		/* output file pointer		*/
   double *		/* vector of 6-axis buoy motion */
);

extern int WriteDynamicSnapshot (
   ResFile,		/* output file pointer	        */
   int *,		/* output information map	*/
   Node *,		/* array of nodes		*/
   int,			/* number of nodes		*/
   int,			/* decimation factor		*/
   int			/* 2D y array flag	*/
);

extern int FakeDynamicSnapshot (
   ResFile,		/* output file pointer	        */
   Node *,		/* array of nodes		*/
   int,			/* number of nodes		*/
   int,			/* decimation factor		*/
   int			/* 2D y array flag	        */
);

extern int WriteDynamicResult (
   ResFile,		/* output file pointer	        */
   int *,		/* output information map  	*/
   int *,		/* array of nodes		*/
   int,			/* number of nodes		*/
   int,			/* decimation factor		*/
   Node *,		/* array of nodes		*/
   int, 
   int			/* 2D y array flag	*/
);

extern int LoadStaticSolution (
   char	*,		/* name of file to load solution from         */
   Node *,		/* array of nodes			      */
   int,			/* number of nodes in problem (sanity check)  */
   int,			/* 2D y array flag			      */
   double       /* theta rotation 2D to 3D */
);

extern void FlushOutput ( );

extern void InitializeProgressFile(Problem *, ResFile, int);
extern void WriteProgress(double, Problem *, ResFile, int);
extern long LoadProgressPoint(char *, Problem *, double, int);

# endif	/* _OUTPUT_H */
