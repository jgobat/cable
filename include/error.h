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
 * File:	error.h							*
 *									*
 * Description:	This file contains the function declarations for the	*
 *		error handling routines.				*
 ************************************************************************/

# ifndef _LOCAL_ERROR_H
# define _LOCAL_ERROR_H
#include "problem.h"

extern void error (char *, ...);
extern void Fatal (char *, ...);
extern void ExitErr (char *, ...);
extern void Exit ( int );
extern void DisplayMessage (char *, ...);

extern void BufferErrors ( int );
extern void SetTempDir ( char * );
extern void QueryX ( char * );
extern void DisplayAuxError ( double );

extern void DisplayInfo (
   int,         /* step number          */
   double,      /* time                 */
   int,         /* iteration            */
   double,      /* error                */
   double       /* relaxation factor    */	
);

extern void DisplayMode (
   int,         /* mode (0 == ascii, 1 = graphical)      */
   int          /* quiet (0 == only errors, 1 = verbose) */
);

extern void DisplayStaticHeader ( );

extern void DisplayDynamicHeader ( );


extern void RefreshDisplay (Analysis *);

extern void SetError ( int );
extern int GetError ( );

#define C_INITIALSOLUTIONFAILED  1
#define C_STATICSOLUTIONFAILED   2
#define C_MAXITERATIONSEXCEEDED  3
#define C_SINGULARITY	         4
#define C_COMPRESSION            5
#define C_NANINF		 6
#define C_ITERATIONSTALLED	 7
#define C_BUOYNOTATSURFACE	 8
#define C_BUOYNOTINWATER         9
#define C_NOSHOOTING		10
#define C_MAXADAPTEXCEEDED      11
#define C_CATENARYSOLUTIONFAILED 12
#define C_SHOOTINGSOLUTIONFAILED 13
#define C_INSUFFICIENTLIFT	14
#define C_NOMORENODES 15

# endif /* _ERROR_H */
