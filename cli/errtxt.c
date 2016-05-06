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
 * File:	error.c							*
 *									*
 * Description:	This file contains the function definitions for the	*
 *		error handling routines.				*
 ************************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <stdarg.h>
# include "error.h"


/************************************************************************
 * Function:	error							*
 *									*
 * Description:	Prints an error message specified as a format string	*
 *		and arguments to standard error.  If the current line	*
 *		exists (is not zero) then current file and line are	*
 *		printed before the error message.			*
 ************************************************************************/

void 
error (char *format, ...)
{
    va_list ap;


    va_start (ap, format);

    fprintf (stderr, "res2mat: ");

    vfprintf (stderr, format, ap);
    fprintf (stderr, "\n");
    va_end (ap);
}


/************************************************************************
 * Function:	Fatal							*
 *									*
 * Description:	Prints an error message specified as a format string	*
 		and arguments to standard error and exits the program.	*
 ************************************************************************/

void 
Fatal (char *format, ...)
{
    va_list ap;


    va_start (ap, format);
    fprintf (stderr, "res2mat: ");
    vfprintf (stderr, format, ap);
    fprintf (stderr, "\n");
    va_end (ap);
    exit (1);
}
