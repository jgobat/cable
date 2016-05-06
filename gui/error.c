/*
    This file is part of WHOI Cable, a program for the static and dynamic
    analysis of oceanographic cable structures.

    Copyright (C) 2008-2016 by Jason Gobat

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
# include <malloc.h>
# include <unistd.h>

# include "localerror.h"
# include "compress.h"
# include "problem.h"
# include "localcontrol.h"
# include <stdarg.h>

extern int    unlink_input;
extern int    static_finished;
extern char  *out_name;
extern char  *in_name;

extern Problem problem;

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
    char buffer [1024];
    va_list ap;

    va_start (ap, format);

    vsprintf (buffer, format, ap);
    ControlMessage(buffer, NULL);
    vfprintf(stderr, format, ap);
    va_end (ap);


    problem.num_errors ++;
}


/************************************************************************
 * Function:	ExitErr							*
 *									*
 * Description:	Prints an error message specified as a format string	*
 * 		and arguments to standard error or to the error dialog  *
 *              and exits the program.					*
 ************************************************************************/

void 
ExitErr (char *format, ...)
{
    char    buffer [1024];
    va_list ap;

    va_start (ap, format);
    vsprintf (buffer, format, ap);
    ControlMessage(buffer, NULL);
    va_end (ap);


    if (unlink_input)
       unlink (in_name);

    abort();
}

/************************************************************************
 * Function:	Fatal							*
 *									*
 * Description:	Prints an error message specified as a format string	*
 		and arguments to standard error and exits the program.	*
 ************************************************************************/

void Fatal (char *format, ...)
{
    va_list ap;


    va_start (ap, format);
    fprintf (stderr, "fatal: ");
    vfprintf (stderr, format, ap);
    fprintf (stderr, "\n");
    va_end (ap);

    if (unlink_input)
       unlink (in_name);

    abort();
}

static char *tmpdir = NULL;

void SetTempDir (name)
   char *name;
{
   tmpdir = name;
}

void 
BufferErrors (int flag)
{
}   

void 
DisplayMode (int mode, int quiet)
{
}

void 
DisplayAuxError (double err)
{
    ControlAuxError (err, NULL);
} 

void DisplayInfo (step, tm, it, err, rlx)
   int		step;
   double	tm;
   int         	it;
   double	err;
   double	rlx;
{
    if (step)
        ControlInfo ((double) step, it, err, NULL);
    else
        ControlInfo (tm, it, err, NULL);

   return;
}

void DisplayStaticHeader ( )
{
    ControlLabels ("step", "iterations", "error", "aux err", NULL);
}

void DisplayDynamicHeader ( )
{
    ControlLabels ("time", "iterations", "error", "aux err", NULL);
}

void DisplayMessage (char *format, ...)
{
   char	   msg [1024];
   va_list ap;

   va_start (ap, format);
   vsprintf (msg, format, ap);
   va_end (ap);
  
   ControlMessage (msg, NULL);
   ControlProcessEvents (NULL);
}

void 
RefreshDisplay (Analysis *a)
{
    ControlDialogInitialize (a, NULL, 0);
    return;
}
   
void 
Exit(int status)
{
   if (unlink_input)
      unlink (in_name);

   if (!static_finished && out_name)
      unlink (out_name);

   // main loop
 
   abort();
}

static int error_code = 0;
 
void SetError(code)
   int  code;
{
   error_code = code;
}

int GetError ( )
{
   int	code;
  
   code = error_code;
   error_code = 0;   

   return code;
}
