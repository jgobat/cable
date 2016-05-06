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
# include <malloc.h>
# include <unistd.h>

# include "error.h"
# include "compress.h"
# include "problem.h"
# include "control.h"
# include <stdarg.h>

extern int    unlink_input;
extern int    static_finished;
extern char  *out_name;
extern char  *in_name;

extern Problem *problem;

# ifdef WINGUI

# include <windows.h>

# define False FALSE
# define True TRUE
# define Boolean WINBOOL

# endif

# ifdef GUI

# include <X11/Intrinsic.h>
# include "OutputDialog.h"
# include "util.h"

extern Widget	 	toplevel;
extern XtAppContext	app_context;
static OutputDialog	error_dialog = NULL;
static String 		buttons [] = {"dismiss"};

# endif

# if defined (GUI) || defined (WINGUI)

static FILE	       *output;
static Boolean   	buffer_errors = False;
static Boolean   	errors = False;

# endif

static int	graphic_display = 0;
static int	display;



/************************************************************************
 * Function:	error							*
 *									*
 * Description:	Prints an error message specified as a format string	*
 *		and arguments to standard error.  If the current line	*
 *		exists (is not zero) then current file and line are	*
 *		printed before the error message.			*
 ************************************************************************/

void error (char *format, ...)
{
# ifdef WINGUI
    char buffer [1024];
# endif
    va_list ap;

# ifdef GUI
    if (error_dialog == NULL && graphic_display)
       error_dialog = OutputDialogCreate (toplevel, "errorDialog",
					  buttons, XtNumber(buttons));
# endif


    va_start (ap, format);

# if !defined (GUI) && !defined (WINGUI)

    if (problem -> line)
	fprintf (stderr, "%s:%d: ", problem -> filename, problem -> line);

    vfprintf (stderr, format, ap);
    fprintf (stderr, "\n");
    va_end (ap);

# else
    if (graphic_display) {
        if (buffer_errors == True) {
            if (problem -> line)
                fprintf (output, "%s:%d: ", problem -> filename, problem -> line); 
            vfprintf (output, format, ap);

            fprintf (output, "\n");

            va_end (ap);
            errors = True;
        } 
        else {
# ifdef GUI
            OutputDialogVprintf (error_dialog, format, ap);
            CenterOnWidget (OutputDialogShell (error_dialog), toplevel, True); 
            WarpToCenter (OutputDialogShell (error_dialog)); 
            OutputDialogSelect (error_dialog, "Error", "dismiss");
# elif defined (WINGUI)
            vsprintf (buffer, format, ap);
            WinErrorDialog(buffer);
# endif
            va_end (ap);
        }
    }
    else {
        if (problem -> line)
            fprintf (stderr, "%s:%d: ", problem -> filename, problem -> line);

        vfprintf (stderr, format, ap);
        fprintf (stderr, "\n");
        va_end (ap);
    }
# endif

    problem -> num_errors ++;
}


/************************************************************************
 * Function:	ExitErr							*
 *									*
 * Description:	Prints an error message specified as a format string	*
 * 		and arguments to standard error or to the error dialog  *
 *              and exits the program.					*
 ************************************************************************/

void ExitErr (char *format, ...)
{
# ifdef WINGUI
    char    buffer [1024];
# endif
    va_list ap;

# ifdef GUI
    if (error_dialog == NULL && graphic_display)
       error_dialog = OutputDialogCreate (toplevel, "errorDialog",
					  buttons, XtNumber(buttons));
# endif


    va_start (ap, format);

# if !defined (GUI) && !defined (WINGUI)

    fprintf (stderr, "fatal: ");
    vfprintf (stderr, format, ap);
    fprintf (stderr, "\n");

# else

    if (graphic_display) {

# ifdef GUI

       OutputDialogVprintf (error_dialog, format, ap);
       CenterOnWidget (OutputDialogShell (error_dialog), toplevel, True); 
       WarpToCenter (OutputDialogShell (error_dialog)); 
       OutputDialogSelect (error_dialog, "Error", "dismiss");

# elif defined WINGUI

       vsprintf (buffer, format, ap);
       WinErrorDialog(buffer);

# endif

    }
    else {
       fprintf (stderr, "fatal: ");
       vfprintf (stderr, format, ap);
       fprintf (stderr, "\n");
    }
# endif

    va_end (ap);

    if (unlink_input)
       unlink (in_name);

    exit (1);
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

    exit (1);
}

#ifdef _CONVEX_SOURCE
# define TEMPNAME(a,b) tmpnam((b))
#else
# define TEMPNAME(a,b) tempnam((a), (b))
#endif

static char *tmpdir = NULL;

void SetTempDir (name)
   char *name;
{
   tmpdir = name;
}

void BufferErrors (flag)
    int flag;
{
# if defined (GUI) || defined (WINGUI)

    char name[] = "cablXXXXXX";

# ifdef GUI
    static OutputDialog	  output_dialog = NULL;
    static String	  buttons [] = {"dismiss"};

    if (output_dialog == NULL)
       output_dialog = OutputDialogCreate (toplevel, "outputDialog",
                                           buttons, XtNumber (buttons));
# endif
    
    if (flag == True && buffer_errors == False) {
        if ((output = fdopen (mkstemp(name), "w")) == NULL) {
            error ("Could not create temporary file for output.");
            return;
        }
        buffer_errors = flag;
        errors = False;

    } 
    else if (flag == False && buffer_errors == True) {
        fclose (output);
        
        if (errors == True) {

# ifdef GUI
            OutputDialogView (output_dialog, name, 10, 60);
            CenterOnScreen (OutputDialogShell (output_dialog), False);
            WarpToCenter (OutputDialogShell (output_dialog));
            OutputDialogSelect (output_dialog, "Error", "dismiss");
# elif defined (WINGUI)
            WinErrorFileDialog (name);
# endif

        }
        
        
        buffer_errors = flag;
        (void) unlink (name);

    }
# endif
}   

void DisplayMode (mode, quiet)
   int 	mode;
   int	quiet;
{
   display = !quiet;
   graphic_display = mode;
}

void DisplayAuxError (err)
   double err;
{
   if (!display)
      return;

# if defined (GUI) || defined (WINGUI)
   if (graphic_display)
      ControlAuxError (err, NULL);
# endif
} 

void DisplayInfo (step, tm, it, err, rlx)
   int		step;
   double	tm;
   int         	it;
   double	err;
   double	rlx;
{
   if (!display)
      return;

# if defined (GUI) || defined (WINGUI)

   if (graphic_display)
      if (step)
         ControlInfo ((double) step, it, err, NULL);
      else
         ControlInfo (tm, it, err, NULL);
   else if (step)
      fprintf (stdout,"%4d      %4d       %11.5g    %11.5g\n", 
               step, it, err, rlx);
   else
      fprintf (stdout,"%7.4f    %4d      %11.5g    %11.5g\n", 
               tm, it, err, rlx);

# else

   if (step)
      fprintf (stdout,"%4d      %4d       %11.5g    %11.5g\n", 
               step, it, err, rlx);
   else
      fprintf (stdout,"%7.4f    %4d      %11.5g    %11.5g\n", 
               tm, it, err, rlx);

# endif

   return;
}

void DisplayStaticHeader ( )
{
   if (!display)
      return;

# if defined (GUI) || defined (WINGUI)
   if (graphic_display)
      ControlLabels ("step", "iterations", "error", "aux err", NULL);
   else
      fprintf (stdout," step    iterations   total error      relaxation\n");
# else
   fprintf (stdout," step    iterations   total error      relaxation\n");
# endif
}

void DisplayDynamicHeader ( )
{
   if (!display)
      return;

# if defined (GUI) || defined (WINGUI)
   if (graphic_display)
      ControlLabels ("time", "iterations", "error", "aux err", NULL);
   else
      fprintf (stdout," time    iteration   total error\n");
# else
   fprintf (stdout," time    iteration   total error\n");
# endif
}

void DisplayMessage (char *format, ...)
{
   char	   msg [1024];
   va_list ap;

   if (!display)
      return;

   va_start (ap, format);
   vsprintf (msg, format, ap);
   va_end (ap);
  
# if defined (GUI) || defined (WINGUI)
   if (graphic_display) {
      ControlMessage (msg, NULL);
      ControlProcessEvents (problem -> solution);
   }
   else
      fprintf (stdout,"--> %s\n", msg);
# else
   fprintf (stdout,"--> %s\n", msg);
# endif
}

void RefreshDisplay (Analysis *a)
{
   if (!display)
      return;

# if defined (GUI) || defined (WINGUI)
   if (graphic_display) 
      ControlDialogInitialize (a, NULL, 0);
# endif

   return;
}
   
void Exit (status)
   int	status;
{
   if (unlink_input)
      unlink (in_name);

   if (!static_finished && out_name)
      unlink (out_name);

# if defined (GUI)
   if (graphic_display && display)
      XtAppMainLoop (app_context);
# elif defined (WINGUI)
   if (graphic_display && display)
      ControlMainLoop();
# endif
 
   exit (status);
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

void ControlDialogQuit(Solution *s)
{
    return;
}

