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

/**************************************************************************** * 
 * File:	WinControl.c
 *
 * Description:	Contains routines to offer graphical solution controls
 *		in a 32-bit Windows (95/NT) environments
 *
 ****************************************************************************/

# ifdef WINGUI

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <sys/types.h>
# include <sys/stat.h>
# include <unistd.h>
# include <fcntl.h>
# include <string.h>
# include <windows.h>
# include "compress.h"
# include "allocate.h"
# include "error.h"
# include "problem.h"
# include "solve.h"


extern int	 unlink_input;
extern char	*in_name;
extern int	 static_finished;
extern char	*out_name;

extern Analysis analysis;

static double	 base_static_relaxation;
static double	 base_outer_relaxation;
static double	 base_dynamic_relaxation;

static double	 base_static_tolerance;
static double	 base_outer_tolerance;
static double	 base_dynamic_tolerance;

static int	 base_static_iterations;
static int	 base_outer_iterations;
static int	 base_dynamic_iterations;

static double	 base_duration;
static double    base_time_step;
static double    base_ramp_time;

static HWND cd = NULL;

extern int CableMain();
// extern int isspace();

static void SetSensitivity (hControld, state)
   HWND	     hControld;
   WINBOOL   state;
{
    EnableWindow (GetDlgItem(hControld, 100), state);
    EnableWindow (GetDlgItem(hControld, 101), state);
    EnableWindow (GetDlgItem(hControld, 102), state);

    EnableWindow (GetDlgItem(hControld, 200), state);
    EnableWindow (GetDlgItem(hControld, 201), state);
    EnableWindow (GetDlgItem(hControld, 202), state);

    EnableWindow (GetDlgItem(hControld, 300), state);
    EnableWindow (GetDlgItem(hControld, 301), state);
    EnableWindow (GetDlgItem(hControld, 302), state);

    EnableWindow (GetDlgItem(hControld, 400), state);
    EnableWindow (GetDlgItem(hControld, 401), state);
    EnableWindow (GetDlgItem(hControld, 402), state);

    EnableWindow (GetDlgItem(hControld, 401), FALSE);

    EnableWindow (GetDlgItem(hControld, 5002), state);
    EnableWindow (GetDlgItem(hControld, 5003), state);

    if (!state)
       SetFocus (GetDlgItem(hControld, 5001));
    else
       SetFocus (GetDlgItem(hControld, 100));

    return;
}

static void Restore (hControld)
   HWND  hControld;
{
   char	 buffer [32];

   analysis.static_relaxation = base_static_relaxation;
   sprintf (buffer, "%g", analysis.static_relaxation);
   SetDlgItemText (hControld, 100, buffer);

   analysis.outer_relaxation = base_outer_relaxation;
   sprintf (buffer, "%g", analysis.outer_relaxation);
   SetDlgItemText (hControld, 101, buffer);

   analysis.dynamic_relaxation = base_dynamic_relaxation;
   sprintf (buffer, "%g", analysis.dynamic_relaxation);
   SetDlgItemText (hControld, 102, buffer);

   analysis.static_tolerance = base_static_tolerance;
   sprintf (buffer, "%g", analysis.static_tolerance);
   SetDlgItemText (hControld, 200, buffer);

   analysis.outer_tolerance = base_outer_tolerance;
   sprintf (buffer, "%g", analysis.outer_tolerance);
   SetDlgItemText (hControld, 201, buffer);

   analysis.dynamic_tolerance = base_dynamic_tolerance;
   sprintf (buffer, "%g", analysis.dynamic_tolerance);
   SetDlgItemText (hControld, 202, buffer);

   analysis.static_it = base_static_iterations;
   sprintf (buffer, "%d", analysis.static_it);
   SetDlgItemText (hControld, 300, buffer);

   analysis.outer_it = base_outer_iterations;
   sprintf (buffer, "%d", analysis.outer_it);
   SetDlgItemText (hControld, 301, buffer);

   analysis.dynamic_it = base_dynamic_iterations;
   sprintf (buffer, "%d", analysis.dynamic_it);
   SetDlgItemText (hControld, 302, buffer);

   analysis.duration = base_duration;
   sprintf (buffer, "%g", analysis.duration);
   SetDlgItemText (hControld, 401, buffer);

   analysis.ramp_time = base_ramp_time;
   sprintf (buffer, "%g", analysis.ramp_time);
   SetDlgItemText (hControld, 402, buffer);

   analysis.dt = base_time_step;
   sprintf (buffer, "%g", analysis.dt);
   SetDlgItemText (hControld, 400, buffer);

   return;
}

static void Change (hControld)
   HWND  hControld;
{
   char		  ptr [32];
   char		  buffer [32];
   double	  value;
   int		  i;

   GetDlgItemText (hControld, 100, ptr, 32);
   value = atof (ptr);
   if (value > 0.0 && value <= 1.0)
      analysis.static_relaxation = value;
   else {
      sprintf (buffer, "%g", analysis.static_relaxation);
      SetDlgItemText (hControld, 100, buffer);
   }

   GetDlgItemText (hControld, 102, ptr, 32);
   value = atof (ptr);
   if (value > 0.0 && value <= 1.0)
      analysis.dynamic_relaxation = value;
   else {
      sprintf (buffer, "%g", analysis.dynamic_relaxation);
      SetDlgItemText (hControld, 102, buffer);
   }

   GetDlgItemText (hControld, 101, ptr, 32);
   value = atof (ptr);
   if (value > 0.0)
      analysis.outer_relaxation = value;
   else {
      sprintf (buffer, "%g", analysis.outer_relaxation);
      SetDlgItemText (hControld, 101, buffer);
   }

   GetDlgItemText (hControld, 300, ptr, 32);
   i = atoi(ptr);
   if (i)
      analysis.static_it = i;
   else {
      sprintf (buffer, "%d", analysis.static_it);
      SetDlgItemText (hControld, 300, buffer);
   }

   GetDlgItemText (hControld, 302, ptr, 32);
   i = atoi(ptr);
   if (i) 
      analysis.dynamic_it = i;
   else {
      sprintf (buffer, "%d", analysis.dynamic_it);
      SetDlgItemText (hControld, 302, buffer);
   }

   GetDlgItemText (hControld, 301, ptr, 32);
   i = atoi(ptr);
   if (i > 0)
      analysis.outer_it = i;
   else {
      sprintf (buffer, "%d", analysis.outer_it);
      SetDlgItemText (hControld, 301, buffer);
   }
      
   GetDlgItemText (hControld, 200, ptr, 32);
   value = atof(ptr);
   if (value)
      analysis. static_tolerance = value;
   else {
      sprintf (buffer, "%g", analysis.static_tolerance);
      SetDlgItemText (hControld, 200, buffer);
   }

   GetDlgItemText (hControld, 202, ptr, 32);
   value = atof(ptr);
   if (value)
      analysis.dynamic_tolerance = value;
   else {
      sprintf (buffer, "%g", analysis.dynamic_tolerance);
      SetDlgItemText (hControld, 202, buffer);
   }

   GetDlgItemText (hControld, 201, ptr, 32);
   value = atof(ptr);
   if (value)
      analysis.outer_tolerance = value;
   else {
      sprintf (buffer, "%g", analysis.outer_tolerance);
      SetDlgItemText (hControld, 201, buffer);
   }

   GetDlgItemText (hControld, 401, ptr, 32);
   value = atof(ptr);
   if (value)
      analysis.duration = value;
   else {
      sprintf (buffer, "%g", analysis.duration);
      SetDlgItemText (hControld, 401, buffer);
   }

   GetDlgItemText (hControld, 402, ptr, 32);
   value = atof(ptr);
   if (value)
      analysis.ramp_time = value;
   else {
      sprintf (buffer, "%g", analysis.ramp_time);
      SetDlgItemText (hControld, 402, buffer);
   }

   GetDlgItemText (hControld, 400, ptr, 32);
   value = atof(ptr);
   if (value)
      analysis.dt = value;
   else {
      sprintf (buffer, "%g", analysis.dt);
      SetDlgItemText (hControld, 400, buffer);
   }

   return;
}


int 
ControlProcessEvents (Solution *s)
{
   MSG	msg;

   if (cd) {
      while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE)) {
         if (!IsDialogMessage (cd, &msg)) {
            TranslateMessage (&msg);
            DispatchMessage (&msg);
         } 
      }
   }

   return;
}

int 
ControlDialogInitialize (Analysis *a, void *controls)
{
    HWND	     hControld;
    char	     buffer [32];

    hControld = cd;

	/*
	 * initialize the controls
	 */

    sprintf (buffer, "%g", a -> static_relaxation);
    SetDlgItemText (hControld, 100, buffer);
    sprintf (buffer, "%g", a -> outer_relaxation);
    SetDlgItemText (hControld, 101, buffer);
    sprintf (buffer, "%g", a -> dynamic_relaxation);
    SetDlgItemText (hControld, 102, buffer);

    sprintf (buffer, "%g", a -> static_tolerance);
    SetDlgItemText (hControld, 200, buffer);
    sprintf (buffer, "%g", a -> outer_tolerance);
    SetDlgItemText (hControld, 201, buffer);
    sprintf (buffer, "%g", a -> dynamic_tolerance);
    SetDlgItemText (hControld, 202, buffer);

    sprintf (buffer, "%d", a -> static_it);
    SetDlgItemText (hControld, 300, buffer);
    sprintf (buffer, "%d", a -> outer_it);
    SetDlgItemText (hControld, 301, buffer);
    sprintf (buffer, "%d", a -> dynamic_it);
    SetDlgItemText (hControld, 302, buffer);

    sprintf (buffer, "%g", a -> dt);
    SetDlgItemText (hControld, 400, buffer);
    sprintf (buffer, "%g", a -> duration);
    SetDlgItemText (hControld, 401, buffer);
    sprintf (buffer, "%g", a -> ramp_time);
    SetDlgItemText (hControld, 402, buffer);

	/*
	 * de-sensitize the controls
	 */

    SetSensitivity (hControld, FALSE);

	/*
	 * store the baseline numbers for easy restoration
	 */

    base_static_relaxation  = a -> static_relaxation;
    base_outer_relaxation   = a -> outer_relaxation;
    base_dynamic_relaxation = a -> dynamic_relaxation;
    base_static_tolerance   = a -> static_tolerance;
    base_dynamic_tolerance  = a -> dynamic_tolerance;
    base_outer_tolerance    = a -> outer_tolerance;
    base_static_iterations  = a -> static_it;
    base_dynamic_iterations = a -> dynamic_it;
    base_outer_iterations   = a -> outer_it;

    base_time_step          = a -> dt;
    base_ramp_time          = a -> ramp_time;
    base_duration           = a -> duration;

    ControlProcessEvents (a -> solution);

    return 0;
}

void ControlLabels (n1, n2, n3, n4)
   char	 *n1;
   char  *n2;
   char  *n3;
   char	 *n4;
{
   SetDlgItemText(cd, 550, n1);
   SetDlgItemText(cd, 551, n2);
   SetDlgItemText(cd, 552, n3);
   SetDlgItemText(cd, 553, n4);


   return;
}

void ControlInfo (tm, it, err)
   double	tm;
   int		it;
   double	err;
{
   char		buffer [32];
   static double prev_tm = -1.0;

   if (tm != prev_tm) {
      sprintf (buffer,"%g", tm);
      SetDlgItemText (cd, 500, buffer);
      prev_tm = tm;
   }

   sprintf (buffer,"%d", it);
   SetDlgItemText (cd, 501, buffer);
   sprintf (buffer,"%g", err);
   SetDlgItemText (cd, 502, buffer);
}

void ControlAuxError(err)
   double err;
{
   char   buffer [32];

   sprintf (buffer,"%g", err);
   SetDlgItemText (cd, 503, buffer);
}

static char *msg_buffer [100];
static int   msg_line;

void ControlMessage (msg)
   char	*msg;
{
   static int           i = 0;
   char                 buffer [64];

   if (msg_buffer [i])
      free(msg_buffer [i]);

   sprintf(buffer,"%02d: %s", i+1, msg);
   msg_buffer [i] = strdup(buffer);
   SetDlgItemText (cd, 504, buffer);

   msg_line = i;

   i ++;
   if (i == 100)
      i = 0;

   return;
}

static void ScrollMessage(dir)
   int   dir;
{
   int   i;

   i = msg_line + dir;

   if (i < 0)
      i = 99;
   else if (i > 99)
      i = 0;

   SetDlgItemText (cd, 504, msg_buffer [i]);

   msg_line = i;

   return;
}

HINSTANCE hInst;

static int pause_flag = 0;

static BOOL WINAPI MainDlgProc(HWND hDlg, UINT msg, 
			       WPARAM wParam, LPARAM lParam)
{
   MSG	 message;
   int   i;
   HICON up_icon;
   HICON down_icon;

   switch( msg ) {

      case WM_INITDIALOG:                 
         up_icon = LoadIcon(hInst, MAKEINTRESOURCE( 1001 ));
         down_icon = LoadIcon(hInst, MAKEINTRESOURCE ( 1002 )); 

         SendDlgItemMessage(hDlg, 510, BM_SETIMAGE, IMAGE_ICON, (LPARAM) up_icon);
         SendDlgItemMessage(hDlg, 511, BM_SETIMAGE, IMAGE_ICON, (LPARAM) down_icon);

         SetDlgItemText(hDlg, 100, "");
         SetDlgItemText(hDlg, 101, "");
         SetDlgItemText(hDlg, 102, "");

         SetDlgItemText(hDlg, 200, "");
         SetDlgItemText(hDlg, 201, "");
         SetDlgItemText(hDlg, 202, "");

         SetDlgItemText(hDlg, 300, "");
         SetDlgItemText(hDlg, 301, "");
         SetDlgItemText(hDlg, 302, "");

         SetDlgItemText(hDlg, 400, "");
         SetDlgItemText(hDlg, 401, "");
         SetDlgItemText(hDlg, 402, "");

         SetDlgItemText(hDlg, 500, "");
         SetDlgItemText(hDlg, 501, "");
         SetDlgItemText(hDlg, 502, "");
         SetDlgItemText(hDlg, 503, "");
         SetDlgItemText(hDlg, 504, "");

         SetDlgItemText(hDlg, 550, "step");
         SetDlgItemText(hDlg, 551, "iterations");
         SetDlgItemText(hDlg, 552, "error");
         SetDlgItemText(hDlg, 553, "aux error");

         SetSensitivity(hDlg, FALSE);

         for (i = 0 ; i < 32 ; i++)
            msg_buffer [i] = NULL;

         break;

      case WM_COMMAND:
         if (wParam == IDCANCEL)   {
            EndDialog( hDlg, TRUE );

            if (unlink_input)
               unlink(in_name);

            if (!static_finished && out_name)
               unlink(out_name);
 
            exit (0);
         }
         else if (wParam == 510) { /* scroll back in message buffer */
            ScrollMessage(-1);
         }
         else if (wParam == 511) { /* scroll forward in message buffer */
            ScrollMessage(1);
         }
         else if (wParam == 5001) { /* pause */
            pause_flag = !pause_flag;

            if (!pause_flag) 
               SetSensitivity (hDlg, FALSE);
            else if (pause_flag) {
               SetSensitivity (hDlg, TRUE);

               while (1) {
                  if (PeekMessage (&message, NULL, 0, 0, PM_REMOVE)) {
                     if (message.message == WM_QUIT)
                        break;
      
                     if (!IsDialogMessage(hDlg, &message)) {
                        TranslateMessage (&message);
                        DispatchMessage (&message);
                     }
                  }
                  else if (!pause_flag) {
                     SetSensitivity (hDlg, FALSE);
                     break;
                  }
               }
            }
         }
         else if (wParam == 5002) { /* update */
            Change(hDlg);
         }
         else if (wParam == 5003) { /* restore */
            Restore(hDlg);
         }

         break;      
   }

   return FALSE;
} 

int APIENTRY WinMain(HINSTANCE hInstance, HINSTANCE hPrev, 
		     LPSTR lpCmd, int nShow)
{
   int    i;
   HWND   hDlg;
   char **argv, **argvlist, *p;
   int    argc, size;
   char   buffer[MAX_PATH];
   WNDCLASS wc;

	/*
	 * put the command line into argc, argv format
	 */

   for (size = 5, p = lpCmd; *p != '\0'; p++) {
      if (isspace(*p)) {
         size++;
	 while (isspace(*p)) 
	    p++;
         

         if (*p == '\0') 
	    break;
      }
   }

   argvlist = (char **) malloc((unsigned) (size * sizeof(char *)));
   argv = argvlist;

    /*
     * Parse the Windows command line string.  If an argument begins with a
     * double quote, then spaces are considered part of the argument until the
     * next double quote.  The argument terminates at the second quote.  Note
     * that this is different from the usual Unix semantics.
     *
     * We actually start to look only after we find a - ... this really
     * only works for progs like cable with clearly - delineated args.
     * The problem is that the shell command from VB seems to stick
     * parts of the the app name (argv [0]) into lpCmd if the app name
     * has spaces in it, even when there are quotes. This seems to work
     * ok from a DOS prompt, but not from Shell(), go figure.
     */


   for (i = 1, p = strchr(lpCmd, '-') ; *p != '\0'; i++) {
      while (isspace(*p)) 
         p++;
	
      if (*p == '\0') 
         break;
	
      if (*p == '"') {
	 p++;
	 argv[i] = p;
	 while ((*p != '\0') && (*p != '"')) 
	    p++;
	 
      } 
      else {
         argv[i] = p;
	 while (*p != '\0' && !isspace(*p)) 
            p++;
      }

      if (*p != '\0') {
         *p = '\0';
         p++;
      }
   }

   GetModuleFileName(NULL, buffer, sizeof(buffer));
   argv [0]= buffer;

   argv[i] = NULL;
   argc = i;

           
   GetClassInfo(NULL, "#32770", &wc);
   wc.hIcon = LoadIcon(hInstance, MAKEINTRESOURCE( 1000 ));
   wc.hInstance = hInstance;
 
   RegisterClass(&wc);
                                        
   hDlg = CreateDialog(hInstance, MAKEINTRESOURCE( 10000 ), NULL, MainDlgProc);

   cd = hDlg;
   hInst = hInstance;

   ControlProcessEvents (NULL);

	/*
	 * transfer control to the real main()
	 */

   CableMain(argc, argv);

   return FALSE;
}

static BOOL WINAPI WinErrorDlgProc(HWND hDlg, UINT msg, 
		                  WPARAM wParam, LPARAM lParam)
{
   switch( msg )
   {
      case WM_INITDIALOG:                
         SetDlgItemText( hDlg, 100, (LPSTR)lParam );   //  Show any passed text
         return( TRUE );
         break;

      case WM_COMMAND:
         if( wParam == IDOK )
         {
            EndDialog( hDlg, ( wParam == IDOK ) );        //  Return TRUE on OK
            return( TRUE );
         }
         break;
   }

   return( FALSE );
}

void WinErrorDialog (buffer)
   char	*buffer;
{
   DialogBoxParam( hInst, MAKEINTRESOURCE( 10001 ), cd,
                  WinErrorDlgProc, (LPARAM) buffer ); 

   return;
}

static void TextSize(file_name, lines, cols, size)
   char	*file_name;
   int  *lines;
   int  *cols;
   int  *size;
{
    int   length;
    int   num_lines;
    int   num_columns;
    int   num_bytes;
    FILE *fp;
    int   n;

    if ((fp = fopen (file_name, "r")) == NULL)
	return;
   
    num_lines = 0;
    num_columns = 0;
    num_bytes = 0;
    
    while ((n = fscanf (fp, " %*[^\n]%n%*[\n]", &length)) != EOF) {
	if (length > num_columns)
	    num_columns = length;
       
        num_bytes += length;
 
	num_lines ++;
        if (feof(fp))
           break;
    }
   
    fclose (fp);
   
    *lines = num_lines;
    *cols = num_columns;
    *size = num_bytes; 
}

static BOOL WINAPI WinErrorFileDlgProc(HWND hDlg, UINT msg, 
		                       WPARAM wParam, LPARAM lParam)
{
   int          sz;
   int		lines;
   int		cols;
   FILE	       *fp;
   char        *buffer;
   char		line [80];
   char	        temp [80];
   int		i;
   int		pos;
   

   switch( msg )
   {
      case WM_INITDIALOG:                
         TextSize((char *) lParam, &lines, &cols, &sz);
        
         fp = fopen((char *) lParam, "r");

         buffer = (char *) malloc(sizeof(buffer) * (sz + lines*2));
       
         pos = 0;
         for (i = 0 ; i < lines ; i++) {
            fgets(line, 80, fp);

            line [strlen(line) - 1] = 0;

            sprintf(&(buffer [pos]), "%s \r\n", line);

            pos = strlen(buffer);
         }

         fclose (fp);

         SetWindowText(GetDlgItem(hDlg, 100), buffer);
         SendMessage (GetDlgItem(hDlg, 100), EM_SETSEL, (WPARAM) -1, (LPARAM) 0);
         HideCaret(GetDlgItem(hDlg, 100));
         

         return( TRUE );
        
         break;

      case WM_COMMAND:
         if( wParam == IDOK )
         {
            EndDialog( hDlg, ( wParam == IDOK ) );        //  Return TRUE on OK
            return( TRUE );
         }
         break;
   }

   return( FALSE );
}

void WinErrorFileDialog (fname)
   char *fname;
{
   DialogBoxParam( hInst, MAKEINTRESOURCE( 10002 ), cd,
                  WinErrorFileDlgProc, (LPARAM) fname ); 


   return;
}

void ControlMainLoop ( )
{
   MSG  message;
/*
   SetFocus (GetDlgItem(cd, IDCANCEL));
*/
   while (1) {
      if (PeekMessage (&message, NULL, 0, 0, PM_REMOVE)) {
         if (message.message == WM_QUIT)
            break;
      
         if (!IsDialogMessage(cd, &message)) {
            TranslateMessage (&message);
            DispatchMessage (&message);
         }
      }
   }

   return;
}

void ControlPlotResults()
{
}

void ControlTabulateResults()
{
}


void
ControlPlotSnaps(Problem *p, Environment *e)
{
}

void
ControlPlotTime(double t, Problem *p, Environment *e)
{
}

# endif /* if defined WINGUI */
