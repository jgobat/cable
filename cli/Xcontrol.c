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

/****************************************************************************
 * 
 * File:	Xcontrol.c
 *
 * Description:	Contains routines to offer graphical solution controls
 *
 ****************************************************************************/


# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "compress.h"
# include "allocate.h"
# include "error.h"
# include "problem.h"
# include "solve.h"
# include "output.h"

# ifdef GUI

# include <X11/Xos.h>
# include <X11/Intrinsic.h>
# include <X11/StringDefs.h>
# include <X11/Shell.h>
# include <X11/Xaw/AsciiText.h>
# include <X11/Xaw/Toggle.h>
# include <X11/Xaw/Command.h>
# include <X11/Xaw/Label.h>
# include <X11/Xaw/Repeater.h>
# include "Layout.h"
# include "TabGroup.h"
# include "util.h"

static Analysis *aptr;

static char *defaults [ ] = {
#include "cable_ad.h"
NULL
};

extern int	 unlink_input;
extern char	*in_name;
extern int	 static_finished;
extern char	*out_name;

extern XtArgVal Float2Arg ( );

Widget 		toplevel;
XtAppContext 	app_context = NULL;

typedef struct control_dialog *ControlDialog;

struct control_dialog {
    Widget  shell;
    Widget  layout;			/*	Layout  layout		*/
    Widget  static_tolerance;	
    Widget  dynamic_tolerance;
    Widget  outer_tolerance;
    Widget  static_relaxation;
    Widget  dynamic_relaxation;
    Widget  outer_relaxation;
    Widget  static_iterations;
    Widget  dynamic_iterations;
    Widget  outer_iterations;
    Widget  time_step;
    Widget  duration;
    Widget  ramp_time;
    Widget  quit;
    Widget  pause;
    Widget  update;
    Widget  restore;
    Widget  step;
    Widget  iter;
    Widget  err1;
    Widget  err2;
    Widget  message;
    Widget  up;
    Widget  down;
    double  base_static_tolerance;	
    double  base_dynamic_tolerance;
    double  base_outer_tolerance;
    double  base_time_step;
    double  base_duration;
    double  base_ramp_time;
    double  base_static_relaxation;
    double  base_outer_relaxation;
    double  base_dynamic_relaxation;
    int     base_static_iterations;
    int     base_dynamic_iterations;
    int     base_outer_iterations;
    char   *buffer [100];
    int     buffer_line;
};

static String labels [ ] = {
    "tolerance","iterations","relaxation",
    "static", "dynamic","outer",
    "dynamic", "time-step", "duration", "ramp-time",
    "step", "iter", "error", "aux err"
};

static String label_names [ ] = {
    "tolerance_name","iterations_name","relaxation_name",
    "static_name", "dynamic_name","outer_name",
    "dynamic1_name", "time_step_name", "duration_name", "ramp_time_name",
    "step_name", "iter_name", "err1_name", "err2_name"
};


/* Resources */

static Pixel highlight;

static char layout_string [ ] =
"vertical { \
   8 \
   horizontal { \
       4 \
       width dynamic_name \
       2 \
       relaxation_name \
       8 \
       (width static_relaxation - width relaxation_name) \
       tolerance_name \
       8 \
       (width static_tolerance - width tolerance_name) \
       iterations_name \
       4 <+inf -100%> \
   } \
   2 \
   horizontal { \
       4 \
       ((width dynamic_name - width static_name) / 2) \
       static_name \
       ((width dynamic_name - width static_name) / 2) \
       2 \
       static_relaxation \
       8 \
       static_tolerance \
       8 \
       static_iterations \
       4 <+inf -100%> \
   } \
   4 \
   horizontal { \
       4 \
       ((width dynamic_name - width outer_name) / 2) \
       outer_name \
       ((width dynamic_name - width outer_name) / 2) \
       2 \
       outer_relaxation \
       8 \
       outer_tolerance \
       8 \
       outer_iterations \
       4 <+inf -100%> \
   } \
   4 \
   horizontal { \
       4 \
       dynamic_name \
       2 \
       dynamic_relaxation \
       8 \
       dynamic_tolerance \
       8 \
       dynamic_iterations \
       4 <+inf -100%> \
   } \
   8 \
   separator1 <+inf -100% *> \
   4 \
   horizontal { \
       4 \
       width dynamic1_name \
       2 \
       time_step_name \
       8 \
       (width time_step - width time_step_name) \
       duration_name \
       8 \
       (width duration - width duration_name) \
       ramp_time_name \
       4 <+inf -100%> \
   } \
   2 \
   horizontal { \
       4 \
       dynamic1_name \
       2 \
       time_step \
       8 \
       duration \
       8 \
       ramp_time \
       4 <+inf -100%> \
   } \
   8 \
   separator3 <+inf -100% *> \
   4 \
   horizontal { \
       4 \
       ((width step - width step_name) / 2) \
       step_name \
       ((width step - width step_name) / 2) \
       ((width dynamic_name + 3*width static_relaxation - 4*width step + 18) / 3)  \
       ((width iter - width iter_name) / 2) \
       iter_name \
       ((width iter - width iter_name) / 2) \
       ((width dynamic_name + 3*width static_relaxation - 4*width step + 18) / 3)  \
       ((width err1 - width err1_name) / 2) \
       err1_name \
       ((width err1 - width err1_name) / 2) \
       ((width dynamic_name + 3*width static_relaxation - 4*width step + 18) / 3)  \
       ((width err2 - width err2_name) / 2) \
       err2_name \
       ((width err2 - width err2_name) / 2) \
    } \
    2 \
    horizontal { \
       4 \
       step \
       ((width dynamic_name + 3*width static_relaxation - 4*width step + 18) / 3)  \
       iter \
       ((width dynamic_name + 3*width static_relaxation - 4*width step + 18) / 3)  \
       err1 \
       ((width dynamic_name + 3*width static_relaxation - 4*width step + 18) / 3)  \
       err2 \
    } \
    4 \
    horizontal { \
       4 \
       message <+inf -100% *> \
       4 \
       down \
       4 \
       up \
       4 \
    } \
    8 \
    separator2 <+inf -100% *> \
    8 \
    horizontal { \
       4 \
       quit \
       8 \
       pause \
       8 \
       update \
       8 \
       restore \
       8 <+inf -100%> \
    } \
    8 \
}";

/* Bitmaps */

#define up_width 12
#define up_height 12
static char up_bits[] = {
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x60, 0x00, 0xf0, 0x00, 0xf8, 0x01,
   0xfc, 0x03, 0xfe, 0x07, 0xfe, 0x07, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};

#define down_width 12
#define down_height 12
static char down_bits[] = {
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xfe, 0x07, 0xfe, 0x07, 0xfc, 0x03,
   0xf8, 0x01, 0xf0, 0x00, 0x60, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};

static Pixmap up_bitmap;
static Pixmap down_bitmap;

static Arg color_args [ ] = {
    {XtNborderColor, (XtArgVal) &highlight},
};

static Arg layout_args [ ] = {
    {XtNlayout, (XtArgVal) NULL},
};

static Arg label_args [ ] = {
    {XtNlabel,	            (XtArgVal) ""},
    {XtNborderWidth,        (XtArgVal) 0},
};

static Arg input_args [ ] = {
    {XtNeditType,    (XtArgVal) XawtextEdit},
    {XtNborderWidth, (XtArgVal) 0},
    {XtNpieceSize,   (XtArgVal) 32},
    {XtNcursorName,  (XtArgVal) "left_ptr"},
    {XtNwidth,       (XtArgVal) 100},
    {XtNdisplayCaret, (XtArgVal) True},
};

static Arg output_args [ ] = {
    {XtNeditType,    (XtArgVal) XawtextRead},
    {XtNborderWidth, (XtArgVal) 0},
    {XtNpieceSize,   (XtArgVal) 32},
    {XtNcursorName,  (XtArgVal) "left_ptr"},
    {XtNwidth,       (XtArgVal) 100},
    {XtNdisplayCaret, (XtArgVal) False},
};

static Arg info_args [ ] = {
    {XtNeditType,    (XtArgVal) XawtextRead},
    {XtNborderWidth, (XtArgVal) 0},
    {XtNpieceSize,   (XtArgVal) 32},
    {XtNcursorName,  (XtArgVal) "left_ptr"},
    {XtNwidth,       (XtArgVal) 100},
    {XtNdisplayCaret, (XtArgVal) False},
};

static Arg core_args [ ] = {
    {XtNwidth,  (XtArgVal) 3},
    {XtNheight, (XtArgVal) 3},
};

static Arg repeater_args [ ] = {
    {XtNbitmap, (XtArgVal) NULL},
};

/* Translation tables */

static String text_table =
"<Btn1Down>:  SetFocus() select-start()";

static XtTranslations text_translations;

static String command_table =
"<Key>space:   AutoRepeat(off) set()\n\
 <KeyUp>space: AutoRepeat(saved) notify() unset()";

static XtTranslations command_translations;

static String toggle_table =
"<Btn1Down>:  SetFocus()\n\
 <Key>space:  toggle() notify()";

static XtTranslations toggle_translations;

static void SetSensitivity (controld, state)
   ControlDialog	controld;
   Boolean		state;
{
    XtSetSensitive (controld -> dynamic_iterations, state);
    XtSetSensitive (controld -> static_iterations, state);
    XtSetSensitive (controld -> outer_iterations, state);
    XtSetSensitive (controld -> dynamic_tolerance, state);
    XtSetSensitive (controld -> static_tolerance, state);
    XtSetSensitive (controld -> outer_tolerance, state);
    XtSetSensitive (controld -> dynamic_relaxation, state);
    XtSetSensitive (controld -> outer_relaxation, state);
    XtSetSensitive (controld -> static_relaxation, state);

    XtSetSensitive (controld -> time_step, state);
    XtSetSensitive (controld -> ramp_time, state);

    XtSetSensitive (controld -> update, state);
    XtSetSensitive (controld -> restore, state);

    return;
}

static void Pause (w, client_data, call_data)
   Widget       w;
   XtPointer    client_data;
   XtPointer    call_data;
{
   Arg		  args [1];
   Boolean	  state;
   XEvent	  event;
   ControlDialog  controld;

   controld = (ControlDialog) client_data;

   XtSetArg (args [0], XtNstate, &state);
   XtGetValues (w, args, 1);

   if (!state) {
      SetSensitivity (controld, True);
      return;
   }

   SetSensitivity (controld, True);
   FlushOutput ( );

   while (1) {

      XtGetValues (w, args, 1);
      if (!state)
         break;

      while (XtAppPending (app_context)) {
         XtAppNextEvent (app_context, &event);
         XtDispatchEvent (&event);
      }
   }

   SetSensitivity (controld, False);
   return;
}

static void Quit (w, client_data, call_data)
   Widget       w;
   XtPointer    client_data;
   XtPointer    call_data;
{
   XtUnmapWidget (toplevel);
   XtDestroyApplicationContext (app_context);

   if (unlink_input)
      unlink(in_name);

   if (!static_finished && out_name)
      unlink(out_name);

   exit (0);
}

static void ControlAction (w, event, params, num_params)
   Widget	w;
   XEvent	*event;
   String	*params;
   Cardinal	*num_params;
{
   Quit (w, NULL, NULL);
}

static void Restore (w, client_data, call_data)
   Widget       w;
   XtPointer    client_data;
   XtPointer    call_data;
{
   ControlDialog  controld;
   char		  buffer [32];

   controld = (ControlDialog) client_data;

   aptr -> static_relaxation = controld -> base_static_relaxation;
   sprintf (buffer, "%g", aptr -> static_relaxation);
   SetTextString (controld -> static_relaxation, buffer);

   aptr -> dynamic_relaxation = controld -> base_dynamic_relaxation;
   sprintf (buffer, "%g", aptr -> dynamic_relaxation);
   SetTextString (controld -> dynamic_relaxation, buffer);

   aptr -> outer_relaxation = controld -> base_outer_relaxation;
   sprintf (buffer, "%g", aptr -> outer_relaxation);
   SetTextString (controld -> outer_relaxation, buffer);

   aptr -> static_tolerance = controld -> base_static_tolerance;
   sprintf (buffer, "%g", aptr -> static_tolerance);
   SetTextString (controld -> static_tolerance, buffer);

   aptr -> dynamic_tolerance = controld -> base_dynamic_tolerance;
   sprintf (buffer, "%g", aptr -> dynamic_tolerance);
   SetTextString (controld -> dynamic_tolerance, buffer);

   aptr -> outer_tolerance = controld -> base_outer_tolerance;
   sprintf (buffer, "%g", aptr -> outer_tolerance);
   SetTextString (controld -> outer_tolerance, buffer);

   aptr -> static_it = controld -> base_static_iterations;
   sprintf (buffer, "%d", aptr -> static_it);
   SetTextString (controld -> static_iterations, buffer);

   aptr -> dynamic_it = controld -> base_dynamic_iterations;
   sprintf (buffer, "%d", aptr -> dynamic_it);
   SetTextString (controld -> dynamic_iterations, buffer);

   aptr -> outer_it = controld -> base_outer_iterations;
   sprintf (buffer, "%d", aptr -> outer_it);
   SetTextString (controld -> outer_iterations, buffer);

   aptr -> duration = controld -> base_duration;
   sprintf (buffer, "%g", aptr -> duration);
   SetTextString (controld -> duration, buffer);

   aptr -> ramp_time = controld -> base_ramp_time;
   sprintf (buffer, "%g", aptr -> ramp_time);
   SetTextString (controld -> ramp_time, buffer);

   aptr -> dt = controld -> base_time_step;
   sprintf (buffer, "%g", aptr -> dt);
   SetTextString (controld -> time_step, buffer);

   return;
}

static void Scroll (w, client_data, call_data)
   Widget       w;
   XtPointer    client_data;
   XtPointer    call_data;
{
   int		  dir;
   int		  i;
   ControlDialog  controld = (ControlDialog) client_data;
  
   if (w == controld -> up)
      dir = 1;
   else
      dir = -1;

   i = controld -> buffer_line + dir;
   
   if (i < 0)
      i = 99;
   else if (i > 99)
      i = 0;

   SetTextString(controld -> message, controld -> buffer [i]);

   controld -> buffer_line = i;

   return;
}

static void Change (w, client_data, call_data)
   Widget       w;
   XtPointer    client_data;
   XtPointer    call_data;
{
   ControlDialog  controld;
   char		 *ptr;
   char		  buffer [32];
   double	  value;
   int		  i;

   controld = (ControlDialog) client_data;

   ptr = GetTextString (controld -> static_relaxation);
   value = atof (ptr);
   if (value > -1.0 && value <= 1.0)
      aptr -> static_relaxation = value;
   else {
      sprintf (buffer, "%g", aptr -> static_relaxation);
      SetTextString (controld -> static_relaxation, buffer);
   }

   ptr = GetTextString (controld -> dynamic_relaxation);
   value = atof (ptr);
   if (value > 0.0 && value <= 1.0)
      aptr -> dynamic_relaxation = value;
   else {
      sprintf (buffer, "%g", aptr -> dynamic_relaxation);
      SetTextString (controld -> dynamic_relaxation, buffer);
   }

   ptr = GetTextString (controld -> outer_relaxation);
   value = atof (ptr);
   if (value > 0.0)
      aptr -> outer_relaxation = value;
   else {
      sprintf (buffer, "%g", aptr -> outer_relaxation);
      SetTextString (controld -> outer_relaxation, buffer);
   }

   ptr = GetTextString (controld -> static_iterations);
   i = atoi(ptr);
   if (i)
      aptr -> static_it = i;
   else {
      sprintf (buffer, "%d", aptr -> static_it);
      SetTextString (controld -> static_iterations, buffer);
   }

   ptr = GetTextString (controld -> dynamic_iterations);
   i = atoi(ptr);
   if (i) 
      aptr -> dynamic_it = i;
   else {
      sprintf (buffer, "%d", aptr -> dynamic_it);
      SetTextString (controld -> dynamic_iterations, buffer);
   }

   ptr = GetTextString (controld -> outer_iterations);
   i = atoi(ptr);
   if (i > 0)
      aptr -> outer_it = i;
   else {
      sprintf (buffer, "%d", aptr -> outer_it);
      SetTextString (controld -> outer_iterations, buffer);
   }
      
   ptr = GetTextString (controld -> static_tolerance);
   value = atof(ptr);
   if (value)
      aptr ->  static_tolerance = value;
   else {
      sprintf (buffer, "%g", aptr -> static_tolerance);
      SetTextString (controld -> static_tolerance, buffer);
   }

   ptr = GetTextString (controld -> dynamic_tolerance);
   value = atof(ptr);
   if (value)
      aptr -> dynamic_tolerance = value;
   else {
      sprintf (buffer, "%g", aptr -> dynamic_tolerance);
      SetTextString (controld -> dynamic_tolerance, buffer);
   }

   ptr = GetTextString (controld -> outer_tolerance);
   value = atof(ptr);
   if (value)
      aptr -> outer_tolerance = value;
   else {
      sprintf (buffer, "%g", aptr -> outer_tolerance);
      SetTextString (controld -> outer_tolerance, buffer);
   }

   ptr = GetTextString (controld -> duration);
   value = atof(ptr);
   if (value)
      aptr -> duration = value;
   else {
      sprintf (buffer, "%g", aptr -> duration);
      SetTextString (controld -> duration, buffer);
   }

   ptr = GetTextString (controld -> ramp_time);
   value = atof(ptr);
   if (value)
      aptr -> ramp_time = value;
   else {
      sprintf (buffer, "%g", aptr -> ramp_time);
      SetTextString (controld -> ramp_time, buffer);
   }

   ptr = GetTextString (controld -> time_step);
   value = atof(ptr);
   if (value)
      aptr -> dt = value;
   else {
      sprintf (buffer, "%g", aptr -> dt);
      SetTextString (controld -> time_step, buffer);
   }

   return;
}

static ControlDialog cd;

int CreateControlDialog (argc, argv, quiet)
   int		*argc;
   char		*argv [];
   int		 quiet;
{
    Window		window;
    
    Cardinal		i;
    Widget		group [16];
    ControlDialog	controld;
    XEvent		event;
    static XtActionsRec actions [ ] = {{"ControlAction", ControlAction}};
# ifdef __CYGWIN32__
    Display	       *display;
# endif
    
    fprintf(stderr,"creating control dialog\n");
	/*
	 * we may need to do this here because XtAppInitialize isn't
	 * nice about how it exits and we want to make sure that
	 * our own Fatal gets called
	 */

# ifdef __CYGWIN32__
    display = XOpenDisplay(NULL);
    if (display == NULL)
       Fatal ("could not open connection to X server");
# endif

    toplevel = XtAppInitialize (&app_context, "Cable", NULL, 0,
                                argc, argv, defaults, NULL, 0);

    if (toplevel == NULL || app_context == NULL)
       Fatal ("could not initialize a shell for X");

    if (quiet)	/* we never actually need the control dialog */
       return 0;

    XtAppAddActions (app_context, actions, 1);
    AddAutoRepeatAction (app_context);

    text_translations = XtParseTranslationTable (text_table);
    command_translations = XtParseTranslationTable (command_table);
    toggle_translations = XtParseTranslationTable (toggle_table);

    controld = XtNew (struct control_dialog);

    controld -> shell = toplevel;

    layout_args [0].value = StringToLayout (toplevel, layout_string);
    controld -> layout   = XtCreateManagedWidget ("layout",
			 layoutWidgetClass, toplevel,
			 layout_args, XtNumber (layout_args));

    controld -> static_relaxation = XtCreateManagedWidget ("static_relaxation",
                         asciiTextWidgetClass, controld -> layout,
                         input_args, XtNumber(input_args));

    controld -> outer_relaxation = XtCreateManagedWidget ("outer_relaxation",
                         asciiTextWidgetClass, controld -> layout,
                         input_args, XtNumber(input_args));

    controld -> dynamic_relaxation =XtCreateManagedWidget ("dynamic_relaxation",
                         asciiTextWidgetClass, controld -> layout,
                         input_args, XtNumber(input_args));

    controld -> static_iterations = XtCreateManagedWidget ("static_iterations", 
                         asciiTextWidgetClass, controld -> layout,
                         input_args, XtNumber(input_args));

    controld -> static_tolerance = XtCreateManagedWidget ("static_tolerance", 
                         asciiTextWidgetClass, controld -> layout,
                         input_args, XtNumber(input_args));

    controld -> dynamic_iterations =XtCreateManagedWidget ("dynamic_iterations",
                         asciiTextWidgetClass, controld -> layout,
                         input_args, XtNumber(input_args));

    controld -> dynamic_tolerance = XtCreateManagedWidget ("dynamic_tolerance", 
                         asciiTextWidgetClass, controld -> layout,
                         input_args, XtNumber(input_args));

    controld -> outer_iterations =XtCreateManagedWidget ("outer_iterations",
                         asciiTextWidgetClass, controld -> layout,
                         input_args, XtNumber(input_args));

    controld -> outer_tolerance = XtCreateManagedWidget ("outer_tolerance", 
                         asciiTextWidgetClass, controld -> layout,
                         input_args, XtNumber(input_args));

    controld -> time_step = XtCreateManagedWidget ("time_step", 
                         asciiTextWidgetClass, controld -> layout,
                         input_args, XtNumber(input_args));

    controld -> ramp_time = XtCreateManagedWidget ("ramp_time", 
                         asciiTextWidgetClass, controld -> layout,
                         input_args, XtNumber(input_args));

    controld -> duration = XtCreateManagedWidget ("duration", 
                         asciiTextWidgetClass, controld -> layout,
                         input_args, XtNumber(input_args));

    controld -> step =XtCreateManagedWidget ("step",
                         asciiTextWidgetClass, controld -> layout,
                         info_args, XtNumber(info_args));

    controld -> iter = XtCreateManagedWidget ("iter", 
                         asciiTextWidgetClass, controld -> layout,
                         info_args, XtNumber(info_args));

    controld -> err1 = XtCreateManagedWidget ("err1", 
                         asciiTextWidgetClass, controld -> layout,
                         info_args, XtNumber(info_args));

    controld -> err2 = XtCreateManagedWidget ("err2", 
                         asciiTextWidgetClass, controld -> layout,
                         info_args, XtNumber(info_args));

    controld -> message = XtCreateManagedWidget ("message", 
                         asciiTextWidgetClass, controld -> layout,
                         output_args, XtNumber(output_args));

    window = RootWindowOfScreen (XtScreen (toplevel));

    up_bitmap = XCreateBitmapFromData (XtDisplay (toplevel), window,
                up_bits, up_width, up_height);

    down_bitmap = XCreateBitmapFromData (XtDisplay (toplevel), window,
                  down_bits, down_width, down_height);

    repeater_args [0].value = (XtArgVal) up_bitmap;

    controld -> up = XtCreateManagedWidget ("up",
                        repeaterWidgetClass, controld -> layout,
                        repeater_args, XtNumber (repeater_args));

    repeater_args [0].value = (XtArgVal) down_bitmap;

    controld -> down = XtCreateManagedWidget ("down",
                        repeaterWidgetClass, controld -> layout,
                        repeater_args, XtNumber (repeater_args));


    controld -> quit = XtCreateManagedWidget ("quit", 
                         commandWidgetClass, controld -> layout,
                         NULL, 0);

    controld -> update = XtCreateManagedWidget ("update", 
                         commandWidgetClass, controld -> layout,
                         NULL, 0);

    controld -> restore = XtCreateManagedWidget ("restore", 
                         commandWidgetClass, controld -> layout,
                         NULL, 0);

    controld -> pause = XtCreateManagedWidget ("pause", 
                         toggleWidgetClass, controld -> layout,
                         NULL, 0);

    for (i = 0 ; i < XtNumber (labels) ; i++) {
        label_args [0].value = (XtArgVal) labels [i];
        XtCreateManagedWidget (label_names [i], labelWidgetClass,
                   controld -> layout, label_args, XtNumber (label_args));
    }

    XtCreateManagedWidget ("separator1", coreWidgetClass,
                        controld -> layout, core_args, XtNumber (core_args));

    XtCreateManagedWidget ("separator2", coreWidgetClass,
                        controld -> layout, core_args, XtNumber (core_args));

    XtCreateManagedWidget ("separator3", coreWidgetClass,
                        controld -> layout, core_args, XtNumber (core_args));

    /* Create a tab group for the solution dialog. */

    i = 0;
    group [i++] = controld -> static_relaxation;
    group [i++] = controld -> outer_relaxation;
    group [i++] = controld -> dynamic_relaxation;
    group [i++] = controld -> static_tolerance;
    group [i++] = controld -> outer_tolerance;
    group [i++] = controld -> dynamic_tolerance;
    group [i++] = controld -> static_iterations;
    group [i++] = controld -> outer_iterations;
    group [i++] = controld -> dynamic_iterations;
    group [i++] = controld -> time_step;
    group [i++] = controld -> duration;
    group [i++] = controld -> ramp_time;
    group [i++] = controld -> quit;
    group [i++] = controld -> pause;
    group [i++] = controld -> update;
    group [i++] = controld -> restore;

    XtGetValues (controld -> layout, color_args, XtNumber (color_args));
    CreateTabGroup (controld -> shell, group, XtNumber (group), 
                    highlight, True);


    	/* 
	 * Add the necessary callbacks and translations
	 */

    XtAddCallback (controld -> up,   XtNcallback, Scroll, (XtPointer) controld);
    XtAddCallback (controld -> down, XtNcallback, Scroll, (XtPointer) controld);
    XtAddCallback (controld -> quit, XtNcallback, Quit,   NULL);
    XtAddCallback (controld -> pause,XtNcallback, Pause,  (XtPointer) controld);
    XtAddCallback (controld -> update, 
                    XtNcallback, Change, (XtPointer) controld);
    XtAddCallback (controld -> restore, 
                    XtNcallback, Restore, (XtPointer) controld);
    
    XtOverrideTranslations (controld -> static_relaxation, text_translations);
    XtOverrideTranslations (controld -> outer_relaxation, text_translations);
    XtOverrideTranslations (controld -> dynamic_relaxation, text_translations);
    XtOverrideTranslations (controld -> static_iterations, text_translations);
    XtOverrideTranslations (controld -> dynamic_iterations, text_translations);
    XtOverrideTranslations (controld -> outer_iterations, text_translations);
    XtOverrideTranslations (controld -> static_tolerance, text_translations);
    XtOverrideTranslations (controld -> dynamic_tolerance, text_translations);
    XtOverrideTranslations (controld -> outer_tolerance, text_translations);
    XtOverrideTranslations (controld -> time_step, text_translations);
    XtOverrideTranslations (controld -> duration, text_translations);
    XtOverrideTranslations (controld -> ramp_time, text_translations);
    XtOverrideTranslations (controld -> quit, command_translations);
    XtOverrideTranslations (controld -> update, command_translations);
    XtOverrideTranslations (controld -> restore, command_translations);
    XtOverrideTranslations (controld -> pause, toggle_translations);

	/*
	 * get it on the screen
	 */

    XtRealizeWidget (toplevel);
    AddDeleteWindowProtocol (toplevel, "ControlAction()");


    while (XtAppPending (app_context)) {
       XtAppNextEvent (app_context, &event);
       XtDispatchEvent (&event);
    }

    cd = controld;

    cd -> buffer_line = 0;

    for (i = 0 ; i < 100 ; i++)
       cd -> buffer [i] = NULL;

	/*
	 * set the duration insensitive
	 */

    XtSetSensitive(controld -> duration, False);

    return 0;
}

int 
ControlDialogInitialize(Analysis *a, void *controls)
{
    ControlDialog	controld;
    char		buffer [32];
    XEvent		event;

    aptr = a;

    controld = cd;

	/*
	 * initialize the controls
	 */

    sprintf (buffer, "%g", a -> static_relaxation);
    SetTextString (controld -> static_relaxation, buffer);
    sprintf (buffer, "%g", a -> dynamic_relaxation);
    SetTextString (controld -> dynamic_relaxation, buffer);
    sprintf (buffer, "%g", a -> outer_relaxation);
    SetTextString (controld -> outer_relaxation, buffer);

    sprintf (buffer, "%g", a -> static_tolerance);
    SetTextString (controld -> static_tolerance, buffer);
    sprintf (buffer, "%g", a -> outer_tolerance);
    SetTextString (controld -> outer_tolerance, buffer);
    sprintf (buffer, "%g", a -> dynamic_tolerance);
    SetTextString (controld -> dynamic_tolerance, buffer);

    sprintf (buffer, "%d", a -> static_it);
    SetTextString (controld -> static_iterations, buffer);
    sprintf (buffer, "%d", a -> outer_it);
    SetTextString (controld -> outer_iterations, buffer);
    sprintf (buffer, "%d", a -> dynamic_it);
    SetTextString (controld -> dynamic_iterations, buffer);

    sprintf (buffer, "%g", a -> dt);
    SetTextString (controld -> time_step, buffer);
    sprintf (buffer, "%g", a -> ramp_time);
    SetTextString (controld -> ramp_time, buffer);
    sprintf (buffer, "%.0f", a -> duration);
    SetTextString (controld -> duration, buffer);

	/*
	 * de-sensitize the controls
	 */

    SetSensitivity (controld, False);

	/*
	 * store the baseline numbers for easy restoration
	 */

    controld -> base_static_relaxation  = a -> static_relaxation;
    controld -> base_outer_relaxation   = a -> outer_relaxation;
    controld -> base_dynamic_relaxation = a -> dynamic_relaxation;
    controld -> base_static_tolerance   = a -> static_tolerance;
    controld -> base_dynamic_tolerance  = a -> dynamic_tolerance;
    controld -> base_outer_tolerance    = a -> outer_tolerance;
    controld -> base_static_iterations  = a -> static_it;
    controld -> base_dynamic_iterations = a -> dynamic_it;
    controld -> base_outer_iterations   = a -> outer_it;

    controld -> base_time_step          = a -> dt;
    controld -> base_ramp_time          = a -> ramp_time;
    controld -> base_duration           = a -> duration;

    SetFocus (controld -> pause);

    while (XtAppPending (app_context)) {
       XtAppNextEvent (app_context, &event);
       XtDispatchEvent (&event);
    }

    return 0;
}

/*
void
ControlPlotResults(Problem *p, Environment *e)
{
}
*/
void 
ControlLabels (char *n1, char *n2, char *n3, char *n4, void *data)
{
   Arg		args [1];
   Widget	w;

   w = XtNameToWidget (cd -> layout, "step_name");

   XtSetArg (args [0], XtNlabel, n1);
   XtSetValues (w, args, 1);

   return;
}

void 
ControlInfo (double tm, int it, double err, void *data)
{
   char		buffer [32];
   static double prev_tm = -1.0;

   if (tm != prev_tm) {
      sprintf (buffer,"%g", tm);
      SetTextString (cd -> step, buffer);
      prev_tm = tm;
   }

   sprintf (buffer,"%d", it);
   SetTextString (cd -> iter, buffer);
   sprintf (buffer,"%g", err);
   SetTextString (cd -> err1, buffer);
}

void 
ControlAuxError(double err, void *data)
{
   char   buffer [32];

   sprintf (buffer,"%g", err);
   SetTextString (cd -> err2, buffer);
}

void 
ControlMessage (char *msg, void *data)
{
   static int		i = 0;
   char			buffer [512];

   if (cd -> buffer [i])
      free(cd -> buffer [i]);

   sprintf(buffer,"%02d: %s", i+1, msg);
   cd -> buffer [i] = strdup(buffer);
   SetTextString (cd -> message, buffer);

   cd -> buffer_line = i;

   i ++;
   if (i == 100)
      i = 0;
  
   return;
}

void 
ControlProcessEvents (Solution *s)
{
   XEvent	event;

   if (app_context) {
      while (XtAppPending (app_context)) {
         XtAppNextEvent (app_context, &event);
         XtDispatchEvent (&event);
      }
   }
}

# endif /* if defined GUI */

void 
ControlPlotResults(Node *node, int nn)
{
}

void 
ControlTabulateResults(Node *node, int nn)
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
