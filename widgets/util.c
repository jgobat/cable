/*
    This file is part of the FElt finite element analysis package.
    Copyright (C) 1993-1997 Jason I. Gobat and Darren C. Atkinson

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

/************************************************************************
 * File:	util.c							*
 *									*
 * Description: This file contains miscellaneous utility functions.	*
 ************************************************************************/

# include <stdio.h>
# include <X11/Intrinsic.h>
# include <X11/StringDefs.h>
# include <X11/Shell.h>
# include <X11/Xaw/AsciiText.h>
# include <X11/Xaw/List.h>
# include <X11/Xaw/MenuButton.h>
# include <X11/Xaw/Viewport.h>
# include "Layout.h"
# include "scroll.h"
# include "util.h"
# include "post.h"


# define min(x,y) ((int) (x) < (int) (y) ? (x) : (y))
# define max(x,y) ((int) (x) > (int) (y) ? (x) : (y))

# if NeedWidePrototypes
# define BOOLEAN int
# define DIMENSION unsigned
# else
# define BOOLEAN Boolean
# define DIMENSION Dimension
# endif

	/*
	 * allow for the _possibility_ of compiling with X11R4
	 */

# ifdef XlibSpecificationRelease
# if XlibSpecificationRelease < 5
typedef char *XPointer;
# endif
# else
typedef char *XPointer;
# endif

static String text_table =
"<Key>:       no-op()\n\
 <BtnDown>:   no-op()\n\
 <BtnUp>:     no-op()\n\
 <BtnMotion>: no-op()";

static XtTranslations text_translations = NULL;


/************************************************************************
 * Function:	StringToLayout						*
 *									*
 * Description:	Converts a value of type XtRString to a value of type	*
 *		XtRLayout.						*
 ************************************************************************/

XtArgVal StringToLayout (widget, string)
    Widget widget;
    String string;
{
    XrmValue from;
    XrmValue to;
    XtArgVal layout;


    from.addr = (XPointer) string;
    from.size = strlen (string) + 1;

    to.addr = (XPointer) &layout;
    to.size = sizeof (layout);

    XtInitializeWidgetClass (layoutWidgetClass);

    if (XtConvertAndStore (widget, XtRString, &from, XtRLayout, &to) == False)
	return (XtArgVal) NULL;

    return layout;
}


/************************************************************************
 * Function:	CenterOnWidget						*
 *									*
 * Description: Sets the x and y resources so that specified shell	*
 *		widget will be centered over the specified widget.  The	*
 *		shell will be centered as closely as possible so that	*
 *		the entire shell is located on the screen.		*
 ************************************************************************/

void CenterOnWidget (shell, center, force)
    Widget  shell;
    Widget  center;
    BOOLEAN force;
{
    Arg       args [4];
    String    geometry;
    Position  x;
    Position  y;
    Position  center_x;
    Position  center_y;
    Dimension center_width;
    Dimension center_height;
    Dimension width;
    Dimension height;
    Dimension bw;


    if (force == False) {
	XtSetArg (args [0], XtNgeometry, &geometry);
	XtGetValues (shell, args, 1);
	if (geometry != NULL)
	    return;
    }


    XtSetArg (args [0], XtNx,      &center_x);
    XtSetArg (args [1], XtNy,      &center_y);
    XtSetArg (args [2], XtNwidth,  &center_width);
    XtSetArg (args [3], XtNheight, &center_height);
    XtGetValues (center, args, 4);

    XtSetArg (args [0], XtNwidth,       &width);
    XtSetArg (args [1], XtNheight,      &height);
    XtSetArg (args [2], XtNborderWidth, &bw);
    XtGetValues (shell, args, 3);

    x = center_x + (center_width - width) / (unsigned) 2;
    y = center_y + (center_height - height) / (unsigned) 2;

    x = max (0, min (x, WidthOfScreen (XtScreen (shell)) - width - 2 * bw));
    y = max (0, min (y, HeightOfScreen (XtScreen (shell)) - height - 2 * bw));

    XtSetArg (args [0], XtNx, x);
    XtSetArg (args [1], XtNy, y);
    XtSetValues (shell, args, 2);
}


/************************************************************************
 * Function:	CenterOnScreen						*
 *									*
 * Description: Sets the x and y resources so that specified shell	*
 *		widget will be centered on the screen.			*
 ************************************************************************/

void CenterOnScreen (shell, force)
    Widget  shell;
    BOOLEAN force;
{
    Arg       args [3];
    String    geometry;
    Position  x;
    Position  y;
    Dimension width;
    Dimension height;
    Dimension bw;


    if (force == False) {
	XtSetArg (args [0], XtNgeometry, &geometry);
	XtGetValues (shell, args, 1);
	if (geometry != NULL)
	    return;
    }


    XtSetArg (args [0], XtNwidth,       &width);
    XtSetArg (args [1], XtNheight,      &height);
    XtSetArg (args [2], XtNborderWidth, &bw);
    XtGetValues (shell, args, 3);

    x = (WidthOfScreen (XtScreen (shell)) - width - 2 * bw) / (unsigned) 2;
    y = (HeightOfScreen (XtScreen (shell)) - height - 2 * bw) / (unsigned) 2;

    XtSetArg (args [0], XtNx, x);
    XtSetArg (args [1], XtNy, y);
    XtSetValues (shell, args, 2);
}


/************************************************************************
 * Function:	AddDeleteWindowProtocol					*
 *									*
 * Description:	Adds the WM_DELETE_WINDOW atom to the display and adds	*
 *		the WM_PROTOCOL action to the specified widget.		*
 ************************************************************************/

void AddDeleteWindowProtocol (shell, action)
    Widget shell;
    String action;
{
    static Atom    atom = None;
    XtTranslations translations;
    char           table [256];


    if (atom == None)
	atom = XInternAtom (XtDisplay (shell), "WM_DELETE_WINDOW", False);

    XSetWMProtocols (XtDisplay (shell), XtWindow (shell), &atom, 1);

    sprintf (table, "<ClientMessage>WM_PROTOCOLS: %s", action);
    translations = XtParseTranslationTable (table);
    XtOverrideTranslations (shell, translations);
}


/************************************************************************
 * Function:	WarpToCenter						*
 *									*
 * Description:	Warps the pointer to the center of the specified	*
 *		widget.							*
 ************************************************************************/

void WarpToCenter (w)
    Widget w;
{
    Arg       args [2];
    Dimension width;
    Dimension height;


    XtSetArg (args [0], XtNwidth,  &width);
    XtSetArg (args [1], XtNheight, &height);
    XtGetValues (w, args, XtNumber (args));

    width /= 2;
    height /= 2;

    XWarpPointer (XtDisplay (w), None, XtWindow (w), 0, 0, 0, 0, width, height);
}


/************************************************************************
 * Function:	ListCursorMove						*
 *									*
 * Description:	An action procedure which moves to the previous or next	*
 *		entry in the list.  The list will scroll as needed.  If	*
 *		no item is currently selected then an initial selection	*
 *		of the first or last visible item will be made.		*
 ************************************************************************/

static void ListCursorMove (w, event, params, num_params)
    Widget    w;
    XEvent   *event;
    String   *params;
    Cardinal *num_params;
{
    Arg			 args [6];
    int			 top_index;
    int			 bottom_index;
    int			 new_index;
    Position		 y;
    Dimension		 height;
    Dimension		 internal_height;
    Dimension		 row_spacing;
    Dimension		 entry_height;
    Widget		 list;
    WidgetList		 children;
    Cardinal		 num_children;
    Cardinal		 i;
    int			 number_strings;
    String		*entries;
    XFontStruct		*font;
    XawListReturnStruct	*info;
    XawListReturnStruct	 call_data;


    /* Locate the list widget. */

    XtSetArg (args [0], XtNchildren,	&children);
    XtSetArg (args [1], XtNnumChildren,	&num_children);
    XtGetValues (w, args, 2);

    for (i = 0; i < num_children; i ++)
	if (XtClass (children [i]) == listWidgetClass)
	    break;

    list = children [i];
    info = XawListShowCurrent (list);


    /* Retrieve the properties of the viewport and list widgets. */

    XtSetArg (args [0], XtNy,              &y);
    XtSetArg (args [1], XtNinternalHeight, &internal_height);
    XtSetArg (args [2], XtNrowSpacing,     &row_spacing);
    XtSetArg (args [3], XtNnumberStrings,  &number_strings);
    XtSetArg (args [4], XtNfont,           &font);
    XtSetArg (args [5], XtNlist,           &entries);
    XtGetValues (list, args, XtNumber (args));

    XtSetArg (args [0], XtNheight, &height);
    XtGetValues (w, args, 1);

    if (number_strings == 0)
	return;


    /* Compute the height in pixels of a list entry. */

    entry_height  = font -> max_bounds.ascent + font -> max_bounds.descent;
    entry_height += row_spacing;


    /* Compute the top and bottom visible list indicies. */

    new_index = -1;
    top_index = (-y + internal_height) / (unsigned) entry_height;
    bottom_index = (-y + internal_height + height) / (unsigned)entry_height - 1;
    if (bottom_index >= number_strings)
	bottom_index = number_strings - 1;


    /* Move to the next entry. */

    if (!strcmp (params [0], "down")) {
	if (info -> list_index == XAW_LIST_NONE) {
	    new_index = top_index;
	    XawListHighlight (list, top_index);

	} else if (info -> list_index < number_strings - 1) {
	    new_index = info -> list_index + 1;
	    XawListHighlight (list, new_index);
	    if (info -> list_index == bottom_index)
		XawViewportSetCoordinates (w, 0, -y + entry_height);
	}


    /* Move to the previous entry. */

    } else {
	if (info -> list_index == XAW_LIST_NONE) {
	    new_index = bottom_index;
	    XawListHighlight (list, bottom_index);

	} else if (info -> list_index > 0) {
	    new_index = info -> list_index - 1;
	    XawListHighlight (list, new_index);
	    if (info -> list_index == top_index)
		XawViewportSetCoordinates (w, 0, -y - entry_height);
	}
    }


    /* If the selected entry has changed then call the callbacks. */

    if (new_index != -1) {
	call_data.string = entries [new_index];
	call_data.list_index = new_index;
	XtCallCallbacks (list, XtNcallback, &call_data);
    }
}


/************************************************************************
 * Function:	ListAddCursorTranslations				*
 *									*
 * Description: Add translations to a "scrolled list" (viewport with	*
 *		list child) to allow the cursor keys to be used to	*
 *		change the selected item.				*
 ************************************************************************/

void ListAddCursorTranslations (viewport)
    Widget viewport;
{
    ListAddCursorAccelerators (viewport, NULL);
}


/************************************************************************
 * Function:	ListAddCursorAccelerators				*
 *									*
 * Description: Adds accelerations or translations from the "scrolled	*
 *		list" widget to the specified widget to allow the	*
 *		cursor keys to be used to change the selected item.	*
 ************************************************************************/

void ListAddCursorAccelerators (viewport, w)
    Widget viewport;
    Widget w;
{
    Arg	   args		  [1];
    static XtAccelerators accelerators;
    static XtTranslations translations;
    static XtAppContext	  app_context = NULL;
    static XtActionsRec	  actions [ ] = {{"ListCursorMove", ListCursorMove}};
    static String	  table = "#override \
			  <Key>Down: ListCursorMove(down)\n\
			  <Key>Up:   ListCursorMove(up)";


    if (app_context == NULL) {
	app_context = XtWidgetToApplicationContext (viewport);
	XtAppAddActions (app_context, actions, XtNumber (actions));
	translations = XtParseTranslationTable (table);
	accelerators = XtParseAcceleratorTable (table);
    }

    if (w == NULL)
	XtOverrideTranslations (viewport, translations);
    else {
	XtSetArg (args [0], XtNaccelerators, accelerators);
	XtSetValues (viewport, args, 1);
	XtInstallAccelerators (w, viewport);
    }
}


/************************************************************************
 * Function:	SetTextString						*
 *									*
 * Description:	Updates the value of a text widget positioning the	*
 *		cursor at the end of the string.			*
 ************************************************************************/

void SetTextString (w, value)
    Widget w;
    String value;
{
    Arg args [1];


    if (value == NULL)
	value = "";

    XawTextDisableRedisplay (w);
    XtSetArg (args [0], XtNstring, value);
    XtSetValues (w, args, 1);
    XawTextSetInsertionPoint (w, strlen (value));
    ScrollToInsertionPoint (w);
    XawTextEnableRedisplay (w);
}


/************************************************************************
 * Function:	GetTextString						*
 *									*
 * Description:	Retrieves the value of a text widget.			*
 ************************************************************************/

String GetTextString (w)
    Widget w;
{
    Arg    args [1];
    String value;


    XtSetArg (args [0], XtNstring, &value);
    XtGetValues (w, args, 1);
    return value;
}


/************************************************************************
 * Function:	GetTextWidth						*
 *									*
 * Description:	Returns the width of a text string.			*
 ************************************************************************/

Cardinal GetTextWidth (font, text, length)
    XFontStruct *font;
    String	 text;
    Cardinal	 length;
{
    int		dr;
    int		far;
    int		fdr;
    XCharStruct	info;


    XTextExtents (font, text, length, &dr, &far, &fdr, &info);
    return info.width;
}


/************************************************************************
 * Function:	SetLabelString						*
 *									*
 * Description:	Sets the label string of a label widget.		*
 ************************************************************************/

void SetLabelString (w, value)
    Widget w;
    String value;
{
    Arg args [1];


    if (value == NULL)
	value = "";

    XtSetArg (args [0], XtNlabel, value);
    XtSetValues (w, args, 1);
}


/************************************************************************
 * Function:	GetLabelString						*
 *									*
 * Description:	Retrieves the label string of a label widget.		*
 ************************************************************************/

String GetLabelString (w)
    Widget w;
{
    Arg    args [1];
    String value;


    XtSetArg (args [0], XtNlabel, &value);
    XtGetValues (w, args, 1);
    return value;
}


/************************************************************************
 * Function:	AutoRepeat						*
 *									*
 * Description:	Sets the keyboard auto repeat mode to the specified	*
 *		value.							*
 ************************************************************************/

static void AutoRepeat (w, event, params, num_params)
    Widget    w;
    XEvent   *event;
    String   *params;
    Cardinal *num_params;
{
# if 0
    static XKeyboardState state;


    if (*num_params == 0 || !strcmp (params [0], "saved")) {
	if (state.global_auto_repeat)
	    XAutoRepeatOn (XtDisplay (w));
	else
	    XAutoRepeatOff (XtDisplay (w));

    } else if (!strcmp (params [0], "on") || !strcmp (params [0], "true")) {
	XGetKeyboardControl (XtDisplay (w), &state);
	XAutoRepeatOn (XtDisplay (w));

    } else if (!strcmp (params [0], "off") || !strcmp (params [0], "false")) {
	XGetKeyboardControl (XtDisplay (w), &state);
	XAutoRepeatOff (XtDisplay (w));
    }

    XSync (XtDisplay (w), 0);	 /* very important */
# endif
}


/************************************************************************
 * Function:	AddAutoRepeatAction					*
 *									*
 * Description:	Adds the AutoRepeat action to the specified application	*
 *		context.						*
 ************************************************************************/

void AddAutoRepeatAction (app_context)
    XtAppContext app_context;
{
    static XtAppContext	context = NULL;
    static XtActionsRec	actions [ ] = {{"AutoRepeat", AutoRepeat}};


    if (context == NULL) {
	context = app_context;
	XtAppAddActions (app_context, actions, XtNumber (actions));
    }
}


/************************************************************************
 * Function:	CreateHelpButton					*
 *									*
 * Description:	Creates a help button which will popup a shell with	*
 *		formatted text when activated by a key or button press.	*
 ************************************************************************/

Widget CreateHelpButton (parent, name)
    Widget parent;
    String name;
{
    Arg    args [3];
    Widget button;
    Widget shell;
    Widget text;


    if (text_translations == NULL)
	text_translations = XtParseTranslationTable (text_table);


    XtSetArg (args [0], XtNmenuName, "shell");

    button = XtCreateManagedWidget (name, menuButtonWidgetClass,
		parent, args, 1);


    XtSetArg (args [0], XtNallowShellResize, True);

    shell = XtCreatePopupShell ("shell", overrideShellWidgetClass,
		button, args, 1);

    AddPostMenuActions (shell);


    XtSetArg (args [0], XtNresize,       XawtextResizeHeight);
    XtSetArg (args [1], XtNwrap,         XawtextWrapWord);
    XtSetArg (args [2], XtNdisplayCaret, False);

    text = XtCreateManagedWidget ("text", asciiTextWidgetClass,
		shell, args, 3);

    XtOverrideTranslations (text, text_translations);

    return button;
}


/************************************************************************
 * Function:	UpdateHelpMessage					*
 *									*
 * Description:	Updates the text widget within the help shell.		*
 ************************************************************************/

void UpdateHelpMessage (button, message, width)
    Widget    button;
    String    message;
    DIMENSION width;
{
    Widget   text;
    Cardinal i;
    Arg      args [1];


    text = XtNameToWidget (button, "shell.text");

    if (width != 0) {
	XtSetArg (args [0], XtNwidth, width);
	XtSetValues (text, args, 1);
    }

    if (message != NULL) {
	SetTextString (text, message);
	for (i = 0; i < 30; i ++)
	    XtCallActionProc (text, "end-of-file", NULL, NULL, 0);
    }
}
