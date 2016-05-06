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
 * File:	post.c							*
 *									*
 * Description:	This file contains the public and private function	*
 *		definitions for the postable menu mechanism.		*
 ************************************************************************/

# include <stdio.h>
# include <X11/keysym.h>
# include <X11/cursorfont.h>
# include <X11/Intrinsic.h>
# include <X11/StringDefs.h>
# include <X11/Xaw/MenuButton.h>
# include <X11/Xaw/SimpleMenu.h>


#define printf do_nothing

void do_nothing ( ) { }


static Position	x1;
static Position	y1;
static Position	x2;
static Position	y2;
static Time	when = 0;
static Boolean	posted = False;
static Boolean	from_inside = False;
static Boolean	from_button = False;
static Cursor	default_cursor;


static String new_table =
"<EnterWindow>: no-op()\n\
 <LeaveWindow>: no-op()";

static XtTranslations new_translations;


static String old_table =
"<EnterWindow>: highlight()\n\
 <LeaveWindow>: unhighlight()";

static XtTranslations old_translations;


static String post_table =
"<Key>:     PostMenu()\n\
 <BtnDown>: PostMenu()\n\
 <BtnUp>:   PostMenu()";

static XtTranslations post_translations;


/************************************************************************
 * Function:	highlight						*
 *									*
 * Description:	Highlights a menu entry given by event location.	*
 ************************************************************************/

static void highlight (w, event)
    Widget  w;
    XEvent *event;
{
    if (XtClass (w) == simpleMenuWidgetClass)
	XtCallActionProc (w, "highlight", event, NULL, 0);
}


/************************************************************************
 * Function:	notify							*
 *									*
 * Description:	Notifies a menu entry given by event location.		*
 ************************************************************************/

static void notify (w, event)
    Widget  w;
    XEvent *event;
{
    if (XtClass (w) == simpleMenuWidgetClass)
	XtCallActionProc (w, "notify", event, NULL, 0);
}


/************************************************************************
 * Function:	unhighlight						*
 *									*
 * Description:	Unhighlights a menu entry given by event location.	*
 ************************************************************************/

static void unhighlight (w, event)
    Widget  w;
    XEvent *event;
{
    if (XtClass (w) == simpleMenuWidgetClass)
	XtCallActionProc (w, "unhighlight", event, NULL, 0);
}


/************************************************************************
 * Function:	move							*
 *									*
 * Description:	Highlights the next or previous menu entry.		*
 ************************************************************************/

static void move (w, sign)
    Widget w;
    int    sign;
{
    Arg        args [3];
    Position   y;
    Position   min;
    Position   max;
    Position   value;
    Dimension  height;
    Dimension  width;
    WidgetList children;
    Cardinal   i;
    Cardinal   num_children;
    Widget     entry;
    Widget     next;
    XEvent     event;


    /* Make sure we have a menu and not something else. */

    if (XtClass (w) != simpleMenuWidgetClass)
	return;


    /* Get information from the menu. */

    XtSetArg (args [0], XtNwidth,       &width);
    XtSetArg (args [1], XtNchildren,    &children);
    XtSetArg (args [2], XtNnumChildren, &num_children);
    XtGetValues (w, args, 3);


    /* Determine the y coordinate of the active entry. */

    entry = XawSimpleMenuGetActiveEntry (w);

    if (entry != NULL) {
	XtSetArg (args [0], XtNy, &value);
	XtGetValues (entry, args, 1);
    } else
	value = -1;



    /* Find the widget with the least y coordinate greater than our own. */

    if (sign >= 0) {
	next = NULL;
	min = 32767;
	for (i = 0; i < num_children; i ++)
	    if (XtIsSensitive (children [i]) && XtIsManaged (children [i])) {
		XtSetArg (args [0], XtNy, &y);
		XtGetValues (children [i], args, 1);
		if (y < min && y > value) {
		    min = y;
		    next = children [i];
		}
	    }

	event.xbutton.y = min;


    /* Find the widget with the greatest y coordinate less than our own. */

    } else {
	next = NULL;
	max = -1;
	for (i = 0; i < num_children; i ++)
	    if (XtIsSensitive (children [i]) && XtIsManaged (children [i])) {
		XtSetArg (args [0], XtNy, &y);
		XtGetValues (children [i], args, 1);
		if (y > max && y < value) {
		    max = y;
		    next = children [i];
		}
	    }

	event.xbutton.y = max;
    }


    /* Nothing to do. */

    if (next == NULL || next == entry)
	return;



    /* Construct a synthetic event and highlight the entry. */

    XtSetArg (args [0], XtNheight, &height);
    XtGetValues (next, args, 1);

    event.type = ButtonPress;
    event.xbutton.x = width / 2;
    event.xbutton.y += height / 2;
    highlight (w, &event);
}


/************************************************************************
 * Function:	grab							*
 *									*
 * Description:	Grabs the keyboard and pointer.				*
 ************************************************************************/

static void grab (w, time)
    Widget w;
    Time   time;
{
    Arg    args [1];
    Widget v;
    Cursor cursor;


    if (XtClass (w) != menuButtonWidgetClass) {
	cursor = None;
	for (v = w; v && !cursor; v = XtParent (v)) {
	    XtSetArg (args [0], XtNcursor, &cursor);
	    XtGetValues (v, args, 1);
	}

	if (cursor == None)
	    cursor = default_cursor;

	XtGrabKeyboard (w, True, GrabModeAsync, GrabModeAsync, time);
	XtGrabPointer  (w, True, ButtonPressMask | ButtonReleaseMask,
			GrabModeAsync, GrabModeAsync, None, cursor, time);
    }
}


/************************************************************************
 * Function:	ungrab							*
 *									*
 * Description:	Ungrabs the keyboard and pointer			*
 ************************************************************************/

static void ungrab (w, time)
    Widget w;
    Time   time;
{
    if (XtClass (w) != menuButtonWidgetClass) {
	XtUngrabKeyboard (w, time);
	XtUngrabPointer  (w, time);
    }
}


/************************************************************************
 * Function:	inside							*
 *									*
 * Description:	Determines if the pointer is inside a rectangle.	*
 ************************************************************************/

static int inside (event)
    XEvent *event;
{
    Position x;
    Position y;


    x = event -> xbutton.x_root;
    y = event -> xbutton.y_root;

    return x >= x1 && x <= x2 && y >= y1 && y <= y2;
}


/************************************************************************
 * Function:	popup							*
 *									*
 * Description:	Pops up the menu.					*
 ************************************************************************/

static void popup (w)
    Widget w;
{
    Arg       args [4];
    Position  x;
    Position  y;
    Dimension width;
    Dimension height;


    /* Make sure we have the menu button and not the menu. */

    if (XtClass (w) != menuButtonWidgetClass)
	return;


    /* Retrieve and save the location of the button. */

    XtSetArg (args [0], XtNx,      &x);
    XtSetArg (args [1], XtNy,      &y);
    XtSetArg (args [2], XtNwidth,  &width);
    XtSetArg (args [3], XtNheight, &height);
    XtGetValues (w, args, 4);

    XtTranslateCoords (XtParent (w), x, y, &x1, &y1);
    XtTranslateCoords (XtParent (w), x + width, y + height, &x2, &y2);


    /* Pop up the menu. */

    XtCallActionProc (w, "PopupMenu", NULL, NULL, 0);
    posted = True;
    from_inside = False;
}


/************************************************************************
 * Function:	popdown							*
 *									*
 * Description:	Pops down the menu.					*
 ************************************************************************/

static void popdown (w)
    Widget w;
{
    if (XtClass (w) == simpleMenuWidgetClass)
	XtOverrideTranslations (w, old_translations);

    XtPopdown (w);
    from_inside = False;
    posted = False;
}


/************************************************************************
 * Function:	nop							*
 *									*
 * Description:	Empty action procedure.					*
 ************************************************************************/

static void nop (w, event, params, num_params)
    Widget    w;
    XEvent   *event;
    String   *params;
    Cardinal *num_params;
{
}


/************************************************************************
 * Function:	action							*
 *									*
 * Description:	Action procedure to handle events occuring on a		*
 *		postable menu.						*
 ************************************************************************/

void action (w, event, params, num_params)
    Widget    w;
    XEvent   *event;
    String   *params;
    Cardinal *num_params;
{
    Time      now;
    KeySym    keysym;
    Modifiers modifiers;


    printf ("%s ", XtName (w));

    switch (event -> type) {
    case KeyPress:
	printf ("key ");
	now = event -> xkey.time;
	keysym = XtGetActionKeysym (event, &modifiers);

	switch (keysym) {
	case XK_Return:
	case XK_space:
	    if (!posted && when != now) {
		printf ("popup");
		grab (w, now);
		popup (w);
		when = now;
	    } else if (posted && when != now) {
		printf ("popdown");
		ungrab (w, now);
		popdown (w);
		notify (w, event);
		unhighlight (w, event);
	    } else if (when == now) {
		printf ("here");
		grab (w, now);
		if (XtClass (w) == simpleMenuWidgetClass)
		    XtOverrideTranslations (w, new_translations);
	    } else
		printf ("%d %d\n", when, now);
	    break;


	case XK_Down:
	    if (posted) {
		move (w, 1);
		break;
	    }


	case XK_Up:
	    if (posted) {
		move (w, -1);
		break;
	    }

	default:

	    if (!posted && when != now) {
		printf ("popup");
		grab (w, now);
		popup (w);
		when = now;
	    } else if (posted && when != now) {
		printf ("popdown");
		ungrab (w, now);
		popdown (w);
		unhighlight (w, event);
	    } else if (when == now) {
		printf ("here");
		grab (w, now);
		if (XtClass (w) == simpleMenuWidgetClass)
		    XtOverrideTranslations (w, new_translations);
	    }
	    break;
	}
	break;


    case ButtonPress:
	printf ("down ");
	now = event -> xbutton.time;

	if (XtClass (w) == menuButtonWidgetClass)
	    from_button = True;

	if (posted && when != now) {
	    printf ("highlight");
	    highlight (w, event);

	} else if (!posted) {
	    printf ("popup");
	    grab (w, now);
	    popup (w);
	    when = now;
	}
	break;


    case ButtonRelease:
	printf ("up ");
	now = event -> xbutton.time;

	if (from_button && !from_inside && inside (event)) {
	    printf ("inside");
	    XtOverrideTranslations (w, new_translations);
	    grab (w, now);
	    when = now;
	    from_inside = True;

	} else if (when != now) {
	    printf ("popdown");
	    ungrab (w, now);
	    popdown (w);
	    notify (w, event);
	    unhighlight (w, event);
	}
	break;
    }

    printf ("\n");
}


/************************************************************************
 * Function:	AddPostMenuActions					*
 *									*
 * Description:	Adds the postable menu actions to the specified widget.	*
 ************************************************************************/

void AddPostMenuActions (w)
    Widget w;
{
    static XtAppContext	app_context = NULL;
    static XtActionsRec	actions [ ] = {{"PostMenu", action}, {"no-op", nop}};


    if (app_context == NULL) {
	app_context = XtWidgetToApplicationContext (w);
	XtAppAddActions (app_context, actions, XtNumber (actions));
	post_translations = XtParseTranslationTable (post_table);
	new_translations = XtParseTranslationTable (new_table);
	old_translations = XtParseTranslationTable (old_table);
	default_cursor = XCreateFontCursor (XtDisplay (w), XC_top_left_arrow);
    }

    XtOverrideTranslations (w, post_translations);
}
