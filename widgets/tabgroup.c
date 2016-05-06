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
 * File:	tabgroup.c						*
 *									*
 * Description:	This file contains the private and public function and	*
 *		type definitions for the tab group mechanism.		*
 ************************************************************************/

# include <stdio.h>
# include <X11/Intrinsic.h>
# include <X11/StringDefs.h>
# include <X11/Xmu/Drawing.h>
# include <X11/Xaw/Text.h>
# include "TabGroup.h"

# ifndef XtNshadowWidth
# define XtNshadowWidth "shadowWidth"
# endif

# if NeedWidePrototypes
# define BOOLEAN int
# else
# define BOOLEAN Boolean
# endif


# define MaxTabGroups 64
# define TranslationTable \
"None<Key>Tab:    ChangeFocus(%d,next)\n\
 !Shift<Key>Tab:  ChangeFocus(%d,prev)\n\
 <ClientMessage>: ChangeFocus(%d,request)"


typedef struct tab_group {
    Widget     ancestor;
    WidgetList members;
    Cardinal   num_members;
    Cardinal   focused;
    Pixel      highlight;
    Pixel      border_color;
    Pixmap     border_pixmap;
    Boolean    clear_caret;
} *TabGroup;

static TabGroup group_list [MaxTabGroups];
static Cardinal num_groups = 0;
static Widget	focused_widget;


/************************************************************************
 * Function:	ChangeFocusProc						*
 *									*
 * Description:	Changes the focused widget in a tab group by restoring	*
 *		the original appearance of the widget, determining the	*
 *		newly focused widget, saving its appearance, and	*
 *		highlighting the widget.				*
 ************************************************************************/

static void ChangeFocusProc (w, event, params, num_params)
    Widget    w;
    XEvent   *event;
    String   *params;
    Cardinal *num_params;
{
    Arg      args [2];
    Cardinal offset;
    Cardinal focused;
    Cardinal group_number;
    TabGroup group;


    sscanf (params [0], "%u", &group_number);
    group = group_list [group_number];
    focused = group -> focused;


    /* Determine if a request was made. */

    if (!strcmp (params [1], "request"))
	if (!strcmp (event -> xclient.data.b, "get")) {
	    focused_widget = group -> members [focused];
	    return;
	}


    /* Restore the original appearance of the widget if necessary. */

    if (focused != group -> num_members) {
	XtSetArg (args [0], XtNborderColor,  group -> border_color);
	XtSetArg (args [1], XtNborderPixmap, group -> border_pixmap);
	XtSetValues (group -> members [focused], args, XtNumber (args));

	if (group -> clear_caret == True) {
	    XtSetArg (args [0], XtNdisplayCaret, False);
	    XtSetValues (group -> members [focused], args, 1);
	}
    }


    /* Determine the index of the newly focused widget. */

    offset = 1;
    if (!strcmp (params [1], "prev"))
	offset = group -> num_members - 1;

    if (!strcmp (params [1], "next"))
	focused = (focused + 1) % group -> num_members;
    else if (!strcmp (params [1], "prev"))
	focused = (group -> num_members + focused - 1) % group -> num_members;
    else
	for (focused = 0; focused < group -> num_members; focused ++)
	    if (w == group -> members [focused])
		break;


    /* Move off the widget if it is not sensitive or managed. */

    while (!XtIsSensitive (group -> members [focused]) ||
	   !XtIsManaged (group -> members [focused]))
	focused = (focused + offset) % group -> num_members;


    /* Save the original appearance. */

    XtSetArg (args [0], XtNborderColor,  &group -> border_color);
    XtSetArg (args [1], XtNborderPixmap, &group -> border_pixmap);
    XtGetValues (group -> members [focused], args, XtNumber (args));


    /* Highlight the widget. */

    XtSetArg (args [0], XtNborderColor,  group -> highlight);
    XtSetArg (args [1], XtNborderPixmap, XtUnspecifiedPixmap);
    XtSetValues (group -> members [focused], args, XtNumber (args));

    if (group -> clear_caret == True) {
	XtSetArg (args [0], XtNdisplayCaret, True);
	XtSetValues (group -> members [focused], args, 1);
    }

    XtSetKeyboardFocus (group -> ancestor, group -> members [focused]);
    group -> focused = focused;
}


/************************************************************************
 * Function:	SetFocusProc						*
 *									*
 * Description:	Sets the focus to the specified widget by simulating	*
 *		an appropriate client message event.				*
 ************************************************************************/

static void SetFocusProc (w, event, params, num_params)
    Widget    w;
    XEvent   *event;
    String   *params;
    Cardinal *num_params;
{
    XEvent synthetic_event;


    synthetic_event.type            = ClientMessage;
    synthetic_event.xclient.display = XtDisplay (w);
    synthetic_event.xclient.window  = XtWindow (w);
    synthetic_event.xclient.format  = 8;

    strcpy (synthetic_event.xclient.data.b, "set");
    XtDispatchEvent (&synthetic_event);
}


/************************************************************************
 * Function:	CreateTabGroup						*
 *									*
 * Description:	Creates a new tab group by adding actions on each	*
 *		member of group and initializes the appearance of each	*
 *		widget.							*
 ************************************************************************/

void CreateTabGroup (ancestor, members, num_members, highlight, clear_caret)
    Widget     ancestor;
    WidgetList members;
    Cardinal   num_members;
    Pixel      highlight;
    BOOLEAN    clear_caret;
{
    int			depth;
    char		table [256];
    Arg			args [5];
    Pixel		border_color;
    Pixel		background;
    Pixel		parent_background;
    Pixmap		border_pixmap;
    Dimension		border_width;
    Dimension		shadow_width;
    Cardinal		i;
    TabGroup		group;
    XtTranslations	translations;
    static XtAppContext	app_context = NULL;
    static XtActionsRec	actions [ ] =
		{{"ChangeFocus", ChangeFocusProc}, {"SetFocus", SetFocusProc}};


    /* Add the actions to the application context if necessary. */

    if (app_context == NULL) {
	app_context = XtWidgetToApplicationContext (ancestor);
	XtAppAddActions (app_context, actions, XtNumber (actions));
    }


    /* Create and initialize the tab group structure. */

    group = XtNew (struct tab_group);
    group -> members = (WidgetList) XtMalloc (sizeof (Widget) * num_members);
    group -> num_members = num_members;
    group -> focused = num_members;
    group -> ancestor = ancestor;
    group -> clear_caret = clear_caret;


    /* Create and parse the translations. */

    sprintf (table, TranslationTable, num_groups, num_groups, num_groups);
    translations = XtParseTranslationTable (table);


    /* Configure each widget. */

    for (i = 0; i < num_members; i ++) {
	group -> members [i] = members [i];
	XtOverrideTranslations (members [i], translations);


	/* Retrieve information about the widget's appearance. */

	shadow_width = 0;
	border_pixmap = XtUnspecifiedPixmap;

	XtSetArg (args [0], XtNdepth,       &depth);
	XtSetArg (args [1], XtNbackground,  &background);
	XtSetArg (args [2], XtNborderColor, &border_color);
	XtSetArg (args [3], XtNborderWidth, &border_width);
	XtSetArg (args [4], XtNshadowWidth, &shadow_width);
	XtGetValues (members [i], args, 5);

	XtSetArg (args [0], XtNbackground, &parent_background);
	XtGetValues (XtParent (members [i]), args, 1);


	/* If the widget is monochrome then we don't have a lot of options.
	   We need to make sure the highlight color is visible.  If the
	   widget is 3d and doesn't have a border then make its border color
	   the same as its parent's background color as below.  Otherwise,
	   give the widget a stippled pixmap as a border so that it can be
	   distinguished from its parent and visible when highlighted. */

	if (depth == 1) {
	    highlight = !background;
	    if (shadow_width != 0 && border_width == 0)
		border_color = parent_background;
	    else
		border_pixmap = XmuCreateStippledPixmap (XtScreen (members [i]),
				!background, background, depth);


	/* The widget is color.  If it has no border then make its border color
	   its parent's background color so it "appears" to have no border. */

	} else if (border_width == 0)
	    border_color = parent_background;


	/* Set the final appearance of the widget. */

	XtSetArg (args [0], XtNborderWidth,  2);
	XtSetArg (args [1], XtNborderColor,  border_color);
	XtSetArg (args [2], XtNborderPixmap, border_pixmap);
	XtSetValues (members [i], args, 3);


	/* Clear the caret if necessary. */

	if (group -> clear_caret == True) {
	    XtSetArg (args [0], XtNdisplayCaret, False);
	    XtSetValues (members [i], args, 1);
	}
    }


    /* Add the group to the list of tab groups. */

    group_list [num_groups ++] = group;
    group -> highlight = highlight;
}


/************************************************************************
 * Function:	SetFocus						*
 *									*
 * Description:	Sets the keyboard focus to the specified widget by	*
 *		calling the SetFocus() action procedure.		*
 ************************************************************************/

void SetFocus (member)
    Widget member;
{
    XtCallActionProc (member, "SetFocus", NULL, NULL, 0);
}


/************************************************************************
 * Function:	GetFocus						*
 *									*
 * Description:	Returns which widget has the keyboard focus.		*
 ************************************************************************/

Widget GetFocus (member)
    Widget member;
{
    XEvent synthetic_event;


    synthetic_event.type            = ClientMessage;
    synthetic_event.xclient.display = XtDisplay (member);
    synthetic_event.xclient.window  = XtWindow (member);
    synthetic_event.xclient.format  = 8;

    strcpy (synthetic_event.xclient.data.b, "get");
    XtDispatchEvent (&synthetic_event);
    return focused_widget;
}


/************************************************************************
 * Function:	HasFocus						*
 *									*
 * Description:	Determines if the specified widget has the keyboard	*
 *		focus.							*
 ************************************************************************/

Boolean HasFocus (member)
    Widget member;
{
    return GetFocus (member) == member ? True : False;
}
