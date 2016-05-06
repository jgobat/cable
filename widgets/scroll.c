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
 * File:	scroll.c						*
 *									*
 * Description:	This file contains the public and private function	*
 *		and variable declarations for the scrollable text	*
 *		mechanism.						*
 ************************************************************************/

# include <stdio.h>
# include <X11/IntrinsicP.h>
# include <X11/StringDefs.h>
# include <X11/Xaw/TextP.h>
# include "scroll.h"
# include "util.h"


static String scroll_table =
"<KeyUp>:           ScrollText()\n\
 <ConfigureNotify>: ScrollText()";

static XtTranslations scroll_translations;


/************************************************************************
 * Function:	ScrollText						*
 *									*
 * Description:	Action procedure interface to ScrollToInsertionPoint().	*
 ************************************************************************/

static void ScrollText (w, event, params, num_params)
    Widget    w;
    XEvent   *event;
    String   *params;
    Cardinal *num_params;
{
    ScrollToInsertionPoint (w);
}


/************************************************************************
 * Function:	AddScrollableTextTranslations				*
 *									*
 * Description:	Adds translations to the specified text widget to	*
 *		enable a primitive character-based scrolling.		*
 ************************************************************************/

void AddScrollableTextTranslations (text)
    Widget text;
{
    static XtAppContext	app_context = NULL;
    static XtActionsRec	actions [ ] = {{"ScrollText", ScrollText}};


    /* Add the actions to the application context if necessary. */

    if (app_context == NULL) {
	app_context = XtWidgetToApplicationContext (text);
	XtAppAddActions (app_context, actions, XtNumber (actions));
	scroll_translations = XtParseTranslationTable (scroll_table);
    }

    XtVaSetValues (text, XtNrightMargin, 0, NULL);
    XtOverrideTranslations (text, scroll_translations);


    /* The SetValues really should do this for us. */

    ((TextWidget) text) -> text.margin.right = 0;
}


/************************************************************************
 * Function:	ScrollToInsertionPoint					*
 *									*
 * Description:	Scrolls the text so the insertion point is always	*
 *		visible.  The scrolling is accomplished by updating the	*
 *		display position of the string and is rather primitive.	*
 ************************************************************************/

void ScrollToInsertionPoint (text)
    Widget text;
{
    int		 ins_pos;
    int		 disp_pos;
    int		 length;
    String	 ptr;
    String	 value;
    Boolean	 caret;
    Position	 l_margin;
    Position	 r_margin;
    Dimension	 width;
    XFontStruct	*font;


    XtVaGetValues (text, XtNdisplayPosition, &disp_pos,
			 XtNinsertPosition,  &ins_pos,
			 XtNdisplayCaret,    &caret,
			 XtNrightMargin,     &r_margin,
			 XtNleftMargin,      &l_margin,
			 XtNstring,	     &value,
			 XtNwidth,	     &width,
			 XtNfont,	     &font, NULL);


    /* Check if the translations have been installed. */

    if (r_margin != 0)
	return;


    /* Scroll the text string to the right. */

    if (ins_pos <= disp_pos) {
	disp_pos = ins_pos - 1;
	if (disp_pos < 0)
	    disp_pos = 0;


    /* Scroll the text string to the left. */

    } else {
	width -= l_margin + r_margin;
	width -= GetTextWidth (font, "W", 1);	/* the widest character? */

	while (1) {
	    ptr = &value [disp_pos];
	    length = 0;

	    while (ptr [length] && GetTextWidth (font, ptr, length) < width)
		length ++;

	    if (ins_pos > disp_pos + length)
		disp_pos ++;
	    else
		break;
	}
    }


    /* Scroll to the new position, disabling the caret temporarily. */

    if (caret == True)
	XtVaSetValues (text, XtNdisplayCaret, False, NULL);

    XtVaSetValues (text, XtNdisplayPosition, disp_pos, NULL);

    if (caret == True)
	XtVaSetValues (text, XtNdisplayCaret, True, NULL);
}
