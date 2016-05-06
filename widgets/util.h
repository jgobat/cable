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
 * File:	util.h							*
 *									*
 * Description:	This file contains the function declarations for the	*
 *		miscellaneous utility functions.			*
 ************************************************************************/

# ifndef _UTIL_H
# define _UTIL_H

extern XtArgVal StringToLayout (
# if NeedFunctionPrototypes
    Widget		/* widget */,
    String		/* string */
# endif
);


extern void CenterOnWidget (
# if NeedFunctionPrototypes
    Widget		/* shell  */,
    Widget		/* center */,
# if NeedWidePrototypes
    int			/* force */
# else
    Boolean		/* force */
# endif
# endif
);


extern void CenterOnScreen (
# if NeedFunctionPrototypes
    Widget		/* shell */,
# if NeedWidePrototypes
    int			/* force */
# else
    Boolean		/* force */
# endif
# endif
);


extern void AddDeleteWindowProtocol (
# if NeedFunctionPrototypes
    Widget		/* shell  */,
    String		/* action */
# endif
);


extern void WarpToCenter (
# if NeedFunctionPrototypes
    Widget		/* w */
# endif
);


extern void ListAddCursorTranslations (
# if NeedFunctionPrototypes
    Widget		/* viewport */
# endif
);


extern void ListAddCursorAccelerators (
# if NeedFunctionPrototypes
    Widget		/* viewport */,
    Widget		/* w        */
# endif
);


extern void SetTextString (
# if NeedFunctionPrototypes
    Widget		/* text  */,
    String		/* value */
# endif
);


extern String GetTextString (
# if NeedFunctionPrototypes
    Widget		/* text */
# endif
);


extern Cardinal GetTextWidth (
# if NeedFunctionPrototypes
    XFontStruct *	/* font   */,
    String		/* text   */,
    Cardinal		/* length */
# endif
);


extern void SetLabelString (
# if NeedFunctionPrototypes
    Widget		/* label */,
    String		/* value */
# endif
);


extern String GetLabelString (
# if NeedFunctionPrototypes
    Widget		/* label */
# endif
);


extern void AddAutoRepeatAction (
# if NeedFunctionPrototypes
    XtAppContext	/* app_context */
# endif
);


extern Widget CreateHelpButton (
# if NeedFunctionPrototypes
    Widget		/* parent */,
    String		/* name   */
# endif
);


extern void UpdateHelpMessage (
# if NeedFunctionPrototypes
    Widget		/* button  */,
    String		/* message */,
# if NeedWidePrototypes
    unsigned		/* width   */
# else
    Dimension		/* width   */
# endif
# endif
);

# endif /* _UTIL_H */
