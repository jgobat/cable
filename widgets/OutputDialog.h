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
 * File:	OutputDialog.h						*
 *									*
 * Description:	This file contains the public function and type		*
 *		declarations for the output dialog box.			*
 ************************************************************************/

# ifndef _OutputDialog_h
# define _OutputDialog_h

# if NeedVarargsPrototypes
# include <stdarg.h>
# else
# include <varargs.h>
# endif

typedef struct output_dialog *OutputDialog;


extern OutputDialog OutputDialogCreate (
# if NeedFunctionPrototypes
    Widget			/* parent      */,
    String			/* name        */,
    String *			/* buttons     */,
    Cardinal			/* num_buttons */
# endif
);


extern String OutputDialogSelect (
# if NeedFunctionPrototypes
    OutputDialog		/* output_dialog */,
    String			/* title	 */,
    String			/* preferred	 */
# endif
);


extern void OutputDialogPopup (
# if NeedFunctionPrototypes
    OutputDialog		/* output_dialog */,
    String			/* title	 */,
    String			/* preferred	 */,
    XtCallbackProc		/* callback	 */,
    XtPointer			/* client_data	 */
# endif
);


extern void OutputDialogPopdown (
# if NeedFunctionPrototypes
    OutputDialog		/* output_dialog */
# endif
);


extern void OutputDialogView (
# if NeedFunctionPrototypes
    OutputDialog		/* output_dialog */,
    String			/* file_name	 */,
    Cardinal			/* max_lines	 */,
    Cardinal			/* max_columns	 */
# endif
);


extern void OutputDialogPrintf (
# if NeedVarargsPrototypes
    OutputDialog		/* output_dialog */,
    String			/* format	 */,
    ...				/* arguments	 */
# endif
);


extern void OutputDialogVprintf (
# if NeedFunctionPrototypes
    OutputDialog		/* output_dialog */,
    String			/* format	 */,
    va_list			/* argument list */
# endif
);


extern Widget OutputDialogShell (
# if NeedFunctionPrototypes
    OutputDialog		/* output_dialog */
# endif
);

# endif /* _OutputDialog_h */
