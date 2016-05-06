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
 * File:	outputdialog.c						*
 *									*
 * Description:	This file contains the public and private function and	*
 *		type definitions for the output dialog box.		*
 ************************************************************************/

# include <stdio.h>
# include <X11/Intrinsic.h>
# include <X11/StringDefs.h>
# include <X11/Shell.h>
# include <X11/Xaw/AsciiText.h>
# include <X11/Xaw/Command.h>
# include "Layout.h"
# include "OutputDialog.h"
# include "TabGroup.h"
# include "util.h"

# if NeedVarargsPrototypes
# include <stdarg.h>
# define Va_start(a,b) va_start(a,b)
# else
# include <varargs.h>
# define Va_start(a,b) va_start(a)
# endif

# define MaxButtons 4
# define Waiting ((String) 1)


struct output_dialog {
    Widget	   shell;
    Widget	   layout;
    Widget	   text;
    Widget	   button [MaxButtons];
    Widget	   dummy;
    Cardinal	   num_buttons;
    XtCallbackProc callback;
    XtPointer	   client_data;
    String	   selected;
};


/* Resources */

static Pixel highlight;

static String layout_string = "\
vertical { \
    4 \
    horizontal { \
	4 \
	text <+inf -100% * +inf -100%> \
	4 \
    } \
    4 \
    horizontal { \
	4 \
	button1 4 button2 4 button3 4 button4 \
	4 \
    } \
    4 \
}";

static Arg color_args [ ] = {
    {XtNborderColor, (XtArgVal) &highlight},
};

static Arg layout_args [ ] = {
    {XtNlayout, (XtArgVal) NULL},
};

static Arg text_args [ ] = {
    {XtNeditType,	(XtArgVal) XawtextRead},
    {XtNwrap,		(XtArgVal) XawtextWrapWord},
    {XtNscrollVertical, (XtArgVal) XawtextScrollWhenNeeded},
    {XtNborderWidth,	(XtArgVal) 0},
    {XtNwidth,		(XtArgVal) 800},
};

static Arg command_args [ ] = {
    {XtNlabel, (XtArgVal) ""},
};


/* Translations */

static String command_table =
"<Key>Return:   AutoRepeat(off) set()\n\
 <KeyUp>Return: AutoRepeat(saved) notify() unset()\n\
 <Key>space:    AutoRepeat(off) set()\n\
 <KeyUp>space:  AutoRepeat(saved) notify() unset()\n\
 <Key>Escape:   OutputDialogDelete()";

static XtTranslations command_translations;


/************************************************************************
 * Function:	Delete							*
 *									*
 * Description:	An action procedure which is called when the dialog is	*
 *		to be deleted; calls the callback of the dummy widget.	*
 ************************************************************************/

static void Delete (w, event, params, num_params)
    Widget    w;
    XEvent   *event;
    String   *params;
    Cardinal *num_params;
{
    while (XtClass (w) != transientShellWidgetClass)
	w = XtParent (w);

    w = XtNameToWidget (w, "layout.dummy");

    XtCallCallbacks (w, XtNcallback, NULL);
}


/************************************************************************
 * Function:	Notify							*
 *									*
 * Description:	Notifies the user through a callback or selection.	*
 ************************************************************************/

static void Notify (w, client_data, call_data)
    Widget    w;
    XtPointer client_data;
    XtPointer call_data;
{
    OutputDialog outputd;
    String	 label;


    outputd = (OutputDialog) client_data;

    label = strcmp (XtName (w), "dummy") ? GetLabelString (w) : NULL;

    if (outputd -> callback != NULL)
	outputd -> callback (w, outputd -> client_data, label);
    else
	outputd -> selected = label;
}


/************************************************************************
 * Function:	OutputDialogCreate					*
 *									*
 * Description:	Creates and returns a new output dialog box.  The	*
 *		dialog can be popped up using OutputDialogSelect()	*
 *		which forces the user to make a selection, or by using	*
 *		OutputialogPopUp() which will call a callback when a	*
 *		button is selected or the dialog is deleted.  The	*
 *		dialog can then be popped down using OutputDialog-	*
 *		Popdown().						*
 ************************************************************************/

OutputDialog OutputDialogCreate (parent, name, buttons, num_buttons)
    Widget   parent;
    String   name;
    String  *buttons;
    Cardinal num_buttons;
{
    Arg			args [1];
    Widget		group [MaxButtons + 1];
    char		buffer [32];
    Cardinal		i;
    OutputDialog	outputd;
    XFontStruct	       *font;
    static XtAppContext	app_context = NULL;
    static XtActionsRec	actions [ ] = {{"OutputDialogDelete", Delete}};


    /* Perform one time initialization. */

    if (app_context == NULL) {
	app_context = XtWidgetToApplicationContext (parent);
	XtAppAddActions (app_context, actions, XtNumber (actions));
	AddAutoRepeatAction (app_context);
	command_translations = XtParseTranslationTable (command_table);
	layout_args [0].value = StringToLayout (parent, layout_string);
    }


    /* Create the output dialog and its widgets. */

    outputd = XtNew (struct output_dialog);

    XtSetArg (args [0], XtNallowShellResize, True);
    outputd -> shell  = XtCreatePopupShell (name,
			transientShellWidgetClass, parent,
			args, 1);

    outputd -> layout = XtCreateManagedWidget ("layout",
			layoutWidgetClass, outputd -> shell,
			layout_args, XtNumber (layout_args));

    outputd -> text   = XtCreateManagedWidget ("text",
			asciiTextWidgetClass, outputd -> layout,
			text_args, XtNumber (text_args));

    outputd -> dummy  = XtCreateWidget ("dummy",
			commandWidgetClass, outputd -> layout,
			NULL, 0);


    XtSetArg (args [0], XtNfont, &font);
    XtGetValues (outputd -> text, args, 1);
    XtSetArg (args [0], XtNfont, font);
    XtSetValues (outputd -> dummy, args, 1);

    XtSetArg (args [0], XtNdisplayCaret, False); 
    XtSetValues (outputd -> text, args, 1);

    if (num_buttons > MaxButtons)
	num_buttons = MaxButtons;

    outputd -> num_buttons = num_buttons;

    for (i = 0; i < num_buttons; i ++) {
	sprintf (buffer, "button%u", i + 1);
	XtSetArg (command_args [0], XtNlabel, buttons [i]);
	outputd -> button [i] = XtCreateManagedWidget (XtNewString (buffer),
				commandWidgetClass, outputd -> layout,
				command_args, XtNumber (command_args));

	group [i] = outputd -> button [i];
    }


    /* Create a tab group for the output dialog. */

/*
    group [i] = outputd -> text;
    XtGetValues (outputd -> layout, color_args, XtNumber (color_args));
    CreateTabGroup (outputd -> shell, group, num_buttons + 1, highlight, True);
*/
    XtGetValues (outputd -> layout, color_args, XtNumber (color_args));
    CreateTabGroup (outputd -> shell, group, num_buttons, highlight, True);
    XtRealizeWidget (outputd -> shell);


    /* Add the translations to each widget. */

    AddDeleteWindowProtocol (outputd -> shell, "OutputDialogDelete()");

    for (i = 0; i < num_buttons; i ++)
	XtOverrideTranslations (outputd -> button [i], command_translations);


    /* Add the necessary callbacks. */

    XtAddCallback (outputd -> dummy, XtNcallback, Notify, outputd);

    for (i = 0; i < num_buttons; i ++)
	XtAddCallback (outputd -> button [i], XtNcallback, Notify, outputd);

    return outputd;
}


/************************************************************************
 * Function:	OutputDialogSelect					*
 *									*
 * Description:	Pops up the output dialog with the specified shell	*
 *		title, sets the input focus to the preferred button,	*
 *		and waits until a selection is made or the dialog is	*
 *		deleted; answer will point to the selection or NULL if	*
 *		the dialog was deleted.					*
 ************************************************************************/

String OutputDialogSelect (outputd, title, preferred)
    OutputDialog outputd;
    String	 title;
    String	 preferred;
{
    XEvent	 event;
    XtAppContext app_context;


    OutputDialogPopup (outputd, title, preferred, NULL, NULL);

    outputd -> selected = Waiting;
    app_context = XtWidgetToApplicationContext (outputd -> shell);

    while (outputd -> selected == Waiting) {
	XtAppNextEvent (app_context, &event);
	XtDispatchEvent (&event);
    }

    OutputDialogPopdown (outputd);

    return outputd -> selected;
}


/************************************************************************
 * Function:	OutputDialogPopup					*
 *									*
 * Description:	Pops up the output dialog with the specified shell	*
 *		title and sets the input focus to the preferred button.	*
 *		If the callback is not NULL then it will be called upon	*
 *		selection or deletion with the specified client	data.	*
 ************************************************************************/

void OutputDialogPopup (outputd, title, preferred, callback, client_data)
    OutputDialog   outputd;
    String	   title;
    String	   preferred;
    XtCallbackProc callback;
    XtPointer	   client_data;
{
    Arg      args [2];
    Cardinal i;


    XtSetArg (args [0], XtNtitle,    title);
    XtSetArg (args [1], XtNiconName, title);
    XtSetValues (outputd -> shell, args, 2);

    outputd -> callback = callback;
    outputd -> client_data = client_data;

    for (i = 0; i < outputd -> num_buttons; i ++)
	if (!strcmp (GetLabelString (outputd -> button [i]), preferred))
	    SetFocus (outputd -> button [i]);

    XtPopup (outputd -> shell, callback == NULL ? XtGrabExclusive : XtGrabNone);
}


/************************************************************************
 * Function:	OutputDialogPopdown					*
 *									*
 * Description:	Pops down the specified output dialog.			*
 ************************************************************************/

void OutputDialogPopdown (outputd)
    OutputDialog outputd;
{
    XtPopdown (outputd -> shell);
}


/************************************************************************
 * Function:	OutputDialogView					*
 *									*
 * Description:	...							*
 ************************************************************************/

void OutputDialogView (outputd, file_name, max_lines, max_columns)
    OutputDialog outputd;
    String	 file_name;
    Cardinal	 max_lines;
    Cardinal	 max_columns;
{
    Arg   args [4];
    Cardinal length;
    Cardinal num_lines;
    Cardinal num_columns;
    XFontStruct *font;
    Dimension width;
    Dimension height;
    FILE *fp;
    int   n;

    if ((fp = fopen (file_name, "r")) == NULL)
	return;
   
    num_lines = 0;
    num_columns = 0;
    
    while ((n = fscanf (fp, " %*[^\n]%n%*[\n]", &length)) != EOF) {
	if (length > num_columns)
	    num_columns = length;
        
	num_lines ++;
        if (feof(fp))
           break;
    }
   
    fclose (fp);
    
    if (max_lines && num_lines > max_lines)
	num_lines = max_lines;

    if (max_columns && num_columns > max_columns)
	num_columns = max_columns;

    XtSetArg (args [0], XtNfont, &font);
    XtGetValues (outputd -> text, args, 1);
    
    width = font -> max_bounds.width;
    height = font -> max_bounds.ascent + font -> max_bounds.descent;

    XtSetArg (args [0], XtNtype,   XawAsciiFile);
    XtSetArg (args [1], XtNstring, file_name);
    XtSetArg (args [2], XtNwidth,  width * (num_columns + 1) + 22);
    XtSetArg (args [3], XtNheight, height * (num_lines + 1));
    XtSetValues (outputd -> text, args, 4);
}


/************************************************************************
 * Function:	OutputDialogPrintf					*
 *									*
 * Description: Uses OutputDialogVprintf() to simulate printf().	*
 ************************************************************************/

# if NeedVarargsPrototypes
void OutputDialogPrintf (OutputDialog outputd, String format, ...)
# else
void OutputDialogPrintf (outputd, format, va_alist)
    OutputDialog outputd;
    String	 format;
    va_dcl
# endif
{
    va_list	     ap;


    Va_start (ap, format);
    OutputDialogVprintf (outputd, format, ap);
    va_end (ap);
}


/************************************************************************
 * Function:	OutputDialogVprintf					*
 *									*
 * Description:	Acts like vprintf() only writes to the text widget.  We	*
 *		need to force the size of the text widget so that the	*
 *		string can be seen.  Ideally, we would query the	*
 *		preferred geometry of the text widget, but alas this	*
 *		doesn't work as the text widget doesn't have a 		*
 *		preferred size since it just assumes that we'll use its	*
 *		built-in scrollbars and thus happily conforms to	*
 *		whatever size it has (yes, I know that resize should	*
 *		work but it doesn't seem to and it doesn't affect the	*
 *		preferred size).  So instead we use the dummy widget to	*
 *		determine the proper size.  This works since command	*
 *		(label) widgets return their preferred size and handle	*
 *		newlines correctly.  Ideally, we would want to use a	*
 *		label widget rather than a text widget, but alas, the	*
 *		text widget allows cutting and pasting and more		*
 *		importantly allows a file to viewed efficiently.	*
 ************************************************************************/

void OutputDialogVprintf (outputd, format, ap)
    OutputDialog outputd;
    String	 format;
    va_list	 ap;
{
    Arg		     args [4];
    char	     buffer [2048];
    XtWidgetGeometry preferred;


    /* Set the label of the dummy widget to the desired string to that
       we can determine its preferred size.  The text widget doesn't
       seem to have a preferred size so we have to do it this way. */

    vsprintf (buffer, format, ap);
    SetLabelString (outputd -> dummy, buffer);
    XtQueryGeometry (outputd -> dummy, NULL, &preferred);

    XtSetArg (args [0], XtNtype,   XawAsciiString);
    XtSetArg (args [1], XtNstring, buffer);
    XtSetArg (args [2], XtNwidth,  preferred.width + 20);
    XtSetArg (args [3], XtNheight, preferred.height + 4);
    XtSetValues (outputd -> text, args, 4);
}


/************************************************************************
 * Function:	OutputDialogShell					*
 *									*
 * Description:	Returns the shell widget of an output dialog.		*
 ************************************************************************/

Widget OutputDialogShell (outputd)
    OutputDialog outputd;
{
    return outputd -> shell;
}
