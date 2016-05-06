/*
    This file is part of WHOI Cable, a program for the static and dynamic
    analysis of oceanographic cable structures.

    Copyright (C) 2008-2016 by Jason Gobat

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

#include <gtk/gtk.h>
#include "duplicates.h"

extern GtkWidget *toplevel;

int gfDuplicatePref;

int
DuplicateDialog(char *name, char *txt)
{
    GtkWidget   *dlg;
    int          ans;

    dlg = gtk_message_dialog_new(GTK_WINDOW(toplevel), 0, GTK_MESSAGE_QUESTION,
                                 GTK_BUTTONS_NONE, 
                                 "%s named '%s' already exists in "
                                 "the global database. Would you like to "
                                 "rename this local version or use the "
                                 "version in in the database?", txt, name);
    gtk_dialog_add_buttons(GTK_DIALOG(dlg), 
                           "Rename", RENAME, 
                           "Rename all", RENAME_ALL, 
                           "Use global", USE_GLOBAL, 
                           "Use global for all", USE_GLOBAL_ALL, NULL);

    ans = gtk_dialog_run (GTK_DIALOG (dlg));
    gtk_widget_destroy (dlg);

    return ans;
}
