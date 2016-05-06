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
#include <stdlib.h>
#include "plot.h"

static gboolean
delete_dialog(GtkWidget *w, gpointer *data)
{
    return TRUE;
}

EntriesDialog *
CreateEntriesDialog(GtkWidget *toplevel, char *title, char **fields)
{
   EntriesDialog *ed;
   GtkWidget *label, *hbox;
   char      **ptr;
   int         i;
   int         num_fields;

   ed = (EntriesDialog *) malloc(sizeof(EntriesDialog));
   
   ed -> dialog = gtk_dialog_new_with_buttons(title,
                                        GTK_WINDOW(toplevel),
                                        GTK_DIALOG_DESTROY_WITH_PARENT,
                                        GTK_STOCK_OK, GTK_RESPONSE_OK,
                                        GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                        NULL);


    num_fields = 0;
    ptr = fields;
    while(*ptr) {
        num_fields ++;
        ptr ++;
    }

    ed -> field = (GtkWidget **) malloc(sizeof(GtkWidget *)*num_fields);

    for (i = 0 ; i < num_fields ; i++) {
        hbox = gtk_hbox_new(FALSE, 0);
        label =  gtk_label_new(fields[i]);
        ed -> field[i] = gtk_entry_new();
        gtk_entry_set_max_length(GTK_ENTRY(ed -> field[i]), 16);
        gtk_container_add(GTK_CONTAINER(hbox), label);
        gtk_container_add(GTK_CONTAINER(hbox), ed -> field[i]);
        gtk_container_add(GTK_CONTAINER(GTK_DIALOG(ed -> dialog)->vbox), hbox);   
    }

    gtk_signal_connect(GTK_OBJECT(ed -> dialog), "destroy", GTK_SIGNAL_FUNC(delete_dialog), NULL);
    return ed;
}

GtkToolItem *
toolbar_insert_stock(GtkWidget *toolbar, const gchar *id, GtkSignalFunc func, gpointer data, int toggle, char *label)
{
    GtkToolItem *button;

    if (!toggle) {
        button = gtk_tool_button_new_from_stock(id);
        gtk_signal_connect(GTK_OBJECT(button), "clicked", func, data);
    }
    else {
        button = gtk_toggle_tool_button_new_from_stock(id);
        gtk_signal_connect(GTK_OBJECT(button), "toggled", func, data);
    }

    if (label) {
        gtk_tool_button_set_label(GTK_TOOL_BUTTON(button), label);
    }

    gtk_widget_show(GTK_WIDGET(button));
    gtk_toolbar_insert(GTK_TOOLBAR(toolbar), button, -1);

    return button;
}

void
GetDisplayObjectGeometry(DisplayObject obj, int *x, int *y, int *width, int *height)
{
    gtk_window_get_size(GTK_WINDOW(obj -> toplevel), width, height);
    gtk_window_get_position(GTK_WINDOW(obj -> toplevel), x, y);
}

void
SetDisplayObjectGeometry(DisplayObject obj, int x, int y, int width, int height){
    gtk_window_move(obj -> toplevel, x, y);
    gtk_window_resize(obj -> toplevel, width, height);
}
