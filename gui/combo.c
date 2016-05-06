/*
    This file is part of WHOI Cable, a program for the static and dynamic
    analysis of oceanographic cable structures.

    Copyright (C) 2008 by Jason Gobat

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

#include <string.h>
#include <stdlib.h>
#include <gtk/gtk.h>

char *
ComboBoxGetText(GtkComboBox *box, int n)
{
    GtkTreeIter iter;
    char        *name;
    GtkListStore *list;

    list = (GtkListStore *) gtk_combo_box_get_model(box);
    if ((n == -1 && gtk_combo_box_get_active_iter(box, &iter)) 
        || (n > -1 && gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(list), &iter, NULL, n))) {
        gtk_tree_model_get(GTK_TREE_MODEL(list), &iter, 0, &name, -1);
    }
    else {
        name = NULL;
    }

    return name;
}

static gboolean
compare(GtkTreeModel *model, GtkTreeIter *iter, char *text)
{
    char         *name;

    gtk_tree_model_get(model, iter, 0, &name, -1);
    if (strcmp(name, text) == 0) {
        free(name);
        return TRUE;
    }

    free(name);
    return FALSE;
}

gboolean
ComboBoxSetText(GtkComboBox *box, char *text, int search_depth)
{
    GtkTreeIter   iter0, iter1;
    GtkTreeModel *model;
    gboolean      found0, found1;

    if (text == NULL) {
        gtk_combo_box_set_active(box, 0);
    }

    model = gtk_combo_box_get_model(GTK_COMBO_BOX(box));
    found0 = gtk_tree_model_get_iter_first(model, &iter0);
    while(found0) {
        if (search_depth == 1) {
            found1 = gtk_tree_model_iter_children(model, &iter1, &iter0);
            while(found1) {
                if (compare(model, &iter1, text)) {
                    gtk_combo_box_set_active_iter(GTK_COMBO_BOX(box), &iter1);
                    return TRUE;
                }

                found1 = gtk_tree_model_iter_next(model, &iter1);
            }
        }
        else {
            if (compare(model, &iter0, text)) {
                gtk_combo_box_set_active_iter(GTK_COMBO_BOX(box), &iter0);
                return TRUE;
            }
        }
        found0 = gtk_tree_model_iter_next(model, &iter0);
    }

    return FALSE;
}

void
ComboBoxClear(GtkWidget *widget, gpointer data)
{
    gtk_combo_box_set_active((GtkComboBox *) data, -1);
}
