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
#include <string.h>
#include <ctype.h>
#include "list.h"

char *
ListToString(GtkWidget *w)
{
    GtkTreeModel    *model;
    GtkTreeView     *tree = (GtkTreeView *) w;
    GtkTreeIter      iter;
    gboolean         status;
    float            x1, x2;
    char             buff[256];
    char             tmp[256];

    buff[0] = 0;

    model = gtk_tree_view_get_model(tree);    
    status = gtk_tree_model_get_iter_first(model, &iter);

    while (status) {
        gtk_tree_model_get(model, &iter, 0, &x1, 1, &x2, -1);
        sprintf(tmp, "(%.2f, %.2f)", x1, x2);
        strcat(buff, tmp);
        status = gtk_tree_model_iter_next(model, &iter);
    }
    
    if (buff[0])
        return strdup(buff);
    else
        return NULL;
}

int *
ListToArray(GtkWidget *w, int *count, int min, int max)
{
    GtkTreeModel    *model;
    GtkTreeView     *tree = (GtkTreeView *) w;
    GtkTreeIter      iter;
    gboolean         status;
    int              x1;
    int              n;
    int             *data;

    n = 0;

    model = gtk_tree_view_get_model(tree);    

    status = gtk_tree_model_get_iter_first(model, &iter);
    while (status) {
        gtk_tree_model_get(model, &iter, 0, &x1, -1);
        if (x1 >= min && x1 <= max) {        
            n ++;
        }

        status = gtk_tree_model_iter_next(model, &iter);
    }
  
    if (n == 0) {
        *count = 0;
        return NULL;
    }

    data = (int *) malloc(sizeof(int) * n);

    n = 0;

    status = gtk_tree_model_get_iter_first(model, &iter);
    while (status) {
        gtk_tree_model_get(model, &iter, 0, &x1, -1);
        if (x1 >= min && x1 <= max) {        
            data[n] = x1;
            n++ ;
        }
        status = gtk_tree_model_iter_next(model, &iter);
    }
 
    *count = n;
    return data;
}

void 
collapse(char *s)
{
    int i,j,n;

    n = strlen(s);
    i = 0;
    for (j = 0 ; j < n ; j++) {
        if (isspace(s[j]))
            continue;

        s[i++] = s[j];
    }
    s[i] = 0;
}

void
StringToList(char *str, GtkWidget *w)
{
    GtkTreeModel     *model;
    GtkTreeView      *tree = (GtkTreeView *) w;
    GtkTreeIter       iter;
    GtkTreeSelection *select;
    char             *ptr1, *ptr2;
    float             x1, x2;

    select = gtk_tree_view_get_selection(tree);
    model = gtk_tree_view_get_model(tree);    
    gtk_list_store_clear(GTK_LIST_STORE(model));

    ptr1 = str;
    collapse(ptr1);
    while ((ptr2 = strsep(&ptr1, "(")) && ptr2 != '\0') {
        if (!ptr1) {
            break;
        }

        if (sscanf(ptr1, "%f,%f)", &x1, &x2) == 2) {
            AddListRow(NULL, (gpointer) w);
            gtk_tree_selection_get_selected(select, &model, &iter);
            gtk_list_store_set((GtkListStore *) model, &iter, 0, x1, 1, x2, -1);
        }
        else
            break;
    }
        
}

void
AddListRow(GtkWidget *w, gpointer data)
{
    GtkTreeView      *tree = (GtkTreeView *) data;
    GtkTreeSelection *select;
    GtkTreeIter       new;
    GtkTreeModel     *model;
    
    select = gtk_tree_view_get_selection(tree);
    model = gtk_tree_view_get_model(tree);
    gtk_list_store_prepend(GTK_LIST_STORE(model), &new);
    gtk_tree_selection_select_iter(select, &new);
//    gtk_list_store_set(GTK_LIST_STORE(model), &new, 0, 0.0, 1, 0.0, -1);
}

static void
DelRow(GtkWidget *w, gpointer data)
{
    GtkTreeView      *tree = (GtkTreeView *) data;
    GtkTreeSelection *select;
    GtkTreeIter       iter;
    GtkTreeModel     *model;

    select = gtk_tree_view_get_selection(tree);
    if (gtk_tree_selection_get_selected(select, &model, &iter)) {
        if (gtk_list_store_remove(GTK_LIST_STORE(model), &iter)) {
            gtk_tree_selection_select_iter(select, &iter);
        }    
    }
}

void
ClearList(GtkWidget *w)
{
    GtkTreeView  *tree = (GtkTreeView *) w;
    GtkTreeModel *model;

    model = gtk_tree_view_get_model(tree);    
    gtk_list_store_clear(GTK_LIST_STORE(model));
}

static void
EditStart(GtkCellRendererText *cell, GtkCellEditable *editable, 
            gchar *path, gpointer data)
{
    fprintf(stderr,"edit started\n");
    gtk_drag_source_unset(GTK_WIDGET(editable));
}

static void
EditStop(GtkCellRendererText *cell, gpointer data)
{
    fprintf(stderr,"edit stopped\n");
}

static void
EditCol(GtkCellRendererText *cell, char *path_string, 
        char *new_text, gpointer data)
{
    ListColumn *col = (ListColumn *) data;
    GtkTreeIter  iter;
    GtkTreePath *path;

    gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(col -> parent), &iter, path_string);
    if (col -> type == G_TYPE_FLOAT) {
        gtk_list_store_set(col -> parent, &iter, col -> number, atof(new_text), -1);
    }
    else if (col -> type == G_TYPE_INT) {
        gtk_list_store_set(col -> parent, &iter, col -> number, atoi(new_text), -1);
    }

    path = gtk_tree_path_new_from_string(path_string), 
    gtk_tree_model_row_changed(GTK_TREE_MODEL(col -> parent), path, &iter);
    gtk_tree_path_free(path);
}

static void
RenderCol(GtkTreeViewColumn *treecol, GtkCellRenderer *renderer,
          GtkTreeModel *model, GtkTreeIter *iter, gpointer data)
{
    ListColumn *col = (ListColumn *) data;
    gfloat  val;
    gint    ival;
    gchar   buf[20];

    if (col -> type == G_TYPE_FLOAT) {
        gtk_tree_model_get(model, iter, col -> number, &val, -1);
        g_snprintf(buf, sizeof(buf), col -> fmt, val);
    }
    else if (col -> type == G_TYPE_INT) {
        gtk_tree_model_get(model, iter, col -> number, &ival, -1);
        g_snprintf(buf, sizeof(buf), col -> fmt, ival);
    }
    
    g_object_set(renderer, "text", buf, NULL);
}

GtkWidget *
MakeList(char **columns, GtkWidget **container, gboolean user_add) 
{
    GtkWidget   *vbox, *button_box;
    GtkWidget   *scroller;
    GtkWidget   *tree;
    GtkCellRenderer *renderer;
    GtkTreeViewColumn *column;
    GtkWidget     *add, *del;
    GtkListStore  *list;
    char         **ptr;
    int            numcol;
    GType         *types;
    ListColumn    *col;
    int            i;

    vbox = gtk_vbox_new(FALSE, 3);

    scroller = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
                                   GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);

    gtk_box_pack_start(GTK_BOX(vbox), scroller, TRUE, TRUE, 3);
    button_box = gtk_hbox_new(FALSE, 3);
    gtk_box_pack_start(GTK_BOX(vbox), button_box, FALSE, FALSE, 3);

    if (user_add) {
        add = gtk_button_new_with_label("+");
        del = gtk_button_new_with_label(" - ");

        gtk_box_pack_start(GTK_BOX(button_box), add, FALSE, FALSE, 3);
        gtk_box_pack_start(GTK_BOX(button_box), del, FALSE, FALSE, 3);
    }
    else {
        add = del = NULL;
    }

    numcol = 0;
    ptr = columns;
    while(*ptr != NULL) {
        numcol ++;
        ptr ++;
    }
    numcol = numcol / 2;
    
    types = (GType *) malloc(sizeof(GType) * numcol);
    for (i = 0 ; i < numcol ; i++) {
        if (strchr(columns[2*i + 1], 's'))
            types[i] = G_TYPE_STRING;
        else if (strchr(columns[2*i + 1], 'f'))
            types[i] = G_TYPE_FLOAT;
        else if (strchr(columns[2*i + 1], 'd'))
            types[i] = G_TYPE_INT;
    }

    list = gtk_list_store_newv(numcol, types);
    gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(list), 
                                         0, GTK_SORT_ASCENDING);
    tree = gtk_tree_view_new_with_model(GTK_TREE_MODEL(list));

    gtk_tree_view_set_reorderable(GTK_TREE_VIEW(tree), FALSE);

    col = (ListColumn *) malloc(sizeof(ListColumn) * numcol);
    for (i = 0 ; i < numcol ; i++) {
        col[i].fmt = columns[2*i + 1];
        col[i].number = i;
        col[i].parent = list;
        col[i].type = types[i];

        renderer = gtk_cell_renderer_text_new();
        column = gtk_tree_view_column_new_with_attributes(columns[2*i], renderer, "text", i, NULL);
        gtk_tree_view_append_column(GTK_TREE_VIEW(tree), column);
        g_object_set(renderer, "editable", TRUE, NULL);
        
        g_signal_connect(renderer, "edited", (GCallback) EditCol, (gpointer) &col[i]);
        g_signal_connect(renderer, "editing-started", (GCallback) EditStart, (gpointer) &col[i]);
        g_signal_connect(renderer, "editing-canceled", (GCallback) EditStop, (gpointer) &col[i]);
        if (types[i] == G_TYPE_FLOAT || types[i] == G_TYPE_INT) {
            gtk_tree_view_column_set_cell_data_func(column, renderer, RenderCol, (gpointer) &col[i], NULL);
        }
    }

    gtk_container_add(GTK_CONTAINER(scroller), tree);

    *container = vbox;

    if (user_add) {
        gtk_signal_connect(GTK_OBJECT(add), "clicked", 
                           GTK_SIGNAL_FUNC(AddListRow), tree);
        gtk_signal_connect(GTK_OBJECT(del), "clicked", 
                           GTK_SIGNAL_FUNC(DelRow), tree);
    }

    return tree;
}
