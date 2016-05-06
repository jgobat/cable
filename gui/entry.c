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

#include <stdlib.h>
#include <gtk/gtk.h>
#include "entry.h"
#include "list.h"

static void
MarkDirty(GtkEditable *editable, gpointer table)
{
    ((EntryTable *) table) -> dirty = TRUE;
}

void
ClearEntryTable(EntryTable *table)
{
    int i;

    for (i = 0 ; i < table -> num_entries ; i++) {
        if (table -> entries[i].type <= STRING_ENTRY) {
            gtk_entry_set_text(GTK_ENTRY(table -> entries[i].w), "");
        }
        else if (table -> entries[i].type == COMBO_ENTRY) {
            gtk_combo_box_set_active(GTK_COMBO_BOX(table -> entries[i].w), 0);
        }
        else if (table -> entries[i].type == LIST_ENTRY) {
	    ClearList(table -> entries[i].w);
        }
    }
       
    return; 
}

static GtkWidget *
MakeWidget(EntryTable *table, int i, GtkWidget **container)
{
    GtkWidget   *w;
    int          j;

    if (table -> entries[i].type <= STRING_ENTRY) {
        w = gtk_entry_new();
        gtk_entry_set_width_chars(GTK_ENTRY(w), 10);
        gtk_signal_connect(GTK_OBJECT(w), "changed",
                           GTK_SIGNAL_FUNC(MarkDirty), table);
        
        *container = w;
    }    
    else if (table -> entries[i].type == COMBO_ENTRY) {
        w = gtk_combo_box_new_text();
        for (j = 0 ; table -> entries[i].options[j] != NULL ; j++) {
            gtk_combo_box_append_text(GTK_COMBO_BOX(w), 
                                      table -> entries[i].options[j]);
        }
        *container = w;
    }
    else if (table -> entries[i].type == LIST_ENTRY) {
        w  = MakeList(table -> entries[i].options, container, TRUE);
    }
    else {
        w  = NULL;
    }

    return w;
}

EntryTable *
MakeEntryTable(EntryTableEntry *entries, char **row_labels, char **col_labels, 
               int num_entries, int num_rows, int num_cols)
{
    EntryTable  *table;
    int          i, j;
    GtkWidget   *label;
    int          row;
    int          col;
    int		 span;
    GtkWidget   *w;

    table = (EntryTable *) malloc(sizeof(EntryTable));

    table -> dirty = FALSE;
    table -> num_entries = num_entries;
    table -> entries = entries;

    table -> tips = gtk_tooltips_new();
    gtk_tooltips_enable(table -> tips);

    if (!row_labels || !col_labels) {
        table -> table = gtk_table_new(num_rows, 2*num_cols, FALSE);

        row = col = 0;
        for (i = 0 ; i < num_entries ; i++)  {
            label = gtk_label_new(entries[i].name);
            gtk_label_set_justify(GTK_LABEL(label), GTK_JUSTIFY_RIGHT);

            table -> entries[i].w = MakeWidget(table, i, &w);

            if (table -> entries[i].w) {
		        if (table -> entries[i].tip) {
                    gtk_tooltips_set_tip(table -> tips, table -> entries[i].w, 
                                         table -> entries[i].tip, NULL);
                    gtk_tooltips_set_tip(table -> tips, label, 
                                        table -> entries[i].tip, NULL);
                }
                if (table -> entries[i].span) {
	                span = table -> entries[i].span;
	            }
                else {
                    span = 1;
                }

                gtk_table_attach(GTK_TABLE(table -> table), label, 
                                 2*col, 2*col + 1, row, row + 1, 0, 0, 0, 0);
                gtk_table_attach(GTK_TABLE(table -> table), w,
                                 2*col + 1, 2*col + 1 + span, row, row + 1, 0, 0, 0, 0);
                col = col + span;
                if (col >= num_cols || entries[i].rowbreak) {
                   row = row + 1;
                   col = 0;
                }
            }
        }
    }
    else {
        table -> table = gtk_table_new(num_rows + 1, num_cols + 1, FALSE);
        for (i = 0 ; i < num_cols ; i++) {
            label = gtk_label_new(col_labels[i]);
            gtk_label_set_justify(GTK_LABEL(label), GTK_JUSTIFY_CENTER);
            gtk_table_attach(GTK_TABLE(table -> table), label, 
                             i+1, i+2, 0, 1, 0, 0, 0, 0); 
        }
        for (i = 0 ; i < num_rows ; i++) {
            label = gtk_label_new(row_labels[i]);
            gtk_label_set_justify(GTK_LABEL(label), GTK_JUSTIFY_RIGHT);
            gtk_table_attach(GTK_TABLE(table -> table), label, 
                             0, 1, i+1, i+2, 0, 0, 0, 0); 
            for (j = 0 ; j < num_cols ; j++) {
                table -> entries[i*num_cols + j].w = 
                                MakeWidget(table, i*num_cols + j, &w);

		        if (table -> entries[i*num_cols + j].w) {
                    if (table -> entries[i*num_cols + j].tip) {
			            gtk_tooltips_set_tip(table -> tips, 
                                   table -> entries[i*num_cols + j].w,
                                   table -> entries[i*num_cols + j].tip, NULL);
                    }
                    gtk_table_attach(GTK_TABLE(table -> table), w,
                                     j+1, j+2, i+1, i+2, 0, 0, 0, 0);
                }
            }
        }
                
    }

    gtk_table_set_row_spacings(GTK_TABLE(table -> table), 3);
    gtk_widget_show_all(table -> table);
    return table;
}



void
SetEntryTableText(EntryTable *table, int entry, char *text)
{
    if (text) 
        gtk_entry_set_text(GTK_ENTRY(table -> entries[entry].w), text);
    else
        gtk_entry_set_text(GTK_ENTRY(table -> entries[entry].w), "");
}
