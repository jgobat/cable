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

#include <stdio.h>
#include <gtk/gtk.h>

void
SetNumericText(GtkWidget *w, char *fmt, double value, gboolean blank_if_zero)
{
    char    buff[20];

    if (value || !blank_if_zero) {
        sprintf(buff, fmt, value);
        gtk_entry_set_text(GTK_ENTRY(w), buff);
    }
}

void
WriteEntry(GtkWidget *w, char *label, FILE *fp)
{
    const gchar    *ptr;

    ptr = gtk_entry_get_text(GTK_ENTRY(w));
    if (ptr && ptr[0]) {
        fprintf(fp, "    %s = %s\n", label, ptr);
    }
}

