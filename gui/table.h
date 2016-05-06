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

#ifndef _TABLE_H
#define _TABLE_H

extern void gtk_table_move_row(GtkTable *, int, int);
extern void gtk_table_insert_row(GtkTable *, int);
extern void gtk_table_delete_row(GtkTable *, int);
extern void gtk_table_rotate_rows(GtkTable *, int, int);
extern int gtk_table_get_child_row(GtkTable *, GtkWidget *);
extern GtkWidget *gtk_table_get_child_at(GtkTable *, int, int);

#endif // _TABLE_H
