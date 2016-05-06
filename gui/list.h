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

#ifndef _LIST_H
#define _LIST_H

typedef struct {
    int           number;
    char         *fmt;
    int           type;
    GtkListStore *parent;
} ListColumn;

extern GtkWidget *MakeList(char **columns, GtkWidget **, gboolean);
extern void AddListRow(GtkWidget *, gpointer);
extern char *ListToString(GtkWidget *w);
extern int *ListToArray(GtkWidget *w, int *, int, int);
extern void StringToList(char *str, GtkWidget *w);
extern void ClearList(GtkWidget *);



#endif // _LIST_H
