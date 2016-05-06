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

#ifndef _ENTRY_H
#define _ENTRY_H

#define NUMBER_ENTRY 0
#define STRING_ENTRY 1
#define COMBO_ENTRY 2
#define LIST_ENTRY 3
#define EMPTY_ENTRY 4

typedef struct {
    char       *name;
    char       *tip;
    int         type;
    int         span;
    char       *options[10];
    int         rowbreak;
    GtkWidget  *w;
} EntryTableEntry;

typedef struct {
    GtkWidget       *table;
    EntryTableEntry *entries;
    int              num_entries;
    int              dirty;
    GtkTooltips     *tips;
} EntryTable;

extern EntryTable *MakeEntryTable(EntryTableEntry *, char **, char **, int, int, int);
extern void ClearEntryTable(EntryTable *);
extern void SetEntryTableText(EntryTable *table, int entry, char *text);


#endif // _ENTRY_H
