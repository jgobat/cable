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

#ifndef _OBJSHEETS_H
#define _OBJSHEETS_H

typedef struct {
    Tree           db_tree;
    GtkTreeStore  *model;
    GtkWidget     *tree;
    EntryTable    *table;
    void          *(*value_function)(CableObject, int);
    int            num_columns;
    GtkTreePath   *curr_path;
    char          *db_file;
    char          *header;
} ObjectSheet;

GtkWidget *BuildMaterialsSheet(void);
GtkWidget *BuildConnectorsSheet(void);
GtkWidget *BuildAnchorsSheet(void);
GtkWidget *BuildBuoysSheet(void);

void WriteBuoys(GSList *, FILE *, int);
void WriteAnchors(GSList *, FILE *, int);
void WriteMaterials(GSList *, FILE *, int);
void WriteConnectors(GSList *, FILE *, int);

extern void WriteDatabases(void);

int ObjectInsert(CableObject, ObjectSheet *);
int ObjectIterator(Item, void *);

char *DatabasePath(char *);
void qprintf(FILE *, ...);

#endif // _OBJSHEETS_H
