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
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#if defined (_WIN32) || defined (__CYGWIN) || defined (_UWIN)
#include <windows.h>
#elif defined (NATIVE_OSX)
#include <mach-o/dyld.h>
#include <libgen.h>
#include <glib/gstdio.h>
#endif
#include "problem.h"
#include "objects.h"
#include "entry.h"
#include "objsheets.h"
#include "combo.h"
#include "list.h"
#include "duplicates.h"
#include "dialog.h"


extern int FileCopy(char *, char*);

extern GtkWidget *toplevel;

ObjectSheet material_sheet;
ObjectSheet connector_sheet;
ObjectSheet anchor_sheet;
ObjectSheet buoy_sheet;

enum {
    MAT_TYPE,
    MAT_MASS,
    MAT_WET,
    MAT_DIAM,
    MAT_LENGTH,
    MAT_CDN, MAT_CDT,
    MAT_AMN, MAT_AMT,
    MAT_EA, MAT_EI, MAT_GJ, 
    MAT_SWL, MAT_YIELD,
    MAT_T, MAT_Te, MAT_Tee,
    MAT_COMMENT,
    NUM_MATERIALS_COLUMNS,
};

static EntryTableEntry materials_entries[] = 
{
    { "type", "choose linear or nonlinear stress-strain relationship", COMBO_ENTRY, 1, {"linear", "nonlinear", NULL} },
    { "mass(kg/m)", "mass per length", NUMBER_ENTRY, 1, { NULL } },
    { "wet(N/m)", "net wet weight per length", NUMBER_ENTRY, 1, { NULL } },
    { "diam(m)", "material diameter", NUMBER_ENTRY, 1,  { NULL } },
    { "length(m)", "length for fixed length objects", NUMBER_ENTRY, 1, { NULL } , 1},
    { "Cdn", "normal drag coefficient", NUMBER_ENTRY, 1, { NULL } },
    { "Cdt", "tangential drag coefficient", NUMBER_ENTRY, 1, { NULL }, 1 },
    { "amn(kg/m)", "normal added mass", NUMBER_ENTRY, 1, { NULL } },
    { "amt(kg/m)", "tangential added mass", NUMBER_ENTRY, 1, { NULL }, 1 },
    { "EA(N)", "axial stiffness", NUMBER_ENTRY, 1, { NULL } },
    { "EI(Nm^2)", "bending stiffness", NUMBER_ENTRY, 1, { NULL } },
    { "GJ(Nm^2)", "torsional stiffness (3D problems only)", NUMBER_ENTRY, 1, { NULL } },
    { "SWL(N)", "safe working load", NUMBER_ENTRY, 1, { NULL } },
    { "yield(N)", "yield or ultimate strength", NUMBER_ENTRY, 1, { NULL }, 1 },
    { "T", "T(e) for nonlinear materials", STRING_ENTRY, 1, { NULL } },
    { "Te", "dT/de for nonlinear materials", STRING_ENTRY, 1, { NULL } },
    { "Tee", "d2T/de2 for nonlinear materials", STRING_ENTRY, 1, { NULL } },
    { "comment", NULL, STRING_ENTRY, 1, { NULL } },
};

enum {
    BUOY_TYPE,
    BUOY_MASS, 
    BUOY_DIAM, 
    BUOY_HEIGHT, 
    BUOY_CDT, 
    BUOY_CDN,
    BUOY_BUOYANCY,
    BUOY_CDW,
    BUOY_SW,
    BUOY_COMMENT,
    BUOY_DIAMETERS,
    NUM_BUOYS_COLUMNS
}; 

static EntryTableEntry buoys_entries[] = 
{
    { "type", "buoy shape or type", 
      COMBO_ENTRY, 1, {"axisymmetric", "sphere", "cylinder", 
                       "capsule", "ship", "platform", NULL} },
    { "mass(kg)", "buoy mass in air", NUMBER_ENTRY, 1, { NULL } },
    { "diam", "basic buoy diameter (use diameters for axisymmetric shapes", 
      NUMBER_ENTRY, 1, { NULL } },
    { "height(m)", "buoy height  (top to bottom length)", 
      NUMBER_ENTRY, 1, { NULL } },
    { "Cdt", "tangential drag coefficient (flow along buoy axis)", 
      NUMBER_ENTRY, 1, { NULL } },
    { "Cdn", "normal or transverse drag coefficient", NUMBER_ENTRY, 1, { NULL } },
    { "buoyancy(N)", "total displacement", NUMBER_ENTRY, 1, { NULL } },
    { "Cdw", "wind drag coefficient", NUMBER_ENTRY, 1, { NULL } },
    { "Sw(m^2)", "projected surface area in wind", NUMBER_ENTRY, 1, { NULL } },
    { "comment", NULL, STRING_ENTRY, 3, { NULL } },
    { "diameters", "shape information for axisymmetric buoys", LIST_ENTRY, 1, {"height", "%.2f", "diam", "%.2f", NULL} },
};

enum {
    CONN_MASS,
    CONN_WET,
    CONN_DIAM,
    CONN_CDT,
    CONN_CDN,
    NUM_CONNECTORS_COLUMNS
};

static EntryTableEntry connectors_entries[] = 
{
    { "mass(kg)", "mass in air", NUMBER_ENTRY, 1, { NULL } },
    { "wet(N)", "net wet weight", NUMBER_ENTRY, 1, { NULL } },
    { "diam(m)", "diameter or effective diameter", NUMBER_ENTRY, 1, { NULL } },
    { "Cdt", "tangential drag coefficient", NUMBER_ENTRY, 1, { NULL } },
    { "Cdn", "normal drag coefficient", NUMBER_ENTRY, 1, { NULL } },
};


enum {
    ANCHOR_MASS,
    ANCHOR_WET,
    ANCHOR_DIAM,
    ANCHOR_CDT,
    ANCHOR_CDN,
    ANCHOR_MU,
    NUM_ANCHORS_COLUMNS
};

static EntryTableEntry anchors_entries[] = 
{
    { "mass(kg)", "mass in air", NUMBER_ENTRY, 1, { NULL } },
    { "wet(N)", "net wet weight", NUMBER_ENTRY, 1, { NULL } },
    { "diam(m)", "diameter or effective diameter", NUMBER_ENTRY, 1, { NULL } },
    { "Cdt", "tangential (along axis) drag coefficient", NUMBER_ENTRY, 1, { NULL } },
    { "Cdn", "normal (transverse) drag coefficient", NUMBER_ENTRY, 1, { NULL } },
    { "mu", "friction coefficient between anchor and bottom", NUMBER_ENTRY, 1, { NULL } },
};

gboolean
whiteSpace(char *str) 
{
    int      has_ws;
    int      has_non_ws;
    char    *ptr;
    int      i;

    has_ws = has_non_ws = 0;
    for (ptr = str ; *ptr ; ptr++) {
        i = isspace(*ptr);
        has_ws += i;
        has_non_ws += !i;
    }

    if (has_ws && has_non_ws)
        return 1;
    else // non whitespace or all whitespaces
        return 0;
}

void
qprintf(FILE *fp, ...)
{
    va_list  ap;
    char    *ptr;
    char     q;
    int      quote_it;

    va_start(ap, fp);    

    while ((ptr = va_arg(ap, char *))) {;
        quote_it = va_arg(ap, int);
        if (quote_it) {
            if (strchr(ptr, '"')) 
                q = '\'';
            else if (strchr(ptr, '\'') || whiteSpace(ptr)) 
                q = '"';
            else 
                q = 0;

            if (q) { fprintf(fp, "%c%s%c", q, ptr, q); }
            else   { fprintf(fp, "%s", ptr); }
        }
        else  {
            fprintf(fp, "%s", ptr); 
        }
    }

    va_end(ap);
}

char *
DatabasePath(char *file)
{
    static char work[256];

#if defined (_WIN32) || defined (__CYGWIN) || defined (_UWIN)
    char *ptr;
    char drive[64], dir[256];

    GetModuleFileName(NULL, work, sizeof(work));
    _splitpath(work, drive, dir, NULL, NULL);

    ptr = dir;
    while (*ptr) {
        if (*ptr == '\\')
            *ptr = '/';
    
        ptr ++;
    }
    sprintf(work, "%s%sdatabase/%s", drive, dir, file);
    // printf("%s\n", work);
#elif defined (NATIVE_OSX)
    static char work2[256];
    uint32_t      sz;
    char         *execpath;
   
    sz = 256;
    _NSGetExecutablePath(work, &sz);
    
    execpath = dirname(work);

    sprintf(work, "%s/Library/Application Support/Cable/", getenv("HOME"));
    if (!g_file_test(work, G_FILE_TEST_IS_DIR)) {
        g_mkdir(work, 0755);
    }
    sprintf(work, "%s/Library/Application Support/Cable/%s", getenv("HOME"), file);
    if (!g_file_test(work, G_FILE_TEST_EXISTS)) {
        sprintf(work2, "%s/../Resources/share/%s", execpath, file);
        if (g_file_test(work2, G_FILE_TEST_EXISTS)) {
            FileCopy(work2, work);
        }            
        else {
        }
    }
    else {
    }
#else
    sprintf(work, "%s/.cable/%s", getenv("HOME"), file);
    // printf("%s\n", work);
#endif
    return g_strdup(work);
}


gboolean
SearchTreeRoot(GtkTreeModel *model, char *search, GtkTreeIter *iter)
{
    gboolean     valid;
    char        *string;

    valid = gtk_tree_model_get_iter_first(model, iter);
    // printf("searching root for %s\n", search);
    while (valid) {
        gtk_tree_model_get(model, iter, 0, &string, -1);
	    // printf("considering %s against %s\n", string, search);
        if (strcmp(string, search) == 0) {
            // free(string);
            return TRUE;
        }
        // free(string);
        valid = gtk_tree_model_iter_next(model, iter);
	    // printf("   not a match, next iter is %s\n", valid ? "valid" : "invalid");
    }
    // printf("%s not found\n", search);

    return FALSE;
    
}

static gboolean
SearchTree(GtkTreeModel *model, char *search, GtkTreeIter *iter1)
{
    GtkTreeIter  iter0;
    gboolean     valid0, valid1;
    char        *string;

    valid0 = gtk_tree_model_get_iter_first(model, &iter0);
    while (valid0) {
        valid1 = gtk_tree_model_iter_children(model, iter1, &iter0); 
        while (valid1) {
            gtk_tree_model_get(model, iter1, 0, &string, -1);
         
            if (strcmp(string, search) == 0) {
                // free(string);
                return TRUE;
            }
            // free(string);
            valid1 = gtk_tree_model_iter_next(model, iter1);
        }
        valid0 = gtk_tree_model_iter_next(model, &iter0);
    }

    return FALSE;
}


char *
BuoyDiameters(Buoy b)
{
    int   j;
    char *buff;
    char  tmp[64];

    buff = (char *) malloc(sizeof(char) * b -> num_diameters * 24);
    buff[0] = 0;

    for (j = 1 ; j <= b -> num_diameters ; j++) {
        sprintf(tmp, "(%.3f, %.3f) ", b -> diameters[j].level, b -> diameters[j].d);
        if (strlen(tmp) + strlen(buff) < b -> num_diameters*24)
            strcat(buff, tmp);
        else
            break;
    }

    return buff;
}

int
ObjectIterator(Item item, void *data)
{
    ObjectSheet *active_sheet = (ObjectSheet *) data;
    CableObject  obj = (CableObject) item;
    gboolean     found;
    GtkTreeIter  iter1, iter2;
    char         buff[20];
    char        *ptr;
    int          i;
    char         name[256];
    int          answer;
    GtkTreePath *path;

    // printf("stuffing %s (%s)\n", obj -> name, obj -> category);

    found = SearchTreeRoot(GTK_TREE_MODEL(active_sheet -> model), 
                           obj -> category, &iter1);
    if (!found) {
	    // printf("adding category %s\n", obj -> category);
        gtk_tree_store_append(active_sheet -> model, &iter1, NULL); 
        gtk_tree_store_set(active_sheet -> model, &iter1, 0, obj -> category, -1);
	    // printf("category added\n");
    }
    
    strncpy(name, obj -> name, 256);
    i = 0;
    while (SearchTree(GTK_TREE_MODEL(active_sheet -> model), name, &iter2)) {
        i ++;
        sprintf(name, "%s_%d", obj -> name, i);
    }
    
    if (i != 0) {
        if (gfDuplicatePref)
            answer = gfDuplicatePref;
        else
            answer = DuplicateDialog(obj -> name, "Object");

        if (answer < 0) {
            answer = gfDuplicatePref = abs(answer);
        }

        if (answer == RENAME) {
	        // printf("found %s, renaming %s\n", obj -> name, name);
            free(obj -> name);
            obj -> name = strdup(name);
        }
        else {
            // printf("ignoring local entry %s\n", obj -> name);
            return 0;
        }
    }    

    // printf("stuffing name %s\n", obj -> name);

    gtk_tree_store_append(active_sheet -> model, &iter2, &iter1);
    gtk_tree_store_set(active_sheet -> model, &iter2, 0, obj -> name, -1);
    for (i = 0 ; i < active_sheet -> num_columns ; i++) {
        ptr = NULL;
        if (active_sheet -> table -> entries[i].type == NUMBER_ENTRY) {
            sprintf(buff, "%g", *((double *) active_sheet -> value_function(obj, i)));
            ptr = buff;
        }
        else if (active_sheet -> table -> entries[i].type == STRING_ENTRY) { 
            ptr = (char *) active_sheet -> value_function(obj, i);
        }
        else if (active_sheet -> table -> entries[i].type == COMBO_ENTRY) {
            ptr = active_sheet -> table -> entries[i].options[*((int *) active_sheet -> value_function(obj, i))];
        }
        else if (active_sheet -> table -> entries[i].type == LIST_ENTRY) {
            // special case so I haven't generalized this functionality yet ...
            if (((Buoy) obj) -> num_diameters) 
                ptr = BuoyDiameters((Buoy) obj);
            else
                ptr = NULL;
        }
        if (ptr) 
            gtk_tree_store_set(active_sheet -> model, &iter2, 
                               i+1, strdup(ptr), -1);
        else
            gtk_tree_store_set(active_sheet -> model, &iter2, 
                               i+1, NULL, -1);
    }

    path = gtk_tree_model_get_path(GTK_TREE_MODEL(active_sheet -> model), &iter1);
    gtk_tree_view_expand_row(GTK_TREE_VIEW(active_sheet -> tree), path, FALSE);
    gtk_tree_path_free(path);

    // printf("done stuffing %s\n", obj -> name);

    return 0;
}

int
ObjectInsert(CableObject obj, ObjectSheet *sheet)
{
    return ObjectIterator((Item) obj, (void *) sheet);
}

static void *
AnchorValue(CableObject obj, int idx)
{
    Anchor a = (Anchor) obj;

    switch(idx) {
    case ANCHOR_MASS: return &a -> m;
    case ANCHOR_WET: return &a -> wet;
    case ANCHOR_DIAM: return &a -> d;
    case ANCHOR_CDT: return &a -> Cdt;
    case ANCHOR_CDN: return &a -> Cdn;
    case ANCHOR_MU: return &a -> mu;
    default:
        printf("unknown anchor index %d\n", idx);
        return NULL;
    }

    return NULL;
}

static void *
ConnectorValue(CableObject obj, int idx)
{
    Connector conn = (Connector) obj;

    switch(idx) {
    case CONN_MASS: return &conn -> m;
    case CONN_WET: return &conn -> wet;
    case CONN_DIAM: return &conn -> d;
    case CONN_CDT: return &conn -> Cdt;
    case CONN_CDN: return &conn -> Cdn;
    default:
        printf("unknown connector index %d\n", idx);
        return NULL;
    }

    return NULL;
}

static void *
BuoyValue(CableObject obj, int idx)
{
    Buoy b = (Buoy) obj;

    switch(idx) {
    case BUOY_TYPE: return &b -> type;
    case BUOY_MASS: return &b -> m; 
    case BUOY_DIAM: return &b -> d; 
    case BUOY_HEIGHT: return &b -> h; 
    case BUOY_CDT: return &b -> Cdt;
    case BUOY_CDN: return &b -> Cdn;
    case BUOY_SW: return &b -> Sw;
    case BUOY_CDW: return &b -> Cdw;
    case BUOY_BUOYANCY: return &b -> buoyancy;
    case BUOY_COMMENT: return b -> comment;
    case BUOY_DIAMETERS:
    default:
        printf("unknown buoy index %d\n", idx);
        return NULL;
    }

    return NULL;
}

static void *
MaterialValue(CableObject obj, int idx)
{
    Material mat = (Material) obj;

    switch(idx) {
    case MAT_TYPE: return &mat -> type;
    case MAT_EA: return &mat -> EA;
    case MAT_EI: return &mat -> EI;
    case MAT_GJ: return &mat -> GJ;
    case MAT_DIAM: return &mat -> d;
    case MAT_WET: return &mat -> wet;
    case MAT_LENGTH: return &mat -> length;
    case MAT_MASS: return &mat -> m;
    case MAT_CDT: return &mat -> Cdt.value;
    case MAT_CDN: return &mat -> Cdn.value;
    case MAT_AMT: return &mat -> amt;
    case MAT_AMN: return &mat -> amn;
    case MAT_YIELD: return &mat -> yield;
    case MAT_SWL: return &mat -> swl;
    case MAT_T:   return mat -> T[0].expr != NULL || mat -> T[0].value != 0 ? 
                         ExpressionString(mat -> T[0]) : NULL;
    case MAT_Te:  return mat -> T[1].expr != NULL || mat -> T[1].value != 0 ? 
                         ExpressionString(mat -> T[1]) : NULL;
    case MAT_Tee: return mat -> T[2].expr != NULL || mat -> T[2].value != 0 ? 
                         ExpressionString(mat -> T[2]) : NULL;
    case MAT_COMMENT: return mat -> comment;
    default:
        printf("unknown material index %d\n", idx);
        return NULL;
    }
}

static void 
PromptDirty(ObjectSheet *sheet)
{
    int          i;
    GtkWidget   *dialog;
    GtkTreeIter  iter;

    if (sheet -> table -> dirty && sheet -> curr_path) { 
        dialog = gtk_message_dialog_new(GTK_WINDOW(toplevel),
                                        GTK_DIALOG_DESTROY_WITH_PARENT,
                                        GTK_MESSAGE_QUESTION,
                                        GTK_BUTTONS_YES_NO,
                                        "Changes to current object not recorded. Record now?");
        if (gtk_dialog_run(GTK_DIALOG (dialog)) == GTK_RESPONSE_YES) {
            gtk_tree_model_get_iter(GTK_TREE_MODEL(sheet -> model), &iter, sheet -> curr_path);
            for (i = 0 ; i < sheet -> num_columns ; i++) {
                if (sheet -> table -> entries[i].type <= STRING_ENTRY) {
                    gtk_tree_store_set(sheet -> model, &iter, i+1, 
                                       gtk_entry_get_text(GTK_ENTRY(sheet -> table -> entries[i].w)), -1);
                }
            }
        }
        gtk_widget_destroy(dialog);
    } 

    return;
}

static void
NewObject(GtkWidget *w, gpointer data)
{
    ObjectSheet *sheet = (ObjectSheet *) data;
    GtkTreeIter child, parent;
    GtkTreeSelection *select;
    GtkTreeModel     *model;

    PromptDirty(sheet);

    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(sheet -> tree));
    if (gtk_tree_selection_get_selected(select, &model, &child)) {
        if (gtk_tree_model_iter_parent(model, &parent, &child) == FALSE) {
            gtk_tree_selection_get_selected(select, &model, &parent);    
        }
        gtk_tree_store_append(GTK_TREE_STORE(model), &child, &parent);
        gtk_tree_store_set(GTK_TREE_STORE(model), &child, 0, "unnamed", -1);
        ClearEntryTable(sheet -> table);

        sheet -> table -> dirty = FALSE;
        if (sheet -> curr_path) {
            gtk_tree_path_free(sheet -> curr_path);
        }
        sheet -> curr_path = gtk_tree_model_get_path(model, &child);
        gtk_tree_selection_select_iter(select, &child);
    }

    return;
}

static void
NewCategory(GtkWidget *w, gpointer data)
{
    ObjectSheet *sheet = (ObjectSheet *) data;
    GtkTreeIter  iter;

    gtk_tree_store_append(sheet -> model, &iter, NULL); 
    gtk_tree_store_set(sheet -> model, &iter, 0, "unnamed", -1);
}

static void
DeleteCableObject(GtkWidget *w, gpointer data)
{
    ObjectSheet *sheet = (ObjectSheet *) data;
    GtkTreeIter iter, parent;
    GtkTreeSelection *select;
    GtkTreeModel     *model;

    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(sheet -> tree));
    if (gtk_tree_selection_get_selected(select, &model, &iter)) {
        if (gtk_tree_model_iter_parent(model, &parent, &iter) == FALSE) {
            if (YesNoDialog("delete entire category?", toplevel) != GTK_RESPONSE_YES)
                return;            
        }
        gtk_tree_store_remove(GTK_TREE_STORE(model), &iter);
    }
}


static void
UpdateObject(GtkWidget *w, gpointer data)
{
    int          i;
    ObjectSheet *sheet = (ObjectSheet *) data;
    GtkTreeIter iter;
    GtkTreeSelection *select;
    GtkTreeModel     *model;

    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(sheet -> tree));
    if (gtk_tree_selection_get_selected(select, &model, &iter)) {
        for (i = 0 ; i < sheet -> num_columns ; i++) {
            if (sheet -> table -> entries[i].type <= STRING_ENTRY) {
                gtk_tree_store_set(GTK_TREE_STORE(model), &iter, i+1, 
                                   gtk_entry_get_text(GTK_ENTRY(sheet -> table -> entries[i].w)), -1);
            }
            else if (sheet -> table -> entries[i].type == LIST_ENTRY) {
                gtk_tree_store_set(GTK_TREE_STORE(model), &iter, i+1, 
                                   ListToString(sheet -> table -> entries[i].w),
                                    -1);
            }
            else if (sheet -> table -> entries[i].type == COMBO_ENTRY) {
                gtk_tree_store_set(GTK_TREE_STORE(model), &iter, i+1, 
                                   gtk_combo_box_get_active_text(GTK_COMBO_BOX(sheet -> table -> entries[i].w)), -1);
            }
        }

        sheet -> table -> dirty = FALSE;
    }

    
    return;
}

static void
FillTable(ObjectSheet *sheet, GtkTreeIter *iter) // CableObject obj)
{
    // char     buff[20];
    char    *ptr;
    int      i;

    for (i = 0 ; i < sheet -> num_columns ; i++) {
        gtk_tree_model_get(GTK_TREE_MODEL(sheet -> model), iter, i + 1, &ptr, -1);
        if (sheet -> table -> entries[i].type <= STRING_ENTRY) {
            SetEntryTableText(sheet -> table, i, ptr);            
        }
        else if (ptr && sheet -> table -> entries[i].type == COMBO_ENTRY) {
            ComboBoxSetText(GTK_COMBO_BOX(sheet -> table -> entries[i].w), ptr, 0);
        }
        else if (sheet -> table -> entries[i].type == LIST_ENTRY) {
	    if (ptr)
            StringToList(ptr, sheet -> table -> entries[i].w);
	    else
		    ClearList(sheet -> table -> entries[i].w);
        }
        g_free(ptr);
    }
}

static void
CopyCableObject(GtkWidget *w, gpointer data)
{
    ObjectSheet *sheet = (ObjectSheet *) data;
    GtkTreeIter iter, iter_new;
    GtkTreeSelection *select;
    GtkTreeModel     *model;
    int               i;
    char             *oldname, *newname, *ptr;

    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(sheet -> tree));
    if (gtk_tree_selection_get_selected(select, &model, &iter)) {

        NewObject(NULL, data);
        fprintf(stderr,"new object created\n");
        // selected object gets changed to newly created object
        gtk_tree_selection_get_selected(select, &model, &iter_new);
        fprintf(stderr,"got new obj ref\n");
        gtk_tree_model_get(GTK_TREE_MODEL(sheet -> model), &iter, 0, &oldname, -1);
        gtk_tree_store_move_after(sheet -> model, &iter_new, &iter);
        fprintf(stderr,"oldname = %s\n", oldname);
        newname = (char *) g_malloc(strlen(oldname) + 6);
        strcpy(newname, oldname);
        strcat(newname, "_copy");
        fprintf(stderr,"newname = %s\n", newname);
        gtk_tree_store_set(sheet -> model, &iter_new, 0, newname, -1);
        for (i = 0 ; i < sheet -> num_columns ; i++) {
            gtk_tree_model_get(GTK_TREE_MODEL(sheet -> model), &iter, i + 1, &ptr, -1);
            gtk_tree_store_set(sheet -> model, &iter_new, i + 1, ptr, -1);
            g_free(ptr);
        }
        FillTable(sheet, &iter_new);
        sheet -> table -> dirty = FALSE;
    }
}

static void
ChangeSelectedObject(GtkTreeSelection *select, gpointer data)
{
    ObjectSheet  *sheet = (ObjectSheet *) data;
    GtkTreeIter   iter;
    GtkTreeModel *model;

    PromptDirty(sheet);

    if (gtk_tree_selection_get_selected(select, &model, &iter)) {
        FillTable(sheet, &iter);

        sheet -> table -> dirty = FALSE;
        if (sheet -> curr_path) {
            gtk_tree_path_free(sheet -> curr_path);
        }
        sheet -> curr_path = gtk_tree_model_get_path(model, &iter);
    }
}

static gboolean
VerifyChildNode(GtkTreeSelection *select, GtkTreeModel *model,
                GtkTreePath *path, gboolean path_selected, gpointer data)
{
    
    if (gtk_tree_path_get_depth(path) < 2) {
        printf("not child node\n");
        return TRUE; // FALSE;
    }

    printf("child node OK\n");
    return TRUE;
}

static void
EditName(GtkCellRendererText *cell, char *path_string, char *new_text, gpointer data)
{
    ObjectSheet *sheet = (ObjectSheet *) data;
    GtkTreeIter  iter;

    printf("editing\n");
    gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(sheet -> model), 
                                        &iter, path_string);
    gtk_tree_store_set(sheet -> model, &iter, 0, new_text, -1);

}

static void
StripUnits(char *src, char *dest)
{
    char *ptr;

    ptr = src;
    while(*ptr) {
        if (*ptr == '(') {
            while (*ptr && *ptr != ')') {
                ptr ++; 
            }   
        }   
        else if (*ptr == ')') {
            ptr ++;
        }
        else {
            *dest = *ptr;
            dest ++;
            ptr ++;
        }
    }

    *dest = 0;
}

static void
WriteObjectEntry(FILE *fp, EntryTableEntry entry, char *ptr)
{

    char    name[256];

    StripUnits(entry.name, name);

    if (ptr && strlen(ptr)) {
        if (entry.type == STRING_ENTRY) 
            qprintf(fp, "        ", 0, name, 1, "=", 0, ptr, 1, "\n", 0, NULL);
        else
            fprintf(fp, "        %s=%s\n", name, ptr);
    }
}

static void
WriteObjectSheet(GtkWidget *w, gpointer data)
{
    ObjectSheet *sheet = (ObjectSheet *) data;
    GtkTreeIter  iter0, iter1;
    int          i;
    gchar       *ptr;
    gchar       *category;
    gboolean     valid1, valid0;
    FILE        *fp;

    fp = fopen(sheet -> db_file, "w");
    if (fp == NULL) {
        fprintf(stderr,"could not open db file %s\n", sheet -> db_file);
        return;
    }
    fprintf(stderr,"opened %s\n", sheet -> db_file);

    fprintf(fp, "%s\n", sheet -> header);

    valid0 = gtk_tree_model_get_iter_first(GTK_TREE_MODEL(sheet -> model), 
                                           &iter0);
    while (valid0) {
        gtk_tree_model_get(GTK_TREE_MODEL(sheet -> model), &iter0, 0, &category, -1);
        valid1 = gtk_tree_model_iter_children(GTK_TREE_MODEL(sheet -> model), 
                                              &iter1, &iter0); 
        while (strcmp(category, "user") && valid1) {
            gtk_tree_model_get(GTK_TREE_MODEL(sheet -> model), &iter1, 0, &ptr, -1);
            qprintf(fp, "    ", 0, ptr, 1, "\n", 0, NULL);
            g_free(ptr);
            qprintf(fp, "        category=", 0, category, 1, "\n", 0, NULL);
            for (i = 0 ; i < sheet -> num_columns ; i++) {
                gtk_tree_model_get(GTK_TREE_MODEL(sheet -> model), &iter1, i + 1, &ptr, -1);
                if (ptr) {
                    WriteObjectEntry(fp, sheet -> table -> entries[i], ptr);
                    g_free(ptr);
                }
            }
         
            valid1 = gtk_tree_model_iter_next(GTK_TREE_MODEL(sheet -> model), &iter1);
        }
        g_free(category);
        valid0 = gtk_tree_model_iter_next(GTK_TREE_MODEL(sheet -> model), &iter0);
    }

    fprintf(fp, "end\n");

    fclose(fp);
}

extern GtkWidget *ImageButton(char *);

GtkWidget *
BuildObjectSheet(ObjectSheet *sheet,
                 EntryTableEntry *entries, int num_entries)
{
    int          i;
    GtkTooltips *tips;
    GtkWidget   *hpaned;
    GtkWidget   *scroller;
    GtkCellRenderer *renderer;
    GtkTreeViewColumn *column;
    GtkTreeSelection *select;
    GtkWidget        *new, *accept, *new_category, *save, *delete, *copy;
    GtkWidget        *hbox_tools, *vbox;
    GType             types[32];

    for (i = 0 ; i < sheet -> num_columns + 1 ; i++) {
        types[i] = G_TYPE_STRING;
    }

    sheet -> model = gtk_tree_store_newv(sheet -> num_columns + 1, types);

    hpaned = gtk_hpaned_new();

    scroller = gtk_scrolled_window_new(NULL, NULL);

    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
                                   GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);

    sheet -> tree = gtk_tree_view_new_with_model(GTK_TREE_MODEL(sheet -> model));
    
    renderer = gtk_cell_renderer_text_new();
    column = gtk_tree_view_column_new_with_attributes("name", renderer, "text", 0, NULL);
    gtk_tree_view_append_column(GTK_TREE_VIEW(sheet -> tree), column);
    g_object_set(renderer, "editable", TRUE, NULL);
    g_signal_connect(renderer, "edited", (GCallback) EditName, sheet);

    for (i = 0 ; i < sheet -> num_columns ; i++) {
        renderer = gtk_cell_renderer_text_new();
        column = gtk_tree_view_column_new_with_attributes(entries[i].name, renderer, "text", i + 1, NULL);
        gtk_tree_view_append_column(GTK_TREE_VIEW(sheet -> tree), column);
    }

    gtk_widget_show(sheet -> tree);

    gtk_container_add(GTK_CONTAINER(scroller), sheet -> tree);
    gtk_paned_pack1(GTK_PANED(hpaned), scroller, TRUE, FALSE); 

    // gtk_tree_view_columns_autosize(GTK_TREE_VIEW(tree));
    gtk_tree_view_set_reorderable(GTK_TREE_VIEW(sheet -> tree), TRUE);
    gtk_tree_view_set_rules_hint(GTK_TREE_VIEW(sheet -> tree), TRUE);

    vbox = gtk_vbox_new(FALSE, 3);
    hbox_tools = gtk_hbox_new(FALSE, 0);
    new = ImageButton(GTK_STOCK_NEW);
    accept = ImageButton(GTK_STOCK_APPLY);
    save = ImageButton(GTK_STOCK_SAVE);
    new_category = ImageButton(GTK_STOCK_OPEN);
    delete = ImageButton(GTK_STOCK_DELETE);
    copy = ImageButton(GTK_STOCK_COPY);

    gtk_box_pack_start(GTK_BOX(hbox_tools), new_category, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox_tools), new, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox_tools), accept, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox_tools), delete, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox_tools), copy, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox_tools), save, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox_tools, FALSE, FALSE, 3);

    gtk_signal_connect(GTK_OBJECT(accept), "clicked", 
                       GTK_SIGNAL_FUNC(UpdateObject), sheet);
    gtk_signal_connect(GTK_OBJECT(delete), "clicked", 
                       GTK_SIGNAL_FUNC(DeleteCableObject), sheet);
    gtk_signal_connect(GTK_OBJECT(new), "clicked", 
                       GTK_SIGNAL_FUNC(NewObject), sheet);
    gtk_signal_connect(GTK_OBJECT(copy), "clicked", 
                       GTK_SIGNAL_FUNC(CopyCableObject), sheet);
    gtk_signal_connect(GTK_OBJECT(new_category), "clicked", 
                       GTK_SIGNAL_FUNC(NewCategory), sheet);
    gtk_signal_connect(GTK_OBJECT(save), "clicked", 
                       GTK_SIGNAL_FUNC(WriteObjectSheet), sheet);

    sheet -> table = MakeEntryTable(entries, NULL, NULL, 
                                    num_entries, ceil(num_entries/3), 3);
    gtk_box_pack_start(GTK_BOX(vbox), sheet -> table -> table, FALSE, FALSE, 3);
    gtk_paned_pack2(GTK_PANED(hpaned), vbox, FALSE, FALSE);

    tips = gtk_tooltips_new();
    gtk_tooltips_set_tip(tips, new, "create new object", NULL);
    gtk_tooltips_set_tip(tips, new_category, "create new category", NULL);
    gtk_tooltips_set_tip(tips, accept, "accept changes to current object", NULL);
    gtk_tooltips_set_tip(tips, delete, "delete current object", NULL);
    gtk_tooltips_set_tip(tips, copy, "create copy of current object", NULL);
    gtk_tooltips_set_tip(tips, save, "save database", NULL);

    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(sheet -> tree));
    gtk_tree_selection_set_mode(select, GTK_SELECTION_SINGLE);
    gtk_tree_selection_set_select_function(select, VerifyChildNode, NULL, NULL);
    g_signal_connect(G_OBJECT (select), "changed",
                     G_CALLBACK(ChangeSelectedObject), sheet);

    sheet -> curr_path = NULL;

    TreeSetAndIterate(sheet -> db_tree, ObjectIterator, (void *) sheet);

    gtk_tree_view_expand_all(GTK_TREE_VIEW(sheet -> tree));

    return hpaned;
}

GtkWidget *
BuildMaterialsSheet(void)
{
    int          i;
    GtkWidget   *hpaned;
    char	*db_path;

    material_sheet.db_tree = TreeCreate(ObjectCompare);

    db_path = DatabasePath("material.db");
    i = ReadDatabaseFile(db_path, material_sheet.db_tree, NULL, NULL, NULL);
    // printf("ReadDatabase returned %d, tree size %d\n", i, TreeSize(material_sheet.db_tree));
    material_sheet.value_function = MaterialValue;
    material_sheet.num_columns = NUM_MATERIALS_COLUMNS;
    material_sheet.db_file = db_path;
    material_sheet.header = "materials";

    hpaned = BuildObjectSheet(&material_sheet, 
                              materials_entries, NUM_MATERIALS_COLUMNS);

    gtk_widget_show_all(hpaned);

    return hpaned;
}

GtkWidget *
BuildBuoysSheet(void)
{
    GtkWidget   *hpaned;
    char        *db_path;

    buoy_sheet.db_tree = TreeCreate(ObjectCompare);

    db_path = DatabasePath("buoy.db");
    ReadDatabaseFile(db_path, NULL, buoy_sheet.db_tree, NULL, NULL);

    buoy_sheet.value_function = BuoyValue;
    buoy_sheet.num_columns = NUM_BUOYS_COLUMNS;
    buoy_sheet.db_file     = db_path;
    buoy_sheet.header = "buoys";

    hpaned = BuildObjectSheet(&buoy_sheet,
                              buoys_entries, NUM_BUOYS_COLUMNS);

    gtk_widget_set_size_request(buoy_sheet.table -> entries[BUOY_DIAMETERS].w, 100, 200);

    gtk_entry_set_width_chars(GTK_ENTRY(buoy_sheet.table -> entries[BUOY_COMMENT].w), 45);

    gtk_widget_show_all(hpaned);

    return hpaned;
}

GtkWidget *
BuildConnectorsSheet(void)
{
    GtkWidget   *hpaned;
    char        *db_path;

    connector_sheet.db_tree = TreeCreate(ObjectCompare);

    db_path = DatabasePath("connect.db");
    ReadDatabaseFile(db_path, NULL, NULL, connector_sheet.db_tree, NULL);

    connector_sheet.value_function = ConnectorValue;
    connector_sheet.num_columns = NUM_CONNECTORS_COLUMNS;
    connector_sheet.db_file     = db_path;
    connector_sheet.header = "connectors";

    hpaned = BuildObjectSheet(&connector_sheet,
                              connectors_entries, NUM_CONNECTORS_COLUMNS);

    gtk_widget_show_all(hpaned);

    return hpaned;
}

GtkWidget *
BuildAnchorsSheet(void)
{
    GtkWidget   *hpaned;
    char 	*db_path;

    anchor_sheet.db_tree = TreeCreate(ObjectCompare);

    db_path = DatabasePath("anchor.db");
    ReadDatabaseFile(db_path, NULL, NULL, NULL, anchor_sheet.db_tree);

    anchor_sheet.value_function = AnchorValue;
    anchor_sheet.num_columns = NUM_ANCHORS_COLUMNS;
    anchor_sheet.db_file     = db_path;
    anchor_sheet.header = "anchors";

    hpaned = BuildObjectSheet(&anchor_sheet,
                              anchors_entries, NUM_ANCHORS_COLUMNS);

    gtk_widget_show_all(hpaned);

    return hpaned;
}


static void
WriteObjects(GSList *list, FILE *fp, ObjectSheet *sheet, int save_all) 
{
    GtkTreeIter  iter, parent;
    gboolean     found;
    int          i;
    gchar       *ptr;

    while(list && list -> data) {
        found = SearchTree(GTK_TREE_MODEL(sheet -> model), list -> data, &iter);
        gtk_tree_model_iter_parent(GTK_TREE_MODEL(sheet -> model), &parent, &iter);
        gtk_tree_model_get(GTK_TREE_MODEL(sheet -> model), &parent, 0, &ptr, -1);

        if (found && (save_all || strcmp(ptr, "user") == 0)) {
            qprintf(fp, "    ", 0, list -> data, 1, "\n", 0, NULL);
            for (i = 0 ; i < sheet -> num_columns ; i++) {
                gtk_tree_model_get(GTK_TREE_MODEL(sheet -> model), &iter, i + 1, &ptr, -1);
                if (ptr && strlen(ptr) > 0) {
                    WriteObjectEntry(fp, sheet -> table -> entries[i], ptr);
                    g_free(ptr);
                }
            }
        }
        else if (found) {
            printf("skipping %s - already in the global database\n", (char *) list -> data);
        }
        else {
            printf("%s not found\n", (char *) list -> data);
        }
        list = list -> next;
    }
}

void
WriteBuoys(GSList *list, FILE *fp, int save_all_objects)
{
    fprintf(fp, "Buoys\n");
    WriteObjects(list, fp, &buoy_sheet, save_all_objects); 
}
void
WriteMaterials(GSList *list, FILE *fp, int save_all_objects)
{
    fprintf(fp, "Materials\n");
    WriteObjects(list, fp, &material_sheet, save_all_objects);
}
void
WriteAnchors(GSList *list, FILE *fp, int save_all_objects)
{
    fprintf(fp, "Anchors\n");
    WriteObjects(list, fp, &anchor_sheet, save_all_objects);
}
void
WriteConnectors(GSList *list, FILE *fp, int save_all_objects)
{
    fprintf(fp, "Connectors\n");
    WriteObjects(list, fp, &connector_sheet, save_all_objects);
}

void
WriteDatabases(void)
{
    WriteObjectSheet(NULL, (gpointer *) &material_sheet);
    WriteObjectSheet(NULL, (gpointer *) &anchor_sheet);
    WriteObjectSheet(NULL, (gpointer *) &connector_sheet);
    WriteObjectSheet(NULL, (gpointer *) &buoy_sheet);
}
