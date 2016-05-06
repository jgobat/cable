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
#include <stdlib.h>
#include <string.h>
#include <gtk/gtk.h>
#include <glib.h>
#include <ctype.h>
#include "problem.h"
#include "text.h"
#include "entry.h"
#include "combo.h"
#include "list.h"
#include "duplicates.h"
#include "objsheets.h"

// the first of these needs to match the first in the table
// definition below

enum {
    ENV_DEPTH,
    ENV_DENSITY,
    ENV_GRAVITY,
    ENV_STIFFNESS,
    ENV_DAMPING,
    ENV_PROFILE,
    ENV_XCURRENT,
    ENV_XCURR_PROFILE,
    ENV_YCURRENT,
    ENV_YCURR_PROFILE,
    ENV_CURRENT_SCALE,
    ENV_XWIND,
    ENV_YWIND,
    ENV_AMP,
    ENV_PERIOD,
    ENV_RANDOM,
    ENV_REGULAR,
    ENV_FORCING,
    NUM_ENVIRONMENT,
};

static GtkListStore *profiles_list;

GtkWidget *envparam[NUM_ENVIRONMENT];
static EntryTable *envtable;

static EntryTableEntry env_entries[] = 
{
    { "water depth",      "water depth (m)", NUMBER_ENTRY, 1, { NULL } },
    { "water density",    "water density (kg/m^3)", NUMBER_ENTRY, 1, { NULL } }, 
    { "gravity",          "acceleration of gravity (typically 9.81) (m/s^2)", NUMBER_ENTRY, 1, { NULL } }, 
    { "bottom stiffness", "spring stiffness of bottom(N/m/m)", NUMBER_ENTRY, 1, { NULL } }, 
    { "bottom damping",   "damping ratio of bottom", NUMBER_ENTRY, 1, { NULL } }, 
};

extern gboolean ComputeTotalsForDepth(GtkWidget *, GdkEventFocus *, gpointer);
extern GtkWidget *ImageButton(char *);

static void
EditName(GtkCellRendererText *cell, char *path_string, char *new_text, gpointer data)
{
    GtkTreeModel *model = (GtkTreeModel *) data;
    GtkTreeIter   iter;

    gtk_tree_model_get_iter_from_string(model, &iter, path_string);
    gtk_list_store_set(GTK_LIST_STORE(model), &iter, 0, new_text, -1);
}

gboolean SearchTreeRoot(GtkTreeModel *, char *, GtkTreeIter *);

static char *
UnusedProfileName(GtkTreeModel *model)
{
    char        buff[16];
    GtkTreeIter iter;
    int         i;

    for (i = 0 ; i < 2000 ; i++) {
        sprintf(buff, "user%04d", i);
        if (SearchTreeRoot(model, buff, &iter) == FALSE) {
            break;
        }
    }

    return g_strdup(buff);
}

static void 
NewProfile(GtkWidget *w, gpointer data)
{
    GtkTreeIter  iter;
    char        *buff;

    buff = UnusedProfileName(GTK_TREE_MODEL(data));
    gtk_list_store_append(GTK_LIST_STORE(data), &iter);
    gtk_list_store_set(GTK_LIST_STORE(data), &iter, 0, buff, 1, "0.0", -1);
}

static void
AddProfile(GtkListStore *list, char *name, char *expr, GtkTreeIter *iter)
{
    char    *ptr;

    ptr = expr;
    while (*ptr) {
        if (*ptr == '{') {
            while (*ptr && *ptr != '}') {
                ptr ++;
            }
        }
        else if (*ptr == '}') {
            ptr ++;
        }
        else if (isspace(*ptr)) {
            ptr ++;
        }
        else
            break;
    }
        
    gtk_list_store_append(list, iter);
    gtk_list_store_set(list, iter, 0, name, 1, ptr, -1);

    return;
}

static void
UpdateProfile(GtkWidget *w, gpointer data)
{
    GtkTreeView      *tree = (GtkTreeView *) data;
    GtkTreeIter       iter;
    GtkTreeSelection *select;
    GtkTreeModel     *model;

    select = gtk_tree_view_get_selection(tree);
    if (gtk_tree_selection_get_selected(select, &model, &iter)) {
        gtk_list_store_set(GTK_LIST_STORE(model), &iter, 1,
                           ListToString(envparam[ENV_PROFILE]), -1); 
    }
}

static void
DeleteProfile(GtkWidget *w, gpointer data)
{
    GtkTreeView      *tree = (GtkTreeView *) data;
    GtkTreeIter       iter;
    GtkTreeSelection *select;
    GtkTreeModel     *model;

    select = gtk_tree_view_get_selection(tree);
    if (gtk_tree_selection_get_selected(select, &model, &iter)) {
        gtk_list_store_remove(GTK_LIST_STORE(model), &iter);
    }

    ClearList(envparam[ENV_PROFILE]);
}

static void
SaveProfileDatabase(GtkWidget *w, gpointer data)
{
    char    *fname;
    FILE    *fp;
    char    *expr, *name;
    GtkTreeIter iter;
    gboolean valid;

    fname = DatabasePath("current.db");
    fp = fopen(fname, "wb");
    g_free(fname);

    if (fp == NULL) {
        return;
    }

    valid = gtk_tree_model_get_iter_first(GTK_TREE_MODEL(profiles_list), &iter);
    while(valid) {
        gtk_tree_model_get(GTK_TREE_MODEL(profiles_list), &iter, 0, &name, 1, &expr, -1);
        fprintf(fp, "{%s} %s\n", name, expr);
        valid = gtk_tree_model_iter_next(GTK_TREE_MODEL(profiles_list), &iter);
    }

    fclose(fp);
}

static int
ReadCurrentDatabase(char *fname, GtkListStore *list)
{
    FILE    *fp;
    char     buffer[2112], expr[2048], name[64];
    GtkTreeIter iter;

    fp = fopen(fname, "rb");
    if (fp == NULL) {
        return 1;
    }

    gtk_list_store_clear(list);

    while(!feof(fp)) {
        if (fgets(buffer, 2112, fp) == NULL) {
            break;
        }
       
        sscanf(buffer, "{%64[^}]} %2048[^\r\n]", name, expr); 
        AddProfile(list, name, expr, &iter);
    }

    fclose(fp);
    return 0;
}

static void
ChangeSelectedProfile(GtkTreeSelection *select, gpointer data)
{
    GtkTreeIter   iter;
    GtkTreeModel *model;
    char         *ptr;
    double        x;
    char         *xptr;

    if (gtk_tree_selection_get_selected(select, &model, &iter)) {
        gtk_tree_model_get(model, &iter, 1, &ptr, -1);        
        x = strtod(ptr, &xptr);
        if (xptr == ptr) { // list(?) as indicated(?) by leading paren
           StringToList(ptr, GTK_WIDGET(data));
        }
        else {
        }
    }
}

static int clearing = 0;

void
ClearProfileBox(GtkWidget *w, gpointer data)
{
    if (!clearing) {
        clearing = 1;
        gtk_combo_box_set_active(GTK_COMBO_BOX(data), -1);
        clearing = 0;
    }
}

void 
ClearCurrentEntry(GtkWidget *w, gpointer data)
{
    GtkTreeIter  iter;
    if (!clearing && gtk_combo_box_get_active_iter(GTK_COMBO_BOX(w), &iter)) {
       clearing = 1;
       gtk_entry_set_text(GTK_ENTRY(data), "");
       clearing = 0;
    }
}

GtkWidget *
BuildEnvironment(void)
{ 
    int          i;
    GtkWidget   *frame;
    GtkWidget        *vbox, *hbox_tools;
    GtkWidget        *new, *accept, *save, *delete;
    GtkWidget        *basic_hbox, *adv_hbox, *hbox;
    GtkWidget        *wave_vbox;
    GtkWidget        *current_vbox;
    EntryTable       *table;
    GtkWidget        *envbook;
    GtkWidget        *current_table;
    GtkWidget        *scroller, *tree;
    GtkTreeSelection *select;
    GtkTooltips      *tips;
    GtkWidget        *wave_random, *wave_regular;
    GtkWidget        *forcing_type;
    GtkWidget        *amplitude, *period;
    char             *fname;
    static char      *current_columns[] = {"depth", "%.0f", "current", "%.2f", NULL};
    GtkCellRenderer *renderer;
    GtkTreeViewColumn *column;

    frame = gtk_frame_new("Environment");
    envbook = gtk_notebook_new();
    gtk_container_add(GTK_CONTAINER(frame), envbook);

    gtk_notebook_set_tab_pos(GTK_NOTEBOOK(envbook), GTK_POS_TOP);

    // tab layout: basic, current profiles, ...

    basic_hbox = gtk_hbox_new(FALSE, 3);
    gtk_notebook_append_page(GTK_NOTEBOOK(envbook), basic_hbox, gtk_label_new("basic"));

    adv_hbox = gtk_hbox_new(FALSE, 3);
    gtk_notebook_append_page(GTK_NOTEBOOK(envbook), adv_hbox, gtk_label_new("current profiles"));

    // current profiles tab 

    profiles_list = gtk_list_store_new(2, G_TYPE_STRING, G_TYPE_STRING);
    scroller = gtk_scrolled_window_new(NULL, NULL);    
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
                                   GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
    tree = gtk_tree_view_new_with_model(GTK_TREE_MODEL(profiles_list));

    renderer = gtk_cell_renderer_text_new();
    column = gtk_tree_view_column_new_with_attributes("name", renderer, "text", 0, NULL);
    gtk_tree_view_append_column(GTK_TREE_VIEW(tree), column);
    g_object_set(renderer, "editable", TRUE, NULL);
    g_signal_connect(renderer, "edited", (GCallback) EditName, profiles_list);

    renderer = gtk_cell_renderer_text_new();
    column = gtk_tree_view_column_new_with_attributes("expression", renderer, "text", 1, NULL);
    gtk_tree_view_append_column(GTK_TREE_VIEW(tree), column);

    gtk_widget_set_size_request(tree, 300, 180);

    gtk_container_add(GTK_CONTAINER(scroller), tree);
    gtk_box_pack_start(GTK_BOX(adv_hbox),  scroller, FALSE, FALSE, 3);

    vbox = gtk_vbox_new(FALSE, 3);
    hbox_tools = gtk_hbox_new(FALSE, 0);
    new = ImageButton(GTK_STOCK_NEW);
    accept = ImageButton(GTK_STOCK_APPLY);
    save = ImageButton(GTK_STOCK_SAVE);
    delete = ImageButton(GTK_STOCK_DELETE); 
    gtk_box_pack_start(GTK_BOX(hbox_tools), new, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox_tools), accept, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox_tools), delete, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox_tools), save, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox_tools, FALSE, FALSE, 3);

    gtk_signal_connect(GTK_OBJECT(new), "clicked",
                       GTK_SIGNAL_FUNC(NewProfile), profiles_list);
    gtk_signal_connect(GTK_OBJECT(accept), "clicked",
                       GTK_SIGNAL_FUNC(UpdateProfile), tree);
    gtk_signal_connect(GTK_OBJECT(delete), "clicked",
                       GTK_SIGNAL_FUNC(DeleteProfile), tree);
    gtk_signal_connect(GTK_OBJECT(save), "clicked",
                       GTK_SIGNAL_FUNC(SaveProfileDatabase), profiles_list);

    tips = gtk_tooltips_new();
    gtk_tooltips_set_tip(tips, new, "create new profile", NULL);
    gtk_tooltips_set_tip(tips, accept, "accept changes to current profile", NULL);
    gtk_tooltips_set_tip(tips, delete, "delete current profile", NULL);
    gtk_tooltips_set_tip(tips, save, "save profile database", NULL);

    envparam[ENV_PROFILE] = MakeList(current_columns, &current_vbox, TRUE);
    gtk_box_pack_start(GTK_BOX(vbox), current_vbox, FALSE, FALSE, 3);
    gtk_widget_set_size_request(envparam[ENV_PROFILE], 100, 180);

    gtk_box_pack_start(GTK_BOX(adv_hbox), vbox, FALSE, FALSE, 3);

    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(tree));
    gtk_tree_selection_set_mode(select, GTK_SELECTION_SINGLE);
    g_signal_connect(G_OBJECT (select), "changed",
                     G_CALLBACK(ChangeSelectedProfile), envparam[ENV_PROFILE]);

    fname = DatabasePath("current.db");
    ReadCurrentDatabase(fname, profiles_list);
    g_free(fname);

    // the basic tab

    wave_vbox = gtk_vbox_new(FALSE, 3);

    envtable = table = MakeEntryTable(env_entries, NULL, NULL, 5, 5, 1);
    for (i = 0 ; i < 5 ; i++) {
        envparam[i] = table -> entries[i].w;
    }

    hbox = gtk_hbox_new(FALSE, 3);
    gtk_box_pack_start(GTK_BOX(hbox), table -> table, FALSE, FALSE, 3);

    current_table = gtk_table_new(5, 3, FALSE);
    gtk_table_attach(GTK_TABLE(current_table), gtk_label_new("constant"),
                     1, 2, 0, 1, 0, 0, 4, 0);
    gtk_table_attach(GTK_TABLE(current_table), gtk_label_new("profile"),
                     2, 3, 0, 1, 0, 0, 4, 0);

    gtk_table_attach(GTK_TABLE(current_table), gtk_label_new("x current"),
                     0, 1, 1, 2, 0, 0, 4, 0);

    envparam[ENV_XCURRENT] = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(envparam[ENV_XCURRENT]), 8);
    gtk_table_attach(GTK_TABLE(current_table), envparam[ENV_XCURRENT], 
                     1, 2, 1, 2, 0, 0, 4, 0);

    envparam[ENV_XCURR_PROFILE] = 
            gtk_combo_box_new_with_model(GTK_TREE_MODEL(profiles_list));

    gtk_table_attach(GTK_TABLE(current_table), envparam[ENV_XCURR_PROFILE], 
                     2, 3, 1, 2, 0, 0, 4, 0);

    renderer = gtk_cell_renderer_text_new();
    gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(envparam[ENV_XCURR_PROFILE]),
                               renderer, FALSE);
    gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(envparam[ENV_XCURR_PROFILE]),
                                   renderer, "text", 0, NULL);

    gtk_signal_connect(GTK_OBJECT(envparam[ENV_XCURRENT]), "changed",
                       GTK_SIGNAL_FUNC(ClearProfileBox), envparam[ENV_XCURR_PROFILE]);
    gtk_signal_connect(GTK_OBJECT(envparam[ENV_XCURR_PROFILE]), "changed", 
                       GTK_SIGNAL_FUNC(ClearCurrentEntry), envparam[ENV_XCURRENT]);

    gtk_table_attach(GTK_TABLE(current_table), gtk_label_new("y current"),
                     0, 1, 2, 3, 0, 0, 4, 0);


    envparam[ENV_YCURRENT] = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(envparam[ENV_YCURRENT]), 8);
    gtk_table_attach(GTK_TABLE(current_table), envparam[ENV_YCURRENT], 
                     1, 2, 2, 3, 0, 0, 4, 0);

    envparam[ENV_YCURR_PROFILE] = 
            gtk_combo_box_new_with_model(GTK_TREE_MODEL(profiles_list));

    renderer = gtk_cell_renderer_text_new();
    gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(envparam[ENV_YCURR_PROFILE]),
                               renderer, FALSE);
    gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(envparam[ENV_YCURR_PROFILE]),
                                   renderer, "text", 0, NULL);

    gtk_table_attach(GTK_TABLE(current_table), envparam[ENV_YCURR_PROFILE], 
                     2, 3, 2, 3, 0, 0, 4, 0);

    gtk_signal_connect(GTK_OBJECT(envparam[ENV_YCURRENT]), "changed",
                       GTK_SIGNAL_FUNC(ClearProfileBox), envparam[ENV_YCURR_PROFILE]);
    gtk_signal_connect(GTK_OBJECT(envparam[ENV_YCURR_PROFILE]), "changed", 
                       GTK_SIGNAL_FUNC(ClearCurrentEntry), envparam[ENV_YCURRENT]);

    envparam[ENV_CURRENT_SCALE] = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(envparam[ENV_CURRENT_SCALE]), 8);
    gtk_table_attach(GTK_TABLE(current_table), envparam[ENV_CURRENT_SCALE], 
                     1, 2, 3, 4, 0, 0, 4, 0);
    gtk_table_attach(GTK_TABLE(current_table), gtk_label_new("scale"),
                     0, 1, 3, 4, 0, 0, 4, 0);

    envparam[ENV_YWIND] = gtk_entry_new();
    envparam[ENV_XWIND] = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(envparam[ENV_YWIND]), 8);
    gtk_entry_set_width_chars(GTK_ENTRY(envparam[ENV_XWIND]), 8);
    gtk_table_attach(GTK_TABLE(current_table), envparam[ENV_XWIND], 
                     1, 2, 4, 5, 0, 0, 4, 0);
    gtk_table_attach(GTK_TABLE(current_table), envparam[ENV_YWIND], 
                     1, 2, 5, 6, 0, 0, 4, 0);
    gtk_table_attach(GTK_TABLE(current_table), gtk_label_new("x wind"),
                     0, 1, 4, 5, 0, 0, 4, 0);
    gtk_table_attach(GTK_TABLE(current_table), gtk_label_new("y wind"),
                     0, 1, 5, 6, 0, 0, 4, 0);
   
    gtk_box_pack_start(GTK_BOX(hbox), current_table, FALSE, FALSE, 3);

 
    gtk_box_pack_start(GTK_BOX(wave_vbox), hbox, FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(wave_vbox), gtk_hseparator_new(), FALSE, FALSE, 3);

    hbox = gtk_hbox_new(FALSE, 3);
    envparam[ENV_FORCING] = forcing_type = gtk_combo_box_new_text();
    gtk_combo_box_append_text(GTK_COMBO_BOX(forcing_type), "wave follower");
    gtk_combo_box_append_text(GTK_COMBO_BOX(forcing_type), "imposed motion");
    gtk_combo_box_append_text(GTK_COMBO_BOX(forcing_type), "morison");
    gtk_combo_box_set_active(GTK_COMBO_BOX(forcing_type), 0);

    envparam[ENV_RANDOM] = wave_random = gtk_radio_button_new_with_label(NULL, "random");
    envparam[ENV_REGULAR] = wave_regular = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(wave_random), "regular");

    gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new("forcing type"), FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(hbox), forcing_type, FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(hbox), wave_random, FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(hbox), wave_regular, FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(wave_vbox), hbox, FALSE, FALSE, 3);

    envparam[ENV_AMP]    = amplitude = gtk_entry_new();
    envparam[ENV_PERIOD] = period    = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(period), 12);
    gtk_entry_set_width_chars(GTK_ENTRY(amplitude), 12);

    hbox = gtk_hbox_new(FALSE, 1);
    gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new("wave amplitude"), FALSE, FALSE, 8);
    gtk_box_pack_start(GTK_BOX(hbox), amplitude, FALSE, FALSE, 1);
    gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new("period"), FALSE, FALSE, 1);
    gtk_box_pack_start(GTK_BOX(hbox), period, FALSE, FALSE, 1);

    gtk_box_pack_start(GTK_BOX(wave_vbox), hbox, FALSE, FALSE, 1);

    gtk_box_pack_start(GTK_BOX(basic_hbox), wave_vbox, FALSE, FALSE, 3);

    gtk_signal_connect(GTK_OBJECT(envparam[ENV_DEPTH]), "focus-out-event",
                       GTK_SIGNAL_FUNC(ComputeTotalsForDepth), NULL);


    gtk_widget_show_all(frame);

    return frame;
}

void 
ClearEnvironment(void)
{
    ClearEntryTable(envtable);
    ClearList(envparam[ENV_PROFILE]);
}

extern int gfDuplicatePref;
extern int DuplicateDialog(char *, char *);

static void
FillCurrent(VarExpr U, GtkWidget *curr, GtkWidget *prof)
{
    GtkTreeIter   iter;
    char          name[256];
    char         *buff;
    int          i;
    int          answer;

    if (!U.expr || IsConstant(U.expr)) {
        SetNumericText(curr, "%.3f", U.value, TRUE);
    }
    else if (U.name) {
        i = 0;
        strncpy(name, U.name, 256);
        while (SearchTreeRoot(GTK_TREE_MODEL(profiles_list), 
                              name, &iter) == TRUE) {
            i ++;
            snprintf(name, 256, "%s_%d", U.name, i);
        }
        
        if (i > 0) {
            if (gfDuplicatePref)
                answer = gfDuplicatePref;
            else
                answer = DuplicateDialog(U.name, "Profile");

            if (answer < 0) {
                answer = gfDuplicatePref = abs(answer);
            }

            if (answer == RENAME) { 
                free(U.name);
                U.name = strdup(name);
                AddProfile(profiles_list, U.name, U.text, &iter);
            }
            else {   // use global
                SearchTreeRoot(GTK_TREE_MODEL(profiles_list), U.name, &iter);
            }
        }
        else {
            AddProfile(profiles_list, U.name, U.text, &iter);
        }

        gtk_combo_box_set_active_iter(GTK_COMBO_BOX(prof), &iter);
    }
    else if (U.expr) {
        buff = UnusedProfileName(GTK_TREE_MODEL(profiles_list));
        AddProfile(profiles_list, buff, U.text, &iter);
        gtk_combo_box_set_active_iter(GTK_COMBO_BOX(prof), &iter);
    }
}

void
FillEnvironment(Environment *env)
{
    int           component;

    SetNumericText(envparam[ENV_DEPTH], "%.0f", env -> depth, TRUE);
    SetNumericText(envparam[ENV_DENSITY], "%.1f", env -> rho, TRUE);
    SetNumericText(envparam[ENV_GRAVITY], "%.2f", env -> gravity, TRUE);
    SetNumericText(envparam[ENV_STIFFNESS], "%.0f", env -> bottom_stiffness, TRUE);
    SetNumericText(envparam[ENV_DAMPING], "%.0f", env -> bottom_damping, TRUE);
    SetNumericText(envparam[ENV_XWIND], "%.1f", env -> y_wind.value, TRUE);
    SetNumericText(envparam[ENV_YWIND], "%.1f", env -> z_wind.value, TRUE);
    SetNumericText(envparam[ENV_CURRENT_SCALE], "%.3f", env -> Uscale, TRUE);

    component = 0;

    switch(env -> forcing) {

    case Velocity:
        gtk_combo_box_set_active(GTK_COMBO_BOX(envparam[ENV_FORCING]), 1);
        component = 1; // input in internal x-direction
        break;
    case WaveFollower:
        gtk_combo_box_set_active(GTK_COMBO_BOX(envparam[ENV_FORCING]), 0);
        component = 2; // wave in internal y-direction
        break;
    case Morison:
        gtk_combo_box_set_active(GTK_COMBO_BOX(envparam[ENV_FORCING]), 2);
        component = 2; // wave in internal y-direction
        break;
    default:
        gtk_combo_box_set_active(GTK_COMBO_BOX(envparam[ENV_FORCING]), -1);
        break;
    }

    if (component) {
        SetNumericText(envparam[ENV_AMP], "%.2f", env -> amplitude[component][0], TRUE);
        SetNumericText(envparam[ENV_PERIOD], "%.2f", env -> period[component][0], TRUE);
    }

    if (env -> input_type == Regular) {
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(envparam[ENV_REGULAR]), 1);
    }
    else {
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(envparam[ENV_RANDOM]), 1);
    }

    FillCurrent(env -> Uy, envparam[ENV_XCURRENT], envparam[ENV_XCURR_PROFILE]);
    FillCurrent(env -> Uz, envparam[ENV_YCURRENT], envparam[ENV_YCURR_PROFILE]);


}

void
WriteEnvironment(FILE *fp) 
{
    char *name, *expr;
    GtkTreeIter iter;

    int forcing;
    char *methods[] = {"wave-follower", "velocity", "morison"};

    fprintf(fp, "Environment\n");
    WriteEntry(envparam[ENV_DEPTH], "depth", fp);
    WriteEntry(envparam[ENV_DENSITY], "density", fp);
    WriteEntry(envparam[ENV_GRAVITY], "gravity", fp);
    WriteEntry(envparam[ENV_STIFFNESS], "bottom-stiffness", fp);
    WriteEntry(envparam[ENV_DAMPING], "bottom-damping", fp);

    if (gtk_combo_box_get_active_iter(GTK_COMBO_BOX(envparam[ENV_XCURR_PROFILE]), &iter)) {
        gtk_tree_model_get(GTK_TREE_MODEL(profiles_list), &iter, 0, &name, 1, &expr, -1);
        fprintf(fp, "    x-current = {\"%s\"} %s\n", name, expr);
    }
    else {
        WriteEntry(envparam[ENV_XCURRENT], "x-current", fp);
    }

    if (gtk_combo_box_get_active_iter(GTK_COMBO_BOX(envparam[ENV_YCURR_PROFILE]), &iter)) {
        gtk_tree_model_get(GTK_TREE_MODEL(profiles_list), &iter, 0, &name, 1, &expr, -1);
        fprintf(fp, "    y-current = {\"%s\"} %s\n", name, expr);
    }
    else {
        WriteEntry(envparam[ENV_YCURRENT], "y-current", fp);
    }

    WriteEntry(envparam[ENV_XWIND], "x-wind", fp);
    WriteEntry(envparam[ENV_YWIND], "y-wind", fp);
    WriteEntry(envparam[ENV_CURRENT_SCALE], "current-scale", fp);

    forcing = gtk_combo_box_get_active(GTK_COMBO_BOX(envparam[ENV_FORCING]));
    if (forcing >= 0 && forcing <= 2) {
        fprintf(fp, "    forcing-method = %s\n", methods[forcing]);
    }
    fprintf(fp, "    input-type = %s\n", 
            gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(envparam[ENV_RANDOM])) ? "random" : "regular");
    if (forcing == 0 || forcing == 2) {
        fprintf(fp, "    x-wave = (%s, %s, 0)\n", 
                gtk_entry_get_text(GTK_ENTRY(envparam[ENV_AMP])), 
                gtk_entry_get_text(GTK_ENTRY(envparam[ENV_PERIOD])));
    }
    else if (forcing == 1) {
        fprintf(fp, "    z-input = (%s, %s)\n", 
                gtk_entry_get_text(GTK_ENTRY(envparam[ENV_AMP])), 
                gtk_entry_get_text(GTK_ENTRY(envparam[ENV_PERIOD])));
    }
}
