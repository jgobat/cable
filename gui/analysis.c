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
#include "entry.h"
#include "problem.h"
#include "list.h"
#include "analysis.h"
#include "text.h"
#include "combo.h"
#include "entry.h"

extern GtkWidget *toplevel;
extern Problem problem;

enum {
    STATIC_RELAXATION,
    STATIC_TOLERANCE,
    STATIC_ITERATIONS,
    OUTER_RELAXATION,
    OUTER_TOLERANCE,
    OUTER_ITERATIONS,
    DYNAMIC_RELAXATION,
    DYNAMIC_TOLERANCE,
    DYNAMIC_ITERATIONS,
    NUM_BASIC_PARAMETERS,
};

static EntryTableEntry param_entries[] =
{
    { NULL, "static relaxation factor (0-1)", NUMBER_ENTRY, 1 },
    { NULL, "tolerance for static convergence", NUMBER_ENTRY, 1 },
    { NULL, "iteration limit in static solver", NUMBER_ENTRY, 1 },
    { NULL,  "outer iteration relaxation factor", NUMBER_ENTRY, 1 },
    { NULL,  "tolerance for outer iteration convergence",NUMBER_ENTRY, 1 },
    { NULL,  "iteration limit for outer iterations",NUMBER_ENTRY, 1 },
    { NULL,  NULL,NUMBER_ENTRY, 1 },
    { NULL,  NULL,NUMBER_ENTRY, 1 },
    { NULL,  NULL,NUMBER_ENTRY, 1 },
};

static EntryTableEntry steps_entries[] =
{
    { NULL, "dynamic solution time step", NUMBER_ENTRY, 1 },
    { NULL, "output time step for time series results at output nodes", NUMBER_ENTRY, 1 },
    { NULL, "output time step for snapshot animation results", NUMBER_ENTRY, 1 },
};

GtkWidget *anparam[NUM_ANALYSIS];
static EntryTable *steps_table, *param_table;

enum {
    InitialCat = 0x01,
    InitialShoot = 0x02,
    InitialFile  = 0x04,
    StaticShoot = 0x08,
    StaticRelax = 0x10,
    StaticAuto = 0x20,
    StaticFile = 0x40,
    Dynamics = 0x80,
};

typedef struct { 
    char *label;
    int   flags;
} SolutionType;

static SolutionType solution_types[] = 
{
    {"catenary initial + static relaxation", InitialCat | StaticRelax },
    {"catenary initial + static relaxation + dynamics", 
                                             InitialCat | StaticRelax | Dynamics },
    {"shooting initial + static relaxation", InitialShoot | StaticRelax }, 
    {"shooting initial + static relaxation + dynamics", 
                                           InitialShoot | StaticRelax | Dynamics },
    {"shooting static", StaticShoot },
    {"shooting static + dynamics", StaticShoot | Dynamics },
    {"initial from file + static relaxation", InitialFile | StaticRelax },
    {"initial from file + static relaxation + dynamics", 
                                            InitialFile | StaticRelax | Dynamics },
    {"static from file + dynamics", StaticFile | Dynamics },
    {"autosolve static", StaticAuto },
    {"autosolve static + dynamics", StaticAuto | Dynamics },
};

static void
ShowLoadFile(GtkComboBox *w, gpointer hbox)
{
    gchar *solution;

    solution = gtk_combo_box_get_active_text(w);
    if (solution == NULL) {
        return;
    }

    if (strstr(solution, "file")) {
        gtk_widget_show(GTK_WIDGET(hbox));
    }
    else {
        gtk_widget_hide(GTK_WIDGET(hbox));
    }

    g_free(solution);
}

static void
BrowseLoadFile(GtkComboBox *w, gpointer data)
{
    GtkEntry *text = (GtkEntry *) data;
    char *filename;
    GtkWidget *dialog;

    dialog = gtk_file_chooser_dialog_new ("load solution from file", 
                      GTK_WINDOW(toplevel),
                      GTK_FILE_CHOOSER_ACTION_OPEN,
                      GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                      GTK_STOCK_OK, GTK_RESPONSE_ACCEPT,
                      NULL);

    if (gtk_dialog_run (GTK_DIALOG (dialog)) == GTK_RESPONSE_ACCEPT) {
        filename = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (dialog));
        gtk_entry_set_text(text, filename);
        g_free (filename);
    }

    gtk_widget_destroy (dialog);
}

GtkWidget *
BuildAnalysisParameters(void)
{
    int          i;
    static char *rows[] = {"static", "outer", "dynamic"};
    static char *cols[] = {"relaxation", "tolerance", "iterations"};
    static char *steps_cols[] = {"solution", "output", "snapshot"};
    static char *steps_rows[] = {"time-step"};
    GtkWidget   *frame;
    EntryTable  *basic, *steps;
    GtkWidget   *vbox, *hbox, *big_hbox, *big_vbox;
    GtkWidget   *nodes_box;
    GtkWidget   *soln_file_text;
    GtkWidget   *soln_file_browse;
    GtkWidget   *duration;
    GtkWidget   *solution_type, *soln_2d, *soln_3d;
    GtkWidget   *title, *type;
    int          num_types;
    GtkWidget   *anbook;
    static char *nodes_cols[] = {"nodes", "%d", NULL};

    // basic layout is a two tab notebook widget. Each tab has a 
    // a main vbox on it

    frame = gtk_frame_new("Model parameters");
    anbook = gtk_notebook_new();
    gtk_container_add(GTK_CONTAINER(frame), anbook);

    gtk_notebook_set_tab_pos(GTK_NOTEBOOK(anbook), GTK_POS_TOP);

    // the basic tab

    big_vbox = gtk_vbox_new(FALSE, 1);
    gtk_notebook_append_page(GTK_NOTEBOOK(anbook), big_vbox, gtk_label_new("basic"));


    big_hbox = gtk_hbox_new(FALSE, 1);
    gtk_box_pack_start(GTK_BOX(big_vbox), big_hbox, FALSE, FALSE, 1);

    vbox = gtk_vbox_new(FALSE, 1); 

    // a 3x3 table for relaxation factors, iteration limits, and tolerances

    basic = param_table = MakeEntryTable(param_entries, rows, cols, 9, 3, 3);
    gtk_box_pack_start(GTK_BOX(vbox), basic -> table, FALSE, FALSE, 1);
    for (i = 0 ; i < 9 ; i++) {
        anparam[i + ANALYSIS_STATIC_RELAX] = basic -> entries[i].w;
    }

    // duration entry and its label get their own hbox

    anparam[ANALYSIS_DURATION]  = duration = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(duration), 15);


    hbox = gtk_hbox_new(FALSE, 1);
    gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new("duration"), FALSE, FALSE, 2);
    gtk_box_pack_start(GTK_BOX(hbox), duration, FALSE, FALSE, 1);

    gtk_box_pack_start(GTK_BOX(vbox), gtk_hseparator_new(), FALSE, FALSE, 1);
   
    steps = steps_table = MakeEntryTable(steps_entries, steps_rows, steps_cols, 3, 1, 3); 
    anparam[ANALYSIS_DT] = steps -> entries[0].w;
    anparam[ANALYSIS_OUTPUT_DT] = steps -> entries[1].w;
    anparam[ANALYSIS_SNAPSHOT_DT] = steps -> entries[2].w;
    // gtk_entry_set_width_chars(GTK_ENTRY(time_step), 12);

    gtk_box_pack_start(GTK_BOX(vbox), steps -> table, FALSE, FALSE, 1);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 1);

    hbox = gtk_hbox_new(FALSE, 1);
    anparam[ANALYSIS_RAMP]  = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(anparam[ANALYSIS_RAMP]), 15);
    gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new("ramp time"), FALSE, FALSE, 2);
    gtk_box_pack_start(GTK_BOX(hbox), anparam[ANALYSIS_RAMP], FALSE, FALSE, 1);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 1);

    gtk_box_pack_start(GTK_BOX(big_hbox), vbox, FALSE, FALSE, 1);
    gtk_box_pack_start(GTK_BOX(big_hbox), gtk_vseparator_new(), FALSE, FALSE, 1);
    anparam[ANALYSIS_NODES] = MakeList(nodes_cols, &nodes_box, TRUE);
    gtk_widget_set_size_request(anparam[ANALYSIS_NODES], 0, 120);

    gtk_box_pack_start(GTK_BOX(big_hbox), nodes_box, FALSE, FALSE, 1);
    vbox = gtk_vbox_new(FALSE, 1);
    anparam[ANALYSIS_TERMINALS] = gtk_check_button_new_with_label("terminals");
    anparam[ANALYSIS_CONNECTORS] = gtk_check_button_new_with_label("connectors");
    anparam[ANALYSIS_FIRST] = gtk_check_button_new_with_label("first");
    anparam[ANALYSIS_LAST] = gtk_check_button_new_with_label("last");
    gtk_box_pack_start(GTK_BOX(vbox), anparam[ANALYSIS_TERMINALS], FALSE, FALSE, 1);
    gtk_box_pack_start(GTK_BOX(vbox), anparam[ANALYSIS_CONNECTORS], FALSE, FALSE, 1);
    gtk_box_pack_start(GTK_BOX(vbox), anparam[ANALYSIS_FIRST], FALSE, FALSE, 1);
    gtk_box_pack_start(GTK_BOX(vbox), anparam[ANALYSIS_LAST], FALSE, FALSE, 1);
    gtk_box_pack_start(GTK_BOX(big_hbox), vbox, FALSE, FALSE, 1);

    gtk_box_pack_start(GTK_BOX(big_vbox), gtk_hseparator_new(), FALSE, FALSE, 1);

    // problem type and title

    hbox = gtk_hbox_new(FALSE, 3);
    anparam[ANALYSIS_TYPE] = type = gtk_combo_box_new_text();
    gtk_combo_box_append_text(GTK_COMBO_BOX(type), "surface");
    gtk_combo_box_append_text(GTK_COMBO_BOX(type), "subsurface");
    gtk_combo_box_append_text(GTK_COMBO_BOX(type), "towing");
    gtk_combo_box_append_text(GTK_COMBO_BOX(type), "deployment");
    gtk_combo_box_append_text(GTK_COMBO_BOX(type), "horizontal");
    gtk_combo_box_set_active(GTK_COMBO_BOX(type), 0);
    anparam[ANALYSIS_TITLE] = title = gtk_entry_new();
    gtk_entry_set_width_chars(GTK_ENTRY(anparam[ANALYSIS_TITLE]), 30);
    gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new("problem type"), FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(hbox), type, FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new("title"), FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(hbox), title, FALSE, FALSE, 3);
    
    gtk_box_pack_start(GTK_BOX(big_vbox), hbox, FALSE, FALSE, 1);

    // a combo box to choose the solution type, radios to choose 2D/3D

    anparam[ANALYSIS_SOLUTION] = solution_type = gtk_combo_box_new_text();
    num_types = sizeof(solution_types) / sizeof(SolutionType);
    for (i = 0 ; i < num_types ; i++) {
        gtk_combo_box_append_text(GTK_COMBO_BOX(solution_type), solution_types[i].label);
    }

    gtk_combo_box_set_active(GTK_COMBO_BOX(solution_type), 0);

    anparam[ANALYSIS_DIMENSION] = soln_2d = gtk_radio_button_new_with_label(NULL, "2D");
    soln_3d = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(soln_2d), "3D");

    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(soln_2d), TRUE);

    hbox = gtk_hbox_new(FALSE, 3);
    gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new("solution type:"), FALSE, FALSE, 3);  
    gtk_box_pack_start(GTK_BOX(hbox), solution_type, FALSE, FALSE, 1);
    gtk_box_pack_start(GTK_BOX(hbox), soln_2d, FALSE, FALSE, 1);
    gtk_box_pack_start(GTK_BOX(hbox), soln_3d, FALSE, FALSE, 1);
    gtk_box_pack_start(GTK_BOX(big_vbox), hbox, FALSE, FALSE, 1);

    // the usually hidden entry and browse button for load solution files

    hbox = gtk_hbox_new(FALSE, 3);
    anparam[ANALYSIS_LOAD_FILE] = soln_file_text = gtk_entry_new();
    soln_file_browse = gtk_button_new_with_label("...");
    gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new("load file"), FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(hbox), soln_file_text, FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(hbox), soln_file_browse, FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(big_vbox), hbox, FALSE, FALSE, 1);
    
    gtk_widget_show_all(big_vbox);
    gtk_widget_hide(hbox);

    gtk_signal_connect(GTK_OBJECT(soln_file_browse), "clicked", GTK_SIGNAL_FUNC(BrowseLoadFile), soln_file_text);
    gtk_signal_connect(GTK_OBJECT(solution_type), "changed", GTK_SIGNAL_FUNC(ShowLoadFile), hbox);

    // the advanced parameters go on a new tab

    vbox = gtk_vbox_new(FALSE, 3);
    gtk_notebook_append_page(GTK_NOTEBOOK(anbook), vbox, gtk_label_new("advanced"));

    hbox = gtk_hbox_new(FALSE, 3);
    anparam[ANALYSIS_INTEGRATION_METHOD] = type = gtk_combo_box_new_text();
    gtk_combo_box_append_text(GTK_COMBO_BOX(type), "spatial");
    gtk_combo_box_append_text(GTK_COMBO_BOX(type), "velocity");
    gtk_combo_box_set_active(GTK_COMBO_BOX(type), 0);
    gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new("dynamic integration"), FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(hbox), type, FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 3);

    hbox = gtk_hbox_new(FALSE, 3);
    anparam[ANALYSIS_OUTER_METHOD] = type = gtk_combo_box_new_text();
    gtk_combo_box_append_text(GTK_COMBO_BOX(type), "bisection");
    gtk_combo_box_append_text(GTK_COMBO_BOX(type), "secant");
    gtk_combo_box_set_active(GTK_COMBO_BOX(type), 0);
    gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new("outer solution method"), FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(hbox), type, FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 3);

    hbox = gtk_hbox_new(FALSE, 3);
    anparam[ANALYSIS_VIVA_ITERATIONS] = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new("VIVA iterations"), 
                       FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(hbox), anparam[ANALYSIS_VIVA_ITERATIONS], 
                       FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 3);
    gtk_widget_show_all(vbox);

    gtk_widget_show(anbook);
    gtk_widget_show(frame);

    return frame;
}

void
FillAnalysis(Analysis *a, Problem *p)
{
    static char *integr [] = {"", "spatial", "temporal"};
    static char *outer [] = {"", "bisection", "secant"};

    SetNumericText(anparam[ANALYSIS_DURATION], "%.0f", a->duration, TRUE);
    SetNumericText(anparam[ANALYSIS_DT], "%.4f", a->dt, TRUE);

    SetNumericText(anparam[ANALYSIS_STATIC_RELAX], "%g", a->static_relaxation, TRUE);
    SetNumericText(anparam[ANALYSIS_STATIC_TOL], "%g", a->static_tolerance, TRUE);
    SetNumericText(anparam[ANALYSIS_STATIC_IT], "%.0f", (double) a->static_it, TRUE);


    SetNumericText(anparam[ANALYSIS_DYNAMIC_RELAX], "%.3f", a->dynamic_relaxation, TRUE);
    SetNumericText(anparam[ANALYSIS_DYNAMIC_TOL], "%g", a->dynamic_tolerance, TRUE);
    SetNumericText(anparam[ANALYSIS_DYNAMIC_IT], "%.0f", (double) a->dynamic_it, TRUE);

    SetNumericText(anparam[ANALYSIS_RAMP], "%g", a->ramp_time, TRUE);

    SetNumericText(anparam[ANALYSIS_OUTER_RELAX], "%.3f", a->outer_relaxation, TRUE);
    SetNumericText(anparam[ANALYSIS_OUTER_TOL], "%g", a->outer_tolerance, TRUE);
    SetNumericText(anparam[ANALYSIS_OUTER_IT], "%.0f", (double) a->outer_it, TRUE);

    SetNumericText(anparam[ANALYSIS_VIVA_ITERATIONS], "%.0f", (double) a->viva_iterations, TRUE);
    ComboBoxSetText(GTK_COMBO_BOX(anparam[ANALYSIS_TYPE]), ProblemTypeName(p -> type), 0);
    ComboBoxSetText(GTK_COMBO_BOX(anparam[ANALYSIS_OUTER_METHOD]), outer[a -> static_outer_method], 0);
    ComboBoxSetText(GTK_COMBO_BOX(anparam[ANALYSIS_INTEGRATION_METHOD]), integr[a -> integration], 0);
    
    if (!p -> title || strlen(p -> title) == 0)
        gtk_entry_set_text(GTK_ENTRY(anparam[ANALYSIS_TITLE]), "untitled");
    else
        gtk_entry_set_text(GTK_ENTRY(anparam[ANALYSIS_TITLE]), p -> title);

    // following is hand maintained along with the solution_types array

    if (a -> static_initial_guess == Catenary && a -> static_solution == Relaxation)
        ComboBoxSetText(GTK_COMBO_BOX(anparam[ANALYSIS_SOLUTION]), solution_types[0].label, 0);
    else if (a -> static_initial_guess == Shooting && a -> static_solution == Relaxation)
        ComboBoxSetText(GTK_COMBO_BOX(anparam[ANALYSIS_SOLUTION]), solution_types[2].label, 0);
    else if (a -> static_solution == Shooting) 
        ComboBoxSetText(GTK_COMBO_BOX(anparam[ANALYSIS_SOLUTION]), solution_types[4].label, 0);
   
    return;
}

char *
AnalysisProblemType(void)
{
    return ComboBoxGetText(GTK_COMBO_BOX(anparam[ANALYSIS_TYPE]), -1);
}

void
WriteAnalysis(FILE *fp)
{
    int i;

    fprintf(fp, "Problem description\n");
    fprintf(fp, "    title = \"%s\"\n", 
            gtk_entry_get_text(GTK_ENTRY(anparam[ANALYSIS_TITLE])));
    fprintf(fp, "    type = %s\n", 
            ComboBoxGetText(GTK_COMBO_BOX(anparam[ANALYSIS_TYPE]), -1));

    fprintf(fp, "\nAnalysis parameters\n");
    WriteEntry(anparam[ANALYSIS_STATIC_RELAX], "static-relaxation", fp);
    WriteEntry(anparam[ANALYSIS_STATIC_IT], "static-iterations", fp);
    WriteEntry(anparam[ANALYSIS_STATIC_TOL], "static-tolerance", fp);
    WriteEntry(anparam[ANALYSIS_OUTER_RELAX], "static-outer-relaxation", fp);
    WriteEntry(anparam[ANALYSIS_OUTER_IT], "static-outer-iterations", fp);
    WriteEntry(anparam[ANALYSIS_OUTER_TOL], "static-outer-tolerance", fp);
    WriteEntry(anparam[ANALYSIS_DYNAMIC_RELAX], "dynamic-relaxation", fp);
    WriteEntry(anparam[ANALYSIS_DYNAMIC_IT], "dynamic-iterations", fp);
    WriteEntry(anparam[ANALYSIS_DYNAMIC_TOL], "dynamic-tolerance", fp);
    WriteEntry(anparam[ANALYSIS_DT], "time-step", fp);
    WriteEntry(anparam[ANALYSIS_DURATION], "duration", fp);
    WriteEntry(anparam[ANALYSIS_RAMP], "ramp-time", fp);
    WriteEntry(anparam[ANALYSIS_VIVA_ITERATIONS], "viva-iterations", fp);
    fprintf(fp, "    static-outer-method = %s\n", 
            ComboBoxGetText(GTK_COMBO_BOX(anparam[ANALYSIS_OUTER_METHOD]), -1));
    fprintf(fp, "    dynamic-integration = %s\n", 
            ComboBoxGetText(GTK_COMBO_BOX(anparam[ANALYSIS_INTEGRATION_METHOD]), -1));

    i = gtk_combo_box_get_active(GTK_COMBO_BOX(anparam[ANALYSIS_SOLUTION]));
    if (solution_types[i].flags & InitialCat)
        fprintf(fp, "    static-initial-guess = catenary\n");
    else if (solution_types[i].flags & InitialShoot)
        fprintf(fp, "    static-initial-guess = shooting\n");

    if (solution_types[i].flags & StaticRelax)
        fprintf(fp, "    static-solution = relaxation\n");
    else if (solution_types[i].flags & StaticShoot)
        fprintf(fp, "    static-solution = shooting\n");
}

void
FillSolutionControl(Solution *s, int num_nodes)
{
    int i;

    s -> twoD = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(anparam[ANALYSIS_DIMENSION]));
    i = gtk_combo_box_get_active(GTK_COMBO_BOX(anparam[ANALYSIS_SOLUTION]));
    s -> static_only = !(solution_types[i].flags & Dynamics);

    s -> output_special = 0;
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(anparam[ANALYSIS_TERMINALS])))
        s -> output_special += OUTPUT_TERMINALS;
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(anparam[ANALYSIS_CONNECTORS])))
        s -> output_special += OUTPUT_CONNECTORS;
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(anparam[ANALYSIS_FIRST])))
        s -> output_special += OUTPUT_FIRST;
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(anparam[ANALYSIS_LAST])))
        s -> output_special += OUTPUT_LAST;

    s -> output_nodes = ListToArray(anparam[ANALYSIS_NODES], 
                                    &s -> num_output_nodes, 1, num_nodes);
    
    s -> initial_file = s -> static_file = NULL;
    if (solution_types[i].flags & InitialFile)
        s -> initial_file = (char *) gtk_entry_get_text(GTK_ENTRY(anparam[ANALYSIS_LOAD_FILE]));
    else if (solution_types[i].flags & StaticFile)     
        s -> static_file = (char *) gtk_entry_get_text(GTK_ENTRY(anparam[ANALYSIS_LOAD_FILE]));
    
    if (solution_types[i].flags & StaticAuto)
        s -> auto_static = 1;
    else
        s -> auto_static = 0;

    s -> output_dt = atof(gtk_entry_get_text(GTK_ENTRY(anparam[ANALYSIS_OUTPUT_DT])));
    s -> snapshot_dt = atof(gtk_entry_get_text(GTK_ENTRY(anparam[ANALYSIS_SNAPSHOT_DT])));
    fprintf(stderr,"snapdhot_dt at fill = %g\n", s -> snapshot_dt);
}

void
ClearAnalysis(void)
{
    ClearEntryTable(param_table);
    ClearEntryTable(steps_table);
    ClearList(anparam[ANALYSIS_NODES]);
    gtk_entry_set_text(GTK_ENTRY(anparam[ANALYSIS_DURATION]), "");
    gtk_entry_set_text(GTK_ENTRY(anparam[ANALYSIS_TITLE]), "");
}
