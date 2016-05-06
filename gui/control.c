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
#include <string.h>
#include <gtk/gtk.h>
#include <glib.h>
#include "localerror.h"
#include "entry.h"
#include "dialog.h"
#include "problem.h"
#include "plot.h"
#include "tension.h"
#include "segments.h"
#include "localcontrol.h"
#include "results.h"
#include "pixbuf2mpeg.h"
#include "text.h"
#include "plotutil.h"

extern void TensionToGlobal(Node node, int twoD, double *Fx, double *Fy, double *Fz);

extern GtkWidget *toplevel;
extern GKeyFile *prefsFile;

static Analysis *baseline = NULL;
static Analysis *aptr = NULL;

static GtkWidget **active_ctlparam = NULL;

typedef struct {
    GtkWidget *toplevel;
    int frameN;
    int interval;
    guint timer;
    GtkWidget *display;
} Playback;
    
static EntryTableEntry relaxation_entries[] =
{
    { NULL, NULL, NUMBER_ENTRY, 1 },
    { NULL, NULL, NUMBER_ENTRY, 1 },
    { NULL, NULL, NUMBER_ENTRY, 1 },
    { NULL, NULL, NUMBER_ENTRY, 1 },
    { NULL, NULL, NUMBER_ENTRY, 1 },
    { NULL, NULL, NUMBER_ENTRY, 1 },
    { NULL, NULL, NUMBER_ENTRY, 1 },
    { NULL, NULL, NUMBER_ENTRY, 1 },
    { NULL, NULL, NUMBER_ENTRY, 1 },
};

static EntryTableEntry dynamic_entries[] =
{
    { NULL, NULL, NUMBER_ENTRY, 1 },
    { NULL, NULL, NUMBER_ENTRY, 1 },
    { NULL, NULL, NUMBER_ENTRY, 1 },
};

static EntryTableEntry runtime_entries[] =
{
    { NULL, NULL, NUMBER_ENTRY, 1 },
    { NULL, NULL, NUMBER_ENTRY, 1 },
    { NULL, NULL, NUMBER_ENTRY, 1 },
};


enum {
    CONTROL_TOPLEVEL,
    CONTROL_STATIC_RELAX,
    CONTROL_STATIC_TOL,
    CONTROL_STATIC_IT,
    CONTROL_OUTER_RELAX,
    CONTROL_OUTER_TOL,
    CONTROL_OUTER_IT,
    CONTROL_DYNAMIC_RELAX,
    CONTROL_DYNAMIC_TOL,
    CONTROL_DYNAMIC_IT,
    CONTROL_DT,
    CONTROL_DURATION,
    CONTROL_RAMP_TIME,
    CONTROL_TIME,
    CONTROL_ITERATION,
    CONTROL_ERROR,
    CONTROL_MESSAGE,
    NUM_CONTROL,
};

int
FileCopy(char *src, char *dest)
{
    FILE *in, *out;
    unsigned char *buff[4096];
    int            n;

    in = fopen(src, "rb");
    out = fopen(dest, "wb");
    if (!in || !out) {
        error("error creating solution copy");
        return 1;
    }

    while(!feof(in)) {
        n = fread(buff, sizeof(unsigned char), 4096, in);    
        if (n) {
            fwrite(buff, sizeof(unsigned char), n, out);
        }
        if (n < 4096) {
            break;
        }
    }

    fclose(in);
    fclose(out);

    return 0;
}

extern int SolutionToMatlab(char *, char *, int, int, int, int, int);
extern char *prev_file;

void 
MatlabSolution(GtkWidget *w, gpointer data) 
{
    Solution *s = (Solution *) data;
    char        *filename;

    filename = FilenameFromDialog("convert results to Matlab", 
                                  gtk_widget_get_toplevel(w),
                                  GTK_FILE_CHOOSER_ACTION_SAVE, &prev_file, 
                                  "*.mat", "Matlab data files", NULL, NULL); 


    if (!filename) {
        return;
    }

    filename = AddFilenameExtension(filename, ".mat");
    SolutionToMatlab(s -> out_name, filename, 0, 0, 1, 1, 0);    
}

void
SaveSolution(GtkWidget *w, gpointer data)
{
    Solution *s = (Solution *) data;
    char        *filename;

    s -> saved = 1;
    filename = FilenameFromDialog("save results", gtk_widget_get_toplevel(w),
                                  GTK_FILE_CHOOSER_ACTION_SAVE, &prev_file, 
                                  "*.crs", "cable results files", NULL, NULL); 

    if (!filename) {
        return;
    }

    filename = AddFilenameExtension(filename, ".crs");
    FileCopy(s -> out_name, filename);    
}

int
CheckSaved(Solution *solution, GtkWidget *win)
{
    GtkWidget   *dialog;
    gint         ans;

    ans = GTK_RESPONSE_YES;
    if (!solution -> saved) {
        dialog = gtk_message_dialog_new(GTK_WINDOW(win),
                                        GTK_DIALOG_MODAL,
                                        GTK_MESSAGE_QUESTION,
                                        GTK_BUTTONS_YES_NO,
                                        "Solution not saved. Save now?");
        gtk_dialog_add_button(GTK_DIALOG(dialog), "Cancel", GTK_RESPONSE_CANCEL);
        ans = gtk_dialog_run(GTK_DIALOG (dialog));
        gtk_widget_destroy(dialog);
        
        if (ans == GTK_RESPONSE_YES) {
            SaveSolution(win, solution);
        }
    }

    return ans;
}


void
QuitSolution(GtkWidget *w, gpointer data)
{
    Solution    *solution = (Solution *) data;

    if (CheckSaved(solution, gtk_widget_get_toplevel(w)) == GTK_RESPONSE_CANCEL)
        return;

    fprintf(stderr, "QuitSolution\n");
    solution -> userQuit = 1;
    if (solution -> solutionComplete) {
        ControlDialogQuit(solution);
    }

    return;
}

gboolean
DeleteSolution(GtkWidget *w, GdkEvent *event, gpointer data)
{
    fprintf(stderr,"DeleteSolution\n");
    QuitSolution(w, data);
    return TRUE;
}

void
ControlDialogQuit(Solution *s)
{
    if (s && s -> out_name && s -> playback == NULL) {
        unlink(s -> out_name);
        s -> out_name = NULL;
    }
    if (s && s -> controls == active_ctlparam) {
        active_ctlparam = NULL;
    }
    // controls[0] is the toplevel win
    if (s -> controls) {
        fprintf(stderr,"destroying solution control\n");
        gtk_widget_destroy(GTK_WIDGET(((GtkWidget **) (s -> controls))[0]));
        s -> controls = NULL;
    }
    if (s && !s -> playback) {
        fprintf(stderr,"freeing solution\n");
        FreeSolution(s);
    }
}

void
PauseSolution(GtkWidget *w, gpointer data)
{
    while(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {  
        ControlProcessEvents((Solution *) data);
    }
}

static void
FillControl(Analysis *a, GtkWidget **ctlparam)
{
    SetNumericText(ctlparam[CONTROL_STATIC_RELAX], "%g", a->static_relaxation, TRUE);
    SetNumericText(ctlparam[CONTROL_STATIC_TOL], "%g", a->static_tolerance, TRUE);
    SetNumericText(ctlparam[CONTROL_STATIC_IT], "%.0f", (double) a->static_it, TRUE);
    
    SetNumericText(ctlparam[CONTROL_DYNAMIC_RELAX], "%.3f", a->dynamic_relaxation, TRUE);
    SetNumericText(ctlparam[CONTROL_DYNAMIC_TOL], "%g", a->dynamic_tolerance, TRUE);
    SetNumericText(ctlparam[CONTROL_DYNAMIC_IT], "%.0f", (double) a->dynamic_it, TRUE);

    SetNumericText(ctlparam[CONTROL_OUTER_RELAX], "%.3f", a->outer_relaxation, TRUE);

    SetNumericText(ctlparam[CONTROL_OUTER_TOL], "%g", a->outer_tolerance, TRUE);
    SetNumericText(ctlparam[CONTROL_OUTER_IT], "%.0f", (double) a->outer_it, TRUE);

    SetNumericText(ctlparam[CONTROL_DURATION], "%g", a->duration, TRUE);
    SetNumericText(ctlparam[CONTROL_DT], "%g", a->dt, TRUE);
    SetNumericText(ctlparam[CONTROL_RAMP_TIME], "%g", a->ramp_time, TRUE);
   
    return; 
}

void
UpdateParameters(GtkWidget *w, gpointer data)
{
    GtkWidget **ctlparam = (GtkWidget **) data;
    double  value;

    value = atof(gtk_entry_get_text(GTK_ENTRY(ctlparam[CONTROL_STATIC_RELAX])));
    if (value > -1.0 && value <= 1.0) 
        aptr -> static_relaxation = value;

    value = atof(gtk_entry_get_text(GTK_ENTRY(ctlparam[CONTROL_STATIC_TOL])));
    if (value > 0)
        aptr->static_tolerance = value;

    value = atof(gtk_entry_get_text(GTK_ENTRY(ctlparam[CONTROL_STATIC_IT])));
    if (value > 0)
        aptr->static_it = value;

    value = atof(gtk_entry_get_text(GTK_ENTRY(ctlparam[CONTROL_DYNAMIC_RELAX])));
    if (value > 0 && value <= 1.0)
        aptr->dynamic_relaxation = value;

    value = atof(gtk_entry_get_text(GTK_ENTRY(ctlparam[CONTROL_DYNAMIC_TOL])));
    if (value > 0)
        aptr->dynamic_tolerance = value;

    value = atof(gtk_entry_get_text(GTK_ENTRY(ctlparam[CONTROL_DYNAMIC_IT])));
    if (value > 0)
        aptr->dynamic_it = value;

    value = atof(gtk_entry_get_text(GTK_ENTRY(ctlparam[CONTROL_OUTER_RELAX])));
    if (value)
        aptr->outer_relaxation = value;

    value = atof(gtk_entry_get_text(GTK_ENTRY(ctlparam[CONTROL_OUTER_TOL])));
    if (value > 0)
        aptr->outer_tolerance = value;

    value = atof(gtk_entry_get_text(GTK_ENTRY(ctlparam[CONTROL_OUTER_IT])));
    if (value > 0)
        aptr->outer_it = value;

    value = atof(gtk_entry_get_text(GTK_ENTRY(ctlparam[CONTROL_DT])));
    if (value > 0)
        aptr->dt = value;

    value = atof(gtk_entry_get_text(GTK_ENTRY(ctlparam[CONTROL_RAMP_TIME])));
    if (value >= 0)
        aptr->ramp_time = value;

    FillControl(aptr, ctlparam);
}


void
RestoreParameters(GtkWidget *w, gpointer data)
{
    if (!baseline || !aptr) {
        return;
    }

    memcpy(aptr, baseline, sizeof(Analysis));
    FillControl(aptr, (GtkWidget **) data);
}

void
PositionResult(Solution *s, GtkWidget *w)
{
    static Solution *prev_s = NULL;
    static int x, y;

    if (s != prev_s) {
        
        if (s -> controls)
            gtk_window_get_position(GTK_WINDOW(((GtkWidget **) (s -> controls))[0]), &x, &y);
        else if (s -> playback)
            gtk_window_get_position(GTK_WINDOW(((Playback *) s -> playback) -> toplevel), &x, &y);
        else {
            x = y = 0;
        }

        prev_s = s;
    }

    x += 20;
    y += 20;
    gtk_window_move(GTK_WINDOW(w), x, y);
}

void
ControlMakePlots(Solution *s, int flags)
{
    double      params[4];
    GError     *err = NULL;
    int         nsamples;

    params[0] = 0;
    params[1] = s -> analysis -> duration;
    params[2] = 0;
    params[3] = 0;

    nsamples = 0;
    if (s -> output_dt) {
        nsamples = (int) ((s -> analysis -> duration + s -> output_dt/2.0) / s -> output_dt) + 1;
    }
 

    if (nsamples && (flags & PLOT_TIME) &&
        g_key_file_get_boolean(prefsFile, "output", "tensionTimeSeries", &err)) {
        printf("making force time plot %d %d\n",
               err == NULL ? 0 : 1, g_key_file_get_boolean(prefsFile, "output", "tensionTimeSeries", &err));

        s -> tsPlot[FORCE] = CreateChartDisplay(s -> problem -> title, 
                                                CHART_XY, 
                                                s -> num_output_nodes, 
                                                "time", "tension", 
                                                NULL, NULL, 
                                                params, 0, nsamples, s);

        PositionResult(s, ((DisplayObject) s -> tsPlot[FORCE]) -> toplevel);
    }
    if (nsamples && (flags & PLOT_TIME) &&
        g_key_file_get_boolean(prefsFile, "output", "positionTimeSeries", &err)) {
        s -> tsPlot[DISPLACEMENT] =
                    CreateChartDisplay(s -> problem -> title,
                                       CHART_XY,
                                       s -> num_output_nodes,
                                       "time", "displacement",
                                       NULL, NULL, params, 0, nsamples, s);
        PositionResult(s, ((DisplayObject) s -> tsPlot[DISPLACEMENT]) -> toplevel);
    }

    if ((flags & PLOT_SPACE) &&
        g_key_file_get_boolean(prefsFile, "output", "positionSnapshots", &err)) {
        s -> snapPlot[DISPLACEMENT] = 
               CreateChartDisplay(s -> problem -> title, 
                                  CHART_ANIMATE, 1, 
                                  "x (m)", "z (m)", 
                                  NULL, NULL, NULL, 0, 
                                  s -> problem -> num_nodes, s);
        PositionResult(s, ((DisplayObject) s -> snapPlot[DISPLACEMENT]) -> toplevel);
    }
        
    if ((flags & PLOT_SPACE) &&
        g_key_file_get_boolean(prefsFile, "output", "tensionSnapshots", &err)) {
        
        s -> snapPlot[FORCE] = 
               CreateChartDisplay(s -> problem -> title, 
                                  CHART_ANIMATE, 1, 
                                  "s (m)", "tension (lbs)", 
                                  NULL, NULL, NULL, 0, 
                                  s -> problem -> num_nodes, s);
        PositionResult(s, ((DisplayObject) s -> snapPlot[FORCE]) -> toplevel);
    }

    if (s -> playback)
        gtk_window_present(GTK_WINDOW(((Playback *) s -> playback) -> toplevel));
}

void 
ControlDialogInitialize(Analysis *a, void *controls, int flags)
{
    GtkWidget **ctlparam = (GtkWidget **) controls;

    aptr = a;

    if (controls) {
        ctlparam = active_ctlparam = controls;
    }
    else {
        ctlparam = active_ctlparam;
    }

    if (baseline) {
        free(baseline);
    }

    baseline = (Analysis *) malloc(sizeof(Analysis));
    memcpy(baseline, a, sizeof(Analysis));

    if (a -> solution -> plotProgress) {
        ControlMakePlots(a -> solution, flags);
    }

    gtk_window_present(GTK_WINDOW(ctlparam[0]));
    FillControl(aptr, ctlparam);
}

void
ControlLabels(char *l1, char *l2, char *l3, char *l4, void *controls)
{
    // GtkWidget **ctlparam = (GtkWidget **) controls;
}

void 
ControlInfo(double tm, int it, double err, void *controls)
{
    GtkWidget    **ctlparam;
    static double  prev_tm = -1.0;

    if (controls) 
        ctlparam = (GtkWidget **) controls;
    else
        ctlparam = active_ctlparam;

    if (tm != prev_tm) {
        SetNumericText(ctlparam[CONTROL_TIME], "%g", tm, FALSE);
        prev_tm = tm;
    }
    
    SetNumericText(ctlparam[CONTROL_ITERATION], "%.0f", (double) it, FALSE);
    SetNumericText(ctlparam[CONTROL_ERROR], "%g", err, FALSE);
}

void
ControlAuxError(double err, void *controls)
{
    // GtkWidget **ctlparam = (GtkWidget **) controls;
}

void
ControlMessage(char *msg, void *controls)
{
    GtkWidget **ctlparam; 
    GtkWidget *dialog;

    if (controls) 
        ctlparam = (GtkWidget **) controls;
    else
        ctlparam = active_ctlparam;

    if (!ctlparam) {
        dialog = gtk_message_dialog_new(NULL, GTK_DIALOG_MODAL, GTK_MESSAGE_ERROR, GTK_BUTTONS_OK, msg);
        gtk_dialog_run(GTK_DIALOG(dialog));
        gtk_widget_destroy(dialog);
        return;
    }

    gtk_text_buffer_insert_at_cursor(GTK_TEXT_BUFFER(gtk_text_view_get_buffer(GTK_TEXT_VIEW(ctlparam[CONTROL_MESSAGE]))), msg, -1);
    gtk_text_buffer_insert_at_cursor(GTK_TEXT_BUFFER(gtk_text_view_get_buffer(GTK_TEXT_VIEW(ctlparam[CONTROL_MESSAGE]))), "\n", 1);
}

int
ControlProcessEvents(Solution *s)
{
    while(gtk_events_pending()) {
        gtk_main_iteration();
    }
    return 0;
}

void
ControlPlotTimeSeries(Result *res, Problem *p, Environment *e)
{
    int i;
    double  *t;

    t = (double *) malloc(sizeof(double) * res -> nsamples);
    for (i = 0 ; i < res -> nsamples ; i++) {
        t[i] = (i - 1)*res -> sample_dt;
    }

    for (i = 1 ; i < MAXOUTPUT ; i++) {
        if (p -> solution -> tsPlot[i]) {
            
            switch(i) {
            case FORCE:
                PlotChartDisplay((DisplayObject) p -> solution -> tsPlot[i],
                                 t, res -> force_t[0]);
                break;
            case DISPLACEMENT:
                PlotChartDisplay((DisplayObject) p -> solution -> tsPlot[i],
                                 t, res -> y_t);
                break;
            }
        }
    }

    free(t);
}

void
ControlPlotTime(double t, Problem *p, Environment *e)
{
    Node            *node; 
    double           y[20];
    int              i, j, k;

    node = p -> node;

    for (i = 1 ; i < MAXOUTPUT ; i++) {
        if (p -> solution -> tsPlot[i] 
            && !(((DisplayObject) (p -> solution -> tsPlot[i])) -> destroy)) {
            for (j = 0 ; j < p -> solution -> num_output_nodes ; j++) {
                k = p -> solution -> output_nodes[j];
                switch(i) {

                case DISPLACEMENT:
                    y[j] = node[k] -> x;
                    break;
                case FORCE:
                    y[j] = Tension(node[k] -> Y[1], node[k] -> material)/4.4482216;
                    break;
                default:
                    break;
                }
            }
            UpdateChartDisplay(p -> solution -> tsPlot[i], t, y, 0);
        }
    }
}

void
ControlPlotSnaps(Problem *p, Environment *e)
{
    Node            *node; 
    int              nn; 
    double          *x;
    double         **y;
    int              i, j;

    node = p -> node;
    nn   = p -> num_nodes;

    for (i = 1 ; i < MAXOUTPUT ; i++) {
        if (p -> solution -> snapPlot[i] 
            && !(((DisplayObject) (p -> solution -> snapPlot[i])) -> destroy)) {
            
            x = ChartX(p -> solution -> snapPlot[i]);
            y = ChartY(p -> solution -> snapPlot[i]);
            for (j = 1 ; j <= nn ; j++) {
                switch(i) {

                case DISPLACEMENT:
                    x[j-1] = node[j] -> y;
                    y[0][j-1] = node[j] -> x;
                    break;
                case FORCE:
                    x[j-1] = node[j] -> s;
                    y[0][j-1] = Tension(node[j] -> Y[1], node[j] -> material)/4.4482216;
                    break;
                default:
                    x[j-1] = node[j] -> s;
                    y[0][j-1] = 0;
                }
            }

            UpdateChartAnimation(p -> solution -> snapPlot[i], NULL);
        }
    }
}

void
TableQuit (GtkWidget *widget, gpointer data)
{
    Solution *s = (Solution *) data;

    if (s && s -> table) {
        gtk_widget_destroy(GTK_WIDGET(s -> table));
        s -> table = NULL;
    }
    else {
        gtk_widget_destroy(widget);
    }
}

gboolean
TableDelete(GtkWidget *w, GdkEvent *event, gpointer data)
{
    TableQuit(w, data);
    return TRUE;
}


static void
ColorSWL(GtkTreeViewColumn *column, GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter, gpointer data)
{
    float   swl;
    gchar   buf[20];

    gtk_tree_model_get(model, iter, (int) data, &swl, -1);

    g_snprintf(buf, sizeof(buf), "%.1f%%", swl*100.0);
    g_object_set(renderer, "text", buf, 
                           "foreground", swl >= 1.0 ? "red" : "black", NULL);
}

double
sf(Node n, double level)
{
    if (level)
        return NodeTension(n)/level;
    else
        return 0.0;
}

static void
SaveTable(gpointer data, guint action, GtkWidget *w)
{
    GtkTreeView *tree = (GtkTreeView *) data;
    GtkTreeModel *model;
    GtkTreePath *path;
    GtkTreeIter  iter1, iter2;
    gboolean     valid1, valid2;
    int          i;
    char        *ptr;
    float        x; 
    char        *filename;
    FILE        *fp;
    GtkWidget   *headers;
    gboolean     incl_headers;
    GList       *cols;
    int          ncols;
       
    headers = gtk_check_button_new_with_label("include headers");

    filename = FilenameFromDialog("save table", gtk_widget_get_toplevel(w),
                                  GTK_FILE_CHOOSER_ACTION_SAVE, NULL, NULL, NULL, headers, &incl_headers);

    if (!filename) {
        return;
    }


    fp = fopen(filename, "w");
    if (!fp) {
        return;
    } 

    cols = gtk_tree_view_get_columns(tree);
    ncols = g_list_length(cols);

    if (incl_headers) {
        while(cols) {
            fprintf(fp, "%8.8s ", gtk_tree_view_column_get_title(GTK_TREE_VIEW_COLUMN(cols -> data))); 
            cols = cols -> next;
        } 
        fprintf(fp, "\n");
    }
    
    model = gtk_tree_view_get_model(tree);
 
    valid1 = gtk_tree_model_get_iter_first(model, &iter1);
    while(valid1) {
        gtk_tree_model_get(model, &iter1, 0, &ptr, -1);
        fprintf(fp, "%8.8s ", ptr);
        for (i = 1 ; i < ncols ; i++) {
            gtk_tree_model_get(model, &iter1, i, &x, -1);
            fprintf(fp,"%8.2f ", x);
        }
        fprintf(fp, "\n");

        path = gtk_tree_model_get_path(model, &iter1);
        
        if (gtk_tree_view_row_expanded(tree, path)) {
            valid2 = gtk_tree_model_iter_children(model, &iter2, &iter1);
            while(valid2) {
                gtk_tree_model_get(model, &iter2, 0, &ptr, -1);
                fprintf(fp, "%8.8s ", ptr);
                for (i = 1 ; i < ncols ; i++) {
                    gtk_tree_model_get(model, &iter2, i, &x, -1);
                    fprintf(fp, "%8.2f ", x);
                }
                fprintf(fp, "\n");

                valid2 = gtk_tree_model_iter_next(model, &iter2);
            }

        } 
        valid1 = gtk_tree_model_iter_next(model, &iter1);
    }

    fclose(fp);
    return;
}

static void
PrintTable(gpointer data, guint action, GtkWidget *w)
{
}

static GtkWidget *
tablemenu(GtkWidget *tree)
{
    GtkItemFactory *item_factory;
    // GtkAccelGroup *accel_group;
    static GtkItemFactoryEntry menu_items[] = {
        { "/_File",         NULL,         NULL,    0, "<Branch>" },
        { "/File/_Print",    "<control>P", PrintTable,    1, "<StockItem>", GTK_STOCK_PRINT },
        { "/File/_Save",    "<control>S", SaveTable,    0, "<StockItem>", GTK_STOCK_SAVE },
    };

    item_factory = gtk_item_factory_new(GTK_TYPE_MENU_BAR, "<table>",
                                        NULL); // accel_group);


    gtk_item_factory_create_items (item_factory,
                                   sizeof(menu_items)/sizeof(menu_items[0]),
                                   menu_items, tree);

    return gtk_item_factory_get_widget (item_factory, "<table>");

}

void
ControlTabulateResults(Problem *p, Environment *e)
{
    Node         *node;
    int           nn;
    char          buff[8];
    double        T, maxT;
    double        Fx, Fy, Fz;
    Node          n, max_n;
    Segment      *seg;
    int          ns;
    int          i;
    Solution    *s;
    GtkWidget   *menubar;
    GtkTreeIter  iter1, iter2;
    GtkWidget   *window;
    GtkWidget   *scroller, *tree;
    GtkCellRenderer *renderer;
    GtkTreeViewColumn *column;
    GtkTreeStore *model;
    GtkWidget    *vbox;
    GError       *err = NULL;
    GType         types[] = {G_TYPE_STRING, G_TYPE_FLOAT, G_TYPE_FLOAT, G_TYPE_FLOAT, G_TYPE_FLOAT, G_TYPE_FLOAT, G_TYPE_FLOAT, G_TYPE_FLOAT, G_TYPE_FLOAT, G_TYPE_FLOAT};
    char *cols[] = {"node", "s", "z", "depth", "x", "Fx", "Fz", "tension (lbs)", "% swl", "% yield"};
#define NCOLS 10

    s = p -> solution;

    if (!g_key_file_get_boolean(prefsFile, "output", "tabulate", &err)) {
        s -> table = NULL;
        return;
    }

    node = p -> node;
    nn   = p -> num_nodes;

    window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(window), p -> title);
    gtk_widget_set_usize(GTK_WIDGET(window), 700, 500);

    gtk_signal_connect(GTK_OBJECT (window), "delete-event",
                       GTK_SIGNAL_FUNC(TableDelete), s);

    model = gtk_tree_store_newv(NCOLS, types);

    scroller = gtk_scrolled_window_new(NULL, NULL);

    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
                                   GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);

    tree = gtk_tree_view_new_with_model(GTK_TREE_MODEL(model));


    for(i = 0 ; i < NCOLS ; i++) {
        renderer = gtk_cell_renderer_text_new();
        column = gtk_tree_view_column_new_with_attributes(cols[i], renderer, "text", i, NULL);
        gtk_tree_view_append_column(GTK_TREE_VIEW(tree), column);
        if (i == 8 || i == 9) {
            gtk_tree_view_column_set_cell_data_func(column, renderer, ColorSWL, (gpointer) i, NULL);
        }
    }

    gtk_tree_view_set_rules_hint(GTK_TREE_VIEW(tree), TRUE);

    vbox = gtk_vbox_new(FALSE, 1);
    gtk_container_add(GTK_CONTAINER(window), vbox);

    menubar = tablemenu(tree);
    gtk_box_pack_start(GTK_BOX(vbox), menubar, FALSE, TRUE, 0);

    gtk_container_add(GTK_CONTAINER(scroller), tree);
    gtk_box_pack_start(GTK_BOX(vbox), scroller, TRUE, TRUE, 0);

    gtk_widget_show_all(window);

    if (p -> terminal[2] -> buoy) {
        T = NodeTension(node[nn]);
        TensionToGlobal(node[nn], TRUE, &Fx, &Fy, &Fz);
        gtk_tree_store_append(model, &iter1, NULL);
        gtk_tree_store_set(model, &iter1, 
                           0, p -> terminal[2] -> buoy -> name, 
                           1, node[nn] -> s,
                           2, node[nn] -> x, 
                           3, e -> depth - node[nn] -> x,
                           4, node[nn] -> y, 
                           5, Fy/4.4482216,
                           6, Fx/4.4482216,
                           7, T/4.4482216, 
                           8, sf(node[nn], node[nn] -> segment -> material -> swl),
                           9, sf(node[nn], node[nn] -> segment -> material -> yield),
                           -1);
    }

    seg = BuildSegmentArray(p, &ns, 0);
    for (i = ns ; i >= 1 ; i--) {
        maxT = 0;
        max_n = 0;
        n = seg[i] -> last;
        while (n && n -> segment == seg[i]) {
            T = NodeTension(n);
            if (T > maxT) {     
                maxT = T;
                max_n = n;
            }
            n = n -> prev;
        }

        n = seg[i] -> last;

        TensionToGlobal(max_n, TRUE, &Fx, &Fy, &Fz);

        gtk_tree_store_append(model, &iter1, NULL);
        gtk_tree_store_set(model, &iter1, 
                           0, seg[i] -> material -> name, 
                           1, n -> s, 
                           2, n -> x, 
                           3, e -> depth - n -> x,
                           4, n -> y, 
                           5, Fy/4.4482216,
                           6, Fx/4.4482216,
                           7, maxT/4.4482216, 
                           8, sf(max_n, n -> segment -> material -> swl),
                           9, sf(max_n, n -> segment -> material -> yield),
                           -1);

        while(n && n -> segment == seg[i]) {
            T = NodeTension(n); 
            TensionToGlobal(n, TRUE, &Fx, &Fy, &Fz);
            sprintf(buff, "%d", n -> number);
            gtk_tree_store_append(model, &iter2, &iter1);
            gtk_tree_store_set(model, &iter2, 
                               0, buff, 
                               1, n -> s,
                               2, n -> x, 
                               3, e -> depth - n -> x,
                               4, n -> y, 
                               5, Fy/4.4482216,
                               6, Fx/4.4482216,
                               7, T/4.4482216, 
                               8, sf(n, n -> segment -> material -> swl),
                               9, sf(n, n -> segment -> material -> yield),
                               -1);
            n = n -> prev;
        }
    }

    if (p -> terminal[1] -> anchor) {
        T = NodeTension(node[1]); 
        TensionToGlobal(node[1], TRUE, &Fx, &Fy, &Fz);
        gtk_tree_store_append(model, &iter1, NULL);
        gtk_tree_store_set(model, &iter1, 
                           0, p -> terminal[1] -> anchor -> name, 
                           1, node[1] -> s,
                           2, node[1] -> x, 
                           3, e -> depth - node[1] -> x, 
                           4, node[1] -> y, 
                           5, Fy/4.4482216,
                           6, Fx/4.4482216,
                           7, T/4.4482216, 
                           8, sf(node[1], node[1] -> segment -> material -> swl),
                           9, sf(node[1], node[1] -> segment -> material -> yield),
                           -1);
    }

    s -> table = (void *) window;        
    PositionResult(s, window);
}

#if 0
void
ControlTabulateResults(Node *node, int nn, Environment *e)
{
    char *cols[] = {"status", "node", "segment", "tension", "z", "x"};
    char           buff[16];
    char         **y;
    int            i, j;
    int            max_T_i;
    double         max_T;
    double         T;
    int            show;
    int            nc = 6;

    y = (char **) malloc(sizeof(char *) * nn * nc);
    bzero(y, nn*nc);

    max_T = 0;
    
    for (i = nn ; i >= 1 ; i--) {
        show = 0;

        j = nn - i;

        T = NodeTension(node[i])/4.4482216;
        if (T > max_T) {
            max_T = T;
            max_T_i = i;
        }

        sprintf(buff, "%d\n", node[i] -> number); 
        y[j*nc + 1] = strdup(buff);

        if (node[i] == node[i] -> segment -> first ||
            node[i] == node[i] -> segment -> last) {
            y[j*nc + 2] = strdup(node[i] -> segment -> material -> name);
            show = 1;
    
            if (node[i] == node[i] -> segment -> last) {
                max_T = 0;
            }
            else {  // must be first (lowest) - time to check max
                if (y[(nn - max_T_i)*nc + 0]) {
                    free(y[(nn - max_T_i)*nc + 0]);
                }
                y[(nn - max_T_i)*nc + 0] = strdup("bold");
            }
        }
        else {
            y[j*nc + 2] = NULL;
        }
   
     
        sprintf(buff, "%11.5f\n", NodeTension(node[i])/4.4482216); 
        y[j*nc + 3] = strdup(buff);

        sprintf(buff, "%11.5f\n", e -> depth - node[i] -> x); 
        y[j*nc + 4] = strdup(buff);

        sprintf(buff, "%11.5f\n", node[i] -> y); 
        y[j*nc + 5] = strdup(buff);

        if (show && !y[i*nc + 0]) { // don't overwrite a bold statement
            y[j*nc + 0] = strdup("show");
        }
        else {
            y[j*nc + 0] = NULL;
        }
    }

    DisplaySheet(cols, 6, nn, y);
}
#endif

void
UpdateSpatialAnimation(DisplayObject obj, Result *res, int var, double t)
{
    int j;
    double *x;
    double **y;
    char    buff[20];

    x = ChartX(obj);
    y = ChartY(obj);

    for (j = 0 ; j < res -> npoints ; j++) {
        switch(var) {

        case DISPLACEMENT:
            x[j] = res -> x[j];
            y[0][j] = res -> y[j];
            break;
        case FORCE:
            x[j-1] = res -> s[j];
            y[0][j-1] = res -> force[0][j]/4.4482216;
            break;
        default:
            x[j-1] = res -> s[j];
            y[0][j-1] = 0;
        }
    }

    sprintf(buff, "t=%g", t);
    UpdateChartAnimation(obj, buff);
}

#ifndef WINDOWS

void
ExportResultMovie(DisplayObject obj)
{
    ChartMovie   mov;
    GdkPixbuf   *pix = NULL;
    char        *filename;
    static char *prev_file = NULL;
    int          i, j;
    Solution    *s;
    int          spatial = 0;
    int          var = 0;
    Result      *res;
    double      *x;
    double     **y;
   
    s = (Solution *) obj -> data;
    res = (Result *) s -> results;

    if (res == NULL || !res -> dynamic) {
        return;
    }

    for (i = 1 ; i < MAXOUTPUT ; i++) {
        if (s -> snapPlot[i] == obj) {
            spatial = 1;
            var = i;
            break;
        }
        else if (s -> tsPlot[i] == obj) {
            spatial = 0;
            var = i;
            break;
        }
    }

    if (!spatial) {
	
        return;
    }

    filename = FilenameFromDialog("Export movie", obj -> toplevel,
                                  GTK_FILE_CHOOSER_ACTION_SAVE,
                                  &prev_file,
                                  "*.mpg", "movie files",
                                  NULL, NULL);
    if (!filename) {
        return;
    }

    filename = AddFilenameExtension(filename, ".mpg");

    x = ChartX(obj);
    y = ChartY(obj);
    for (j = 0 ; j < res -> nsteps && !obj -> destroy ; j++) {
        ReadResultSnapshot(res, j, 1, 1, 0);
        // printf("updating plot\n"); fflush(stdout);
        UpdateSpatialAnimation(obj, res, var, j*res -> snap_dt);

        // printf("getting pixbuf\n"); fflush(stdout);
        pix = ChartPixbuf(obj, pix);

        if (j == 0) {
            // printf("initializing movie\n"); fflush(stdout);
            InitializeMovie(&mov, filename, 
                            gdk_pixbuf_get_width(pix), 
                            gdk_pixbuf_get_height(pix));
        }

        // printf("adding movie\n"); fflush(stdout);
        AddFrameToMovie(&mov, pix, 3);
    }

    CloseMovie(&mov);

    gdk_pixbuf_unref(pix);

    return;
}

#endif

gboolean
ControlUpdatePlayback(gpointer data)
{
    Solution *s = (Solution *) data;
    Playback *p = (Playback *) s -> playback;
    Result   *res;
    int       i;
    char      buff[32];
    double    t;

    res = (Result *) s -> results;

    p -> frameN ++;
    if (p -> frameN == res -> nsteps) {
        p -> frameN = 0;
    }

    t = p -> frameN * res -> snap_dt;
    snprintf(buff, 32, "%g", t);
    gtk_label_set_text(GTK_LABEL(p -> display), buff);

    ReadResultSnapshot(res, p -> frameN, 1, 1, 0);

    for (i = 1 ; i < MAXOUTPUT ; i++) {
        if (s -> snapPlot[i] 
            && !(((DisplayObject) (s -> snapPlot[i])) -> destroy)) {

            UpdateSpatialAnimation(s -> snapPlot[i], res, i, t);
        }
    }

    return TRUE;
}

void
frameFwd(GtkWidget *widget, gpointer data)
{
    Solution *s = (Solution *) data;
    Playback *p = (Playback *) s -> playback;

    if (p -> timer)
        return;

    ControlUpdatePlayback(s);
}

void
frameRev(GtkWidget *widget, gpointer data)
{
    Solution *s = (Solution *) data;
    Playback *p = (Playback *) s -> playback;

    if (p -> timer)
        return;

    p -> frameN -= 2;
    if (p -> frameN < -1)
        p -> frameN = ((Result *) (s -> results)) -> nsteps - 2;

    ControlUpdatePlayback(s);
}

extern void quitPlot(GtkWidget *, gpointer);

void
destroyTimer(int timer)
{
    GSource *src;

    if (timer) {
        src = g_main_context_find_source_by_id(NULL, timer);
        if (src && !g_source_is_destroyed(src)) {
            g_source_destroy(src);
            while(!g_source_is_destroyed(src)) { } ;
        }
    }
}

static void
closePlots(Solution *s)
{
    int i;

    if (s -> table) {
        TableQuit(s -> table, s);
    }

    for (i = 1 ; i < MAXOUTPUT ; i++) {
        if (s -> tsPlot[i] && ((DisplayObject) s -> tsPlot[i]) -> destroy == 0) {
            quitPlot(NULL, s -> tsPlot[i]);
            // gtk_widget_destroy(GTK_WIDGET(((DisplayObject) (s -> tsPlot[i])) -> toplevel));
        }
        if (s -> snapPlot[i] && ((DisplayObject) s -> snapPlot[i]) -> destroy == 0) {
            quitPlot(NULL, s -> snapPlot[i]);
            // gtk_widget_destroy(GTK_WIDGET(((DisplayObject) (s -> snapPlot[i])) -> toplevel));
        }
    }

}

void
closePlayback(GtkWidget *widget, gpointer data)
{
    Solution *s = (Solution *) data;
    Playback *p = (Playback *) s -> playback;

    fprintf(stderr,"closePlayback\n");

    if (p -> timer) {
        destroyTimer(p -> timer);
        p -> timer = 0;
    }

    closePlots(s); 

    // scenario 1 - we have already closed the solution control dialog
    // or it never existed so we need to do all the cleanup here
    if (s -> controls == NULL) {
        if (s -> out_name) {
            fprintf(stderr, "playback closing - deleting output file\n");
            unlink(s -> out_name);
        }
        fprintf(stderr, "playback closing - freeing solution\n");
        FreeSolution(s);
    }
    // scenario 2 - solution control still exists so we might still
    // need to save, etc. -> only flag the playback as closed
    // so that the solution control close function knows it is 
    // responsible for cleanup
    else {
        s -> playback = NULL;
    }

    g_object_unref(widget);
    gtk_widget_destroy(gtk_widget_get_toplevel(widget));

}

gboolean
deletePlayback(GtkWidget *widget, GdkEvent *event, gpointer data)
{
    fprintf(stderr,"deletePlayback\n");
    closePlayback(widget, data);
    return TRUE;
}


void
playPlayback(GtkWidget *widget, gpointer data)
{
    Solution *s = (Solution *) data;
    Playback *p = (Playback *) s -> playback;

    if (gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(widget))) {
        printf("starting playback\n");
        p -> timer = g_timeout_add(p -> interval, ControlUpdatePlayback, s);
    }
    else {
        if (p -> timer) {
            destroyTimer(p -> timer);
            p -> timer = 0;
        }
    }
}

void
speedIncr(GtkWidget *widget, gpointer data) 
{
    Solution *s = (Solution *) data;
    Playback *p = (Playback *) s -> playback;

    if (p -> interval == 10)
        return;

    if (p -> interval <= 100)
        p -> interval -= 10; 
    else
        p -> interval -= 100;

    if (p -> timer) {
        destroyTimer(p -> timer);
        p -> timer = g_timeout_add(p -> interval, ControlUpdatePlayback, s);
    }
}

void
speedDecr(GtkWidget *widget, gpointer data)
{
    Solution *s = (Solution *) data;
    Playback *p = (Playback *) s -> playback;

    if (p -> interval == 1000)
        return;

    if (p -> interval >= 100)
        p -> interval += 100;
    else
        p -> interval += 10;

    if (p -> timer) {
        destroyTimer(p -> timer);
        p -> timer = g_timeout_add(p -> interval, ControlUpdatePlayback, s);
    }
}



GtkWidget *
BuildPlayback(Solution *s)
{
    GtkWidget *win;
    char       buffer[128];
    GtkWidget *toolbar;
    Playback  *p;
    GtkWidget *hbox;
    PangoFontDescription *desc; 
    Result    *results;

    results = (Result *) (s -> results);

    if (results -> dynamic)
        sprintf(buffer, "Playback control - %s", s -> problem -> title);
    else
        sprintf(buffer, "Display control - %s", s -> problem -> title);

    p = (Playback *) malloc(sizeof(Playback));
    p -> interval = 100;
    p -> frameN = 0;
    p -> timer = 0;

    s -> playback = p;
 
    win = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    p -> toplevel = win;
    if (results -> dynamic)
        gtk_window_set_default_size(GTK_WINDOW(win), 360, 32);
    else
        gtk_window_set_default_size(GTK_WINDOW(win), 240, 32);

    gtk_window_set_title(GTK_WINDOW(win), buffer);

    toolbar = gtk_toolbar_new();
    gtk_toolbar_set_style(GTK_TOOLBAR(toolbar), GTK_TOOLBAR_ICONS);
    hbox = gtk_hbox_new(FALSE, 5);

    toolbar_insert_stock(toolbar, GTK_STOCK_CLOSE,
                         GTK_SIGNAL_FUNC(closePlayback), s, 0, NULL);
    if (results -> dynamic) {
        toolbar_insert_stock(toolbar, GTK_STOCK_MEDIA_REWIND,
                         GTK_SIGNAL_FUNC(speedDecr), s, 0, NULL);
        toolbar_insert_stock(toolbar, GTK_STOCK_MEDIA_PREVIOUS,
                         GTK_SIGNAL_FUNC(frameRev), s, 0, NULL);
        toolbar_insert_stock(toolbar, GTK_STOCK_MEDIA_PLAY,
                         GTK_SIGNAL_FUNC(playPlayback), s, 1, NULL);
        toolbar_insert_stock(toolbar, GTK_STOCK_MEDIA_NEXT,
                         GTK_SIGNAL_FUNC(frameFwd), s, 0, NULL);
        toolbar_insert_stock(toolbar, GTK_STOCK_MEDIA_FORWARD,
                         GTK_SIGNAL_FUNC(speedIncr), s, 0, NULL);

        p -> display = gtk_label_new("0.0");

        desc = pango_font_description_from_string("Courier 24");
        gtk_widget_modify_font(p -> display, desc);
        pango_font_description_free(desc);

        gtk_widget_set_size_request(p -> display, 120, 32);
    }
    else {
        p -> display = gtk_label_new("static only");

        desc = pango_font_description_from_string("Courier 16");
        gtk_widget_modify_font(p -> display, desc);
        pango_font_description_free(desc);

        gtk_widget_set_size_request(p -> display, 160, 32);
    }

    gtk_signal_connect(GTK_OBJECT (win), "delete-event",
                       GTK_SIGNAL_FUNC(deletePlayback), s);

    gtk_box_pack_start(GTK_BOX(hbox), toolbar, TRUE, TRUE, 1);
    gtk_box_pack_start(GTK_BOX(hbox), p -> display, FALSE, FALSE, 1);

    gtk_container_add(GTK_CONTAINER(win), hbox);

    gtk_widget_show_all(win);

    return win;

}



GtkWidget **
BuildControl(Solution *s)
{
    int           i;
    GtkWidget   **ctlparam;
    static char  *rows[] = {"static", "outer", "dynamic"};
    static char  *cols[] = {"relaxation", "tolerance", "iterations"};
    static char  *dynamic_cols[] = {"time-step", "duration", "ramp-time"};
    static char  *dynamic_rows[] = {"dynamic"};
    static char  *runtime_cols[] = {"time", "iter", "error"};
    static char  *runtime_rows[] = {"progress"};
    EntryTable   *basic, *dynamic, *runtime;
    GtkWidget    *vbox, *hbox;
    GtkWidget    *ctlwin;
    GtkWidget    *scroller;
    GtkWidget    *quit, *pause, *update, *restore, *save, *matlab;

    ctlparam = (GtkWidget **) g_malloc(sizeof(GtkWidget *) * NUM_CONTROL);

    ctlwin = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(ctlwin), "WHOI Cable Solution Control");

    ctlparam[CONTROL_TOPLEVEL] = ctlwin;

    // the basic tab

    vbox = gtk_vbox_new(FALSE, 3);
    gtk_container_set_border_width(GTK_CONTAINER(vbox),5);
    gtk_container_add(GTK_CONTAINER(ctlwin), vbox);
    gtk_widget_show(vbox);


    // a 3x3 table for relaxation factors, iteration limits, and tolerances

    basic = MakeEntryTable(relaxation_entries, rows, cols, 9, 3, 3);
    gtk_box_pack_start(GTK_BOX(vbox), basic -> table, FALSE, FALSE, 1);
    for (i = 0 ; i < 9 ; i++) {
        ctlparam[i + CONTROL_STATIC_RELAX] = basic -> entries[i].w;
    }

    gtk_box_pack_start(GTK_BOX(vbox), gtk_hseparator_new(), FALSE, FALSE, 0);

    // a 1x3 table for dynamic timing

    dynamic = MakeEntryTable(dynamic_entries, dynamic_rows, dynamic_cols, 3, 1, 3);
    gtk_box_pack_start(GTK_BOX(vbox), dynamic -> table, FALSE, FALSE, 1);
    for (i = 0 ; i < 3 ; i++) {
        ctlparam[i + CONTROL_DT] = dynamic -> entries[i].w;
    }

    gtk_box_pack_start(GTK_BOX(vbox), gtk_hseparator_new(), FALSE, FALSE, 0);

    // a 1x3 table for runtime status

    runtime = MakeEntryTable(runtime_entries, runtime_rows, runtime_cols, 3, 1, 3);
    gtk_box_pack_start(GTK_BOX(vbox), runtime -> table, FALSE, FALSE, 1);
    for (i = 0 ; i < 3 ; i++) {
        ctlparam[i + CONTROL_TIME] = runtime -> entries[i].w;
    }


    gtk_box_pack_start(GTK_BOX(vbox), gtk_hseparator_new(), FALSE, FALSE, 0);

    scroller = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
                                   GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);

    ctlparam[CONTROL_MESSAGE] = gtk_text_view_new();
    gtk_container_add(GTK_CONTAINER(scroller), ctlparam[CONTROL_MESSAGE]);
    gtk_box_pack_start(GTK_BOX(vbox), scroller, FALSE, FALSE, 3);

    gtk_box_pack_start(GTK_BOX(vbox), gtk_hseparator_new(), FALSE, FALSE, 0);

    hbox = gtk_hbox_new(FALSE, 3);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 3);
    
    quit = gtk_button_new_with_label("quit");
    pause = gtk_toggle_button_new_with_label("pause");
    update = gtk_button_new_with_label("update");
    restore = gtk_button_new_with_label("restore");
    save = gtk_button_new_with_label("save");
    matlab = gtk_button_new_with_label("matlab");
    
    gtk_box_pack_start(GTK_BOX(hbox), quit, FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(hbox), pause, FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(hbox), update, FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(hbox), restore, FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(hbox), save, FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(hbox), matlab, FALSE, FALSE, 3);

    gtk_widget_show_all(vbox);

    gtk_signal_connect(GTK_OBJECT(quit), "clicked", GTK_SIGNAL_FUNC(QuitSolution), s);
    gtk_signal_connect(GTK_OBJECT(pause), "clicked", GTK_SIGNAL_FUNC(PauseSolution), s);
    gtk_signal_connect(GTK_OBJECT(update), "clicked", GTK_SIGNAL_FUNC(UpdateParameters), ctlparam);
    gtk_signal_connect(GTK_OBJECT(restore), "clicked", GTK_SIGNAL_FUNC(RestoreParameters), ctlparam);
    gtk_signal_connect(GTK_OBJECT(save),   "clicked", GTK_SIGNAL_FUNC(SaveSolution), s);
    gtk_signal_connect(GTK_OBJECT(matlab), "clicked", GTK_SIGNAL_FUNC(MatlabSolution), s);

    gtk_signal_connect(GTK_OBJECT (ctlwin), "delete-event",
                       GTK_SIGNAL_FUNC(DeleteSolution), s);

    active_ctlparam = ctlparam;
   
    return ctlparam;
}


