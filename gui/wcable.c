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
#include <gtk/gtk.h>
#include <glib.h>
#include <string.h>
#include <math.h>
#include <libgen.h>
#include <unistd.h>
#ifdef WINDOWS
#include <windows.h>
#endif
#ifdef NATIVE_OSX
#include <mach-o/dyld.h>
#include <libgen.h>
#include <glib/gstdio.h>
#endif
#include <sys/stat.h>
#include "svn_version.h"
#include "table.h"
#include "problem.h"
#include "localcontrol.h"
#include "entry.h"
#include "objsheets.h"
#include "dialog.h"
#include "text.h"
#include "combo.h"
#include "environment.h"
#include "analysis.h"
#include "segments.h"
#include "duplicates.h"
#include "results.h"
#include "compress.h"
#include "plot.h"
#include "localerror.h"
#ifdef NATIVE_OSX
#include <ige-mac-menu.h>
#endif


#define DEFAULT_PRECISION 3
#define DEFAULT_SPACE 8

typedef struct {
    GtkWidget   *vbox;
    GtkWidget   *x, *y, *z;
    GtkWidget   *select;
    GtkWidget   *type;
    GtkWidget   *release;
    GtkWidget   *safety;
    GtkWidget   *friction;
    GtkWidget   *buoy;
    GtkWidget   *anchor;
    int          object;
} TerminalTable;

TerminalTable terminal_1;
TerminalTable terminal_2;
#define TERMINAL_ANCHOR 1
#define TERMINAL_BUOY 2

extern ObjectSheet material_sheet;
extern ObjectSheet connector_sheet;
extern ObjectSheet anchor_sheet;
extern ObjectSheet buoy_sheet;



enum {
    LAYOUT = 0,
    MATERIALS = 1,
    CONNECTORS = 2,
    BUOYS = 3,
    ANCHORS = 4,
    NUM_SHEETS,
};

static char *sheet_title[] = {
    "layout",
    "materials",
    "attachments",
    "buoys",
    "anchors",
};

enum {
    NUMBER = 0,
    CONNECTOR = 1,
    MATERIAL,
    LENGTH,
    NODES,
    ATTACH,
    RUNNING_DEPTH,
    TOTAL_LENGTH,
    TOTAL_NODES,
    NUM_LAYOUT_COLUMNS,
};

static char *layout_headers[] = {
    "#", 
    "pin",
    "material",
    "length",
    "nodes",
    "",
    "DEPTH",
    "ALTIT",
    "NODES",
};


typedef struct {
    GtkWidget   *labelNumber;
    GtkWidget   *comboMaterial;
    GtkWidget   *checkPin;
    GtkWidget   *entryLength;
    GtkWidget   *spinNodes;
    GtkWidget   *labelDepth;
    GtkWidget   *labelAltit;
    GtkWidget   *labelNodes;
    GtkWidget   *expanderAttach;
    GtkTable    *tableAttachments; // maintain table and list in parallel
    GList       *attachmentList;
} guiSegment;

typedef struct {
    GtkWidget  *comboAttachment;
    GtkTable   *tableNodes;
    GList      *nodeList; // maintain table and list in parallel
    guiSegment *parent;
} guiAttachment;

typedef struct {
    GtkWidget *spinNode;
    GtkWidget *labelDepth;
    GtkWidget *labelAltit;
    guiAttachment *parent;
} guiAttachmentNode;

static GtkTable *layout_table; // maintain in parallel with segmentList;
static GList *segmentList; 
static int active_row = -1;

char *current_filename = NULL;

GtkWidget *toplevel;
GtkWidget *main_vbox;
GtkWidget *notebook;
GtkWidget *layout_hpane = NULL;

GKeyFile *prefsFile;
static GQueue *recentFiles;
static GtkWidget *recentWidgets[4];

static GtkWidget *solve_item;


static gboolean bufferTotalChanges = FALSE;

void
SetCurrentFile(char *name)
{
    if (current_filename) {
        free(current_filename);
    }
    current_filename = strdup(name);
    gtk_window_set_title(GTK_WINDOW(toplevel), basename(current_filename));
}

void
LoadPrefs(void)
{
    char *fname;
   
    prefsFile = g_key_file_new();
 
    fname = DatabasePath("prefs.ini");
    g_key_file_load_from_file(prefsFile, fname, G_KEY_FILE_KEEP_COMMENTS, NULL);
    if (!g_key_file_has_group(prefsFile, "output")) {
        g_key_file_set_string(prefsFile, "output", "toggles", "on");
    }
}

void SavePrefs(void)
{
    char    *fname;
    FILE    *fp;

    fname = DatabasePath("prefs.ini");
     
    fp = fopen(fname, "w");
    fprintf(fp, "%s", g_key_file_to_data(prefsFile, NULL, NULL));
    fclose(fp);
}

void
quit ()
{
  SavePrefs();
  gtk_main_quit();
}

static void
SetRowColor(gint row, GdkColor *color)
{
    guiSegment  *s;

    s = g_list_nth_data(segmentList, row - 1);
    
    if (s) {
        gtk_widget_modify_base(s -> labelNumber, GTK_STATE_NORMAL, color);
        gtk_widget_modify_base(s -> entryLength, GTK_STATE_NORMAL, color);
        gtk_widget_modify_base(s -> spinNodes, GTK_STATE_NORMAL, color);
    }
}

static gboolean
ActivateRow(GtkWidget *widget, GdkEventFocus *event, gpointer data)
{
    GdkColor     yellow, white;
    int          this_row;

    gdk_color_parse("yellow", &yellow);
    gdk_color_parse("white", &white);

    this_row = g_list_index(segmentList, data) + 1;  

    if (active_row > -1 && active_row != this_row) {
        SetRowColor(active_row, &white);
    }

    active_row = this_row;
    SetRowColor(this_row, &yellow);

    return FALSE;
}

static double
CurrentDepthEntry(void)
{
    double          depth;
    const gchar     *ptr;
    extern GtkWidget *envparam[];

    ptr = gtk_entry_get_text(GTK_ENTRY(envparam[0]));
    if (ptr && ptr[0]) {
        depth = atof(ptr);
        return atof(ptr);
    }
    else {
        depth = 0;
    }

    return depth;
}

static void ChangeAttachmentNode(GtkWidget *, gpointer);

static gboolean
ComputeTotals(GtkWidget *widget, GdkEventFocus *event, gpointer data)
{
    double       n;
    float        envDepth;
    double       length, thislen, y;
    int          nodecount;
    int          segment_nodes = 0;
    double       prev_seg_nodes; 
    char         buffer[32];
    char        *problem_type;
    char        *markup;
    gboolean     topdown;
    GtkAdjustment *adj;
    GList       *s, *a, *nodelist;
    guiSegment  *gs;
    guiAttachment *ga;
    guiAttachmentNode *gn;
    
     
    if (bufferTotalChanges)
        return FALSE;

    envDepth = (double) CurrentDepthEntry(); // atof(gtk_entry_get_text(GTK_ENTRY(envparam[0])));

    problem_type = AnalysisProblemType();
    if (strcmp(problem_type, "surface") == 0 || strcmp(problem_type, "towing"))
        topdown = TRUE;
    else
        topdown = FALSE;

    // make a bottom-up pass

    length = 0.0;
    nodecount = 0;
  
    s = g_list_last(segmentList);
 
    while(s) {
        gs = (guiSegment *) s -> data;
        
        markup = g_markup_printf_escaped("<span weight=\"bold\">%03d</span>", g_list_index(segmentList, gs) + 1);
        gtk_label_set_markup(GTK_LABEL(gs -> labelNumber), markup);
        g_free(markup);

        thislen = atof(gtk_entry_get_text(GTK_ENTRY(gs -> entryLength)));
        length += thislen;
        segment_nodes = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(gs -> spinNodes));

        nodecount += segment_nodes;

        snprintf(buffer, 10, "%6.1f", length);
        gtk_label_set_text(GTK_LABEL(gs -> labelAltit), buffer);
/*
        if (envDepth) {
            snprintf(buffer, 10, "%6.1f", envDepth - length);
            gtk_label_set_text(GTK_LABEL(gs -> labelDepth), buffer);
        }
        else {
            gtk_label_set_text(GTK_LABEL(gs -> labelDepth), "0.0");
        }
*/
        snprintf(buffer, 10, "%6d", nodecount);
        gtk_label_set_text(GTK_LABEL(gs -> labelNodes), buffer);

        s = s -> prev;
    }

    // make a top-down pass

    s = g_list_first(segmentList);

    if (topdown || !envDepth) { 
        length = 0.0;
    }
    else {
        length = envDepth - length;
    }
 
    while(s) {
        gs = (guiSegment *) s -> data;

        snprintf(buffer, 10, "%6.1f", length);
        gtk_label_set_text(GTK_LABEL(gs -> labelDepth), buffer);
        thislen = atof(gtk_entry_get_text(GTK_ENTRY(gs -> entryLength)));
        length += thislen;

        segment_nodes = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(gs -> spinNodes));

        // depth and altit are set now so ChangeAttachmentNode works

        a = g_list_first(gs -> attachmentList);
        while(a) {
            ga = (guiAttachment *) a -> data;
            nodelist = g_list_first(ga -> nodeList);
            while(nodelist) {
                gn = (guiAttachmentNode *) nodelist -> data;
                n = gtk_spin_button_get_value(GTK_SPIN_BUTTON(gn -> spinNode));
                adj = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(gn -> spinNode));
                prev_seg_nodes = adj -> upper;
                if (prev_seg_nodes != segment_nodes) {
                    adj -> upper = segment_nodes;
                    gtk_adjustment_changed(adj);
                    y = (n - 1.0)*thislen / (prev_seg_nodes - 1.0);
                    n = rint(y*(segment_nodes - 1.0)/thislen + 1); 
                    gtk_spin_button_set_value(GTK_SPIN_BUTTON(gn -> spinNode), 
                                              (double) n);
                }

                ChangeAttachmentNode(gn -> spinNode, gn);

                nodelist = nodelist -> next;
            }
            a = a -> next;
        }

        s = s -> next;
    }

    return FALSE;
}

gboolean
ComputeTotalsForDepth(GtkWidget *widget, GdkEventFocus *event, gpointer data)
{
    return ComputeTotals(NULL, NULL, NULL);
}

static void
ComputeNodes(GtkWidget *widget, gpointer data)
{
    ComputeTotals(widget, NULL, data);
    return;
}

guiSegment *BuildSegmentRow(int);

static guiSegment *
AddNewSegment(int row, gboolean run_totals)
{
    guiSegment *gs;

    gtk_table_insert_row(layout_table, row);
    gs = BuildSegmentRow(row);

    if (run_totals) {
        ComputeTotals(NULL, NULL, NULL);
    }

    if (row <= active_row) {
        ActivateRow(NULL, NULL, g_list_nth_data(segmentList, active_row));
    }

    return gs;
}

static void
AddNewSegmentActive(GtkWidget *widget, gpointer data)
{
    AddNewSegment(active_row + (int) data, TRUE);
};


GtkWidget *
ImageButton(char *stock)
{
    GtkWidget   *image;
    GtkWidget   *button;

    button = gtk_button_new();
    image = gtk_image_new_from_stock(stock, GTK_ICON_SIZE_SMALL_TOOLBAR);
    gtk_container_add(GTK_CONTAINER(button), image);
    gtk_widget_show(button);
    gtk_widget_show(image);

    return button;
}

static void
FreeSegmentMemory(Segment seg)
{
    int j;

    for (j = 1 ; j <= seg -> num_attach ; j++) {
        if (seg -> attach[j].num_nodes) {
            seg -> attach[j].nodes ++; 
            free(seg -> attach[j].nodes);
        }
    }
    if (seg -> num_attach) {
        seg -> attach ++; free(seg -> attach);
    }
}

static guiSegment *
DecodeSegmentBox(int row, Segment seg, int fill_attachments)
{
    gchar       *str;
    GtkTreeIter iter;
    guiSegment *s;
    guiAttachment *ga;
    guiAttachmentNode *gn;
    GList *alist, *nlist;
    int    j, k;

    s = g_list_nth_data(segmentList, row - 1);
 
    seg -> length = atof(gtk_entry_get_text(GTK_ENTRY(s -> entryLength)));
    seg -> dist[1].percent = 1.0;
    seg -> num_nodes = seg -> dist[1].nodes = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(s -> spinNodes));
    seg -> num_dist = 1;
    if (gtk_combo_box_get_active_iter(GTK_COMBO_BOX(s -> comboMaterial), &iter)) {
        gtk_tree_model_get(GTK_TREE_MODEL(material_sheet.model), &iter, 0, &str, -1);
        seg -> material = (Material) str;
    }
    else {
        seg -> material = NULL;
    }

    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(s -> checkPin))) {
        seg -> connector = NULL; // (void *) 1;
        seg -> connection = Pinned;
    }
    else {
        seg -> connector = NULL;
    }

    if (fill_attachments) {
        seg -> num_attach = 0;
        // g_list_length(s -> attachmentList);
        seg -> attach = NULL; 
        // (attachment *) malloc(sizeof(attachment)*seg -> num_attach); seg -> attach --;
        for (j = 1, alist = g_list_first(s -> attachmentList); alist ; j++) {
            ga = (guiAttachment *) alist -> data;
            if (gtk_combo_box_get_active_iter(GTK_COMBO_BOX(ga -> comboAttachment), &iter)) {
                gtk_tree_model_get(GTK_TREE_MODEL(connector_sheet.model), &iter, 0, &str, -1);
            }
            else {
                str = NULL;
            }

            if (g_list_length(ga -> nodeList) && str) {
                if (seg -> attach) {
                    seg -> attach ++;
                }
                seg -> num_attach ++;
                seg -> attach = (attachment *) realloc(seg -> attach, sizeof(attachment)*seg -> num_attach); seg -> attach --;
                seg -> attach[seg -> num_attach].num_nodes = g_list_length(ga -> nodeList);
                seg -> attach[seg -> num_attach].nodes = (int *) malloc(sizeof(int) * seg -> attach[seg -> num_attach].num_nodes); seg -> attach[seg -> num_attach].nodes --; // unit offset
                seg -> attach[seg -> num_attach].object = (Connector) str;
                for (k = 1, nlist = g_list_first(ga -> nodeList) ; nlist ; k++) {
                    gn = (guiAttachmentNode *) nlist -> data;
                    seg -> attach[seg -> num_attach].nodes[k] = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(gn -> spinNode));
                    nlist = nlist -> next;
                }
            }    
            alist = alist -> next;
        }
    }

    return s;
}


static void
MaterialChanged(GtkWidget *widget, gpointer data)
{
    gchar       *str;
    GtkTreeIter iter;
   
    if (gtk_combo_box_get_active_iter(GTK_COMBO_BOX(widget), &iter)) {
        gtk_tree_model_get(GTK_TREE_MODEL(material_sheet.model), &iter, 5, &str, -1);
        if (atof(str) > 0) {
            gtk_entry_set_text(GTK_ENTRY(data), str);
        }
    }
}

static void
ComboSensitivity(GtkCellLayout *layout, GtkCellRenderer *cell,
                 GtkTreeModel *model, GtkTreeIter *iter, gpointer data)
{
    gboolean sensitive;

     sensitive = !gtk_tree_model_iter_has_child(model, iter);
     g_object_set(cell, "sensitive", sensitive, NULL);
}

static void
DeleteAttachmentNode(GtkWidget *widget, gpointer data)
{
    guiAttachmentNode *n = (guiAttachmentNode *) data;

    gtk_table_delete_row(n -> parent -> tableNodes, g_list_index(n -> parent -> nodeList, n) + 1);
    n -> parent -> nodeList = g_list_remove(n -> parent -> nodeList, n);
    free(n);
}

static void 
DeleteAttachment(GtkWidget *widget, gpointer data)
{
    guiAttachment *a = (guiAttachment *) data;
    GList        *n;
   
    n = g_list_first(a -> nodeList);
    while(n) {
        DeleteAttachmentNode(NULL, n -> data);
        n = n -> next;
    } 
    gtk_table_delete_row(a -> parent -> tableAttachments, g_list_index(a -> parent -> attachmentList, a) + 1);
    a -> parent -> attachmentList = g_list_remove(a -> parent -> attachmentList, a);
    if (g_list_length(a -> parent -> attachmentList) == 0) {
        gtk_widget_modify_base(a -> parent -> expanderAttach, GTK_STATE_NORMAL, NULL);
        gtk_widget_modify_text(a -> parent -> expanderAttach, GTK_STATE_NORMAL, NULL);
        gtk_expander_set_label(GTK_EXPANDER(a -> parent -> expanderAttach), "_Attach");
    }

    free(a);
}

static void 
DeleteSegment(int row)
{
    guiSegment *s;
    GList      *a;

    s = (guiSegment *) g_list_nth_data(segmentList, row - 1);

    a = g_list_first(s -> attachmentList);
    while(a) {
        DeleteAttachment(NULL, a -> data);
        a = a -> next;
    } 

    gtk_table_delete_row(GTK_TABLE(layout_table), row);
    segmentList = g_list_remove(segmentList, s);
    free(s);

    if (row == active_row) {
        active_row = -1;
        if (layout_table -> nrows == 1) {
            s = BuildSegmentRow(1);
            row = 1;
        }
        else if (row > 1) {
            row = row - 1;
        }

        ActivateRow(NULL, NULL, g_list_nth_data(segmentList, row-1));
    }

    ComputeTotals(NULL, NULL, NULL);
}

static GtkTargetEntry targetEntry[1] = {{"segDrag", GTK_TARGET_SAME_APP, 1}};

static GtkWidget *dragSource = NULL;

static void
SourceBeginDrag(GtkWidget *w, GdkDragContext *dc, gpointer data)
{
    dragSource = w;
}

static void FillSegmentBox(guiSegment *, Segment, gboolean);

static void
DestEndDrag(GtkWidget *w, GdkDragContext *dc, 
            gint x, gint y, guint t, gpointer data)
{
    int dest_row, source_row; 
    guiSegment *sourceSeg, *destSeg;
    struct segment seg;

    if (!dragSource) {
        return;
    }

    sourceSeg = g_object_get_data(G_OBJECT(dragSource), "segment");
    destSeg   = g_object_get_data(G_OBJECT(w), "segment");
    source_row = g_list_index(segmentList, sourceSeg) + 1;
    dest_row = g_list_index(segmentList, destSeg) + 1;

    DecodeSegmentBox(source_row, &seg, TRUE);
    destSeg = AddNewSegment(dest_row + 1, FALSE);
    FillSegmentBox(destSeg, &seg, FALSE);
    FreeSegmentMemory(&seg);

    source_row = g_list_index(segmentList, sourceSeg) + 1;
    DeleteSegment(source_row);

    if (dest_row > source_row)
        ActivateRow(NULL, NULL, g_list_nth_data(segmentList, dest_row - 1));
    else if (dest_row < source_row)
        ActivateRow(NULL, NULL, g_list_nth_data(segmentList, dest_row));
  
    dragSource = NULL;
}

static void
ChangeAttachmentNode(GtkWidget *widget, gpointer data)
{
    guiAttachmentNode *n = (guiAttachmentNode *) data;
    char       buffer[16];
    float      node;
    double     A_top, A, D_top, D;
    struct segment seg;
    int        row;
    char      *problem_type;
    gboolean   topdown;

    problem_type = AnalysisProblemType();
    if (strcmp(problem_type, "surface") == 0 || strcmp(problem_type, "towing"))
        topdown = TRUE;

    row = g_list_index(segmentList, n -> parent -> parent);
    DecodeSegmentBox(row + 1, &seg, FALSE );
 
    node = gtk_spin_button_get_value(GTK_SPIN_BUTTON(widget));

    A_top = atof(gtk_label_get_text(GTK_LABEL(n -> parent -> parent -> labelAltit)));
    D_top = atof(gtk_label_get_text(GTK_LABEL(n -> parent -> parent -> labelDepth)));

    if (D_top >= 0) {
        D = D_top + (seg.num_nodes - node)/(seg.num_nodes - 1)*seg.length;
     
        sprintf(buffer, "%6.1f\n", D);
        gtk_label_set_text(GTK_LABEL(n -> labelDepth), buffer);
    }
    else {
        gtk_label_set_text(GTK_LABEL(n -> labelDepth), "-");
    }

    A = (node - 1.0)*seg.length/(seg.num_nodes - 1) + (A_top - seg.length);
    sprintf(buffer, "%6.1f\n", A);
    gtk_label_set_text(GTK_LABEL(n -> labelAltit), buffer);
}

static guiAttachmentNode *
AddAttachmentNode(GtkWidget *widget, gpointer data)
{
    guiAttachmentNode *n;
    guiAttachment *a = (guiAttachment *) data;
    GtkWidget     *delete;
    int            num_rows;
    struct segment seg;
    GtkAdjustment *adj;

    // Do this decode first before we've added this node to the list.
    // Otherwise the decoder will try to fill out the segment's 
    // attachment data for the very node we're trying to add
    DecodeSegmentBox(g_list_index(segmentList, a -> parent) + 1, &seg, FALSE);

    n = (guiAttachmentNode *) malloc(sizeof(guiAttachmentNode));

    n -> parent = a;
    a -> nodeList = g_list_append(a -> nodeList, n);

    adj = (GtkAdjustment *) gtk_adjustment_new(1, 1, seg.num_nodes, 1, 10, 0); 

    delete = ImageButton(GTK_STOCK_DELETE);

    n -> spinNode  = gtk_spin_button_new(adj, 1.0, 0);
    n -> labelDepth = gtk_label_new("-");
    n -> labelAltit = gtk_label_new("-");
    
    num_rows = a -> tableNodes -> nrows;
    gtk_table_resize(a -> tableNodes, num_rows + 1, 4);

    gtk_table_attach(a -> tableNodes, delete, 
                     0, 1, num_rows, num_rows + 1, 0, 0, 2, 0);
    gtk_table_attach(a -> tableNodes, n -> spinNode,   
                     1, 2, num_rows, num_rows + 1, 0, 0, 2, 0);
    gtk_table_attach(a -> tableNodes, n -> labelDepth,  
                     2, 3, num_rows, num_rows + 1, 0, 0, 2, 0);
    gtk_table_attach(a -> tableNodes, n -> labelAltit,  
                     3, 4, num_rows, num_rows + 1, 0, 0, 2, 0);

    gtk_signal_connect(GTK_OBJECT(delete), "clicked",
                       GTK_SIGNAL_FUNC(DeleteAttachmentNode), n);

    gtk_signal_connect(GTK_OBJECT(n -> spinNode), "value-changed",
                       GTK_SIGNAL_FUNC(ChangeAttachmentNode), (gpointer) n);

    gtk_widget_show(delete);
    gtk_widget_show(n -> spinNode);
    gtk_widget_show(n -> labelDepth);
    gtk_widget_show(n -> labelAltit);

    ChangeAttachmentNode(n -> spinNode, n);
    return n;
}

static guiAttachment *
AddAttachment(GtkWidget *widget, gpointer data)
{
    guiSegment   *s = (guiSegment *) data;
    GtkCellRenderer *renderer;
    GtkWidget     *delete;
    GtkWidget     *nodes;
    GtkWidget     *nodeAdd;
    int            num_rows;
    guiAttachment *a;
    GtkWidget     *align;
    GdkColor       red;

    gdk_color_parse("red", &red);

    a = (guiAttachment *) malloc(sizeof(guiAttachment));
    a -> parent = s; 
    s -> attachmentList = g_list_append(s -> attachmentList, a);
    a -> nodeList = NULL;

    a -> comboAttachment = gtk_combo_box_new_with_model(GTK_TREE_MODEL(connector_sheet.model));
    gtk_combo_box_set_wrap_width(GTK_COMBO_BOX(a -> comboAttachment), 1);

    renderer = gtk_cell_renderer_text_new();
    gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(a -> comboAttachment), 
                               renderer, FALSE);
    gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(a -> comboAttachment), 
                                   renderer, "text", 0, NULL);
    gtk_cell_layout_set_cell_data_func(GTK_CELL_LAYOUT(a -> comboAttachment), 
                                       GTK_CELL_RENDERER(renderer), 
                                       ComboSensitivity, NULL, NULL);

    delete = ImageButton(GTK_STOCK_DELETE);

    num_rows = s -> tableAttachments -> nrows;
    gtk_table_resize(s -> tableAttachments, num_rows + 1, 3);

    align = gtk_alignment_new(0, 0, 0, 0);
    gtk_container_add(GTK_CONTAINER(align), delete);
    gtk_table_attach(s -> tableAttachments, align,
                     0, 1, num_rows, num_rows + 1, 
                     0, GTK_EXPAND | GTK_FILL, 2, 0);
    gtk_widget_show_all(align);

    align = gtk_alignment_new(0, 0, 0, 0);
    gtk_container_add(GTK_CONTAINER(align), a -> comboAttachment);
    gtk_table_attach(s -> tableAttachments, align, 
                     1, 2, num_rows, num_rows + 1, 
                     GTK_EXPAND | GTK_FILL, GTK_EXPAND | GTK_FILL, 2, 0);
    gtk_widget_show_all(align);

    nodes = gtk_expander_new_with_mnemonic("Nodes");
    align = gtk_alignment_new(0, 0, 0, 0);
    gtk_container_add(GTK_CONTAINER(align), nodes);
    gtk_table_attach(s -> tableAttachments, align,
                     2, 3, num_rows, num_rows + 1, 
                     GTK_EXPAND | GTK_FILL, GTK_EXPAND | GTK_FILL, 2, 0);
    gtk_widget_show_all(align);

    a -> tableNodes = (GtkTable *) gtk_table_new(1, 4, FALSE);
    nodeAdd = gtk_button_new_with_label("Add");
    gtk_table_attach(GTK_TABLE(a -> tableNodes), nodeAdd, 
                     0, 1, 0, 1, 0, 0, 4, 0);
    gtk_table_attach(GTK_TABLE(a -> tableNodes), gtk_label_new("node"),
                     1, 2, 0, 1, 0, 0, 2, 0);
    gtk_table_attach(GTK_TABLE(a -> tableNodes), gtk_label_new("depth"),
                     2, 3, 0, 1, 0, 0, 2, 0);
    gtk_table_attach(GTK_TABLE(a -> tableNodes), gtk_label_new("altit"),
                     3, 4, 0, 1, 0, 0, 2, 0);
    gtk_container_add(GTK_CONTAINER(nodes), GTK_WIDGET(a -> tableNodes));
    gtk_signal_connect(GTK_OBJECT(nodeAdd), "clicked",
                       GTK_SIGNAL_FUNC(AddAttachmentNode), a);


    gtk_signal_connect(GTK_OBJECT(delete), "clicked",
                       GTK_SIGNAL_FUNC(DeleteAttachment), a);
    
    gtk_widget_show_all(nodes);

    gtk_widget_modify_base(s -> expanderAttach, GTK_STATE_NORMAL, &red);
    gtk_widget_modify_text(s -> expanderAttach, GTK_STATE_NORMAL, &red);
    gtk_expander_set_label(GTK_EXPANDER(s -> expanderAttach), "<span foreground=\"red\">_Attach</span>");

    return a;
}

static int
compare (const void *p1, const void *p2)
{
    int n1 = (int) p1;
    int n2 = (int) p2;

    if (n1 < n2)
        return -1;
    else if (n1 == n2)
        return 0;
    else
        return 1;
}

static void
FillSegmentBox(guiSegment *gs, Segment seg, gboolean material_is_complete)
{
    guiAttachmentNode *gn;
    guiAttachment *ga;
    char           buff[20];
    int            j, k;

    sprintf(buff, "%.1f", seg -> length);
    gtk_entry_set_text(GTK_ENTRY(gs -> entryLength), buff);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(gs -> spinNodes), 
                              (double) seg -> dist[1].nodes);
    if (seg -> material) {
        if (material_is_complete)
            ComboBoxSetText(GTK_COMBO_BOX(gs -> comboMaterial), 
                            seg -> material -> name, 1);
        else
            ComboBoxSetText(GTK_COMBO_BOX(gs -> comboMaterial), 
                            (char *) seg -> material, 1);
    }
    if (seg -> connector || seg -> connection == Pinned) {
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(gs -> checkPin), TRUE);
    }
    else {
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(gs -> checkPin), FALSE);
    }

    for (k = 1 ; k <= seg -> num_attach ; k++) {
        ga = AddAttachment(NULL, gs);    
        if (material_is_complete) 
            ComboBoxSetText(GTK_COMBO_BOX(ga -> comboAttachment),
                            seg -> attach[k].object -> name, 1);
        else
            ComboBoxSetText(GTK_COMBO_BOX(ga -> comboAttachment),
                            (char *) seg -> attach[k].object, 1);
        qsort(&(seg -> attach[k].nodes[1]), seg -> attach[k].num_nodes, sizeof(int), compare);
        for (j = seg -> attach[k].num_nodes; j >= 1; j--) {
            gn = AddAttachmentNode(NULL, ga); 
            gtk_spin_button_set_value(GTK_SPIN_BUTTON(gn -> spinNode), 
                                      (double) seg -> attach[k].nodes[j]);
            ChangeAttachmentNode(gn -> spinNode, gn);
        }
    }
    
    return;
}

guiSegment *
BuildSegmentRow(int row)
{
    int              i;
    guiSegment      *s;
    GtkWidget       *align;
    GtkAdjustment   *adj;
    GtkWidget       *attachAdd;
    GtkCellRenderer *renderer;
    GtkRequisition   req;
    GtkWidget       *drags[3];
       

    s = (guiSegment *) malloc(sizeof(guiSegment));

    s -> attachmentList = NULL;

    if (row == 0) 
        segmentList = g_list_prepend(segmentList, s);
    else
        segmentList = g_list_insert(segmentList, s, row - 1);

    adj = (GtkAdjustment *) gtk_adjustment_new(10.0, 3.0, 10000.0, 
                                               1.0, 10.0, 0);
    
    s -> labelNumber = gtk_label_new("-");
    align = gtk_alignment_new(0, 0, 0, 0);
    gtk_container_add(GTK_CONTAINER(align), s -> labelNumber);
    gtk_table_attach(GTK_TABLE(layout_table), align,
                     NUMBER, NUMBER + 1, row, row + 1,
                     0, GTK_EXPAND | GTK_FILL, 4, 0);

    s -> comboMaterial = gtk_combo_box_new_with_model(GTK_TREE_MODEL(material_sheet.model));
    gtk_combo_box_set_wrap_width(GTK_COMBO_BOX(s -> comboMaterial), 1);
    align = gtk_alignment_new(0, 0, 0, 0);
    gtk_container_add(GTK_CONTAINER(align), s -> comboMaterial);
    gtk_table_attach(GTK_TABLE(layout_table), align,
                     MATERIAL, MATERIAL + 1, row, row + 1,
                     GTK_EXPAND | GTK_FILL, GTK_EXPAND | GTK_FILL, 0, 0);

    renderer = gtk_cell_renderer_text_new();
    gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(s -> comboMaterial), 
                               renderer, FALSE);
    gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(s -> comboMaterial), 
                                   renderer, "text", 0, NULL);
    gtk_cell_layout_set_cell_data_func(GTK_CELL_LAYOUT(s -> comboMaterial), 
                                       GTK_CELL_RENDERER(renderer), 
                                       ComboSensitivity, NULL, NULL);

    s -> checkPin = gtk_check_button_new();
    align = gtk_alignment_new(0, 0, 0, 0);
    gtk_container_add(GTK_CONTAINER(align), s -> checkPin);
    gtk_table_attach(GTK_TABLE(layout_table), align,
                     CONNECTOR, CONNECTOR + 1, row, row + 1,
                     0, GTK_EXPAND | GTK_FILL, 0, 0);

    s -> spinNodes = gtk_spin_button_new(adj, 1.0, 0);
    align = gtk_alignment_new(0, 0, 0, 0);
    gtk_container_add(GTK_CONTAINER(align), s -> spinNodes);
    gtk_table_attach(GTK_TABLE(layout_table), align,
                     NODES, NODES + 1, row, row + 1,
                     GTK_EXPAND | GTK_FILL, GTK_EXPAND | GTK_FILL, 0, 0);

    s -> entryLength = gtk_entry_new();
    align = gtk_alignment_new(0, 0, 0, 0);
    gtk_container_add(GTK_CONTAINER(align), s -> entryLength);
    gtk_table_attach(GTK_TABLE(layout_table), align,
                     LENGTH, LENGTH + 1, row, row + 1,
                     GTK_EXPAND | GTK_FILL, GTK_EXPAND | GTK_FILL, 0, 0);

    s -> expanderAttach = gtk_expander_new_with_mnemonic("_Attach");
    gtk_expander_set_use_markup(GTK_EXPANDER(s -> expanderAttach), TRUE);
    align = gtk_alignment_new(0, 0, 0, 0);
    gtk_container_add(GTK_CONTAINER(align), s -> expanderAttach);
    gtk_table_attach(GTK_TABLE(layout_table), align,
                     ATTACH, ATTACH + 1, row, row + 1,
                     GTK_EXPAND | GTK_FILL, GTK_EXPAND | GTK_FILL, 0, 0);

    s -> tableAttachments = GTK_TABLE(gtk_table_new(1, 3, FALSE));
    attachAdd = gtk_button_new_with_label("Add");
    gtk_table_attach(GTK_TABLE(s -> tableAttachments), attachAdd, 
                     0, 1, 0, 1, 0, 0, 4, 0);
    gtk_table_attach(GTK_TABLE(s -> tableAttachments), 
                     gtk_label_new("attachment"),
                     1, 2, 0, 1, 0, 0, 2, 0);
    gtk_container_add(GTK_CONTAINER(s -> expanderAttach), 
                      GTK_WIDGET(s -> tableAttachments));
    gtk_signal_connect(GTK_OBJECT(attachAdd), "clicked",
                       GTK_SIGNAL_FUNC(AddAttachment), s);

    s -> labelDepth = gtk_label_new("-");
    align = gtk_alignment_new(0, 0, 0, 0);
    gtk_container_add(GTK_CONTAINER(align), s -> labelDepth);
    gtk_table_attach(GTK_TABLE(layout_table), align,
                     RUNNING_DEPTH, RUNNING_DEPTH + 1, row, row + 1,
                     GTK_EXPAND | GTK_FILL, GTK_EXPAND | GTK_FILL, 0, 0);

    s -> labelAltit = gtk_label_new("-");
    align = gtk_alignment_new(0, 0, 0, 0);
    gtk_container_add(GTK_CONTAINER(align), s -> labelAltit);
    gtk_table_attach(GTK_TABLE(layout_table), align,
                     TOTAL_LENGTH, TOTAL_LENGTH + 1, row, row + 1,
                     GTK_EXPAND | GTK_FILL, GTK_EXPAND | GTK_FILL, 0, 0);

    s -> labelNodes = gtk_label_new("-");
    align = gtk_alignment_new(0, 0, 0, 0);
    gtk_container_add(GTK_CONTAINER(align), s -> labelNodes);
    gtk_table_attach(GTK_TABLE(layout_table), align,
                     TOTAL_NODES, TOTAL_NODES + 1, row, row + 1,
                     GTK_EXPAND | GTK_FILL, GTK_EXPAND | GTK_FILL, 0, 0);

    gtk_entry_set_width_chars(GTK_ENTRY(s -> entryLength), 8);
    gtk_label_set_max_width_chars(GTK_LABEL(s -> labelAltit), 8);
    gtk_label_set_max_width_chars(GTK_LABEL(s -> labelNodes), 5);
    gtk_label_set_width_chars(GTK_LABEL(s -> labelAltit), 8);
    gtk_label_set_width_chars(GTK_LABEL(s -> labelNodes), 5);
    gtk_label_set_justify(GTK_LABEL(s -> labelAltit), GTK_JUSTIFY_LEFT);

    gtk_widget_size_request(s -> comboMaterial, &req);
    gtk_cell_renderer_set_fixed_size(renderer, req.width, -1);

    gtk_signal_connect(GTK_OBJECT(s -> labelNumber), "focus-in-event",
                       GTK_SIGNAL_FUNC(ActivateRow), s);

    gtk_signal_connect(GTK_OBJECT(s -> spinNodes), "focus-in-event",
                       GTK_SIGNAL_FUNC(ActivateRow), s);
    gtk_signal_connect(GTK_OBJECT(s -> spinNodes), "value-changed",
                       GTK_SIGNAL_FUNC(ComputeNodes), (gpointer) 1);

    gtk_signal_connect(GTK_OBJECT(s -> entryLength), "focus-in-event",
                       GTK_SIGNAL_FUNC(ActivateRow), s);
    gtk_signal_connect(GTK_OBJECT(s -> entryLength), "focus-out-event",
                       GTK_SIGNAL_FUNC(ComputeTotals), (gpointer) 1);

    gtk_combo_box_set_focus_on_click(GTK_COMBO_BOX(s -> comboMaterial), TRUE);
    gtk_signal_connect(GTK_OBJECT(s -> comboMaterial), "focus",
                       GTK_SIGNAL_FUNC(ActivateRow), s);
    gtk_signal_connect(GTK_OBJECT(s -> comboMaterial), "changed",
                       GTK_SIGNAL_FUNC(MaterialChanged), s -> entryLength);


    gtk_widget_show_all(GTK_WIDGET(layout_table));
    // gtk_widget_show(s -> labelNumber);
    // gtk_widget_show(s -> comboMaterial);
    // gtk_widget_show(s -> checkPin);
    // gtk_widget_show(s -> entryLength);
   //  gtk_widget_show(s -> spinNodes);
   //  gtk_widget_show_all(attach);
   //  gtk_widget_show(s -> labelDepth);
   //  gtk_widget_show(s -> labelAltit);
   //  gtk_widget_show(s -> labelNodes);

    gtk_label_set_selectable(GTK_LABEL(s -> labelNumber), TRUE);

    drags[0] = s -> entryLength;
    drags[1] = s -> spinNodes;
    drags[2] = s -> labelNumber;

    for (i = 0 ; i < 3 ; i++) {
        gtk_drag_source_set(drags[i], GDK_BUTTON1_MASK, 
                            targetEntry, 1, GDK_ACTION_MOVE);
        gtk_drag_source_set_icon_stock(drags[i], GTK_STOCK_PASTE);
        gtk_drag_dest_set(drags[i], GTK_DEST_DEFAULT_ALL, 
                          targetEntry, 1, GDK_ACTION_MOVE);

        gtk_signal_connect(GTK_OBJECT(drags[i]), "drag-begin",
                           GTK_SIGNAL_FUNC(SourceBeginDrag), NULL);
        gtk_signal_connect(GTK_OBJECT(drags[i]), "drag-drop", 
                           GTK_SIGNAL_FUNC(DestEndDrag), NULL);

        g_object_set_data(G_OBJECT(drags[i]), "segment", s);
    }

    return s;
}

static void
HideTerminalOptions(GtkToggleButton *w, gpointer table)
{
    if (gtk_toggle_button_get_active(w))
        gtk_widget_show_all(GTK_WIDGET(table));
    else
        gtk_widget_hide_all(GTK_WIDGET(table));
}

static void
AnchorSelector(GtkToggleButton *w, gpointer table)
{
    if (gtk_toggle_button_get_active(w)) {
        ((TerminalTable *) table) -> object = TERMINAL_ANCHOR;
        gtk_combo_box_set_model(GTK_COMBO_BOX(((TerminalTable *) table) -> select), 
                                GTK_TREE_MODEL(anchor_sheet.model));
    }
}

static void
BuoySelector(GtkToggleButton *w, gpointer table)
{
    if (gtk_toggle_button_get_active(w)) {
        ((TerminalTable *) table) -> object = TERMINAL_BUOY;
        gtk_combo_box_set_model(GTK_COMBO_BOX(((TerminalTable *) table) -> select), 
                                GTK_TREE_MODEL(buoy_sheet.model));
    }
}

char *
WriteTerminal(TerminalTable *table, FILE *fp)
{
    char    *type;
    const gchar    *x, *y, *z;
    char    *ret;


    ret = ComboBoxGetText(GTK_COMBO_BOX(table -> select), -1);
    if (!ret) {
        return NULL;
    }

    fprintf(fp, "    terminal = {\n");

    if (table -> object == TERMINAL_ANCHOR) {    
        qprintf(fp, "       anchor = ", 0, ret, 1, "\n", 0, NULL);
    }
    else {
        qprintf(fp, "       buoy = ", 0, ret, 1, "\n", 0, NULL);
    }

    WriteEntry(table -> release, "release-time", fp);
    WriteEntry(table -> safety, "safety", fp);
    WriteEntry(table -> friction, "mu", fp);
    
    // type = ComboBoxGetText(GTK_COMBO_BOX(table -> type), -1);    
    type = gtk_combo_box_get_active_text(GTK_COMBO_BOX(table -> type));
    if (type && type[0] != 0) {
        x = gtk_entry_get_text(GTK_ENTRY(table -> x));
        y = gtk_entry_get_text(GTK_ENTRY(table -> y));
        z = gtk_entry_get_text(GTK_ENTRY(table -> z));
        if (strcmp(type, "speed") == 0) {
            fprintf(fp, "       x-speed = %s\n", x);
            fprintf(fp, "       y-speed = %s\n", y);
            fprintf(fp, "       z-speed = %s\n", z);
        }
        else if (strcmp(type, "force") == 0) {
            fprintf(fp, "       x-force = %g\n", atof(x));
            fprintf(fp, "       y-force = %g\n", atof(y));
            fprintf(fp, "       z-force = %g\n", atof(z));
        }
        else if (strcmp(type, "thrust") == 0) {
            fprintf(fp, "       x-thrust = %s\n", x);
            fprintf(fp, "       y-thrust = %s\n", y);
            fprintf(fp, "       z-thrust = %s\n", z);
        }
        else if (strcmp(type, "position") == 0) {
            fprintf(fp, "       x = %g\n", atof(x));
            fprintf(fp, "       y = %g\n", atof(y));
            fprintf(fp, "       z = %g\n", atof(z));
        }

        g_free(type);
    }
    fprintf(fp, "    }\n");

    return ret;
}

void
BuildTerminalTable(TerminalTable *table, int default_radio, int number)
{
    GtkCellRenderer *renderer;
    GtkWidget   *hbox1;
    GtkWidget   *adv_table;
    GtkWidget   *more;
    char         buff[12];

    table -> vbox = gtk_vbox_new(FALSE, 1);

    table -> buoy = gtk_radio_button_new_with_label(NULL, "buoy");
    table -> anchor = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(table -> buoy), "anchor");

    if (default_radio == TERMINAL_BUOY) {
        table -> select = gtk_combo_box_new_with_model(GTK_TREE_MODEL(buoy_sheet.model));
        table -> object = TERMINAL_BUOY;
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(table -> buoy), TRUE);
    }
    else {
        table -> select = gtk_combo_box_new_with_model(GTK_TREE_MODEL(anchor_sheet.model));
        table  -> object = TERMINAL_ANCHOR;
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(table -> anchor), TRUE);
    }

    gtk_combo_box_set_wrap_width(GTK_COMBO_BOX(table -> select), 1);

    gtk_signal_connect(GTK_OBJECT(table -> buoy), "toggled", GTK_SIGNAL_FUNC(BuoySelector), table);
    gtk_signal_connect(GTK_OBJECT(table -> anchor), "toggled", GTK_SIGNAL_FUNC(AnchorSelector), table);

    hbox1 = gtk_hbox_new(FALSE, 1);

    more = gtk_toggle_button_new_with_label("...");
    sprintf(buff, "terminal %d\n", number);
    gtk_box_pack_start(GTK_BOX(hbox1), gtk_label_new(buff), FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(hbox1), more, FALSE, FALSE, 1);
    gtk_box_pack_start(GTK_BOX(hbox1), table -> buoy, FALSE, FALSE, 1);
    gtk_box_pack_start(GTK_BOX(hbox1), table -> anchor, FALSE, FALSE, 1);
    gtk_box_pack_start(GTK_BOX(hbox1), table -> select, FALSE, FALSE, 1);

    gtk_box_pack_start(GTK_BOX(table -> vbox), hbox1, FALSE, FALSE, 1);

    renderer = gtk_cell_renderer_text_new();
    gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(table -> select), 
                               renderer, FALSE);
    gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(table -> select), 
                                   renderer, "text", 0, NULL);
    gtk_cell_layout_set_cell_data_func(GTK_CELL_LAYOUT(table -> select), 
                                       GTK_CELL_RENDERER(renderer), 
                                       ComboSensitivity, NULL, NULL);

    table -> x = gtk_entry_new();
    table -> y = gtk_entry_new();
    table -> z = gtk_entry_new();

    table -> type = gtk_combo_box_new_text();
    gtk_combo_box_append_text(GTK_COMBO_BOX(table -> type), "speed");
    gtk_combo_box_append_text(GTK_COMBO_BOX(table -> type), "force");
    gtk_combo_box_append_text(GTK_COMBO_BOX(table -> type), "thrust");
    gtk_combo_box_append_text(GTK_COMBO_BOX(table -> type), "position");

    table -> release = gtk_entry_new();
    table -> safety = gtk_entry_new();
    table -> friction = gtk_entry_new();

    adv_table = gtk_table_new(4, 3, FALSE);

    gtk_table_attach(GTK_TABLE(adv_table), table -> type, 0,1,0,1,0,0,0,0);
    gtk_table_attach(GTK_TABLE(adv_table), gtk_label_new("x"), 
                     1,2,0,1,0,0,2,2);
    gtk_table_attach(GTK_TABLE(adv_table), gtk_label_new("y"), 
                     1,2,1,2,0,0,2,2);
    gtk_table_attach(GTK_TABLE(adv_table), gtk_label_new("z"), 
                     1,2,2,3,0,0,2,2);
    gtk_table_attach(GTK_TABLE(adv_table), table -> x, 2,3,0,1,0,0,0,2);
    gtk_table_attach(GTK_TABLE(adv_table), table -> y, 2,3,1,2,0,0,0,2);
    gtk_table_attach(GTK_TABLE(adv_table), table -> z, 2,3,2,3,0,0,0,2);

    gtk_table_attach(GTK_TABLE(adv_table), gtk_label_new("release-time"),
                     0,2,3,4,0,0,2,2);
    gtk_table_attach(GTK_TABLE(adv_table), table -> release, 2,3,3,4,0,0,0,2);
    gtk_table_attach(GTK_TABLE(adv_table), gtk_label_new("friction-coeff"),
                     0,2,4,5,0,0,2,2);
    gtk_table_attach(GTK_TABLE(adv_table), table -> friction, 2,3,4,5,0,0,0,2);
    gtk_table_attach(GTK_TABLE(adv_table), gtk_label_new("safety-factor"),
                     0,2,5,6,0,0,2,2);
    gtk_table_attach(GTK_TABLE(adv_table), table -> safety, 2,3,5,6,0,0,0,2);
    
    

    gtk_box_pack_start(GTK_BOX(table -> vbox), adv_table, FALSE, FALSE, 1);

    gtk_widget_show_all(table -> vbox);
    gtk_widget_hide_all(adv_table);

    gtk_signal_connect(GTK_OBJECT(more), "toggled", 
                       GTK_SIGNAL_FUNC(HideTerminalOptions), adv_table);
    return;
}

void 
ClearLayout(gboolean rebuild_segment_one)
{
    guiSegment *s;


    active_row = -1; // so we don't try to replace it
    while (layout_table -> nrows > 1) {
        DeleteSegment(layout_table -> nrows - 1);
    }

    g_list_free(segmentList);
    segmentList = NULL;

    active_row = -1;
    if (rebuild_segment_one) {
        s = BuildSegmentRow(1);
        ActivateRow(NULL, NULL, s);
    }
}

static void
ClearLayoutUser(GtkWidget *w, gpointer data)
{
    GtkWidget   *dialog;

    dialog = gtk_message_dialog_new(GTK_WINDOW(toplevel), GTK_DIALOG_MODAL, 
                                    GTK_MESSAGE_QUESTION, GTK_BUTTONS_YES_NO, 
                                    "clear entire layout?");

    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_YES) {
        ClearLayout(TRUE);
        ComputeTotals(NULL, NULL, NULL);
    }

    gtk_widget_destroy(dialog);
}

static struct segment segCopy;
static int clipboardFull = 0;

static void
CopyActiveSegment(GtkWidget *w, gpointer data)
{
    if (active_row > -1) {
        DecodeSegmentBox(active_row, &segCopy, TRUE);
        clipboardFull = 1;
    }
}

static void
PasteBelowActiveSegment(GtkWidget *w, gpointer data)
{
    guiSegment *gs;

    if (clipboardFull) {
        gs = AddNewSegment(active_row + 1, TRUE);
        FillSegmentBox(gs, &segCopy, FALSE);
        FreeSegmentMemory(&segCopy);
    }
}

static void
DeleteActiveSegment(GtkWidget *widget, gpointer data)
{
    if (active_row == -1) {
        return;
    }

    clipboardFull = 1;
    DecodeSegmentBox(active_row, &segCopy, TRUE);
    DeleteSegment(active_row);
}

GtkWidget *
BuildLayoutSheet(void)
{
    int          i;
    GtkWidget   *hpaned;
    GtkWidget   *hbox_tools, *scroller, *vbox1, *vbox2;
    GtkWidget   *new_above, *new_below, *delete, *new, *paste, *copy;
    GtkWidget   *vbox_ctrls, *frame_analysis, *frame_environment;
    GtkWidget   *sep;    
    GtkTooltips *tips;
    GtkWidget   *lengthLabel, *nodesLabel, *numberLabel;
    guiSegment  *gs;

    hpaned = gtk_hpaned_new();
    gtk_widget_show(hpaned);

    tips = gtk_tooltips_new();
    gtk_tooltips_enable(tips);

    // right-hand pane - analysis and environment

    vbox_ctrls = gtk_vbox_new(FALSE, 2);
    frame_analysis = BuildAnalysisParameters();
    gtk_box_pack_start(GTK_BOX(vbox_ctrls), frame_analysis, FALSE, FALSE, 0);
    frame_environment = BuildEnvironment();
    gtk_box_pack_start(GTK_BOX(vbox_ctrls), frame_environment, FALSE, FALSE, 0);
    gtk_widget_show(vbox_ctrls);
    gtk_paned_pack2(GTK_PANED(hpaned), vbox_ctrls, TRUE, TRUE);

    // left-hand pane - layout controls

    vbox1 = gtk_vbox_new(FALSE, 1);
    gtk_widget_show(vbox1);

    gtk_paned_pack1(GTK_PANED(hpaned), vbox1, FALSE, FALSE);

    hbox_tools = gtk_hbox_new(FALSE, 1);
    gtk_box_pack_start(GTK_BOX(vbox1), hbox_tools, FALSE, FALSE, 0);

    new_above = ImageButton(GTK_STOCK_GO_UP);
    new_below = ImageButton(GTK_STOCK_GO_DOWN);
    delete = ImageButton(GTK_STOCK_CUT);
    copy = ImageButton(GTK_STOCK_COPY);
    paste = ImageButton(GTK_STOCK_PASTE);
    new = ImageButton(GTK_STOCK_DELETE);
    gtk_box_pack_start(GTK_BOX(hbox_tools), new_above, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox_tools), new_below, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox_tools), delete, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox_tools), copy, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox_tools), paste, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox_tools), new, FALSE, FALSE, 0);
    gtk_widget_show(hbox_tools);

    gtk_tooltips_set_tip(tips, new_above, "add segment above current segment", NULL);
    gtk_tooltips_set_tip(tips, new_below, "add segment below current segment", NULL);
    gtk_tooltips_set_tip(tips, delete, "cut current segment", NULL);
    gtk_tooltips_set_tip(tips, new, "clear layout (delete all segments)", NULL);
    gtk_tooltips_set_tip(tips, copy, "copy segment to clipboard", NULL);
    gtk_tooltips_set_tip(tips, paste, "paste segment from clipboard below active segment", NULL);

    vbox2 = gtk_vbox_new(FALSE, 1);

    segmentList = NULL;

    layout_table = GTK_TABLE(gtk_table_new(2, 7, FALSE));
    for (i = 0 ; i < NUM_LAYOUT_COLUMNS ; i++) {
        gtk_table_attach(layout_table,
                         gtk_label_new(layout_headers[i]),
                         i, i+1, 0, 1,
                         0, 0, 0, 0);
    }

    lengthLabel = gtk_table_get_child_at(layout_table, 0, LENGTH);
    nodesLabel = gtk_table_get_child_at(layout_table, 0, NODES);
    numberLabel = gtk_table_get_child_at(layout_table, 0, NUMBER);

    gtk_drag_dest_set(numberLabel, GTK_DEST_DEFAULT_ALL, 
                      targetEntry, 1, GDK_ACTION_MOVE);
    gtk_signal_connect(GTK_OBJECT(numberLabel), "drag-drop", 
                       GTK_SIGNAL_FUNC(DestEndDrag), NULL);
    gtk_drag_dest_set(lengthLabel, GTK_DEST_DEFAULT_ALL, 
                      targetEntry, 1, GDK_ACTION_MOVE);
    gtk_signal_connect(GTK_OBJECT(lengthLabel), "drag-drop", 
                       GTK_SIGNAL_FUNC(DestEndDrag), NULL);
    gtk_drag_dest_set(nodesLabel, GTK_DEST_DEFAULT_ALL, 
                      targetEntry, 1, GDK_ACTION_MOVE);
    gtk_signal_connect(GTK_OBJECT(nodesLabel), "drag-drop", 
                       GTK_SIGNAL_FUNC(DestEndDrag), NULL);

    gs = BuildSegmentRow(1);

    BuildTerminalTable(&terminal_1, TERMINAL_ANCHOR, 1);
    BuildTerminalTable(&terminal_2, TERMINAL_BUOY, 2);
    gtk_box_pack_start(GTK_BOX(vbox2), terminal_2.vbox, FALSE, FALSE, 0);
    sep = gtk_hseparator_new(); gtk_widget_show(sep);
    gtk_box_pack_start(GTK_BOX(vbox2), sep, FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(vbox2), GTK_WIDGET(layout_table), FALSE, FALSE, 0);
    sep = gtk_hseparator_new(); gtk_widget_show(sep);
    gtk_box_pack_start(GTK_BOX(vbox2), sep, FALSE, FALSE, 3);
    gtk_box_pack_start(GTK_BOX(vbox2), terminal_1.vbox, FALSE, FALSE, 0);
    gtk_widget_show_all(GTK_WIDGET(layout_table));
    gtk_widget_show(vbox2);

    scroller = gtk_scrolled_window_new(NULL, NULL);
    gtk_box_pack_start(GTK_BOX(vbox1), scroller, TRUE, TRUE, 0);

    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller), 
                                   GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scroller), vbox2);
    gtk_widget_show(scroller);

    gtk_signal_connect(GTK_OBJECT(new_above), "clicked",
                       GTK_SIGNAL_FUNC(AddNewSegmentActive), (gpointer) 0);
    gtk_signal_connect(GTK_OBJECT(new_below), "clicked",
                       GTK_SIGNAL_FUNC(AddNewSegmentActive), (gpointer) 1);
    gtk_signal_connect(GTK_OBJECT(delete), "clicked",
                       GTK_SIGNAL_FUNC(DeleteActiveSegment), NULL);
    gtk_signal_connect(GTK_OBJECT(new), "clicked",
                       GTK_SIGNAL_FUNC(ClearLayoutUser), NULL);
    gtk_signal_connect(GTK_OBJECT(copy), "clicked",
                       GTK_SIGNAL_FUNC(CopyActiveSegment), NULL);
    gtk_signal_connect(GTK_OBJECT(paste), "clicked",
                       GTK_SIGNAL_FUNC(PasteBelowActiveSegment), NULL);

    ActivateRow(NULL, NULL, gs);

    ComputeTotals(NULL, NULL, NULL);


    return hpaned;
}

static void
NewModel(gpointer data, guint action, GtkWidget *w)
{
    char *path;
    Problem p;
    Analysis a;
    Environment e;
    struct stat buf;

    ClearEnvironment();
    ClearAnalysis();

    // clear the old model ...
    ClearLayout(TRUE);
    ComputeTotals(NULL, NULL, NULL);

    SetCurrentFile("new problem");

    path = DatabasePath("defaults.cbl");
    if (stat(path, &buf)) {
        return;
    }

    fprintf(stderr,"reading defaults file %s\n", path);
    p.buoy_tree = p.material_tree = p.anchor_tree = p.connector_tree = NULL;
    if (ReadModelFile(path, &p, &a, &e, 1)) {
        fprintf(stderr,"did not get defaults\n");
        return;
    }

    // FillIn ??
    FillEnvironment(&e);
    FillAnalysis(&a, &p);
}

static int
ShowOutputOptions(gpointer data, guint action, GtkWidget *w)
{
    GtkWidget   *win, *vbox, *check[6];
    GError      *err = NULL;
    gint         result;
    gboolean     state;
    int          i;
    char        *labels[] = {
        "always show this dialog for new or loaded solutions",
        "tabulated results",
        "position time series",
        "position snapshots",
        "tension time series",
        "tension snapshots",
    };
    char        *keys[] = {
        "alwaysShowOutputOptions",
        "tabulate",
        "positionTimeSeries",
        "positionSnapshots",
        "tensionTimeSeries",
        "tensionSnapshots",
    };
    gboolean     defaults[] = {
        1, 1, 0, 1, 1, 0,
    };

    win = gtk_dialog_new_with_buttons("output options", 
                                      GTK_WINDOW(toplevel),
                                      GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
                                      GTK_STOCK_OK,
                                      GTK_RESPONSE_ACCEPT,
                                      GTK_STOCK_CANCEL,
                                      GTK_RESPONSE_REJECT,
                                      NULL);

    vbox = GTK_DIALOG(win) -> vbox; 

    for (i = 0 ; i < 6 ; i++) {
        check[i] = gtk_check_button_new_with_label(labels[i]);
        gtk_box_pack_start(GTK_BOX(vbox), check[i], FALSE, TRUE, 0);

        state = g_key_file_get_boolean(prefsFile, "output", keys[i], &err);
        if (err) {
            state = defaults[i];
            err = NULL;
        }
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(check[i]), state);
    }
    gtk_widget_show_all(vbox);
    result = gtk_dialog_run(GTK_DIALOG(win));

    if (result == GTK_RESPONSE_ACCEPT) {
        for (i = 0 ; i < 6 ; i++) {
            g_key_file_set_boolean(prefsFile, "output", keys[i], 
                     gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(check[i])));
        }
    }

    gtk_widget_destroy(win);
    
    return (result == GTK_RESPONSE_ACCEPT ? 0 : 1);
}

static void
FillTerminalTable(Terminal terminal, TerminalTable *table)
{

    if (terminal -> anchor) {
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(table -> anchor), TRUE);
        ComboBoxSetText(GTK_COMBO_BOX(table -> select), 
                        terminal -> anchor -> name, 1);
    }
    else if (terminal -> buoy) {
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(table -> buoy), TRUE);
        ComboBoxSetText(GTK_COMBO_BOX(table -> select), 
                        terminal -> buoy -> name, 1);
    }

    if (terminal -> xspeed.expr || terminal -> xspeed.value
        || terminal -> yspeed.expr || terminal -> yspeed.value 
        || terminal -> zspeed.expr || terminal -> zspeed.value) {

        printf("setting speed\n");
        ComboBoxSetText(GTK_COMBO_BOX(table -> type), "speed", 1);
        gtk_combo_box_set_active(GTK_COMBO_BOX(table -> type), 0);
        if (terminal -> yspeed.expr)
            gtk_entry_set_text(GTK_ENTRY(table -> x), terminal -> yspeed.text);
        else  
            SetNumericText(table -> x, "%g", terminal -> yspeed.value, FALSE);

        if (terminal -> zspeed.expr)
            gtk_entry_set_text(GTK_ENTRY(table -> y), terminal -> zspeed.text);
        else  
            SetNumericText(table -> y, "%g", terminal -> zspeed.value, FALSE);

        if (terminal -> xspeed.expr)
            gtk_entry_set_text(GTK_ENTRY(table -> z), terminal -> xspeed.text);
        else  
            SetNumericText(table -> z, "%g", terminal -> xspeed.value, FALSE);
    }
    else if (terminal -> xforce || terminal -> yforce || terminal -> zforce) {
        ComboBoxSetText(GTK_COMBO_BOX(table -> type), "force", 1);
        gtk_combo_box_set_active(GTK_COMBO_BOX(table -> type), 1);
        SetNumericText(table -> x, "%g", terminal -> yforce, FALSE);
        SetNumericText(table -> y, "%g", terminal -> zforce, FALSE);
        SetNumericText(table -> z, "%g", terminal -> xforce, FALSE);
    }
    else if (terminal -> xthrust.expr || terminal -> xthrust.value
            || terminal -> ythrust.expr || terminal -> ythrust.value 
            || terminal -> zthrust.expr || terminal -> zthrust.value) {
        ComboBoxSetText(GTK_COMBO_BOX(table -> type), "thrust", 1);
        gtk_combo_box_set_active(GTK_COMBO_BOX(table -> type), 2);
        if (terminal -> ythrust.expr)
            gtk_entry_set_text(GTK_ENTRY(table -> x), terminal -> ythrust.text);
        else  
            SetNumericText(table -> x, "%g", terminal -> ythrust.value, FALSE);

        if (terminal -> zthrust.expr)
            gtk_entry_set_text(GTK_ENTRY(table -> y), terminal -> zthrust.text);
        else  
            SetNumericText(table -> y, "%g", terminal -> zthrust.value, FALSE);

        if (terminal -> xthrust.expr)
            gtk_entry_set_text(GTK_ENTRY(table -> z), terminal -> xthrust.text);
        else  
            SetNumericText(table -> z, "%g", terminal -> xthrust.value, FALSE);
    }
    else {
        ComboBoxSetText(GTK_COMBO_BOX(table -> type), "position", 1);
        gtk_combo_box_set_active(GTK_COMBO_BOX(table -> type), 3);
        SetNumericText(table -> x, "%g", terminal -> y, FALSE);
        SetNumericText(table -> y, "%g", terminal -> z, FALSE);
        SetNumericText(table -> z, "%g", terminal -> x, FALSE);
    }
}

static void
PushRecent(char *name)
{
    int i, n;
    GtkWidget   *label;

    if (name == NULL) {
        return;
    }

    if (g_queue_find_custom(recentFiles, name, strcmp))
        return;

    g_queue_push_head(recentFiles, strdup(name));

    n = g_queue_get_length(recentFiles);
    for (i = 0 ; i < n && i < 4 ; i++) {
        name = g_queue_peek_nth(recentFiles, i);
        label = gtk_bin_get_child(GTK_BIN(recentWidgets[i]));
        gtk_label_set_text(GTK_LABEL(label), basename(name));
        gtk_widget_show(recentWidgets[i]);
    }
}

void
LatchRecent()
{
    int      i,n;
    char     key[12];
    char    *name;

    n = g_queue_get_length(recentFiles);
    for (i = 0 ; i < n && i < 4 ; i++) {
        name = g_queue_peek_nth(recentFiles, i);
        sprintf(key, "recent%d", i + 1);
        g_key_file_set_string(prefsFile, "files", key, name);
    }
}

extern void ControlPlotTimeSeries(Result *, Problem *, Environment *);
extern void BuildPlayback(Solution *);

char *prev_file = NULL;

static void
OpenCableFile(char *filename, int results_mode)
{
    Segment   *seg;
    int        ns, i, j;
    Problem   *gui_problem;
    Analysis  *gui_analysis;
    Environment *gui_environment;
    Solution    *solution;
    char      *cabname;
    Result    *res = NULL;
    ResFile    in;
    GError    *err = NULL;
    int        fid;
    gboolean   show_output_options;
    guiSegment *gs = NULL;
    char       path[PATH_MAX];

    gui_problem = (Problem *) calloc(1, sizeof(Problem));
    gui_analysis = (Analysis *) calloc(1, sizeof(Analysis));
    gui_environment = (Environment *) calloc(1, sizeof(Environment));
    gui_problem -> material_tree = NULL;
    gui_problem -> anchor_tree = NULL;
    gui_problem -> connector_tree = NULL;
    gui_problem -> buoy_tree = NULL;

    // gui_problem.material_tree = material_sheet.db_tree;
    // gui_problem.anchor_tree = anchor_sheet.db_tree;
    // gui_problem.connector_tree = connector_sheet.db_tree;
    // gui_problem.buoy_tree = buoy_sheet.db_tree;
    if (results_mode) {
        in = res_open(filename, "rb");
        if (in == NULL) {
            error("OpenCableFile(1): Unable to open %s", filename);
            return;
        }
        res = ReadResultsFile(in, 0, 1, 0, 0);

        if (!res) {
            res_close(in);
            return;
        }

        fid = g_file_open_tmp("wcabXXXXXX", &cabname, &err); 
        if (fid < 0) {
            error("error creating temp file (OpenCableFile cab)");
            res_close(in);
            return;
        }
        else {
            write(fid, res -> input, strlen(res -> input));
            close(fid);
        }

        if (ReadModelFile(cabname, gui_problem, gui_analysis, gui_environment, 0)) {
            free(cabname);
            res_close(in);
            return;
        }
        free(cabname);
    }
    else {
        if (ReadModelFile(filename, gui_problem, gui_analysis, gui_environment, 0)) {
            free(gui_problem);
            free(gui_analysis);
            free(gui_environment);
            return;
        }
    }
    
    FillInEnvironment(gui_environment);
    FillInAnalysis(gui_problem, gui_analysis);
    // FillInObjectProperties(gui_problem, gui_analysis, gui_environment);

    bufferTotalChanges = TRUE;

    FillEnvironment(gui_environment);
    FillAnalysis(gui_analysis, gui_problem);

    // clear the old model ...
    ClearLayout(FALSE);

    // depending on the rigor with which we fill environment and
    // analysis we might need to clear those as well ...

    gfDuplicatePref = 0; // reset 

    TreeSetAndIterate(gui_problem -> material_tree, 
                      ObjectIterator, (void *) &material_sheet);
    TreeSetAndIterate(gui_problem -> buoy_tree, 
                      ObjectIterator, (void *) &buoy_sheet);
    TreeSetAndIterate(gui_problem -> connector_tree, 
                      ObjectIterator, (void *) &connector_sheet);
    TreeSetAndIterate(gui_problem -> anchor_tree, 
                      ObjectIterator, (void *) &anchor_sheet);

    seg = BuildSegmentArray(gui_problem, &ns, 0); // seg is unit offset array

    // FillSegmentBox(active_row, seg[1]);
    for (i = 1, j = ns ; i <= ns ; i++, j--) {
        gs = AddNewSegment(i, FALSE);
        FillSegmentBox(gs, seg[j], TRUE);
    }
    
    ActivateRow(NULL, NULL, gs);

    FillTerminalTable(gui_problem -> terminal[1], &terminal_1);
    FillTerminalTable(gui_problem -> terminal[2], &terminal_2);


    SetNumericText(terminal_1.friction, "%.2f", gui_problem -> terminal[1] -> friction, TRUE);
    SetNumericText(terminal_1.safety, "%.2f", gui_problem -> terminal[1] -> safety, TRUE);

    seg ++; free(seg);

    bufferTotalChanges = FALSE;
    ComputeTotals(NULL, NULL, NULL);


    SetCurrentFile(filename);

#ifdef WINDOWS
    GetFullPathNameA(filename, PATH_MAX, path, &cabname);
#else
    realpath(dirname(filename), path);
    strncat(path, "/", PATH_MAX);
    strncat(path, basename(filename), PATH_MAX);
#endif
    PushRecent(path);
    LatchRecent();

    if (prev_file) {
        free(prev_file);
    }
    prev_file = strdup(dirname(path));
 
    if (results_mode) {

        err = NULL;
        show_output_options = g_key_file_get_boolean(prefsFile, "output", "alwaysShowOutputOptions", &err);

        if (show_output_options || err) {
            if (ShowOutputOptions(NULL, 0, NULL)) {
                FreeProblem(gui_problem);
                FreeEnvironment(gui_environment);
                FreeAnalysis(gui_analysis);
                return;
            }
            
            err = NULL;
        }

        solution = (Solution *) malloc(sizeof(Solution));

        solution -> in_name = NULL;
        solution -> out_name = NULL;
        solution -> static_file = NULL;
        solution -> initial_file = NULL;
        solution -> dynstat_name = NULL;
        solution -> tmpdir = NULL;

        solution -> rotation = 0;

        solution -> progress_dt = 0;
        solution -> restart_t = 0;
        solution -> restart_file = NULL;
        solution -> progress_file = NULL;

        solution -> segment_dt = res -> seg_dt;
        solution -> buoy_dt =  res -> buoy_dt;
        solution -> ext_dt = res -> ext_dt;
        solution -> analysis = gui_analysis;
        solution -> problem = gui_problem;
        solution -> environment = gui_environment;
        solution -> output_dt = res -> sample_dt;
        solution -> snapshot_dt = res -> snap_dt;
        solution -> num_output_nodes = res -> num_output_nodes;
        if (res -> num_output_nodes) {
            solution -> output_nodes = (int *) malloc(sizeof(int) * res -> num_output_nodes);
            memcpy(solution -> output_nodes, 
                   res -> output_nodes, res -> num_output_nodes * sizeof(int));
        }
        else
            solution -> output_nodes = NULL;

        for (i = 1 ; i < MAXOUTPUT ; i++) {
            solution -> snapPlot[i] = NULL;
            solution -> tsPlot[i] = NULL;
        }

        solution -> table = NULL;
        solution -> controls = NULL;
        solution -> results = res;
        gui_problem -> solution = solution;
        gui_analysis -> solution = solution;
        gui_environment -> solution = solution;

        if (res -> dynamic) {
            solution -> static_only = 0;
            BuildPlayback(solution);
        }
        else {
            solution -> static_only = 1;
            solution -> playback = NULL;
        }

        ResultsToProblem(res, gui_problem, gui_environment, TRUE);
        ControlTabulateResults(gui_problem, gui_environment);

        ControlMakePlots(solution, 
                         res -> dynamic ? (PLOT_SPACE | PLOT_TIME) : PLOT_SPACE);
        ControlPlotSnaps(gui_problem, gui_environment);
        if (res -> num_output_nodes && res -> sample_dt)
            ControlPlotTimeSeries(res, gui_problem, gui_environment);
    }
    else {

        FreeProblem(gui_problem);
        FreeEnvironment(gui_environment);
        FreeAnalysis(gui_analysis);
    }

    return;
}


static void
OpenModel(gpointer data, guint action, GtkWidget *w)
{
    char      *filename;
    
    if (w) {
        filename = FilenameFromDialog("open model", toplevel,    
                                      GTK_FILE_CHOOSER_ACTION_OPEN,
                                      &prev_file, 
                                      "*.cbl", "cable input files",
                                      NULL, NULL);

    }
    else {
        filename = (char *) data;
    }

    if (!filename) {
        return;
    }

    OpenCableFile(filename, FALSE);
}

static void
OpenResults(gpointer data, guint action, GtkWidget *w)
{
    char      *resname;

    resname = FilenameFromDialog("open results", toplevel,
                                GTK_FILE_CHOOSER_ACTION_OPEN, &prev_file,
                                "*.crs", "cable results files", 
                                NULL, NULL);
    if (!resname) {
        return;
    }

    OpenCableFile(resname, TRUE);
}


static GSList *
insert_unique(GSList *list, char *elt)
{
    if (elt && !g_slist_find_custom(list, elt, strcmp)) {
        return g_slist_append(list, elt);
    }

    return list;
}


char *
GUItoModel(char *given_name, int *numnodes, int save_all_objects)
{
    char            basis[32] = "/tmp/cabXXXXXX";
    int             i, j, k;
    char           *file;
    char           *name;
    struct segment  seg;
    FILE           *fp;
    GSList         *materials = NULL;
    GSList         *connectors = NULL;
    GSList         *anchors = NULL;
    GSList         *buoys = NULL;
    int             nn;
    guiSegment     *gs;

    if (given_name) {
        file = given_name;
    }   
    else {
        file = strdup(basis);
        mktemp(file);
    }

    fp = fopen(file, "w");
    if (!fp) {
	    error("could not open file %s for writing", file);
        return NULL;
    }

    WriteAnalysis(fp);
    fprintf(fp, "\n");
    WriteEnvironment(fp);


    fprintf(fp, "\nLayout\n");
    
    name = WriteTerminal(&terminal_1, fp);
    if (terminal_1.object == TERMINAL_ANCHOR) {
        anchors = insert_unique(anchors, name);
    }
    else {
        buoys = insert_unique(buoys, name);
    }

    nn = 0;
    for (i = layout_table -> nrows - 1 ; i >= 1 ; i--) {
        gs = DecodeSegmentBox(i, &seg, TRUE);
        if (seg.material == NULL || seg.length == 0 || seg.dist[1].nodes == 0) {
            continue;
            // skip silently?
        }
        fprintf(fp, "    segment = {\n");
        fprintf(fp, "        length = %g\n", seg.length);
        qprintf(fp, "        material = ", 0, (char *) seg.material, 1, "\n", 0, NULL);
        fprintf(fp, "        nodes = (%d, %g)\n", seg.dist[1].nodes, seg.dist[1].percent);

        for (j = 1 ; j <= seg.num_attach ; j++) {
            connectors = insert_unique(connectors, (char *) seg.attach[j].object);
            if (j == 1)  {
                qprintf(fp, "        attachments = ", 0, 
                            (char *) seg.attach[j].object, 1, " : (", 0, NULL);
            }
            else {
                qprintf(fp, "                      ", 0, 
                            (char *) seg.attach[j].object, 1, " : ( ", 0, NULL);
            }
        
            
            for (k = 1 ; k <= seg.attach[j].num_nodes ; k++) {
                fprintf(fp, "%d ", seg.attach[j].nodes[k]); 
            }
            fprintf(fp, ")\n");
        }
      
        fprintf(fp, "    }\n"); 

        nn += seg.dist[1].nodes;

        if (seg.connection == Pinned) { //  == (void *) 1) {
            fprintf(fp, "    connector = pinned\n");
        }
        else if (seg.connector) { // can't happen in current GUI
            qprintf(fp, "    connector = ", 0, (char *) seg.connector, 1, "\n", 0, NULL);
        }

        materials = insert_unique(materials, (char *) seg.material);
        if (seg.connector && seg.connector != (void *) 1)
            connectors = insert_unique(connectors, (char *) seg.connector);

        FreeSegmentMemory(&seg);
    }

    name = WriteTerminal(&terminal_2, fp);
    if (terminal_2.object == TERMINAL_ANCHOR) {
        anchors = insert_unique(anchors, name);
    }
    else {
        buoys = insert_unique(buoys, name);
    }

    fprintf(fp, "\n");
    WriteBuoys(buoys, fp, save_all_objects);           fprintf(fp, "\n");
    WriteAnchors(anchors, fp, save_all_objects);       fprintf(fp, "\n");
    WriteConnectors(connectors, fp, save_all_objects); fprintf(fp, "\n");
    WriteMaterials(materials, fp, save_all_objects);   fprintf(fp, "\n");
 

    g_slist_free(buoys);
    g_slist_free(materials);
    g_slist_free(connectors);
    g_slist_free(anchors);

    fprintf(fp, "\nEnd\n");

    fclose(fp);

    if (numnodes) 
        *numnodes = nn;
    
    return file;
}

static void
OpenDatabaseFile(gpointer data, guint action, GtkWidget *w)
{
    char      *dbname;

    dbname = FilenameFromDialog("open database", toplevel,
                                GTK_FILE_CHOOSER_ACTION_OPEN, &prev_file,
                                "*.db", "cable db files",
                                NULL, NULL);
    if (!dbname) {
        return;
    }

    ReadDatabaseFile(dbname, NULL, NULL, NULL, NULL);
}

static void
SaveDatabases(gpointer data, guint action, GtkWidget *w)
{
    WriteDatabases();
}


static void
SaveModel(gpointer data, guint action, GtkWidget *w)
{
    char      *filename; 
    
    if (w) {
        if (action || !current_filename) {
            filename = FilenameFromDialog("save model", toplevel, 
                                          GTK_FILE_CHOOSER_ACTION_SAVE, 
                                          &prev_file, 
                                          "*.cbl", "cable input files",
                                          NULL, NULL);
        }
        else {
            filename = strdup(current_filename);
        }
    }
    else {
        filename = (char *) data;
    }

    if (!filename) {
        return;
    }

    filename = AddFilenameExtension(filename, ".cbl");
 
    GUItoModel(filename, NULL, 1);
    SetCurrentFile(filename);
}

extern GtkWidget **BuildControl(Solution *);
extern int SolutionDriver(Solution *, int);

static void
SolveModel(gpointer data, guint action, GtkWidget *w)
{
   ResFile    in;
   int        i, fid, nn;
   Solution  *solution;
   GtkWidget **ctls;
   GError    *err = NULL;
   gboolean   show_output_options;

   solution = (Solution *) malloc(sizeof(Solution));
   solution -> out_name = NULL;

   if (action) { // model summary - we specify text file to export to
        solution -> out_name = FilenameFromDialog("export model summary", toplevel,
                                                  GTK_FILE_CHOOSER_ACTION_SAVE,
                                                  &prev_file, 
                                                  NULL, NULL, NULL, NULL);
        if (!solution -> out_name) {
            free(solution);
            return;
        }
   }
   else {
       fid = g_file_open_tmp("wcrsXXXXXX", &(solution -> out_name), &err); 
       if (fid < 0) {
           free(solution);
           error("error creating temp file %s", solution -> out_name);       
           return;
       }
       else 
           close(fid);
   }

   fid = g_file_open_tmp("wcabXXXXXX", &(solution -> in_name), &err); 
   if (fid < 0) {
       error("error creating temp file %s", solution -> in_name);
       free(solution);
       return;
   }
   close(fid);

   err = NULL;

   if (!action) {
       show_output_options = g_key_file_get_boolean(prefsFile, "output", 
                                            "alwaysShowOutputOptions", &err);
       if (show_output_options || err) {
           if (ShowOutputOptions(NULL, 0, NULL)) {
               if (solution -> out_name) free(solution -> out_name);
               free(solution);
               return;
           }
           err = NULL;
       }
   }

   GUItoModel(solution -> in_name, &nn, 1);
 
   for (i = 1 ; i < MAXOUTPUT ; i++) {
       solution -> output_map[i] = 0;
   }
   solution -> output_map[DISPLACEMENT] = 1;
   solution -> output_map[FORCE] = 1;

   FillSolutionControl(solution, nn);

   solution -> dynstat_name = NULL;
   solution -> table_name = NULL;
   solution -> buoy_dt = 0.0;
   solution -> segment_dt = 0.0;
   solution -> ext_dt = 0.0;
   solution -> tmpdir = NULL;
   solution -> debug_input = 0;
   solution -> bill_of_materials = action;
   solution -> refine = 0;
   solution -> simple = 0;
   solution -> quit = 1;
   solution -> quiet = 0;
   solution -> unlink_input = 0; // take care of this here
   solution -> X = 1;
   solution -> decimate = 1;
   solution -> saved = 0;
   solution -> userQuit = 0;
   solution -> solutionComplete = 0;
   solution -> rotation = 0.0;

   solution -> progress_dt = 0;
   solution -> restart_t = 0;
   solution -> progress_file = NULL;
   solution -> restart_file = NULL;

   solution -> table = NULL;
   solution -> results = NULL;
   for (i = 1 ; i < MAXOUTPUT ; i++) {
       solution -> snapPlot[i] = NULL;
       solution -> tsPlot[i] = NULL;
   }

   solution -> problem = (Problem *) malloc(sizeof(Problem));
   solution -> analysis = (Analysis *) malloc(sizeof(Analysis));
   solution -> environment = (Environment *) malloc(sizeof(Environment));

   solution -> problem -> material_tree = NULL;
   solution -> problem -> anchor_tree = NULL;
   solution -> problem -> buoy_tree = NULL;
   solution -> problem -> connector_tree = NULL;
/*
   solution -> problem -> solution = solution;
   solution -> analysis -> solution = solution;
   solution -> environment -> solution = solution;
*/
   if (solution -> static_only) {
       solution -> plotProgress = 1;
   }

   solution -> plotProgress = 1;

   solution -> playback = NULL; 

   if (!action) {
       ctls = BuildControl(solution);
       gtk_widget_show(ctls[0]);
       solution -> controls = ctls;
   }
   else {
       solution -> controls = NULL;
   }

   if (!action) {
       gtk_widget_set (solve_item, "sensitive", FALSE, NULL);
   }

   SolutionDriver(solution, 0);
   unlink(solution -> in_name);
   solution -> solutionComplete = 1;
   if (!solution -> static_only) {
       for (i = 1 ; i < MAXOUTPUT ; i++) {
           if (solution -> snapPlot[i]) {
               in = res_open(solution -> out_name, "rb");
               if (in) {
                   solution -> results = ReadResultsFile(in, 0, 1, 0, 0);
                   if (solution -> results) {
                       BuildPlayback(solution);
                   }
                   else {
                       res_close(in);
                   }
               }
               // only need to do the above if we found at least one snap plot
               break;
           }
       }
   }
 
   if (!action) {
       gtk_widget_set (solve_item, "sensitive", TRUE, NULL);
   }

   if (action) {
        FreeSolution(solution);
   } 
   
   return;
}


static void
About(gpointer data, guint action, GtkWidget *w)
{
    GtkWidget   *dialog;

    dialog = gtk_message_dialog_new(GTK_WINDOW(toplevel), GTK_DIALOG_MODAL, 
                                    GTK_MESSAGE_INFO, GTK_BUTTONS_OK, 
                                    "WHOI Cable v3.08\ncompiled %s\nbuild %s", 
                                    TIMESTAMP, SVN_Version);

    gtk_dialog_run(GTK_DIALOG(dialog));
    gtk_widget_destroy(dialog);
}


static void
OpenRecent(gpointer data, guint action, GtkWidget *w)
{
    char *name;

    name = g_queue_peek_nth(recentFiles, action - 1);
    if (name) {
        OpenCableFile(name, strstr(name, "crs") ? 1 : 0);
    }
}

static GtkWidget *
makemenu(void)
{
    GtkItemFactory *item_factory;
    GtkAccelGroup *accel_group;

    static GtkItemFactoryEntry menu_items[] = {
        { "/_File",         NULL,         NULL,    0, "<Branch>",  },
        { "/File/_New",    "<control>N", NewModel,    0, "<StockItem>", GTK_STOCK_NEW },
        { "/File/_Open",    "<control>O", OpenModel,    1, "<StockItem>", GTK_STOCK_OPEN },
        { "/File/_Save",    "<control>S", SaveModel,    0, "<StockItem>", GTK_STOCK_SAVE },
        { "/File/Save _As", NULL,         SaveModel,    1, "<StockItem>", GTK_STOCK_SAVE_AS },
        { "/File/Open _Results",    "<control>R", OpenResults,    0, "<StockItem>", GTK_STOCK_OPEN },
        { "/File/Open database file", NULL, OpenDatabaseFile, 0, "<Item>", },
        { "/File/Save databases", NULL, SaveDatabases, 0, "<Item>", },
        { "/File/sep1",     NULL,         NULL,    0, "<Separator>" },
        { "/File/Preferences", NULL, NULL, 0, "<Item>", },
        { "/File/sep2",     NULL,         NULL,    0, "<Separator>" },
        { "/File/recent1",  NULL,       OpenRecent, 1, "<Item>" }, 
        { "/File/recent2",  NULL,       OpenRecent, 2, "<Item>" }, 
        { "/File/recent3",  NULL,       OpenRecent, 3, "<Item>" }, 
        { "/File/recent4",  NULL,       OpenRecent, 4, "<Item>" }, 
        { "/File/sep3",     NULL,         NULL,    0, "<Separator>" },
        { "/File/_Quit",    "<CTRL>Q", quit, 0, "<StockItem>", GTK_STOCK_QUIT },
        { "/_Model",      NULL,         NULL,             0, "<Branch>" },
        { "/Model/_Solve",  NULL,         SolveModel,           0, "<Item>",  },
        { "/Model/_Results", NULL,  ShowOutputOptions, 0, "<Item>", }, 
        { "/Model/_Export summary",  NULL,     SolveModel,           1, "<Item>",  },
        { "/_Help",         NULL,         NULL,           0, "<Branch>" },
        { "/_Help/About",   NULL,         About,           0, "<Item>" },
    };


    // accel_group = gtk_accel_group_new ();

    item_factory = gtk_item_factory_new(GTK_TYPE_MENU_BAR, "<main>",
                                        NULL); // accel_group);

    gtk_item_factory_create_items (item_factory, 
                                   sizeof(menu_items)/sizeof(menu_items[0]), 
                                   menu_items, NULL);

    // gtk_window_add_accel_group (GTK_WINDOW (window), accel_group);

    recentWidgets[0] = gtk_item_factory_get_widget(item_factory, "/File/recent1");
    recentWidgets[1] = gtk_item_factory_get_widget(item_factory, "/File/recent2");
    recentWidgets[2] = gtk_item_factory_get_widget(item_factory, "/File/recent3");
    recentWidgets[3] = gtk_item_factory_get_widget(item_factory, "/File/recent4");

    gtk_widget_hide(recentWidgets[0]);
    gtk_widget_hide(recentWidgets[1]);
    gtk_widget_hide(recentWidgets[2]);
    gtk_widget_hide(recentWidgets[3]);

    solve_item = gtk_item_factory_get_widget(item_factory, "/Model/Solve");
    return gtk_item_factory_get_widget (item_factory, "<main>");
}

void
cableLog(const gchar *log_domain, GLogLevelFlags log_level,
         const gchar *message, gpointer unused_data)
{
    fprintf(stderr,"log message: %s\n", message);
}

int 
main(int argc, char *argv[]) 
{
    char recent[9];
    GtkWidget *w;
    GtkWidget *menubar;
    GdkColormap *colormap;
    gint   i, n;
    char **ptr;
    GError *err;
    char   filetype;
    char   docdir[256];

    g_log_set_always_fatal(G_LOG_LEVEL_ERROR);

    ptr = (char **) malloc(sizeof(char *) * argc);

    ptr[0] = argv[0]; 
    n = 1;
    for (i = 1 ; i < argc ; i++) {
        if (strncmp(argv[i], "-psn", 4)) {
            ptr[n++] = argv[i];
        }
    }
             
    argc = n;
    argv = ptr;

    gtk_init(&argc, &argv);
    // g_log_set_default_handler(cableLog, NULL);

    colormap = gdk_colormap_get_system();

    toplevel = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(toplevel), "WHOI Cable");
#ifdef NATIVE_OSX
    gtk_widget_set_usize(GTK_WIDGET(toplevel), 1120, 675);
#else
    gtk_widget_set_usize(GTK_WIDGET(toplevel), 1000, 630);
#endif

    gtk_signal_connect(GTK_OBJECT (toplevel), "destroy",
                       GTK_SIGNAL_FUNC (quit), NULL); 


    main_vbox = gtk_vbox_new(FALSE,1);
    gtk_container_set_border_width(GTK_CONTAINER(main_vbox),0); 
    gtk_container_add(GTK_CONTAINER(toplevel), main_vbox);
    gtk_widget_show(main_vbox);

    menubar = makemenu();
#ifdef NATIVE_OSX
    // ige_mac_menu_install_key_handler (); /* only for version 0.8 and later */
    ige_mac_menu_set_menu_bar (GTK_MENU_SHELL (menubar));
#else
    gtk_widget_show(menubar);
    gtk_box_pack_start(GTK_BOX(main_vbox), menubar, FALSE, TRUE, 0);
#endif

    notebook = gtk_notebook_new();
    gtk_notebook_set_tab_pos(GTK_NOTEBOOK(notebook), GTK_POS_TOP);
    gtk_box_pack_start(GTK_BOX(main_vbox), notebook, TRUE, TRUE, 0);
    gtk_widget_show(notebook);

    // we need to get some prefs stuff out of the way so that
    // cpp is setup before we try to build the database pages

    LoadPrefs();
    err = NULL;

#ifdef NATIVE_OSX
    // we don't have a real setup script under OSX so we have to figure out
    // pathnames here
    {
        char *cpp;
        struct stat statbuf;
        static char work2[256], work[256];
        unsigned int  sz;
        char         *execpath;
    
        // even under OSX someone might have set the prefs.ini key
        cpp = g_key_file_get_string(prefsFile, "parser", "cpp", &err);
        if (!cpp) { // more likely
            printf("no cpp prefs found\n");
            if (stat("/usr/bin/cpp", &statbuf) == 0) {
                printf("/usr/bin/cpp exists\n");
                cpp = NULL; // will default to /usr/bin/cpp in model/problem.c
            }
            else if (stat("/usr/local/bin/cpp", &statbuf) == 0) {
                printf("/usr/bin/cpp exists\n");
                cpp = "/usr/local/bin/cpp";
            }
            else {
                printf("searching for local Cable.app/cpp\n");
                sz = 256;
                _NSGetExecutablePath(work, &sz);
                execpath = dirname(work);
                sprintf(work2, "%s/../Resources/bin/cpp", execpath);
                if (stat(work2, &statbuf) == 0) {
                    cpp = work2;
                }
            }
        }
        SetupCpp(cpp, NULL, NULL);
    }
#else
    // the windows install script puts a pointer to cpp into the prefs.ini
    // file. If that setting exists SetupCpp will honor it - otherwise
    // it will use its regular fallback paths (no CPP under windows,
    // /usr/bin/cpp everywhere else
    SetupCpp(g_key_file_get_string(prefsFile, "parser", "cpp", &err), NULL, NULL);
#endif

    gtk_notebook_append_page(GTK_NOTEBOOK(notebook),
                             BuildMaterialsSheet(),
                             gtk_label_new(sheet_title[MATERIALS]));

    gtk_notebook_append_page(GTK_NOTEBOOK(notebook),
                             BuildConnectorsSheet(),
                             gtk_label_new(sheet_title[CONNECTORS]));

    gtk_notebook_append_page(GTK_NOTEBOOK(notebook),
                             BuildBuoysSheet(),
                             gtk_label_new(sheet_title[BUOYS]));

    gtk_notebook_append_page(GTK_NOTEBOOK(notebook),
                             BuildAnchorsSheet(),
                             gtk_label_new(sheet_title[ANCHORS]));

    gtk_notebook_prepend_page(GTK_NOTEBOOK(notebook),
                             (layout_hpane = BuildLayoutSheet()),
                             gtk_label_new(sheet_title[LAYOUT]));


    gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook), LAYOUT);

    gtk_widget_show(toplevel);

    w = gtk_table_get_child_at(layout_table, 0, layout_table -> ncols - 1);
    gtk_paned_set_position(GTK_PANED(layout_hpane), w -> allocation.x + w -> allocation.width + 80);

    // need to flush events so the GUI is fully up before   
    // we try to load any files
    while(gtk_events_pending()) {
        gtk_main_iteration();
    } 

    
    recentFiles = g_queue_new();
    
    for (i = 4 ; i >= 1 ; i--) {
        err = NULL;
        sprintf(recent, "recent%d", i);
        PushRecent(g_key_file_get_string(prefsFile, "files", recent, &err));
    }
    LatchRecent();

    if (argc == 2) {
        filetype = CheckMagic(argv[1]);
        if (filetype == -1) {
            SetCurrentFile(argv[1]); 
            NewModel(NULL, 0, NULL);
        }
        else {
            OpenCableFile(argv[1], filetype);
            // g_str_has_suffix(argv[1], ".crs") ? TRUE : FALSE);
        }
    }
    else {
        NewModel(NULL, 0, NULL);
#if defined (WINDOWS) || defined(NATIVE_OSX)
#ifdef WINDOWS
        GetEnvironmentVariable("userprofile", docdir, 256);        
#else
        sprintf(docdir, "%s/Documents", getenv("HOME"));
#endif
        if (g_file_test(docdir, G_FILE_TEST_IS_DIR)) {
            prev_file = strdup(docdir);
        }
#endif
    }

    gtk_main();

    return(0);
}
