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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gtk/gtk.h>
#include <gdk/gdk.h>
#include <gtk/gtkstock.h>
#include <gtkextra/gtkplot.h>
#include <gtkextra/gtkplotdata.h>
#include <gtkextra/gtkplotcanvas.h>
#include <gtkextra/gtkplotcanvasplot.h>
#include <gtkextra/gtkplotcanvastext.h>
#include <gtkextra/gtkplotps.h>
#include <gtkextra/gtkplotprint.h>
#include "plot.h"
#include "plotutil.h"



typedef struct {
    GtkWidget   *toplevel;
    int          destroy;
    void        *data;
    gdouble  *px;
    gdouble **py;
    gdouble  *pda;
    gchar   **labels;
    int       npoints;
    int       maxpoints;
    int       curr_idx;
    int       head_idx;
    int       running;
    GtkPlotData **dataset;
    GtkPlotData **currpt;
    GtkPlot     *plot;
    GtkWidget   *canvas;
    GtkWidget   *menu;
    int         style;
    int         num_sets;
    int         legend_available;
    int         legend_drawn;
    int         autoscroll_x;
    EntriesDialog *limdlg;
    GdkPixbuf *pixbuf;  // master pixbuf
    double     pixlim[4]; // lon, lat limits of master pixbuf
    double      peakring[10];
    int         peakdex;
    int         autoscale_x;
    int         autoscale_y;
    GtkWidget   *t_label;
} ChartDisplay;
#ifdef DEBUG
gint32 timer;
#endif

	/*
	 * the following two routines are based on the heuristics used in
	 * xmgr to figure out axis extreme and tick spacing ...
	 */

static double 
NiceNumber (double x, int round_mode)
{
   double	order;
   double	fraction;
   double	y;

   order = floor (log10(fabs(x)));
   fraction = fabs(x) / pow(10.0, order);

   if (round_mode) {
      if (fraction < 1.5)
         y = 1.0;
      else if (fraction < 3.0)
         y = 2.0;
      else if (fraction < 7.0)
         y = 5.0;
      else
         y = 10.0;
   }
   else if (fraction <= 1.0001)
      y = 1.0;
   else if (fraction <= 2.0001)
      y = 2.0;
   else if (fraction <= 5.0001)
      y = 5.0;
   else
      y = 10.0;


   if (x > 0)
      return y*pow(10.0, order);
   else
      return -y*pow(10.0, order); 

}

static int set_xrange(ChartDisplay *, double, double);
static int set_yrange(ChartDisplay *, double, double);

static void
render_background_image(ChartDisplay *cd)
{
    GdkPixmap *pixmap;
    GdkColormap *colormap;
    GdkPixbuf *pixbuf;
    GdkBitmap *mask;
    double lon_per_pix, lat_per_pix;
    int imwidth, imheight;
    int x0, y0;
    double xmin, xmax, ymin, ymax;

    if (!cd -> pixbuf) {
        return;
    }

    imwidth = gdk_pixbuf_get_width(cd -> pixbuf);
    imheight = gdk_pixbuf_get_height(cd -> pixbuf);

   fprintf(stderr,"master image width = %d, height = %d\n", imwidth, imheight);

    gtk_plot_get_xrange(cd -> plot, &xmin, &xmax);
    gtk_plot_get_yrange(cd -> plot, &ymin, &ymax);
    fprintf(stderr,"axes limits: %g %g %g %g\n", xmin, xmax, ymin, ymax);

    // we have to do a bounds check here as we might get called
    // to render before the display limits get updated (the update
    // limits dialog is not modal)

    if (xmin < cd -> pixlim[0])
        xmin = cd -> pixlim[0];
    if (xmax > cd -> pixlim[1])
        xmax = cd -> pixlim[1];
    if (ymin < cd -> pixlim[2])
        ymin = cd -> pixlim[2];
    if (ymax > cd -> pixlim[3])
        ymax = cd -> pixlim[3];

    // fprintf(stderr,"ymin, ymax = %g, %g\n", ymin, ymax);
    set_xrange(cd, xmin, xmax); // will return if unchanged
    set_yrange(cd, ymin, ymax); 

    fprintf(stderr, "axes limits: %g %g %g %g\n", xmin, xmax, ymin, ymax);
    fprintf(stderr, "pixlim: %g %g %g %g\n", cd->pixlim[0], cd->pixlim[1], cd->pixlim[2], cd->pixlim[3]);

    lon_per_pix = (cd -> pixlim[1] - cd -> pixlim[0]) / (double) imwidth; 
    lat_per_pix = (cd -> pixlim[3] - cd -> pixlim[2]) / (double) imheight; 

    x0 = (xmin - cd -> pixlim[0])/lon_per_pix;
    y0 = (cd -> pixlim[3] - ymax)/lat_per_pix;

    imwidth = (xmax - xmin)/lon_per_pix;
    imheight = (ymax - ymin)/lat_per_pix; 

    fprintf(stderr,"x0 = %d, y0 = %d\n", x0, y0);
    fprintf(stderr,"width = %d, height = %d\n", imwidth, imheight);

    pixbuf = gdk_pixbuf_new_subpixbuf(cd -> pixbuf, x0, y0, imwidth, imheight);
    gdk_pixbuf_render_pixmap_and_mask (pixbuf, &pixmap, &mask, 0);
    gdk_drawable_get_size(pixmap, &imwidth, &imheight);
    gtk_plot_set_background_pixmap(cd -> plot, pixmap);
    gtk_plot_canvas_set_size(GTK_PLOT_CANVAS(cd -> canvas), imwidth, imheight);
    return;
}

static void
setup_background_image(ChartDisplay *cd, char *fname)
{
    double pixlim[4];
    int    n;
   
    n = sscanf(fname, "%*[/A-Za-z0-9]_%lf_%lf_%lf_%lf.", 
               &pixlim[0], &pixlim[1], &pixlim[2], &pixlim[3]);
    
    if (n == 4) {
        cd -> pixlim[0] = pixlim[0];
        cd -> pixlim[1] = pixlim[1];
        cd -> pixlim[2] = pixlim[2];
        cd -> pixlim[3] = pixlim[3];
        fprintf(stderr,"setting master image limits from filename\n");
        fprintf(stderr,"%g %g %g %g\n", pixlim[0], pixlim[1], pixlim[2], pixlim[3]);
    }
    else {
        fprintf(stderr,"could not parse image extents from filename %s\n", fname);
        cd -> pixbuf = NULL;
        return;
    }
                 
    if (cd -> pixbuf) {
        g_object_unref(G_OBJECT(cd -> pixbuf));
    }
    cd -> pixbuf = gdk_pixbuf_new_from_file(fname, NULL);
    if (cd -> pixbuf) {
        g_object_ref(G_OBJECT(cd -> pixbuf));
    }
    return;
}

void
quitPlot(GtkWidget *button, gpointer data)
{
    ChartDisplay *cd = (ChartDisplay *) data;
    int           j;

#ifdef DEBUG
    gtk_timeout_remove(timer);
    gtk_main_quit();
#endif
    if (cd -> destroy == 0) {
        g_free(cd -> px);

        for (j = 0 ; j < cd -> num_sets ; j++)
            g_free(cd -> py[j]);

        g_free(cd -> py);
    }

    cd -> destroy = 1;
    gtk_widget_destroy(cd -> toplevel);
}

gboolean
deletePlot(GtkWidget *w, GdkEvent *event, gpointer data)
{ 
    quitPlot(w, data);
    return TRUE;
}

static void
clear(GtkWidget *button, gpointer data)
{
    ChartDisplay *cd = (ChartDisplay *) data;
    int i;
    int state;

    state = cd -> running;
    cd -> running = 0;

    cd -> npoints = 0;
    cd -> head_idx = cd -> curr_idx = 0;
    for (i = 0 ; i < cd -> num_sets ; i++) {
        gtk_plot_data_set_numpoints(cd -> dataset[i], cd -> npoints);
        gtk_plot_data_draw_points(cd -> dataset[i], 1);
    }

    gtk_plot_refresh(cd -> plot, NULL);

    gtk_plot_canvas_paint(GTK_PLOT_CANVAS(cd -> canvas));
    gtk_widget_queue_draw(GTK_WIDGET(cd -> canvas));


    cd -> running = state;
}

static void
my_pause(GtkWidget *button, gpointer data)
{
    ((ChartDisplay *) data) -> running = !gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(button));
}

static void popup_set_limits(GtkWidget *, gpointer *);

static void 
do_openbg(GtkWidget *dialog, gint response, gpointer *data)
{
    ChartDisplay *cd = (ChartDisplay *) data;
    char *filename;

    if (response == GTK_RESPONSE_DELETE_EVENT) 
        return;
    else if (response == GTK_RESPONSE_CANCEL) {
        gtk_widget_destroy(dialog);
        return;
    }

    filename = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (dialog));
    setup_background_image(cd, filename);
    if (cd -> pixbuf) {
        popup_set_limits(NULL, data);
        render_background_image(cd);
    }

    g_free(filename);

    gtk_widget_destroy (dialog);
}

static void
openbg(GtkWidget *button, gpointer data)
{
    ChartDisplay *cd = (ChartDisplay *) data;
    GtkWidget *dialog;

    dialog = gtk_file_chooser_dialog_new ("Open Background Image",
                      GTK_WINDOW(cd -> toplevel),
                      GTK_FILE_CHOOSER_ACTION_OPEN,
                      GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                      GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
                      NULL);

    // gtk_file_chooser_set_filename(GTK_FILE_CHOOSER (dialog), gtk_window_get_title(GTK_WINDOW(cd -> toplevel)));

    gtk_signal_connect(GTK_OBJECT(dialog), "response", GTK_SIGNAL_FUNC(do_openbg), cd);

    gtk_widget_show_all(dialog);
    return;
}

static void
toggle_legend(GtkWidget *item, gpointer data)
{
    ChartDisplay *cd = (ChartDisplay *) data;

    printf("toggling\n");

    if (cd -> legend_available && cd -> legend_drawn) {
        gtk_plot_hide_legends(cd -> plot);
        // gtk_plot_legends_set_attributes(cd -> plot, NULL, 0, NULL, NULL);
       // gtk_plot_set_legends_border(cd -> plot, GTK_PLOT_BORDER_NONE, 0);
        cd -> legend_drawn = 0;
    }
    else if (cd -> legend_available) {
        gtk_plot_show_legends(cd -> plot);
        // gtk_plot_legends_set_attributes(cd -> plot, NULL, 12, NULL, NULL);
        // gtk_plot_set_legends_border(cd -> plot, 2, 3);
        cd -> legend_drawn = 1;
    }
    gtk_plot_canvas_paint(GTK_PLOT_CANVAS(cd -> canvas));
    gtk_widget_queue_draw(GTK_WIDGET(cd -> canvas));
}

static int 
set_xrange(ChartDisplay *cd, double new_min, double new_max)
{
    double  xmin, xmax;
    int     precision;

    gtk_plot_get_xrange(cd -> plot, &xmin, &xmax);
    if (xmin == new_min && xmax == new_max) {
        return 0;
    }

    if (cd -> pixbuf) {
        if (new_min < cd -> pixlim[0])
            new_min = cd -> pixlim[0];
        if (new_max > cd -> pixlim[1])
            new_max = cd -> pixlim[1];
    }

    xmin = new_min;
    xmax = new_max;
 
    gtk_plot_set_xrange(cd -> plot, new_min, new_max);
    gtk_plot_set_ticks(cd -> plot, GTK_PLOT_AXIS_X, NiceNumber((xmax - xmin)/8.0, 1), 1);
    precision = (int) log10((xmax - xmin)/80);
    precision = precision >= 0 ? 0 : abs(precision);
    gtk_plot_axis_set_labels_style(gtk_plot_get_axis(cd -> plot, GTK_PLOT_AXIS_BOTTOM), GTK_PLOT_LABEL_FLOAT, precision);


    return 1; 
}

static int
set_yrange(ChartDisplay *cd, double new_min, double new_max)
{
    double  ymin, ymax;
    int     precision;
 
    gtk_plot_get_yrange(cd -> plot, &ymin, &ymax);
    if (ymin == new_min && ymax == new_max) {
        return 0;
    }

    if (cd -> pixbuf) {
        fprintf(stderr, "pixlim = %g %g %g %g\n", cd->pixlim[0], cd->pixlim[1], cd->pixlim[2], cd->pixlim[3]);
        
        if (new_max > cd -> pixlim[3])
            new_max = cd -> pixlim[3];
        if (new_min < cd -> pixlim[2])
            new_min = cd -> pixlim[2];
    }
 
    ymin = new_min;
    ymax = new_max;
 
    gtk_plot_set_yrange(cd -> plot, new_min, new_max);
    gtk_plot_set_ticks(cd -> plot, GTK_PLOT_AXIS_Y, NiceNumber((ymax - ymin)/8.0, 1), 1);
    precision = (int) log10((ymax - ymin)/80);
    precision = precision >= 0 ? 0 : abs(precision);
    gtk_plot_axis_set_labels_style(gtk_plot_get_axis(cd -> plot, GTK_PLOT_AXIS_LEFT), GTK_PLOT_LABEL_FLOAT, precision);

    return 1; 
}

static int
set_zrange(ChartDisplay *cd, double new_min, double new_max)
{
    double  zmin, zmax;

    zmin = cd->dataset[0]->gradient->ticks.min;
    zmax = cd->dataset[0]->gradient->ticks.max;

    if (zmin == new_min && zmax == new_max) {
        return 0;
    }
 
    gtk_plot_data_set_gradient(cd -> dataset[0], new_min, new_max, 3, 0);
    gtk_plot_data_set_gradient_size(cd -> dataset[0], cd -> plot -> internal_allocation.width);
    gtk_plot_data_move_gradient(cd -> dataset[0], 0, 1.13);

    return 1;
}

static gboolean
canvas_click(GtkWidget *widget, GdkEventButton *event, gpointer data)
{
    ChartDisplay *cd = (ChartDisplay *) data;
    GtkAllocation alloc;
    int plot_l, plot_r, plot_t, plot_b;
    int grad_l, grad_r, grad_t, grad_b;
    int x, y, b;
    float xzone, yzone;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    float dir;

    plot_l = cd -> plot -> internal_allocation.x;
    plot_r = plot_l + cd -> plot -> internal_allocation.width;
    plot_t = cd -> plot -> internal_allocation.y;
    plot_b = plot_t + cd -> plot -> internal_allocation.height;
   
    xzone = 0.2*(plot_r - plot_l); 
    yzone = 0.2*(plot_b - plot_t);

    if (cd -> style == CHART_XYZ) {
        alloc = gtk_plot_data_get_gradient_allocation(cd -> dataset[0]);
        grad_l = alloc.x;
        grad_r = grad_l + alloc.width;
        grad_t = alloc.y;
        grad_b = grad_t + alloc.height;
        zmin = cd->dataset[0]->gradient->ticks.min;
        zmax = cd->dataset[0]->gradient->ticks.max;
    }
    else {
        grad_t = grad_b = 10000;
        grad_l = grad_r = 0;
    }

    x = event -> x;
    y = event -> y;
    b = event -> button; 

    if (x > plot_l && x < plot_r && y > plot_t && y < plot_b && b == 3) { 
        printf("here\n");
        gtk_menu_popup (GTK_MENU(cd -> menu), NULL, NULL, NULL, NULL, 
                        event -> button, event -> time);
        return TRUE;
    }

    if (b == 1)
        dir = 1.0;
    else if (b == 3)
        dir = -1.0;
    else
        goto DONE;

    gtk_plot_get_xrange(cd -> plot, &xmin , &xmax);
    gtk_plot_get_yrange(cd -> plot, &ymin , &ymax);

    if (x < plot_l && y < plot_t + yzone && y > plot_t - yzone) {
        set_yrange(cd, ymin, ymax + dir*(ymax - ymin)*0.05);
    }
    else if (x < plot_l && y < plot_b + yzone && y > plot_b - yzone) {
        set_yrange(cd, ymin + dir*(ymax - ymin)*0.05, ymax);
    }
    else if (y > plot_b && y < grad_t && x > plot_l - xzone && x < plot_l + xzone) {
        set_xrange(cd, xmin + dir*(xmax - xmin)*0.05, xmax);
    }
    else if (y > plot_b && y < grad_t && x > plot_r - xzone && x < plot_r + xzone) {
        set_xrange(cd, xmin, xmax + dir*(xmax - xmin)*0.05);
    }
    else if (y > grad_t && y < grad_b && x > grad_l - xzone && x < grad_l + xzone) {
        set_zrange(cd, zmin + dir*0.05*(zmax - zmin), zmax);
    }
    else if (y > grad_t && y < grad_b && x > grad_r - xzone && x < grad_r + xzone) {
        set_zrange(cd, zmin, zmax + dir*0.05*(zmax - zmin));
    }
    else {
        goto DONE;
    }

    gtk_plot_canvas_paint(GTK_PLOT_CANVAS(cd -> canvas));
    gtk_widget_queue_draw(GTK_WIDGET(cd -> canvas));
    if (cd -> pixbuf) {
        render_background_image(cd);
    }

DONE:
    return TRUE;
}

static void
resize(GtkWidget *canvas, gpointer allocation, gpointer data)
{
   gdouble h, w;
   ChartDisplay *cd = (ChartDisplay *) data;

   // printf("%d %d\n", cd -> canvas -> allocation.width, cd -> canvas -> allocation.height);
   gtk_plot_canvas_set_size(GTK_PLOT_CANVAS(cd -> canvas), cd -> canvas -> allocation.width, canvas -> allocation.height);
 // gtk_plot_data_set_gradient_size(cd -> dataset, cd -> canvas -> allocation.width);
    if (cd -> style == CHART_XYZ)
        gtk_plot_data_set_gradient_size(cd -> dataset[0], cd -> plot -> internal_allocation.width);

   gtk_plot_canvas_paint(GTK_PLOT_CANVAS(cd -> canvas));
   gtk_widget_queue_draw(GTK_WIDGET(cd -> canvas));
}

static void
play(GtkWidget *button, gpointer data)
{
    ((ChartDisplay *) data) -> running = 1;
}

static void
draw_page(GtkPrintOperation *operation, GtkPrintContext *context, int page_nr, gpointer data) 
{
    ChartDisplay    *cd = (ChartDisplay *) data;
    GtkPrintSettings *settings;
    int        width, height;
    GdkPixbuf   *pix;
    cairo_t     *cr;
    cairo_surface_t *surf;
    cairo_matrix_t trans;
    GtkPaperSize *size;
    double        paper_width, paper_height;
    GtkPageSetup *setup;

    setup = gtk_print_operation_get_default_page_setup(operation);
    paper_width = gtk_page_setup_get_paper_width(setup, GTK_UNIT_INCH);
    paper_height = gtk_page_setup_get_paper_height(setup, GTK_UNIT_INCH);

    fprintf(stderr,"w = %g, h = %g\n", paper_width, paper_height);

    width = GTK_PLOT_CANVAS(cd -> canvas) -> pixmap_width; 
    height = GTK_PLOT_CANVAS(cd -> canvas) -> pixmap_height;

    pix = gdk_pixbuf_get_from_drawable(NULL, GTK_PLOT_CANVAS(cd -> canvas) -> pixmap, NULL, 0, 0, 0, 0, width, height);

    gdk_pixbuf_save(pix, "foo.png", "png", NULL, NULL);


    surf = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
    cr = gtk_print_context_get_cairo_context(context);

    if ((double) width / (double) height > paper_width/paper_height) {
        cairo_scale(cr, paper_width/width, paper_width/width);
    }
    else {
        cairo_scale(cr, paper_height/height, paper_height/height);
    }
 
    gdk_cairo_set_source_pixbuf(cr, pix, 0.0, 0.0);

    // cairo_matrix_init_scale(&trans, 800.0/width, 600.0/height);

    cairo_paint(cr);
    g_object_unref(G_OBJECT(pix));
}

#ifndef WINDOWS

extern void ExportResultMovie(DisplayObject);

static void
movie(GtkWidget *button, gpointer data)
{
    ExportResultMovie((DisplayObject) data);
    return;
}

#endif

static void
print(GtkWidget *button, gpointer data)
{
    ChartDisplay *cd = (ChartDisplay *) data;
    GtkWidget *dialog;
    GtkPageSetup *setup;

    GtkPrintOperation *op;
    GtkPrintOperationResult res; 

    setup = gtk_page_setup_new();
    op = gtk_print_operation_new (); 
    // gtk_print_operation_set_print_settings (op, settings); 
 
///    gtk_print_run_page_setup_dialog(cd ->  toplevel, setup, op);
    gtk_print_operation_set_default_page_setup(op, setup);

    gtk_print_operation_set_n_pages (op, 1); 
    gtk_print_operation_set_unit (op, GTK_UNIT_INCH); 
    // gtk_print_operation_set_use_full_page(op, TRUE);
    g_signal_connect (op, "draw_page", G_CALLBACK (draw_page), cd); 
    gtk_print_operation_set_default_page_setup(op, gtk_page_setup_new());
    gtk_print_operation_set_print_settings(op, gtk_print_settings_new());
    res = gtk_print_operation_run (op, GTK_PRINT_OPERATION_ACTION_PRINT_DIALOG, 
                                   NULL, NULL); 


    printf("assigned draw_page\n");
     
#if 0
    dialog = gtk_file_chooser_dialog_new ("Print to file",
                      GTK_WINDOW(cd -> toplevel),
                      GTK_FILE_CHOOSER_ACTION_SAVE,
                      GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                      GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT,
                      NULL);

    gtk_file_chooser_set_filename(GTK_FILE_CHOOSER (dialog), gtk_window_get_title(GTK_WINDOW(cd -> toplevel)));

    if (gtk_dialog_run(GTK_DIALOG (dialog)) == GTK_RESPONSE_ACCEPT) {
        char *filename;

        filename = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (dialog));
        // gtk_plot_canvas_export_ps(GTK_PLOT_CANVAS(((ChartDisplay *) data) -> canvas), 
          //                     filename, 0, 0, GTK_PLOT_LETTER);
        gtk_plot_export_ps(cd -> plot, filename, GTK_PLOT_PORTRAIT, TRUE, GTK_PLOT_LETTER);

        g_free(filename);
    }

    gtk_widget_destroy (dialog);
#endif
}

static void
scroll(ChartDisplay *cd, int dir)
{
    gdouble xmin, xmax;

    if (cd -> running) {
        return;
    }

    gtk_plot_get_xrange(cd -> plot, &xmin , &xmax);
    gtk_plot_set_xrange(cd -> plot, xmin + dir*(xmax - xmin)*0.2, 
                                    xmax + dir*(xmax - xmin)*0.2);
    gtk_plot_canvas_paint(GTK_PLOT_CANVAS(cd -> canvas));
    gtk_widget_queue_draw(GTK_WIDGET(cd -> canvas));
}


static void
scrollback(GtkWidget *button, gpointer data)
{
    scroll((ChartDisplay *) data, -1);
}

static void
scrollfwd(GtkWidget *button, gpointer data)
{
    scroll((ChartDisplay *) data, 1);
}




#define HALFBUFF 20000
#define FULLBUFF 40000

void
GetChartDisplayLimits(DisplayObject obj, double *limits)
{
    ChartDisplay *cd = (ChartDisplay *) obj;

    gtk_plot_get_xrange(cd -> plot, &limits[0] , &limits[1]);
    gtk_plot_get_yrange(cd -> plot, &limits[2] , &limits[3]);
    if (cd -> style == CHART_XYZ) {
        limits[4] = cd->dataset[0]->gradient->ticks.min;
        limits[5] = cd->dataset[0]->gradient->ticks.max;
    }

    return;
}

static gchar **
CopyLabels(gchar **labels, int n)
{
    int     i;
    gchar **copy;

    copy = (gchar **) g_malloc(sizeof(gchar *) * n);

    for (i = 0 ; i < n ; i++) {
        if (labels[i])
            copy[i] = strdup(labels[i]);
        else
            copy[i] = NULL;
    }

    return copy;
}

static double
max(double *x, int n)
{
    double  y;
    int     i;

    y = fabs(x[0]);
    for (i = 1 ; i < n ; i++) {
        if (fabs(x[i]) > y)
            y = fabs(x[i]);
    }

    return y;
}

double *
ChartX(DisplayObject obj)
{
    ChartDisplay *cd = (ChartDisplay *) obj;
    return cd -> px;
}

double **
ChartY(DisplayObject obj)
{
    ChartDisplay *cd = (ChartDisplay *) obj;
    return cd -> py;
}


GdkPixbuf *
ChartPixbuf(DisplayObject obj, GdkPixbuf *dest)
{
    ChartDisplay *cd = (ChartDisplay *) obj;
    int           width, height;
    GdkPixbuf    *pix;

    width = GTK_PLOT_CANVAS(cd -> canvas) -> pixmap_width; 
    height = GTK_PLOT_CANVAS(cd -> canvas) -> pixmap_height;

    pix = gdk_pixbuf_get_from_drawable(dest, GTK_PLOT_CANVAS(cd -> canvas) -> pixmap, NULL, 0, 0, 0, 0, width, height);
    return pix;
}

void 
UpdateChartAnimation(DisplayObject obj, char *label)
{
    GtkPlotAxis *axis;
    ChartDisplay *cd = (ChartDisplay *) obj;
    int           i, j;
    double        xmin, xmax, ymin, ymax;
    
    for (i = 0 ; i < cd -> num_sets ; i++) {
        gtk_plot_data_set_numpoints(cd -> dataset[i], cd -> npoints);
        gtk_plot_data_set_x(cd -> dataset[i], cd -> px);
        gtk_plot_data_set_y(cd -> dataset[i], cd -> py[i]);
    }

    for (i = 0 ; i < cd -> num_sets ; i++) {
        gtk_plot_data_draw_points(cd -> dataset[i], 1);
    }

    if (cd -> autoscale_x || cd -> autoscale_y)  {
        gtk_plot_get_xrange(cd -> plot, &xmin, &xmax);
        gtk_plot_get_yrange(cd -> plot, &ymin, &ymax);
        if (xmin == xmax) {
            xmin = xmax = cd -> px[0];
        }
        if (ymin == ymax) {
            ymin = ymax = cd -> py[0][0];
        }

        for (i = 0 ; i < cd -> npoints ; i++) {
            if (cd -> px[i] > xmax)
                xmax = cd -> px[i];
            else if (cd -> px[i] < xmin)
                xmin = cd -> px[i];
            
            for (j = 0 ; j < cd -> num_sets ; j++) {
                if (cd -> py[j][i] > ymax)
                    ymax = cd -> py[j][i];
                else if (cd -> py[j][i] < ymin)
                    ymin = cd -> py[j][i];
            }
        }

        set_xrange(cd, xmin, xmax);
        set_yrange(cd, ymin, ymax);
        
        // gtk_plot_autoscale(cd -> plot);
    }


    if (label) {
        axis = gtk_plot_get_axis(cd -> plot, GTK_PLOT_AXIS_BOTTOM);
        gtk_plot_axis_set_title(axis, label);
        fprintf(stderr,"%s\n", label);
    }
    if (label && cd -> t_label) {
        gtk_label_set_text(GTK_LABEL(cd -> t_label), label);
        fprintf(stderr,"%s\n", label);
    }

    gtk_plot_refresh(cd -> plot, NULL);

    gtk_plot_canvas_paint(GTK_PLOT_CANVAS(cd -> canvas));
    gtk_widget_queue_draw(GTK_WIDGET(cd -> canvas));
}

void
PlotChartDisplay(DisplayObject obj, gdouble *x, gdouble **y)
{
    ChartDisplay *cd = (ChartDisplay *) obj;
    double  ymin, ymax, xmin, xmax;
    int     j, i;
    int     updated;

    for (j = 0 ; j < cd -> maxpoints ; j++) {
        cd -> px[j] = x[j];
        for (i = 0 ; i < cd -> num_sets ; i++) {
            cd -> py[i][j] = y[i][j];
        }
    }

    for (i = 0 ; i < cd -> num_sets ; i++) { 
        gtk_plot_data_set_numpoints(cd -> dataset[i], cd -> maxpoints);
        gtk_plot_data_set_x(cd -> dataset[i], cd -> px);
        gtk_plot_data_set_y(cd -> dataset[i], cd -> py[i]);
        gtk_plot_data_draw_points(cd -> dataset[i], 1);
    }

    updated = 0;

    if (cd -> autoscale_y) {
        gtk_plot_get_yrange(cd -> plot, &ymin, &ymax);
        if (ymin == ymax) {
            ymin = ymax = cd -> py[0][0];
        }
        for (i = 0 ; i < cd -> maxpoints ; i++) {
            for (j = 0 ; j < cd -> num_sets ; j++) {
                if (cd -> py[j][i] > ymax)
                    ymax = cd -> py[j][i];
                else if (cd -> py[j][i] < ymin)
                    ymin = cd -> py[j][i];
            }
        }

        updated = set_yrange(cd, ymin, ymax);
    }


    gtk_plot_refresh(cd -> plot, NULL);
    
    if (updated) {
        gtk_plot_canvas_paint(GTK_PLOT_CANVAS(cd -> canvas));
        gtk_widget_queue_draw(GTK_WIDGET(cd -> canvas));
    }

    return;
}

void
UpdateChartDisplay(DisplayObject obj, gdouble x, gdouble *y, gdouble z)
{
    int     updated;
    ChartDisplay *cd = (ChartDisplay *) obj;
//    int     mirror_idx
    int i, j;
    double  ymin, ymax, xmin, xmax;
    int     label_it = 0;
    double  lastpt1, lastpt2;
    GtkPlotCanvasChild *child;
    gchar **labels;
    char    buffer[32];

/*
    if (cd -> curr_idx >= HALFBUFF) 
        mirror_idx = cd -> curr_idx - HALFBUFF;
    else
        mirror_idx = cd -> curr_idx + HALFBUFF;
*/

    for (i = 0 ; i < cd -> num_sets ; i++) {
        // cd -> py[i][cd -> curr_idx] = cd -> py[i][mirror_idx] = y[i];
        cd -> py[i][cd -> curr_idx] = y[i];
    }

    // cd -> px[cd -> curr_idx] = cd -> px[mirror_idx] = x;
    cd -> px[cd -> curr_idx] = x;

    if (cd -> style == CHART_XYZ) {
        // cd -> pda[cd -> curr_idx] = cd -> pda[mirror_idx] = z;
        cd -> pda[cd -> curr_idx] = z;
    }


    cd -> curr_idx++;
    cd -> npoints ++;
   
 
    for (i = 0 ; i < cd -> num_sets ; i++) {
        gtk_plot_data_set_numpoints(cd -> dataset[i], cd -> npoints);
        gtk_plot_data_set_x(cd -> dataset[i], cd -> px); 
        gtk_plot_data_set_y(cd -> dataset[i], cd -> py[i]);
        if (cd -> style == CHART_XYZ && i == 0) {
            gtk_plot_data_set_da(cd -> dataset[0], cd -> pda); 
        }
    }


    if(!cd -> running)
        return;

    if (cd -> autoscroll_x) {
        gtk_plot_get_xrange(cd -> plot, &xmin , &xmax);

        if(x > xmax){
        // gtk_plot_set_xrange(cd -> plot, xmin + 5. , xmax + 5.);
            if (xmax - xmin <= 1.0)
                gtk_plot_set_xrange(cd -> plot, x + 0.25 - (xmax - xmin),  x + 0.25);
            else
                gtk_plot_set_xrange(cd -> plot, x + 0.5 - (xmax - xmin),  x + 0.5);
        // gtk_plot_set_xrange(cd -> plot, xmin + 0.25*(xmax - xmin), xmax + 0.25*(xmax - xmin));
        // gtk_plot_set_xrange(cd -> plot, x - (xmax - xmin), x);
            gtk_plot_canvas_paint(GTK_PLOT_CANVAS(cd -> canvas));
            gtk_widget_queue_draw(GTK_WIDGET(cd -> canvas));
        } 
    }

    for(i = 0 ; i < cd -> num_sets ; i++) {
        gtk_plot_data_draw_points(cd -> dataset[i], 1);
        if (cd -> currpt) {
            gtk_plot_data_draw_points(cd -> currpt[i], 1);
        }
    }

    updated = 0;
    if (cd -> autoscale_y) {
        gtk_plot_get_yrange(cd -> plot, &ymin, &ymax);
        if (ymin == ymax) {
            ymin = ymax = cd -> py[0][0];
        }
        for (i = 0 ; i < cd -> npoints ; i++) {
            for (j = 0 ; j < cd -> num_sets ; j++) {
                if (cd -> py[j][i] > ymax)
                    ymax = cd -> py[j][i];
                else if (cd -> py[j][i] < ymin)
                    ymin = cd -> py[j][i];
            }
        }

        updated = set_yrange(cd, ymin, ymax);
    }

    gtk_plot_refresh(cd -> plot, NULL);
    
    // no need to do this unless we're a scatter or profile where we
    // constantly move a single cross-hair to identify the most recent point

    if (cd -> currpt || updated) {
        gtk_plot_canvas_paint(GTK_PLOT_CANVAS(cd -> canvas));
        gtk_widget_queue_draw(GTK_WIDGET(cd -> canvas));
    }

    return;
}


static double
increment(double min, double max)
{
    double  range, d;
    int     numticks = 8;
    double  nice_min, nice_max;

    range = NiceNumber((max - min), 0);
    d = NiceNumber(range / (numticks + 1.0), 1);

    nice_min = floor(min / d) * d;
    nice_max = ceil(max / d) * d;

    return d;
}

static void
toggle_set(GtkWidget *w, gpointer data)
{
   ChartDisplay *cd = (ChartDisplay *) data;
   int           set_number;

   if (sscanf(gtk_label_get_text(GTK_LABEL(GTK_BIN(w) -> child)), 
              "toggle set %d", &set_number) != 1) {
        return;
   }

   set_number --; // zero-offset
   if (set_number < 0 || set_number > cd -> num_sets -1) {
        return;
   }

   if (!GTK_WIDGET_VISIBLE(GTK_WIDGET(cd -> dataset[set_number]))) {
      gtk_widget_show(GTK_WIDGET(cd -> dataset[set_number]));
      if (cd -> currpt && cd -> currpt[set_number])
          gtk_widget_show(GTK_WIDGET(cd -> currpt[set_number]));
   }
   else {
      gtk_widget_hide(GTK_WIDGET(cd -> dataset[set_number]));
      if (cd -> currpt && cd -> currpt[set_number])
          gtk_widget_hide(GTK_WIDGET(cd -> currpt[set_number]));
   }

   gtk_plot_canvas_paint(GTK_PLOT_CANVAS(cd -> canvas));
   gtk_widget_queue_draw(GTK_WIDGET(cd -> canvas));
}

static void do_set_limits(GtkWidget *, gint, gpointer *);

static void
makelimdlg(ChartDisplay *cd)
{
    char buff[64];
    char *limit_fields[] = {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax", NULL};

    sprintf(buff, "%s limits", gtk_window_get_title(GTK_WINDOW(cd -> toplevel)));
    if (cd -> style == CHART_XYZ) {
        cd -> limdlg = CreateEntriesDialog(cd -> toplevel, buff, limit_fields);
    }
    else {
        limit_fields[4] = NULL;        
        cd -> limdlg = CreateEntriesDialog(cd -> toplevel, buff, limit_fields);
    }

    gtk_signal_connect(GTK_OBJECT(cd -> limdlg -> dialog), "response", GTK_SIGNAL_FUNC(do_set_limits), cd);
}

static void
do_set_limits(GtkWidget *w, gint response, gpointer *data)
{
    ChartDisplay *cd = (ChartDisplay *) data;
    double        limits[6];
    char         *ptr1, *ptr2;
    int           i, nlim;
    int           errors, changes;
    GdkColor      red, black;

    gdk_color_parse("red", &red);
    gdk_color_parse("black", &black);

    if (response == GTK_RESPONSE_DELETE_EVENT) {
        cd -> limdlg = NULL;
        return;
    }

    if (response == GTK_RESPONSE_CANCEL) {
        gtk_widget_hide(cd -> limdlg -> dialog);
        return;
    }

    nlim = cd -> style == CHART_XYZ ? 6 : 4;
    errors = 0;
    for (i = 0 ; i < nlim ; i++) {
        ptr1 = (char *) gtk_entry_get_text(GTK_ENTRY(cd -> limdlg -> field[i]));
        limits[i] = strtod(ptr1, &ptr2);
        if (ptr1 == ptr2 || *ptr2 != 0) {
            gtk_widget_modify_text(cd -> limdlg -> field[i], GTK_STATE_NORMAL, &red);
            errors ++;
            
        }   
    }
 
    if (!errors) { 
        if (limits[0] >= limits[1]) {
            gtk_widget_modify_text(cd -> limdlg -> field[0], GTK_STATE_NORMAL, &red);
            gtk_widget_modify_text(cd -> limdlg -> field[1], GTK_STATE_NORMAL, &red);
            errors ++;
        }
        if (limits[2] >= limits[3]) {
            gtk_widget_modify_text(cd -> limdlg -> field[2], GTK_STATE_NORMAL, &red);
            gtk_widget_modify_text(cd -> limdlg -> field[3], GTK_STATE_NORMAL, &red);
            errors ++;
        }
        if (cd -> style == CHART_XYZ && limits[4] >= limits[5]) {
            gtk_widget_modify_text(cd -> limdlg -> field[4], GTK_STATE_NORMAL, &red);
            gtk_widget_modify_text(cd -> limdlg -> field[5], GTK_STATE_NORMAL, &red);
            errors ++;
        }
    } 

    if (errors) {
        return;
    }

    changes = set_xrange(cd, limits[0], limits[1]);
    changes += set_yrange(cd, limits[2], limits[3]);

    if (cd -> style == CHART_XYZ) {
        changes += set_zrange(cd, limits[4], limits[5]);

    }

    if (changes) {
        gtk_plot_canvas_paint(GTK_PLOT_CANVAS(cd -> canvas));
        gtk_widget_queue_draw(GTK_WIDGET(cd -> canvas));
    }
    
    for (i = 0 ; i < nlim ; i++) {
        gtk_widget_modify_text(cd -> limdlg -> field[i], GTK_STATE_NORMAL, &black);
    }
 
    gtk_widget_hide(cd -> limdlg -> dialog);
}

static void
popup_set_limits(GtkWidget *w, gpointer *data)
{
    ChartDisplay *cd = (ChartDisplay *) data;
    double        limits[6];
    char          buff[32];
    int           nlim, i;

    if (cd -> limdlg == NULL) {
        makelimdlg(cd);
    }

    nlim = 4;
    gtk_plot_get_xrange(cd -> plot, &limits[0] , &limits[1]);
    gtk_plot_get_yrange(cd -> plot, &limits[2] , &limits[3]);

    if (cd -> style == CHART_XYZ) {
        limits[4] = cd->dataset[0]->gradient->ticks.min;
        limits[5] = cd->dataset[0]->gradient->ticks.max;
        nlim = 6;
    }

    for (i = 0 ; i < nlim ; i++) {
        sprintf(buff, "%g", limits[i]); 
        gtk_entry_set_text(GTK_ENTRY(cd -> limdlg -> field[i]), buff);
    }
    gtk_widget_show_all(cd -> limdlg -> dialog);

    return;
}

static GtkWidget *
makemenu(ChartDisplay *cd, char **legends, int num_sets)
{
    GtkWidget   *menu;
    GtkWidget   *menuitem;
    int          i;
    char         buff[64]; 

    menu = gtk_menu_new();

    if (cd -> legend_available) {
        menuitem = gtk_menu_item_new_with_label("toggle legend");
        gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);
        gtk_signal_connect (GTK_OBJECT(menuitem), "activate",
                            GTK_SIGNAL_FUNC (toggle_legend), cd);
        gtk_widget_show(menuitem);
    }

    for (i = 0 ; i < num_sets ; i++) {
        sprintf(buff, "toggle set %d", i+1);
        menuitem = gtk_menu_item_new_with_label(buff);
        gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);
        gtk_signal_connect (GTK_OBJECT(menuitem), "activate",
                            GTK_SIGNAL_FUNC(toggle_set), cd);
        gtk_widget_show(menuitem);
    }

    menuitem = gtk_menu_item_new_with_label("set limits");
    gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);
    gtk_signal_connect (GTK_OBJECT(menuitem), "activate",
                        GTK_SIGNAL_FUNC(popup_set_limits), cd);
    gtk_widget_show(menuitem);

    return menu;
}

DisplayObject
CreateChartDisplay(char *title, int style, int num_sets, char *xlabel, char *ylabel, char  **legends, char *bg_image, double *params, int autoscroll, int num_points, void *userData)
{
    int i;
    GtkWidget *vbox1, *hbox;
    GtkWidget *toolbar;
    GtkWidget *button;
    GtkPlotAxis *gradient;
    GtkPlotAxis *axis;
    GtkPlotCanvasChild *child;
    GdkColor color, grad_min, grad_max;
    ChartDisplay  *cd;
    char *set_colors[] = {"red", "blue", "green", "cyan", "magenta", "orange"};
    char  buff[64];
    double ymin, ymax;

    cd = (ChartDisplay *) malloc(sizeof(ChartDisplay));
    cd -> autoscroll_x = autoscroll;
    cd -> num_sets = num_sets;
    cd -> style = style; 
    cd -> data = userData;

    if (num_points) {
        cd -> maxpoints = num_points;
    }
    else {
        cd -> maxpoints = FULLBUFF;
    }

    cd -> px = (gdouble *) g_malloc(cd -> maxpoints*sizeof(gdouble));

    cd -> py = (gdouble **) g_malloc(sizeof(gdouble *) * num_sets);
    for (i = 0 ; i < num_sets ; i++)
        cd -> py[i] = (gdouble *) g_malloc(cd -> maxpoints*sizeof(gdouble));

    if (style == CHART_XYZ)
        cd -> pda = (gdouble *) g_malloc(cd -> maxpoints*sizeof(gdouble));
    else
        cd -> pda = NULL;

    if (style == CHART_PEAK) {
        cd -> labels = (gchar **) g_malloc(cd -> maxpoints*sizeof(gchar *));
        cd -> peakdex = 0;
        for (i = 0 ; i < 10 ; i++) {
            cd -> peakring[i] = 0;
        }
    }
    else
        cd -> labels = NULL;

    cd -> head_idx = cd -> curr_idx = 0;
    if (style == CHART_ANIMATE) {
        cd -> npoints = cd -> maxpoints;
    }
    else {
        cd -> npoints = 0;
    }

    cd -> destroy = 0;
    cd -> toplevel = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(cd ->  toplevel), title);
    // gtk_widget_set_usize(cd ->  toplevel, 600, 550);
    // SIZE
    gtk_window_set_default_size(GTK_WINDOW(cd -> toplevel), 600, 500);
    gtk_container_border_width(GTK_CONTAINER(cd -> toplevel), 0);

    gtk_widget_show(cd -> toplevel);

    gtk_signal_connect (GTK_OBJECT(cd -> toplevel), "delete-event", 
                     GTK_SIGNAL_FUNC (deletePlot), cd);

    vbox1=gtk_vbox_new(FALSE,0);
    gtk_container_add(GTK_CONTAINER(cd -> toplevel),vbox1);
    gtk_widget_show(vbox1);

    // SIZE
    cd -> canvas = gtk_plot_canvas_new(600, 468, 1.);

    GTK_PLOT_CANVAS_UNSET_FLAGS(GTK_PLOT_CANVAS(cd -> canvas), GTK_PLOT_CANVAS_DND_FLAGS);

    gdk_color_parse("light grey", &color);
    gdk_color_alloc(gtk_widget_get_colormap(cd -> canvas), &color);
    gtk_plot_canvas_set_background(GTK_PLOT_CANVAS(cd -> canvas), &color);

    gtk_widget_show(cd -> canvas);

    cd -> plot = GTK_PLOT(gtk_plot_new_with_size(NULL, 1.0, 1.0));

    toolbar = gtk_toolbar_new();
    gtk_toolbar_set_style(GTK_TOOLBAR(toolbar), GTK_TOOLBAR_ICONS);
/*
    toolbar_insert_stock(toolbar, GTK_STOCK_MEDIA_PAUSE, 
                         GTK_SIGNAL_FUNC(my_pause), cd, 1, NULL);
    toolbar_insert_stock(toolbar, GTK_STOCK_MEDIA_REWIND, 
                         GTK_SIGNAL_FUNC(scrollback), cd, 0, NULL);
    toolbar_insert_stock(toolbar, GTK_STOCK_MEDIA_FORWARD, 
                         GTK_SIGNAL_FUNC(scrollfwd), cd, 0, NULL);
*/
    toolbar_insert_stock(toolbar, GTK_STOCK_ZOOM_FIT,
                         GTK_SIGNAL_FUNC(popup_set_limits), cd, 0, NULL);
    toolbar_insert_stock(toolbar, GTK_STOCK_PRINT, 
                         GTK_SIGNAL_FUNC(print), cd, 0, NULL);

#ifndef WINDOWS
    if (style == CHART_ANIMATE) {
        toolbar_insert_stock(toolbar, GTK_STOCK_EXECUTE,
                             GTK_SIGNAL_FUNC(movie), cd, 0, NULL);
    }
#endif

/*
    toolbar_insert_stock(toolbar, GTK_STOCK_DELETE, 
                         GTK_SIGNAL_FUNC(clear), cd, 0, NULL);
*/
    if (bg_image) {
        toolbar_insert_stock(toolbar, GTK_STOCK_OPEN, 
                             GTK_SIGNAL_FUNC(openbg), cd, 0, NULL);
    }
    // toolbar_insert_stock(toolbar, GTK_STOCK_QUIT, 
    //                      GTK_SIGNAL_FUNC(quitPlot), cd, 0, NULL);

    hbox=gtk_hbox_new(FALSE,4);
    gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(toolbar), TRUE, TRUE, 0);

    if (style == CHART_ANIMATE) { 
        cd -> t_label = gtk_label_new("          ");
        gtk_widget_set_size_request(GTK_WIDGET(cd -> t_label), 100, 40);
        gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(cd -> t_label), FALSE, FALSE, 0);
    }
    else {
        cd -> t_label = NULL;
    }

    gtk_widget_show_all(GTK_WIDGET(hbox));
    
    gtk_box_pack_start(GTK_BOX(vbox1), GTK_WIDGET(hbox), FALSE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(vbox1), cd -> canvas, TRUE, TRUE, 0);

    gdk_color_parse("light yellow", &color);
    gdk_color_alloc(gtk_widget_get_colormap(GTK_WIDGET(cd -> plot)), &color);
    gtk_plot_set_background(cd -> plot, &color);

    gdk_color_parse("white", &color);
    gdk_color_alloc(gtk_widget_get_colormap(cd -> canvas), &color);
 // gtk_plot_legends_set_attributes(cd -> plot,
 //                                 NULL, 0, NULL, &color);

    cd -> pixbuf = NULL;
    if (bg_image) {
        setup_background_image(cd, bg_image); // do this so cd -> pixlim 
                                              // are set so that the range
                                              // sets below can honor them
    }

    cd -> autoscale_x = cd -> autoscale_y = 1;

    if (params) {
        if (params[1] > params[0]) {
            set_xrange(cd, params[0], params[1]);
            cd -> autoscale_x = 0;
        }

        if (params[3] > params[2]) {
            set_yrange(cd, params[2], params[3]);
            cd -> autoscale_y = 0;
        }
    }

    if (cd -> pixbuf) {
        render_background_image(cd);
    }

    gtk_plot_get_yrange(cd -> plot, &ymin, &ymax);
    
    gtk_plot_axis_set_visible(gtk_plot_get_axis(cd -> plot, GTK_PLOT_AXIS_TOP), FALSE);
    gtk_plot_axis_set_visible(gtk_plot_get_axis(cd -> plot, GTK_PLOT_AXIS_RIGHT), FALSE);
    gtk_plot_grids_set_visible(cd -> plot, TRUE, TRUE, TRUE, TRUE);

    child = gtk_plot_canvas_plot_new(cd -> plot);
    // SIZE
    gtk_plot_canvas_put_child(GTK_PLOT_CANVAS(cd -> canvas), child, .15, .08, .95, 0.9); // 0.77
    gtk_widget_show(GTK_WIDGET(cd -> plot));

    gtk_plot_axis_hide_title(gtk_plot_get_axis(cd -> plot, GTK_PLOT_AXIS_TOP));
    gtk_plot_axis_hide_title(gtk_plot_get_axis(cd -> plot, GTK_PLOT_AXIS_RIGHT));
    if (ylabel) {
        axis = gtk_plot_get_axis(cd -> plot, GTK_PLOT_AXIS_LEFT);
        gtk_plot_axis_set_title(axis, ylabel);
        gtk_plot_axis_move_title(axis, axis->title.angle, axis->title.x+0.035, axis->title.y);
    }

    if (xlabel) {
        axis = gtk_plot_get_axis(cd -> plot, GTK_PLOT_AXIS_BOTTOM);
        gtk_plot_axis_set_title(axis, xlabel);
        gtk_plot_axis_move_title(axis, axis->title.angle, axis->title.x, axis->title.y - 0.03);
    }

    gtk_plot_set_legends_border(cd -> plot, 2, 3);
    gtk_plot_legends_move(cd -> plot, .60, .10);
/*
    if ((num_sets == 1 && style == CHART_XYZ) || legends == NULL) {
        gtk_plot_legends_set_attributes(cd -> plot, NULL, 0, NULL, NULL);
        gtk_plot_set_legends_border(cd -> plot, GTK_PLOT_BORDER_NONE, 0);
    }
*/
    cd -> legend_drawn = 0;

    gtk_plot_get_yrange(cd -> plot, &ymin, &ymax);

    cd -> dataset = (GtkPlotData **) malloc(sizeof(GtkPlotData *) * num_sets);
    for (i = 0 ; i < num_sets ; i++) {
        cd -> dataset[i] = GTK_PLOT_DATA(gtk_plot_data_new());
        gtk_plot_add_data(cd -> plot, cd -> dataset[i]);
        gtk_widget_show(GTK_WIDGET(cd -> dataset[i]));

        if ((i == 0 && style == CHART_XYZ) || legends == NULL || legends[i] == NULL) {
            gtk_plot_data_hide_legend(cd -> dataset[i]);
        }
        else {
            gtk_plot_data_set_legend(cd -> dataset[i], legends[i]);
            cd -> legend_drawn = 1;
        } 

        if (style == CHART_PEAK && i == 0) {
            gtk_plot_data_show_labels(cd -> dataset[i], TRUE);
            gtk_plot_data_labels_set_attributes(cd -> dataset[i], "Times-Roman", 9, 0, NULL, NULL);
            gtk_plot_data_set_labels(cd -> dataset[i], NULL);
        }
    }

    if (!cd -> legend_drawn) {
        gtk_plot_legends_set_attributes(cd -> plot, NULL, 0, NULL, NULL);
        gtk_plot_set_legends_border(cd -> plot, GTK_PLOT_BORDER_NONE, 0);
    }

    cd -> legend_available = cd -> legend_drawn;

    cd -> menu = makemenu(cd, legends, num_sets);

    // setup the current point marker (has to have its own single point
    // dataset unfortunately)

    gtk_plot_get_yrange(cd -> plot, &ymin, &ymax);

    if (style == CHART_SCATTER || style == CHART_PROFILE) {
        gdk_color_parse("black", &color);
        gdk_color_alloc(gdk_colormap_get_system(), &color);
        cd -> currpt = (GtkPlotData **) malloc(sizeof(GtkPlotData *) * num_sets);
        for (i = 0 ; i < num_sets ; i++) {
            cd -> currpt[i] = GTK_PLOT_DATA(gtk_plot_data_new());
            gtk_plot_add_data(cd -> plot, cd -> currpt[i]);
            gtk_widget_show(GTK_WIDGET(cd -> currpt[i]));
            gtk_plot_data_hide_legend(cd -> currpt[i]);
            gtk_plot_data_set_line_attributes(cd -> currpt[i],
                                           GTK_PLOT_LINE_NONE,
                                           0, 0, 1, &color);
            gtk_plot_data_set_symbol(cd -> currpt[i],
                                  GTK_PLOT_SYMBOL_CROSS,
                                  GTK_PLOT_SYMBOL_EMPTY,
                                  8, 2, &color, &color);
            gtk_plot_data_set_numpoints(cd -> currpt[i], 1);
        }
    }
    else {
        cd -> currpt = NULL;
    }


    gtk_plot_get_yrange(cd -> plot, &ymin, &ymax);

    for (i = 0 ; i < num_sets ; i++) {
        gdk_color_parse(set_colors[i % 6], &color);
        gdk_color_alloc(gdk_colormap_get_system(), &color);

        if (i == 0 && style == CHART_XYZ) {
            gtk_plot_data_set_symbol(cd -> dataset[0],
                             GTK_PLOT_SYMBOL_CIRCLE,
                             GTK_PLOT_SYMBOL_FILLED,
                             4, 1, &color, &color);
        }

        if (style == CHART_SCATTER) {
            gtk_plot_data_set_line_attributes(cd -> dataset[i],
                                              GTK_PLOT_LINE_NONE,
                                              0, 0, 1, &color);
            gtk_plot_data_set_symbol(cd -> dataset[i],
                                     GTK_PLOT_SYMBOL_CIRCLE,
                                     GTK_PLOT_SYMBOL_FILLED,
                                     4, 1, &color, &color);
        }
        else {
            gtk_plot_data_set_line_attributes(cd -> dataset[i],
                                             GTK_PLOT_LINE_SOLID,
                                             0, 0, 1, &color);
        }
    }

    gtk_plot_canvas_paint(GTK_PLOT_CANVAS(cd -> canvas));
    gtk_plot_canvas_refresh(GTK_PLOT_CANVAS(cd -> canvas));

    gtk_plot_clip_data(cd -> plot, TRUE);

    // gtk_plot_canvas_set_size(GTK_PLOT_CANVAS(cd -> canvas), 600, 500);

    gtk_widget_set_size_request(GTK_WIDGET(cd -> canvas), -1, -1);
    gtk_widget_set_size_request(GTK_WIDGET(cd -> plot), -1, -1);
    gtk_widget_set_size_request(GTK_WIDGET(vbox1), -1, -1);
    gtk_window_set_default_size(GTK_WINDOW(cd -> toplevel), 300, 300);
    // gtk_widget_show(cd -> toplevel);

    makelimdlg(cd);


    if (cd -> pixbuf) {
        render_background_image(cd);
    }

    gtk_signal_connect(GTK_OBJECT(cd -> canvas), "size-allocate",
                       GTK_SIGNAL_FUNC(resize), cd);

    gtk_signal_connect(GTK_OBJECT(cd -> canvas), "button_press_event",
                       G_CALLBACK(canvas_click), cd);
    {
        gint w, h;
//        gtk_window_move(GTK_WINDOW(cd -> toplevel), 0, 400);
//        gtk_window_get_size(GTK_WINDOW(cd -> toplevel), &w, &h);
        // printf("before resize %d %d\n", w, h);
 //       gtk_window_resize(GTK_WINDOW(cd -> toplevel), 600, 500);
  //      gtk_window_get_size(GTK_WINDOW(cd -> toplevel), &w, &h);
        // printf("after resize %d %d\n", w, h);
    }
    
    cd -> running = 1;
        
    return (DisplayObject) cd;
}
