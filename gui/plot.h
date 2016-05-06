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

#ifndef _DISPLAY_H
#define _DISPLAY_H

typedef struct {
    void *toplevel;
    int   destroy;
    void *data;
} *DisplayObject;

typedef struct {
    GtkWidget *dialog;
    GtkWidget **field;
} EntriesDialog;

EntriesDialog *CreateEntriesDialog(GtkWidget *top, char *title, char **fields);

typedef void *StatusItem;

#define CHART_XY 1
#define CHART_XYZ 2
#define CHART_SCATTER 3
#define CHART_IMAGE 4
#define CHART_PROFILE 5
#define CHART_PEAK 6
#define CHART_ANIMATE 7

void SetDisplayObjectGeometry(DisplayObject obj, int x, int y, int width, int height);
void GetDisplayObjectGeometry(DisplayObject obj, int *x, int *y, int *width, int *height);

DisplayObject CreateChartDisplay(char *title, int style, int numsets, char *xlabel, char *ylabel, char **legend, char *background_image, double *params, int autoscroll, int num_points, void *data);
void PlotChartDisplay(DisplayObject obj, double *x, double **y);
void UpdateChartDisplay(DisplayObject obj, double x, double *y, double z);
void UpdateChartAnimation(DisplayObject obj, char *);
double *ChartX(DisplayObject obj);
double **ChartY(DisplayObject obj);
void GetChartDisplayLimits(DisplayObject obj, double *limits);
GdkPixbuf *ChartPixbuf(DisplayObject obj, GdkPixbuf *dest);


DisplayObject CreateElapsedClockDisplay(char *title, char *fontname);
DisplayObject CreateUTCClockDisplay(char *title, char *fontname);
void UpdateElapsedClockDisplay(DisplayObject obj, double hours);
void UpdateUTCClockDisplay(DisplayObject obj, double yday);
char *GetClockDisplayFont(DisplayObject obj);

DisplayObject CreateTextDisplay(char *title, char **formats, int nformats, char *fontname);
void UpdateTextDisplay(DisplayObject obj, double x[][6]);
char *GetTextDisplayFont(DisplayObject obj);

#define TRIAXUS_E 1
#define TRIAXUS_D 2
#define TRISOARUS 3

DisplayObject CreateVehicleDisplay(char *title, int type);
void UpdateVehicleDisplay(DisplayObject obj, double p, double r, double y, double pflap, double rflap, double yflap);

DisplayObject CreateStatusDisplay(char *title, GtkWidget *);
void UpdateDisplayItem(StatusItem item, int state);
StatusItem CreateStatusItem(DisplayObject status_display, char *name);

DisplayObject CreateWaterfallDisplay(char *title, char *xlabel, char *ylabel, char *zlabel, char *clabel, double *params);

#endif // _DISPLAY_H
