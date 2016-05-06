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

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gtk/gtk.h>

#include <ctype.h>

#ifdef WINDOWS
char *
strcasestr (char *haystack, char *needle)
{
    char *p, *startn = 0, *np = 0;

    for (p = haystack; *p; p++) {
        if (np) {
            if (toupper(*p) == toupper(*np)) {
                if (!*++np)
                    return startn;
            } else
                np = 0;
        } else if (toupper(*p) == toupper(*needle)) {
            np = needle + 1;
            startn = p;
        }
    }

    return 0;
}
#endif

char *
AddFilenameExtension(char *filename, char *ext)
{
    char *newname;

    if (strcasestr(filename, ext)) {
        return filename;
    };

    newname = (char *) malloc(sizeof(char) * (strlen(filename) + strlen(ext) + 1));
    strcpy(newname, filename);
    strcat(newname, ext);
    free(filename);
    return newname; 
}

char *
FilenameFromDialog(char *title, GtkWidget *parent, 
                   GtkFileChooserAction action, char **prev_file, 
                   char *ext, char *descr,
                   GtkWidget *extra, gboolean *state)
{
    char    cwd[256];
    GtkWidget *dialog;
    char      *filename;
    static char *prev_prev_file = NULL;
    GtkFileFilter *filt;

    dialog = gtk_file_chooser_dialog_new (title,
                  GTK_WINDOW(parent),
                  action,
                  GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                  GTK_STOCK_OK, GTK_RESPONSE_ACCEPT,
                  NULL);

    if (ext && descr) {
        filt = gtk_file_filter_new();
        gtk_file_filter_set_name(filt, "all files");
        gtk_file_filter_add_pattern(filt, "*");
        gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(dialog), filt);

        filt = gtk_file_filter_new();
        gtk_file_filter_set_name(filt, descr);
        gtk_file_filter_add_pattern(filt, ext);
        gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(dialog), filt);
        gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(dialog), filt);
    }

    if (extra) {
        gtk_file_chooser_set_extra_widget(GTK_FILE_CHOOSER(dialog), extra);
        gtk_widget_show(extra);
    }

    if (prev_file && *prev_file) {
        fprintf(stderr,"set 1 %s\n", *prev_file);
        gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), *prev_file);  
    }
    else if (prev_prev_file) {
        fprintf(stderr,"set 2\n");
        gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), prev_prev_file);  
    }
    else {
        getcwd(cwd, 256);
        fprintf(stderr,"setting current to cwd %s\n", cwd);

        gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), cwd);
    }

    filename = NULL;
    if (gtk_dialog_run (GTK_DIALOG (dialog)) == GTK_RESPONSE_ACCEPT) {
        filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER (dialog));
    }

    if (prev_file) {
        if (*prev_file) {
            g_free(*prev_file);
        }
        *prev_file = gtk_file_chooser_get_current_folder(GTK_FILE_CHOOSER(dialog));
        prev_prev_file = *prev_file;
    }

    if (extra && state) {
        *state = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(extra));
    }

    gtk_widget_destroy(dialog);
    return filename;
}

int
YesNoDialog(char *msg, GtkWidget *parent)
{
    GtkWidget   *dialog;
    int          resp;

    dialog = gtk_message_dialog_new(GTK_WINDOW(parent),
                                        GTK_DIALOG_MODAL,
                                        GTK_MESSAGE_QUESTION,
                                        GTK_BUTTONS_YES_NO,
                                        msg);
    resp = gtk_dialog_run(GTK_DIALOG (dialog));
    gtk_widget_destroy(dialog);
    return resp;
}

