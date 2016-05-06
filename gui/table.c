/*
 * Routines to insert and delete rows in a GtkTable.
 * The table is assumed to have no children that span more than one
 * row: odd displays may result if this is not true.
 *
 * Rows are numbered by the value of the top_attach field;
 * a table with 2 rows has rows numbered 0 and 1, etc.
 */

#include <gtk/gtk.h>

static GtkTableChild *
collect_children(GtkTable *table, int row, int *nchp)
{
    GtkTableChild *row_children;
    GtkTableChild *table_child;
    int nch;
    int i;
    GList *list;

    row_children = g_new0(GtkTableChild, table->ncols);
    nch = 0;

    /* locate all of the relevant children, and add an extra reference
       to them */
    for (list = table->children; list; list = list->next) {
        table_child = list->data;

        if(table_child->top_attach == row) {
            if(nch == table->ncols) {
                g_error("collect_children: too many children on row %d\n", row);
                return NULL; // abort();
            }
            row_children[nch++] = *table_child;
            gtk_widget_ref(table_child->widget);
        }
    }
    if(nchp)
        *nchp = nch;
    return row_children;
}

static void
reattach_children(GtkTable *table, GtkTableChild *row_children, int nch,
          int row) 
{
    int i;
    GtkAttachOptions xoptions;
    GtkAttachOptions yoptions;

    /* add in new position and unref */
    for(i = 0; i < nch; i++) {
        xoptions = row_children[i].xexpand 
            | row_children[i].xshrink<<1
            | row_children[i].xfill<<2;
        yoptions = row_children[i].yexpand 
            | row_children[i].yshrink<<1
            | row_children[i].yfill<<2;

        gtk_table_attach(table, row_children[i].widget,
                 row_children[i].left_attach,
                 row_children[i].right_attach,
                 row,
                 row + (row_children[i].bottom_attach -
                       row_children[i].top_attach),
                 xoptions,
                 yoptions,
                 row_children[i].xpadding,
                 row_children[i].ypadding);

        gtk_widget_unref(row_children[i].widget);
    }
}



/*
 * move a table row to another row position.
 */
void
gtk_table_move_row(GtkTable *table, int from, int to)
{
    GtkTableChild *row_children;
    int nch;
    int i;
    GtkTableChild *table_child;
    GtkAttachOptions xoptions;
    GtkAttachOptions yoptions;
    GList *list;

/*  printf("move row %d to %d: ", from, to); */
    row_children = g_new(GtkTableChild, table->ncols);
    nch = 0;

    /* locate all of the relevant children, and add an extra reference
       to them */
    for (list = table->children; list; list = list->next) {
        table_child = list->data;

        if(table_child->top_attach == from) {
            if(nch == table->ncols) {
                g_error("gtk_table_move_row: too many children on row %d\n", from);
                return; // abort();
            }
            row_children[nch++] = *table_child;
            gtk_widget_ref(table_child->widget);
        }
    }

    /* remove from old position */
    for(i = 0; i < nch; i++) {
        gtk_container_remove(GTK_CONTAINER(table), row_children[i].widget);
    }
    /* add in new position and unref */
    for(i = 0; i < nch; i++) {
        xoptions = row_children[i].xexpand 
            | row_children[i].xshrink<<1
            | row_children[i].xfill<<2;
        yoptions = row_children[i].yexpand 
            | row_children[i].yshrink<<1
            | row_children[i].yfill<<2;

        gtk_table_attach(table, row_children[i].widget,
                 row_children[i].left_attach,
                 row_children[i].right_attach,
                 to,
                 to + (row_children[i].bottom_attach -
                       row_children[i].top_attach),
                 xoptions,
                 yoptions,
                 row_children[i].xpadding,
                 row_children[i].ypadding);

        gtk_widget_unref(row_children[i].widget);
    }
/*  printf("%d widgets moved\n", nch); */
    
    g_free(row_children);
}

/*
 * Insert a new empty row into a GtkTable just before the specified row,
 * moving subsequent rows down.
 *  - resize the table one larger
 *  - move rows down to put the hole where it belongs
 */
void
gtk_table_insert_row(GtkTable *table, int row)
{
    int old_nrows, old_ncols;
    int i;

    old_nrows = table->nrows;
    old_ncols = table->ncols;

    gtk_table_resize(table, old_nrows+1, old_ncols);
    
    for(i = old_nrows; i > row; i--) {
        gtk_table_move_row(table, i-1, i);
    }
    
}

/*
 * Delete the indicated row of a GtkTable, moving subsequent rows up.
 *  - remove the specified row
 *  - move other rows
 *  - resize the table one smaller
 */
void
gtk_table_delete_row(GtkTable *table, int row)
{
    GtkTableChild *table_child;
    GList *list;
    GList *del_list = NULL;
    int old_nrows, old_ncols;
    int i;

    old_nrows = table->nrows;
    old_ncols = table->ncols;

    /* remove widgets from the old row 
     can't delete while walking the list; must make a temporary list
    */
    for (list = table->children; list; list = list->next) {
        table_child = list->data;
        if(table_child->top_attach == row) {
            del_list = g_list_append(del_list, table_child->widget);
        }
    }

    while(del_list) {
        GtkWidget *child;
        child = del_list->data;
        /* printf("remove child 0x%lx\n", child); */
        gtk_container_remove(GTK_CONTAINER(table), child);
        del_list = g_list_remove(del_list, child);
    }

    for(i = row; i < old_nrows-1; i++) {
        gtk_table_move_row(table, i+1, i);
    }

    gtk_table_resize(table, old_nrows-1, old_ncols);
}

/*
 * Reorder a table, moving row FROM to position TO, and shifting
 * all rows in between up or down one row.
 * given this arrangement:
0 A
1 B
2 C
3 D
rotate_rows(, 2, 0) produces this result:  ("up" case)
0 C
1 A
2 B
3 D
while applying rotate_rows(, 1, 3) to the original produces: ("down" case)
0 A
1 C
2 D
3 B
 */
void
gtk_table_rotate_rows(GtkTable *table, int from, int to)
{
    int down;  /* does the named "jumping" row move up or down? */
    GtkTableChild *moving_children;
    int nch, i;

    if(to == from) return;
    else if(to > from) down = 1;
    else down = 0;

    /* collect the children in the row that's jumping */
    moving_children = collect_children(table, from, &nch);
    /* remove them from their old row */
    for(i = 0; i < nch; i++) {
        gtk_container_remove(GTK_CONTAINER(table), moving_children[i].widget);
    }

    /* now move the other ones one position.
    * if the "jumping" row moves down, thers move up, etc. */
    if(down) {
        for(i = from; i < to; i++) {
            gtk_table_move_row(table, i+1, i);
        }
    } else { /* up*/
        for(i = from; i > to; i--) {
            gtk_table_move_row(table, i-1, i);
        }
    }

    /* finally, put those children back down in the vacated target row */
    reattach_children(table, moving_children, nch, to);
    g_free(moving_children);
    
}

/*
 * given a pointer to a child widget of a table, 
 * return the (top_attach) row number it is currently on.
 */
int
gtk_table_get_child_row(GtkTable *table, GtkWidget *child)
{
    GtkTableChild *table_child;
    GList *list;
    int i;

    for (list = table->children; list; list = list->next) {
        table_child = list->data;
        
        if(table_child->widget == child)
            return table_child->top_attach;
    }
    return -1;
}

GtkWidget *
gtk_table_get_child_at(GtkTable *table, int row, int col)
{
    GtkTableChild *table_child;
    GList *list;
    int i;

    for (list = table->children; list; list = list->next) {
        table_child = list->data;
        
        if(table_child->top_attach == row && table_child->left_attach == col)
            return table_child->widget;
    }
    return NULL;
}
