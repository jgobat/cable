/*
    This file is part of WHOI Cable, a program for the static and dynamic
    analysis of oceanographic cable structures.

    Copyright (C) 1997-2016 by Woods Hole Oceanographic Institution (WHOI)
    and Jason Gobat

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

/****************************************************************************
 *
 * File:        bom.c
 *
 * Description: generate materials bills from WHOI Cable input
 *
 * History:
 *
 ****************************************************************************/

# include <stdio.h>
# include <math.h>
# include <malloc.h>
# include <stdlib.h>
# include <string.h>
# include <ctype.h>
# include "compress.h"
# include "problem.h"
# include "tension.h"
# include "output.h"
# include "error.h"
# include "solve.h"
# include "segments.h"

static int nseg;
static Segment *seg;

void
rtrim(char *x)
{
    char *tail;
    if (x == NULL || x[0] == 0 || strlen(x) < 2)
        return;

    tail = &(x [strlen(x) - 1]);
    while (isspace(*tail) && tail != x)
        tail --;

    *(tail + 1) = 0;
    return;
}

void 
SummarizeWeights (FILE *fp, Problem *p, Segment *seg, int ns, double *w0, double *EA, double *length)
{
   double	weight;
   double	ws;
   double	segw;
   int		i;
   int		j;
   double	k_equiv, k_total; 
   char     buff[24];

   weight = p -> terminal[2] -> buoy -> w - p -> terminal[2] -> buoy -> buoyancy;
   if (length)
       *length = 0.0;
   k_total = 0.0;
  
   fprintf(fp, "\n\nBuoyancy/weight summary - from buoy downwards, weights below named piece\n");
   if (p -> type == Surface) {
       weight = 0;
       fprintf(fp, "WARNING: this summation is not an accurate proxy for tension\n");
   }

   fprintf(fp, "\n%23s: %g lbs\n", p -> terminal[2] -> buoy -> name, weight/4.4482216);
 
   for (i = ns ; i >= 1 ; i--) {
      segw = 0.0;

      ws = seg [i] -> length * seg [i] -> material -> wet;
      segw += ws;

      for (j = 1 ; j <= seg [i] -> num_attach ; j++) {
         segw += seg [i] -> attach [j].object -> wet 
                 * seg [i] -> attach [j].num_nodes;         
      }

      weight += segw;

      sprintf(buff, "%6.1f m %-16s", 
              seg[i] -> length, seg[i] -> material -> name);
      rtrim(buff);   
      fprintf(fp, "%23s: %-g lbs\n", buff, weight/4.4482216);
 
      if (seg [i] -> connector != NULL) {
         weight += seg [i] -> connector -> wet;
         
         fprintf(fp, "%23s: %-g lbs\n", seg[i] -> connector -> name,
                 weight/4.4482216);
      }

/*
      if (weight <= 0.0)
         weight = 0.0;
*/
      k_equiv = seg [i] -> material -> EA / seg [i] -> length;
      k_total += 1.0/k_equiv;

      if (length) {
         *length += seg [i] -> length;
      }
   }

   if (EA) {
       k_total = 1.0 / k_total;
       *EA = k_total * (*length);
   }

   if (w0 && length) {
       *w0 = weight / *length;
   }

   return;
}

static int
AddMaterial(Item item, void *call_data)
{
    FILE     *out_fp = (FILE *) call_data;
    Material mat = (Material) item;
    int      i;
    double   length, mass, wet;
    int      header;

    header = 0;
    length = mass = wet = 0;
    for (i = 1 ; i <= nseg ; i++) {
        if (seg[i] -> material == mat) {
            if (!header) {
                header = 1;
                fprintf(out_fp, "%s lengths: %.1f", mat -> name, seg[i] -> length);
            }                
            else {
                fprintf(out_fp, ", %.1f", seg[i] -> length);
            }

            wet    += seg[i] -> length * mat -> wet;
            mass   += seg[i] -> length * mat -> m;
            length += seg[i] -> length;
        }
    }
    if (header) {
        fprintf(out_fp, "\n   totals: length     %.1f\n", length);
        fprintf(out_fp, "           air mass   %.1f\n", mass);
        fprintf(out_fp, "           wet weight %.1f lbs\n", wet/4.4482216);
    }

    return 0;
}

static int
AddConnector(Item item, void *call_data)
{
    FILE     *out_fp = (FILE *) call_data;
    Connector conn = (Connector) item;
    int i, j;
    int count;

    count = 0;
    for (i = 1 ; i <= nseg ; i++) {
        if (seg[i] -> connector == conn) {
            count ++;
        }
        for (j = 1 ; j <= seg[i] -> num_attach ; j++) {
            if (seg[i] -> attach[j].object == conn) {
                count += seg[i] -> attach[j].num_nodes;
            }
        }
    }

    if (count) {
        fprintf(out_fp, "%s:\n", conn -> name);
        fprintf(out_fp, "   quantity:   %d\n", count);
        fprintf(out_fp, "   wet weight: %.1f\n", conn -> wet * count);
        fprintf(out_fp, "   air mass:   %.1f\n", conn -> m * count);
    }

    return 0;
}

void
LayoutSummary(FILE *fp, Problem *p, Environment *e, Segment *seg, int nseg)
{
    int     i, j, k;
    double  length;
    int     nodecount, *nodes;
    int    *nodelist;
    double *lengthlist;

    nodelist = (int *) malloc(sizeof(int) * nseg); nodelist --;
    lengthlist = (double *) malloc(sizeof(double) * nseg); lengthlist --;
    nodes = (int *) malloc(sizeof(int) * nseg); nodes --;

    length    = 0.0;
    nodecount = 0;

    for (i = 1; i <= nseg ; i++) {
        length += seg[i] -> length;
        nodes[i] = 0;
        for (j = 1 ; j <= seg[i] -> num_dist ; j++) {
            nodes[i] += seg[i] -> dist[j].nodes;
        }
        nodecount += nodes[i];

        nodelist[i] = nodecount;
        lengthlist[i] = length;
    }

    fprintf(fp, "\nTop terminal is %s\n", 
            p -> terminal[2] -> anchor ? p -> terminal[2] -> anchor -> name :
            p -> terminal[2] -> buoy -> name);
            
    fprintf(fp, "\n     Segment               length  nodes  DEPTH   ALTIT   NODES\n");
    fprintf(fp, "--------------------------------------------------------------\n");
    for (i = nseg, k = 1 ; i >= 1 ; i--, k++) {
        if (i < nseg && seg[i] -> connector)  {
            fprintf(fp, "     +%s\n", seg[i] -> connector -> name);   
        }
        else if (i < nseg && seg[i] -> connection == Pinned) {
            fprintf(fp, "     +pinned\n");
        }
        fprintf(fp, "%04d %-20s %7.1f %5d %7.1f %7.1f %6d\n", nseg - i + 1,
                seg[i] -> material -> name, seg[i] -> length,
                nodes[i], e -> depth - lengthlist[i], lengthlist[i], nodelist[i]); 
    }

    fprintf(fp, "\nBottom terminal is %s\n\n", 
            p -> terminal[1] -> anchor ? p -> terminal[1] -> anchor -> name :
            p -> terminal[1] -> buoy -> name);
            

    nodelist ++;
    lengthlist ++;

    free(nodelist);
    free(lengthlist);
}

int 
MaterialsBill(FILE *fp, Problem *p, Environment *e) 
{
    if (!fp) {
        fp = stdout;
    }

    // BuildSegmentArray(p, &nseg, 1); // include branches
    seg = p -> segment;
    nseg = p -> num_segments;

    LayoutSummary(fp, p, e, seg, nseg);

    TreeSetAndIterate(p -> material_tree, AddMaterial, fp);
    TreeSetAndIterate(p -> connector_tree, AddConnector, fp);

    if (p -> type == Subsurface || p -> type == Surface) {
        SummarizeWeights(fp, p, seg, nseg, NULL, NULL, NULL);
    }

    return 0;
}

