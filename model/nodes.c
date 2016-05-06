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
 * File:	nodes.c
 *
 * Description: routines to create the array of nodes from the segment tree
 *		
 ****************************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "compress.h"
# include "problem.h"
# include "allocate.h"
# include "error.h"

extern Analysis *analysis;

static Node	*n;	/* array of node pointers	*/
static int	 nn;	/* total number of nodes	*/
static int	 cn;	/* current node number		*/
static Node 	 prev_main_node;

static int 
CountNodes (Item item, void *data)
{
   Segment	s;
   int		i;

   s = (Segment) item;

   for (i = 1 ; i <= s -> num_dist ; i++)
      nn += s -> dist [i].nodes;

   nn -= (s -> num_dist - 1);

   if (s -> num_top_dist) {
      for (i = 1 ; i <= s -> num_top_dist ; i++)
         nn += s -> top_dist [i].nodes;

      nn -= s -> num_top_dist;
   }

   if (s -> num_bottom_dist) {
      for (i = 1 ; i <= s -> num_bottom_dist ; i++)
         nn += s -> bottom_dist [i].nodes;

      nn -= s -> num_bottom_dist;
   }

   return 0;
}

static void 
HookupAttachments(Segment s, int start_node)
{
   int	i, j;

   for (i = 1 ; i <= s -> num_attach ; i++) {
      for (j = 1 ; j <= s -> attach [i].num_nodes ; j++) {
         if (s -> attach [i].nodes [j] > s -> num_nodes ||
             s -> attach [i].nodes [j] < 1)
            ExitErr ("segment %d has no node %d for an attachment", 
                   s -> number, s -> attach [i].nodes [j]);

         n [start_node - 1 + s -> attach [i].nodes [j]] -> attachment = 
            s -> attach [i].object;
      }
   }

   return;
}

static void 
ProcessSegment(Segment s)
{
   int	  i, j;
   double ds;
   int	  count;
   int    start_node;
   double p;
   int	  first_active, last_active;

   ds = 0.0;

   start_node = cn + 1;
   s -> num_nodes = 0;
   s -> first = n [start_node];

   if (s -> branch) {
      if (s == s -> branch -> segment [1])
         n [start_node] -> s = n [cn] -> s; // no dimensionality to of  
   					                        // connectors in the branch 
				                            // direction, so ignore ds  
				 	                        // of the mainline node     
      else
         n [start_node] -> s = n [cn] -> s + n [cn] -> ds;
   }
   else if (prev_main_node)
      n [start_node] -> s = prev_main_node -> s
                            + prev_main_node -> ds;
  
	// add nodes spooled at bottom of segment

   if (s -> num_bottom_dist) {

      s -> bottom_wet = s -> bottom_length * s -> material -> wet;
      s -> bottom_spooled = s -> bottom_length;

      p = 0.0;
      for (i = 1 ; i <= s -> num_bottom_dist ; i++) {

         s -> num_nodes += s -> bottom_dist [i].nodes;

         if (i < s -> num_bottom_dist) 
            count = s -> bottom_dist [i].nodes - 1.0;
         else
            count = s -> bottom_dist [i].nodes;

         count = s -> bottom_dist [i].nodes - 1.0;

         ds = s -> bottom_dist [i].percent * s -> bottom_length / (s -> bottom_dist [i].nodes - 1.0);

         for (j = 1 ; j <= count ; j++) {
            cn ++;
            n [cn] -> number       = cn;
            n [cn] -> segment    = s;
            n [cn] -> material   = s -> material;
            n [cn] -> attachment = NULL;
            n [cn] -> ds = 0;
            n [cn] -> ds0         = ds;
            n [cn] -> active     = 0;
            n [cn] -> position   = Cable;
            if (cn > start_node)
               n [cn] -> s          = n [cn - 1] -> s + ds;
         }

         p += s -> bottom_dist [i].percent;
      }

      if (p < 1.0 - 1.0e-6 || p > 1.0 + 1.0e-6)
         ExitErr ("bottom spool nodes for segment %d only cover %g of spooled length", s -> number, p);

      // n [cn] -> ds = 0;			
      n [cn + 1] -> s = n [cn] -> s + n [cn] -> ds;	
   }

	// add nodes on active segment

   first_active = cn + 1;

   p = 0.0;
   for (i = 1 ; i <= s -> num_dist ; i++) {

      s -> num_nodes += s -> dist [i].nodes;

      if (i < s -> num_dist) 
         count = s -> dist [i].nodes - 1.0;
      else
         count = s -> dist [i].nodes;

      ds = s -> dist [i].percent * s -> length / (s -> dist [i].nodes - 1.0);

      for (j = 1 ; j <= count ; j++) {
         cn ++;
         n [cn] -> number       = cn;
         n [cn] -> segment    = s;
         n [cn] -> material   = s -> material;
         n [cn] -> attachment = NULL;
         n [cn] -> ds = n [cn] -> ds0         = ds;
         n [cn] -> active     = 1;
         n [cn] -> position   = Cable;
         if (cn > first_active)
            n [cn] -> s      = n [cn - 1] -> s + ds;
      }

      p += s -> dist [i].percent;
   }

   if (p < 1.0 - 1.0e-6 || p > 1.0 + 1.0e-6)
      ExitErr ("nodes for segment %d only cover %g of length", s -> number, p);

   last_active = cn;

	// add nodes on spool at top of segment

   if (s -> num_top_dist) {

      s -> top_wet = s -> top_length * s -> material -> wet;
      s -> top_spooled = s -> top_length;

      p = 0.0;
      for (i = 1 ; i <= s -> num_top_dist ; i++) {

         s -> num_nodes += s -> top_dist [i].nodes;
         if (i < s -> num_top_dist) 
            count = s -> top_dist [i].nodes - 1.0;
         else
            count = s -> top_dist [i].nodes;

         count = s -> top_dist [i].nodes - 1.0;

         ds = s -> top_dist [i].percent * s -> top_length / (s -> top_dist [i].nodes - 1.0);

         for (j = 1 ; j <= count ; j++) {
            cn ++;
            n [cn] -> number     = cn;
            n [cn] -> segment    = s;
            n [cn] -> material   = s -> material;
            n [cn] -> attachment = NULL;
            n [cn] -> ds = 0;
            n [cn] -> ds0         = ds;
            n [cn] -> active     = 0;
            n [cn] -> position   = Cable;
            if (cn > last_active + 1)
               n [cn] -> s       = n [cn - 1] -> s + ds;
            else 
               n [cn] -> s       = n [last_active] -> s;

            n [cn] -> s       = n [cn - 1] -> s + ds;
         }

         p += s -> top_dist [i].percent;
      }

      if (p < 1.0 - 1.0e-6 || p > 1.0 + 1.0e-6)
         ExitErr ("top spool nodes for segment %d only cover %g of spooled length", s -> number, p);

      n [cn] -> ds = 0;			
      n [last_active] -> ds0 = ds;
   }

   if (s -> connector && s -> num_top_dist == 0)
      n [last_active] -> ds       = s -> connector -> length;
   else
      n [last_active] -> ds       = 0.0;

   s -> last = n [cn];
   s -> last_active = n [last_active];
   s -> first_active = n [first_active];

	/*
	 * Decide how to flag the first and last nodes on this segment.
	 * We can't rely on checking branch -> first,last for branched 
	 * segments because those don't get set until we build the branch 
	 * array.  We also need to determine if this branch 
	 * branches _to_ a connector (i.e., the branch forms a
 	 * loop).  If it does then we need to turn a mainline node 
	 * into a junction.  We do these steps in ProcessJunctions
	 */


   if (s -> branch) {
      if (s == s -> branch -> segment [s -> branch -> num_segment]) 
         n [last_active] -> position = BranchTerminal;
      else 
         n [last_active] -> position = Cable;

      if (s == s -> branch -> segment [1]) 
         n [first_active] -> position = BranchStart;
      else 
         n [first_active] -> position = Connection;
   }
   else {
      if (cn == nn)
         n [last_active] -> position = TopBoundary;
      else 
         n [last_active] -> position = Cable;

      if (s -> number == 1)
         n [first_active] -> position = BottomBoundary;
      else 
         n [first_active] -> position = Connection;
   }

   HookupAttachments(s, start_node);

   return;
}

static int 
ProcessMainSegment (Item item, void *data)
{
   Segment	      s = (Segment) item;
   int		      i, j;

	/*
	 * save the node number of the first node on this segment	
	 */

   ProcessSegment(s);

   prev_main_node = s -> last;

   for (i = 1 ; i <= s -> num_branch_to ; i++) {
      for (j = 1 ; j <= s -> branch_to [i] -> num_segment ; j++) 
         ProcessSegment (s -> branch_to [i] -> segment [j]);

      s -> branch_to [i] -> first = s -> branch_to [i] -> segment [1] -> first;
      s -> branch_to [i] -> last = s -> branch_to [i] -> segment [s -> branch_to [i] -> num_segment] -> last;
   }

   return 0;
}

static int 
ProcessJunctions (Item item, void *data)
{
   Problem *problem  = (Problem *) data;
   int	  i, idx;
   Node   nd;
   Branch b;

   b = (Branch) item;

	/*
	 * this is code to resolve looped junctions
	 */

   if (b -> terminal -> loop_main_node) { 
      idx = (int) b -> terminal -> loop_main_node;
 
      if (idx < 1 || idx > nn || n [idx] -> segment -> branch)
         ExitErr("looped branches must connect into nodes on main line");

      if (n [idx] != n [idx] -> segment -> last || !n [idx] -> segment -> connector)
         ExitErr("looped branches must join into a connector node");

      b -> terminal -> loop_main_node = n [idx];

      n [idx] -> next -> position = Junction;
      n [idx] -> next -> segment -> junction.num_nodes ++;

      i = n [idx] -> next -> segment -> junction.num_nodes;
      n [idx] -> next -> segment -> junction.node [i] = 
                                  b -> segment [b -> num_segment] -> last;

      if (i > problem -> junction_size)
         problem -> junction_size = i;
   }

	/*
	 * this is code to resolve regular junctions
	 */

   nd = b -> segment_from -> last -> next;

   nd -> position = Junction;
   nd -> segment -> junction.num_nodes ++;

   i = nd -> segment -> junction.num_nodes;
   nd -> segment -> junction.node [i] =  b -> segment [1] -> first;

   if (i > problem -> junction_size)
      problem -> junction_size = i;

   return 0;
}

void 
MakeNextPrevOutput(Problem *problem, Node *node, int num_nodes)
{
   int		i, j;
   Node		prev;
   Node		prev_main;
   int		main_number;
   int		branch_number;
   Node		a;
 
   prev = prev_main = NULL;

   main_number = 0;  
   branch_number = 0;
  
   for (i = 1 ; i <= num_nodes ; i++) {
      node[i] -> next = NULL;
      node[i] -> prev = NULL;
      
      if (!node[i] -> segment -> branch) {	 // mainline
	     node[i] -> prev = prev_main;
         if (prev_main)
	        prev_main -> next = node[i];

	     prev = prev_main = node[i];
         
         main_number ++;
         node[i] -> output_number = main_number;
      }
      else {
         if (prev -> segment -> branch == node[i] -> segment -> branch) {
            node[i] -> prev= prev;
            prev -> next = node[i];
         }
	     else
	        node[i] -> prev = NULL;

	     prev = node[i];

         branch_number ++;
         node[i] -> output_number = branch_number;
      }           
   }

   for (j = 1; j <= problem -> num_branch ; j++) {
      a = problem -> branch[j] -> first;

      while (a) {
         a -> output_number += main_number;
         a = a -> next;
      }
   }

   for (i = 1 ; i <= num_nodes ; i++)    {
      // fprintf(stderr,"%d < %d > %d\n", node[i] -> prev ? node[i] -> prev -> number : 0, node[i] -> number, node[i] -> next ? node[i] -> next -> number : 0);
   }

   return;
}

int 
MakeNextPrevActive(Problem *problem, Node *node, int num_nodes, Node* active)
{
   int		n, i;
   Node		prev;
   Node		prev_main;
  
   prev = prev_main = NULL;
  
   n = 0;
   for (i = 1 ; i <= num_nodes ; i++)
      if (node [i] -> active) {
         active [++ n] = node [i];
         active [n] -> active_number = n;
         active [n] -> next_active = NULL;
         active [n] -> prev_active = NULL;

         if (!active[n] -> segment -> branch) {	 // mainline
            active[n] -> prev_active = prev_main;
            if (prev_main)
               prev_main -> next_active = active[n];

            prev = prev_main = active[n];
         }
         else {
            if (prev -> segment -> branch == active[n] -> segment -> branch) {
               active [n] -> prev_active = prev;
               prev -> next_active = active [n];

               // this is a bit inefficient but will eventually be right
               active[n] -> segment -> branch -> terminal -> node = active[n];
            }
            else
               active [n] -> prev_active = NULL;

            prev = active[n];
         }           
      }
      else
         node [i] -> active_number = 0;

   problem -> terminal[1] -> node = active[1];
   problem -> terminal[2] -> node = active[n];

#if 0
   fprintf(stderr,"num_active = %d\n", n);

   for (i = 1 ; i <= num_nodes ; i++)    {
     fprintf(stderr,"%d active=%d, active_num=%d, type=%d, nextAct=%d, prevAct=%d\n",
             i, node[i] -> active, node[i] -> active_number, 
             node[i] -> position, node[i] -> next_active ? node[i] -> next_active -> number : -1, node[i] -> prev_active ? node[i] -> prev_active -> number : -1);
 

      // fprintf(stderr,"%d < %d > %d\n", active[i] -> prev_active ? active[i] -> prev_active -> number : 0, active[i] -> number, active[i] -> next_active ? active[i] -> next_active -> number : 0);
   }

   printf("ACTIVE = %d\n", n);
   for (i = 1 ; i <= n ; i++)
     fprintf(stderr,"%d %d active=%d, active_num=%d, type=%d, nextAct=%d, prevAct=%d, ds=%g, ds0=%g\n",
             i, active[i] -> number,  active[i] -> active, active[i] -> active_number, 
             active[i] -> position, active[i] -> next_active ? active[i] -> next_active -> number : -1, active[i] -> prev_active ? active[i] -> prev_active -> number : -1, active[i] -> ds, active[i] -> ds0);
#endif

   return n;
}


void 
CreateNodeArray (Problem *p, Node **node, int *number_of_nodes, 
                 Node **active, int *number_of_active_nodes)
{
   int		 i, j;
   Node		*a;

	/*
	 * figure out how many nodes there are in the problem
	 */

   nn = 0;
   TreeSetIterator (p -> segment_tree, CountNodes, NULL);
   TreeIterate (p -> segment_tree);

   if (nn <= 0) {
      error("nothing to do - no nodes defined in problem");
      return;
   }
   
   TreeSetIterator (p -> branch_segment_tree, CountNodes, NULL);
   TreeIterate (p -> branch_segment_tree);

	/*
	 * create the array of node pointers
	 */

   n = Allocate (Node, nn);
   if (n == NULL)
      ExitErr ("could not allocate memory for node array");

   UnitOffset (n);

   a = Allocate(Node, nn);
   if (a == NULL)
      ExitErr ("could not allocate memory for active node array");

   UnitOffset(a);
 
	/*
	 * create the actual node structures at each point in the array
	 */

   for (i = 1 ; i <= nn ; i++) {
      n [i] = AllocNew (struct node);

      n [i] -> next = NULL;
      n [i] -> prev = NULL;
      n [i] -> next_active = NULL;
      n [i] -> prev_active = NULL;
      n [i] -> number = n [i] -> active_number = n [i] -> output_number = 0;
      n [i] -> segment = NULL;
      n [i] -> material = NULL;
      n [i] -> position = Cable;
      n [i] -> attachment = NULL;
      n [i] -> s = 0.0;
      n [i] -> ds = 0.0;
      n [i] -> ds0 = 0.0;
      n [i] -> active = 0;
   
      n[i] -> Ys = (double *) malloc(sizeof(double) * 10); n[i] -> Ys --;
      n[i] -> Y = (double *) malloc(sizeof(double) * 13); n[i] -> Y --;
      n[i] -> Y_o = (double *) malloc(sizeof(double) * 13); n[i] -> Y_o --;      
      n[i] -> Y_f = (double *) malloc(sizeof(double) * 13); n[i] -> Y_f --;
      n[i] -> Y_o_f = (double *) malloc(sizeof(double) * 13); n[i] -> Y_o_f --;      
      n[i] -> Yd_o = (double *) malloc(sizeof(double) * 13); n[i] -> Yd_o --;
      n[i] -> Yd = (double *) malloc(sizeof(double) * 13); n[i] -> Yd --;
   
      for (j = 1 ; j <= 10 ; j++)
         n [i] -> Ys [j] = 0.0;

      for (j = 1 ; j <= 13 ; j++) {
         n [i] -> Y [j] = 0.0;
         n [i] -> Yd [j] = 0.0;
         n [i] -> Y_o [j] = 0.0;
         n [i] -> Yd_o [j] = 0.0;
         n [i] -> Y_f [j] = 0.0;
         n [i] -> Y_o_f [j] = 0.0;
      }

      n [i] -> x = 0.0;
      n [i] -> y = 0.0;
      n [i] -> z = 0.0;
      n [i] -> x_o = 0.0;
      n [i] -> y_o = 0.0;
      n [i] -> z_o = 0.0;
      n [i] -> xdot_o = n[i] -> xdot = 0.0;
      n [i] -> ydot_o = n[i] -> ydot = 0.0;
      n [i] -> zdot_o = n[i] -> zdot = 0.0;
 
      n [i] -> pay = 0.0;
      n [i] -> pay_o = 0.0;

      n [i] -> x_f = 0.0;
      n [i] -> y_f = 0.0;
      n [i] -> z_f = 0.0;
      n [i] -> x_o_f = 0.0;
      n [i] -> y_o_f = 0.0;
      n [i] -> z_o_f = 0.0;
      n [i] -> xdot_o_f = n[i] -> xdot_f = 0.0;
      n [i] -> ydot_o_f = n[i] -> ydot_f = 0.0;
      n [i] -> zdot_o_f = n[i] -> zdot_f = 0.0;

      n [i] -> pay_f = 0.0;
      n [i] -> pay_o_f = 0.0;

      n [i] -> lift = 0.0;
      n [i] -> drag = 0.0;
   }
 
	/*
	 * walk the segment tree once more to actually fill the
	 * node structures out
	 */

   cn = 0;

   n [1] -> s = 0.0;
   prev_main_node = NULL;
   TreeSetIterator (p -> segment_tree, ProcessMainSegment, NULL);
   TreeIterate (p -> segment_tree);

   if (number_of_nodes)
      *number_of_nodes = nn;

   *node = n;

   *number_of_active_nodes = MakeNextPrevActive(p, n, nn, a);
   *active = a;

   MakeNextPrevOutput(p, n, nn);

   TreeSetIterator (p -> branch_tree, ProcessJunctions, p);
   TreeIterate (p -> branch_tree);


   for (i = 1 ; i <= nn ; i++) {
      // n [i] -> ds0 = n [i] -> ds;
      n [i] -> pay = 0.0;
      n [i] -> pay_o = 0.0;
 
//      fprintf(stderr,"n = %d, l = %d, a = %d, active = %d, s = %g, ds = %g, pos = %d\n", n [i] -> number, n [i] -> local_number, n [i] -> active_number, n [i] -> active, n [i] -> s, n [i] -> ds, n [i] -> position);        
   }
 
   return;
}


