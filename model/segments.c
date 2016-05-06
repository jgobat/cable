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
 * File:        segments.c
 *
 * Description: routines for building arrays from trees
 *
 * History:
 *
 ****************************************************************************/

# include <stdio.h>
# include <math.h>
# include <malloc.h>
# include <stdlib.h>
# include "compress.h"
# include "problem.h"
# include "Tree.h"
# include "tension.h"


static int       count;
static Segment	*array = NULL;
static Branch	*br_array = NULL;

static int 
AddSegment(Item item, void *data)
{
   Segment	s = (Segment) item;

   array [count ++] = s; 
   
   return 0;
}

static int 
AddBranch(Item item, void *data)
{
   Branch	b = (Branch) item;

    
   // b -> first = b -> segment [1] -> first;
   // b -> last  = b -> segment [b -> num_segment] -> last;

   br_array [count ++] = b; 
   
   return 0;
}

Segment *
BuildSegmentArray (Problem *p, int *ns, int branches)
{
   *ns = TreeSize (p -> segment_tree);
   if (branches)
      *ns += TreeSize(p -> branch_segment_tree);

   if (*ns == 0)
      return NULL;

/* this is a bad idea - calling function should be responsible

   if (array) {
      array ++;
      free (array);
   }
*/

   array = (Segment *) malloc(sizeof(Segment) * *ns);
   array --;

   count = 1;

   TreeSetIterator (p -> segment_tree, AddSegment, NULL);
   TreeIterate (p -> segment_tree);

   if (branches) {
      TreeSetIterator (p -> branch_segment_tree, AddSegment, NULL);
      TreeIterate (p -> branch_segment_tree);
   }

   return array; 
}

Branch *
BuildBranchArray (Problem *p, int *nb)
{

   *nb = TreeSize (p -> branch_tree);
   if (*nb == 0)
      return NULL;

   if (br_array) {
      br_array ++;
      free (br_array);
   }

   br_array = (Branch *) malloc(sizeof(Branch) * *nb);
   br_array --;

   count = 1;

   TreeSetIterator (p -> branch_tree, AddBranch, NULL);
   TreeIterate (p -> branch_tree);

   return br_array; 
}

void CopyNodeData(Node dest, Node src, int ne)
{
   int	k;

   for (k = 1 ; k <= ne ; k++) {
      dest  -> Y[k] = src -> Y[k];
      dest  -> Yd[k] = src -> Yd[k];
      dest  -> Y_o[k] = src -> Y_o[k];
      dest  -> Y_f[k] = src -> Y_f[k];
      dest  -> Y_o_f[k] = src -> Y_o_f[k];
      dest  -> Yd_o[k] = src -> Yd_o[k];
   }
   dest  -> x = src -> x;
   dest  -> y = src -> y;
   dest  -> z = src -> z;

   dest  -> x_o = src -> x_o;
   dest  -> y_o = src -> y_o;
   dest  -> z_o = src -> z_o;

   dest  -> x_f = src -> x_f;
   dest  -> y_f = src -> y_f;
   dest  -> z_f = src -> z_f;

   dest  -> x_o_f = src -> x_o_f;
   dest  -> y_o_f = src -> y_o_f;
   dest  -> z_o_f = src -> z_o_f;

   dest  -> xdot_o = src -> xdot_o;
   dest  -> ydot_o = src -> ydot_o;
   dest  -> zdot_o = src -> zdot_o;

   dest  -> xdot = src -> xdot;
   dest  -> ydot = src -> ydot;
   dest  -> zdot = src -> zdot;
 
   dest  -> xdot_o_f = src -> xdot_o_f;
   dest  -> ydot_o_f = src -> ydot_o_f;
   dest  -> zdot_o_f = src -> zdot_o_f;

   dest  -> xdot_f = src -> xdot_f;
   dest  -> ydot_f = src -> ydot_f;
   dest  -> zdot_f = src -> zdot_f;
 
   dest -> pay_o = src -> pay_o;
   dest -> pay = src -> pay;

   dest -> pay_o_f = src -> pay_o_f;
   dest -> pay_f = src -> pay_f;

   return;
}

int 
ProcessSpools(Problem *p, Node *node, int num_nodes, 
              Node *active, int num_active,
              double t, double dt, int ne, int twoD)
{
   Segment	    *seg;
   int		     nseg;
   int			     i, j;
   double		     T;
   double		     added, speed;
   double	         prev_speed, prev_t, prev_dt;
   int			     first_spooled, last_spooled;
   int			     first_unspooled, last_unspooled;
   EquationType		 pos;
   int			     active_changed;
   double            sf;

   seg  = p -> segment;
   nseg = p -> num_segments;
 
   active_changed = 0;

   for (i = 1 ; i <= nseg ; i++) {
      prev_t = seg[i] -> prev_t;
      prev_dt = seg[i] -> prev_dt;

      if (seg [i] -> num_bottom_dist) {
         printf("processing bottom spool segment %d\n", i);
                 
    	 added = 0.0;
         if (t < prev_t) {
            if (seg [i] -> bottom_pay.expr) {
               T = Tension(seg[i] -> first_active -> Y_o[1], seg [i] -> material);
               prev_speed = EvalCode(seg [i] -> bottom_pay.expr, 
                                     seg[i] -> first_active, prev_t, T,
                                     0, 0, 0, OLDNODEDATA);
            }
            else 
               prev_speed = seg [i] -> bottom_pay.value;
         }
         else
            prev_speed = 0.0;

         if (seg [i] -> bottom_pay.expr) {
            T = Tension(seg [i] -> first_active -> Y[1], seg [i] -> material);
            speed = EvalCode(seg[i] -> bottom_pay.expr, 
                             seg[i] -> first_active, t, T, 
                             0, 0, 0, CURRNODEDATA);
         }
         else  {
            T = Tension(seg [i] -> first_active -> Y[1], seg [i] -> material);
            speed = seg [i] -> bottom_pay.value;
         }
         printf("t = %g, T_bot = %g, payspeed = %g\n", t, T, speed);

	    // we only need to act if speed is non-zero
 	 
         if (speed || prev_speed) {
            added = speed * dt - prev_speed * prev_dt;
            sf = 1.0 + seg[i]->first_active->Y[1];
            if (added > seg[i] -> bottom_spooled*sf) { // can't add more than we have left on spool
               added = seg[i] -> bottom_spooled*sf;
               speed = added / dt;
            }
            else if (added/sf < (seg[i] -> bottom_spooled + seg[i]->top_spooled) - (seg[i] -> length + seg[i] -> top_length + seg[i] -> bottom_length)) {
               added = ((seg[i] -> bottom_spooled + seg[i]->top_spooled) - (seg[i] -> length + seg[i] -> top_length + seg[i] -> bottom_length))*sf;
               speed = added / dt;
            }
 
            seg[i] -> bottom_spooled -= added/sf;
            seg[i] -> bottom_wet = seg[i] -> bottom_spooled * seg[i] -> material -> wet;

            if (added > 0) {

               if (seg[i] -> first_active -> ds == seg[i] -> first_active -> ds0)
                  first_spooled = seg [i] -> first_active -> number - 1; 
               else
                  first_spooled = seg [i] -> first_active -> number; 

               last_spooled = seg [i] -> first -> number; 
               if (first_spooled < last_spooled)
                  first_spooled = last_spooled;

               printf("t = %g, bottomAdd = %g, f = %d, l = %d\n", 
                      t, added, first_spooled, last_spooled); 

	   // first we add any "whole" segments we can

               j = first_spooled;
               while (j >= last_spooled) {
                  sf = 1.0 + node[j] -> Y[1];
                  if (node[j] -> ds + added/sf < node [j] -> ds0)
                     break;

                  if (!node [j] -> active) {
                     printf("activating %d\n", node [j] -> number);
                     CopyNodeData(node [j], seg[i] -> first_active, ne);
                     node [j] -> active = 1;
                     active_changed = 1;
                  }
                  printf("setting %d ds to ds0 (%g)\n", 
                         node[j] -> number, node[j] -> ds0);
                  added -= (node [j] -> ds0 - node[j] -> ds)*sf;
                  node [j] -> ds = node [j] -> ds0;
                  // node [j] -> Y [u_idx] = speed;

                  if (added == 0 || j == last_spooled)
                     break;

                  j --;
               }

               sf = 1.0 + node[j] -> Y[1];

		// take care of any fractional amount by adjusting ds

               if (j > 0) {
                  if (!node [j] -> active) {
                     CopyNodeData(node[j], seg[i] -> first_active, ne);
                     node [j] -> active = 1;
                     active_changed = 1;
                     node [j] -> ds += added/sf;
                     printf("activating %d with ds = %g\n", node [j] -> number, node[j] -> ds);
	                 // above is because the following else if does not
		             // work the first time we activate (ds == ds0 at first) 
                  }
                  else if (node [j] -> ds < node [j] -> ds0) {  
                     node [j] -> ds += added/sf;
                     printf("changing %d ds to %g (ds0 = %g)\n", node [j] -> number, node [j] -> ds, node[j] -> ds0);
                  }
               }

		// change the segment -> first_active if necessary
 
               if (node [j] != seg[i] -> first_active) {
                  pos = seg [i] -> first_active -> position;
                  seg [i] -> first_active -> position = Cable;
/*
                  next = seg[i] -> first_active -> next;
                  this = seg[i] -> first_active;
                  while(this != seg[i] -> last_active) {
                     this -> Y_o[3] = next -> Y_o[3];
                     this -> Y[3] = next -> Y[3];
                     this = this -> next;
                     next = this -> next;
                  }
*/
                    
#if 0
                  for (k = 1 ; k <= 6 ; k++) {
                  seg[i] -> first_active -> Y_o[k] = seg[i] -> first_active -> next -> Y_o[k] + (seg[i] -> first_active -> next -> Y_o[k] - seg[i] -> first_active -> next -> next -> Y_o[k])/seg[i] -> first_active -> next -> ds*seg[i] -> first_active -> ds;
                  printf("%d: %f %f %f\n", k,  
                         seg[i] -> first_active -> Y_o[k],
                         seg[i] -> first_active -> next -> Y_o[k],
                         seg[i] -> first_active -> next -> next -> Y_o[k]);
                  seg[i] -> first_active -> Y[k] = seg[i] -> first_active -> next -> Y[k] + (seg[i] -> first_active -> next -> Y[k] - seg[i] -> first_active -> next -> next -> Y[k])/seg[i] -> first_active -> next -> ds*seg[i] -> first_active -> ds;
                  }

                  // seg [i] -> first_active -> Y_o[3] = seg[i] -> first_active -> next -> Y_o[3];
                  // seg [i] -> first_active -> Y[3] = seg[i] -> first_active -> next -> Y[3];
#endif
                  seg [i] -> first_active -> pay   = 0.0;
                  seg [i] -> first_active -> pay_o = 0.0;
                  seg [i] -> first_active = node [j];
                  seg [i] -> first_active -> position = pos;
               }
            }
            else if (added < 0) {	// bottom_added is negative

               if (seg[i] -> first_active -> ds == 0)
                  first_unspooled = seg [i] -> first_active -> number + 1; 
               else
                  first_unspooled = seg [i] -> first_active -> number; 

               last_unspooled = seg [i] -> last_active -> number; 

               // if (first_spooled < last_spooled) // do we need something like this????
                 // first_spooled = last_spooled; 

               printf("t = %g, bottomAdd = %g, f = %d, l = %d\n", 
                       t, added, first_unspooled, last_unspooled); 

		// first we remove any "whole" segments we can

               j = first_unspooled;
               while (j <= last_unspooled) {
                  if (-added < node [j] -> ds*(1.0 + node [j] -> Y[1]))
                     break;

                  if (node [j] -> active) {
                     printf("deactivating %d\n", node [j] -> number);
                     // CopyNodeData(node [j], seg[i] -> first_active, ne);
                     node [j] -> active = 0;
                     active_changed = 1;
                  }
                  // node [j] -> ds = node [j] -> ds0;
                  added += node [j] -> ds*(1.0 + node [j] -> Y[1]);
                  node[j] -> ds = 0;
                  // node [j] -> Y [u_idx] = 0;


                  j ++;
                  if (added == 0)
                     break;
               }

               if (j > seg[i] -> last -> number)
                  j = seg[i] -> last -> number;

		       // take care of any fractional amount by adjusting ds

               if (added) {
                  node [j] -> ds += added/(1.0 + node[j] -> Y[1]);
                  printf("changing %d ds to %g (ds0 = %g)\n", node [j] -> number, node [j] -> ds, node[j] -> ds0);
               }

		// change the segment -> first_active if necessary
 
               if (node [j] != seg[i] -> first_active) {  // was != Jan 16 2015
                  pos = seg [i] -> first_active -> position;
                  seg [i] -> first_active -> position = Cable;
                  seg [i] -> first_active -> pay   = 0.0;
                  seg [i] -> first_active -> pay_o = 0.0;
                  CopyNodeData(node [j], seg[i] -> first_active, ne);

                  seg [i] -> first_active = node [j];
                  seg [i] -> first_active -> position = pos;
               }
            }
         }

         if (added) {
            seg[i] -> first_active -> pay = speed;
            seg[i] -> bottom_pay_speed = speed;
         }
         else {
            seg[i] -> first_active -> pay = 0.0;
            seg[i] -> bottom_pay_speed = 0;
         }
      }

      if (seg [i] -> num_top_dist) {
         printf("processing top spool segment %d\n", i);

         added = 0.0;
         if (t < prev_t) {
            if (seg [i] -> top_pay.expr) {
               T = Tension(seg [i] -> last_active -> Y_o[1], seg [i] -> material);
               prev_speed = EvalCode(seg [i] -> top_pay.expr, 
                                     seg[i] -> last_active, prev_t, T,
                                     0, 0, 0, OLDNODEDATA);
            }
            else 
               prev_speed = seg [i] -> top_pay.value;
         }
         else
            prev_speed = 0.0;

         if (seg [i] -> top_pay.expr) {
            T = Tension(seg [i] -> last_active -> Y[1], seg [i] -> material);
            speed = EvalCode(seg [i] -> top_pay.expr, 
                             seg[i] -> last_active, t, T, 
                             0, 0, 0, CURRNODEDATA);
                // seg[i] -> last_active -> Y[u_idx], 0.0, 0.0, 0.0, 0.0);
            printf("t = %g, T_top = %g, payspeed = %g\n", t, T, speed);
         }
         else  {
            speed = seg [i] -> top_pay.value;
         }

         if (speed || prev_speed) {
            added = speed * dt - prev_speed * prev_dt;

            sf = 1.0 + seg[i] -> last_active -> Y[1];
            if (added > seg[i] -> top_spooled*sf) { // can't add more than we have left on spool
               added = seg[i] -> top_spooled*sf;
               speed = added / dt;
            }
            else if (added/sf < (seg[i] -> bottom_spooled + seg[i]->top_spooled) - (seg[i] -> length + seg[i] -> top_length + seg[i] -> bottom_length)) {
               added = ((seg[i] -> bottom_spooled + seg[i]->top_spooled) - (seg[i] -> length + seg[i] -> top_length + seg[i] -> bottom_length))*sf;
               speed = added / dt;
            }

            seg[i] -> top_spooled -= added/sf;
            seg[i] -> top_wet = seg[i] -> top_spooled * seg[i] -> material -> wet;
            if (added > 0) {
               if (node [seg[i] -> last_active -> number - 1] -> ds == node[seg[i] -> last_active -> number - 1] -> ds0)
                  first_spooled = seg [i] -> last_active -> number + 1; 
               else
                  first_spooled = seg [i] -> last_active -> number; 

               last_spooled = seg [i] -> last -> number; 

               if (first_spooled > last_spooled)
                  first_spooled = last_spooled;

               printf("t = %g, topAdd = %g, f = %d, l = %d\n", 
                       t, added, first_spooled, last_spooled); 

               j = first_spooled;
               while (j <= last_spooled) {
                  sf = 1.0 + node[j-1] -> Y[1];
                  if (added/sf + node[j-1] -> ds < node [j-1] -> ds0)
                     break;

                  if (!node [j] -> active) {
                     printf("activating %d\n", node [j] -> number);
                     CopyNodeData(node[j], seg[i] -> last_active, ne);
                     node [j] -> active = 1;
                     active_changed = 1;
                  }
                  printf("setting %d ds to ds0 (%g)\n", 
                         node[j-1] -> number, node[j-1] -> ds0);
              
                  added -= (node[j-1]->ds0 - node [j-1] -> ds)*sf;
                  node [j-1] -> ds = node [j-1] -> ds0;

                  if (added == 0 || j == last_spooled)
                     break;

                  j++;
               }

               if (added) {
                  if (!node [j] -> active) {
                     printf("activating %d\n", node [j] -> number);
                     CopyNodeData(node[j], seg[i] -> last_active, ne);
                     node [j] -> active = 1;
                     active_changed = 1;
                     node [j] -> ds = 0;
                  }

	              // the following always works because ds != ds0 to start 
                  if (node [j-1] -> ds < node [j-1] -> ds0)  {
                     node [j-1] -> ds += added/(1.0 + node[j-1] -> Y[1]); 
                     printf("changing %d ds to %g (ds0 = %g)\n", 
                            node [j-1] -> number, node [j-1] -> ds, 
                            node[j-1] -> ds0);
                  }
               }
            }
            else if (added < 0) {	// top_added negative
               if (seg[i] -> last_active -> prev -> ds == 0)
                  last_unspooled = seg [i] -> last_active -> prev -> number - 1; 
               else
                  last_unspooled = seg [i] -> last_active -> prev -> number;
/*
               if (seg[i] -> last_active -> prev -> ds == 0)
                  last_unspooled = seg [i] -> last_active -> prev -> number;
               else
                  last_unspooled = seg [i] -> last_active -> number; 
*/
               first_unspooled = seg [i] -> first_active -> number; 

               printf("t = %g, topAdd = %g, f = %d, l = %d\n", 
                       t, added, first_unspooled, last_unspooled); 

               j = last_unspooled;
               while (j >= first_unspooled) {
                  if (-added < node [j] -> ds*(1.0 + node[j]->Y[1]))
                     break;

                  if (node [j+1] -> active) {
                     printf("deactivating %d\n", node [j+1] -> number);
                     node [j+1] -> active = 0;
                     active_changed = 1;
                  }
              
                  added += node [j] -> ds*(1.0 + node[j]->Y[1]);
                  node [j] -> ds = 0; // node [j-1] -> ds0;

                  j--;
                  if (added == 0 || j == first_unspooled)
                     break;

               }

               if (added) {
                  node [j] -> ds += added/(1.0 + node[j]->Y[1]); 
                  printf("changing %d ds to %g (ds0 = %g)\n", 
                         node [j] -> number, node [j] -> ds, node[j] -> ds0);
               }
               j ++;
            }

            if (node [j] != seg[i] -> last_active) {
               pos = seg [i] -> last_active -> position;
               seg [i] -> last_active -> position = Cable;
               seg [i] -> last_active -> pay = 0.0;
               seg [i] -> last_active -> pay_o = 0.0;
               CopyNodeData(node[j], seg[i] -> last_active, ne);
               seg [i] -> last_active = node [j];
               seg [i] -> last_active -> position = pos;
               seg [i] -> last_active -> ds = 0.0; 
            }
         }

         if (added) {
            seg [i] -> last_active -> pay = -speed;
            seg [i] -> top_pay_speed = -speed;
         }
         else {
            seg [i] -> last_active -> pay = 0.0;
            seg [i] -> top_pay_speed = 0.0;
         }
      }

      seg[i] -> prev_t = t;
      seg[i] -> prev_dt = dt;
   } // end loop over segments

 
   if (active_changed) {
      printf("num_active0 = %d\n", num_active);
      num_active = MakeNextPrevActive(p, node, num_nodes, active);
      printf("num_active1 = %d\n", num_active);
   }

   return num_active;
}
