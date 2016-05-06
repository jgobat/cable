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
 * File:	driver.c
 *
 * Description: the driver module for the cable simulation application
 *
 * History:
 *		
 ****************************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <malloc.h>
# include <unistd.h>
# include <math.h>
# include <sys/stat.h>
# include "compress.h"
# include "problem.h"
# include "options.h"
# include "error.h"
# include "Tree.h"
# include "output.h"
# include "solve.h"
# include "allocate.h"
# include "control.h"
# include "segments.h"
# include "tension.h"

# define VERSION "3.12"

Problem     *problem;
Environment *environment;
Analysis    *analysis;
Debug         debug;

int		 unlink_input;
int		 static_finished = 0; 
char		*in_name;
char		*out_name = NULL;

extern int AutoGenerateMesh(Problem *, Node *);
extern void MaterialsBill(FILE *, Problem *, Environment *);

extern void TensionToGlobal(Node node, int twoD, double *Fx, double *Fy, double *Fz);

static double
sf(Node n, double level)
{
    if (level)
        return NodeTension(n)/level;
    else
        return 0.0;
}


int
SolutionDriver(Solution *s, int iterations)
{
   int      dynamic;
   int      i, n;
   int		nn;
   int		nan;
   int	 	num_added;
   Node	     *node;
   Node	     *active;
   int		(*dynamic_solve) ( );
   int      (*static_solve) ( );
   ResFile      out;
   ResFile	dynstat;
   int		sing;
   int		status;
   int terminals, first, last, connectors;
   FILE     *fp;
   char     *text;
   struct stat  st;

   
   unlink_input = s -> unlink_input;
   in_name = s -> in_name;
   out_name = s -> out_name;

	/*
	 * set the temporary directory in case we need it for error
 	 * buffering
	 */

   if (s -> tmpdir)
      SetTempDir (s -> tmpdir);

   problem = s -> problem;
   environment = s -> environment;
   analysis = s -> analysis;

   // set the parent solution structure on the "children" so we can 
   // always get back to it
   problem -> solution = environment -> solution = analysis -> solution = s;
   s -> userQuit = 0;

	/*
	 * read the input file and dump debugging or BOM output if requested
	 */

# if defined (GUI) || defined (WINGUI)
   if (s -> X) 
      BufferErrors (1);

   n = ReadModelFile (s -> in_name, problem, analysis, environment, 0);

   if (s -> X) 
      BufferErrors (0);

   if (n) 
      return 1;

# else
   n = ReadModelFile (s -> in_name, problem, analysis, environment, 0);
   if (n) 
      return 1;
# endif

   FillInEnvironment(environment);
   FillInAnalysis(problem, analysis);
   FillInObjectProperties(problem, environment);

//   FillInEnvironment(environment);
//   printf("pre-fill: %g\n", analysis -> tolerance);
//   FillInAnalysis(problem, analysis);
//   FillInObjectProperties(problem, environment);

   if (s-> debug_input) {
      fDumpModelFile(stdout, problem, analysis, environment);
      Exit(0);
   }
   else if (s -> bill_of_materials) {
      if (out_name) {
         fp = fopen(out_name, "w");        
      }
      else
         fp = NULL;

      MaterialsBill(fp, problem, environment); 
      if (fp)
         fclose(fp);

      return 0;
   } 


	/*
	 * polish off the problem definition work that isn't
	 * actually done by the parser -- build the array of
	 * nodes from the segment tree and calculcate derived
	 * material and body properties
	 */


   // problem -> branch  = BuildBranchArray(problem, &(problem -> num_branch));
   // CreateNodeArray (&node, &nn, &active, &nan);

   node = problem -> node;
   nn = problem -> num_nodes;
   nan = problem -> num_active;
   active = problem -> active;

   if (nn == 0 || nan == 0) {
       return 1;
   }

   if (s -> num_output_nodes) {
      for (i = 0 ; i < s -> num_output_nodes ; i++) {
         n = s -> output_nodes [i];
         
         if (n < 1 || n > nn)
            error("listed node %d not defined in current problem", n);
      }
   }
   else
      s -> output_nodes = NULL;

   first      = s -> output_special & OUTPUT_FIRST;
   last       = s -> output_special & OUTPUT_LAST;
   terminals  = s -> output_special & OUTPUT_TERMINALS;
   connectors = s -> output_special & OUTPUT_CONNECTORS;

   if ((first && last) || terminals) {
      s -> num_output_nodes += 2;
      Reallocate(s -> output_nodes, int, s -> num_output_nodes);
      s -> output_nodes [s -> num_output_nodes - 2] = 1;
      s -> output_nodes [s -> num_output_nodes - 1] = nn;
   }
   else if (first) {
      s -> num_output_nodes += 1;
      Reallocate(s -> output_nodes, int, s -> num_output_nodes);
      s -> output_nodes [s -> num_output_nodes - 1] = 1;
   }   
   else if (last) {
      s -> num_output_nodes += 1;
      Reallocate(s -> output_nodes, int, s -> num_output_nodes);
      s -> output_nodes [s -> num_output_nodes - 1] = nn;
   }   

   if (connectors || terminals) {
      num_added = 0;
      
      for (i = 1 ; i <= nn ; i++) {
         if (connectors && node [i] -> segment -> last -> number == i 
                        && node [i] -> segment -> connector) {
            num_added += 2;
         }
         else if  (terminals && node [i] -> position == BranchTerminal) {
            num_added ++;
         }
      } 

      if (num_added) {
         s -> output_nodes = Reallocate(s -> output_nodes, int, (s -> num_output_nodes + num_added));

         num_added = 0;
          
         for (i = 1 ; i <= nn ; i++) {
            if (connectors && node [i] -> segment -> last -> number == i 
                           && node [i] -> segment -> connector) {
               s -> output_nodes [num_added + s -> num_output_nodes] = i;
               s -> output_nodes [num_added + s -> num_output_nodes + 1] = i+1;
               num_added += 2;
            }
	        else if (terminals && node [i] -> position == BranchTerminal) {
               s -> output_nodes [num_added + s -> num_output_nodes] = i;
               num_added ++;
            }
         }

         s -> num_output_nodes += num_added;
      }
   }

# if defined (GUI) || defined (WINGUI)
   if (s -> X)
      BufferErrors (1);
# endif

   n = CheckTypeParameters ( );
   n += CheckStaticParameters ( );

   if (!s -> static_only) {
      n += CheckDynamicParameters ( );
      n += CheckEnvironmentParameters ( );
   }
 
   n += CheckMaterialProperties (s -> twoD);
   n += CheckBuoyProperties ( );

   n += CheckBranchTerminalProperties ( );

   if (n)
      error ("%d errors found in problem description", n);

# if defined (GUI) || defined (WINGUI)
   if (s -> X)
      BufferErrors (0);
# endif

   if (n)
      return 1;

   if (s -> output_nodes == NULL && s -> output_dt && !s -> static_only) {
      error("you gave a sampling rate but no nodes (-nodes) to sample");
      return 1;
   }
   if (s -> output_nodes != NULL && !s -> output_dt && !s -> static_only) {
      error ("you gave nodes to sample but no sampling rate (-sample)");
      return 1;
   }
   if (s -> output_nodes == NULL && !s -> snapshot_dt && !s -> static_only) {
      error ("no nodes to sample or snapshot rate given");
      return 1;
   }

   if (!s -> static_only && s -> output_dt && fabs(s -> output_dt - (int)((s -> output_dt + 1e-6)/analysis -> dt)*analysis -> dt) > 1e-6) {
      error("sample rate must be an integer multiple of problem time step");
      return 1;
   }

   if (!s -> static_only && s -> snapshot_dt && fabs(s -> snapshot_dt - (int)((s -> snapshot_dt + 1e-6)/analysis -> dt)*analysis -> dt) > 1e-6) {
      error ("snap_dt must be an integer multiple of problem time step");
      return 1;
   }
   if (s -> decimate < 0) {
      error("decimate must be an integer >= 1");
      return 1;
   }
   if (s -> decimate > 1 && (nn - 1) % s -> decimate > 0)  {
      error("decimate must be specified to include nodes 1 and n");
      return 1;
   }
  
# if defined (GUI) || defined (WINGUI)
   if (s -> X && !s -> quiet) {
      fprintf(stderr,"initializing control\n");
      ControlDialogInitialize (analysis, s -> controls, s -> static_only ? PLOT_SPACE : PLOT_SPACE | PLOT_TIME);
   }
# endif
 
   // read in a copy of the input file text to stash into the results file

   if (stat(s -> in_name, &st) == 0) {
      text = (char *) malloc(sizeof(char) * (st.st_size + 1));
      fp = fopen(s -> in_name, "rb");
      fread(text, sizeof(char), st.st_size, fp);
      text[st.st_size] = 0;
      fclose(fp);
   }
   else {
      text = NULL;
   }

   // initialize the results file

   out = res_open (s -> out_name, "wb");
   if (out == NULL)  {
      error("could not open output file %s for writing", s -> out_name);
      return 1;
   }


   if (iterations) {
      debug.status           = 1;
      debug.sample_it        = s -> output_dt;
      debug.snap_it          = s -> snapshot_dt;
      debug.output_map       = s -> output_map;
      debug.num_output_nodes = s -> num_output_nodes;
      debug.output_nodes     = s -> output_nodes; 
      debug.out              = out;
      debug.decimate	     = s -> decimate;
      s -> static_only = 0;
   }
   else
      debug.status = 0;

   if (s -> static_only || environment -> forcing != LAMP) 
      s -> buoy_dt = 0.0;

   if (s -> static_only)
      s -> ext_dt = 0;

   dynamic = 0;
   if (!s -> static_only)  {
      dynamic |= 1;
      dynamic |= 16; // all files now have start t in them
   }
   if (s -> buoy_dt)
      dynamic |= 2;
   if (s -> segment_dt)
      dynamic |= 4;
   if (s -> ext_dt)
      dynamic |= 8;

   InitializeResultsFile(out, nn, problem -> title, s -> output_map, 
                         s -> decimate, dynamic, s -> twoD, text);
   if (debug.status)
      s -> static_only = 1;

   if (text) {
      free(text);
   }

	/*
	 * solve the static problem and just bail out afterwards
	 * if that is all we are doing for now
	 */

   if (s -> twoD) {
      if (analysis -> static_solution == Shooting)
         static_solve = ShootStaticProblem2D;
      else
         static_solve = SolveStaticProblem2D;

      dynamic_solve = SolveDynamicProblem2D;
   }
   else {
      if (analysis -> static_solution == Shooting) 
         static_solve = ShootStaticProblem3D;
      else
         static_solve = SolveStaticProblem3D;
   
      dynamic_solve = SolveDynamicProblem3D;
   }
 

   // static option 1: load static from file

   if (s -> static_file) {
      status = LoadStaticSolution (s -> static_file, node, nn, s -> twoD, s -> rotation);

      if (status) {
         res_close(out);
         unlink(s -> out_name);
         error("could not load static solution from file");
         return 1;
      }

      if (s -> refine && s -> twoD && analysis -> mesh_amplify) {
         i = AutoGenerateMesh (problem, node);
         if (!i)
            status = static_solve (0, node, nn, active, nan, out, s -> output_map);
      }

      if (!status)
         WriteStaticSolution (node, nn, out, s -> output_map, s -> decimate, s -> twoD);
   }

   // static option 2: automatic static solve

   else if (s -> auto_static) {
      status = AutoStaticSolve(node, nn, active, nan, s -> twoD, s -> initial_file);

      if (!status)
         WriteStaticSolution(node, nn, out, s -> output_map, s -> decimate, s -> twoD);
   }

   // static option 3 (most common): solve for static solution

   else {

      // option 3a: get initial guess from a file

      if (s -> initial_file) {
         status = LoadStaticSolution(s -> initial_file, node, nn, s -> twoD, s -> rotation);

         if (status) {
            res_close(out);
            unlink(s -> out_name);
            error("could not load initial guess static solution from file");
            return 1;
         }
      }	
 
      if (static_solve) {
         status = static_solve(s -> initial_file ? 1 : 0, node, nn, 
                               active, nan, out, s -> output_map);
         fprintf(stderr,"static solved\n");
      }

      if (s -> userQuit) {
         ControlDialogQuit(s);
         return 1;
      }

      // option 3b: if we got a solution we can optionally regrid
      // and then resolve

      if (!status && analysis -> mesh_amplify && s -> twoD) {
         i = AutoGenerateMesh (problem, node);
         if (!i) {
            status = static_solve (NULL, node, nn, active, nan, out, s -> output_map);
            if (s -> userQuit) {
               ControlDialogQuit(s);
               return 1;
            }
         }
      }

      if (!debug.status) {
         if (status) {
            error("could not get static solution");
            return 1;
         }

         WriteStaticSolution (node, nn, 
		      	              out, s -> output_map, s -> decimate, s -> twoD);

        if (s -> table_name) {
            FILE *fp;
            int   i, ns;
            double T, maxT, Fx, Fy, Fz, Wmin, Sn, phi;
            Segment *seg;
            Node n, max_n;

            fp= fopen(s -> table_name, "w");
            if (!fp) {
                fprintf(stderr,"file open error table\n");
            } 
            T = NodeTension(node[nn]);
            TensionToGlobal(node[nn], 1, &Fx, &Fy, &Fz);
            fprintf(fp, "%s %f %f %f %f %f %f 0\n",
                    problem -> terminal[2] -> buoy -> name, 
                    node[nn] -> s,
                    environment -> depth - node[nn] -> x,
                    node[nn] -> Y[3]*180/M_PI,
                    T/4.4482216, 
                    sf(node[nn], node[nn] -> segment -> material -> swl),
                    sf(node[nn], node[nn] -> segment -> material -> yield));


            // build local version of seg array without branches
            seg = BuildSegmentArray(problem, &ns, 0);
            for (i = ns ; i >= 1 ; i--) {
                maxT = 0;
                n = seg[i] -> last;
                max_n = n;
                while (n && n -> segment == seg[i]) {
                    T = NodeTension(n);
                    if (T > maxT) {     
                        maxT = T;
                        max_n = n;
                    }
                    n = n -> prev;
                }

                n = seg[i] -> last;

                TensionToGlobal(max_n, 1, &Fx, &Fy, &Fz);
                fprintf(fp, "%s %f %f %f %f %f %f 0\n",
                        seg[i] -> material -> name, 
                        n -> s, 
                        environment -> depth - n -> x,
                        n -> Y[3]*180/M_PI,
                        maxT/4.4482216, 
                        sf(max_n, n -> segment -> material -> swl),
                        sf(max_n, n -> segment -> material -> yield));

            }

            if (problem -> terminal[1] -> anchor) {
                T = NodeTension(node[1]); 
                TensionToGlobal(node[1], 1, &Fx, &Fy, &Fz);
                phi = node[1] -> Y[3];
                Sn  = node[1] -> Y[2];
                Wmin = ((T*sin(phi) + Sn*cos(phi))*problem -> terminal[1]->safety/problem -> terminal[1]->friction +  T*cos(phi) - Sn*sin(phi))/4.4482216;
                fprintf(fp, "%s %f %f %f %f %f %f %f\n",
                    problem -> terminal[1] -> anchor -> name, 
                    node[1] -> s,
                    environment -> depth - node[1] -> x, 
                    node[1] -> Y[3]*180/M_PI,
                    T/4.4482216, 
                    sf(node[1], node[1] -> segment -> material -> swl),
                    sf(node[1], node[1] -> segment -> material -> yield), Wmin);

                
            }

            fclose(fp);
        }
      }
   }
   static_finished = 1;
   if (status) {
      fprintf(stderr, "numerical error - deleting static result file\n");
      res_close(out);
      unlink(s -> out_name);

      if (!s -> quit)
         Exit (1);
      else
         return 1;
   }

   if (s -> simple) 
      SimpleDynamics(active, nan);
   
   if (s -> static_only) {
      res_close (out);
      DisplayMessage ("Static solution complete");  

#if (defined GUI || defined WINGUI)
      ControlPlotSnaps(problem, environment);
      ControlTabulateResults(problem, environment);
#endif

      if (!s -> quit)
         Exit (0);
      else
         return 0;
   }
   
   if (s -> dynstat_name) 
      problem -> dynstat = 1;

#if (defined GUI || defined WINGUI)
   ControlTabulateResults(problem, environment);
#endif

   sing = dynamic_solve (node, nn, active, nan, out, s -> output_map, 
                         s -> output_nodes, s -> num_output_nodes, 
			             s -> output_dt, s -> snapshot_dt, 
                         s -> segment_dt, s -> buoy_dt, s -> ext_dt, 
                         s -> decimate, 
                         s -> progress_file, s -> progress_dt, 
                         s -> restart_file, s -> restart_t);

   res_close (out);
   if (s -> userQuit) {
      ControlDialogQuit(s);
      return 1;
   }

   if (s -> dynstat_name) {
      dynstat = res_open (s -> dynstat_name, "wb");
      if (dynstat == NULL)  {
         error("could not open dynstat file %s for writing", s -> dynstat_name);
         return 1;
      }

      InitializeResultsFile (dynstat, nn, problem -> title, s -> output_map, 1, 0, s -> twoD, NULL);

      WriteStaticSolution (node, nn, dynstat, s -> output_map, 1, s -> twoD);
      res_close(dynstat);
   }

   DisplayMessage ("Dynamic solution complete");

   if (!s -> quit)
      Exit(0);
 
   return 0;		
}
