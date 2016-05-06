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
# include "compress.h"
# include "problem.h"
# include "options.h"
# include "error.h"
# include "Tree.h"
# include "output.h"
# include "solve.h"
# include "allocate.h"
# include "control.h"
# include "results.h"

# define VERSION "3.12"

extern int SolutionDriver(Solution *s, int iterations);

int  output_map [MAXOUTPUT];

extern int unlink_input;
extern char *in_name;
			    
static int	default_static	    = 0;
static int 	default_motion      = 1;
static int 	default_vel         = 1;
static int 	default_force       = 1;
static int 	default_moment      = 1;
static int	default_euler	    = 1;
static double	default_sample	    = 0.0;
static double   default_snap_dt     = 0.0;
static double   default_seg_dt      = 0.0;
static double   default_buoy_dt     = 0.0;
static double   default_ext_dt     = 0.0;
static int	default_decimate    = 1;
static int	default_debug 	    = 0;
static int	default_twoD	    = 1;
static int	default_X	    = 0;
static int	default_quiet	    = 0;
static int	default_quit	    = 0;
static int      default_iterations  = 0;
static int      default_refine      = 0;
static int	default_auto	    = 0;
static  double default_progress_dt = 0.0;
static  double default_restart_t  = 0.0;


# ifndef WINGUI

static char *usage_msg = "\
usage: cable [-in model_description_file] [-out output_file] [options]\n\
   -static       only do the static solution (default=off)\n\
   -simple       follow the static solution with simple dynamics (default=off)\n\
   -load ...     name of file to get static solution from (default=none)\n\
   -initial ...  load static solution in file as initial guess (default=none)\n\
   -twoD	 use a two-dimensional solver (default=on)\n\
   -nodes ...    node numbers to give time histories for\n\
   -first        include first node in list for node time histories\n\
   -last         include last node in list for node time histories\n\
   -terminals    include nodes at terminals in list for node time histories\n\
   -connectors   include nodes w/connectors in list for node time histories\n\
   -sample ...   time step to sample output data (default=problem dt)\n\
   -snap_dt ...  time step to dump configuration snapshots (default=0.0)\n\
   -seg_dt ...   time step for segment data (pay/length/spooled) (default=0.0)\n\
   -buoy_dt ...  time step for saving LAMP 6-axis buoy motion (default=0.0)\n\
   -ext_dt ...   time step for saving external forces (default=0.0)\n\
   -decimate ... nodal increment to using in writing spatial results (default=1)\n\
   -dynstat ...  write final dynamic state as a static solution (default=none)\n\
   -progress ... progress file written for later restarts (default=none)\n\
   -progress_dt ... time step for saving progress (default=0.0)\n\
   -restart ...   progress file read for restart (default=none)\n\
   -restart_t ... time value to load from progress file for restart\n\
   -motion       include motion results in output file (default=on)\n\
   -vel          include velocity results in output file (default=on)\n\
   -force        include force results in output file (default=on)\n\
   -moment       include bending moment results in output file (default=on)\n\
   -euler	 include euler parameter results in output file (default=on)\n\
   -refine       use loaded static solution as mesh to refine (default=off)\n\
   -X            use graphical display and control (default=off)\n\
   -debug        interpret input and exit\n\
   -bom          generate a bill of materials and exit\n\
   -version      display version information and exit\n\
   -quiet        suppress all display information (for speed) (default=off)\n\
   -quit         quit when solution is complete (default=off)\n\
   -cpp		 name of C pre-processor to use (default=/lib/cpp)\n\
   -nocpp        do not run input file through C pre-processor (default=off)\n\
   -I		 add a directory to the C pre-processor include path\n\
   -D            define a macro for the C pre-processor\n\
   -U            undefine a macro for the C pre-processor\n\
";

# endif

static int usage(msg)
   char	*msg;
{
# ifdef WINGUI
   char  buffer [1024];

   sprintf (buffer,"invokation error: %s", msg);
   WinErrorDialog (buffer);

   if (unlink_input)
      unlink(in_name);

   exit (0);

# else

   fprintf (stderr,"cable: %s\n\n", msg);
   fputs(usage_msg, stderr);
   exit (0);

# endif

}

# ifdef WINGUI
int 
CableMain (argc, argv)
# else
int 
main (argc, argv)
# endif
   int  argc;
   char *argv [ ];
{
   Solution     solution;
   static char *pipe_name = "-";
   char	       *load_name;
   char	       *initial_name;
   char	       *dynstat_name;
   char	       *tmpdir;
   char        *out_name;
   char        *table_name;
   char        *in_name;
   char        *progress_name;
   char        *restart_name;
   int		motion;
   int		vel;
   int		force;
   int		moment;
   int		euler;
   int		static_only;
   int		simple;
   int      bom;
   int		debug_input;
   int		first;
   int		last;
   int		terminals;
   int		connectors; 
   int		X;
   int		quiet;
   int		quit;
   int		twoD;
   int	       *output_nodes;
   int          num_output_nodes;
   int		iterations;
   int		auto_static;
   int		refine;
   double	sample;
   double	snap_dt;
   double   seg_dt;
   double	buoy_dt;
   double   ext_dt;
   double   progress_dt, restart_t;
   double   rot;
   int	 	decimate;
   int	        n;
   int		i;
   int		version;
   int      unlink_input;

	/*
	 * parse the cpp specific options and remove them
	 * from the command line
	 */

    if (ParseCppOptions (&argc, argv)) 
        usage("error in cpp specific options");

	/*
	 * fetch the parameters
	 */

   unlink_input = 0;
   n = GetBooleanOption (argc, argv, "unlink", &unlink_input);
   if (n < 0)
      usage ("error in unlink parameter");

   version = 0;
   n = GetBooleanOption (argc, argv, "version", &version);
   if (n < 0)
      usage ("error in version parameter");
 
   if (version) {
      fprintf (stderr,"WHOI Cable ver. %s (built %s, %s)\n", VERSION, __DATE__, __TIME__);
      exit (0);
   }

   n = GetSoloStringOption(argc, argv, "in", &in_name);
   if (n < 0)
      usage("error in input file parameter (-in)");
   else if (n == 0)
      in_name = pipe_name;

   debug_input = default_debug;  
   n = GetBooleanOption (argc, argv, "debug", &debug_input);
   if (n < 0)
      usage ("error in debug parameter");

   bom = 0;
   n = GetBooleanOption (argc, argv, "bom", &bom);
   if (n < 0)
      usage ("error in bom parameter");

   out_name = NULL;
   n = GetSoloStringOption(argc, argv, "out", &out_name);
   if (n != 1 && !debug_input && !bom)
      usage("must specify an output file (-out)");

   table_name = NULL;
   n = GetSoloStringOption(argc, argv, "table", &table_name);

   if (out_name && strcmp(out_name, in_name) == 0)
      usage("input and output files must have different names");

   load_name = NULL;
   n = GetSoloStringOption(argc, argv, "load", &load_name);
   if (n < 0)
      usage ("error in load file option");

   if (load_name && strcmp(out_name, load_name) == 0)
      usage("output and load files must have different names");

   initial_name = NULL;
   n = GetSoloStringOption(argc, argv, "initial", &initial_name);
   if (n < 0)
      usage ("error in initial file option");

   if (load_name && strcmp(out_name, load_name) == 0)
      usage("output and load files must have different names");

   dynstat_name = NULL;
   n = GetSoloStringOption(argc, argv, "dynstat", &dynstat_name);
   if (n < 0)
      usage ("error in dynstat file option");

   if (dynstat_name && strcmp(out_name, dynstat_name) == 0)
      usage("output and dynstat files must have different names");

   progress_name = NULL;
   n = GetSoloStringOption(argc, argv, "progress", &progress_name);
   if (n < 0)
      usage ("error in progress file option");

   if (progress_name && strcmp(out_name, progress_name) == 0)
      usage("output and progress files must have different names");

   restart_name = NULL;
   n = GetSoloStringOption(argc, argv, "restart", &restart_name);
   if (n < 0)
      usage ("error in restart file option");

   if (restart_name && strcmp(restart_name, progress_name) == 0)
      usage("restart and progress files must have different names");

   static_only = default_static;
   n = GetBooleanOption (argc, argv, "static", &static_only);
   if (n < 0)
      usage ("error in static parameter");

   simple = 0;
   n = GetBooleanOption (argc, argv, "simple", &simple);
   if (n < 0)
      usage ("error in simple parameter");

   refine = default_refine;
   n = GetBooleanOption (argc, argv, "refine", &refine);
   if (n < 0)
      usage ("error in refine parameter");

   auto_static = default_auto;
   n = GetBooleanOption (argc, argv, "auto", &auto_static);
   if (n < 0)
      usage ("error in auto parameter");

   motion = default_motion;
   n = GetBooleanOption (argc, argv, "motion", &motion);
   if (n < 0)
      usage ("error in motion parameter");

   vel = default_vel;
   n = GetBooleanOption (argc, argv, "vel", &vel);
   if (n < 0)
      usage ("error in vel parameter");

   force = default_force;
   n = GetBooleanOption (argc, argv, "force", &force);
   if (n < 0)
      usage ("error in force parameter");

   moment = default_moment;
   n = GetBooleanOption (argc, argv, "moment", &moment);
   if (n < 0)
      usage ("error in moment parameter");

   euler = default_euler;
   n = GetBooleanOption (argc, argv, "euler", &euler);
   if (n < 0)
      usage ("error in euler parameter");

   iterations = default_iterations;  
   n = GetBooleanOption (argc, argv, "iterations", &iterations);
   if (n < 0)
      usage ("error in iterations parameter");

   tmpdir = NULL;
   n = GetSoloStringOption (argc, argv, "tmpdir", &tmpdir);
   if (n < 0)
      usage ("error in tmpdir parameter");

# ifndef WINGUI
   X = default_X;  
   n = GetBooleanOption (argc, argv, "X", &X);
   if (n < 0)
      usage ("error in X parameter");
# else
   X = 1;
# endif

   quiet = default_quiet;  
   n = GetBooleanOption (argc, argv, "quiet", &quiet);
   if (n < 0)
      usage ("error in quiet parameter");

   quit = default_quit;  
   n = GetBooleanOption (argc, argv, "quit", &quit);
   if (n < 0)
      usage ("error in quit parameter");

   twoD = default_twoD;  
   n = GetBooleanOption (argc, argv, "twoD", &twoD);
   if (n < 0)
      usage ("error in twoD parameter");

   sample = default_sample;
   n = GetSoloDoubleOption (argc, argv, "sample", &sample);
   if (n < 0)
      usage ("error in sample parameter");

   seg_dt = default_seg_dt;
   n = GetSoloDoubleOption (argc, argv, "seg_dt", &seg_dt);
   if (n < 0)
      usage ("error in seg_dt parameter");
 
   snap_dt = default_snap_dt;
   n = GetSoloDoubleOption (argc, argv, "snap_dt", &snap_dt);
   if (n < 0)
      usage ("error in snap_dt parameter");
 
   progress_dt = default_progress_dt;
   n = GetSoloDoubleOption (argc, argv, "progress_dt", &progress_dt);
   if (n < 0)
      usage ("error in progress_dt parameter");
 
   restart_t = default_restart_t;
   n = GetSoloDoubleOption (argc, argv, "restart_t", &restart_t);
   if (n < 0)
      usage ("error in restart_t parameter");
 
   rot = 0.0;
   n = GetSoloDoubleOption (argc, argv, "rot", &rot);
   if (n < 0)
      usage ("error in rot parameter");
 
   buoy_dt = default_buoy_dt;
   n = GetSoloDoubleOption (argc, argv, "buoy_dt", &buoy_dt);
   if (n < 0)
      usage ("error in buoy_dt parameter");

   ext_dt = default_ext_dt;
   n = GetSoloDoubleOption (argc, argv, "ext_dt", &ext_dt);
   if (n < 0)
      usage ("error in ext_dt parameter");

   decimate = default_decimate;
   n = GetSoloIntegerOption (argc, argv, "decimate", &decimate);
   if (n < 0)
      usage ("error in decimate parameter");
 
   output_nodes = NULL;
   num_output_nodes = GetIntegerOption (argc, argv, "nodes", &output_nodes);
   if (num_output_nodes < 0)
      usage ("error in nodes parameter");
 
   first = 0;
   n = GetBooleanOption (argc, argv, "first", &first);
   if (n < 0)
      usage ("error in first parameter");
 
   terminals = 0;
   n = GetBooleanOption (argc, argv, "terminals", &terminals);
   if (n < 0)
      usage ("error in terminals parameter");
 
   last = 0;
   n = GetBooleanOption (argc, argv, "last", &last);
   if (n < 0)
      usage ("error in last parameter");
 
   connectors = 0;
   n = GetBooleanOption (argc, argv, "connectors", &connectors);
   if (n < 0)
      usage ("error in connectors parameter");

# ifdef GUI
   if (X)  {
      fprintf(stderr,"calling create\n");
      CreateControlDialog (&argc, argv, quiet);
   } 
   DisplayMode (X, quiet);
# elif defined (WINGUI)
   DisplayMode (X, quiet);
# else
   DisplayMode (0, quiet);
# endif

   if (ArgsUsed () != argc - 1) 
      usage ("incorrect parameter syntax");

	/*
	 * make a little map of what kinds of output variables
	 * are being requested 
	 */

   for (i = 1 ; i < MAXOUTPUT ; i++)
      solution.output_map [i] = 0;

   solution.output_map [MOTION] = motion;
   solution.output_map [VEL]    = vel;
   solution.output_map [FORCE]  = force;
   solution.output_map [MOMENT] = moment;
   solution.output_map [EULER]  = euler;

   solution.rotation = rot;
   solution.in_name = in_name;
   solution.table_name = table_name;
   solution.out_name = out_name; 
   solution.dynstat_name = dynstat_name;
   solution.progress_file = progress_name;
   solution.restart_file  = restart_name;
   solution.progress_dt = progress_dt;
   solution.restart_t   = restart_t;
   solution.static_only = static_only;
   solution.twoD   = twoD;
   solution.static_file = load_name;
   solution.initial_file = initial_name;
   solution.output_dt = sample;
   solution.snapshot_dt = snap_dt;
   solution.segment_dt = seg_dt;
   solution.buoy_dt = buoy_dt;
   solution.ext_dt  = ext_dt;
   solution.tmpdir = tmpdir;
   solution.debug_input = debug_input;
   solution.bill_of_materials = bom;
   solution.refine = refine;
   solution.simple = simple;
   solution.auto_static = auto_static;
   solution.quit = quit;
   solution.quiet = quiet;
   solution.unlink_input = unlink_input;
   solution.X = X;
   solution.decimate = decimate;

   solution.userQuit = 0;
  
   solution.output_nodes = output_nodes;
   solution.num_output_nodes = num_output_nodes;  
   solution.output_special = first*OUTPUT_FIRST + last*OUTPUT_LAST 
                             + terminals*OUTPUT_TERMINALS 
                             + connectors*OUTPUT_CONNECTORS;

   solution.problem = (Problem *) calloc(1, sizeof(Problem));
   solution.environment = (Environment *) calloc(1, sizeof(Environment));
   solution.analysis = (Analysis *) calloc(1, sizeof(Analysis));
   
   return SolutionDriver(&solution, 0);
}
