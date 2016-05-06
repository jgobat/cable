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

/************************************************************************
 * File:	problem.c						*
 *									*
 * Description:	This file contains the public and private function	*
 *		definitions for the creation of a problem instance.	*
 ************************************************************************/

# include <stdio.h>
# include <string.h>
# include <unistd.h>
# include <errno.h>
# include <sys/stat.h>
# include "error.h"
# include "compress.h"
# include "problem.h"
# include "allocate.h"
# include "objects.h"
# include "segments.h"

# define streq(a,b)	!strcmp(a,b)
# define strneq(a,b)	strcmp(a,b)

// # ifndef LIBDIR
// # define LIBDIR "/usr/local/lib"
// # define LIBDIR  NULL
// # endif

#ifdef WINDOWS
#define popen _popen
#define pclose _pclose
#endif

# ifndef CPP

# ifdef WINDOWS
# define CPP NULL
# elif defined __CYGWIN32__
# define CPP "/usr/bin/cpp"
# elif defined __APPLE__
# define CPP "/usr/bin/cpp"
# else
# define CPP "/usr/bin/cpp"
# endif

# endif

static char *cpp = CPP;
static char  cpp_command [2048];


ParserControl ctl;

/************************************************************************
 * Function:	segment_cmp						           
 *									
 * Description:	Compares two segments by number.			
 ************************************************************************/

static int segment_cmp (s1, s2)
    Item s1;
    Item s2;
{
    return ((Segment) s1) -> number - ((Segment) s2) -> number;
}

/************************************************************************
 * Function:	branch_cmp						
 *									
 * Description:	Compares two branches by number.			
 ************************************************************************/

static int branch_cmp (b1, b2)
    Item b1;
    Item b2;
{
    return ((Branch) b1) -> number - ((Branch) b2) -> number;
}

/************************************************************************
 * Function:	ObjectCompare  						                    *
 *									                                    *
 * Description:	Compares two materials by category and then name		*
 ************************************************************************/

int 
ObjectCompare (Item item1, Item item2)
{
    CableObject      obj1 = (CableObject) item1;
    CableObject      obj2 = (CableObject) item2;
    int         i;

    // category sort doesn't work - given a red-black tree we can
    // only sort on one key
/*
    if (obj1 -> category && obj2 -> category 
        && (i = strcmp(obj1 -> category, obj2 -> category))) {
        if (i == 0) printf("matched on category\n");
        return i;
    }
*/
    // must have the same category - sort by name

    i = strcmp(obj1 -> name, obj2 -> name);
    return i;
}

/************************************************************************
 * Function:	resolve_segment	
 *								
 * Descriptin:	Resolve the names of objects for a segment.		
 ************************************************************************/

static int 
resolve_segment (Item item, void *call_data)
{
    Problem *p = (Problem *) call_data;
    struct material    m;
    struct connector   c;
    unsigned	       number;
    Tree	       tree;
    Segment	       segment;
    int		       i;
    int            n;


    segment = (Segment) item;
    number = segment -> number;

    /* Resolve the material. */

    tree = p -> material_tree;
    m.name = (char *) segment -> material;

    n = 0;

    if (m.name) {
	    segment -> material = (Material) TreeSearch (tree, &m);
       
	    if (!segment -> material) {
	        error ("segment %u uses undefined material %s", number, m.name);
            n ++;
        } 
	    Deallocate (m.name);
    } 
    else 
        error ("segment %u has no material assigned", number);


    if (!segment -> length && segment -> material -> length) {
        segment -> length = segment -> material -> length;

        if (!segment -> num_dist) {
           segment -> num_dist = 1;
           segment -> dist [1].nodes = 3;
	       segment -> dist [1].percent = 1.0;
	    }
    }

    if (!segment -> length) {
        error ("segment %u has no length assigned", number);
        n ++;
    }

    if (!segment -> num_dist) {
        error ("segment %u has no nodes defined", number);
        n ++;
    }


    /* Resolve any attachments */

    tree = p -> connector_tree;

    for (i = 1 ; i <= segment -> num_attach ; i++) {
       c.name = (char *) segment -> attach [i].object;

       if (c.name) {
          segment -> attach [i].object = (Connector) TreeSearch(tree, &c);

          if (!segment -> attach [i].object) {
             error ("layout uses undefined attachment %s", c.name);
             n ++;
          }
          Deallocate (c.name);
       }
       else
          segment -> attach [i].object = NULL;
    }


    /* Resolve the optional connector. */

    tree = p -> connector_tree;
    c.name = (char *) segment -> connector;

    if (c.name) {
	    segment -> connector = (Connector) TreeSearch (tree, &c);

	    if (!segment -> connector)  {
	        error ("layout uses undefined connector %s", c.name);
            n ++;
        }
	    Deallocate (c.name);
    } 
    else
        segment -> connector = NULL; 

    if (segment -> branch) 
        tree = p -> branch_segment_tree;
    else
        tree = p -> segment_tree;

    segment -> next = TreeSuccessor(tree, segment); 
    segment -> prev = TreePredecessor(tree, segment); 

    if (segment -> branch) {
       if (segment -> next && segment -> next -> branch != segment -> branch)
          segment -> next = NULL;

       if (segment -> prev && segment -> prev -> branch != segment -> branch)
          segment -> prev = NULL;
    }

    return n;
}

/************************************************************************
 * Function:	resolve_branch						*
 *									*
 * Description:	Resolve the names of objects assigned to the 		*
 *		terminal of a branch					*
 ************************************************************************/

static int 
resolve_branch (Item item, void *call_data)
{
    Problem        *p = (Problem *) call_data;
    Branch		    br;
    struct anchor	a;
    struct buoy		b;
    Terminal		t;

    br = (Branch) item;

    t = br -> terminal;

    if (t -> anchor) {
       a.name = (char *) t -> anchor;
       t -> anchor = (Anchor) TreeSearch (p -> anchor_tree, &a);

       if (!t -> anchor)  {
          error ("layout uses undefined anchor %s", a.name);
          return 1;
       }

       Deallocate (a.name);
    }
    else if (t -> buoy) {
       b.name = (char *) t -> buoy;

       t -> buoy = (Buoy) TreeSearch (p -> buoy_tree, &b);

       if (!t -> buoy) {
          error ("layout uses undefined buoy %s", b.name);
          return 1;
       }

       Deallocate (b.name);
    }

    return 0;
}


/************************************************************************
 * Function:	resolve_names						*
 *									*
 * Description:	Resolves the names and numbers of objects making sure	*
 *		all that are specified are actually defined.  We've	*
 *		actually cheated and stored the name or number of the	*
 *		object instead of a pointer to the object in the 	*
 *		structures. 		 				*
 ************************************************************************/

static int
resolve_names (Problem *p)
{
    struct buoy		b;
    struct anchor	a;
    int			i;
    Terminal		t;
    int             n;


	/*
	 * resolve the bouy and anchor names
	 */
   
    n = 0;
    for (i = 1 ; i <= 2 ; i++) {
       t = p -> terminal [i];
       if (!t) {
          error("terminal %d missing", i);
          n++;
          continue;
       }
       if (t -> anchor && t -> buoy) {
          error ("terminal %d has both an anchor (%lx) and a buoy (%lx) defined", i, t -> anchor, t -> buoy);
          n ++;
          continue;
       }
       if (!t -> anchor && !t -> buoy) {
          error ("terminal %d must have an anchor or a buoy defined", i);
          n ++;
          continue;
       }

       if (t -> anchor) {
          a.name = (char *) t -> anchor;
          t -> anchor = (Anchor) TreeSearch (p -> anchor_tree, &a);

	      if (!t -> anchor) {
	          error ("layout uses undefined anchor %s", a.name);
              n ++;
         }
 	     Deallocate (a.name);
      }
      else {
          b.name = (char *) t -> buoy;

          t -> buoy = (Buoy) TreeSearch (p -> buoy_tree, &b);

  	      if (!t -> buoy) {
	         error ("layout uses undefined buoy %s", b.name);
             n ++;
          }

	      Deallocate (b.name);
       }
    }

    if (TreeSize (p -> segment_tree)) {
	    n += TreeSetAndIterate (p -> segment_tree, resolve_segment, p);
    }

    if (TreeSize (p -> branch_segment_tree)) {
	    n += TreeSetAndIterate (p -> branch_segment_tree, resolve_segment, p);
    }

    if (TreeSize (p -> branch_tree)) {
        n += TreeSetAndIterate (p -> branch_tree, resolve_branch, p);
    }

    return n;
}

int
ReadDatabaseFile (char *filename, Tree mat, Tree buo, Tree con, Tree anc)
{
    char    buffer[256];
    char   *plural;
    FILE   *input;
    Analysis a;
    Environment e;
    Problem  p;
    int      n;

    if (!filename || (n = access (filename, R_OK))) {
	    error ("database(1): Unable to open %s (error %d)", filename ? filename : "NULL");
        return 1;
    }

    // cpp = CPP;
/*
    if (cpp) 
        printf("cpp = %s\n", cpp);
    else
        printf("CPP undefined\n");
*/

    if (cpp != NULL) {
#ifdef WINDOWS
        sprintf (buffer, "\"\"%s\" \"%s\"\"", cpp, filename);
#else
        sprintf (buffer, "'%s' < '%s'", cpp, filename);
#endif
        if (!(input = popen (buffer, "r"))) {
	        error ("database(2): Unable to execute %s", cpp);
	        return 1;
        }
    } 
    else if (!(input = fopen (filename, "r"))) {
	    error ("database(3): Unable to open %s", filename);
	    return 1;
	}

    // printf("reading database file %s\n", filename);

    // save pointers to the global trees

    p.num_errors = 0;
    p.line = 1;
    p.filename = filename;
    p.material_tree = mat;
    p.connector_tree = con;
    p.buoy_tree = buo;
    p.anchor_tree = anc;

    ctl.problem = &p;
    ctl.analysis = &a;
    ctl.environment = &e;
    CableInitLexer(input);
    Cable_parse ( );

    if (cpp) {
        pclose (input);
    }
    else { // if (input != stdin) { // impure??
        fclose (input);
    }

    if (p.num_errors) {
        plural = p.num_errors != 1 ? "errors" : "error";
        error ("%u %s found in input", p.num_errors, plural);
        return p.num_errors;
    }

    return 0;
}

/************************************************************************
 * Function:	ReadModelFile						
 *									
 * Description:	Reads a model file using the preprocessor if desired.  	
 *		A filename of "-" indicates standard input (can only be	
 *		used initially) and a NULL filename indicates no file	
 *		(an empty problem is created).				
 ************************************************************************/

int 
ReadModelFile (char *filename, Problem *p, Analysis *a, Environment *e, int parse_only)
{
    char  buffer [2048];
    char *plural;
    FILE *input;
    int   i,j,n;

    // cpp = CPP;

    /* Open the file and send it through the preprocessor. */

    if (filename) {
# ifndef DOS
	    if (cpp != NULL) {
	        if (streq (filename, "-"))
		        sprintf (buffer, "\"%s\"", cpp);
	        else {
		        if (access (filename, R_OK)) {
		            error ("model(1): Unable to open %s", filename);
		            return 1;
		        }
#ifdef WINDOWS
		        sprintf (buffer, "\"%s \"%s\"\"", cpp_command, filename);
#else
		        sprintf (buffer, "%s < '%s'", cpp_command, filename);
#endif
	        }

	        if (!(input = popen (buffer, "r"))) {
		        error ("model(2): Unable to execute %s", cpp);
		        return 1;
	        }

	    } 
        else
# endif
	    {
	        if (streq (filename, "-"))
		        input = NULL; // stdin; // impure??
	        else if (!(input = fopen (filename, "r"))) {
		        error ("model(3): Unable to open %s", filename);
		        return 1;
	        }
	    }
    }
    else
       input = NULL;


    /* Initialize the problem instance. */

    p->input      = NULL;
    p->title	     = strdup ("");
    p->type	     = General;
    p->dynstat	     = 0;
    p->num_errors	     = 0;
    p->line	     = 1;
    p->segment_tree     = TreeCreate (segment_cmp);
    p->branch_tree	     = TreeCreate (branch_cmp);
    p->branch_segment_tree	= TreeCreate (segment_cmp);

    if (!p -> buoy_tree)
        p->buoy_tree	     = TreeCreate (ObjectCompare);
    if (!p -> material_tree)
        p->material_tree    = TreeCreate (ObjectCompare);
    if (!p -> anchor_tree)
        p->anchor_tree      = TreeCreate (ObjectCompare);
    if (!p -> connector_tree)
        p->connector_tree   = TreeCreate (ObjectCompare);
    p->terminal [1]     = NULL;
    p->terminal [2]     = NULL;
    p->junction_size    = 0;
    p->num_branch = 0;
    p->branch = NULL;
    p->node = NULL;
    p->num_nodes = 0;
    p -> segment = NULL;
    p -> num_segments = 0;
    p -> twoD = 0;
    p -> dynamic = 0;

    if (filename) {
	    p->filename = streq (filename, "-") ? "stdin" : filename;
    }
    else
	    p->filename = "";


    /* Initialize the analysis structure. */

    a->maxit		= 0;
    a->outer_it		    = 0;
    a->static_it		= 0;
    a->dynamic_it		= 0;
    a->shooting_it      = 0;

    a->relaxation 	= 0;
    a->dynamic_relaxation = 0.0;
    a->static_relaxation  = 0.0;
    a->outer_relaxation = 0.0;

    a->tolerance		= 0.0;
    a->outer_tolerance = a->static_tolerance = a->dynamic_tolerance = 0.0;

    a->duration		= 0.0;
    a->dt			= 0.0;
    a->ramp_time		= 0.0;

    a->current_steps	= 0;

    a->mesh_smooth	= 0.0;
    a->mesh_amplify	= 0.0;

    a->relax_up		= 1.02;
    a->relax_down		= 1.1;
    a->stall_limit        = 200;

    a->alpha_k		= 0.0;
    a->alpha_m		= 0.0;
    a->gamma		= 0.0;		
    a->sp_rho		= -0.5;		/* typically useful value */

    a->integration	= Spatial;
    a->static_initial_guess = Catenary;
    a->static_solution	  = Relaxation;
    a->static_outer_method = Bisection;

    a->adapt_factor = 10;
    a->adapt_levels = 4;

    a -> viva_dt = 0.0;
    a -> viva_decimate = 0;
    a -> viva_iterations = 0;

    a -> dynamic_var_smooth.value = 1;

    /* Initialize the environment structure. */

    e->input_type		= Regular;
    e->forcing			= Velocity;

    for (i = 1 ; i <= 3 ; i++) {
       e->num_components [i] = 0;

       for (j = 0 ; j < 10 ; j++) {
          e->amplitude [i][j]	= 0.0;
          e->period [i][j]	= 0.0;
          e->phase [i][j]	= 0.0;
          e->omega [i][j]	= 0.0;
          e->k [i][j]		= 0.0;
       }
    }

    e->rho			= 0.0;
    e->surface = e->depth = 0.0;
    e->gravity			= 0.0;
    e->bottom_stiffness	= 0.0;
    e->bottom_damping		= 0.0;
    e->bottom_friction		= 0.0;

    e->fill_density		= 0.0;

    e->bottom_elevation.expr	= NULL;
    e->bottom_elevation.text	= NULL;
    e->bottom_elevation.value	= 0.0;

    e->Uz.expr			= NULL;
    e->Uz.text			= NULL;
    e->Uz.value		= 0.0;
    e->Uy.expr			= NULL;
    e->Uy.text			= NULL;
    e->Uy.value		= 0.0;
    e->Ux.expr			= NULL;
    e->Ux.text			= NULL;
    e->Ux.value		= 0.0;

    e->Uscale			= 1.0;
    e->Uxscale			= 1.0;
    e->Uyscale			= 1.0;
    e->Uzscale			= 1.0;

    e->current_rotation = 0.0;

    e->Uz_mod.expr		= NULL;
    e->Uz_mod.text		= NULL;
    e->Uz_mod.value		= 1.0;
    e->Uy_mod.expr		= NULL;
    e->Uy_mod.text		= NULL;
    e->Uy_mod.value		= 1.0;
    e->Ux_mod.expr		= NULL;
    e->Ux_mod.text		= NULL;
    e->Ux_mod.value		= 1.0;

    e->z_wind.expr		= NULL;
    e->z_wind.text		= NULL;
    e->z_wind.value		= 0.0;
    e->y_wind.expr		= NULL;
    e->y_wind.text		= NULL;
    e->y_wind.value		= 0.0;

    e->wave_file		= NULL;
    e->velocity_file		= NULL;
    e->current_file		= NULL;

    /* Parse the input and resolve the names. */

    ctl.problem    = p;
    ctl.analysis   = a;
    ctl.environment = e;

    if (input) {
	    CableInitLexer(input);
	    Cable_parse ( );
	    p->line = 0;

# ifdef DOS
        fclose (input);
# else
	    if (cpp)
	        pclose (input);
    	else // if (input != stdin)  // impure??
	        fclose (input);
# endif

        if (!parse_only) {
	        if (p->num_errors) {
	            plural = p->num_errors != 1 ? "errors" : "error";
	            error ("%u %s found in input", p->num_errors, plural);
	            return p->num_errors;
	        }

	        if ((n = resolve_names (p))) {
                error("error in name resolution");
                return n;
            }

            p -> branch  = BuildBranchArray(p, &(p -> num_branch));
            CreateNodeArray (p, &(p -> node), &(p -> num_nodes), 
                             &(p -> active), &(p -> num_active));

            p -> segment = BuildSegmentArray(p, &(p -> num_segments), 1);
        }
    }

    return 0;
}

int
SetupCpp(char *cppCmd, char **includes, char **defines)
{
    struct stat     statbuf;

    if (cppCmd) {
        cpp = cppCmd;
    }

    if (cpp && stat(cpp, &statbuf)) { 
        error("program %s does not exist, pre-processing disabled", cpp);
        cpp = NULL;
    }

    if (cpp) {
        snprintf(cpp_command, sizeof(cpp_command), "\"%s\"", cpp);
        while (includes && *includes) {
            strcat(cpp_command, " -I");
            strcat(cpp_command, *includes);
            includes ++;
        }
        while (defines) {
            strcat(cpp_command, " -D");
            strcat(cpp_command, *defines);
            defines ++;
        }
    }

    return 0;
}

/************************************************************************
 * Function:	ParseCppOptions						*
 *									*
 * Description:	Parses and removes the preprocesor options from the	*
 *		command line arguments.					*
 ************************************************************************/

int 
ParseCppOptions (argc, argv)
    int  *argc;
    char *argv [ ];
{
    int   i;
    int   j;
    char *arg;
    char  cpp_args [2048];

    cpp = CPP;

    j = 1;
    cpp_args [0] = 0;

    for (i = 1; i < *argc; i ++)
	if ((arg = argv [i]) [0] != '-') {
	    argv [j ++] = arg;
	} else if (streq (arg, "-nocpp")) {
	    cpp = NULL;
	} else if (streq (arg, "-cpp")) {
	    if (++ i == *argc)
		return 1;
	    cpp = argv [i];
	} else if (arg [1] == 'D' || arg [1] == 'U' || arg [1] == 'I') {
	    strcat (cpp_args, " '");
	    strcat (cpp_args, arg);
	    strcat (cpp_args, "'");
	} else
	    argv [j ++] = arg;

#ifdef LIBDIR
    if (cpp != NULL)
	    sprintf (cpp_command, "\"%s\" -I%s %s", cpp, (char *) LIBDIR, cpp_args);
#else
    if (cpp != NULL)
	    sprintf (cpp_command, "\"%s\" %s", cpp, cpp_args);
#endif

    argv [*argc = j] = NULL;
    return 0;
}

void
FreeNode(Node n)
{
    ZeroOffset(n -> Ys); Deallocate(n -> Ys);
    ZeroOffset(n -> Y); Deallocate(n -> Y);
    ZeroOffset(n -> Y_o); Deallocate(n -> Y_o);
    ZeroOffset(n -> Yd_o); Deallocate(n -> Yd_o);
    ZeroOffset(n -> Yd); Deallocate(n -> Yd);
    ZeroOffset(n -> Y_f); Deallocate(n -> Y_f);
    ZeroOffset(n -> Y_o_f); Deallocate(n -> Y_o_f);

    Deallocate(n);  
}

void
FreeProblem(Problem *p)
{
    int i;

    if (!p)
        return;

    Deallocate(p -> title);
    Deallocate(p -> input);
    // Deallocate(p -> filename); // not allocated anymore
    
    DestroyTerminal(p -> terminal[1]);
    DestroyTerminal(p -> terminal[2]);

    TreeDestroy(p -> segment_tree);
    TreeDestroy(p -> branch_tree);
    TreeDestroy(p -> branch_segment_tree);
    TreeDestroy(p -> buoy_tree);
    TreeDestroy(p -> material_tree);
    TreeDestroy(p -> connector_tree);
    TreeDestroy(p -> anchor_tree);
    // TreeDestroy(p -> layout_tree);
    
    for (i = 1 ; i <= p -> num_branch ; i++) {
        DestroyBranch(p -> branch[i]);
    }
    ZeroOffset(p -> branch); Deallocate(p -> branch);

    for (i = 1 ; i <= p -> num_nodes ; i++) {
        FreeNode(p -> node[i]); 
    }

    ZeroOffset(p -> node); Deallocate(p -> node);

    Deallocate(p);
}

void
FreeAnalysis(Analysis *a)
{
    Deallocate(a);
}

void
FreeEnvironment(Environment *e)
{
    if (!e) 
        return;

    Deallocate(e -> velocity_file);
    Deallocate(e -> wave_file);
    Deallocate(e -> current_file);

    FreeCode(e -> bottom_elevation.expr);
    Deallocate(e -> bottom_elevation.text);
    FreeCode(e -> z_wind.expr);
    Deallocate(e -> z_wind.text);
    FreeCode(e -> y_wind.expr);
    Deallocate(e -> y_wind.text);
    FreeCode(e -> Ux.expr);
    Deallocate(e -> Ux.text);
    FreeCode(e -> Uy.expr);
    Deallocate(e -> Uy.text);
    FreeCode(e -> Uz.expr);
    Deallocate(e -> Uz.text);
    FreeCode(e -> Uz_mod.expr);
    Deallocate(e -> Uz_mod.text);
    FreeCode(e -> Uy_mod.expr);
    Deallocate(e -> Uy_mod.text);
    FreeCode(e -> Ux_mod.expr);
    Deallocate(e -> Ux_mod.text);
   
    Deallocate(e); 
}

void
FreeSolution(Solution *s)
{
    FreeProblem(s -> problem);
    FreeEnvironment(s -> environment);
    FreeAnalysis(s -> analysis);

    Deallocate(s -> in_name);
    Deallocate(s -> out_name);
    Deallocate(s -> static_file);
    Deallocate(s -> initial_file);
    Deallocate(s -> dynstat_name);

    Deallocate(s -> output_nodes);
    Deallocate(s -> tmpdir);

    // s -> results
    // s -> playback


    Deallocate(s);
}

