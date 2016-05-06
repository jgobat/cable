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
 * File:	file.c							*
 *									*
 * Description:	This file contains the public and private function	*
 *		definitions for writing model files and objects.	*
 ************************************************************************/

# include <stdio.h>
# include <ctype.h>
# include <string.h>
# include "compress.h"
# include "problem.h"
# include "error.h"

/* Nasty macros for printing things just oh so right. */

# define cfprintf(fp,f,x) \
	if (x) fprintf (fp, f, x)

# define PrintHeader(fp) \
	if (!printed_header ++) fprintf (fp, "\ncanvas configuration\n"); \
	if (last_section) {fprintf (fp, "\n"); last_section = 0;} \
	this_section = 1

# define PrintBoolean(fp,n,x) \
	if (x != UnspecifiedValue) {\
	    PrintHeader(fp); \
	    fprintf (fp, "%s%s=", this_line ++ ? " " : "", n); \
	    fprintf (fp, "%s", x ? "true" : "false"); \
	}

# define PrintFloat(fp,n,x) \
	if (x != UnspecifiedValue) {\
	    PrintHeader(fp); \
	    fprintf (fp, "%s%s=%g", this_line ++ ? " " : "", n, x); \
	}

# define PrintInteger(fp,n,x) \
	if (x != UnspecifiedValue) {\
	    PrintHeader(fp); \
	    fprintf (fp, "%s%s=%d", this_line ++ ? " " : "", n, x); \
	}

# define PrintString(fp,n,x) \
	if (x != NULL) {\
	    PrintHeader(fp); \
	    fprintf (fp, "%s%s=%s", this_line ++ ? " " : "", n, Quote (x)); \
	}

# define InitFormat(fp) \
	printed_header = this_line = last_section = 0

# define StartSection(fp) \
	if (!last_section) last_section = this_section; \
	this_section = 0

# define EndSection(fp) \

# define PrintNewline(fp) \
	if (this_line) {fprintf (fp, "\n"); this_line = 0;}


typedef struct {
    FILE *fp;
    int   only_show_used;
} output;


/************************************************************************
 * Function:	Quote							*
 *									*
 * Description:	Quotes a string if necessary.				*
 ************************************************************************/

static char *Quote (s)
    char *s;
{
    char	c;
    char       *ptr;
    static char buffer [256];

    
    if (s == NULL || strcmp (s, "") == 0)
       return "\"\"";

    for (ptr = s; (c = *ptr); ptr ++)
	if (!isalpha (c) && c != '_' && (isdigit (c) ? ptr == s : 1)) {
	    sprintf (buffer, "\"%s\"", s);
	    return buffer;
	}

    return s;
}

static void 
WriteTerminal (Terminal t, output *out)
{

    fprintf (out->fp, "    terminal = {\n");
    if (t -> anchor) {
        fprintf (out->fp, "        anchor = %s\n", t -> anchor -> name);
        t -> anchor -> aux = (char * ) 1;
    }

    cfprintf (out->fp, "        z = %g\n", t -> x); 
    cfprintf (out->fp, "        x = %g\n", t -> y); 
    cfprintf (out->fp, "        y = %g\n", t -> z); 

    else if (t -> buoy) {
        fprintf (out->fp, "        buoy = %s\n", t -> buoy -> name);
        t -> buoy -> aux = (char *) 1;
    }

    if (t -> release >= 0.0)
       fprintf (out->fp, "        release-time = %g\n", t -> release);

    cfprintf (out->fp, "        z-force = %g\n", t -> xforce);
    cfprintf (out->fp, "        x-force = %g\n", t -> yforce);
    cfprintf (out->fp, "        y-force = %g\n", t -> zforce);

    if (t -> xthrust.expr)
       fprintf (out->fp, "        z-thrust = %s\n", 
                t -> xthrust.text);
    else 
       cfprintf (out->fp, "        z-thrust = %g\n", 
                 t -> xthrust.value);
    if (t -> ythrust.expr)
       fprintf (out->fp, "        x-thrust = %s\n", 
                t -> ythrust.text);
    else 
       cfprintf (out->fp, "        x-thrust = %g\n", 
                 t -> ythrust.value);
    if (t -> zthrust.expr)
       fprintf (out->fp, "        y-thrust = %s\n", 
                t -> zthrust.text);
    else 
       cfprintf (out->fp, "        y-thrust = %g\n", 
                 t -> zthrust.value);

    if (t -> xspeed.expr)
       fprintf (out->fp, "        z-speed = %s\n", 
                t -> xspeed.text);
    else 
       cfprintf (out->fp, "        z-speed = %g\n", 
                 t -> xspeed.value);
    if (t -> yspeed.expr)
       fprintf (out->fp, "        x-speed = %s\n", 
                t -> yspeed.text);
    else 
       cfprintf (out->fp, "        x-speed = %g\n", 
                 t -> yspeed.value);
    if (t -> zspeed.expr)
       fprintf (out->fp, "        y-speed = %s\n", 
                t -> zspeed.text);
    else 
       cfprintf (out->fp, "        y-speed = %g\n", 
                 t -> zspeed.value);

    fprintf (out->fp, "    }\n");
}


int CountNodes (Segment s)
{
   int		i;
   int		nn;

   nn = 0;

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

   return nn;
}

/************************************************************************
 * Function:	WriteSegment						*
 *									*
 * Description:	Writes an segment to the specified stream.		*
 ************************************************************************/

static int nodes = 0;
static double main_s = 0.0;
static int segments = 0;

static int 
WriteSegment (Item item, void *call_data)
{
    output *out = (output *) call_data;
    Segment	segment;
    int		i, j;
    Branch	br;
   
    segment = (Segment) item;

    main_s += segment -> length;
    nodes += CountNodes(segment);

    segments ++;

    fprintf (out->fp, "    segment = { /* %d */\n", segments);
    fprintf (out->fp, "        length   = %g\n", segment -> length);
    fprintf (out->fp, "        material = %s\n", segment -> material -> name);

    if (segment -> num_dist) {
       fprintf (out->fp, "        nodes    ="); 
       for (i = 1 ; i <= segment -> num_dist ; i++)
          fprintf (out->fp, " (%d, %g)", segment -> dist [i].nodes,
                                    segment -> dist [i].percent);

       fprintf (out->fp, "\n");
    }

    if (segment -> num_attach) {
        fprintf (out->fp, "        attachments = ");
        for (i = 1 ; i <= segment -> num_attach ; i++) {
           if (i > 1)
              fprintf (out->fp, "                      ");

           fprintf (out->fp, "%s : (", segment -> attach [i].object -> name);
           segment -> attach[1].object -> aux = (char *) 1;
           for (j = 1 ; j <= segment -> attach [i].num_nodes ; j++) {
              fprintf (out->fp,"%d", segment -> attach [i].nodes [j]);
              if (j < segment -> attach [i].num_nodes)
                 fprintf (out->fp, ", ");
           }
           fprintf (out->fp, ")\n");
        }
    }
    if (segment -> num_bottom_dist) {
        fprintf (out->fp, "        bottom-spool-nodes = ");
        for (i = 1 ; i <= segment -> num_bottom_dist ; i++)   
            fprintf(out->fp, "(%d, %g)", segment -> bottom_dist[i].nodes,
                                         segment -> bottom_dist[i].percent);
        fprintf(out->fp, "\n");
        cfprintf (out->fp, "        bottom-spool-length = %g\n", segment -> bottom_length);
        if (segment -> bottom_pay.expr) 
            fprintf (out->fp, "        bottom-spool-pay = %s\n", segment -> bottom_pay.text);
        else if (segment -> bottom_pay.value) 
            cfprintf (out->fp, "        bottom-spool-pay = %f\n", segment -> bottom_pay.value);
    }
    if (segment -> num_top_dist) {
        fprintf (out->fp, "        top-spool-nodes = ");
        for (i = 1 ; i <= segment -> num_top_dist ; i++)   
            fprintf(out->fp, "(%d, %g)", segment -> top_dist[i].nodes,
                                         segment -> top_dist[i].percent);
        fprintf(out->fp, "\n");
        cfprintf (out->fp, "        top-spool-length = %g\n", segment -> top_length);
        if (segment -> top_pay.expr) 
            fprintf (out->fp, "        top-spool-pay = %s\n", segment -> top_pay.text);
        else if (segment -> top_pay.value) 
            cfprintf (out->fp, "        top-spool-pay = %f\n", segment -> top_pay.value);
    }
    fprintf (out->fp, "    } /* %d, %g */\n", nodes, main_s);
 
    if (segment -> connector) {
       if (segment -> connector_x
           || segment -> connector_y
           || segment -> connector_z
           || segment -> connector_xforce
           || segment -> connector_yforce
           || segment -> connector_zforce
           || segment -> connector_xthrust.expr
           || segment -> connector_ythrust.expr
           || segment -> connector_zthrust.expr
           || segment -> connector_xthrust.value
           || segment -> connector_ythrust.value
           || segment -> connector_zthrust.value) {

          fprintf(out -> fp, "    connector = {\n");
          fprintf(out -> fp, "        body = %s\n", segment -> connector -> name);
	      cfprintf (out->fp, "        x = %g\n", segment -> connector_y);
	      cfprintf (out->fp, "        y = %g\n", segment -> connector_z);
	      cfprintf (out->fp, "        z = %g\n", segment -> connector_x);
	      cfprintf (out->fp, "        x-force = %g\n", segment -> connector_yforce);
	      cfprintf (out->fp, "        y-force = %g\n", segment -> connector_zforce);
	      cfprintf (out->fp, "        z-force = %g\n", segment -> connector_xforce);
          if (segment -> connector_xthrust.expr)
	        fprintf (out->fp, "        z-thrust = %s\n", segment -> connector_xthrust.text);
        else
	        cfprintf (out->fp, "        z-thrust = %g\n", segment -> connector_xthrust.value);
          if (segment -> connector_ythrust.expr)
	        fprintf (out->fp, "        x-thrust = %s\n", segment -> connector_ythrust.text);
        else
	        cfprintf (out->fp, "        x-thrust = %g\n", segment -> connector_ythrust.value);
          if (segment -> connector_zthrust.expr)
	        fprintf (out->fp, "        y-thrust = %s\n", segment -> connector_zthrust.text);
        else
	        cfprintf (out->fp, "        y-thrust = %g\n", segment -> connector_zthrust.value);

           fprintf(out -> fp, "    }\n");
       }
       else {
          fprintf (out->fp, "    connector = %s\n", segment -> connector -> name);
       }
       segment -> connector -> aux = (char *) 1;
    }

    if (segment -> num_branch_to) {
       for (i = 1 ; i <= segment -> num_branch_to ; i++) {
          br = segment -> branch_to [i];

          fprintf (out->fp, "    branch = { \n");

          for (j = 1 ; j <= br -> num_segment ; j++)
             WriteSegment((Item) br -> segment [j], out);

          WriteTerminal(br -> terminal, out);
          fprintf (out->fp, "    }\n");
       }
    }

    segment -> material -> aux = (char *) 1;

    return 0;
}


/************************************************************************
 * Function:	WriteMaterial						*
 *									*
 * Description:	Writes a material to the specified stream if the aux	*
 *		pointer matches the mark flag.  The pointer is then	*
 *		cleared.						*
 ************************************************************************/

static int 
WriteMaterial (Item item, void *call_data)
{
    output *out = (output *) call_data;
    Material material;


    material = (Material) item;

    if (material -> aux || out -> only_show_used == 0) {
	    fprintf (out->fp, "  %s\n", Quote (material -> name));

	    cfprintf (out->fp, "    EA    = %g\n", material -> EA);
	    cfprintf (out->fp, "    EI    = %g\n", material -> EI);
	    cfprintf (out->fp, "    GJ    = %g\n", material -> GJ);
	    cfprintf (out->fp, "    d     = %g\n", material -> d);
	    cfprintf (out->fp, "    A     = %g\n", material -> A);
	    cfprintf (out->fp, "    m     = %g\n", material -> m);
	    cfprintf (out->fp, "    w     = %g\n", material -> w);
	    cfprintf (out->fp, "    wet   = %g\n", material -> wet);
	    cfprintf (out->fp, "    amn   = %g\n", material -> amn);
	    cfprintf (out->fp, "    amt   = %g\n", material -> amt);
	    cfprintf (out->fp, "    Cmn   = %g\n", material -> Cmn);
	    cfprintf (out->fp, "    Cmt   = %g\n", material -> Cmt);
	    cfprintf (out->fp, "    Can   = %g\n", material -> Can);
	    cfprintf (out->fp, "    Cat   = %g\n", material -> Cat);
        cfprintf (out->fp, "    bt    = %g\n", material -> bt);
        cfprintf (out->fp, "    bn    = %g\n", material -> bn);

        if (material -> Cdn.expr)
	        fprintf (out->fp, "    Cdn   = %s\n", material -> Cdn.text);
        else
	        fprintf (out->fp, "    Cdn   = %g\n", material -> Cdn.value);

        if (material -> Cdt.expr)
	        fprintf (out->fp, "    Cdt   = %s\n", material -> Cdt.text);
        else
	        fprintf (out->fp, "    Cdt   = %g\n", material -> Cdt.value);

        if (material -> color)
            fprintf (out->fp, "    color = %s\n", material -> color);

        if (material -> type == Nonlinear) {
           if (material -> T [0].expr)
              fprintf (out->fp, "    T = %s\n", material -> T [0].text);
           if (material -> T [1].expr)
              fprintf (out->fp, "    Te = %s\n", material -> T [1].text);
           if (material -> T [2].expr)
              fprintf (out->fp, "    Tee = %s\n", material -> T [2].text);
       }
    }

    material -> aux = NULL;
    return 0;
}


/************************************************************************
 * Function:	WriteBuoy						*
 *									*
 * Description:	Writes a buoy to the specified stream if the aux	*
 *		pointer matches the mark flag.  The pointer is then	*
 *		cleared.						*
 ************************************************************************/

static int 
WriteBuoy (Item item, void *call_data)
{
    output *out = (output *) call_data;
    int	  i;
    Buoy  buoy;
    static char *names [] = {"axisymmetric", "sphere", "cylinder", "capsule", "ship", "platform"};


    buoy = (Buoy) item;

    if (buoy -> aux || out -> only_show_used == 0) {
	    fprintf (out->fp, "  %s\n", Quote (buoy -> name));

        if (buoy -> type)
           fprintf (out->fp, "    type     = %s\n", names [buoy -> type]);

        cfprintf (out->fp, "    d        = %g\n", buoy -> d);
        cfprintf (out->fp, "    h        = %g\n", buoy -> h);
        cfprintf (out->fp, "    m        = %g\n", buoy -> m);
        cfprintf (out->fp, "    cg       = %g\n", buoy -> cg);
        cfprintf (out->fp, "    buoyancy = %g\n", buoy -> buoyancy);
        cfprintf (out->fp, "    Cdt      = %g\n",  buoy -> Cdt);
        cfprintf (out->fp, "    Cdn      = %g\n",  buoy -> Cdn);
        cfprintf (out->fp, "    Cdw      = %g\n",  buoy -> Cdw);
        cfprintf (out->fp, "    Sw       = %g\n",  buoy -> Sw);

        if (buoy -> diameters) {
           fprintf (out->fp, "    diameters = ");

           for (i = 1 ; i <= buoy -> num_diameters ; i++)
              fprintf (out->fp, "(%g, %g) ", buoy -> diameters [i].level,
                                        buoy -> diameters [i].d);

           fprintf (out->fp,"\n");
        }

        if (buoy -> color)
            fprintf (out->fp, "    color    = %s\n", buoy -> color);
    }

    buoy -> aux = NULL;
    return 0;
}

/************************************************************************
 * Function:	WriteAnchor						*
 *									*
 * Description:	Writes an anchor to the specified stream if the aux	*
 *		pointer matches the mark flag.  The pointer is then	*
 *		cleared.						*
 ************************************************************************/

static int 
WriteAnchor (Item item, void *call_data)
{
    output *out = (output *) call_data;
    Anchor anchor;


    anchor = (Anchor) item;

    if (anchor -> aux || out->only_show_used == 0) {
	    fprintf (out->fp, "  %s\n", Quote (anchor -> name));

        if (anchor -> color)
            fprintf (out->fp, "    color = %s\n", anchor -> color);
    }

    anchor -> aux = NULL;
    return 0;
}

/************************************************************************
 * Function:	WriteConnector						*
 *									*
 * Description:	Writes a connector to the specified stream if the aux	*
 *		pointer matches the mark flag.  The pointer is then	*
 *		cleared.						*
 ************************************************************************/

static int 
WriteConnector (Item item, void *call_data)
{
    output *out = (output *) call_data;
    Connector connector;


    connector = (Connector) item;

    if (connector -> aux || out -> only_show_used == 0) {
	    fprintf (out->fp, "  %s\n", Quote (connector -> name));

        cfprintf (out->fp, "    m        = %g\n", connector -> m);
        cfprintf (out->fp, "    d        = %g\n", connector -> d);
        cfprintf (out->fp, "    wet      = %g\n", connector -> wet);
        cfprintf (out->fp, "    Cdt      = %g\n", connector -> Cdt);
        cfprintf (out->fp, "    Cdn      = %g\n", connector -> Cdn);
       
        if (connector -> color)
            fprintf (out->fp, "    color    = %s\n", connector -> color);
    }

    connector -> aux = NULL;
    return 0;
}


/************************************************************************
 * Function:	WriteAnalysisParameters					*
 *									*
 * Description:	Writes out the analysis parameters section.		*
 ************************************************************************/

static void 
WriteAnalysisParameters (Analysis *analysis, output *out)
{
    static char *integr [] = {"", "spatial", "temporal"};
    static char *solution [] = {"", "relaxation", "catenary", "shooting"};

    cfprintf (out->fp, "    max-iterations      = %d\n", analysis->maxit);
    cfprintf (out->fp, "    relaxation          = %g\n", analysis->relaxation);
    cfprintf (out->fp, "    tolerance           = %g\n", analysis->tolerance);

    fprintf (out->fp,  "    static-solution      = %s\n", 
             solution [analysis->static_solution]);
    fprintf (out->fp,  "    static-initial-guess = %s\n", 
             solution [analysis->static_initial_guess]);

    cfprintf (out->fp, "    static-iterations    = %d\n", analysis->static_it);
    cfprintf (out->fp, "    static-relaxation    = %g\n", analysis->static_relaxation);
    cfprintf (out->fp, "    static-tolerance     = %g\n", analysis->static_tolerance);
    cfprintf (out->fp, "    current-steps        = %d\n", analysis->current_steps);
    cfprintf (out->fp, "    relax-adapt-up       = %g\n", analysis->relax_up);
    cfprintf (out->fp, "    relax-adapt-down     = %g\n", analysis->relax_down);
    cfprintf (out->fp, "    mesh-amplification   = %g\n", analysis->mesh_amplify);

    cfprintf (out->fp, "    static-outer-iterations = %d\n", analysis->outer_it);
    cfprintf (out->fp, "    static-outer-relaxation = %g\n", analysis->outer_relaxation);
    cfprintf (out->fp, "    static-outer-tolerance  = %g\n", analysis->outer_tolerance);

    cfprintf (out->fp, "    dynamic-iterations  = %d\n", analysis->dynamic_it);
    cfprintf (out->fp, "    dynamic-relaxation  = %g\n", analysis->dynamic_relaxation);
    cfprintf (out->fp, "    dynamic-tolerance   = %g\n", analysis->dynamic_tolerance);
    cfprintf (out->fp, "    time-step           = %g\n", analysis->dt);
    cfprintf (out->fp, "    duration            = %g\n", analysis->duration);
    cfprintf (out->fp, "    ramp-time           = %g\n", analysis->ramp_time);
    if (analysis->duration) {
       fprintf (out->fp, "    dynamic-integration = %s\n", integr [analysis->integration]);
       cfprintf (out->fp, "    dynamic-rho         = %g\n", analysis->sp_rho);
       cfprintf (out->fp, "    dynamic-alpha-k     = %g\n", analysis->alpha_k);
       cfprintf (out->fp, "    dynamic-alpha-m     = %g\n", analysis->alpha_m);
       cfprintf (out->fp, "    dynamic-gamma       = %g\n", analysis->gamma);
    }

}

static void 
WriteEnvironment (Environment *environment, output *out)
{
    static char *methods [] = {"", "velocity", "wave-follower", "force", "spar", "froude-krylov", "morison"};
    static char *types [] = {"", "regular", "random"};
    static char *input_names [] = {"", "z-input", "x-input", "y-input"};
    static char *wave_names [] = {"", NULL, "x-wave ", "y-wave "};
    char **names;
    int    i, j;

    fprintf  (out->fp, "    forcing-method     = %s\n", methods [environment->forcing]);
    fprintf  (out->fp, "    input-type         = %s\n", types [environment->input_type]);
    cfprintf (out->fp, "    rho                = %g\n", environment->rho);
    cfprintf (out->fp, "    gravity            = %g\n", environment->gravity);
    cfprintf (out->fp, "    depth              = %g\n", environment->depth);
    cfprintf (out->fp, "    bottom-stiffness   = %g\n", environment->bottom_stiffness);
    cfprintf (out->fp, "    bottom-damping     = %g\n", environment->bottom_damping);
    cfprintf (out->fp, "    bottom-friction    = %g\n", environment->bottom_friction);
    if (environment->bottom_elevation.expr) 
       fprintf (out->fp, "    bottom-elevation    = %s\n", 
		environment->bottom_elevation.text);
    else 
       cfprintf (out->fp, "    bottom-elevation    = %g\n", 
		 environment->bottom_elevation.value);
   
    if (environment->forcing == Velocity || environment->forcing == Force)
       names = input_names;
    else
       names = wave_names;
 
    for (j = 1 ; j <= 3 ; j ++) {
       if (!names [j])
          continue;

       if (environment->num_components [j]) { 
          fprintf (out->fp, "    %s            = ", names [j]);

          for (i = 0 ; i < environment->num_components [j] ; i++) {
             fprintf (out->fp, "(%g,%g,%g) ", 
                      environment->amplitude [j][i], 
                      environment->period [j][i], 
                      environment->phase [j][i]);
          }
          fprintf (out->fp, "\n");
       }
    }

    fprintf(out->fp, "    current-scale      = %f\n", environment->Uscale);
    cfprintf(out->fp, "    current-rotation  = %f\n", environment->current_rotation);
    cfprintf(out->fp, "    x-current-scale    = %f\n", environment->Uyscale);
    cfprintf(out->fp, "    y-current-scale    = %f\n", environment->Uzscale);
    cfprintf(out->fp, "    z-current-scale    = %f\n", environment->Uxscale);
    if (environment -> current_file) {
       fprintf (out->fp, "    current-file = %s\n", environment -> current_file);
    }
    if (environment->Uz.expr)
       fprintf (out->fp, "    y-current           = %s\n", environment->Uz.text);
    else
       cfprintf (out->fp, "    y-current           = %g\n", environment->Uz.value);

    if (environment->Uy.expr)
       fprintf (out->fp, "    x-current           = %s\n", environment->Uy.text);
    else
       cfprintf (out->fp, "    x-current           = %g\n", environment->Uy.value);

    if (environment->Ux.expr)
       fprintf (out->fp, "    z-current           = %s\n", environment->Ux.text);
    else
       cfprintf (out->fp, "    z-current           = %g\n", environment->Ux.value);

    if (environment->Uz_mod.expr)
       fprintf (out->fp, "    y-current-modulation = %s\n", environment->Uz_mod.text);
    else if (environment->Uz_mod.value != 1)
       fprintf (out->fp, "    y-current-modulation = %g\n", environment->Uz_mod.value);

    if (environment->Uy_mod.expr)
       fprintf (out->fp, "    x-current-modulation = %s\n", environment->Uy_mod.text);
    else if (environment->Uy_mod.value != 1)
       fprintf (out->fp, "    x-current-modulation = %g\n", environment->Uy_mod.value);

    if (environment->Ux_mod.expr)
       fprintf (out->fp, "    z-current-modulation = %s\n", environment->Ux_mod.text);
    else if (environment->Ux_mod.value != 1)
       fprintf (out->fp, "    z-current-modulation = %g\n", environment->Ux_mod.value);

    if (environment->z_wind.expr)
       fprintf (out->fp, "    y-wind              = %s\n", environment->z_wind.text);
    else
       cfprintf (out->fp, "    y-wind              = %g\n", environment->z_wind.value);

    if (environment->y_wind.expr)
       fprintf (out->fp, "    x-wind              = %s\n", environment->y_wind.text);
    else
       cfprintf (out->fp, "    x-wind              = %g\n", environment->y_wind.value);
}

char *
ProblemTypeName(int problem_type)
{
    static char	*type_names [] = {"", "general", "surface", "subsurface", "horizontal", "towing", "drifter", "deployment", "riser", "webster", "horizontal-drifter"};

    return type_names[problem_type];
}

/************************************************************************
 * Function:	WriteFile						*
 *									*
 * Description:	Writes a buoy file.  A filename of "-" indicates	*
 *		standard output.					*
 ************************************************************************/

static int 
WriteFile (output *out, Problem *problem, Analysis *analysis, Environment *environment)
{

    /* Write the problem description section. */

    fprintf (out -> fp, "problem description\n");
    if (problem->title != NULL && strcmp (problem->title, ""))   
       fprintf (out -> fp, "    title = %s\n", Quote (problem->title));

    fprintf (out -> fp,"     type = %s\n", ProblemTypeName(problem->type));


    /* Write the analysis parameters section */

    fprintf (out -> fp, "\nanalysis parameters\n");
    WriteAnalysisParameters (analysis, out);	

    fprintf (out -> fp, "\nenvironment\n");
    WriteEnvironment (environment, out);


    /* Write the layout section marking referenced objects. */

    fprintf (out -> fp, "\nlayout\n");

    WriteTerminal (problem->terminal [1], out);


    nodes = 0;
    main_s = 0.0;
    TreeSetIterator (problem->segment_tree, WriteSegment, (void *) out);
    TreeIterate (problem->segment_tree);


    WriteTerminal (problem->terminal [2], out);

    /* Write the materials section. */

    if (TreeSize (problem->material_tree) > 0) {
	    fprintf (out -> fp, "\nmaterials\n");
	    TreeSetIterator (problem->material_tree, WriteMaterial, (void *) out);
	    TreeIterate (problem->material_tree);
    }


    /* Write the buoys section. */

    if (TreeSize (problem->buoy_tree) > 0) {
	    fprintf (out -> fp, "\nbuoys\n");
	    TreeSetIterator (problem->buoy_tree, WriteBuoy, (void *) out);
	    TreeIterate (problem->buoy_tree);
    }


    /* Write the anchors section. */

    if (TreeSize (problem->anchor_tree) > 0) {
	    fprintf (out -> fp, "\nanchors\n");
	    TreeSetIterator (problem->anchor_tree, WriteAnchor, (void *) out);
	    TreeIterate (problem->anchor_tree);
    }


    /* Write the connectors section. */

    if (TreeSize (problem->connector_tree) > 0) {
	    fprintf (out -> fp, "\nconnectors\n");
	    TreeSetIterator (problem->connector_tree, WriteConnector, (void *) out);
	    TreeIterate (problem->connector_tree);
    }

    fprintf (out -> fp, "\nend\n");

    return 0;
}

static FILE *
OpenFile (char *filename)
{
    FILE *fp;

    if (strcmp (filename, "-")) {
	    if (!(fp = fopen (filename, "w"))) {
	        error ("Unable to open %s", filename);
	        return NULL;
	    }
    } 
    else
#ifdef WINDOWS
	    fp = NULL;
#else
	    fp = stdout;
#endif

    return fp;
}

static void 
CloseFile (FILE *fp)
{
#ifndef WINDOWS
    if (fp != stdout)
	    fclose (fp);
#endif

    return;
}

/************************************************************************
 * Function:	WriteModelFile						*
 *									*
 * Description:	Writes a model file -- only referenced objects will be	*
 *		written.						*
 ************************************************************************/

int 
WriteModelFile (char *filename, Problem *problem, Analysis *analysis, Environment *environment)
{
    output  out;

    if ((out.fp = OpenFile (filename)) == NULL)
       return 1;

    out.only_show_used = 1;

    if (WriteFile(&out, problem, analysis, environment))
       return 1;

    CloseFile (out.fp);

    return 0;
}


/************************************************************************
 * Function:	DumpModelFile						
 *	
 * Description:	Dumps a Model file -- referenced and unreferenced	
 *		objects will be written.				
 ************************************************************************/

int 
DumpModelFile (char *filename, Problem *problem, Analysis *analysis, Environment *environment)
{
    output out;

    if ((out.fp = OpenFile(filename)) == NULL)
       return 1;

    out.only_show_used = 0;

    if (WriteFile(&out, problem, analysis, environment))
       return 1;

    CloseFile (out.fp);

    return 0;
}

/************************************************************************
 * Function:	fWriteModelFile						*
 *									*
 * Description:	Writes a model file -- only referenced objects will be	*
 *		written.						*
 ************************************************************************/

int 
fWriteModelFile (FILE *stream, Problem *problem, Analysis *analysis, Environment *environment)
{
    output out;

    out.fp = stream;
    out.only_show_used = 1;

    if (WriteFile(&out, problem, analysis, environment))
       return 1;

    return 0;
}


/************************************************************************
 * Function:	fDumpModelFile						*
 *									*
 * Description:	Dumps a model file -- referenced and unreferenced	*
 *		objects will be written.				*
 ************************************************************************/

int
fDumpModelFile (FILE *stream, Problem *problem, Analysis *analysis, Environment *environment)
{
    output out;

    out.fp = stream;
    out.only_show_used = 0;

    if (WriteFile(&out, problem, analysis, environment))
       return 1;

    return 0;
}
