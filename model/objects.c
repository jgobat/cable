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
 * File:	objects.c						*
 *									*
 * Description:	This file contains the function definitions for		*
 *		operations on the various objects.			*
 ************************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "compress.h"
# include "problem.h"
# include "error.h"
# include "objects.h"
# include "allocate.h"


/************************************************************************
 * Function:	 CreateBranch						*
 ************************************************************************/

Branch CreateBranch (number)
    unsigned   number;
{
    Branch	branch;
    int		i;

    if (!(branch = AllocNew (struct branch)))
	Fatal ("unable to allocate memory for new branch");

    branch -> aux = NULL;
    branch -> number = number;

    branch -> first = 0;
    branch -> last  = 0;
    branch -> num_segment = 0;

    branch -> terminal     = NULL;
    branch -> segment_from = NULL;

    branch -> w0 = 0.0;
    branch -> length = 0.0;

    for (i = 0 ; i < 10 ; i++)
       branch -> segment [i] = NULL;
   
    return branch;
}


/************************************************************************
 * Function:	 DestroyBranch						*
 ************************************************************************/

void DestroyBranch (branch)
    Branch	branch;
{
    if (branch) {
	    Deallocate (branch -> aux);
	    Deallocate (branch);
    }

    return; 
}


/************************************************************************
 * Function:	 CreateSegment						*
 ************************************************************************/

Segment CreateSegment (number)
    unsigned   number;
{
    Segment	segment;
    int		i;

    if (!(segment = AllocNew (struct segment)))
	    Fatal ("unable to allocate memory for new segment");

    segment -> aux = NULL;
    segment -> name = NULL;
    segment -> number = number;
    segment -> num_nodes = 0;

    segment -> num_dist = 0;
    segment -> num_attach = 0;
    segment -> attach = NULL;

    segment -> junction.num_nodes = 0;

    segment -> material  = NULL;
    segment -> connector = NULL;
    segment -> connection = Spliced;

    segment -> length   = 0.0;

    segment -> branch = NULL;

    for (i = 0 ; i < 10 ; i++)
       segment -> branch_to [i] = NULL;

    segment -> num_branch_to = 0;

    segment -> top_length = 0.0;
    segment -> bottom_length = 0.0;

    segment -> top_wet = 0.0;
    segment -> bottom_wet = 0.0;

    segment -> top_added = 0.0;
    segment -> bottom_added = 0.0;
  
    segment -> num_top_dist = 0; 
    segment -> num_bottom_dist = 0; 

    segment -> top_pay.expr = NULL;
    segment -> top_pay.text = NULL;
    segment -> top_pay.value = 0.0;

    segment -> bottom_pay.expr = NULL;
    segment -> bottom_pay.text = NULL;
    segment -> bottom_pay.value = 0.0;

    segment -> next = NULL;
    segment -> prev = NULL;

    segment -> top_pay_speed = segment -> top_pay_speed_o = 0;
    segment -> bottom_pay_speed = segment -> bottom_pay_speed_o = 0;

    segment -> prev_t = segment -> prev_dt = 0;

    segment -> connector_x = 0;
    segment -> connector_y = 0;
    segment -> connector_z = 0;

    segment -> connector_xforce = 0;
    segment -> connector_yforce = 0;
    segment -> connector_zforce = 0;

    segment -> connector_xthrust.expr = NULL;
    segment -> connector_xthrust.text = NULL;
    segment -> connector_xthrust.value = 0;

    segment -> connector_ythrust.expr = NULL;
    segment -> connector_ythrust.text = NULL;
    segment -> connector_ythrust.value = 0;

    segment -> connector_zthrust.expr = NULL;
    segment -> connector_zthrust.text = NULL;
    segment -> connector_zthrust.value = 0;

    return segment;
}


/************************************************************************
 * Function:	 DestroySegment						*
 ************************************************************************/

void DestroySegment (segment)
    Segment segment;
{
    if (segment) {
	    Deallocate (segment -> aux);
	    Deallocate (segment -> name);
	    Deallocate (segment);
    }
}


/************************************************************************
 * Function:	 CreateBuoy						*
 ************************************************************************/

Buoy CreateBuoy (name)
    char *name;
{
    Buoy  buoy;


    if (!(buoy = AllocNew (struct buoy)))
	Fatal ("unable to allocate memory for new buoy");

    buoy -> comment = NULL;

    buoy -> aux = NULL;
    buoy -> color = NULL;
    buoy -> name = name;
    buoy -> type = 0;		/* undefined	*/
    buoy -> category = strdup("user");

    buoy -> cg       = 0.0;
    buoy -> h        = 0.0;
    buoy -> d        = 0.0;
    buoy -> m        = 0.0;
    buoy -> am       = 0.0;
    buoy -> Cdn      = 0.0;
    buoy -> Cdt      = 0.0;
    buoy -> buoyancy = 0.0;
    buoy -> draft    = 0.0;
    buoy -> S	     = 0.0;
    buoy -> w        = 0.0;
    buoy -> diameters = NULL;
    buoy -> num_diameters = 0;
    buoy -> max_draft  = 0.0;
    buoy -> min_draft  = 0.0;

    buoy -> Mmma  = 0.0;
    buoy -> Marma = 0.0;
    buoy -> Mdr   = 0.0; 

    buoy -> Cdw = 0.0;
    buoy -> Sw  = 0.0; 

    buoy -> Cl = 0.0;

    return buoy; 
}


/************************************************************************
 * Function:	 DestroyBuoy						*
 ************************************************************************/

void DestroyBuoy (buoy)
    Buoy  buoy;
{
    if (buoy) {
	    Deallocate (buoy -> name);
	    Deallocate (buoy -> comment);
	    Deallocate (buoy -> color);
	    Deallocate (buoy -> aux);
        Deallocate (buoy -> category);
	    Deallocate (buoy);
    }
}

/************************************************************************
 * Function:	 CreateConnector					*
 ************************************************************************/

Connector CreateConnector (name)
    char *name;
{
    Connector  connector;


    if (!(connector = AllocNew (struct connector)))
	    Fatal ("unable to allocate memory for new connector");

    connector -> aux = NULL;
    connector -> color = NULL;
    connector -> name = name;
    connector -> category = strdup("user");
    
    connector -> m          = 0;
    connector -> wet   	    = 0;
    connector -> d	    = 0;
    connector -> length	    = 0;
    connector -> Cdt        = 0;
    connector -> Cdn        = 0;
    connector -> am	    = 0;
    connector -> comment    = NULL;

    return connector; 
}


/************************************************************************
 * Function:	 DestroyConnector					*
 ************************************************************************/

void DestroyConnector (connector)
    Connector  connector;
{
    if (connector) {
	    Deallocate (connector -> comment);
	    Deallocate (connector -> name);
	    Deallocate (connector -> color);
	    Deallocate (connector -> aux);
        Deallocate (connector -> category);
	    Deallocate (connector);
    }
}

/************************************************************************
 * Function:	 CreateAnchor						*
 ************************************************************************/

Anchor CreateAnchor (name)
    char *name;
{
    Anchor  anchor;


    if (!(anchor = AllocNew (struct anchor)))
	    Fatal ("unable to allocate memory for new anchor");

    anchor -> aux = NULL;
    anchor -> color = NULL;
    anchor -> name = name;
    anchor -> category = strdup("user");
    anchor -> m   = 0.0;
    anchor -> am  = 0.0;
    anchor -> wet = 0.0;
    anchor -> mu  = 0.0;
    anchor -> Cdn = 0.0;
    anchor -> Cdt = 0.0;
    anchor -> S   = 0.0;
    anchor -> d   = 0.0;
    anchor -> comment = NULL;

    return anchor; 
}

/************************************************************************
 * Function:	 CreateTerminal						*
 ************************************************************************/

Terminal CreateTerminal ( )
{
    Terminal	terminal;

    if (!(terminal = AllocNew (struct terminal)))
       Fatal ("unable to allocate memory for new terminal");

    terminal -> anchor = NULL;
    terminal -> buoy = NULL;
    terminal -> loop_main_node = NULL;
    terminal -> node = NULL;

    terminal -> release = -1.0;

    terminal -> friction = 1.0;
    terminal -> safety = 2.0;

    terminal -> x = 0.0;
    terminal -> y = 0.0;
    terminal -> z = 0.0;

    terminal -> xforce = 0.0;
    terminal -> yforce = 0.0;
    terminal -> zforce = 0.0;

    terminal -> initial_xforce = 0.0;
    terminal -> initial_yforce = 0.0;
    terminal -> initial_zforce = 0.0;

    terminal -> tension = 0.0;

    terminal -> xspeed.value = 0.0;
    terminal -> xspeed.expr = NULL;
    terminal -> xspeed.text = NULL;

    terminal -> yspeed.value = 0.0;
    terminal -> yspeed.expr = NULL;
    terminal -> yspeed.text = NULL;

    terminal -> zspeed.value = 0.0;
    terminal -> zspeed.expr = NULL;
    terminal -> zspeed.text = NULL;

    terminal -> xthrust.value = 0.0;
    terminal -> xthrust.expr = NULL;
    terminal -> xthrust.text = NULL;

    terminal -> ythrust.value = 0.0;
    terminal -> ythrust.expr = NULL;
    terminal -> ythrust.text = NULL;

    terminal -> zthrust.value = 0.0;
    terminal -> zthrust.expr = NULL;
    terminal -> zthrust.text = NULL;

    terminal -> profile.value = 0.0;
    terminal -> profile.expr  = NULL;
    terminal -> profile.text  = NULL;

    terminal -> profile_H = 0.0;
    terminal -> profile_turn = 0.0;
    terminal -> profile_m = 0.0;

    terminal -> flap_file = NULL;

    terminal -> Ki = terminal -> Kd = terminal -> Kp = 0.0;

    return terminal;
}

/************************************************************************
 * Function:	 DestroyTerminal					*
 ************************************************************************/

void DestroyTerminal (terminal)
    Terminal  terminal;
{
    if (terminal) {
        FreeCode (terminal -> xspeed.expr);
        Deallocate (terminal -> xspeed.text);

        FreeCode (terminal -> yspeed.expr);
        Deallocate (terminal -> yspeed.text);

        FreeCode (terminal -> zspeed.expr);
        Deallocate (terminal -> zspeed.text);

        FreeCode (terminal -> xthrust.expr);
        Deallocate (terminal -> xthrust.text);

        FreeCode (terminal -> ythrust.expr);
        Deallocate (terminal -> ythrust.text);

        FreeCode (terminal -> zthrust.expr);
        Deallocate (terminal -> zthrust.text);

        FreeCode (terminal -> profile.expr);
        Deallocate (terminal -> profile.text);

	    if (terminal -> flap_file)
           free(terminal -> flap_file);

	    Deallocate (terminal);
    }
}

/************************************************************************
 * Function:	 DestroyAnchor						*
 ************************************************************************/

void DestroyAnchor (anchor)
    Anchor  anchor;
{
    if (anchor) {
	    Deallocate (anchor -> comment);
	    Deallocate (anchor -> name);
	    Deallocate (anchor -> color);
	    Deallocate (anchor -> aux);
        Deallocate (anchor -> category);
	    Deallocate (anchor);
    }
}


/************************************************************************
 * Function:	 CreateMaterial						*
 *									*
 * Description:	 CreateMaterial creates and initializes a new material	*
 *		 structure.  The name is assigned (not copied) and the	*
 *		 fields are initialized to zero.			*
 ************************************************************************/

Material CreateMaterial (name)
    char *name;
{
    Material material;
    int	     i;

    if (!(material = AllocNew (struct material)))
	Fatal ("unable to allocate memory for new material");

    material -> aux   = NULL;
    material -> color = NULL;
    material -> name  = name;
    material -> category = strdup("user");
    material -> EA    = 0;
    material -> EI    = 0;
    material -> GJ    = 0;
    material -> E     = 0;
    material -> I     = 0;
    material -> G     = 0;
    material -> J     = 0;
    material -> nu    = 0;
    material -> A     = 0;
    material -> d     = 0;
    material -> id    = 0;
    material -> bt     = 0;
    material -> bn     = 0;
    material -> m     = 0;
    material -> w     = 0;
    material -> length = 0;
    material -> rho   = 0;
    material -> wet   = 0;
    material -> Cmn   = 0.0;
    material -> Cmt   = 0.0;
    material -> Can   = 0.0;
    material -> Cat   = 0.0;
    material -> amn   = 0.0;
    material -> amt   = 0.0;
    material -> type  = Linear;
    material -> comment = NULL;
    material -> swl = material -> yield = 0.0;

    material -> Cdn.text = NULL;
    material -> Cdn.expr = NULL;
    material -> Cdn.value = 0;

    material -> Cdt.text = NULL;
    material -> Cdt.expr = NULL;
    material -> Cdt.value = 0;


    for (i = 0 ; i < 3 ; i++) {
       material -> T [i].text     = NULL;
       material -> T [i].expr     = NULL;
       material -> T [i].value    = 0.0;
    }

    return material;
}


/************************************************************************
 * Function:	 DestroyMaterial					*
 *									*
 * Description:	 DestroyMaterial deallocates a material structure.  The	*
 *		 name and auxillary structure are also deallocated.	*
 ************************************************************************/

void DestroyMaterial (material)
    Material material;
{
    int i;

    if (material) {
  	    Deallocate (material -> comment);
	    Deallocate (material -> name);
	    Deallocate (material -> color);
	    Deallocate (material -> aux);
        Deallocate (material -> category);

        FreeCode(material -> Cdn.expr);
        Deallocate(material -> Cdn.text);

        FreeCode(material -> Cdt.expr);
        Deallocate(material -> Cdt.text);

        for (i = 0 ; i < 3 ; i++) {
           FreeCode(material -> T [i].expr);
           Deallocate(material -> T [i].expr);
        }

	    Deallocate (material);
    }
}

char *
ExpressionString(VarExpr e)
{
    char    buff[20];

    if (e.expr)
        return e.text;
    else {
        snprintf(buff, 20, "%g", e.value);
        return strdup(buff);
    }        
}

/************************************************************************
 * Function:	AssignElevation 					*
 *									*
 * Description:	Assigns the bottom elevation expression			*
 ************************************************************************/

void AssignElevation (Environment *environment, Code expr, char *text)
{
    Deallocate (environment -> bottom_elevation.text);
    FreeCode (environment -> bottom_elevation.expr);

    environment -> bottom_elevation.value = EvalCode (expr, NULL, 0, 0, 
                                                      0, 0, 0, CURRNODEDATA);
    environment -> bottom_elevation.expr  = IsConstant (expr) ? NULL : CopyCode (expr);
    environment -> bottom_elevation.text  = text ? strdup (text) : NULL;

   return;
}

void 
AssignSmoothing(Analysis *a, Code expr, char *text)
{
    Deallocate(a -> dynamic_var_smooth.text);
    Deallocate(a -> dynamic_var_smooth.expr);
    a -> dynamic_var_smooth.value = EvalCode(expr, NULL, 0, 0, 0, 0, 0, CURRNODEDATA);
    a -> dynamic_var_smooth.expr = IsConstant(expr) ? NULL : CopyCode(expr);
    a -> dynamic_var_smooth.text = text ? strdup(text) : NULL;
}

/************************************************************************
 * Function:	AssignCurrentModulation					*
 *									*
 * Description:	Assigns a current modulation given as a piece of code	*
 ************************************************************************/

void 
AssignCurrentModulation (Environment *environment, 
                         char dof, Code expr, char *text)
{
    if (dof == 'z') {
        Deallocate (environment -> Uz.text);
        FreeCode (environment -> Uz.expr);

        environment -> Uz_mod.value = EvalCode (expr, NULL, 0, 0, 
                                                      0, 0, 0, CURRNODEDATA);
        environment -> Uz_mod.expr  = IsConstant (expr) ? NULL : CopyCode (expr);
        environment -> Uz_mod.text  = text ? strdup (text) : NULL;
    }
    else if (dof == 'y') {
        Deallocate (environment -> Uy.text);
        FreeCode (environment -> Uy.expr);

        environment -> Uy_mod.value = EvalCode (expr, NULL, 0, 0, 
                                                      0, 0, 0, CURRNODEDATA);
        environment -> Uy_mod.expr  = IsConstant (expr) ? NULL : CopyCode (expr);
        environment -> Uy_mod.text  = text ? strdup (text) : NULL;
    }
    else if (dof == 'x') {
        Deallocate (environment -> Ux.text);
        FreeCode (environment -> Ux.expr);

        environment -> Ux_mod.value = EvalCode (expr, NULL, 0, 0,
                                                      0, 0, 0, CURRNODEDATA);
        environment -> Ux_mod.expr  = IsConstant (expr) ? NULL : CopyCode (expr);
        environment -> Ux_mod.text  = text ? strdup (text) : NULL;
    }
}


/************************************************************************
 * Function:	AssignCurrent						*
 *									*
 * Description:	Assigns a current given as a piece of codea		*
 ************************************************************************/

void 
AssignCurrent (Environment *environment, char dof, Code expr, char *text)
{
    if (dof == 'z') {
        Deallocate (environment -> Uz.text);
        FreeCode (environment -> Uz.expr);

        environment -> Uz.value = EvalCode (expr, NULL, 0, 0,
                                                      0, 0, 0, CURRNODEDATA);
        environment -> Uz.expr  = IsConstant (expr) ? NULL : CopyCode (expr);
        environment -> Uz.text  = text ? strdup (text) : NULL;
    }
    else if (dof == 'y') {
        Deallocate (environment -> Uy.text);
        FreeCode (environment -> Uy.expr);
        
        environment -> Uy.value = EvalCode (expr, NULL, 0, 0,
                                                      0, 0, 0, CURRNODEDATA);
        environment -> Uy.expr  = IsConstant (expr) ? NULL : CopyCode (expr);
        environment -> Uy.text  = text ? strdup (text) : NULL;
    }
    else if (dof == 'x') {
        Deallocate (environment -> Ux.text);
        FreeCode (environment -> Ux.expr);

        environment -> Ux.value = EvalCode (expr, NULL, 0, 0,
                                                      0, 0, 0, CURRNODEDATA);
        environment -> Ux.expr  = IsConstant (expr) ? NULL : CopyCode (expr);
        environment -> Ux.text  = text ? strdup (text) : NULL;
    }
}

/************************************************************************
 * Function:	AssignWind						*
 *									*
 * Description:	Assigns a wind speed given as a piece of codea		*
 ************************************************************************/

void 
AssignWind (Environment *environment, char dof, Code expr, char *text)
{
    if (dof == 'z') {
        Deallocate (environment -> z_wind.text);
        FreeCode (environment -> z_wind.expr);

        environment -> z_wind.value = EvalCode (expr, NULL, 0, 0,
                                                      0, 0, 0, CURRNODEDATA);
        environment -> z_wind.expr  = IsConstant (expr) ? NULL : CopyCode (expr);
        environment -> z_wind.text  = text ? strdup (text) : NULL;
    }
    else if (dof == 'y') {
        Deallocate (environment -> y_wind.text);
        FreeCode (environment -> y_wind.expr);

        environment -> y_wind.value = EvalCode (expr, NULL, 0, 0,
                                                      0, 0, 0, CURRNODEDATA);
        environment -> y_wind.expr  = IsConstant (expr) ? NULL : CopyCode (expr);
        environment -> y_wind.text  = text ? strdup (text) : NULL;
    }
}

/************************************************************************
 * Function:	AssignSpeed						*
 *									*
 * Description:	Assigns a speed given as a piece of code		*
 ************************************************************************/

void AssignSpeed (dof, term, expr, text)
    char      dof;
    Terminal  term;
    Code      expr;
    char     *text;
{
    if (dof == 'x') {
        Deallocate (term -> xspeed.text);
        FreeCode (term -> xspeed.expr);

        term -> xspeed.value = EvalCode (expr, NULL, 0, 0,
                                                      0, 0, 0, CURRNODEDATA);
        term -> xspeed.expr  = IsConstant (expr) ? NULL : CopyCode (expr);
        term -> xspeed.text  = text ? strdup (text) : NULL;
    }
    else if (dof == 'y') {
        Deallocate (term -> yspeed.text);
        FreeCode (term -> yspeed.expr);

        term -> yspeed.value = EvalCode (expr, NULL, 0, 0,
                                                      0, 0, 0, CURRNODEDATA);
        term -> yspeed.expr  = IsConstant (expr) ? NULL : CopyCode (expr);
        term -> yspeed.text  = text ? strdup (text) : NULL;
    }
    else if (dof == 'z') {
        Deallocate (term -> zspeed.text);
        FreeCode (term -> zspeed.expr);

        term -> zspeed.value = EvalCode (expr, NULL, 0, 0,
                                                      0, 0, 0, CURRNODEDATA);
        term -> zspeed.expr  = IsConstant (expr) ? NULL : CopyCode (expr);
        term -> zspeed.text  = text ? strdup (text) : NULL;
    }
}

/************************************************************************
 * Function:	AssignConnectorThrust						
 *									
 * Description:	Assigns a thrust given as a piece of code		
 ************************************************************************/

void AssignConnectorThrust (dof, seg, expr, text)
    char      dof;
    Segment   seg;
    Code      expr;
    char     *text;
{
    if (dof == 'x') {
        Deallocate (seg -> connector_xthrust.text);
        FreeCode (seg -> connector_xthrust.expr);

        seg -> connector_xthrust.value = EvalCode (expr, NULL, 0, 0,
                                                      0, 0, 0, CURRNODEDATA);
        seg -> connector_xthrust.expr  = IsConstant (expr) ? NULL : CopyCode (expr);
        seg -> connector_xthrust.text  = text ? strdup (text) : NULL;
    }
    else if (dof == 'y') {
        Deallocate (seg -> connector_ythrust.text);
        FreeCode (seg -> connector_ythrust.expr);

        seg -> connector_ythrust.value = EvalCode (expr, NULL, 0, 0,
                                                      0, 0, 0, CURRNODEDATA);
        seg -> connector_ythrust.expr  = IsConstant (expr) ? NULL : CopyCode (expr);
        seg -> connector_ythrust.text  = text ? strdup (text) : NULL;
    }
    else if (dof == 'z') {
        Deallocate (seg -> connector_zthrust.text);
        FreeCode (seg -> connector_zthrust.expr);

        seg -> connector_zthrust.value = EvalCode (expr, NULL, 0, 0,
                                                      0, 0, 0, CURRNODEDATA);
        seg -> connector_zthrust.expr  = IsConstant (expr) ? NULL : CopyCode (expr);
        seg -> connector_zthrust.text  = text ? strdup (text) : NULL;
    }
}

/************************************************************************
 * Function:	AssignThrust						*
 *									*
 * Description:	Assigns a thrust given as a piece of code		*
 ************************************************************************/

void AssignThrust (dof, term, expr, text)
    char      dof;
    Terminal  term;
    Code      expr;
    char     *text;
{
    if (dof == 'x') {
        Deallocate (term -> xthrust.text);
        FreeCode (term -> xthrust.expr);
        
        term -> xthrust.value = EvalCode (expr, NULL, 0, 0,
                                                      0, 0, 0, CURRNODEDATA);
        term -> xthrust.expr  = IsConstant (expr) ? NULL : CopyCode (expr);
        term -> xthrust.text  = text ? strdup (text) : NULL;
    }
    else if (dof == 'y') {
        Deallocate (term -> ythrust.text);
        FreeCode (term -> ythrust.expr);

        term -> ythrust.value = EvalCode (expr, NULL, 0, 0,
                                                      0, 0, 0, CURRNODEDATA);
        term -> ythrust.expr  = IsConstant (expr) ? NULL : CopyCode (expr);
        term -> ythrust.text  = text ? strdup (text) : NULL;
    }
    else if (dof == 'z') {
        Deallocate (term -> zthrust.text);
        FreeCode (term -> zthrust.expr);

        term -> zthrust.value = EvalCode (expr, NULL, 0, 0,
                                                      0, 0, 0, CURRNODEDATA);
        term -> zthrust.expr  = IsConstant (expr) ? NULL : CopyCode (expr);
        term -> zthrust.text  = text ? strdup (text) : NULL;
    }
}

/************************************************************************
 * Function:	AssignTension						*
 *									*
 * Description:	Assigns a tension function given as a piece of code	*
 ************************************************************************/

void AssignTension (deriv, material, expr, text)
    int	      deriv;
    Material  material;
    Code      expr;
    char     *text;
{
     Deallocate (material -> T [deriv].text);
     FreeCode (material -> T [deriv].expr);

     material -> T [deriv].value = EvalCode (expr, NULL, 0, 0,
                                                      0, 0, 0, CURRNODEDATA);
     material -> T [deriv].expr  = IsConstant (expr) ? NULL : CopyCode (expr);
     material -> T [deriv].text  = text ? strdup (text) : NULL;
}

void AssignDrag(int dof, Material material, Code expr, char *text)
{
    if (dof == 0) {
 	    Deallocate(material -> Cdt.text);
 	    Deallocate(material -> Cdt.expr);

	    material -> Cdt.value = EvalCode (expr, NULL, 0, 0,
                                                      0, 0, 0, CURRNODEDATA);
	    material -> Cdt.expr = IsConstant(expr) ? NULL : CopyCode(expr);
	    material -> Cdt.text = text ? strdup(text) : NULL;
    } 
    else if (dof == 1) {
 	    Deallocate(material -> Cdn.text);
 	    Deallocate(material -> Cdn.expr);

	    material -> Cdn.value = EvalCode (expr, NULL, 0, 0,
                                                      0, 0, 0, CURRNODEDATA);
	    material -> Cdn.expr = IsConstant(expr) ? NULL : CopyCode(expr);
	    material -> Cdn.text = text ? strdup(text) : NULL;
    } 
}

void AssignSegmentPay (s, top, expr, text)
    Segment   s;
    int	      top;
    Code      expr;
    char     *text;
{
    VarExpr   *e;

    
    if (top) 
       e = &(s -> top_pay);
    else
       e = &(s -> bottom_pay);

    Deallocate (e -> text);
    FreeCode (e -> expr);

    e -> value = EvalCode (expr, NULL, 0, 0,
                           0, 0, 0, CURRNODEDATA);
    e -> expr  = IsConstant (expr) ? NULL : CopyCode (expr);
    e -> text  = text ? strdup (text) : NULL;
}

void AssignProfile (term, expr, text)
    Terminal  term;
    Code      expr;
    char     *text;
{
    Deallocate (term -> profile.text);
    FreeCode (term -> profile.expr);

    term -> profile.value = EvalCode (expr, NULL, 0, 0,
                                      0, 0, 0, CURRNODEDATA);
    term -> profile.expr  = IsConstant (expr) ? NULL : CopyCode (expr);
    term -> profile.text  = text ? strdup (text) : NULL;
}
