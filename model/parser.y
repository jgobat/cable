%{
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
 * File:	parser.y						
 *									
 * Description:	This file contains the yacc specification for buoy	
 ************************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include "code.h"
# include "error.h"
# include "problem.h"
# include "objects.h"
# include "compress.h"
# include "allocate.h"
# define lint

extern ParserControl ctl;

extern void yyerror ( );
extern int  yylex  ( );

/* Stuff to handle variable expressions */

static int z_pos;

static char *varExName;

/* Current objects (inherited attributes). */

static Item		  found;		/* current object	*/
static Buoy		  buoy;			/* current buoy		*/
static Material		  material;		/* current material	*/
static Anchor		  anchor;		/* current anchor	*/
static Connector	  connector;	/* current connector 	*/
static Segment		  segment;		/* current segment	*/
static Terminal		  terminal;		/* current terminal	*/
static Branch		  branch;		/* current branch	*/

static int		  num_terminals;	/* terminal count	*/
static int		  num_segs;		/* segment count	*/

static Segment		  segment_from;
static int		  num_branches;
static int		  num_branch_segment;

	/*
	 * for sort of matlab like specification of regularly
	 * spaced attachments
	 */

static int attach_start;
static int attach_step;
static int attach_stop;

/* Dummy strucutures.  If an object is defined twice, the second definition
   is illegal and the current object is set to a dummy object so that
   we can get through the rest of the definition for this unused second
   instance without a fuss. */

static struct buoy	  dummy_buoy;		/* dummy buoy		   */
static struct anchor 	  dummy_anchor;		/* dummy anchor		   */
static struct material	  dummy_material;	/* dummy material	   */
static struct connector	  dummy_connector;	/* dummy connector	   */


/* Temporary arrays. */

static distribution       dist_array [1024];    
static distribution      *dist_ptr;             
static DiameterPair       pair_array [1024];    /* temporary pair array    */
static DiameterPair      *pair_ptr;             /* pointer into array      */
static attachment	  attach_array [1024];
static attachment	 *attach_ptr;
static int		  int_array [1024];
static int		 *int_ptr;
static wave		  wave_array [10];
static wave		 *wave_ptr;



/* Discrete expression variables. */

static int		  table_error = 0;	/* error indicator	   */
static double		  last_depth = 0;	/* last depth coordinate   */
static double		 *table = NULL;		/* table of values	   */
static unsigned		  table_count = 0;	/* count of values	   */
static unsigned		  table_size = 0;	/* size of table	   */

%}


%union {
    int           i;
    double        d;
    DiameterPair  p;
    distribution  r;
    wave          w;
    attachment    a;
    char         *s;
}

%right	'?' ':'
%left	OR
%left	AND
%left	'|'
%left	'^'
%left	'&'
%left	EQUALS NEQUAL
%left	'<' '>' LT_EQ GT_EQ
%left	LSHIFT RSHIFT
%left	'+' '-'
%left	'*' '/' '%'
%right	UNARY '!' '~'


%token	NAME INTEGER DOUBLE BUOY_TYPE MATERIAL_TYPE 
%token  VAR_t VAR_s VAR_n

%token  VAR_p VAR_H VAR_e VAR_x VAR_y VAR_z VAR_T VAR_u VAR_v VAR_w
%token  VAR_Fx VAR_Fy VAR_Fz VAR_Ux VAR_Vy VAR_Wz
%token  VAR_Om2 VAR_Om1 VAR_Om3 VAR_Sn VAR_Sb VAR_phi
%token  VAR_B0 VAR_B1 VAR_B2 VAR_B3

%token  VAR_p_f VAR_H_f VAR_e_f VAR_x_f VAR_y_f VAR_z_f 
%token  VAR_T_f VAR_u_f VAR_v_f VAR_w_f
%token  VAR_Fx_f VAR_Fy_f VAR_Fz_f VAR_Ux_f VAR_Vy_f VAR_Wz_f
%token  VAR_Om2_f VAR_Om1_f VAR_Om3_f VAR_Sn_f VAR_Sb_f VAR_phi_f
%token  VAR_B0_f VAR_B1_f VAR_B2_f VAR_B3_f

%token  CONST_PI
%token  FORCING_METHOD INPUT_TYPE PROBLEM_TYPE INTEGRATION_METHOD
%token  CONNECTION_TYPE
%token  SOLUTION_METHOD
%token  OUTER_METHOD

%token  SEG FIRSTACTIVE LASTACTIVE FIRST LAST FILTER
%token	SIN COS TAN POW EXP LOG LOG10 SQRT HYPOT FLOOR CEIL FMOD FABS
%token  SINH COSH TANH ACOS ASIN ATAN2

%token	PROBLEM MATERIALS ANCHORS CONNECTORS 
%token  LAYOUT ENVIRONMENT ANALYSIS END 
%token  BUOYS

%token	TITLE_EQ TYPE_EQ

%token	BUOYANCY_EQ CG_EQ H_EQ CDT_EQ CDN_EQ DRAFT_EQ 
%token  DIAMETERS_EQ
%token  CL_EQ CDW_EQ SW_EQ CMN_EQ CMT_EQ AMN_EQ AMT_EQ CAT_EQ CAN_EQ
%token  MU_EQ 
%token  INITIAL_XFORCE_EQ INITIAL_YFORCE_EQ INITIAL_ZFORCE_EQ
%token  SAFETY_EQ TENSION_EQ XFORCE_EQ YFORCE_EQ ZFORCE_EQ 
%token	XSPEED_EQ YSPEED_EQ ZSPEED_EQ RELEASE_TIME_EQ PAYRATE_EQ 
%token	PROFILE_EQ PROFILE_TURN_DEPTH_EQ PROFILE_RANGE_EQ PROFILE_SLOPE_EQ
%token  FLAP_FILE_EQ
%token  KI_EQ KD_EQ KP_EQ 
%token  XTHRUST_EQ YTHRUST_EQ ZTHRUST_EQ
%token  BUOYANCY_HEAD_EQ

%token	CONNECTOR_EQ BODY_EQ SEGMENT_EQ ANCHOR_EQ NODE_EQ TERMINAL_EQ BRANCH_EQ 
%token	MATERIAL_EQ LENGTH_EQ NODES_EQ ATTACHMENTS_EQ NAME_EQ
%token  TOP_PAY_EQ BOTTOM_PAY_EQ TOP_LENGTH_EQ BOTTOM_LENGTH_EQ
%token  TOP_NODES_EQ BOTTOM_NODES_EQ
%token  BUOY_EQ

%token  COMMENT_EQ CATEGORY_EQ
%token	M_EQ BT_EQ BN_EQ EA_EQ D_EQ COLOR_EQ EI_EQ GJ_EQ A_EQ 
%token  WET_EQ W_EQ AM_EQ
%token  B_EQ SWL_EQ YIELD_EQ
%token  X_EQ Y_EQ Z_EQ
%token  T_EQ TE_EQ TEE_EQ
%token  E_EQ G_EQ I_EQ J_EQ NU_EQ ID_EQ

%token	TOLERANCE_EQ MAX_ITERATIONS_EQ RELAXATION_EQ RAMP_TIME_EQ
%token  DYNAMIC_TOLERANCE_EQ STATIC_TOLERANCE_EQ STATIC_OUTER_TOLERANCE_EQ
%token  DYNAMIC_RELAXATION_EQ STATIC_RELAXATION_EQ OUTER_RELAXATION_EQ
%token  ALPHA_K_EQ ALPHA_M_EQ DYN_RHO_EQ GAMMA_EQ
%token  STATIC_ITERATIONS_EQ DYNAMIC_ITERATIONS_EQ STATIC_OUTER_ITERATIONS_EQ
%token  STATIC_OUTER_METHOD_EQ
%token  SHOOTING_ITERATIONS_EQ
%token	DURATION_EQ TIME_STEP_EQ SAMPLE_RATE_EQ OUTPUT_NODES_EQ
%token  TIME_STEP_ADAPT_EQ TIME_STEP_LEVELS_EQ
%token  CURRENT_STEPS_EQ DYNAMIC_INTEGRATION_EQ
%token  STATIC_SOLUTION_EQ STATIC_INITIAL_GUESS_EQ
%token  MESH_SMOOTH_EQ MESH_AMPLIFICATION_EQ
%token  RELAX_ADAPT_UP_EQ RELAX_ADAPT_DOWN_EQ RELAX_STALL_LIMIT_EQ
%token  VIVA_TIME_STEP_EQ VIVA_DECIMATE_EQ VIVA_ITERATIONS_EQ
%token  DYNAMIC_VAR_SMOOTHING_EQ

%token	ZCURRENT_EQ YCURRENT_EQ XCURRENT_EQ	XINPUT_EQ YINPUT_EQ ZINPUT_EQ 
%token  ZCURRENT_MODULATION_EQ XCURRENT_MODULATION_EQ YCURRENT_MODULATION_EQ
%token  XCURRENT_PROFILE_EQ YCURRENT_PROFILE_EQ CURRENT_SCALE_EQ
%token  CURRENT_ROTATION_EQ
%token  CURRENT_X_SCALE_EQ CURRENT_Y_SCALE_EQ CURRENT_Z_SCALE_EQ
%token	XWAVE_EQ YWAVE_EQ YWIND_EQ XWIND_EQ
%token  DEPTH_EQ GRAVITY_EQ RHO_EQ 
%token	WAVE_FILE_EQ VELOCITY_FILE_EQ CURRENT_FILE_EQ
%token  INPUT_TYPE_EQ FORCING_METHOD_EQ
%token  BOTTOM_STIFFNESS_EQ BOTTOM_DAMPING_EQ BOTTOM_ELEVATION_EQ
%token  BOTTOM_FRICTION_EQ
%token  FILL_FLUID_DENSITY_EQ

%type	<i> INTEGER MATERIAL_TYPE BUOY_TYPE CONNECTION_TYPE
%type   <i> FORCING_METHOD INPUT_TYPE PROBLEM_TYPE INTEGRATION_METHOD
%type   <i> SOLUTION_METHOD
%type   <i> OUTER_METHOD
%type	<d> DOUBLE constant_expression
%type	<s> NAME
%type   <p> buoy_diameter
%type   <r> distribution_pair
%type   <w> wave_triple
%type   <a> attachment
%type	<i> expression function or_action and_action if_action else_action 
            
%%

specification
	: initialize problem_description section_list END
	;


initialize
	: /* empty */
	    {
 		num_segs = 0;
                num_terminals = 0;
                num_branches = 0;
		num_branch_segment = 0;
		branch = NULL;
		segment_from = NULL;
	    }
	;


/* Problem description */

problem_description
	: PROBLEM problem_parameter_list
	| /* empty */
	;


problem_parameter_list
	: problem_parameter_list problem_parameter
	| /* empty */
	;


problem_parameter
	: TITLE_EQ NAME
	    {
		Deallocate (ctl.problem -> title);
		ctl.problem -> title = $2;
	    }

        | TYPE_EQ PROBLEM_TYPE
	    {
		ctl.problem -> type = $2;
	    }

	| error
	;


/* Sections */

section_list
	: section_list section
	| /* empty */
	;


section
	: buoy_section
        | anchor_section
        | connector_section
	| material_section
	| layout_section
	| analysis_section
	| environment_section
	| END
	;


/* Buoy section */

buoy_section
	: BUOYS buoy_definition_list
	;


buoy_definition_list
	: buoy_definition_list buoy_definition
	| /* empty */
	;


buoy_definition
	: buoy_name buoy_parameter_list
	;


buoy_name
	: NAME
	    {
		    buoy = CreateBuoy ($1);
		    found = TreeInsert (ctl.problem -> buoy_tree, buoy);

		    if (found != (Item) buoy) {
		        error ("buoy %d is previously defined", $1);
		        DestroyBuoy (buoy);
		        buoy = &dummy_buoy;
		    }
	    }
	;

buoy_parameter_list
	: buoy_parameter_list buoy_parameter
	| /* empty */
	;


buoy_parameter
	: COLOR_EQ NAME
	    {
		    Deallocate (buoy  -> color);
            buoy -> color = $2;
	    }

	| CATEGORY_EQ NAME
	    {
            Deallocate (buoy -> category);
		    buoy -> category = $2;
	    }

	| COMMENT_EQ NAME
	    {
            Deallocate (buoy -> comment);
		    buoy -> comment = $2;
	    }

   | TYPE_EQ BUOY_TYPE
        {
            buoy -> type = $2;
            // fprintf(stderr,"set buoy type = %d\n", buoy -> type);
        }

	| D_EQ constant_expression
	    {
		buoy -> d = $2;
	    }

        | DRAFT_EQ constant_expression
            {
                buoy -> draft = $2;
            }

	| CG_EQ constant_expression
	    {
		buoy -> cg = $2;
	    }

	| H_EQ constant_expression
	    {
		buoy -> h = $2;
	    }

	| BUOYANCY_EQ constant_expression
	    {
		buoy -> buoyancy = $2;
	    }

	| CL_EQ constant_expression
	    {
		buoy -> Cl = $2;
	    }

        | CDW_EQ constant_expression
            {
                buoy -> Cdw = $2;
            }

        | SW_EQ constant_expression
            {
                buoy -> Sw = $2;
            }

        | CDN_EQ constant_expression
            {
                buoy -> Cdn = $2;
            }

        | CDT_EQ constant_expression
            {
                buoy -> Cdt = $2;
            }

        | M_EQ constant_expression
            {
                buoy -> m = $2;
            }

        | AM_EQ constant_expression
            {
                buoy -> am = $2;
            }

        | B_EQ constant_expression
            {
                buoy -> b = $2;
            }

        | DIAMETERS_EQ buoy_diameter_list
            {
                unsigned i;
                unsigned size;

                if (buoy == &dummy_buoy)
                    break;

                size = pair_ptr - pair_array;

                if (!(buoy -> diameters = Allocate (DiameterPair, size)))
                    Fatal ("unable to allocate memory for diameter pairs");

                UnitOffset (buoy -> diameters);
                buoy -> num_diameters = size;

                for (i = 1; i <= size; i ++)
                    buoy -> diameters [i] = pair_array [i - 1];
            }

        | error
        ;

buoy_diameter_list
        : buoy_diameter_list ',' buoy_diameter
            {
                *pair_ptr ++ = $3;
            }

        | buoy_diameter_list buoy_diameter
            {
                *pair_ptr ++ = $2;
            }

        | buoy_diameter
            {
                pair_ptr = pair_array;
                *pair_ptr ++ = $1;
            }
        ;

buoy_diameter
        : '(' constant_expression ',' constant_expression ')'
            {
                $$.level = $2;
                $$.d     = $4;
            }
        ;


/* Anchor section */

anchor_section
	: ANCHORS anchor_definition_list
	;


anchor_definition_list
	: anchor_definition_list anchor_definition
	| /* empty */
	;


anchor_definition
	: anchor_name anchor_parameter_list
	;


anchor_name
	: NAME
	    {
		anchor = CreateAnchor ($1);
		found = TreeInsert (ctl.problem -> anchor_tree, anchor);

		if (found != (Item) anchor) {
		    error ("anchor %s is previously defined", $1);
		    DestroyAnchor (anchor);
		    anchor = &dummy_anchor;
		}
	    }
	;

anchor_parameter_list
	: anchor_parameter_list anchor_parameter
	| /* empty */
	;


anchor_parameter
	: COLOR_EQ NAME
	    {
    		Deallocate (anchor -> color);
            anchor -> color = $2;
	    }

	| CATEGORY_EQ NAME
	    {
            Deallocate (anchor -> category);
		    anchor -> category = $2;
	    }

	| COMMENT_EQ NAME
	    {
            Deallocate (anchor -> comment);
    		anchor -> comment = $2;
	    }
    | CDT_EQ constant_expression
        {
            anchor -> Cdt = $2;
        }
        | CDN_EQ constant_expression
            {
                anchor -> Cdn = $2;
            }
        | WET_EQ constant_expression
            {
                anchor -> wet = $2;
            }
        | D_EQ constant_expression
            {
                anchor -> d = $2;
            }
        | M_EQ constant_expression
            {
                anchor -> m = $2;
            }
	| MU_EQ constant_expression
	    {
		anchor -> mu = $2;
 	    }

	| error
	;


/* Connector section */

connector_section
	: CONNECTORS connector_definition_list
	;


connector_definition_list
	: connector_definition_list connector_definition
	| /* empty */
	;


connector_definition
	: connector_name connector_parameter_list
	;


connector_name
	: NAME
	    {
		connector = CreateConnector ($1);
		found = TreeInsert (ctl.problem -> connector_tree, connector);

		if (found != (Item) connector) {
		    error ("connector %s is previously defined", $1);
		    DestroyConnector (connector);
		    connector = &dummy_connector;
		}
	    }
	;

connector_parameter_list
	: connector_parameter_list connector_parameter
	| /* empty */
	;


connector_parameter
	: COLOR_EQ NAME
	    {
		    Deallocate (connector -> color);
            connector -> color = $2;
	    }

	| CATEGORY_EQ NAME
	    {
            Deallocate (connector -> category);
		    connector -> category = $2;
	    }

	| CDT_EQ constant_expression
	    {
		connector -> Cdt = $2;
	    }

	| CDN_EQ constant_expression
	    {
		connector -> Cdn = $2;
	    }

	| WET_EQ constant_expression
	    {
		connector -> wet = $2;
	    }

	| D_EQ constant_expression
	    {
		connector -> d = $2;
	    }

	| LENGTH_EQ constant_expression
	    {
		connector -> length = $2;
	    }

	| M_EQ constant_expression
	    {
		connector -> m = $2;
	    }

	| AM_EQ constant_expression
	    {
		connector -> am = $2;
	    }

	| COMMENT_EQ NAME
	    {
            Deallocate (connector -> comment);
		    connector -> comment = $2;
	    }

	| error
	;


/* Material section */

material_section
	: MATERIALS material_definition_list
    {
        // printf("inside materials\n");
    }
	;


material_definition_list
	: material_definition_list material_definition
	| /* empty */
	;


material_definition
	: material_name material_parameter_list
    {
        // printf("material definition\n");
    }
	;


material_name
	: NAME
	    {
		material = CreateMaterial ($1);
        // printf("reading material %s\n", $1);
		found = TreeInsert (ctl.problem -> material_tree, material);
                
		    if (found != (Item) material) {
		        error ("material %s is previously defined", $1);
		        DestroyMaterial (material);
		        material = &dummy_material;
		    }
	    }
	;

material_parameter_list
	: material_parameter_list material_parameter
	| /* empty */
	;


material_parameter
	: COLOR_EQ NAME
	    {
		    Deallocate (material -> color);
            material -> color = $2;
	    }

	| CATEGORY_EQ NAME
	    {
            Deallocate (material -> category);
		    material -> category = $2;
	    }

	| EA_EQ constant_expression
	    {
		    material -> EA = $2;
	    }

	| EI_EQ constant_expression
	    {
		    material -> EI = $2;
	    }

	| GJ_EQ constant_expression
	    {
		    material -> GJ = $2;
	    }

	| D_EQ constant_expression
	    {
		    material -> d = $2;
	    }

	| E_EQ constant_expression
	    {
		    material -> E = $2;
	    }

	| G_EQ constant_expression
	    {
		    material -> G = $2;
	    }

	| J_EQ constant_expression
	    {
		    material -> J = $2;
	    }

	| I_EQ constant_expression
	    {
		    material -> I = $2;
	    }

	| NU_EQ constant_expression
	    {
		    material -> nu = $2;
	    }


	| ID_EQ constant_expression
	    {
		    material -> id = $2;
	    }

	| A_EQ constant_expression
	    {
		    material -> A = $2;
	    }

	| W_EQ constant_expression
	    {
		    material -> w = $2;
	    }

	| WET_EQ constant_expression
	    {
		    material -> wet = $2;
	    }

	| LENGTH_EQ constant_expression
	    {
		    material -> length = $2;
	    }

	| AM_EQ constant_expression
	    {
		    material -> amn = $2;
	    }

    | AMN_EQ constant_expression
        {
		    material -> amn = $2;
	    }

    | AMT_EQ constant_expression
        {
		    material -> amt = $2;
	    }

	| CDN_EQ  enable_copy variable_expression
	    {
		    AssignDrag(1, material, InCore, copy_input (0));
	    }


	| CDT_EQ  enable_copy variable_expression
	    {
		    AssignDrag(0, material, InCore, copy_input (0));
	    }

    | CMN_EQ constant_expression
            {
		material -> Cmn = $2;
	    }

    | CMT_EQ constant_expression
        {
		    material -> Cmt = $2;
	    }

    | CAN_EQ constant_expression
        {
		    material -> Can = $2;
	    }

    | CAT_EQ constant_expression
        {
		    material -> Cat = $2;
	    }

	| M_EQ constant_expression
	    {
            material -> m = $2;
	    }

    | BN_EQ constant_expression
        {
            material -> bn = $2;
        }

    | BT_EQ constant_expression
        {
            material -> bt = $2;
        }

	| T_EQ  enable_copy variable_expression
            {
                AssignTension (0, material, InCore, copy_input (0));
            }

	| TE_EQ  enable_copy variable_expression
            {
                AssignTension (1, material, InCore, copy_input (0));
            }

	| TEE_EQ  enable_copy variable_expression
            {
                AssignTension (2, material, InCore, copy_input (0));
            }

    | TYPE_EQ MATERIAL_TYPE
	    {
	        material -> type = $2;
        }

	| YIELD_EQ constant_expression
	    {
		    material -> yield = $2;
	    }

	| SWL_EQ constant_expression
	    {
		    material -> swl = $2;
	    }

	| COMMENT_EQ NAME
	    {
            Deallocate (material -> comment);
		    material -> comment = $2;
	    }

	| error
	;


/* Layout section */

layout_section
	: LAYOUT layout_terminal layout_segment_connector_list layout_terminal
	;

layout_terminal
	: terminal_start '{' terminal_parameter_list '}'
        ;

terminal_start
	: TERMINAL_EQ
            {
                 terminal = CreateTerminal ( );
                 if (branch) {
                    branch -> terminal = terminal;
                    branch -> num_segment = num_branch_segment;
                    branch = NULL;
                 }
		 else {
                    num_terminals ++;
                    ctl.problem -> terminal [num_terminals] = terminal;
		 }
            }
        | error
 	;

terminal_parameter_list
	: terminal_parameter_list terminal_parameter
	| /* empty */
	;

terminal_parameter
	: BUOY_EQ NAME
            {   
                terminal -> buoy = (Buoy) $2;
                
            }  

	| ANCHOR_EQ NAME
	    {
		terminal -> anchor = (Anchor) $2;
	    }

	| NODE_EQ INTEGER
	    {
		terminal -> loop_main_node = (Node) $2;
	    }

	| RELEASE_TIME_EQ constant_expression
	    {
		terminal -> release = $2;
	    }

	| X_EQ constant_expression
	    {
		/* user X is internal Y */
		terminal -> y = $2;
	    }

	| Y_EQ constant_expression
	    {
		/* user Y is internal Z */
		terminal -> z = $2;
	    }

	| Z_EQ constant_expression
	    {
		/* user Z is internal X */
		terminal -> x = $2;
	    }

	| SAFETY_EQ constant_expression
  	    {
		terminal -> safety = $2;
	    }

	| MU_EQ constant_expression
  	    {
		terminal -> friction = $2;
	    }

	| TENSION_EQ constant_expression
  	    {
		terminal -> tension = $2;
	    }

	| INITIAL_XFORCE_EQ constant_expression
	    {
		/* user X is internal Y */
		terminal -> initial_yforce = $2;
	    }

	| INITIAL_YFORCE_EQ constant_expression
	    {
		/* user Y is internal Z */
		terminal -> initial_zforce = $2;
	    }

	| INITIAL_ZFORCE_EQ constant_expression
	    {
		/* user Z is internal X */
		terminal -> initial_xforce = $2;
	    }


	| XFORCE_EQ constant_expression
	    {
		/* user X is internal Y */
		terminal -> yforce = $2;
	    }

	| YFORCE_EQ constant_expression
	    {
		/* user Y is internal Z */
		terminal -> zforce = $2;
	    }

	| ZFORCE_EQ constant_expression
	    {
		/* user Z is internal X */
		terminal -> xforce = $2;
	    }

	| XSPEED_EQ  enable_copy variable_expression
            {
		/* user X is internal Y */
                AssignSpeed ('y', terminal, InCore, copy_input (0));
            }

	| YSPEED_EQ  enable_copy variable_expression
            {
		/* user Y is internal Z */
                AssignSpeed ('z', terminal, InCore, copy_input (0));  
            }

	| ZSPEED_EQ  enable_copy variable_expression
            {
		/* user Z is internal X */
                AssignSpeed ('x', terminal, InCore, copy_input (0));  
            }

	| XTHRUST_EQ  enable_copy variable_expression
            {
		        /* user X is internal Y */
                AssignThrust ('y', terminal, InCore, copy_input (0));
            }

	| YTHRUST_EQ  enable_copy variable_expression
            {
		        /* user Y is internal Z */
                AssignThrust ('z', terminal, InCore, copy_input (0));  
            }

	| ZTHRUST_EQ {z_pos = 1;} enable_copy variable_expression
            {
		        /* user Z is internal X */
                AssignThrust ('x', terminal, InCore, copy_input (0));  
                z_pos = 0;
            }
            

        | FLAP_FILE_EQ NAME
            {
	 	Deallocate (terminal -> flap_file);
		terminal -> flap_file = $2;
	    }

        | PROFILE_TURN_DEPTH_EQ constant_expression
	    {
		terminal -> profile_turn = $2;
            }

        | PROFILE_RANGE_EQ constant_expression
	    {
		terminal -> profile_H = $2;
            }

        | PROFILE_SLOPE_EQ constant_expression
	    {
		terminal -> profile_m = $2;
            }

	| PROFILE_EQ  enable_copy variable_expression
	    {
		AssignProfile(terminal, InCore, copy_input(0));
	    }
	| KI_EQ constant_expression
	    {
		terminal -> Ki = $2;
	    }
	| KD_EQ constant_expression
	    {
		terminal -> Kd = $2;
	    }
	| KP_EQ constant_expression
	    {
		terminal -> Kp = $2;
	    }

	| error
	;

layout_segment_connector_list
	: layout_segment_connector_list layout_segment_connector
	| /* empty */
	;

layout_segment_connector
	: layout_segment '{' segment_parameter_list '}' layout_connector layout_branch_list
	;

layout_branch_list
	: layout_branch_list layout_branch
	| /* empty */
	;

branch_segment_connector_list
	: branch_segment_connector_list branch_segment_connector
	| /* empty */
	;

branch_segment_connector
	: layout_segment '{' segment_parameter_list '}' layout_connector 
	;
	
layout_branch
        : branch_start '{' branch_segment_connector_list layout_terminal '}'
	;

branch_start
	: BRANCH_EQ
	    {
                num_branches ++;
		num_branch_segment = 0;

                branch = CreateBranch (num_branches);
                TreeInsert (ctl.problem -> branch_tree, branch);

                branch -> segment_from = segment_from;
                segment_from -> num_branch_to ++;
                segment_from -> branch_to [segment_from -> num_branch_to] = branch;
            }
	| error
        ;

layout_segment
	: SEGMENT_EQ
	    {
                 num_segs ++;
		         segment = CreateSegment (num_segs);

                 if (!branch) {
                    segment_from = segment;
                    TreeInsert (ctl.problem -> segment_tree, segment);
                 }
                 else {
                    TreeInsert (ctl.problem -> branch_segment_tree, segment);
                     
                    num_branch_segment ++;
		            branch -> segment [num_branch_segment] = segment;

                    segment -> branch = branch;
		         }

	    }
	| error
	;


segment_parameter_list
	: segment_parameter_list segment_parameter
	| /* empty */
	;

segment_parameter
	: LENGTH_EQ constant_expression
	    {
		segment -> length = $2;
	    }

	| MATERIAL_EQ NAME
	    {
		segment -> material = (Material) $2;
	    }

	| NAME_EQ NAME
	    {
        Deallocate(segment -> name);
		segment -> name = $2;
	    }

	| NODES_EQ distribution_pair_list
	    {
		    unsigned i;
 		    unsigned size;
            unsigned count;

		    size = dist_ptr - dist_array;

            segment -> num_dist = size;
 
            count = 0;
            for (i = 1 ; i <= size ; i++) {
               segment -> dist [i].nodes = dist_array [i - 1].nodes;
               segment -> dist [i].percent = dist_array [i - 1].percent;

               count += dist_array [i - 1].nodes;
		    }
            if (count <= 1) {
               error ("segment must have at least 2 nodes");
               break;
            } 
	    }

	| ATTACHMENTS_EQ attachment_list
	    {
		unsigned	i;
 		unsigned	size;

		size = attach_ptr - attach_array;

		segment -> num_attach = size;
		segment -> attach = Allocate(attachment, size);
                UnitOffset (segment -> attach);

                for (i = 1 ; i <= size ; i++) 
      	           segment -> attach [i] = attach_array [i - 1];
	    }

	| TOP_PAY_EQ  enable_copy variable_expression
	    {
                AssignSegmentPay (segment, 1, InCore, copy_input (0));  
	    }

	| BOTTOM_PAY_EQ  enable_copy variable_expression
	    {
                AssignSegmentPay (segment, 0, InCore, copy_input (0));  
	    }

	| TOP_LENGTH_EQ constant_expression
	    {
		        segment -> top_length = $2;
	    }

	| BOTTOM_LENGTH_EQ constant_expression
	    {
		        segment -> bottom_length = $2;
	    }

	| TOP_NODES_EQ distribution_pair_list
	    {
		unsigned i;
 		unsigned size;
                unsigned count;

		size = dist_ptr - dist_array;

                segment -> num_top_dist = size;
 
                count = 0;
                for (i = 1 ; i <= size ; i++) {
                   segment -> top_dist [i].nodes = dist_array [i - 1].nodes;
                   segment -> top_dist [i].percent = dist_array [i - 1].percent;

                   count += dist_array [i - 1].nodes;
		        }
                if (count <= 1) {
                   error ("top spool of segment must have at least 2 nodes");
                   break;
                } 
	    }

	| BOTTOM_NODES_EQ distribution_pair_list
	    {
		        unsigned i;
 		        unsigned size;
                unsigned count;

		        size = dist_ptr - dist_array;

                segment -> num_bottom_dist = size;
 
                count = 0;
                for (i = 1 ; i <= size ; i++) {
                   segment -> bottom_dist [i].nodes = dist_array [i - 1].nodes;
                   segment -> bottom_dist [i].percent = dist_array [i - 1].percent;

                   count += dist_array [i - 1].nodes;
		        }
                if (count <= 1) {
                   error ("bottom spool of segment must have at least 2 nodes");
                   break;
                } 
	    }

	| error
        ;
   
layout_connector
	: CONNECTOR_EQ NAME
	    {
		    segment -> connector = (Connector) $2;	
	    }
    | CONNECTOR_EQ CONNECTION_TYPE
        {
            segment -> connection = $2;
            segment -> connector = NULL;
        }
    | CONNECTOR_EQ '{' segment_connector_parameter_list '}'
    ;
	| /* empty */
	;

segment_connector_parameter_list
	: segment_connector_parameter_list segment_connector_parameter
	| /* empty */
	;

segment_connector_parameter
    : BODY_EQ NAME
        {
            segment -> connector = (Connector) $2;
        }     
	| XTHRUST_EQ  enable_copy variable_expression
            {
		        /* user X is internal Y */
                AssignConnectorThrust ('y', segment, InCore, copy_input (0));
            }

	| YTHRUST_EQ  enable_copy variable_expression
            {
		        /* user Y is internal Z */
                AssignConnectorThrust ('z', segment, InCore, copy_input (0));  
            }

	| ZTHRUST_EQ  enable_copy variable_expression
            {
		        /* user Z is internal X */
                AssignConnectorThrust ('x', segment, InCore, copy_input (0));  
            }
    | X_EQ constant_expression
            {
                segment -> connector_y = $2;
            }
    | Y_EQ constant_expression
            {
                segment -> connector_z = $2;
            }
    | Z_EQ constant_expression
            {
                segment -> connector_x = $2;
            }
    | INITIAL_XFORCE_EQ constant_expression
            {
                segment -> connector_yforce = $2;
            }
    | INITIAL_YFORCE_EQ constant_expression
            {
                segment -> connector_zforce = $2;
            }
    | INITIAL_ZFORCE_EQ constant_expression
            {
                segment -> connector_xforce = $2;
            }
    ;
     

attachment_list
	: attachment_list ',' attachment
	    {
		*attach_ptr ++ = $3;
            }

	| attachment_list attachment
	    {
		*attach_ptr ++ = $2;
	    }	

	| attachment
	    {
                attach_ptr = attach_array;
		*attach_ptr ++ = $1;
            }
	;

attachment	
	: NAME ':' '(' attachment_node_list ')'
	    {
		unsigned	i;
 		unsigned	size;

		size = int_ptr - int_array;

                $$.num_nodes = size;

                if (!($$.nodes = Allocate(int, size)))
                   Fatal ("unable to allocate memory for attachment nodes");

                UnitOffset ($$.nodes);

                for (i = 1 ; i <= size ; i++)
                   $$.nodes [i] = int_array [i - 1];

		$$.object = (Connector) $1;
	    }
        | NAME ':' '[' attachment_node_array ']'
            {
                unsigned        i;
                unsigned        size;

                size = (attach_stop - attach_start)/attach_step + 1;
                $$.num_nodes = size;

                if (!($$.nodes = Allocate(int, size)))
                   Fatal ("unable to allocate memory for attachment nodes");

                UnitOffset ($$.nodes);

                for (i = 0 ; i < size ; i++)
                   $$.nodes [i+1] = attach_start + i*attach_step;

                $$.object = (Connector) $1;
            }
        ;

attachment_node_array
        : INTEGER ',' INTEGER ',' INTEGER
            {
                attach_start = $1;
                attach_step  = $3;
                attach_stop  = $5;
            }
        ;

attachment_node_list
	: attachment_node_list ',' INTEGER
	    {
		*int_ptr ++ = $3;
	    }

	| attachment_node_list INTEGER
	    {
		*int_ptr ++ = $2;
	    }

	| /* empty */
	    {
		int_ptr = int_array;
	    }
	;


distribution_pair_list
        : distribution_pair_list ',' distribution_pair
            {
                *dist_ptr ++ = $3;
            }

        | distribution_pair_list distribution_pair
            {
                *dist_ptr ++ = $2;
            }

        | distribution_pair
            {
                dist_ptr = dist_array;
                *dist_ptr ++ = $1;
            }
        ;

distribution_pair
        : '(' INTEGER ',' constant_expression ')'
            {
                $$.nodes   = $2;
                $$.percent = $4;
            }
        ;


/* Environment section */

environment_section
	: ENVIRONMENT environment_parameter_list
	;


environment_parameter_list
	: environment_parameter_list environment_parameter
	| /* empty */
	;


environment_parameter
	: INPUT_TYPE_EQ INPUT_TYPE
	    {
		ctl.environment -> input_type = $2;
	    }

	| FORCING_METHOD_EQ FORCING_METHOD
	    {
		ctl.environment -> forcing = $2;
	    }

        | RHO_EQ constant_expression
	    {
		ctl.environment -> rho = $2;
	    }

        | DEPTH_EQ constant_expression
            {
                ctl.environment -> depth = $2;
            }
 
        | GRAVITY_EQ constant_expression
            {
                ctl.environment -> gravity = $2;
            }
 
        | FILL_FLUID_DENSITY_EQ constant_expression
            {
                ctl.environment -> fill_density = $2;
            }
 
        | BOTTOM_STIFFNESS_EQ constant_expression
	    {
		ctl.environment -> bottom_stiffness = $2;
 	    }

        | BOTTOM_FRICTION_EQ constant_expression
	    {
		ctl.environment -> bottom_friction = $2;
 	    }

 	| BOTTOM_DAMPING_EQ constant_expression
	    {
		ctl.environment -> bottom_damping = $2;
	    }

 	| BOTTOM_ELEVATION_EQ  enable_copy variable_expression
	    {
 		AssignElevation(ctl.environment, InCore, copy_input(0));
	    }

        | VELOCITY_FILE_EQ NAME
	    {
	 	Deallocate (ctl.environment -> velocity_file);
    		ctl.environment -> velocity_file = $2;
	    }

    | CURRENT_SCALE_EQ constant_expression
        {
        ctl.environment -> Uscale = $2;
        }

    | CURRENT_ROTATION_EQ constant_expression
        {
        ctl.environment -> current_rotation = $2;
        }

    | CURRENT_X_SCALE_EQ constant_expression
        {
 		/* user X is internal Y */
        ctl.environment -> Uyscale = $2;
        }

    | CURRENT_Y_SCALE_EQ constant_expression
        {
 		/* user Y is internal Z */
        ctl.environment -> Uzscale = $2;
        }

    | CURRENT_Z_SCALE_EQ constant_expression
        {
 		/* user Z is internal X */
        ctl.environment -> Uxscale = $2;
        }

	| CURRENT_FILE_EQ NAME
	    {
		Deallocate(ctl.environment -> current_file);
		ctl.environment -> current_file = $2;
	    }

        | WAVE_FILE_EQ NAME
	    {
	 	Deallocate (ctl.environment -> wave_file);
    		ctl.environment -> wave_file = $2;
	    }

	| XWAVE_EQ wave_triple_list		/* user X is internal Y */
	    {
                unsigned i;
                unsigned size;

                size = wave_ptr - wave_array;

                ctl.environment -> num_components [2] = size;
 
                for (i = 0 ; i < size ; i++) {
                   ctl.environment -> amplitude [2][i] = wave_array [i].amplitude;
                   ctl.environment -> period [2][i] = wave_array [i].period;
                   ctl.environment -> phase [2][i] = wave_array [i].phase;
                }
	    }

	| YWAVE_EQ  wave_triple_list		/* user Y is internal Z */
	    {
                unsigned i;
                unsigned size;

                size = wave_ptr - wave_array;

		ctl.environment -> num_components [3] = size;
 
                for (i = 0 ; i < size ; i++) {
                   ctl.environment -> amplitude [3][i] = wave_array [i].amplitude;
                   ctl.environment -> period [3][i] = wave_array [i].period;
                   ctl.environment -> phase [3][i] = wave_array [i].phase;
                }
	    }

	| XINPUT_EQ wave_triple_list		/* user X is internal Y */
	    {
                unsigned i;
                unsigned size;

                size = wave_ptr - wave_array;

                ctl.environment -> num_components [2] = size;

 
                for (i = 0 ; i < size ; i++) {
                   ctl.environment -> amplitude [2][i] = wave_array [i].amplitude;
                   ctl.environment -> period [2][i] = wave_array [i].period;
                   ctl.environment -> phase [2][i] = wave_array [i].phase;
                }
	    }

	| YINPUT_EQ  wave_triple_list	/* user Y is internal Z */
	    {
                unsigned i;
                unsigned size;

                size = wave_ptr - wave_array;

                ctl.environment -> num_components [3] = size;
 
                for (i = 0 ; i < size ; i++) {
                   ctl.environment -> amplitude [3][i] = wave_array [i].amplitude;
                   ctl.environment -> period [3][i] = wave_array [i].period;
                   ctl.environment -> phase [3][i] = wave_array [i].phase;
                }
	    }

	| ZINPUT_EQ  wave_triple_list	/* user Z is internal X */
	    {
                unsigned i;
                unsigned size;

                size = wave_ptr - wave_array;

                ctl.environment -> num_components [1] = size;
 
                for (i = 0 ; i < size ; i++) {
                   ctl.environment -> amplitude [1][i] = wave_array [i].amplitude;
                   ctl.environment -> period [1][i] = wave_array [i].period;
                   ctl.environment -> phase [1][i] = wave_array [i].phase;
                }
	    }

        | YWIND_EQ  enable_copy variable_expression
            {
 		    /* user Y is internal Z */
		    AssignWind (ctl.environment, 'z', InCore, copy_input(0));
            }

        | XWIND_EQ  enable_copy variable_expression
            {
 		    /* user X is internal y */
		    AssignWind (ctl.environment, 'y', InCore, copy_input(0));
            }

	| ZCURRENT_EQ  enable_copy named_variable_expression
	    {
		        /* user Z is internal X */
                ctl.environment -> Ux.name = varExName;
		        AssignCurrent (ctl.environment, 'x', InCore, copy_input (0));
	    }

	| YCURRENT_EQ  enable_copy named_variable_expression
	    {
		        /* user Y is internal Z */
                ctl.environment -> Uz.name = varExName;
		        AssignCurrent(ctl.environment, 'z', InCore, copy_input (0));
	    }

	| XCURRENT_EQ  enable_copy named_variable_expression
	    {
		        /* user X is internal Y */
                ctl.environment -> Uy.name = varExName;
		        AssignCurrent(ctl.environment, 'y', InCore, copy_input (0));
	    }

	| ZCURRENT_MODULATION_EQ  enable_copy variable_expression
	    {
		/* user Z is internal X */
		AssignCurrentModulation (ctl.environment, 'x', InCore, copy_input (0));
	    }

	| YCURRENT_MODULATION_EQ  enable_copy variable_expression
	    {
		/* user Y is internal Z */
		AssignCurrentModulation (ctl.environment, 'z', InCore, copy_input (0));
	    }

	| XCURRENT_MODULATION_EQ  enable_copy variable_expression
	    {
		/* user X is internal Y */
		AssignCurrentModulation (ctl.environment, 'y', InCore, copy_input (0));
	    }

	| error
	;


wave_triple_list
        : wave_triple_list ',' wave_triple 
            {
                *wave_ptr ++ = $3;
            }

        | wave_triple_list wave_triple
            {
                *wave_ptr ++ = $2;
            }

        | wave_triple
            {
                wave_ptr = wave_array;
                *wave_ptr ++ = $1;
            }
        ;

wave_triple
        : '(' constant_expression ',' constant_expression ',' constant_expression ')'
            {
                $$.amplitude = $2;
                $$.period    = $4;
                $$.phase     = $6;
            }
        ;


/* Analysis section */

analysis_section
	: ANALYSIS analysis_parameter_list
	;


analysis_parameter_list
	: analysis_parameter_list analysis_parameter
	| /* empty */
	;


analysis_parameter
	: TOLERANCE_EQ constant_expression
    {
	    ctl.analysis -> tolerance = $2;
    }

    | STATIC_TOLERANCE_EQ constant_expression
    {
	    ctl.analysis -> static_tolerance = $2;
    }

    | STATIC_OUTER_TOLERANCE_EQ constant_expression
    {
	    ctl.analysis -> outer_tolerance = $2;
    }

    | DYNAMIC_TOLERANCE_EQ constant_expression
    {
	    ctl.analysis -> dynamic_tolerance = $2;
	}

    | DYNAMIC_VAR_SMOOTHING_EQ enable_copy variable_expression
    {
        AssignSmoothing(ctl.analysis, InCore, copy_input(0));
    }

	| CURRENT_STEPS_EQ INTEGER
    {
		ctl.analysis -> current_steps = $2;
    }

	| MAX_ITERATIONS_EQ INTEGER
    {
		ctl.analysis -> maxit = $2;
    }

    | STATIC_ITERATIONS_EQ INTEGER
    {
		ctl.analysis -> static_it = $2;
    }

        | SHOOTING_ITERATIONS_EQ INTEGER
 	    {
		ctl.analysis -> shooting_it = $2;
	    }

        | STATIC_OUTER_ITERATIONS_EQ INTEGER
 	    {
		ctl.analysis -> outer_it = $2;
	    }

    | STATIC_OUTER_METHOD_EQ OUTER_METHOD
        {
        ctl.analysis -> static_outer_method = $2;
        }

	| STATIC_SOLUTION_EQ SOLUTION_METHOD
	    {
		ctl.analysis -> static_solution = $2;
	    }

	| STATIC_INITIAL_GUESS_EQ SOLUTION_METHOD
	    {
		ctl.analysis -> static_initial_guess = $2;
	    }

        | DYNAMIC_ITERATIONS_EQ INTEGER
 	    {
		ctl.analysis -> dynamic_it = $2;
	    }

	| RELAX_ADAPT_UP_EQ constant_expression
	    {
		ctl.analysis -> relax_up = $2;
	    }

	| RELAX_ADAPT_DOWN_EQ constant_expression
	    {
		ctl.analysis -> relax_down = $2;
	    }

	| RELAX_STALL_LIMIT_EQ INTEGER
	    {
		ctl.analysis -> stall_limit = $2;
	    }

	| MESH_SMOOTH_EQ constant_expression
	    {
		ctl.analysis -> mesh_smooth = $2;
       	    }

	| MESH_AMPLIFICATION_EQ constant_expression
	    {
		ctl.analysis -> mesh_amplify = $2;
	    }

	| ALPHA_K_EQ constant_expression
	    {
		ctl.analysis -> alpha_k = $2;
	    }

	| ALPHA_M_EQ constant_expression
	    {
		ctl.analysis -> alpha_m = $2;
	    }

	| DYN_RHO_EQ constant_expression
	    {
		ctl.analysis -> sp_rho = $2;
	    }

	| GAMMA_EQ constant_expression
	    {
		ctl.analysis -> gamma = $2;
	    }

	| RELAXATION_EQ constant_expression
	    {
		ctl.analysis -> relaxation = $2;
	    }

	| DYNAMIC_RELAXATION_EQ constant_expression
	    {
		ctl.analysis -> dynamic_relaxation = $2;
	    }

	| STATIC_RELAXATION_EQ constant_expression
	    {
		ctl.analysis -> static_relaxation = $2;
	    }

	| OUTER_RELAXATION_EQ constant_expression
	    {
		ctl.analysis -> outer_relaxation = $2;
	    }

	| DURATION_EQ constant_expression
	    {
		ctl.analysis -> duration = $2;
	    }

	| TIME_STEP_EQ constant_expression
	    {
		ctl.analysis -> dt = $2;
	    }

	| TIME_STEP_ADAPT_EQ INTEGER
	    {
		ctl.analysis -> adapt_factor = $2;
	    }

	| TIME_STEP_LEVELS_EQ INTEGER
	    {
		ctl.analysis -> adapt_levels = $2;
	    }

	| DYNAMIC_INTEGRATION_EQ INTEGRATION_METHOD
	    {
		ctl.analysis -> integration = $2;
	    }

	| RAMP_TIME_EQ constant_expression
	    {
		ctl.analysis -> ramp_time = $2;
	    }

	| VIVA_TIME_STEP_EQ constant_expression
	    {
		ctl.analysis -> viva_dt = $2;
	    }

	| VIVA_ITERATIONS_EQ INTEGER
	    {
		ctl.analysis -> viva_iterations = $2;
	    }


	| VIVA_DECIMATE_EQ INTEGER
	    {
		ctl.analysis -> viva_decimate = $2;
	    }

	| error
	;


/* Expressions */

named_variable_expression
    : variable_expression
        {
            varExName = NULL;
        }

    | '{' NAME '}' variable_expression
        {
            varExName = $2;
        }
    ;

variable_expression
	: expression
	    {
		EmitCode (HaltOp);
		SetIP (0);
	    }

	| discrete_pair_list
	    {
		if (table_error)
		    EmitCode (PushOp, 0.0);
		else
		    EmitCode (TableOp, table, table_count);

		EmitCode (HaltOp);
		table_count = 0;
		table_error = 0;
		last_depth = 0;
		SetIP (0);
	    }

	| discrete_pair_list '+'
	    {
		if (table_error)
		    EmitCode (PushOp, 0.0);
		else
		    EmitCode (CycleOp, table, table_count);

		EmitCode (HaltOp);
		table_count = 0;
		table_error = 0;
		last_depth = 0;
		SetIP (0);
	    }
	;


discrete_pair_list
	: discrete_pair_list ',' discrete_pair
	| discrete_pair_list discrete_pair
	| discrete_pair
	;


discrete_pair
	: '(' constant_expression ',' constant_expression ')'
	    {
		if ($2 < last_depth) {
		    error ("point not in nondecreasing order");
		    table_error = 1;
		    break;
		}

		if (table_count == table_size) {
		    table_size = table_size ? table_size << 1 : 8;
		    if (!Reallocate (table, double, table_size))
			Fatal ("unable to expand table");
		}

		table [table_count ++] = last_depth = $2;
		table [table_count ++] = $4;
	    }
	;


enable_copy
        : /* empty */
            {
                copy_input (1);
            }
        ;


constant_expression
	: expression
	    {
		EmitCode (HaltOp);
		SetIP (0);
		$$ = EvalCode (InCore, NULL, 0.0, 0.0, 0.0, 0.0, 0.0, CURRNODEDATA);
	    }
	;


expression
	: expression '?' if_action expression ':' else_action expression
	    {
		int ip = GetIP ( );
		SetIP (ip - $7 - 2);
		EmitCode (JmpOp, $7);
		SetIP (GetIP ( ) - $4 - 4);
		EmitCode (JzOp, $4 + 2);
		SetIP (ip);
		$$ = $1 + $3 + $4 + $6 + $7;
	    }

	| expression OR or_action expression
	    {
		int ip = GetIP ( );
		SetIP (ip - $4 - 3);
		EmitCode (JnzOp, $4 + 1);
		SetIP (ip);
		EmitCode (TestOp);
		$$ = $1 + $3 + $4 + 1;
	    }

	| expression AND and_action expression
	    {
		int ip = GetIP ( );
		SetIP (ip - $4 - 3);
		EmitCode (JzOp, $4 + 1);
		SetIP (ip);
		EmitCode (TestOp);
		$$ = $1 + $3 + $4 + 1;
	    }

	| expression '|' expression
	    {
		EmitCode (OrOp);
		$$ = $1 + 1 + $3;
	    }

	| expression '^' expression
	    {
		EmitCode (XorOp);
		$$ = $1 + 1 + $3;
	    }

	| expression '&' expression
	    {
		EmitCode (AndOp);
		$$ = $1 + 1 + $3;
	    }

	| expression EQUALS expression
	    {
		EmitCode (EqOp);
		$$ = $1 + 1 + $3;
	    }

	| expression NEQUAL expression
	    {
		    EmitCode (NeqOp);
		    $$ = $1 + 1 + $3;
	    }

	| expression '<' expression
	    {
		    EmitCode (LtOp);
		    $$ = $1 + 1 + $3;
	    }

	| expression '>' expression
	    {
		    EmitCode (GtOp);
		    $$ = $1 + 1 + $3;
	    }

	| expression LT_EQ expression
	    {
		    EmitCode (LteqOp);
		    $$ = $1 + 1 + $3;
	    }

	| expression GT_EQ expression
	    {
		    EmitCode (GteqOp);
		    $$ = $1 + 1 + $3;
	    }

	| expression LSHIFT expression
	    {
		    EmitCode (LsftOp);
		    $$ = $1 + 1 + $3;
	    }

	| expression RSHIFT expression
	    {
		    EmitCode (RsftOp);
		    $$ = $1 + 1 + $3;
	    }

	| expression '+' expression
	    {
		    EmitCode (AddOp);
		    $$ = $1 + 1 + $3;
	    }

	| expression '-' expression
	    {
		    EmitCode (SubOp);
		    $$ = $1 + 1 + $3;
	    }

	| expression '*' expression
	    {
		    EmitCode (MulOp);
		    $$ = $1 + 1 + $3;
	    }

	| expression '/' expression
	    {
		    EmitCode (DivOp);
		    $$ = $1 + 1 + $3;
	    }

	| expression '%' expression
	    {
		    EmitCode (ModOp);
		    $$ = $1 + 1 + $3;
	    }

	| '+' expression		%prec UNARY
	    {
		    $$ = $2;
	    }

	| '-' expression		%prec UNARY
	    {
		    EmitCode (NegOp);
		    $$ = 1 + $2;
	    }

	| '!' expression
	    {
		    EmitCode (NotOp);
		    $$ = 1 + $2;
	    }

	| '~' expression
	    {
		    EmitCode (InvOp);
		    $$ = 1 + $2;
	    }

	| '(' expression ')'
	    {
		    $$ = $2;
	    }

	| INTEGER
	    {
		    EmitCode (PushOp, (double) $1);
		    $$ = 2;
	    }

	| DOUBLE
	    {
		    EmitCode (PushOp, $1);
		    $$ = 2;
	    }

    | CONST_PI
	    {
		    EmitCode (PushOp, M_PI);
	        $$ = 2;
	    }


	| function
	;

function
    : SIN '(' expression ')'
	    {
		EmitCode (SinOp);
		$$ = $3 + 1;
	    }

	| COS '(' expression ')'
	    {
		EmitCode (CosOp);
		$$ = $3 + 1;
	    }

	| TAN '(' expression ')'
	    {
		EmitCode (TanOp);
		$$ = $3 + 1;
	    }

	| TANH '(' expression ')'
	    {
		EmitCode (TanhOp);
		$$ = $3 + 1;
	    }

	| SINH '(' expression ')'
	    {
		EmitCode (SinhOp);
		$$ = $3 + 1;
	    }

	| COSH '(' expression ')'
	    {
		EmitCode (CoshOp);
		$$ = $3 + 1;
	    }

	| ASIN '(' expression ')'
	    {
		EmitCode (AsinOp);
		$$ = $3 + 1;
	    }

	| ACOS '(' expression ')'
	    {
		EmitCode (AcosOp);
		$$ = $3 + 1;
	    }

    | ATAN2 '(' expression ',' expression ')'
        {
        EmitCode(Atan2Op);
        $$ = $3 + $5 + 1;
        }

	| POW '(' expression ',' expression ')'
	    {
		EmitCode (PowOp);
		$$ = $3 + $5 + 1;
	    }

	| EXP '(' expression ')'
	    {
		EmitCode (ExpOp);
		$$ = $3 + 1;
	    }

	| LOG '(' expression ')'
	    {
		EmitCode (LnOp);
		$$ = $3 + 1;
	    }

	| LOG10 '(' expression ')'
	    {
		EmitCode (LogOp);
		$$ = $3 + 1;
	    }

	| SQRT '(' expression ')'
	    {
		EmitCode (SqrtOp);
		$$ = $3 + 1;
	    }

	| HYPOT '(' expression ',' expression ')'
	    {
		EmitCode (HypotOp);
		$$ = $3 + $5 + 1;
	    }

	| FLOOR '(' expression ')'
	    {
		EmitCode (FloorOp);
		$$ = $3 + 1;
	    }

	| CEIL '(' expression ')'
	    {
		EmitCode (CeilOp);
		$$ = $3 + 1;
	    }

	| FMOD '(' expression ',' expression ')'
	    {
		EmitCode (FmodOp);
		$$ = $3 + $5 + 1;
	    }

	| FABS '(' expression ')'
	    {
		EmitCode (FabsOp);
		$$ = $3 + 1;
	    }
    | VAR_t 
        {
                EmitCode(VartOp);
		        $$ = 1;
        }

    | VAR_n
        {
                EmitCode(VarNOp);
                $$ = 1;
        }
    | VAR_s
        {       EmitCode(PushOp, 0.0);
                EmitCode(VarSOp);
                $$ = 3;
        }
    | VAR_s '(' expression ')'
        {
                EmitCode(VarSOp);
		        $$ = $3 + 1;
        }

    | VAR_e '(' expression ')'
        {
                EmitCode(VarEOp);
		        $$ = $3 + 1;
        }
    | VAR_e 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarEOp);
		        $$ = 3;
        }

    | VAR_e_f '(' expression ')'
        {
                EmitCode(VarEOp | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_e_f
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarEOp | FILTERNODEDATA);
		        $$ = 3;
        }

    | VAR_Sn '(' expression ')'
        {
                EmitCode(VarSnOp);
		        $$ = $3 + 1;
        }
    | VAR_Sn 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarSnOp);
		        $$ = 3;
        }
    | VAR_Sn_f '(' expression ')'
        {
                EmitCode(VarSnOp | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_Sn_f 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarSnOp | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_Sb '(' expression ')'
        {
                EmitCode(VarSbOp);
		        $$ = $3 + 1;
        }
    | VAR_Sb 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarSbOp);
		        $$ = 3;
        }
    | VAR_Sb_f '(' expression ')'
        {
                EmitCode(VarSbOp | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_Sb_f 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarSbOp | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_u '(' expression ')'
        {
                EmitCode(VarUOp);
		        $$ = $3 + 1;
        }
    | VAR_u 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarUOp);
		        $$ = 3;
        }
    | VAR_u_f '(' expression ')'
        {
                EmitCode(VarUOp | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_u_f 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarUOp | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_v '(' expression ')'
        {
                EmitCode(VarVOp);
		        $$ = $3 + 1;
        }
    | VAR_v 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarVOp);
		        $$ = 3;
        }
    | VAR_v_f '(' expression ')'
        {
                EmitCode(VarVOp | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_v_f 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarVOp | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_w '(' expression ')'
        {
                EmitCode(VarWOp);
		        $$ = $3 + 1;
        }
    | VAR_w 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarWOp);
		        $$ = 3;
        }
    | VAR_w_f '(' expression ')'
        {
                EmitCode(VarWOp | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_w_f 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarWOp | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_B0
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarB0Op);
		        $$ = 3;
        }
    | VAR_B0 '(' expression ')'
        {
                EmitCode(VarB0Op);
		        $$ = $3 + 1;
        }
    | VAR_B0_f
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarB0Op | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_B0_f '(' expression ')'
        {
                EmitCode(VarB0Op | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_B1
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarB1Op);
		        $$ = 3;
        }
    | VAR_B1 '(' expression ')'
        {
                EmitCode(VarB1Op);
		        $$ = $3 + 1;
        }
    | VAR_B1_f
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarB1Op | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_B1_f '(' expression ')'
        {
                EmitCode(VarB1Op | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_B2
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarB2Op);
		        $$ = 3;
        }
    | VAR_B2 '(' expression ')'
        {
                EmitCode(VarB2Op);
		        $$ = $3 + 1;
        }
    | VAR_B2_f
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarB2Op | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_B2_f '(' expression ')'
        {
                EmitCode(VarB2Op | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_B3
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarB3Op);
		        $$ = 3;
        }
    | VAR_B3 '(' expression ')'
        {
                EmitCode(VarB3Op);
		        $$ = $3 + 1;
        }
    | VAR_B3_f
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarB3Op | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_B3_f '(' expression ')'
        {
                EmitCode(VarB3Op | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_Om1
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarOm1Op);
		        $$ = 3;
        }
    | VAR_Om1 '(' expression ')'
        {
                EmitCode(VarOm1Op);
		        $$ = $3 + 1;
        }
    | VAR_Om1_f
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarOm1Op | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_Om1_f '(' expression ')'
        {
                EmitCode(VarOm1Op | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_Om2
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarOm2Op);
		        $$ = 3;
        }
    | VAR_Om2 '(' expression ')'
        {
                EmitCode(VarOm2Op);
		        $$ = $3 + 1;
        }
    | VAR_Om2_f
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarOm2Op | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_Om2_f '(' expression ')'
        {
                EmitCode(VarOm2Op | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_Om3
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarOm3Op);
		        $$ = 3;
        }
    | VAR_Om3 '(' expression ')'
        {
                EmitCode(VarOm3Op);
		        $$ = $3 + 1;
        }
    | VAR_Om3_f
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarOm3Op | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_Om3_f '(' expression ')'
        {
                EmitCode(VarOm3Op | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_phi
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarPhiOp);
		        $$ = 3;
        }
    | VAR_phi '(' expression ')'
        {
                EmitCode(VarPhiOp);
		        $$ = $3 + 1;
        }
    | VAR_phi_f
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarPhiOp | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_phi_f '(' expression ')'
        {
                EmitCode(VarPhiOp | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_T '(' expression ')'
        {
                EmitCode(VarTOp);
		        $$ = $3 + 1;
        }
    | VAR_T 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarTOp);
		        $$ = 3;
        }
    | VAR_T_f '(' expression ')'
        {
                EmitCode(VarTOp | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_T_f 
        {
                printf("VAR_T_f\n");
                EmitCode(PushOp, 0.0);
                EmitCode(VarTOp | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_x '(' expression ')'
        {
                EmitCode(VarXOp);
		        $$ = $3 + 1;
        }
    | VAR_x 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarXOp);
		        $$ = 3;
        }
    | VAR_x_f '(' expression ')'
        {
                EmitCode(VarXOp | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_x_f 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarXOp | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_y '(' expression ')'
        {
                EmitCode(VarYOp);
		        $$ = $3 + 1;
        }
    | VAR_y 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarYOp);
		        $$ = 3;
        }
    | VAR_y_f '(' expression ')'
        {
                EmitCode(VarYOp | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_y_f 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarYOp | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_z '(' expression ')'
        {
                EmitCode(VarZOp);
		        $$ = $3 + 1;
        }
    | VAR_z 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarZOp);
		        $$ = 3;
        }
    | VAR_z_f '(' expression ')'
        {
                EmitCode(VarZOp | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_z_f 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarZOp | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_Ux '(' expression ')'
        {
                EmitCode(VarUxOp);
		        $$ = $3 + 1;
        }
    | VAR_Ux 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarUxOp);
		        $$ = 3;
        }
    | VAR_Ux_f '(' expression ')'
        {
                EmitCode(VarUxOp | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_Ux_f 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarUxOp | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_Vy '(' expression ')'
        {
                EmitCode(VarVyOp);
		        $$ = $3 + 1;
        }
    | VAR_Vy 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarVyOp);
		        $$ = 3;
        }
    | VAR_Vy_f '(' expression ')'
        {
                EmitCode(VarVyOp | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_Vy_f 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarVyOp | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_Wz '(' expression ')'
        {
                EmitCode(VarWzOp);
		        $$ = $3 + 1;
        }
    | VAR_Wz 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarWzOp);
		        $$ = 3;
        }
    | VAR_Wz_f '(' expression ')'
        {
                EmitCode(VarWzOp | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_Wz_f 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarWzOp | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_Fx '(' expression ')'
        {
                EmitCode(VarFxOp);
		        $$ = $3 + 1;
        }
    | VAR_Fx 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarFxOp);
		        $$ = 3;
        }
    | VAR_Fx_f '(' expression ')'
        {
                EmitCode(VarFxOp | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_Fx_f 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarFxOp | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_Fy '(' expression ')'
        {
                EmitCode(VarFyOp);
		        $$ = $3 + 1;
        }
    | VAR_Fy 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarFyOp);
		        $$ = 3;
        }
    | VAR_Fy_f '(' expression ')'
        {
                EmitCode(VarFyOp | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_Fy_f
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarFyOp | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_Fz '(' expression ')'
        {
                EmitCode(VarFzOp);
		        $$ = $3 + 1;
        }
    | VAR_Fz 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarFzOp);
		        $$ = 3;
        }
    | VAR_Fz_f '(' expression ')'
        {
                EmitCode(VarFzOp | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_Fz_f 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarFzOp | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_H '(' expression ')'
        {
                EmitCode(VarHOp);
		        $$ = $3 + 1;
        }
    | VAR_H 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarHOp);
		        $$ = 3;
        }
    | VAR_H_f '(' expression ')'
        {
                EmitCode(VarHOp | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_H_f 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarHOp | FILTERNODEDATA);
		        $$ = 3;
        }
    | VAR_p '(' expression ')'
        {
                EmitCode(VarPOp);
		        $$ = $3 + 1;
        }
    | VAR_p 
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarPOp);
		        $$ = 3;
        }
    | VAR_p_f '(' expression ')'
        {
                EmitCode(VarPOp | FILTERNODEDATA);
		        $$ = $3 + 1;
        }
    | VAR_p_f
        {
                EmitCode(PushOp, 0.0);
                EmitCode(VarPOp | FILTERNODEDATA);
		        $$ = 3;
        }
	| FIRST '(' expression ')'
        {
        EmitCode(FirstOp);
        $$ = $3 + 1;
        }
        
    | FIRSTACTIVE '(' expression ')'
        {
        EmitCode(FirstActiveOp);
        $$ = $3 + 1;
        }
        
	| LAST '(' expression ')'
        {
        EmitCode(LastOp);
        $$ = $3 + 1;
        }
        
    | LASTACTIVE '(' expression ')'
        {
        EmitCode(LastActiveOp);
        $$ = $3 + 1;
        }
    ;
        


if_action
	: /* empty */
	    {
		EmitCode (JzOp, 0);
		$$ = 2;
	    }
	;


else_action
	: /* empty */
	    {
		EmitCode (JmpOp, 0);
		$$ = 2;
	    }
	;


or_action
	: /* empty */
	    {
		EmitCode (CopyOp);
		EmitCode (JnzOp, 0);
		EmitCode (PopOp);
		$$ = 4;
	    }
	;


and_action
	: /* empty */
	    {
		EmitCode (CopyOp);
		EmitCode (JzOp, 0);
		EmitCode (PopOp);
		$$ = 4;
	    }
	;

%%
/*
# ifdef YYBYACC
char *buoy_suppress_warnings_from_gcc = yysccsid;
# endif
*/
