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
 * File:	lexer.l							*
 ************************************************************************/

# include <stdio.h>
# include <string.h>
# include <unistd.h>
# include "error.h"
# include "compress.h"
# include "problem.h"
# include "allocate.h"

extern ParserControl ctl;

#include "parser.h"

# undef atof
extern double atof ( );


static char *filename;
static void  update_line(void);

static int   copy_flag = 0;
static char  copied_text [65536];
static char *copy_ptr = copied_text;

# define YY_USER_ACTION \
	if (copy_flag) { \
	    char *ptr = yytext; \
	    while ((*copy_ptr = *ptr ++)) copy_ptr ++; \
	}

# undef  YY_INPUT
# define YY_INPUT(buf,result,max_size) \
	if (ctl.problem -> input) \
	    result = (*buf = *(ctl.problem -> input) ++) ? 1 : (ctl.problem -> input --, 0); \
	else \
	    if ((result = read (fileno (yyin), (char *) buf, max_size)) < 0) \
		YY_FATAL_ERROR ("read() in flex scanner failed");
%}

digit	[0-9]
exp	[eE][+-]?{digit}+
/* punct	[/.:] */
punct	[/.:]
letter	[a-zA-Z_]
/* id	{letter}({letter}|{digit}|{punct})* */
id	{letter}({letter}|{digit}|{punct})*
eq	[ \t]*=

%x	comment
%%



%{
/* Section headings */
%}

problem\ description	{return PROBLEM;}
analysis\ parameters	{return ANALYSIS;}
environment		{return ENVIRONMENT;}
buoys			{return BUOYS;}
platforms		{return BUOYS;}
anchors			{return ANCHORS;}
connectors		{return CONNECTORS;}
materials		{return MATERIALS;}
layout			{return LAYOUT;}
end			{return END;}


%{
/* Reserved words (in expression or right hand side of assignment) */
%}

pi				{return CONST_PI;}
t				{return VAR_t;}	

e				{return VAR_e;}
Sn              {return VAR_Sn;}
Sb              {return VAR_Sb;}
F  				{return VAR_T;}
u				{return VAR_u;}
v				{return VAR_v;}
w				{return VAR_w;}
B0              {return VAR_B0;}
B1              {return VAR_B1;}
B2              {return VAR_B2;}
B3              {return VAR_B3;}
Om1              {return VAR_Om1;}
Om2              {return VAR_Om2;}
Om3              {return VAR_Om3;}
phi             {return VAR_phi;}
x				{return VAR_x;}
y				{return VAR_y;}
z				{return VAR_z;}
ux				{return VAR_Ux;}
vy				{return VAR_Vy;}
wz				{return VAR_Wz;}
Fx              {return VAR_Fx;}
Fy              {return VAR_Fy;}
Fz              {return VAR_Fz;}
H				{return VAR_H;}
P               {return VAR_p;}
s               {return VAR_s;}
n               {return VAR_n;}

e_				{return VAR_e_f;}
Sn_             {return VAR_Sn_f;}
Sb_             {return VAR_Sb_f;}
F_  			{return VAR_T_f;}
u_			    {return VAR_u_f;}
v_				{return VAR_v_f;}
w_				{return VAR_w_f;}
B0_             {return VAR_B0_f;}
B1_             {return VAR_B1_f;}
B2_             {return VAR_B2_f;}
B3_             {return VAR_B3_f;}
Om1_            {return VAR_Om1_f;}
Om2_            {return VAR_Om2_f;}
Om3_            {return VAR_Om3_f;}
phi_            {return VAR_phi_f;}
x_				{return VAR_x_f;}
y_				{return VAR_y_f;}
z_				{return VAR_z_f;}
ux_				{return VAR_Ux_f;}
vy_				{return VAR_Vy_f;}
wz_				{return VAR_Wz_f;}
Fx_             {return VAR_Fx_f;}
Fy_             {return VAR_Fy_f;}
Fz_             {return VAR_Fz_f;}
H_				{return VAR_H_f;}
P_              {return VAR_p_f;}

first           {return FIRST;}
last            {return LAST;}
firstActive     {return FIRSTACTIVE;}
lastActive      {return LASTACTIVE;}
filter          {return FILTER;}

sin				{return SIN;}
cos				{return COS;}
tan				{return TAN;}
sinh				{return SINH;}
cosh				{return COSH;}
tanh				{return TANH;}
pow				{return POW;}
exp				{return EXP;}
log				{return LOG;}
log10				{return LOG10;}
sqrt				{return SQRT;}
hypot				{return HYPOT;}
floor				{return FLOOR;}
ceil				{return CEIL;}
fmod				{return FMOD;}
fabs				{return FABS;}
asin                {return ASIN;}
acos                {return ACOS;}
atan2               {return ATAN2;}

general				{Cable_lval.i = General;    return PROBLEM_TYPE;}
surface				{Cable_lval.i = Surface;    return PROBLEM_TYPE;}
subsurface			{Cable_lval.i = Subsurface; return PROBLEM_TYPE;}
horizontal			{Cable_lval.i = Horizontal; return PROBLEM_TYPE;}
horizontal-drifter  {Cable_lval.i = HorizontalDrifter; return PROBLEM_TYPE;}
towing				{Cable_lval.i = Towing;     return PROBLEM_TYPE;}
drifter				{Cable_lval.i = Drifter;    return PROBLEM_TYPE;}
deployment          {Cable_lval.i = Deployment; return PROBLEM_TYPE;}
riser				{Cable_lval.i = Riser;      return PROBLEM_TYPE;}
webster				{Cable_lval.i = Webster;    return PROBLEM_TYPE;}
tensioned			{Cable_lval.i = Webster;	return PROBLEM_TYPE;}
positioned			{Cable_lval.i = Horizontal;	return PROBLEM_TYPE;}

spatial			{Cable_lval.i = Spatial;      return INTEGRATION_METHOD;}
temporal		{Cable_lval.i = Temporal;     return INTEGRATION_METHOD;}

bisection       {Cable_lval.i = Bisection;    return OUTER_METHOD;}
secant          {Cable_lval.i = Secant;    return OUTER_METHOD;}
fixed          {Cable_lval.i = Fixed;    return OUTER_METHOD;}
variable          {Cable_lval.i = Variable;    return OUTER_METHOD;}

relaxation		{Cable_lval.i = Relaxation;	return SOLUTION_METHOD;}
catenary		{Cable_lval.i = Catenary; 	return SOLUTION_METHOD;}
shooting		{Cable_lval.i = Shooting;	return SOLUTION_METHOD;}

velocity		{Cable_lval.i = Velocity;     return FORCING_METHOD;}
lamp			{Cable_lval.i = LAMP;     	  return FORCING_METHOD;}
morison			{Cable_lval.i = Morison;      return FORCING_METHOD;}
froude-krylov		{Cable_lval.i = FroudeKrylov; return FORCING_METHOD;}
wave-follower		{Cable_lval.i = WaveFollower; return FORCING_METHOD;}
spar			{Cable_lval.i = Spar;         return FORCING_METHOD;}
force			{Cable_lval.i = Force;	  return FORCING_METHOD;}

regular			{Cable_lval.i = Regular; return INPUT_TYPE;}
random			{Cable_lval.i = Random; return INPUT_TYPE;}

sphere          		{Cable_lval.i = Sphere; return BUOY_TYPE;}
cylinder          		{Cable_lval.i = Cylinder; return BUOY_TYPE;}
capsule          		{Cable_lval.i = Capsule; return BUOY_TYPE;}
ship  	        		{Cable_lval.i = Ship; return BUOY_TYPE;}
platform         		{Cable_lval.i = Platform; return BUOY_TYPE;}
axisymmetric			{Cable_lval.i = Axisymmetric; return BUOY_TYPE;}

linear				{Cable_lval.i = Linear; return MATERIAL_TYPE;}
nonlinear			{Cable_lval.i = Nonlinear; return MATERIAL_TYPE;}

spliced             {Cable_lval.i = Spliced; return CONNECTION_TYPE;}
pinned              {Cable_lval.i = Pinned; return CONNECTION_TYPE;}

%{
/* Assignments */
%}

title{eq}			{return TITLE_EQ;}
color{eq}			{return COLOR_EQ;}

buoyancy{eq}			{return BUOYANCY_EQ;}
cg{eq}				{return CG_EQ;}
h{eq}				{return H_EQ;}
height{eq}          {return H_EQ;}
diameters{eq}			{return DIAMETERS_EQ;}
draft{eq}			{return DRAFT_EQ;}
type{eq}			{return TYPE_EQ;}
cl{eq}			        {return CL_EQ;}
cdt{eq}				{return CDT_EQ;}
cdn{eq}				{return CDN_EQ;}
cdw{eq}				{return CDW_EQ;}
sw{eq}				{return SW_EQ;}
m{eq}				{return M_EQ;}
mass{eq}            {return M_EQ;}
mu{eq}				{return MU_EQ;}
safety{eq}          {return SAFETY_EQ;}
tension{eq}			{return TENSION_EQ;}
initial-x-force{eq}			{return INITIAL_XFORCE_EQ;}
initial-y-force{eq}			{return INITIAL_YFORCE_EQ;}
initial-z-force{eq}			{return INITIAL_ZFORCE_EQ;}
x-force{eq}			{return XFORCE_EQ;}
y-force{eq}			{return YFORCE_EQ;}
z-force{eq}			{return ZFORCE_EQ;}
x-speed{eq}			{return XSPEED_EQ;}
y-speed{eq}			{return YSPEED_EQ;}
z-speed{eq}			{return ZSPEED_EQ;}
x-thrust{eq}			{return XTHRUST_EQ;}
y-thrust{eq}			{return YTHRUST_EQ;}
z-thrust{eq}			{return ZTHRUST_EQ;}
profile{eq}			{return PROFILE_EQ;}
profile-Ki{eq}			{return KI_EQ;}
profile-Kd{eq}			{return KD_EQ;}
profile-Kp{eq}			{return KP_EQ;}
profile-slope{eq}		{return PROFILE_SLOPE_EQ;}
profile-range{eq}		{return PROFILE_RANGE_EQ;}
profile-turn-depth{eq}		{return PROFILE_TURN_DEPTH_EQ;}
flap-file{eq}			{return FLAP_FILE_EQ;}
release-time{eq}		{return RELEASE_TIME_EQ;}
buoyancy-head{eq}		{return BUOYANCY_HEAD_EQ;}
	
buoy{eq}			{return BUOY_EQ;}
platform{eq}			{return BUOY_EQ;}
anchor{eq}			{return ANCHOR_EQ;}
node{eq}			{return NODE_EQ;}
connector{eq}			{return CONNECTOR_EQ;}
body{eq}            {return BODY_EQ;}
segment{eq}			{return SEGMENT_EQ;}
nodes{eq}			{return NODES_EQ;}
material{eq}			{return MATERIAL_EQ;}
length{eq}			{return LENGTH_EQ;}
top-spool-pay{eq}		{return TOP_PAY_EQ;}
bottom-spool-pay{eq}		{return BOTTOM_PAY_EQ;}
top-spool-length{eq}		{return TOP_LENGTH_EQ;}
bottom-spool-length{eq}		{return BOTTOM_LENGTH_EQ;}
top-spool-nodes{eq}		{return TOP_NODES_EQ;}
bottom-spool-nodes{eq}		{return BOTTOM_NODES_EQ;}
attachments{eq}			{return ATTACHMENTS_EQ;}
terminal{eq}			{return TERMINAL_EQ;}
branch{eq}			{return BRANCH_EQ;}

category{eq}        {return CATEGORY_EQ;}
comment{eq}			{return COMMENT_EQ;}
ea{eq}				{return EA_EQ;}
ei{eq}				{return EI_EQ;}
gj{eq}				{return GJ_EQ;}
d{eq}				{return D_EQ;}
diam{eq}            {return D_EQ;}
b{eq}				{return B_EQ;}
bt{eq}				{return BT_EQ;}
bn{eq}				{return BN_EQ;}
a{eq}				{return A_EQ;}
w{eq}				{return W_EQ;}
wet{eq}				{return WET_EQ;}
am{eq}				{return AM_EQ;}
amn{eq}				{return AMN_EQ;}
amt{eq}				{return AMT_EQ;}
cmn{eq}				{return CMN_EQ;}
cmt{eq}				{return CMT_EQ;}
can{eq}				{return CAN_EQ;}
cat{eq}				{return CAT_EQ;}
x{eq}				{return X_EQ;}
y{eq}				{return Y_EQ;}
z{eq}				{return Z_EQ;}
t{eq}				{return T_EQ;}
te{eq}				{return TE_EQ;}
tee{eq}				{return TEE_EQ;}
e{eq}				{return E_EQ;}
i{eq}				{return I_EQ;}
j{eq}				{return J_EQ;}
g{eq}				{return G_EQ;}
id{eq}				{return ID_EQ;}
nu{eq}				{return NU_EQ;}
swl{eq}             {return SWL_EQ;}
yield{eq}           {return YIELD_EQ;}

max-iterations{eq}		{return MAX_ITERATIONS_EQ;}
static-iterations{eq}		{return STATIC_ITERATIONS_EQ;}
dynamic-iterations{eq}		{return DYNAMIC_ITERATIONS_EQ;}
static-outer-iterations{eq}	{return STATIC_OUTER_ITERATIONS_EQ;}
shooting-iterations{eq}		{return SHOOTING_ITERATIONS_EQ;}
static-outer-method{eq}     {return STATIC_OUTER_METHOD_EQ;}

dynamic-var-smoothing{eq} {return DYNAMIC_VAR_SMOOTHING_EQ;}

duration{eq}			{return DURATION_EQ;}
ramp-time{eq}			{return RAMP_TIME_EQ;}
dynamic-integration{eq}		{return DYNAMIC_INTEGRATION_EQ;}
time-step{eq}			{return TIME_STEP_EQ;}
time-step-adapt-factor{eq}	{return TIME_STEP_ADAPT_EQ;}
time-step-adapt-levels{eq}	{return TIME_STEP_LEVELS_EQ;}

viva-time-step{eq}      {return VIVA_TIME_STEP_EQ;}
viva-decimate{eq}      {return VIVA_DECIMATE_EQ;}
viva-iterations{eq}     {return VIVA_ITERATIONS_EQ;}

current-steps{eq}		{return CURRENT_STEPS_EQ;}
static-solution{eq}		{return STATIC_SOLUTION_EQ;}
static-initial-guess{eq}        {return STATIC_INITIAL_GUESS_EQ;}

tolerance{eq}			{return TOLERANCE_EQ;}
static-tolerance{eq}		{return STATIC_TOLERANCE_EQ;}
dynamic-tolerance{eq}		{return DYNAMIC_TOLERANCE_EQ;}
static-outer-tolerance{eq}	{return STATIC_OUTER_TOLERANCE_EQ;}

relaxation{eq}			{return RELAXATION_EQ;}
dynamic-relaxation{eq}		{return DYNAMIC_RELAXATION_EQ;}
static-relaxation{eq}		{return STATIC_RELAXATION_EQ;}
static-outer-relaxation{eq}	{return OUTER_RELAXATION_EQ;}

mesh-smoothing-length{eq}	{return MESH_SMOOTH_EQ;}
mesh-amplification{eq}		{return MESH_AMPLIFICATION_EQ;}

relax-adapt-up{eq}		{return RELAX_ADAPT_UP_EQ;}
relax-adapt-down{eq}	 	{return RELAX_ADAPT_DOWN_EQ;}
relax-stall-limit{eq}		{return RELAX_STALL_LIMIT_EQ;}

dynamic-alpha-k{eq}		{return ALPHA_K_EQ;}
dynamic-alpha-m{eq}		{return ALPHA_M_EQ;}
dynamic-gamma{eq}		{return GAMMA_EQ;}
dynamic-rho{eq}			{return DYN_RHO_EQ;}
dynamic-lambda{eq}		{return DYN_RHO_EQ;}

wave-file{eq}			{return WAVE_FILE_EQ;}
force-file{eq}			{return VELOCITY_FILE_EQ;}
velocity-file{eq}		{return VELOCITY_FILE_EQ;}
current-file{eq}		{return CURRENT_FILE_EQ;}
z-current{eq}			{return ZCURRENT_EQ;}
y-current{eq}			{return YCURRENT_EQ;}
x-current{eq}			{return XCURRENT_EQ;}
y-current-profile{eq}			{return YCURRENT_PROFILE_EQ;}
x-current-profile{eq}			{return XCURRENT_PROFILE_EQ;}
x-current-modulation{eq}	{return XCURRENT_MODULATION_EQ;}
y-current-modulation{eq}	{return YCURRENT_MODULATION_EQ;}
z-current-modulation{eq}	{return XCURRENT_MODULATION_EQ;}
current-scale{eq}           {return CURRENT_SCALE_EQ;}
current-rotation{eq}        {return CURRENT_ROTATION_EQ;}
x-current-scale{eq}           {return CURRENT_X_SCALE_EQ;}
y-current-scale{eq}           {return CURRENT_Y_SCALE_EQ;}
z-current-scale{eq}           {return CURRENT_Z_SCALE_EQ;}
forcing-method{eq}		{return FORCING_METHOD_EQ;}
input-type{eq}			{return INPUT_TYPE_EQ;}
x-input{eq}			{return XINPUT_EQ;}
y-input{eq}			{return YINPUT_EQ;}
z-input{eq}			{return ZINPUT_EQ;}
x-wave{eq}			{return XWAVE_EQ;}
y-wave{eq}			{return YWAVE_EQ;}
rho{eq}				{return RHO_EQ;}
density{eq}         {return RHO_EQ;}
depth{eq}			{return DEPTH_EQ;}
gravity{eq}			{return GRAVITY_EQ;}
x-wind{eq}			{return XWIND_EQ;}
y-wind{eq}			{return YWIND_EQ;}
bottom-stiffness{eq}		{return BOTTOM_STIFFNESS_EQ;}
bottom-damping{eq}		{return BOTTOM_DAMPING_EQ;}
bottom-elevation{eq}  		{return BOTTOM_ELEVATION_EQ;}
bottom-friction{eq}		{return BOTTOM_FRICTION_EQ;}
fill-fluid-density{eq}		{return FILL_FLUID_DENSITY_EQ;}


%{
/* Operators */
%}

"&&"				{return AND;}
"||"				{return OR;}
"=="				{return EQUALS;}
"!="				{return NEQUAL;}
"<="				{return LT_EQ;}
">="				{return GT_EQ;}
"<<"				{return LSHIFT;}
">>"				{return RSHIFT;}


%{
/* Numbers and names */
%}

{digit}+			{Cable_lval.i = atoi (yytext); return INTEGER;}
{digit}+{exp}			{Cable_lval.d = atof (yytext); return DOUBLE;}
{digit}+\.{digit}*{exp}?	{Cable_lval.d = atof (yytext); return DOUBLE;}
{digit}*\.{digit}+{exp}?	{Cable_lval.d = atof (yytext); return DOUBLE;}

{id}				{Cable_lval.s = strdup (yytext); return NAME;}
\"[^"]*\"			{Cable_lval.s = strdup (yytext + 1);
				     Cable_lval.s [yyleng - 2] = 0; return NAME;}
\'[^']*\'			{Cable_lval.s = strdup (yytext + 1);
				     Cable_lval.s [yyleng - 2] = 0; return NAME;}


%{
/* Comments, white space, and remaining preprocessor lines (the comment
   rules are straight from the flex documentation). */
%}

#[^\n]*$			{update_line ( );}
[ \t]+				{}
\n				{ctl.problem -> line ++;}
"/*"				{BEGIN comment;}

\r                      {}
<comment>[^*\n]*		{}
<comment>[^*\n]*\n		{ctl.problem -> line ++;}
<comment>"*"+[^*/\n]*		{}
<comment>"*"+[^*/\n]*\n		{ctl.problem -> line ++;}
<comment>"*"+"/"		{BEGIN INITIAL;}

.				{return *yytext;}

%%

/************************************************************************
 * Function:	yyerror							*
 *									*
 * Description:	Writes a syntax error message reporting the unexpected	*
 *		token only if input if being read from a file.		*
 ************************************************************************/

void Cable_error (msg)
    char *msg;
{
    extern int yychar;


    if (!ctl.problem -> input)
	if (yytext) {
	    error ("parse error before `%s'", yytext);
    }
	else {
	    error ("parse error before end of input");
    }
    else
	   ctl.problem -> input -= yyleng;
}


/************************************************************************
 * Function:	init_lexer						*
 *									*
 * Description:	Initializes the lexer for a new file.			*
 ************************************************************************/

void 
CableInitLexer(FILE *fp)
{
    static int first_time = 1;

    if (first_time) {
	    first_time = 0;
	    yyin = fp;
    } 
    else {
	    yyrestart (fp);
    }

    if (fp) {
	    ctl.problem -> input = NULL;
	    filename = ctl.problem -> filename;
    }
}


/************************************************************************
 * Function:	update_line						*
 *									*
 * Description:	Updates the current line and file from a preprocessor	*
 *		directive.						*
 ************************************************************************/

static void update_line ( )
{
     char	*ptr1;
     char	*ptr2;
     int         line;
     static char include_filename [256];


     ptr2 = include_filename;
     line = strtol (yytext + 1, &ptr1, 10) - 1;

     if (ptr1) {
	     ctl.problem -> line = line;

	    while (*ptr1 && *ptr1 ++ != '"');

	    if (*ptr1 && ctl.problem -> filename != filename)  {
	        ctl.problem -> filename = filename;
        }

	    while (*ptr1 && *ptr1 != '"')
	        *ptr2 ++ = *ptr1 ++;

	    *ptr2 = 0;
	    if (ptr2 != include_filename && ctl.problem -> filename != include_filename) {
	        ctl.problem -> filename = include_filename;
        }
    }
}


/************************************************************************
 * Function:	copy_input						*
 *									*
 * Description:	Specifies whether yytext should be copied into a local	*
 *		buffer.  The text is put into a canonical form in which	*
 *		leading and trailing white space is eliminated and all	*
 *		other sequences of white space are replaced by a space.	*
 ************************************************************************/

char *copy_input (flag)
    int flag;
{
    char       *src;
    char       *dest;
    static char canonical_text [65536];


    if (!(copy_flag = flag)) {
	    *(copy_ptr - yyleng) = 0;
	    dest = canonical_text;

	    src = strtok (copied_text, " \t\n");
	    while (src) {
	        while ((*dest ++ = *src ++));
	        src = strtok (NULL, " \t\n");
	        *(dest - 1) = ' ';
	    }

	    *(dest - 1) = 0;
	    *(copy_ptr = copied_text) = 0;
	    return canonical_text;
    }

    // copy_ptr = copied_text;
    return copied_text;
}

int Cable_wrap ( )
{
    if (0) unput (0);
# if 0
# ifdef YY_SC_TO_UI	/* eliminates warnings with gcc -Wall */
    if (0) {
	yy_push_state (yy_top_state ( ));
	yy_pop_state ( );
	yyunput (0, NULL);
    }
# endif
# endif
    return 1;
}
