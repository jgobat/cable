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
 * File:	model.h							*
 *									*
 * Description:	This file contains the primary type and function	*
 *		declarations for the mooring modeling library		*
 ************************************************************************/

# ifndef _MODEL_H
# define _MODEL_H
# include <stdio.h>


# define TINY 1.0e-60
# define UnspecifiedValue (-99999999)


/* Problem types */

typedef enum {
   General	= 1,
   Surface	= 2,
   Subsurface	= 3,
   Horizontal   = 4,
   Towing       = 5,
   Drifter      = 6,
   Deployment   = 7,
   Riser	= 8,
   Webster	= 9,
   HorizontalDrifter = 10,
} ProblemType;


/* Position integration schemes */

typedef enum {
   Spatial  = 1,
   Temporal = 2
} Integration;


/* Various kinds of solution algorithms */

typedef enum {
   Relaxation = 1,
   Catenary = 2,
   Shooting = 3
} SolutionMethod;

typedef enum {
   Bisection = 1,
   Secant = 2,
   Fixed = 3,
   Variable = 4,
} OuterMethod;

/* Forcing methods */

typedef enum {
   Velocity	= 1,		/* 1 - 3 require no buoy knowledge */
   WaveFollower = 2,
   Force	= 3,
   Spar  	= 4,
   FroudeKrylov = 5,
   Morison	= 6,
   LAMP		= 7,
} ForcingMethod;


/* Forcing input types */

typedef enum {
   Regular  = 1,
   Random  = 2 
} InputType;


/* Buoy types */

typedef enum {
    Axisymmetric = 0,		/* defined by diameters array	*/
    Sphere     = 1,		/* d in diameter		*/
    Cylinder   = 2, 		/* d around and h tall		*/
    Capsule    = 3,		/* d high and h long		*/
    Ship       = 4,		/* totally unspecified for now  */
    Platform   = 5, 		/* totally unspecified for now  */
} BuoyType;


/* Material types */

typedef enum {
   Linear	= 0,
   Nonlinear = 1
} MaterialType;


/* The types of equations that the solver can request from the
   Jacobian building routines -- also used for node position defintion */

typedef enum {
   Cable          = 1,
   BottomBoundary = 2,
   TopBoundary    = 3,
   BranchStart    = 4, 
   BranchTerminal = 5,
   Connection     = 6,
   Junction       = 7
} EquationType;

/* Variable expression */

typedef struct {
    double value;		/* value at time zero */
    Code   expr;		/* expression code    */
    char  *text;		/* text of expression */
    char  *name;
} VarExpr;

/* A pair structure to hold the buoy diameter information */

typedef struct pair {
    double   level;		/* height along the buoy for this diameter */
    double   d;                 /* diameter at this level                  */
} DiameterPair;

typedef struct cable_object {
    char *aux;
    char *name;
    char *category;
    char *color;
    char	*comment;
} *CableObject;

/* A material */

typedef struct material {
    char   	*aux;			/* auxillary data pointer        */
    char   	*name;			/* name of material              */
    char    *category;
    char   	*color;			/* name of color to use          */
    char	*comment;
    double   d;			/* diameter			 */
    VarExpr  Cdt;			/* tangential drag coefficient	 */
    VarExpr  Cdn;			/* normal drag coefficient	 */
    double   Cmn;			/* normal added mass coefficient */
    double	 Cmt;			/* tangential added mass coeff   */
    double   Can;			/* normal added mass coefficient */
    double	 Cat;			/* tangential added mass coeff   */
    double   bt;			/* longitudinal damping constant */
    double   bn;			/* transverse damping constant   */
    double   w;			/* dry weight per length         */
    double   wet;			/* wet weight per length         */
    double	 rV;			/* displaced mass per length	 */
    double   A;			/* cross-sectional area	         */
    double   EA;			/* derived stiffness		 */
    double   EI;			/* derived flexural stiffness	 */
    double   GJ;			/* derived torsional stiffness   */
    double	 E;			/* Young's modulus		*/
    double	 I;			/* area moment of inertia	*/
    double	 J;			/* polar moment of inertia	*/
    double	 G;			/* shear modulus		*/
    double	 nu;			/* Poisson's ratio		*/
    double	 id;			/* inside diameter		*/
    double   m;			/* mass per length	        */
    double	 length;		/* optional fixed length	*/
    double   amt;			/* normal added mass/length     */
    double	 amn;			/* tangential added mass/length */
    double	 rho;			/* material volumetric density  */
    double   mud_b;			/* damping on bottom		*/
    VarExpr	 T [3];			/* T(e), T'(e), T''(e) functions */
    double   swl;           // maximum safe working load
    double   yield;         // yield or failure strength
    MaterialType type;			/* type of constitutive relation */
} *Material;


/* A buoy */

typedef struct buoy {
    char	 *aux;
    char	 *name;
    char     *category;
    char	 *color;
    char	 *comment;
    BuoyType	  type;
    double	  S;			/* projected drag area		*/
    double	  cg;			/* location of cg above base	*/
    double	  d;			/* diameter			*/
    double	  h;			/* overall height		*/
    double	  m;			/* mass				*/
    double	  w;			/* weight in air                */
    double	  draft;		/* buoy static draft		*/
    double	  max_draft;		/* maximum draft of buoy	*/
    double	  min_draft;		/* minimum draft of buoy	*/
    double	  buoyancy;		/* total available buoyancy	*/
    double	  Cl;			/* maximum net lift coeff	*/
    double	  Cdt;			/* tangential drag coefficient	*/
    double	  Cdn;			/* normal drag coefficient	*/
    double        Cdw;			/* wind drag coefficient	*/
    double	  Sw;			/* surface area exposed to wind */
    double	  am;			/* Radiation added mass		*/
    double	  b;			/* Radiation damping		*/
    double	  Mdr;			/* Morison total drag coeff     */
    double	  Mmma;			/* Morison mass + added mass	*/
    double	  Marma;		/* Morison added + ghost mass   */ 
    DiameterPair *diameters;		/* body of revolution info	*/
    int		  num_diameters;
} *Buoy;


/* An anchor */

typedef struct anchor {
    char	 *aux;
    char	 *name;
    char     *category;
    char	 *color;
    char	 *comment;
    double	  wet;
    double	  m;
    double	  am;
    double	  Cdn; 
    double	  Cdt;
    double	  S;
    double	  d;
    double	  mu;
} *Anchor;


/* A terminal at the beginning or end of a line */

typedef struct terminal {
    Anchor	  anchor;
    Buoy	  buoy;
    struct node	 *node;		  // node at terminal
    struct node  *loop_main_node;	// node on mainline for looped branch
    double	  x;
    double	  y;
    double	  z;
    double	  xforce;
    double	  yforce;
    double	  zforce;
    double	  initial_xforce;
    double	  initial_yforce;
    double	  initial_zforce;
    double	  tension;
    double	  phi;
    VarExpr	  xspeed;
    VarExpr	  yspeed;
    VarExpr	  zspeed;
    VarExpr	  xthrust;
    VarExpr	  ythrust;
    VarExpr	  zthrust;
    VarExpr	  profile;
    double	  profile_m;
    double	  profile_H;
    double	  profile_turn;
    double	  Ki;
    double	  Kd;
    double	  Kp;
    double	  release;
    double	  buoyancy_head;
    char	 *flap_file;
    double    friction;
    double    safety;
} *Terminal;


/* A connector (massive rigid object between segments) */

typedef struct connector {
    char	 *aux;
    char	 *name;
    char     *category;
    char	 *color;
    char	 *comment;
    double	  m;			/* mass				 */
    double	  am;			/* added mass			 */
    double	  Cdt;			/* tangential drag coefficient	 */
    double	  Cdn;			/* normal drag coefficient	 */
    double	  wet;			/* displacement			 */
    double	  d;			/* characteristic drag diameter  */
    double	  length;		/* length to give dimensionality */
} *Connector;

typedef struct {
   int		num_nodes;
   struct node *node[10];
} junction;
  
typedef struct {
   Connector	object;
   int		num_nodes;
   int	       *nodes;
} attachment; 

typedef struct {
   int		nodes; 
   double	percent;
} distribution;

	/*
	 * we only need this for the parser
	 */

typedef struct {
   double amplitude;
   double period;
   double phase;
} wave;

typedef enum {
   Spliced = 0,
   Pinned = 1, 
} ConnectionType;

/* A cable (a continuous length of of a single material type) */

typedef struct segment {
    char 	   *aux;
    char       *name;
    unsigned	    number;	/* segment number		   */
    struct segment *next;
    struct segment *prev;
    struct node	   *first;
    struct node    *last;
    struct node    *first_active;
    struct node    *last_active;
    unsigned	    num_nodes;	/* total number of nodes	   */
    unsigned	    num_dist;	/* number of distribution pairs	   */
    distribution    dist [10];	/* percent length of nodes	   */
    unsigned	    num_top_dist;
    unsigned	    num_bottom_dist;
    distribution    top_dist [10];
    distribution    bottom_dist [10];
    Material	    material;	/* material type for this cable    */
    double	        length;	    /* natural length of the segment   */
    Connector	    connector;	/* the connector above the segment */
    ConnectionType  connection; /* connection type when connector not spec'd */
    double          connector_x;
    double          connector_y;
    double          connector_z;
    double          connector_xforce;
    double          connector_yforce;
    double          connector_zforce;
    VarExpr         connector_xthrust;
    VarExpr         connector_ythrust;
    VarExpr         connector_zthrust;
    attachment     *attach;	/* array of attachment structures  */
    unsigned	    num_attach;	/* number of attachment structures */
    struct branch  *branch;	/* parent branch 		   */
    struct branch  *branch_to [10]; /* branches leaving from this segment   */
    unsigned        num_branch_to;   /* number of branches from this segment */
    junction	    junction;
    VarExpr	    top_pay;
    VarExpr	    bottom_pay;
    double	    top_wet;
    double	    bottom_wet;
    double	    top_length;
    double	    bottom_length;
    double	    top_spooled;
    double	    bottom_spooled;
    double	    top_added;
    double	    bottom_added;
    double      top_pay_speed;
    double      bottom_pay_speed;
    double      top_pay_speed_o;
    double      bottom_pay_speed_o;
    double      prev_t;
    double      prev_dt;
} *Segment;

/* A branch - a series of segments w/terminal at end that splits off
   from mainline trunk of segments */

typedef struct branch {
    char  	*aux;
    Segment      segment_from;
    Segment      segment [10];
    int		 num_segment;
    struct node *first;
    struct node *last;
    int		 number;	/* branch number			    */
    Terminal	 terminal;	/* terminal at end of branch		    */
    double	 w0;		/* equiv weight for catenary calculations   */
    double	 EA;		/* equiv axial stiffness for catenary calc  */
    double	 length;	/* total length for catenary calculations   */
} *Branch;

/* A node on a segment (a single material point) */

typedef struct node {
    char	*aux;
    unsigned	 number;
    unsigned     active_number;		
    unsigned	 output_number;		/* output order number		*/
    struct node *next;
    struct node *prev;
    struct node *next_active;
    struct node *prev_active;
    Segment	     segment;		/* the parent segment	         */
    Material	 material;		/* alias to segment's material   */
    EquationType position;		
    Connector	 attachment;	/* optional attached mass        */
    double	 s;
    double	 ds;
    double	 ds0;
    int		 active;
    double	*Ys;
    double	*Y;
    double  *Yd;
    double	*Y_o;
    double	*Yd_o;
    double	*Y_f;
    double	*Y_o_f;

    // x,y,z,xdot, etc. vars need to maintain this order for
    // the progress saver to work as written
    double	 x;
    double	 y;
    double	 z;
    double   x_o;
    double	 y_o;
    double	 z_o;

    double	 x_f;
    double	 y_f;
    double	 z_f;
    double   x_o_f;
    double	 y_o_f;
    double	 z_o_f;

    double	 xdot;
    double	 ydot;
    double	 zdot;
    double	 xdot_o;
    double	 ydot_o;
    double	 zdot_o;

    double	 xdot_f;
    double	 ydot_f;
    double	 zdot_f;
    double	 xdot_o_f;
    double	 ydot_o_f;
    double	 zdot_o_f;

    double	 pay;
    double	 pay_o;
    double	 pay_f;
    double	 pay_o_f;

    double   drag;
    double   lift;
} *Node;
 
	/*
	 * routines in expressions.c
	 */

extern void Current (	
   double,		/* current time		  	*/
   double,      /* depth */
   double,
   double,
   double *,		/* X return current		*/
   double *,		/* Y return current		*/
   double *		/* Z return current		*/

);

extern void DragCoeff (
   Node, 
   Material,
   double, // t
   double, // u_rel
   double, // v_rel
   double *, // Cdt
   double *  // Cdn
);

extern void Speed (
   double,		/* time				*/
   Terminal,		/* terminal point structure	*/
   double *,		/* X speed			*/
   double *, 		/* Y speed			*/
   double *		/* Z speed			*/
);

extern void ConnectorThrust (
   double,		/* time				*/
   Segment,
   Node,
   double *,		/* X thrust			*/
   double *, 		/* Y thrust			*/
   double *		/* Z thrust			*/
);

extern void Thrust (
   double,		/* time				*/
   Terminal,		/* terminal point structure	*/
   double *,		/* X thrust			*/
   double *, 		/* Y thrust			*/
   double *		/* Z thrust			*/
);

extern double Profile (
   double 
);


extern double Bottom (
   double,      // horiz  coord
   double,      // 3D horiz coord
   double       // time t
);

extern void WindDrag (
   double,		/* time			*/
   Buoy,		/* buoy			*/
   double *,		/* 2D force return 	*/
   double *		/* 3D force return	*/
);

extern double ProjectedArea (
   Buoy,
   double
);

extern double WaterplaneArea (
   Buoy,
   double
);


extern void WaveParticleMotion (
   double,		/* time			*/
   double,		/* x			*/
   double,		/* y			*/
   double,		/* z			*/
   double *,		/* u return		*/
   double *,		/* v return		*/
   double *,		/* w return		*/
   double *,		/* ud return		*/
   double *,		/* vd return		*/
   double *		/* wd return		*/
);

extern void WaveParticleVelocity (
   double,		/* time			*/
   double,		/* x			*/
   double,		/* y			*/
   double,		/* z			*/
   double *,		/* u return		*/
   double *,		/* v return		*/
   double *		/* w return		*/
);

extern void WaveParticleAcceleration (
   double,		/* time			*/
   double,		/* x			*/
   double,		/* y			*/
   double,		/* z			*/
   double *,		/* u return		*/
   double *,		/* v return		*/
   double *		/* w return		*/
);

extern void WaveSurfaceVelocity (
   double,		/* time			*/
   double,		/* x			*/
   double,		/* y			*/
   double,		/* z			*/
   double *,		/* u return		*/
   double *,		/* v return		*/
   double *		/* w return		*/
);

extern void InputVelocity (
   double,		/* time			*/
   double,		/* x			*/
   double,		/* y			*/
   double,		/* z			*/
   double *,		/* u return		*/
   double *,		/* v return		*/
   double *		/* w return		*/
);

extern void FroudeKrylovCoefficients (
   double,		/* time			*/
   Buoy,		/* buoy			*/
   double,		/* vertical position	*/
   double,		/* horizontal position  */
   double,		/* 3D horiz position    */
   double *,		/* Fex return		*/
   double *,		/* wave damping return  */
   double *		/* net buoyancy return  */
);

extern void FlapForces (
   double,		// t
   double,		// dt
   double,		// curr_x
   double,		// prev_x
   double,		// curr_y
   double		// prev_y
);

extern int TowDepthControl (
   double,		// t
   double,		// dt
   double,		// curr_x
   double,		// prev_x
   double,		// curr_y
   double,		// prev_y
   double		// length
);

extern int CheckTypeParameters ( );
extern int CheckStaticParameters ( );
extern int CheckDynamicParameters ( );
extern int CheckEnvironmentParameters ( );
extern int CheckMaterialProperties (int);
extern int CheckBuoyProperties ( );
extern int CheckBranchTerminalProperties ( );

# endif /* _MODEL_H */
