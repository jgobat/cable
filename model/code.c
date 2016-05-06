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
 * File:	code.c							*
 *									*
 * Description:	This file contains the public and private function and	*
 *		type definitions for the stack machine code.		*
 ************************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "code.h"
# include "allocate.h"
# include "model.h"
# include "problem.h"
# include "tension.h"
# include "rotate.h"
# include <stdarg.h>
# define VA_START(a,b)      va_start(a,b)

# define MaxStackDepth	4096
# define MaxCodeSize	4096

# define push(x)	(* ++ sp = (x))
# define pop()		(* sp --)
# define top()		(* sp)

# define deg2rad(x)	((x) * M_PI / 180.0)


# define None	 0
# define Integer 1
# define Double  2
# define Array	 3

static struct {
    char *opcode;
    int   arg_type;
} data [ ] = {
    {"jmp",   Integer},
    {"jnz",   Integer},
    {"jz",    Integer},
    {"push",  Double},
    {"pop",   None},
    {"copy",  None},
    {"test",  None},
    {"neg",   None},
    {"not",   None},
    {"inv",   None},
    {"mul",   None},
    {"div",   None},
    {"mod",   None},
    {"add",   None},
    {"sub",   None},
    {"lsft",  None},
    {"rsht",  None},
    {"lt",    None},
    {"gt",    None},
    {"lteq",  None},
    {"gteq",  None},
    {"eq",    None},
    {"neq",   None},
    {"and",   None},
    {"xor",   None},
    {"or",    None},
    {"vart",  None},
    {"varn",  None},
    {"vars",  None},
    {"vare",  None},
    {"varSn",  None},
    {"varSb",  None},
    {"varu",  None},
    {"varv",  None},
    {"varw",  None},
    {"varB0",  None},
    {"varB1",  None},
    {"varB2",  None},
    {"varB3",  None},
    {"varOm1",  None},
    {"varOm2",  None},
    {"varOm3",  None},
    {"varphi", None},
    {"varF",  None},
    {"varx",  None},
    {"vary",  None},
    {"varz",  None},
    {"varUx",  None},
    {"varVy",  None},
    {"varWz",  None},
    {"varFx",  None},
    {"varFy",  None},
    {"varFz",  None},
    {"varh",   None},
    {"varp",  None},
    {"first", None},
    {"last", None},
    {"firstactive", None},
    {"lastactive", None},
    {"sin",   None},
    {"cos",   None},
    {"tan",   None},
    {"sinh",  None},
    {"cosh",  None},
    {"tanh",  None},
    {"asin",  None},
    {"acos",  None},
    {"atan2", None},
    {"exp",   None},
    {"ln",    None},
    {"log",   None},
    {"pow",   None},
    {"sqrt",  None},
    {"hypot", None},
    {"floor", None},
    {"ceil",  None},
    {"fmod",  None},
    {"fabs",  None},
    {"table", Array},
    {"cycle", Array},
    {"halt",  None},
};


typedef union instruction {
    int    op;
    int    offset;
    double arg;
} Instruction;


static Instruction in_core [MaxCodeSize];
static Code	   ip = in_core;

static double	   stack [MaxStackDepth];
static double	  *sp;

Code InCore = in_core;


/************************************************************************
 * Function:	EmitCode						*
 *									*
 * Description:	Adds an instruction to the current piece of code.	*
 ************************************************************************/

void 
EmitCode (int op, ...)
{
    va_list ap;
    int     i;
    int     length;
    double *array;


    ip ++ -> op = op;

    op = op & 0x3ff;

    switch (data [op].arg_type) {
    case Integer:
	VA_START (ap, op);
	ip ++ -> offset = va_arg (ap, int);
	va_end (ap);
	break;

    case Double:
	VA_START (ap, op);
	ip ++ -> arg = va_arg (ap, double);
	va_end (ap);
	break;

    case Array:
	VA_START (ap, op);
	array = va_arg (ap, double *);
	ip ++ -> offset = length = va_arg (ap, int);
	va_end (ap);
	for (i = 0; i < length; i ++)
	    ip ++ -> arg = array [i];
	break;
    }
}


/************************************************************************
 * Function:	CopyCode						*
 *									*
 * Description:	Copies a piece of code.					*
 ************************************************************************/

Code CopyCode (code)
    Code code;
{
    Code     pc;
    Code     ptr;
    Code     copy;
    Opcode   op;
    unsigned size;


    size = 0;
    if (!(pc = code))
	return NULL;

    while (1) {
	size ++;
	if ((op = pc ++ -> op) == HaltOp)
	    break;
	if (data [op].arg_type == Array) {
	    size += pc -> offset;
	    pc += pc -> offset;
	}
	if (data [op].arg_type != None) {
	    size ++;
	    pc ++;
	}
    }

    if (!(copy = Allocate (Instruction, size)))
	return NULL;

    pc = code;
    ptr = copy;

    while (size --)
	*ptr ++ = *pc ++;

    return copy;
}


/************************************************************************
 * Function:	FreeCode						*
 *									*
 * Description:	Deallocates a copied program.				*
 ************************************************************************/

void FreeCode (code)
    Code code;
{
    if (code != InCore)
	Deallocate (code);
}


/************************************************************************
 * Function:	SetIP							*
 *									*
 * Description:	Sets the address of the instruction pointer.		*
 ************************************************************************/

void SetIP (new_ip)
    int new_ip;
{
    ip = InCore + new_ip;
}


/************************************************************************
 * Function:	GetIP							*
 *									*
 * Description:	Returns the address of the instruction pointer.		*
 ************************************************************************/

int GetIP ( )
{
    return ip - InCore;
}


/************************************************************************
 * Function:	EvalTable						*
 *									*
 * Description:	Evaluates a table or cycle opcode.			*
 ************************************************************************/

static double EvalTable (array, length, depth, flag)
    double *array;
    int     length;
    double  depth;
    int     flag;
{
    int    i1;
    int    i2;
    double max;
    double t1;
    double t2;
    double x1;
    double x2;


    if (length == 2)
	return depth == array [0] ? array [1] : 0;

    if (flag) {
	max = array [length - 2];
	while (depth > max)
	    depth -= max;
    }

    for (i1 = 0; i1 < length; i1 += 2)
	if (depth <= array [i1])
	    break;

    if (i1 == length) {
	i1 = i1 - 2;
	i2 = i1 - 2;
    } else if (i1)
	i2 = i1 - 2;
    else
	i2 = i1 + 2;

    t1 = array [i1];
    t2 = array [i2];
    x1 = array [i1 + 1];
    x2 = array [i2 + 1];

    /* fprintf (stderr,"(%g,%g) %g (%g,%g)\n", t1, x1, depth, t2, x2); */
    return t1 != t2 ? (x2 - x1) / (t2 - t1) * (depth - t1) + x1 : x2;
}

static void 
rotate(int twoD, int dynamic, int which, Opcode op, Node nd, double *Fx, double *Fy, double *Fz)
{
    int  off;

    if (twoD && dynamic)
        off = 5;
    else if (twoD)
        off = 4;
    else if (!twoD && dynamic)
        off = 7;
    else 
        off = 4;

    if (which == OLDNODEDATA && !(op & FILTERNODEDATA))
        RotateToFixed(Tension(nd -> Y_o[1], nd -> material),
                      nd -> Y_o[2],
                      twoD ? 0 : nd -> Y_o[3],
                      nd -> Y_o[off],
                      twoD ? 0 : nd -> Y_o[off + 1],
                      twoD ? 0 : nd -> Y_o[off + 2],
                      twoD ? 0 : nd -> Y_o[off + 3],
                      Fx, Fy, Fz, twoD);
    else if (which == OLDNODEDATA && (op & FILTERNODEDATA))
        RotateToFixed(Tension(nd -> Y_o_f[1], nd -> material),
                      nd -> Y_o_f[2],
                      twoD ? 0 : nd -> Y_o_f[3],
                      nd -> Y_o_f[off],
                      twoD ? 0 : nd -> Y_o_f[off + 1],
                      twoD ? 0 : nd -> Y_o_f[off + 2],
                      twoD ? 0 : nd -> Y_o_f[off + 3],
                      Fx, Fy, Fz, twoD);
    else if (op & FILTERNODEDATA)
        RotateToFixed(Tension(nd -> Y_f[1], nd -> material),
                      nd -> Y_f[2],
                      twoD ? 0 : nd -> Y_f[3],
                      nd -> Y_f[off],
                      twoD ? 0 : nd -> Y_f[off + 1],
                      twoD ? 0 : nd -> Y_f[off + 2],
                      twoD ? 0 : nd -> Y_f[off + 3],
                      Fx, Fy, Fz, twoD);
    else
        RotateToFixed(Tension(nd -> Y[1], nd -> material),
                      nd -> Y[2],
                      twoD ? 0 : nd -> Y[3],
                      nd -> Y[off],
                      twoD ? 0 : nd -> Y[off + 1],
                      twoD ? 0 : nd -> Y[off + 2],
                      twoD ? 0 : nd -> Y[off + 3],
                      Fx, Fy, Fz, twoD);
}

extern Problem *problem;
extern Environment *environment;

double
selectY(int x, Node n, Node *nodes, int nn, int which, Opcode op, int i)
{
    // printf("op = %d, x = %d, w = %d\n", op, x, which);
    if (op & FILTERNODEDATA) {
        if (x == 0) 
            return (which == OLDNODEDATA ? n -> Y_o_f[i] : n -> Y_f[i]);
        else if (x >= 1 && x <= nn) 
            return (which == OLDNODEDATA ? nodes[x] -> Y_o_f[i] : nodes[x] -> Y_f[i]);
        else 
            return 0.0;
    }
    else {
        if (x == 0) 
            return (which == OLDNODEDATA ? n -> Y_o[i] : n -> Y[i]);
        else if (x >= 1 && x <= nn) 
            return (which == OLDNODEDATA ? nodes[x] -> Y_o[i] : nodes[x] -> Y[i]);
        else 
            return 0.0;
    }
}

double
EvalCode(Code code, Node n, 
         double t, double var, 
         double Rx, double Ry, double Rz, int which)
{
    int     x;
    int     y;
    double  a;
    double  b;
    Code    pc;
    Node   *nodes;
    Segment *seg;
    int      nn, ns;
    int      twoD, dynamic;
    int      rot;
    double   Fx, Fy, Fz;
    double   e;
    Node     nd;
    Opcode   op;

    rot = 0;

    if (problem) {
        nn    = problem -> num_nodes;
        ns    = problem -> num_segments;
        nodes = problem -> node;
        seg   = problem -> segment;
        twoD  = problem -> twoD;
        dynamic = problem -> dynamic; 
    }
    else {
        nn = ns = 0;
        nodes = NULL;
        seg = NULL;
        twoD = 1;
        dynamic = 0;
    }

    sp = stack;
    if (!(pc = code))
	return 0;

    while (1) {
    
    op = pc ++ -> op;
    // printf("op = %d, masked = %d\n", op, op & 0x3ff);

	switch (op & 0x3ff) {
	case JmpOp:
	    y = pc ++ -> offset;
	    pc += y;
	    break;

	case JnzOp:
	    x = pop ( );
	    y = pc ++ -> offset;
	    if (x) pc += y;
	    break;

	case JzOp:
	    x = pop ( );
	    y = pc ++ -> offset;
	    if (!x) pc += y;
	    break;

	case PushOp:
        // printf("pushing %g\n", pc -> arg);
	    push (pc ++ -> arg);
	    break;

	case PopOp:
	    a = pop ( );
	    break;

	case CopyOp:
	    a = top ( );
	    push (a);
	    break;

	case TestOp:
	    a = pop ( );
	    push (a != 0);
	    break;

	case NegOp:
	    a = pop ( );
	    push (-a);
	    break;

	case NotOp:
	    a = pop ( );
	    push (!a);
	    break;

	case InvOp:
	    x = pop ( );
	    push (~x);
	    break;

	case MulOp:
	    b = pop ( );
	    a = pop ( );
        // printf("a=%g, b=%g\n", a, b);
	    push (a * b);
	    break;

	case DivOp:
	    b = pop ( );
	    a = pop ( );
	    push (b ? a / b : 0);
	    break;

	case ModOp:
	    y = pop ( );
	    x = pop ( );
	    push (y ? x % y : 0);
	    break;

	case AddOp:
	    b = pop ( );
	    a = pop ( );
	    push (a + b);
	    break;

	case SubOp:
	    b = pop ( );
	    a = pop ( );
	    push (a - b);
	    break;

	case LsftOp:
	    y = pop ( );
	    x = pop ( );
	    push (x << y);
	    break;

	case RsftOp:
	    y = pop ( );
	    x = pop ( );
	    push (x >> y);
	    break;

	case LtOp:
	    b = pop ( );
	    a = pop ( );
	    push (a < b);
	    break;

	case GtOp:
	    b = pop ( );
	    a = pop ( );
	    push (a > b);
	    break;

	case LteqOp:
	    b = pop ( );
	    a = pop ( );
	    push (a <= b);
	    break;

	case GteqOp:
	    b = pop ( );
	    a = pop ( );
	    push (a >= b);
	    break;

	case EqOp:
	    b = pop ( );
	    a = pop ( );
	    push (a == b);
	    break;

	case NeqOp:
	    b = pop ( );
	    a = pop ( );
	    push (a != b);
	    break;

	case AndOp:
	    y = pop ( );
	    x = pop ( );
	    push (x & y);
	    break;

	case XorOp:
	    y = pop ( );
	    x = pop ( );
	    push (x ^ y);
	    break;

	case OrOp:
	    y = pop ( );
	    x = pop ( );
	    push (x | y);
	    break;

    case VartOp:
        push(t);
        break;

    case VarNOp:
        if (n == NULL)
            push(0);
        else 
            push(n -> number);
        break;

    case VarSOp:
        x = pop( );
        if (n == NULL || nodes == NULL || x < 0 || x > nn) {
            push(var);
        }
        else if (x >= 1 && x <= nn) {
            push(nodes[x] -> s);
        }
        else  {
            push(n -> s);
        }
        break;

    case VarEOp:
        x = pop( );
        if (n == NULL || nodes == NULL)
            push(var);
        else {
            push(selectY(x, n, nodes, nn, which, op, 1));
        }
	    break;

    case VarSnOp:
        x = pop( );
        if (n == NULL || nodes == NULL) 
            push(0.0);
        else {
            push(selectY(x, n, nodes, nn, which, op, 2));
        }
	    break;

    case VarSbOp:
        x = pop( );
        if (twoD || n == NULL || nodes == NULL) 
            push(0.0);
        else  {
            push(selectY(x, n, nodes, nn, which, op, 3));
        }
	    break;

    case VarUOp:
        x = pop( ); 
        if (n == NULL || nodes == NULL) 
            push(Rx);
        else if (!dynamic)
            push(0.0);
        else 
            push(selectY(x, n, nodes, nn, which, op, twoD ? 3 : 4));
        break;

    case VarVOp:
        x = pop( );
        if (n == NULL || nodes == NULL) 
            push(Ry);
        else if (!dynamic)
            push(0.0);
        else 
            push(selectY(x, n, nodes, nn, which, op, twoD ? 4 : 5));
        break;

    case VarWOp:
        x = pop( );
        if (n == NULL || nodes == NULL) 
            push(Rz);
        else if (!dynamic || twoD)
            push(0.0);
        else 
            push(selectY(x, n, nodes, nn, which, op, 6));
        break;

    case VarB0Op:
        x = pop( );
        if (twoD || n == NULL || nodes == NULL)
            push(0.0);
        else {
            push(selectY(x, n, nodes, nn, which, op, dynamic ? 7 : 4));
        }
        break;

    case VarB1Op:
        x = pop( );
        if (twoD || n == NULL || nodes == NULL)
            push(0.0);
        else {
            push(selectY(x, n, nodes, nn, which, op, dynamic ? 8 : 5));
        }
        break;

    case VarB2Op:
        x = pop( );
        if (twoD || n == NULL || nodes == NULL)
            push(0.0);
        else {
            push(selectY(x, n, nodes, nn, which, op, dynamic ? 9 : 6));
        }
        break;

    case VarB3Op:
        x = pop( );
        if (twoD || n == NULL || nodes == NULL)
            push(0.0);
        else {
            push(selectY(x, n, nodes, nn, which, op, dynamic ? 10 : 7));
        }
        break;

    case VarOm1Op:
        x = pop( );
        if (twoD || n == NULL || nodes == NULL)
            push(0.0);
        else {
            push(selectY(x, n, nodes, nn, which, op, dynamic ? 11 : 8));
        }
        break;

    case VarOm2Op:
        x = pop( );
        if (twoD || n == NULL || nodes == NULL)
            push(0.0);
        else {
            push(selectY(x, n, nodes, nn, which, op, dynamic ? 12 : 9));
        }
        break;

    case VarOm3Op:
        x = pop( );
        if (n == NULL || nodes == NULL)
            push(0.0);
        else {
            if (twoD) 
                push(selectY(x, n, nodes, nn, which, op, dynamic ? 6 : 4));
            else 
                push(selectY(x, n, nodes, nn, which, op, dynamic ? 13 : 10));
        }
        break;

    case VarPhiOp:
        x = pop( );
        if (!twoD || n == NULL || nodes == NULL)
            push(0.0);
        else {
            push(selectY(x, n, nodes, nn, which, op, dynamic ? 5 : 3));
        }
        break;

	case VarTOp:
        x = pop( );
        if (n == NULL || nodes == NULL)
            push(var);
        else {
            if (x >= 0 && x <= nn) {
                e = selectY(x, n, nodes, nn, which, op, 1);
                push(Tension(e, n -> material));
            }
            else  {
                push(0.0);
            }
        }
	    break;

	case VarXOp:
        x = pop( );
        if (n == NULL || nodes == NULL) {
            push(Ry);
        }
        else {
            if (op & FILTERNODEDATA) {
                if (x == 0) push(which == OLDNODEDATA ? n -> y_o_f : n -> y_f);
                else if (x >= 1 && x <= nn) push(which == OLDNODEDATA ? nodes[x] -> y_o_f : nodes[x] -> y_f);
                else push(0.0);
            }
            else {
                if (x == 0) push(which == OLDNODEDATA ? n -> y_o : n -> y);
                else if (x >= 1 && x <= nn) push(which == OLDNODEDATA ? nodes[x] -> y_o : nodes[x] -> y);
                else push(0.0);
            }
        }
	    break;

	case VarYOp:
        x = pop( );
        if (twoD) 
            push(0.0);
        else if (n == NULL || nodes == NULL)  {
            push(Rz);
        }
        else {
            if (op & FILTERNODEDATA) {
                if (x == 0) push(which == OLDNODEDATA ? n -> z_o_f : n -> z_f);
                else if (x >= 1 && x <= nn) push(which == OLDNODEDATA ? nodes[x] -> z_o_f : nodes[x] -> z_f);
                else push(0.0);
            }
            else {
                if (x == 0) push(which == OLDNODEDATA ? n -> z_o : n -> z);
                else if (x >= 1 && x <= nn) push(which == OLDNODEDATA ? nodes[x] -> z_o : nodes[x] -> z);
                else push(0.0);
            }
        }
	    break;

	case VarZOp:
        x = pop( );
        if (n == NULL || nodes == NULL) {
            push(Rx);
        }
        else {
            if (op & FILTERNODEDATA) {
                if (x == 0) push(which == OLDNODEDATA ? n -> x_o_f : n -> x_f);
                else if (x >= 1 && x <= nn) push(which == OLDNODEDATA ? nodes[x] -> x_o_f : nodes[x] -> x_f);
                else push(0.0);
            }
            else {
                if (x == 0) push(which == OLDNODEDATA ? n -> x_o : n -> x);
                else if (x >= 1 && x <= nn) push(which == OLDNODEDATA ? nodes[x] -> x_o : nodes[x] -> x);
                else push(0.0);
            }
        }
	    break;

	case VarUxOp:
        x = pop( );
        if (n == NULL || nodes == NULL)  {
            push(0.0);
        }
        else {
            if (op & FILTERNODEDATA) {
                if (x == 0) push(which == OLDNODEDATA ? n -> ydot_o_f : n -> ydot_f);
                else if (x >= 1 && x <= nn) push(which == OLDNODEDATA ? nodes[x] -> ydot_o_f : nodes[x] -> ydot_f);
                else push(0.0);
            }
            else {
                if (x == 0) push(which == OLDNODEDATA ? n -> ydot_o : n -> ydot);
                else if (x >= 1 && x <= nn) push(which == OLDNODEDATA ? nodes[x] -> ydot_o : nodes[x] -> ydot);
                else push(0.0);
            }
        }
	    break;

	case VarVyOp:
        x = pop( );
        if (twoD || n == NULL || nodes == NULL) 
            push(0.0);
        else {
            if (op & FILTERNODEDATA) {
                if (x == 0) push(which == OLDNODEDATA ? n -> zdot_o_f : n -> zdot_f);
                else if (x >= 1 && x <= nn) push(which == OLDNODEDATA ? nodes[x] -> zdot_o_f : nodes[x] -> zdot_f);
                else push(0.0);
            }
            else {
                if (x == 0) push(which == OLDNODEDATA ? n -> zdot_o : n -> zdot);
                else if (x >= 1 && x <= nn) push(which == OLDNODEDATA ? nodes[x] -> zdot_o : nodes[x] -> zdot);
                else push(0.0);
            }
        }
	    break;

	case VarWzOp:
        x = pop( );
        if (n == NULL || nodes == NULL)
            push(0.0);
        else {
            if (op & FILTERNODEDATA) {
                if (x == 0) push(which == OLDNODEDATA ? n -> xdot_o_f : n -> xdot_f);
                else if (x >= 1 && x <= nn) push(which == OLDNODEDATA ? nodes[x] -> xdot_o_f : nodes[x] -> xdot_f);
                else push(0.0);
            }
            else {
                if (x == 0) push(which == OLDNODEDATA ? n -> xdot_o : n -> xdot);
                else if (x >= 1 && x <= nn) push(which == OLDNODEDATA ? nodes[x] -> xdot_o : nodes[x] -> xdot);
                else push(0.0);
            }
        }
	    break;

	case VarFxOp:
        x = pop( );
        if (n == NULL || nodes == NULL) 
            push(0.0);
        else {
            if (x == 0) nd = n;
            else if (x >= 1 && x <= nn) nd = nodes[x];
            else nd = NULL;

            if (nd && !rot) {
                rotate(twoD, dynamic, which, op, nd, &Fx, &Fy, &Fz);
                rot = 1;
            }
            if (nd) push(Fy);
            else    push(0.0);
        }
	    break;

	case VarFyOp:
        x = pop( );
        if (twoD || n == NULL || nodes == NULL) 
            push(0.0);
        else {
            if (x == 0) nd = n;
            else if (x >= 1 && x <= nn) nd = nodes[x];
            else nd = NULL;

            if (nd && !rot) {
                rotate(twoD, dynamic, which, op, nd, &Fx, &Fy, &Fz);
                rot = 1;
            }
            if (nd) push(Fz);
            else    push(0.0);
        }
	    break;

	case VarFzOp:
        x = pop( );
        if (n == NULL || nodes == NULL) 
            push(0.0);
        else {
            if (x == 0) nd = n;
            else if (x >= 1 && x <= nn) nd = nodes[x];
            else nd = NULL;

            if (nd && !rot) {
                rotate(twoD, dynamic, which, op, nd, &Fx, &Fy, &Fz);
                rot = 1;
            }
            if (nd) push(Fx);
            else    push(0.0);
        }
	    break;

    case VarHOp:
        x = pop( );
        if (environment && n == NULL) {
            push(environment -> surface - Rx);
        }
        else if (environment && n && nodes) {
            if (which == OLDNODEDATA && x >= 0 && x <= nn) {
                if (op & FILTERNODEDATA) { 
                    if (x == 0) push(environment -> surface - n -> x_o_f);
                    else        push(environment -> surface - nodes[x] -> x_o_f);
                }
                else {
                    if (x == 0) push(environment -> surface - n -> x_o);
                    else        push(environment -> surface - nodes[x] -> x_o);
                }
            }
            else if (x >= 0 && x <= nn) {
                if (op & FILTERNODEDATA) {
                    if (x == 0) push(environment -> surface - n -> x_f);
                    else        push(environment -> surface - nodes[x] -> x_f);
                }
                else {
                    if (x == 0) push(environment -> surface - n -> x);
                    else        push(environment -> surface - nodes[x] -> x);
                }
            }
            else
                push(0.0);
        }
        else 
            push(0.0);

        break; 
       
    case VarPOp:
        x = pop( );
        if (n == NULL || nodes == NULL)
            push(0.0);
        else {
            if (op & FILTERNODEDATA) {
                if (x == 0) push(which == OLDNODEDATA ? n -> pay_o_f : n -> pay_f);
                else if (x >= 1 && x <= nn) push(which == OLDNODEDATA ? nodes[x] -> pay_o_f : nodes[x] -> pay_f);
                else push(0.0);
            }
            else {
                if (x == 0) push(which == OLDNODEDATA ? n -> pay_o : n -> pay);
                else if (x >= 1 && x <= nn) push(which == OLDNODEDATA ? nodes[x] -> pay_o : nodes[x] -> pay);
                else push(0.0);
            }
        }
        break; 
 
    case FirstOp:
        x = pop( );
        if (x >= 1 && x <= ns) push((double) seg[x] -> first -> number);
        else push(-1.0);
        break;

    case LastOp:
        x = pop( );
        if (x >= 1 && x <= ns) push((double) seg[x] -> last -> number);
        else push(-1.0);
        break;

    case LastActiveOp:
        x = pop( );
        if (x >= 1 && x <= ns) push((double) seg[x] -> last_active -> number);
        else push(-1.0);
        break;

    case FirstActiveOp:
        x = pop( );
        if (x >= 1 && x <= ns) push((double) seg[x] -> first_active -> number);
        else push(-1.0);
        break;

	case SinOp:
	    a = pop ( );
	    push (sin (a));
	    break;

	case CosOp:
	    a = pop ( );
	    push (cos (a));
	    break;

	case TanOp:
	    a = pop ( );
	    push (tan (a));
	    break;

	case SinhOp:
	    a = pop ( );
	    push (sinh (a));
	    break;

	case CoshOp:
	    a = pop ( );
	    push (cosh (a));
	    break;

	case TanhOp:
	    a = pop ( );
	    push (tanh (a));
	    break;

	case Atan2Op:
	    b = pop ( );
	    a = pop ( );
	    push (atan2(a,b));
	    break;

	case AcosOp:
	    a = pop ( );
	    push (a >= -1 && a <= 1 ? acos(a) : 0);
	    break;

	case AsinOp:
	    a = pop ( );
	    push (a >= -1 && a <= 1 ? asin(a) : 0);
	    break;

	case ExpOp:
	    a = pop ( );
	    push (exp (a));
	    break;

	case LnOp:
	    a = pop ( );
	    push (a > 0 ? log (a) : 0);
	    break;

	case LogOp:
	    a = pop ( );
	    push (a > 0 ? log10 (a) : 0);
	    break;

	case PowOp:
	    b = pop ( );
	    a = pop ( );
	    push (a >= 0 || b == (int) b ? pow (a, b) : 0);
	    break;

	case SqrtOp:
	    a = pop ( );
	    push (a >= 0 ? sqrt (a) : 0);
	    break;

	case HypotOp:
	    b = pop ( );
	    a = pop ( );
	    push (hypot (a, b));
	    break;

	case FloorOp:
	    a = pop ( );
	    push (floor (a));
	    break;

	case CeilOp:
	    a = pop ( );
	    push (ceil (a));
	    break;

	case FmodOp:
	    b = pop ( );
	    a = pop ( );
	    push (b ? fmod (a, b) : 0);
	    break;

	case FabsOp:
	    a = pop ( );
	    push (fabs (a));
	    break;

	case TableOp:
	    y = pc ++ -> offset;
	    push (EvalTable (pc, y, var, 0));
	    pc += y;
	    break;

	case CycleOp:
	    y = pc ++ -> offset;
	    push (EvalTable (pc, y, var, 1));
	    pc += y;
	    break;

	case HaltOp:
	    return pop ( );
	}
    
    }
}


/************************************************************************
 * Function:	DebugCode						*
 *									*
 * Description:	Print a piece of stack code as instructions.		*
 ************************************************************************/

void DebugCode (code)
    Code code;
{
    int     i;
    int     x;
    Opcode  op;
    Code    pc;


    if (!(pc = code))
	return;

    while (1) {
	op = pc -> op;
	printf ("%lx\t%s", (long) (pc ++), data [op].opcode);

	switch (data [op].arg_type) {
	case Integer:
	    x = pc ++ -> offset;
	    printf ("\t%lx\n", (long) (pc + x));
	    break;

	case Double:
	    printf ("\t%g\n", pc ++ -> arg);
	    break;

	case Array:
	    x = pc ++ -> offset;
	    printf ("\t");
	    for (i = 0; i < x; i ++)
		printf ("%g ", pc ++ -> arg);
	    printf ("\n");
	    break;

	default:
	    printf ("\n");
	}

	if (op == HaltOp)
	    return;
    }
}


/************************************************************************
 * Function:	IsConstant						*
 *									*
 * Description:	Determines if a piece of code is constant.		*
 ************************************************************************/

int IsConstant (code)
    Code code;
{
    Opcode op;
    Code   pc;

    op = 0;	/* gcc -Wall */

    if (!(pc = code))
	    return 1;

    while (1) {
        op = pc ++ -> op;
        op &= 0x3ff;
        if ((op >= VartOp && op <= VarHOp)
            || op == TableOp
            || op == CycleOp) 
            
            return 0;

        if (op == HaltOp)
            return 1;

	    if (data [op].arg_type != None)
		    pc ++;
    }
}

/************************************************************************
 * Function:	IsTable						            
 *									
 * Description:	Determines if a piece of code is constant.		
 ************************************************************************/

int 
IsTable (Code code)
{
    Opcode op;
    Code   pc;

    op = 0;	/* gcc -Wall */

    if (!(pc = code)) {
	    return 0;
    }

    while (1) {
	    switch (op = pc ++ -> op) {
	    case TableOp:
	        return 1;

        case HaltOp:
	        return 0;

	    default:
	        if (data [op].arg_type != None)
		        pc ++;
        }
    }
}

int
DumpCodeTable(Code code, double **table)
{
    Opcode  op;
    Code    pc;
    int     n;

    op = 0;	/* gcc -Wall */

    if (!(pc = code)) {
	    return 1;
    }

    n = 0;
    while (1) {
	    switch (op = pc ++ -> op) {
	    case TableOp:
            n = pc ++ -> offset;
            *table = (double *) pc;
	        return n;

	    default:
	        if (data [op].arg_type != None)
		        pc ++;
        }
    }

    *table = NULL;
    return 0;
}
