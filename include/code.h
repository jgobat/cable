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
 * File:	code.h							*
 *									*
 * Description:	This file contains the public function and type		*
 *		declarations for the stack machine code.		*
 ************************************************************************/

# ifndef _CODE_H
# define _CODE_H

typedef union instruction *Code;

#define OLDNODEDATA 2
#define CURRNODEDATA 1

#define FILTERNODEDATA 0x400

typedef enum {
    JmpOp,		/* unconditional jump	    */
    JnzOp,		/* jump if not zero	    */
    JzOp,		/* jump if zero		    */
    PushOp,		/* push value		    */
    PopOp,		/* pop top of stack	    */
    CopyOp,		/* copy top of stack	    */
    TestOp,		/* test top of stack	    */
    NegOp,		/* unary negation	    */
    NotOp,		/* logical negation	    */
    InvOp,		/* bitwise negation	    */
    MulOp,		/* multiplication	    */
    DivOp,		/* division		    */
    ModOp,		/* modulo		    */
    AddOp,		/* addition		    */
    SubOp,		/* subtraction		    */
    LsftOp,		/* left shift		    */
    RsftOp,		/* right shift		    */
    LtOp,		/* less than		    */
    GtOp,		/* greater than		    */
    LteqOp,		/* less than or equal	    */
    GteqOp,		/* greater than or equal    */
    EqOp,		/* equality		    */
    NeqOp,		/* inequality		    */
    AndOp,		/* bitwise and		    */
    XorOp,		/* bitwise xor		    */
    OrOp,		/* bitwise or		    */
    VartOp,
    VarNOp,
    VarSOp,
    VarEOp,
    VarSnOp,
    VarSbOp,
    VarUOp,
    VarVOp,
    VarWOp,
    VarB0Op,
    VarB1Op,
    VarB2Op,
    VarB3Op,
    VarOm1Op,
    VarOm2Op,
    VarOm3Op,
    VarPhiOp,
    VarTOp,		
    VarXOp,		
    VarYOp,		
    VarZOp,		
    VarUxOp,
    VarVyOp,
    VarWzOp,
    VarFxOp,
    VarFyOp,
    VarFzOp,
    VarHOp, 
    VarPOp,
    FirstOp,
    LastOp,
    FirstActiveOp,
    LastActiveOp,
    SinOp,		/* sin function		    */
    CosOp,		/* cos function		    */
    TanOp,		/* tan function		    */
    SinhOp,		/* sinh function	    */
    CoshOp,		/* cosh function	    */
    TanhOp,		/* tanh function	    */
    AsinOp,
    AcosOp,
    Atan2Op,
    ExpOp,		/* exp function		    */
    LnOp,		/* log function		    */
    LogOp,		/* log10 function	    */
    PowOp,		/* pow function		    */
    SqrtOp,		/* sqrt function	    */
    HypotOp,		/* hypot function	    */
    FloorOp,		/* floor function	    */
    CeilOp,		/* ceil function	    */
    FmodOp,		/* fmod function	    */
    FabsOp,		/* fabs function	    */
    TableOp,		/* table of values	    */
    CycleOp,		/* circular table of values */
    HaltOp		/* halt execution	    */
} Opcode;

struct node;

extern Code   InCore;
extern void   EmitCode	   (int, ...);
extern Code   CopyCode     (Code);
extern void   FreeCode     (Code);
extern double EvalCode(Code code, 
                       struct node *n, 
                       double t, double var, 
                       double x, double y, double z, int whichdata);

extern void   DebugCode    (Code);
extern int    CompileCode  (char *);
extern int    IsConstant   (Code);
extern void   SetIP	   (int);
extern int    GetIP	   (void);
extern double exptod	   (char *, char **);
extern int DumpCodeTable(Code, double **);
extern int IsTable(Code);

# endif /* _CODE_H */
