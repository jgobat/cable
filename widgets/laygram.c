#ifndef lint
static char LayYYsccsid[] = "@(#)yaccpar	1.8 (Berkeley) 01/20/90";
#endif
#define YYBYACC 1
#line 2 "laygram.y"
#include    <X11/Xlib.h>
#include    <X11/Xresource.h>
#include    <stdio.h>
#include    <X11/IntrinsicP.h>
#include    <X11/cursorfont.h>
#include    <X11/StringDefs.h>

#include    <X11/Xmu/Misc.h>
#include    <X11/Xmu/Converters.h>
#include    "LayoutP.h"

static LayoutPtr    *dest;

#line 17 "laygram.y"
typedef union {
    int		    ival;
    XrmQuark	    qval;
    BoxPtr	    bval;
    BoxParamsPtr    pval;
    GlueRec	    gval;
    LayoutDirection lval;
    ExprPtr	    eval;
    Operator	    oval;
} YYSTYPE;
#line 31 "y.tab.c"
#define OC 257
#define CC 258
#define OA 259
#define CA 260
#define OP 261
#define CP 262
#define NAME 263
#define NUMBER 264
#define INFINITY 265
#define VERTICAL 266
#define HORIZONTAL 267
#define EQUAL 268
#define DOLLAR 269
#define PLUS 270
#define MINUS 271
#define TIMES 272
#define DIVIDE 273
#define PERCENTOF 274
#define PERCENT 275
#define WIDTH 276
#define HEIGHT 277
#define UMINUS 278
#define UPLUS 279
#define YYERRCODE 256
short LayYYlhs[] = {                                        -1,
    0,    1,    1,    1,    1,    3,    2,    2,    4,    4,
    5,    5,    7,    7,    8,    8,    6,    6,    6,   10,
   10,   10,   11,   11,   11,   11,   11,   11,   12,   12,
   12,   12,   12,   12,   12,   12,    9,    9,
};
short LayYYlen[] = {                                         2,
    1,    2,    2,    3,    1,    4,    2,    1,    7,    0,
    4,    0,    2,    0,    2,    0,    2,    1,    1,    2,
    2,    1,    2,    2,    3,    2,    1,    2,    3,    3,
    3,    3,    3,    2,    2,    1,    1,    1,
};
short LayYYdefred[] = {                                      0,
   37,   38,    0,    1,    0,    0,    0,    0,   27,    0,
    0,    0,    0,    0,    0,    0,    5,    0,    0,    0,
    0,    0,    0,    0,    0,    2,   28,    0,    0,   23,
   24,    7,    6,    0,    3,   26,   35,   34,   25,    0,
    0,    0,    0,    0,    0,    0,    4,    0,    0,    0,
   31,   32,   33,   19,   13,    0,    0,    0,    0,   17,
   15,    0,   11,    0,    0,    9,
};
short LayYYdgoto[] = {                                       3,
   15,   16,   17,   26,   35,   55,   46,   58,    5,   18,
   22,   23,
};
short LayYYsindex[] = {                                   -256,
    0,    0,    0,    0, -240, -132, -128, -255,    0, -232,
 -227, -227, -218, -206, -132, -250,    0, -200, -210, -128,
 -128, -210, -104, -201, -114,    0,    0, -210, -210,    0,
    0,    0,    0, -201,    0,    0,    0,    0,    0, -128,
 -128, -128, -128, -128, -249, -194,    0, -194, -202, -202,
    0,    0,    0,    0,    0, -213, -249, -193, -175,    0,
    0, -201,    0, -194, -171,    0,
};
short LayYYrindex[] = {                                      0,
    0,    0,    0,    0,    0,    0,    0, -166,    0,    0,
    0,    0,    0,    0, -162,    0,    0, -149, -223,    0,
    0,  -98,    0, -148,    0,    0,    0, -203, -183,    0,
    0,    0,    0, -253,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0, -170,    0, -154, -180, -163,
    0,    0,    0,    0,    0, -239,    0,    0,    0,    0,
    0, -253,    0, -154,    0,    0,
};
short LayYYgindex[] = {                                      0,
    0,   98,  116,    0,    0,   62,  -33,  -34,    0,  100,
   -6,  -18,
};
#define YYTABLESIZE 176
short LayYYtable[] = {                                      19,
   48,   37,   38,   24,   28,   29,   14,   33,   19,    1,
    2,    7,   25,   59,    9,   54,    6,   14,   19,   10,
   18,   49,   50,   51,   52,   53,   13,   14,   64,   65,
   27,   18,   18,    7,   22,   22,    9,   22,   56,   22,
   22,   10,   22,   22,   30,   22,   22,   22,   13,   14,
   56,   60,   22,   22,   21,   21,   31,   21,   34,   21,
   21,   36,   21,   21,   36,   21,   21,   21,   45,   42,
   43,   44,   21,   21,   20,   20,   57,   20,   62,   20,
   20,   29,   20,   20,   63,   20,   20,   20,   66,   29,
   29,   10,   20,   20,   10,    8,   10,   10,   30,   10,
   10,   16,   10,   10,   10,   16,   30,   30,   12,   10,
   10,   12,   32,   12,   12,    4,   12,   12,   61,   12,
   12,   12,   14,   14,   47,    0,   12,   12,    7,    0,
    8,    9,    7,    1,    2,    9,   10,   11,   12,    0,
   10,   20,   21,   13,   14,    0,    7,   13,   14,    9,
    0,    0,    0,    0,   10,   11,   12,   39,    0,    0,
    0,   13,   14,   36,    0,   40,   41,   42,   43,   44,
    0,   36,   36,   36,   36,   36,
};
short LayYYcheck[] = {                                       6,
   34,   20,   21,  259,   11,   12,  260,  258,   15,  266,
  267,  261,  268,   48,  264,  265,  257,  271,   25,  269,
  260,   40,   41,   42,   43,   44,  276,  277,   62,   64,
  263,  271,  272,  261,  258,  259,  264,  261,   45,  263,
  264,  269,  266,  267,  263,  269,  270,  271,  276,  277,
   57,  265,  276,  277,  258,  259,  263,  261,  259,  263,
  264,  275,  266,  267,  275,  269,  270,  271,  270,  272,
  273,  274,  276,  277,  258,  259,  271,  261,  272,  263,
  264,  262,  266,  267,  260,  269,  270,  271,  260,  270,
  271,  258,  276,  277,  261,  258,  263,  264,  262,  266,
  267,  272,  269,  270,  271,  260,  270,  271,  258,  276,
  277,  261,   15,  263,  264,    0,  266,  267,   57,  269,
  270,  271,  271,  272,   25,   -1,  276,  277,  261,   -1,
  263,  264,  261,  266,  267,  264,  269,  270,  271,   -1,
  269,  270,  271,  276,  277,   -1,  261,  276,  277,  264,
   -1,   -1,   -1,   -1,  269,  270,  271,  262,   -1,   -1,
   -1,  276,  277,  262,   -1,  270,  271,  272,  273,  274,
   -1,  270,  271,  272,  273,  274,
};
#define YYFINAL 3
#ifndef YYDEBUG
#define YYDEBUG 0
#endif
#define YYMAXTOKEN 279
#if YYDEBUG
char *LayYYname[] = {
"end-of-file",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"OC","CC","OA","CA","OP","CP",
"NAME","NUMBER","INFINITY","VERTICAL","HORIZONTAL","EQUAL","DOLLAR","PLUS",
"MINUS","TIMES","DIVIDE","PERCENTOF","PERCENT","WIDTH","HEIGHT","UMINUS",
"UPLUS",
};
char *LayYYrule[] = {
"$accept : layout",
"layout : compositebox",
"box : NAME bothparams",
"box : signedExpr oneparams",
"box : NAME EQUAL signedExpr",
"box : compositebox",
"compositebox : orientation OC boxes CC",
"boxes : box boxes",
"boxes : box",
"bothparams : OA opStretch opShrink TIMES opStretch opShrink CA",
"bothparams :",
"oneparams : OA opStretch opShrink CA",
"oneparams :",
"opStretch : PLUS glue",
"opStretch :",
"opShrink : MINUS glue",
"opShrink :",
"glue : simpleExpr INFINITY",
"glue : simpleExpr",
"glue : INFINITY",
"signedExpr : MINUS simpleExpr",
"signedExpr : PLUS simpleExpr",
"signedExpr : simpleExpr",
"simpleExpr : WIDTH NAME",
"simpleExpr : HEIGHT NAME",
"simpleExpr : OP expr CP",
"simpleExpr : simpleExpr PERCENT",
"simpleExpr : NUMBER",
"simpleExpr : DOLLAR NAME",
"expr : expr PLUS expr",
"expr : expr MINUS expr",
"expr : expr TIMES expr",
"expr : expr DIVIDE expr",
"expr : expr PERCENTOF expr",
"expr : MINUS expr",
"expr : PLUS expr",
"expr : simpleExpr",
"orientation : VERTICAL",
"orientation : HORIZONTAL",
};
#endif
#define LayYYclearin (LayYYchar=(-1))
#define LayYYerrok (LayYYerrflag=0)
#ifdef YYSTACKSIZE
#ifndef YYMAXDEPTH
#define YYMAXDEPTH YYSTACKSIZE
#endif
#else
#ifdef YYMAXDEPTH
#define YYSTACKSIZE YYMAXDEPTH
#else
#define YYSTACKSIZE 500
#define YYMAXDEPTH 500
#endif
#endif
int LayYYdebug;
int LayYYnerrs;
int LayYYerrflag;
int LayYYchar;
short *LayYYssp;
YYSTYPE *LayYYvsp;
YYSTYPE LayYYval;
YYSTYPE LayYYlval;
short LayYYss[YYSTACKSIZE];
YYSTYPE LayYYvs[YYSTACKSIZE];
#define LayYYstacksize YYSTACKSIZE
#line 253 "laygram.y"

LayYYwrap ()
{
    return 1;
}

LayYYsetdest (c)
    LayoutPtr	*c;
{
    dest = c;
}
#line 241 "y.tab.c"
#define YYABORT goto LayYYabort
#define YYACCEPT goto LayYYaccept
#define YYERROR goto LayYYerrlab
int
LayYYparse()
{
    register int LayYYm, LayYYn, LayYYstate;
#if YYDEBUG
    register char *LayYYs;
    extern char *getenv();

    if (LayYYs = getenv("YYDEBUG"))
    {
        LayYYn = *LayYYs;
        if (LayYYn >= '0' && LayYYn <= '9')
            LayYYdebug = LayYYn - '0';
    }
#endif

    LayYYnerrs = 0;
    LayYYerrflag = 0;
    LayYYchar = (-1);

    LayYYssp = LayYYss;
    LayYYvsp = LayYYvs;
    *LayYYssp = LayYYstate = 0;

LayYYloop:
    if (LayYYn = LayYYdefred[LayYYstate]) goto LayYYreduce;
    if (LayYYchar < 0)
    {
        if ((LayYYchar = LayYYlex()) < 0) LayYYchar = 0;
#if YYDEBUG
        if (LayYYdebug)
        {
            LayYYs = 0;
            if (LayYYchar <= YYMAXTOKEN) LayYYs = LayYYname[LayYYchar];
            if (!LayYYs) LayYYs = "illegal-symbol";
            printf("LayYYdebug: state %d, reading %d (%s)\n", LayYYstate,
                    LayYYchar, LayYYs);
        }
#endif
    }
    if ((LayYYn = LayYYsindex[LayYYstate]) && (LayYYn += LayYYchar) >= 0 &&
            LayYYn <= YYTABLESIZE && LayYYcheck[LayYYn] == LayYYchar)
    {
#if YYDEBUG
        if (LayYYdebug)
            printf("LayYYdebug: state %d, shifting to state %d\n",
                    LayYYstate, LayYYtable[LayYYn]);
#endif
        if (LayYYssp >= LayYYss + LayYYstacksize - 1)
        {
            goto LayYYoverflow;
        }
        *++LayYYssp = LayYYstate = LayYYtable[LayYYn];
        *++LayYYvsp = LayYYlval;
        LayYYchar = (-1);
        if (LayYYerrflag > 0)  --LayYYerrflag;
        goto LayYYloop;
    }
    if ((LayYYn = LayYYrindex[LayYYstate]) && (LayYYn += LayYYchar) >= 0 &&
            LayYYn <= YYTABLESIZE && LayYYcheck[LayYYn] == LayYYchar)
    {
        LayYYn = LayYYtable[LayYYn];
        goto LayYYreduce;
    }
    if (LayYYerrflag) goto LayYYinrecovery;
#ifdef lint
    goto LayYYnewerror;
#endif
LayYYnewerror:
    LayYYerror("syntax error");
#ifdef lint
    goto LayYYerrlab;
#endif
LayYYerrlab:
    ++LayYYnerrs;
LayYYinrecovery:
    if (LayYYerrflag < 3)
    {
        LayYYerrflag = 3;
        for (;;)
        {
            if ((LayYYn = LayYYsindex[*LayYYssp]) && (LayYYn += YYERRCODE) >= 0 &&
                    LayYYn <= YYTABLESIZE && LayYYcheck[LayYYn] == YYERRCODE)
            {
#if YYDEBUG
                if (LayYYdebug)
                    printf("LayYYdebug: state %d, error recovery shifting\
 to state %d\n", *LayYYssp, LayYYtable[LayYYn]);
#endif
                if (LayYYssp >= LayYYss + LayYYstacksize - 1)
                {
                    goto LayYYoverflow;
                }
                *++LayYYssp = LayYYstate = LayYYtable[LayYYn];
                *++LayYYvsp = LayYYlval;
                goto LayYYloop;
            }
            else
            {
#if YYDEBUG
                if (LayYYdebug)
                    printf("LayYYdebug: error recovery discarding state %d\n",
                            *LayYYssp);
#endif
                if (LayYYssp <= LayYYss) goto LayYYabort;
                --LayYYssp;
                --LayYYvsp;
            }
        }
    }
    else
    {
        if (LayYYchar == 0) goto LayYYabort;
#if YYDEBUG
        if (LayYYdebug)
        {
            LayYYs = 0;
            if (LayYYchar <= YYMAXTOKEN) LayYYs = LayYYname[LayYYchar];
            if (!LayYYs) LayYYs = "illegal-symbol";
            printf("LayYYdebug: state %d, error recovery discards token %d (%s)\n",
                    LayYYstate, LayYYchar, LayYYs);
        }
#endif
        LayYYchar = (-1);
        goto LayYYloop;
    }
LayYYreduce:
#if YYDEBUG
    if (LayYYdebug)
        printf("LayYYdebug: state %d, reducing by rule %d (%s)\n",
                LayYYstate, LayYYn, LayYYrule[LayYYn]);
#endif
    LayYYm = LayYYlen[LayYYn];
    LayYYval = LayYYvsp[1-LayYYm];
    switch (LayYYn)
    {
case 1:
#line 50 "laygram.y"
{ *dest = LayYYvsp[0].bval; }
break;
case 2:
#line 53 "laygram.y"
{
			BoxPtr	box = New(LBoxRec);
			box->nextSibling = 0;
			box->type = WidgetBox;
			box->params = *LayYYvsp[0].pval;
			Dispose (LayYYvsp[0].pval);
			box->u.widget.quark = LayYYvsp[-1].qval;
			LayYYval.bval = box;
		    }
break;
case 3:
#line 63 "laygram.y"
{
			BoxPtr	box = New(LBoxRec);
			box->nextSibling = 0;
			box->type = GlueBox;
			box->params = *LayYYvsp[0].pval;
			Dispose (LayYYvsp[0].pval);
			box->u.glue.expr = LayYYvsp[-1].eval;
			LayYYval.bval = box;
		    }
break;
case 4:
#line 73 "laygram.y"
{
			BoxPtr	box = New(LBoxRec);
			box->nextSibling = 0;
			box->type = VariableBox;
			box->u.variable.quark = LayYYvsp[-2].qval;
			box->u.variable.expr = LayYYvsp[0].eval;
			LayYYval.bval = box;
		    }
break;
case 5:
#line 82 "laygram.y"
{
			LayYYval.bval = LayYYvsp[0].bval;
		    }
break;
case 6:
#line 87 "laygram.y"
{
			BoxPtr	box = New(LBoxRec);
			BoxPtr	child;

			box->nextSibling = 0;
			box->parent = 0;
			box->type = BoxBox;
			box->u.box.dir = LayYYvsp[-3].lval;
			box->u.box.firstChild = LayYYvsp[-1].bval;
			for (child = LayYYvsp[-1].bval; child; child = child->nextSibling) 
			{
			    if (child->type == GlueBox) 
			    {
				child->params.stretch[!LayYYvsp[-3].lval].expr = 0;
				child->params.shrink[!LayYYvsp[-3].lval].expr = 0;
				child->params.stretch[!LayYYvsp[-3].lval].order = 100000;
				child->params.shrink[!LayYYvsp[-3].lval].order = 100000;
				child->params.stretch[!LayYYvsp[-3].lval].value = 1;
				child->params.shrink[!LayYYvsp[-3].lval].value = 1;
			    }
			    child->parent = box;
			}
			LayYYval.bval = box;
		    }
break;
case 7:
#line 113 "laygram.y"
{ 
			LayYYvsp[-1].bval->nextSibling = LayYYvsp[0].bval;
			LayYYval.bval = LayYYvsp[-1].bval;
		    }
break;
case 8:
#line 118 "laygram.y"
{	LayYYval.bval = LayYYvsp[0].bval; }
break;
case 9:
#line 121 "laygram.y"
{	
			BoxParamsPtr	p = New(BoxParamsRec);
			
			p->stretch[LayoutHorizontal] = LayYYvsp[-5].gval;
			p->shrink[LayoutHorizontal] = LayYYvsp[-4].gval;
			p->stretch[LayoutVertical] = LayYYvsp[-2].gval;
			p->shrink[LayoutVertical] = LayYYvsp[-1].gval;
			LayYYval.pval = p;
		    }
break;
case 10:
#line 131 "laygram.y"
{	
			BoxParamsPtr	p = New(BoxParamsRec);
			
			ZeroGlue (p->stretch[LayoutHorizontal]);
			ZeroGlue (p->shrink[LayoutHorizontal]);
			ZeroGlue (p->stretch[LayoutVertical]);
			ZeroGlue (p->shrink[LayoutVertical]);
			LayYYval.pval = p;
		    }
break;
case 11:
#line 142 "laygram.y"
{	
			BoxParamsPtr	p = New(BoxParamsRec);
			
			p->stretch[LayoutHorizontal] = LayYYvsp[-2].gval;
			p->shrink[LayoutHorizontal] = LayYYvsp[-1].gval;
			p->stretch[LayoutVertical] = LayYYvsp[-2].gval;
			p->shrink[LayoutVertical] = LayYYvsp[-1].gval;
			LayYYval.pval = p;
		    }
break;
case 12:
#line 152 "laygram.y"
{	
			BoxParamsPtr	p = New(BoxParamsRec);
			
			ZeroGlue (p->stretch[LayoutHorizontal]);
			ZeroGlue (p->shrink[LayoutHorizontal]);
			ZeroGlue (p->stretch[LayoutVertical]);
			ZeroGlue (p->shrink[LayoutVertical]);
			LayYYval.pval = p;
		    }
break;
case 13:
#line 163 "laygram.y"
{ LayYYval.gval = LayYYvsp[0].gval; }
break;
case 14:
#line 165 "laygram.y"
{ ZeroGlue (LayYYval.gval); }
break;
case 15:
#line 168 "laygram.y"
{ LayYYval.gval = LayYYvsp[0].gval; }
break;
case 16:
#line 170 "laygram.y"
{ ZeroGlue (LayYYval.gval); }
break;
case 17:
#line 173 "laygram.y"
{ LayYYval.gval.order = LayYYvsp[0].ival; LayYYval.gval.expr = LayYYvsp[-1].eval; }
break;
case 18:
#line 175 "laygram.y"
{ LayYYval.gval.order = 0; LayYYval.gval.expr = LayYYvsp[0].eval; }
break;
case 19:
#line 177 "laygram.y"
{ LayYYval.gval.order = LayYYvsp[0].ival; LayYYval.gval.expr = 0; LayYYval.gval.value = 1; }
break;
case 20:
#line 180 "laygram.y"
{
			LayYYval.eval = New(ExprRec);
			LayYYval.eval->type = Unary;
			LayYYval.eval->u.unary.op = LayYYvsp[-1].oval;
			LayYYval.eval->u.unary.down = LayYYvsp[0].eval;
		    }
break;
case 21:
#line 187 "laygram.y"
{ LayYYval.eval = LayYYvsp[0].eval; }
break;
case 23:
#line 191 "laygram.y"
{	LayYYval.eval = New(ExprRec);
			LayYYval.eval->type = Width;
			LayYYval.eval->u.width = LayYYvsp[0].qval;
		    }
break;
case 24:
#line 196 "laygram.y"
{	LayYYval.eval = New(ExprRec);
			LayYYval.eval->type = Height;
			LayYYval.eval->u.height = LayYYvsp[0].qval;
		    }
break;
case 25:
#line 201 "laygram.y"
{ LayYYval.eval = LayYYvsp[-1].eval; }
break;
case 26:
#line 203 "laygram.y"
{
			LayYYval.eval = New(ExprRec);
			LayYYval.eval->type = Unary;
			LayYYval.eval->u.unary.op = LayYYvsp[0].oval;
			LayYYval.eval->u.unary.down = LayYYvsp[-1].eval;
		    }
break;
case 27:
#line 210 "laygram.y"
{	LayYYval.eval = New(ExprRec);
			LayYYval.eval->type = Constant;
			LayYYval.eval->u.constant = LayYYvsp[0].ival;
		    }
break;
case 28:
#line 215 "laygram.y"
{	LayYYval.eval = New(ExprRec);
			LayYYval.eval->type = Variable;
			LayYYval.eval->u.variable = LayYYvsp[0].qval;
		    }
break;
case 29:
#line 221 "laygram.y"
{ binary: ;
			LayYYval.eval = New(ExprRec);
			LayYYval.eval->type = Binary;
			LayYYval.eval->u.binary.op = LayYYvsp[-1].oval;
			LayYYval.eval->u.binary.left = LayYYvsp[-2].eval;
			LayYYval.eval->u.binary.right = LayYYvsp[0].eval;
		    }
break;
case 30:
#line 229 "laygram.y"
{ goto binary; }
break;
case 31:
#line 231 "laygram.y"
{ goto binary; }
break;
case 32:
#line 233 "laygram.y"
{ goto binary; }
break;
case 33:
#line 235 "laygram.y"
{ goto binary; }
break;
case 34:
#line 237 "laygram.y"
{ unary: ;
			LayYYval.eval = New(ExprRec);
			LayYYval.eval->type = Unary;
			LayYYval.eval->u.unary.op = LayYYvsp[-1].oval;
			LayYYval.eval->u.unary.down = LayYYvsp[0].eval;
		    }
break;
case 35:
#line 244 "laygram.y"
{ LayYYval.eval = LayYYvsp[0].eval; }
break;
case 37:
#line 248 "laygram.y"
{   LayYYval.lval = LayoutVertical; }
break;
case 38:
#line 250 "laygram.y"
{   LayYYval.lval = LayoutHorizontal; }
break;
#line 641 "y.tab.c"
    }
    LayYYssp -= LayYYm;
    LayYYstate = *LayYYssp;
    LayYYvsp -= LayYYm;
    LayYYm = LayYYlhs[LayYYn];
    if (LayYYstate == 0 && LayYYm == 0)
    {
#if YYDEBUG
        if (LayYYdebug)
            printf("LayYYdebug: after reduction, shifting from state 0 to\
 state %d\n", YYFINAL);
#endif
        LayYYstate = YYFINAL;
        *++LayYYssp = YYFINAL;
        *++LayYYvsp = LayYYval;
        if (LayYYchar < 0)
        {
            if ((LayYYchar = LayYYlex()) < 0) LayYYchar = 0;
#if YYDEBUG
            if (LayYYdebug)
            {
                LayYYs = 0;
                if (LayYYchar <= YYMAXTOKEN) LayYYs = LayYYname[LayYYchar];
                if (!LayYYs) LayYYs = "illegal-symbol";
                printf("LayYYdebug: state %d, reading %d (%s)\n",
                        YYFINAL, LayYYchar, LayYYs);
            }
#endif
        }
        if (LayYYchar == 0) goto LayYYaccept;
        goto LayYYloop;
    }
    if ((LayYYn = LayYYgindex[LayYYm]) && (LayYYn += LayYYstate) >= 0 &&
            LayYYn <= YYTABLESIZE && LayYYcheck[LayYYn] == LayYYstate)
        LayYYstate = LayYYtable[LayYYn];
    else
        LayYYstate = LayYYdgoto[LayYYm];
#if YYDEBUG
    if (LayYYdebug)
        printf("LayYYdebug: after reduction, shifting from state %d \
to state %d\n", *LayYYssp, LayYYstate);
#endif
    if (LayYYssp >= LayYYss + LayYYstacksize - 1)
    {
        goto LayYYoverflow;
    }
    *++LayYYssp = LayYYstate;
    *++LayYYvsp = LayYYval;
    goto LayYYloop;
LayYYoverflow:
    LayYYerror("yacc stack overflow");
LayYYabort:
    return (1);
LayYYaccept:
    return (0);
}
