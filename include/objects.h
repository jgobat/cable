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
 * File:	objects.h						*
 *									*
 * Description:	This file contains the function declarations for the	*
 *		operations on the various objects.			*
 ************************************************************************/

# ifndef _OBJECTS_H
# define _OBJECTS_H
# include "model.h"
# include "code.h"
# include "Tree.h"

extern char   *ExpressionString(VarExpr);
extern void	   AssignCurrent(Environment *, char, Code, char *);
extern void	   AssignSpeed(char, Terminal, Code, char *);
extern void	   AssignThrust(char, Terminal, Code, char *);
extern void	   AssignConnectorThrust(char, Segment, Code, char *);
extern void	   AssignPay(Terminal, Code, char *);
extern void	   AssignSegmentPay(Segment, int, Code, char *);
extern void	   AssignProfile(Terminal, Code, char *);
extern void    AssignTension(int, Material, Code, char *);
extern void    AssignElevation(Environment *, Code, char *);
extern void    AssignSmoothing(Analysis *, Code, char *);
extern void    AssignWind(Environment *, char, Code, char *);
extern void	   AssignElevation(Environment *, Code, char *);
extern void	   AssignCurrentModulation(Environment *, char, Code, char *);
extern void	   AssignDrag(int, Material, Code, char *);

extern Material    CreateMaterial(char *);
extern void        DestroyMaterial(Material);

extern Buoy        CreateBuoy(char *);
extern void        DestroyBuoy(Buoy);

extern Anchor      CreateAnchor(char *);
extern void        DestroyAnchor(Anchor);

extern Connector   CreateConnector(char *);
extern void        DestroyConnector(Connector);

extern Branch  CreateBranch(unsigned);
extern void	   DestroyBranch(Branch);

extern Segment     CreateSegment(unsigned);
extern void        DestroySegment(Segment);

extern Terminal    CreateTerminal(void);
extern void        DestroyTerminal(Terminal);

extern Node        CreateNode(Segment);
extern void        DestroyNode(Node);

# endif /* _OBJECTS_H */
