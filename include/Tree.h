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
 * File:	Tree.h							*
 *									*
 * Description:	This file contains the public function and type		*
 *		declarations for the red-black trees.			*
 ************************************************************************/

# ifndef _Tree_h
# define _Tree_h
# ifndef _Item_h
# include "Item.h"
# endif

typedef struct tree *Tree;

extern Tree TreeCreate        (ItemComparator);
extern int  TreeDestroy       (Tree);
extern int  TreeIterate       (Tree);
extern int  TreeSize          (Tree);

extern int  TreePreorder      (Tree);
extern int  TreeInorder       (Tree);
extern int  TreePostorder     (Tree);

extern Item TreeInsert        (Tree, Item);
extern Item TreeDelete        (Tree, Item);
extern Item TreeSearch        (Tree, Item);

extern Item TreeMinimum       (Tree);
extern Item TreeMaximum       (Tree);
extern Item TreePredecessor   (Tree, Item);
extern Item TreeSuccessor     (Tree, Item);

extern int  TreeSetIterator   (Tree, ItemIterator, void *);
extern int  TreeSetAndIterate (Tree, ItemIterator, void *);
extern int  TreeSetDestructor (Tree, ItemDestructor);
extern int  TreeSetDuplicator (Tree, ItemDuplicator);
extern void TreeMerge(Tree into, Tree from);

# endif /* _Tree_h */
