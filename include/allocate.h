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
 * File:	allocate.h						*
 *									*
 * Description:	This file contains useful macros for memory allocation.	*
 ************************************************************************/

# ifndef _ALLOCATE_H
# define _ALLOCATE_H
# include "malloc.h"

# define AllocNew(type) \
	(type *) malloc (sizeof (type))

# define Allocate(type,number) \
	(type *) malloc (sizeof (type) * (number))

# define Reallocate(ptr,type,number) \
	((ptr) = (ptr) ? (type *) realloc ((ptr), sizeof (type) * (number)) : \
			  (type *) malloc (sizeof (type) * (number)))

# define Deallocate(ptr) \
	if (ptr) free ((char *) (ptr))

# define UnitOffset(ptr) \
	if (ptr) (ptr) --

# define ZeroOffset(ptr) \
	if (ptr) (ptr) ++

# endif /* _ALLOCATE_H */
