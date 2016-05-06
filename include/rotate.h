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

# ifndef _ROTATE_H
# define _ROTATE_H

extern void RotateToFixed (
   double, 
   double,
   double,
   double,
   double,
   double,
   double,
   double *,
   double *,
   double *,
   int
);

extern void RotateToLocal (
   double, 
   double,
   double,
   double,
   double,
   double,
   double,
   double *,
   double *,
   double *
);

# endif /* _ROTATE_H */
