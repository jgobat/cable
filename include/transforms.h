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
 * File:	transforms.h						*
 *									*
 * Description:	This file contains useful macros for local -- global	*
 *		coordinate transformations in 3D			*
 ************************************************************************/

# ifndef _TRANSFORMS_H
# define _TRANSFORMS_H

# define XComponent(t,n,b,B0,B1,B2,B3) \
         ((B0*B0 + B1*B1 - B2*B2 - B3*B3)*t \
          + 2.0*(B1*B2 - B0*B3)*n + 2.0*(B1*B3 + B0*B2)*b)

# define dXdB0(t,n,b,B0,B1,B2,B3) \
         (2.0*(B0*t - B3*n + B2*b))
# define dXdB1(t,n,b,B0,B1,B2,B3) \
         (2.0*(B1*t + B2*n + B3*b))
# define dXdB2(t,n,b,B0,B1,B2,B3) \
         (2.0*(-B2*t + B1*n + B0*b))
# define dXdB3(t,n,b,B0,B1,B2,B3) \
         (2.0*(-B3*t - B0*n + B1*b))
# define dXdt(B0,B1,B2,B3) \
         (B0*B0 + B1*B1 - B2*B2 - B3*B3)
# define dXdn(B0,B1,B2,B3) \
         (2.0*(B1*B2 - B0*B3))
# define dXdb(B0,B1,B2,B3) \
         (2.0*(B1*B3 + B0*B2))

# define XtComponent(t,B0,B1,B2,B3) \
         ((B0*B0 + B1*B1 - B2*B2 - B3*B3)*t)

# define YComponent(t,n,b,B0,B1,B2,B3) \
         (2.0*(B1*B2 + B0*B3)*t \
          + (B0*B0 - B1*B1 + B2*B2 - B3*B3)*n \
          + 2.0*(B2*B3 - B0*B1)*b)

# define dYdB0(t,n,b,B0,B1,B2,B3) \
         (2.0*(B3*t + B0*n - B1*b))
# define dYdB1(t,n,b,B0,B1,B2,B3) \
         (2.0*(B2*t - B1*n - B0*b)) 
# define dYdB2(t,n,b,B0,B1,B2,B3) \
         (2.0*(B1*t + B2*n + B3*b))
# define dYdB3(t,n,b,B0,B1,B2,B3) \
         (2.0*(B0*t - B3*n + B2*b))
# define dYdt(B0,B1,B2,B3) \
         (2.0*(B1*B2 + B0*B3))
# define dYdn(B0,B1,B2,B3) \
         (B0*B0 - B1*B1 + B2*B2 - B3*B3)
# define dYdb(B0,B1,B2,B3) \
	 (2.0*(B2*B3 - B0*B1))

# define YtComponent(t,B0,B1,B2,B3) \
         (2.0*(B1*B2 + B0*B3)*t)

# define ZComponent(t,n,b,B0,B1,B2,B3) \
         (2.0*(B1*B3 - B0*B2)*t + 2.0*(B2*B3 + B0*B1)*n \
          + (B0*B0 - B1*B1 - B2*B2 + B3*B3)*b)

# define dZdB0(t,n,b,B0,B1,B2,B3) \
         (2.0*(-B2*t + B1*n + B0*b))
# define dZdB1(t,n,b,B0,B1,B2,B3) \
         (2.0*(B3*t + B0*n - B1*b))
# define dZdB2(t,n,b,B0,B1,B2,B3) \
         (2.0*(-B0*t + B3*n - B2*b))
# define dZdB3(t,n,b,B0,B1,B2,B3) \
         (2.0*(B1*t + B2*n + B3*b))
# define dZdt(B0,B1,B2,B3) \
         (2.0*(B1*B3 - B0*B2))
# define dZdn(B0,B1,B2,B3) \
         (2.0*(B2*B3 + B0*B1))
# define dZdb(B0,B1,B2,B3) \
         (B0*B0 - B1*B1 - B2*B2 + B3*B3)

# define ZtComponent(t,B0,B1,B2,B3) \
         (-2.0*(B0*B2 - B1*B3)*t)

# define tComponent(X,Y,Z,B0,B1,B2,B3) \
         ((B0*B0 + B1*B1 - B2*B2 - B3*B3)*X \
          + 2.0*(B0*B3 + B1*B2)*Y + 2.0*(B1*B3 - B0*B2)*Z)
# define dtdB0(X,Y,Z,B0,B1,B2,B3) \
         (2.0*(B0*X + B3*Y - B2*Z))
# define dtdB1(X,Y,Z,B0,B1,B2,B3) \
         (2.0*(B1*X + B2*Y + B3*Z))
# define dtdB2(X,Y,Z,B0,B1,B2,B3) \
         (2.0*(-B2*X + B1*Y - B0*Z))
# define dtdB3(X,Y,Z,B0,B1,B2,B3) \
         (2.0*(-B3*X + B0*Y + B1*Z))

# define nComponent(X,Y,Z,B0,B1,B2,B3) \
         (2.0*(B1*B2 - B0*B3)*X \
          + (B0*B0 - B1*B1 + B2*B2 - B3*B3)*Y \
          + 2.0*(B0*B1 + B2*B3)*Z)
# define dndB0(X,Y,Z,B0,B1,B2,B3) \
         (2.0*(-B3*X + B0*Y + B1*Z))
# define dndB1(X,Y,Z,B0,B1,B2,B3) \
         (2.0*(B2*X - B1*Y + B0*Z))
# define dndB2(X,Y,Z,B0,B1,B2,B3) \
         (2.0*(B1*X + B2*Y + B3*Z))
# define dndB3(X,Y,Z,B0,B1,B2,B3) \
         (2.0*(-B0*X - B3*Y + B2*Z))

# define bComponent(X,Y,Z,B0,B1,B2,B3) \
         (2.0*(B0*B2 + B1*B3)*X + 2.0*(B2*B3 - B0*B1)*Y \
          + (B0*B0 - B1*B1 - B2*B2 + B3*B3)*Z)
# define dbdB0(X,Y,Z,B0,B1,B2,B3) \
         (2.0*(B2*X - B1*Y + B0*Z))
# define dbdB1(X,Y,Z,B0,B1,B2,B3) \
         (2.0*(B3*X - B0*Y - B1*Z))
# define dbdB2(X,Y,Z,B0,B1,B2,B3) \
         (2.0*(B0*X + B3*Y - B2*Z))
# define dbdB3(X,Y,Z,B0,B1,B2,B3) \
         (2.0*(B1*X + B2*Y + B3*Z))

# endif /* _TRANSFORMS_H */
