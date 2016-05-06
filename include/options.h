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

# ifndef _OPTIONS_H
# define _OPTIONS_H

int ArgsUsed ( void );

int GetIntegerOption (
   int,
   char	**,
   char	*,
   int	**
);

int GetDoubleOption (
   int,
   char	**,
   char	*,
   double **
);

int GetFloatOption (
   int,
   char	**,
   char	*,
   float **
);

int GetBooleanOption (
   int,
   char **,
   char *,
   int *
);

int GetStringOption (
   int,
   char	**,
   char	*,
   char	***
);

int GetSoloIntegerOption (
   int,
   char	**,
   char	*,
   int	*
);

int GetSoloDoubleOption (
   int,
   char	**,
   char	*,
   double *
);

int GetSoloFloatOption (
   int,
   char	**,
   char	*,
   float *
);

int GetSoloStringOption (
   int,
   char	**,
   char	*,
   char	**
);

# endif
