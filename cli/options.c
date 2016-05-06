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

/****************************************************************************
 *
 * File:        options.c
 *
 * Description: routines to read and parse command line options
 *
 * History:
 *
 ****************************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <ctype.h>
# include "allocate.h"
# include "options.h"
# include "error.h"


# define streq(a,b)	 !strcmp(a,b)

static int args_used = 0;

static int 
GetArgNumber (int argc, char *argv[], char *name)
{
   char	buffer [256];
   int	i;
   int  n;

   sprintf (buffer,"-%s", name);      

   n = -1;
   for (i = 1 ; i < argc ; i++) 
      if (streq (argv [i], buffer)) {
         n = i;
         break;
      }
 
   return n;
}

static int 
GetBooleanArgNumber (int argc, char *argv[], char *name, int *state)
{
   char	buffer1 [256];
   char	buffer2 [256];
   int	i;
   int  n;

   sprintf (buffer1, "-%s", name);      
   sprintf (buffer2, "+%s", name);      

   n = -1;
   for (i = 1 ; i < argc ; i++) {
      if (streq (argv [i], buffer1)) {
         n = i;
         *state = 1;
         break;
      }
      else if (streq (argv [i], buffer2)) {
         n = i;
         *state = 0;
         break;
      } 
   }

   return n;
}

static int 
CountArgs (int argc, char *argv[], int n)
{
   int	 i;
   int	 c;

   c = 0;

   for (i = n ; i < argc ; i++) 
      if ((argv [i][0] == '-' || argv [i][0] == '+') && 
          !isdigit(argv [i][1]) && argv [i][1] != '.') {

         c = i - n;
         break;
      }

   if (i == argc)
      c = argc - n;

   return c;
}

int ArgsUsed ( )
{
   return args_used;
}

int GetIntegerOption(argc, argv, name, opt)
   int	 argc;
   char	*argv [ ];
   char	*name;
   int	**opt;
{
   int		n;
   int		c;
   int		i;

   n = GetArgNumber (argc, argv, name);
   if (n == -1) 
      return 0;

   c = CountArgs(argc, argv, n+1);
   if (c == 0)
      return -1;

   *opt = Allocate(int, c);
   if (*opt == NULL)
      Fatal("could not allocate memory for argument options");

   for (i = 0  ; i < c ; i++) 
      (*opt) [i] = atoi (argv [i + n + 1]);

   args_used = args_used + c + 1;

   return c;
}

int GetDoubleOption(argc, argv, name, opt)
   int	     argc;
   char	    *argv [ ];
   char	    *name;
   double  **opt;
{
   int		n;
   int		c;
   int		i;

   n = GetArgNumber (argc, argv, name);
   if (n == -1) 
      return 0;

   c = CountArgs(argc, argv, n+1);
   if (c == 0)
      return -1;

   *opt = Allocate(double, c);
   if (*opt == NULL)
      Fatal("could not allocate memory for argument options");

   for (i = 0  ; i < c ; i++) 
      (*opt) [i] = atof (argv [i + n + 1]);

   args_used = args_used + c + 1;

   return c;
}

int GetFloatOption(argc, argv, name, opt)
   int	     argc;
   char	    *argv [ ];
   char	    *name;
   float   **opt;
{
   int		n;
   int		c;
   int		i;

   n = GetArgNumber (argc, argv, name);
   if (n == -1) 
      return 0;

   c = CountArgs(argc, argv, n+1);
   if (c == 0)
      return -1;

   *opt = Allocate(float, c);
   if (*opt == NULL)
      Fatal("could not allocate memory for argument options");

   for (i = 0  ; i < c ; i++) 
      (*opt) [i] = atof (argv [i + n + 1]);

   args_used = args_used + c + 1;

   return c;
}

int GetBooleanOption(argc, argv, name, opt)
   int	 argc;
   char	*argv [ ];
   char	*name;
   int	*opt;
{
   int		n;
   int		c;
   int		state;

   n = GetBooleanArgNumber (argc, argv, name, &state);
   if (n == -1) 
      return 0;

   c = CountArgs(argc, argv, n+1);
   if (c > 0)
      return -1;

   *opt = state;

   args_used = args_used + c + 1;

   return 1;
}

int GetStringOption(argc, argv, name, opt)
   int	 argc;
   char	*argv [ ];
   char	*name;
   char	***opt;
{
   int		n;
   int		c;
   int		i;
 
   n = GetArgNumber (argc, argv, name);
   if (n == -1) 
      return 0;
 
   c = CountArgs(argc, argv, n+1);
   if (c == 0)
      return -1;

   *opt = Allocate(char *, c);
   if (*opt == NULL)
      Fatal("could not allocate memory for argument options");

   for (i = 0  ; i < c ; i++)
      (*opt) [i] = strdup (argv [i + n + 1]);

   args_used = args_used + c + 1;

   return c;
}

int GetSoloIntegerOption(argc, argv, name, opt)
   int	 argc;
   char	*argv [ ];
   char	*name;
   int	*opt;
{
   int		n;
   int		c;

   n = GetArgNumber (argc, argv, name);
   if (n == -1) 
      return 0;

   c = CountArgs(argc, argv, n+1);
   if (c != 1)
      return -1;

   *opt = atoi (argv [n + 1]);

   args_used = args_used + c + 1;

   return c;
}

int GetSoloFloatOption(argc, argv, name, opt)
   int	    argc;
   char	   *argv [ ];
   char	   *name;
   float   *opt;
{
   int		n;
   int		c;

   n = GetArgNumber (argc, argv, name);
   if (n == -1) 
      return 0;

   c = CountArgs(argc, argv, n+1);
   if (c != 1)
      return -1;

   *opt = atof (argv [n + 1]);

   args_used = args_used + c + 1;

   return c;
}

int GetSoloDoubleOption(argc, argv, name, opt)
   int	    argc;
   char	   *argv [ ];
   char	   *name;
   double  *opt;
{
   int		n;
   int		c;

   n = GetArgNumber (argc, argv, name);
   if (n == -1) 
      return 0;

   c = CountArgs(argc, argv, n+1);
   if (c != 1)
      return -1;

   *opt = atof (argv [n + 1]);

   args_used = args_used + c + 1;

   return c;
}

int GetSoloStringOption(argc, argv, name, opt)
   int	 argc;
   char	*argv [ ];
   char	*name;
   char	**opt;
{
   int		n;
   int		c;
 
   n = GetArgNumber (argc, argv, name);
   if (n == -1) 
      return 0;
 
   c = CountArgs(argc, argv, n+1);
   if (c != 1)
      return -1;

   *opt = strdup (argv [n + 1]);

   args_used = args_used + c + 1;

   return c;
}
