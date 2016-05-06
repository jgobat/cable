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
 * File:        res2mat.c
 *
 * Description: The driver (which is really all there is) for a program
 *              that reads cable binary results files and writes
 *		a binary Matlab (v4) file containing all of the results.
 *
 * History:
 *
 ****************************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <malloc.h>
# include <string.h>
# include "options.h"

extern int SolutionToMatlab(char *resname, char *matname, int global, int lbs, int totals, int use_file_twoD, int requested_twoD);

void *problem = NULL;
void *environment = NULL;

void Basename(char *);

static int	default_global = 0;
static int 	default_twoD = 1;

static char *usage_msg = "\
usage: res2mat [-in model_results_file] [-out matlab_output_file]\n\
   -global       rotate results into global (XYZ) coordinates (default=off)\n\
   -totals       present results as dynamic + static value (default=off)\n\
   -lbs          convert force units from N to lbf (default=off)\n\
   -ft           convert length units from m to ft (default=off)\n\
   -twoD         use 2D transformation formulas and skip 3D info (default=on)\n\
";

static int usage(msg)
   char *msg;
{
   fprintf (stderr,"res2mat: %s\n\n", msg);
   fputs(usage_msg, stderr);
   exit (1);
}

int main (argc, argv)
   int	 argc;
   char *argv [ ];
{
   char	       **in_name;
   char	         in_buffer [256];
   char		 out_buffer [256];
   char		*out_name;
   int		 input;
   int		 num_in;
   int		 global;
   int		 twoD, twoD_spec;
   int		 totals;
   int		 lbs;
   int	     ft;
   int       n;

   num_in = GetStringOption(argc, argv, "in", &in_name);
   if (num_in <= 0)
      usage ("must specify an input file (-in)");

   out_name = NULL;
   n = GetSoloStringOption(argc, argv, "out", &out_name);
   if (n < 0)
      usage ("error in output file parameter (-out)");

   if (n == 1 && num_in > 1)
      usage ("cannot specify an output name with multiple input files");

   if (num_in == 1 && n == 1 && strcmp(in_name [0], out_name) == 0)
      usage ("input and output file names must be different");

   global = default_global;
   n = GetBooleanOption (argc, argv, "global", &global);
   if (n < 0)
      usage ("error in global boolean option");

   twoD = default_twoD;
   n = GetBooleanOption (argc, argv, "twoD", &twoD);
   if (n < 0)
      usage ("error in twoD boolean option");

   twoD_spec = n;

   lbs = 0;
   n = GetBooleanOption (argc, argv, "lbs", &lbs);
   if (n < 0)      
      usage ("error in lbs boolean option");

   ft = 0;
   n = GetBooleanOption (argc, argv, "ft", &ft);   
   if (n < 0)
      usage ("error in ft boolean option");

   totals = 0;
   n = GetBooleanOption (argc, argv, "totals", &totals);
   if (n < 0)
      usage ("error in totals boolean option");

   for (input = 0 ; input < num_in ; input++) {

      if (out_name) {
         sprintf(in_buffer, "%s", in_name [0]);
         sprintf(out_buffer, "%s", out_name);
      }
      else {
         sprintf(in_buffer, "%s", in_name [input]);
         Basename(in_name [input]);
         if (in_name [input]) {
            sprintf(out_buffer, "%s.mat", in_name [input]);

            fprintf (stderr,"in:  %s\n", in_buffer);
            fprintf (stderr,"out: %s\n", out_buffer);
         }
         else {
            fprintf (stderr,"res2mat: could not form output name\n");
            exit (1);
         }
      }

      SolutionToMatlab(in_buffer, out_buffer,  global, lbs, totals, !twoD_spec, twoD);
   }
 
   exit (0);
}
