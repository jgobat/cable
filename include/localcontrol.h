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
 * File:	control.h						*
 *									*
 * Description:	This file contains the public function and type		*
 *		declarations for the GUI control code (both Win and X)  *
 ************************************************************************/

#ifndef _CABLE_CONTROL_H
#define _CABLE_CONTROL_H

enum {
   PLOT_NONE,
   PLOT_SPACE = 1,
   PLOT_TIME = 2
};

extern int ControlProcessEvents (Solution *);

extern void ControlMessage ( char *, void * );

extern void ControlInfo (
   double,			/* time		*/
   int,				/* step		*/
   double,			/* err		*/
   void *           // controls
);

extern void ControlAuxError (
   double,
   void *           // controls
);

extern void ControlLabels (
   char *,			/* step label 	*/
   char *,			
   char *,
   char *,
   void *           // controls
);

extern void ControlDialogInitialize (Analysis *, void *, int);

extern void ControlMakePlots(Solution *, int);

extern int CreateControlDialog (
   int *,			/* argc		*/
   char **,			/* argv 	*/
   int				/* quiet flag   */
);

extern void ControlMainLoop ( );

extern void ControlTabulateResults(Problem *, Environment *);
extern void ControlPlotResults(Problem *, Environment *);
extern void ControlPlotTime(double, Problem *, Environment *);
extern void ControlPlotSnaps(Problem *, Environment *);
extern void ControlDialogQuit(Solution *);

# endif /* _CONTROL_H */
