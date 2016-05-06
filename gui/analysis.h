/*
    This file is part of WHOI Cable, a program for the static and dynamic
    analysis of oceanographic cable structures.

    Copyright (C) 2008-2016 by Jason Gobat

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

#ifndef _ANALYSIS_H
#define _ANALYSIS_H

enum {
    ANALYSIS_TITLE,
    ANALYSIS_TYPE,
    ANALYSIS_STATIC_RELAX,
    ANALYSIS_STATIC_TOL,
    ANALYSIS_STATIC_IT,
    ANALYSIS_OUTER_RELAX,
    ANALYSIS_OUTER_TOL,
    ANALYSIS_OUTER_IT,
    ANALYSIS_DYNAMIC_RELAX,
    ANALYSIS_DYNAMIC_TOL,
    ANALYSIS_DYNAMIC_IT,
    ANALYSIS_DURATION,
    ANALYSIS_DT,
    ANALYSIS_OUTPUT_DT,
    ANALYSIS_SNAPSHOT_DT,
    ANALYSIS_RAMP,
    ANALYSIS_SOLUTION,
    ANALYSIS_DIMENSION,
    ANALYSIS_LOAD_FILE,
    ANALYSIS_NODES,
    ANALYSIS_TERMINALS,
    ANALYSIS_CONNECTORS,
    ANALYSIS_FIRST,
    ANALYSIS_LAST,
    ANALYSIS_OUTER_METHOD,
    ANALYSIS_INTEGRATION_METHOD,
    ANALYSIS_VIVA_ITERATIONS,
    NUM_ANALYSIS,
};

extern GtkWidget *anparam[NUM_ANALYSIS];

extern GtkWidget *BuildAnalysisParameters(void);
extern void FillAnalysis(Analysis *, Problem *);
extern void ClearAnalysis(void);
extern void FillSolutionControl(Solution *, int);
extern char *AnalysisProblemType(void);
extern void WriteAnalysis(FILE *);

#endif // _ANALYSIS_H
