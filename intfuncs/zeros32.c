/************************************************************
 * Allocate and zero space for an int32 matrix of arbitrary
 * size and dimensionality. Functionally equivalent to zeros()
 *
 * SYNTAX: a = zeros32(dims)
 *
 * AUTHOR : Mike Tyszka, Ph.D.
 * PLACE  : Caltech BIC
 * DATES  : 09/28/2001 JMT Adapt from zeros8.c
 *
 * Copyright 2001-2005 California Institute of Technology.
 * All rights reserved.
 *
 * This file is part of FrogSpawn.
 *
 * FrogSpawn is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * FrogSpawn is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with FrogSpawn; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 ************************************************************/

#include <stdio.h>
#include <math.h>
#include <mex.h>

#define DEBUG 0

#define A_MAT plhs[0]
#define D_MAT prhs[0]

/************************************************************
 * MAIN ENTRY POINT TO zeros8()
 ************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int *Adims;
  int nAdim, d;
  int *a;
  double *dims;

  if (DEBUG) mexPrintf("Starting int32 allocator\n");
  
  /* Check for proper number of arguments */
  if (nrhs > 1 || nlhs > 1)
	mexErrMsgTxt("SYNTAX: a = zeros(dims)");
  
  /* Get allocation dimensions information */
  nAdim = mxGetNumberOfElements(D_MAT); /* Dimensionality of target */
  dims = mxGetPr(D_MAT); /* Dimensions of target */
  Adims = mxCalloc(nAdim, sizeof(int));
  for (d = 0; d < nAdim; d++) Adims[d] = (int)dims[d];

  /* Optional report */
  if (DEBUG) {
    mexPrintf("Allocating int32 space : %d dimensions\n", nAdim);
    for (d = 0; d < nAdim; d++) mexPrintf("%d ", Adims[d]);
    mexPrintf("\n");
  }

  /* Allocate int32 space */
  A_MAT = mxCreateNumericArray(nAdim, Adims, mxINT32_CLASS, mxREAL);
  
  /* Clean up */
  mxFree(Adims);

  if (DEBUG) mexPrintf("Done\n");
}
