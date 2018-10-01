/************************************************************
 * Allocate and zero space for a uint8 matrix of arbitrary
 * size and dimensionality. Functionally equivalent to zeros()
 *
 * SYNTAX: a = zeros8(dims)
 *
 * AUTHOR : Mike Tyszka, Ph.D.
 * PLACE  : Caltech BIC
 * DATES  : 08/15/2001 JMT Use tricubic8.c as a template
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

typedef unsigned char uint8;

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
  uint8 *a;
  double *dims;

  if (DEBUG) mexPrintf("Starting uint8 allocator\n");
  
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
    mexPrintf("Allocating uint8 space : %d dimensions\n", nAdim);
    for (d = 0; d < nAdim; d++) mexPrintf("%d ", Adims[d]);
    mexPrintf("\n");
  }

  /* Allocate uint8 space */
  A_MAT = mxCreateNumericArray(nAdim, Adims, mxUINT8_CLASS, mxREAL);
  
  /* Clean up */
  mxFree(Adims);

  if (DEBUG) mexPrintf("Done\n");
}
