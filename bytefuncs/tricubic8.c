/************************************************************
 * Interpolate values from a 3D uint8 dataset using a fast
 * tricubic algorithm.
 *
 * SYNTAX: SI = tricubic8(S, xs, ys, zs)
 *
 * AUTHOR : Mike Tyszka, Ph.D.
 * PLACE  : Caltech BIC
 * DATES  : 08/15/2001 JMT Use mex_unwrap3d.c as a template
 * REFS   : Tricubic interpolation from Graphics Gems V
 *          Chapter III.3 Louis K Arata, pp 107-110
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

#define SI_MAT   plhs[0]

#define S_MAT    prhs[0]
#define XS_MAT   prhs[1]
#define YS_MAT   prhs[2]
#define ZS_MAT   prhs[3]

#define LOC3D(x,y,z) ((x) + nx * ((y) + ny * (z)))
#define xloop for(x=0;x<nx;x++)
#define yloop for(y=0;y<ny;y++)
#define zloop for(z=0;z<nz;z++)
#define iloop for(i=0;i<volsize;i++)

int inbounds(int,		     /* x */
	         int,		     /* y */
	         int,		     /* z */
	         int,		     /* nx */
	         int,		     /* ny */
	         int);		     /* nz */

static void toxyz(int, int *, int *, int *, int, int, int);

double Tricubic(uint8 *, double, double, double, int [3]);
double Trilinear(uint8 *, double, double, double, int [3]);

/************************************************************
 * MAIN ENTRY POINT TO MEX_Unwrap3D()
 ************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int p;
  const int *Sdim, *Xdim, *Ydim, *Zdim;
  int nSdim, nXdim, nYdim, nZdim;
  int nxp, nyp, nzp, np;
  uint8 *s;
  double *si, *xs, *ys, *zs;

  if (DEBUG) mexPrintf("Starting tricubic interpolator\n");
  
  /* Check for proper number of arguments */
  if (nrhs != 4 || nlhs > 1)
	mexErrMsgTxt("SYNTAX: Si = tricubic8(S, xs, ys, zs)");
  
  /* Check that the data is uint8 */
  if (!mxIsClass(S_MAT,"uint8"))
    mexErrMsgTxt("Data must be unsigned 8-bit integers");
    
  /* Check that coordinate vectors are doubles */
  if (!mxIsDouble(XS_MAT) | !mxIsDouble(YS_MAT) | !mxIsDouble(ZS_MAT))
    mexErrMsgTxt("Coordinate vectors must be doubles");

  /* Get matrix dimensions */
  Sdim  = mxGetDimensions(S_MAT);
  Xdim  = mxGetDimensions(XS_MAT);
  Ydim  = mxGetDimensions(YS_MAT);
  Zdim  = mxGetDimensions(ZS_MAT);
  
  /* Get dimensionality of data and coordinate vectors */
  nSdim = mxGetNumberOfDimensions(S_MAT);
  nXdim = mxGetNumberOfDimensions(XS_MAT);
  nYdim = mxGetNumberOfDimensions(YS_MAT);
  nZdim = mxGetNumberOfDimensions(ZS_MAT);

  if (nXdim != nYdim || nXdim != nZdim)
    mexErrMsgTxt("Coordinate matrices must have the same number of dimensions");
  if (nSdim != 3)
    mexErrMsgTxt("Input data must be 3D");
    
  /* Determine the number of coordinates to interpolate */
  nxp = mxGetNumberOfElements(XS_MAT);
  nyp = mxGetNumberOfElements(YS_MAT);
  nzp = mxGetNumberOfElements(ZS_MAT);
  if (nxp != nyp || nxp != nzp)
    mexErrMsgTxt("Coordinate matrices must have the same number of elements");
  np = nxp;

  /* Create real matrices for the interpolated data */
  if (DEBUG) mexPrintf("Making space for interpolated data\n");
  SI_MAT = mxCreateNumericArray(nXdim, Xdim, mxDOUBLE_CLASS, mxREAL);

  /* Get the data pointers */
  s   = (uint8 *)mxGetPr(S_MAT);
  si  = mxGetPr(SI_MAT);
  xs  = mxGetPr(XS_MAT);
  ys  = mxGetPr(YS_MAT);
  zs  = mxGetPr(ZS_MAT);

  /* Loop over all interpolation coordinates */
  if (DEBUG) mexPrintf("Interpolating data at coordinate points\n");
  for (p = 0; p < np; p++) {
    si[p] = Tricubic(s, xs[p], ys[p], zs[p], Sdim);
  }
  
  if (DEBUG) mexPrintf("Done\n");
}

/************************************************************
 * Check whether (x,y,z) is in bounds
 ************************************************************/
int inbounds(int x, int y, int z, int nx, int ny, int nz)
{
  return (x >= 0 && x < nx &&
	  y >= 0 && y < ny &&
	  z >= 0 && z < nz);
}

/************************************************************
 * Convert a vector location into (x,y,z) coordinates
 ************************************************************/
static void toxyz(int loc, int *x, int *y, int *z, int nx, int ny, int nz)
{
  int imsize = nx * ny;
  *z = loc / imsize; loc = loc - *z * imsize;
  *y = loc / nx; loc = loc - *y * nx;
  *x = loc;
}

/**************************************************
 * Tricubic interpolation from Graphics Gems V
 * Chapter III.3 Louis K Arata, pp 107-110
 *
 * This implementation interpolates uint8 data
 * to double values from coordinates passed as
 * doubles (sample units)
 **************************************************/

double Tricubic(uint8 *s, double px, double py, double pz, int *dims)
{
   int x, y, z;
   int i, j, k;
   double dx, dy, dz;
   uint8 *pv;
   double u[4], v[4], w[4];
   double r[4], q[4];
   double vox = 0.0;
   double dx2, dx3;
   double dy2, dy3;
   double dz2, dz3;
   int xDim = dims[0];
   int yDim = dims[1];
   int zDim = dims[2];
   int xyDim = xDim * yDim;

   x = (int)px; y = (int)py; z = (int)pz;

   // Check whether tricubic can be calculated at this location
   if (x <= 1 || x >= (xDim-3) ||
       y <= 1 || y >= (yDim-3) ||
       z <= 1 || z >= (zDim-3)) {
      // Try trilinear instead
      return Trilinear(s, px, py, pz, dims);
   }

   dx = px - (double)x; dy = py - (double)y; dz = pz - (double)z;

   pv = s + (x-1) + (y-1) * xDim + (z-1) * xyDim;

   dx2 = dx * dx; dx3 = dx2 * dx;
   dy2 = dy * dy; dy3 = dy2 * dy;
   dz2 = dz * dz; dz3 = dz2 * dz;

   // Factors for Catmull-Rom tricubic interpolation
   u[0] = -0.5 * dx3 + dx2 - 0.5 * dx;
   u[1] = 1.5 * dx3 - 2.5 * dx2 + 1;
   u[2] = -1.5 * dx3 + 2.0 * dx2 + 0.5 * dx;
   u[3] = 0.5 * dx3 - 0.5 * dx2;
   
   v[0] = -0.5 * dy3 + dy2 - 0.5 * dy;
   v[1] = 1.5 * dy3 - 2.5 * dy2 + 1;
   v[2] = -1.5 * dy3 + 2.0 * dy2 + 0.5 * dy;
   v[3] = 0.5 * dy3 - 0.5 * dy2;
   
   w[0] = -0.5 * dz3 + dz2 - 0.5 * dz;
   w[1] = 1.5 * dz3 - 2.5 * dz2 + 1;
   w[2] = -1.5 * dz3 + 2.0 * dz2 + 0.5 * dz;
   w[3] = 0.5 * dz3 - 0.5 * dz2;

   // Stacked weighted sum
   for (k = 0; k < 4; k++) {
     q[k] = 0.0;
     for (j = 0; j < 4; j++) {
	   r[j] = 0.0;
	   for (i = 0; i < 4; i++) {
	     r[j] += u[i] * (double)(*pv);
	     pv++;
	   }
	   q[k] += v[j] * r[j];
	   pv += xDim - 4;
     }
     vox += w[k] * q[k];
     pv += xyDim - 4 * xDim;
   }
   
   return vox;
}

double Trilinear(uint8 *s, double px, double py, double pz, int dims[3])
{
   int x, y, z;
   int i, j, k;
   double dx, dy, dz;
   uint8 *pv;
   double u[4], v[4], w[4];
   double r[4], q[4];
   double f = 0.0;
   double f000, f001, f010, f011, f100, f101, f110, f111;
   double fx00, fx01, fx10, fx11;
   double fxy0, fxy1;
   double dx2, dx3;
   double dy2, dy3;
   double dz2, dz3;
   int xDim = dims[0];
   int yDim = dims[1];
   int zDim = dims[2];
   int xyDim = xDim * yDim;

   x = (int)px; y = (int)py; z = (int)pz;

   // Check whether trilinear can be calculated at this location
   if (x <= 1 || x >= (xDim-1) ||
       y <= 1 || y >= (yDim-1) ||
       z <= 1 || z >= (zDim-1)) {

      return 0.0;
   }

   dx = px - (double)x; dy = py - (double)y; dz = pz - (double)z;

   // Pointer to lowest corner
   pv = s + (x-1) + (y-1) * xDim + (z-1) * xyDim;

   // Extract samples at corners (fzyx)
   f000 = (double)*pv;
   f001 = (double)*(pv+1);
   f010 = (double)*(pv+xDim);
   f011 = (double)*(pv+xDim+1);
   f100 = (double)*(pv+xyDim);
   f101 = (double)*(pv+xyDim+1);
   f110 = (double)*(pv+xyDim+xDim);
   f111 = (double)*(pv+xyDim+xDim+1);
   
   // Interpolate 4 x edges -> square
   
   fx00 = dx * f001 + (1-dx) * f000;
   fx01 = dx * f011 + (1-dx) * f010;
   fx10 = dx * f101 + (1-dx) * f100;
   fx11 = dx * f111 + (1-dx) * f110;
   
   // Interpolate 2 y edges -> line
   
   fxy0 = dy * fx01 + (1-dy) * fx00;     
   fxy1 = dy * fx11 + (1-dy) * fx10;
   
   // Interpolate 1 z edge -> trilinear interpolated value
   
   f = dz * fxy1 + (1-dz) * fxy0;
   
   return f;
}
