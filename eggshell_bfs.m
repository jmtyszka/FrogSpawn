function [map,maplegend] = eggshell_bfs(s,vsize,r_v,lims,ds,CofM,np,pm)
% [map,maplegend] = eggshell_bfs(s,vsize,r_v,lims,ds,CofM,np,pm)
%
% Extract planes through a given set of  of samples from a uint8 egg
% dataset with a predetermined center of mass and north pole location.
%
% ARGS:
% s         = 3D egg data
% vsize     = voxel dimensions (eg in um)
% r_v       = radius vector (same units as vsize)
% lims      = latitude and longitude limits to extract [LAmin LAmax LOmin LOmax]
% ds        = degrees/sample for shell (isotropic)
% CofM      = egg center of mass in sample units
% np        = [lat lon] of north pole in sampling frame [optional]
% pm        = longitude of prime meridian in sampling frame [optional]
%
% RETURNS:
% map       = 3D egg shell segment
% maplegend = map legend for plotting
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 08/15/2001 JMT Extract from egg.m
%          08/16/2001 JMT Add south pole location transform
%          08/21/2001 JMT Correct polar rotation transform
%                     Rotate egg NP to correct lat, then rotate to correct long
%          09/10/2001 JMT Switch to mapping toolbox rotatem function
%                     Fix code for use of map origin, not north pole
%          01/17/2002 JMT Change to spherical harmonic shell definition
%          01/23/2002 JMT Change to lat and lon vectors to define shell segment
%          01/31/2002 JMT Return to degs/samp and no angular vectors
%                     Add support for r vector and sub shells
%          03/16/2004 JMT Add vsize argument
%
% Copyright 2001-2005 California Institute of Technology.
% All rights reserved.
%
% This file is part of FrogSpawn.
%
% FrogSpawn is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% FrogSpawn is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with FrogSpawn; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% Grab dataset dimensions
[nx,ny,nz] = size(s);

% Default arguments
if nargin < 2 vsize = [1 1 1]; end
if nargin < 3 r_v = nx/2 * mean(vsize); end
if nargin < 4 lims = [-90 90 -180 180]; end
if nargin < 5 ds = 1; end
if nargin < 6 CofM = [nx/2 ny/2 nz/2] .* vsize; end
if nargin < 7 np = []; end
if nargin < 8 pm = 0; end

% Create latitude and longitude vectors
% Remember to include offset for supplied prime meridian
la_v = lims(1):ds:lims(2);
lo_v = (lims(3):ds:lims(4)) + pm;

% Calculate declination, elevation and azimuth vectors from latitude and longitude vector
el_v = la_v * pi / 180;
az_v = lo_v * pi / 180;
de_v = pi/2 - el_v;

% Create a coordinate mesh
[el,az,r] = ndgrid(el_v,az_v,r_v);

% Grab mesh dimensions
[nel,naz,nr] = size(el);

% Make space for extracted shell segment
map = zeros(nel, naz, nr);

% Convert to cartesian coordinates (r is still in um at this stage)
[xs_e,ys_e,zs_e] = sph2cart(az,el,r);

% Scale from distance units to voxels
xs_e = xs_e / vsize(1);
ys_e = ys_e / vsize(2);
zs_e = zs_e / vsize(3);

% Don't rotate if no north pole is specified (for test shells, etc)
if isempty(np)
    
  xs = xs_e;
  ys = ys_e;
  zs = zs_e;
    
else
  
  % Calculate rotations required to place supplied north pole at [90 0] in sampling frame
  dth = np(2);
  dph = 90 - np(1);
  
  % Convert to radians
  dth = dth * pi / 180;
  dph = dph * pi / 180;
  
  % Rotate sampling grid in cartesian space using rotation matrices about z and y
  ct = cos(dth); st = sin(dth);
  cp = cos(dph); sp = sin(dph);
  Rz = [ct -st 0; st ct 0; 0 0 1];
  Ry = [cp 0 sp; 0 1 0; -sp 0 cp];
  
  % Create a 3 x nsamp sample point matrix for application of rotations
  xyz_e = [xs_e(:)';ys_e(:)';zs_e(:)'];
  
  % Apply rotations transforming egg frame to sampling frame
  xyz = Rz * Ry * xyz_e;
  
  % Extract and reshape x, y and z coordinate vectors
  xs = reshape(xyz(1,:), [nel naz nr]);
  ys = reshape(xyz(2,:), [nel naz nr]);
  zs = reshape(xyz(3,:), [nel naz nr]);
  
end

% Recenter coordinate space (remember to scale CofM from distance units to voxels)
xs = xs + CofM(1) / vsize(1);
ys = ys + CofM(2) / vsize(2);
zs = zs + CofM(3) / vsize(3);

% Extract a spherical shell from the original data
map = tricubic8(s,xs,ys,zs);

% Create a legend for use by the Matlab Mapping toolbox
maplegend = [1/ds lims(2) lims(3)];