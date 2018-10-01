function [lplanes,xs_e,ys_e,zs_e] = eggplane_bfs(s,vsize,theta,r0,CofM,np,pm,mode)
% [lplanes,xs_e,ys_e,zs_e] = eggplane_bfs(s,vsize,theta,r0,CofM,np,pm,mode)
%
% Extract planes through a given set of latitudes or longitudes from a uint8 egg
% dataset with a predetermined center of mass and north pole location.
% In longitude mode, sagittal planes passing through the given
% longitudes at the equator are extracted from the SIM data.
%
% All coordinates (xc,yc,zc,r0,CofM) are in um.
%
% ARGS:
% s         = 3D egg data
% vsize     = voxel dimensions (um)
% theta     = vector of latitudes or longitudes [degs]
% r0        = best fit sphere radius (um)
% CofM      = egg center of mass (um)
% np        = [lat lon] of north pole in sampling frame [optional]
% pm        = longitude (degs) of prime meridian in egg frame [optional]
% mode      = 'lat','lon' or 'fan' mode
%
% RETURNS:
% lplanes        = stack of extracted planes
% xs_e,ys_e,zs_e = egg frame coordinate grids for extracted planes
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 02/05/2002 JMT Adapt from eggshell_bfs.m
%          07/25/2003 JMT Update with pm argument
%          03/16/2004 JMT Replace voxel center grids with voxel size
%                         Recode axis extrema calculations from vsize
%                         Add fan mode for longitudinal plane extraction
%
% Copyright 2002-2005 California Institute of Technology.
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
nth = length(theta);

% Find extrema of each dimension relative to center of mass
minx = 0.5 * vsize(1) - CofM(1); maxx = (nx - 0.5) * vsize(1) - CofM(1);
miny = 0.5 * vsize(2) - CofM(2); maxy = (ny - 0.5) * vsize(2) - CofM(2);
minz = 0.5 * vsize(3) - CofM(3); maxz = (nz - 0.5) * vsize(3) - CofM(3);

% Axis vectors
xv = linspace(minx,maxx,nx);
yv = linspace(miny,maxy,ny);
zv = linspace(minz,maxz,nz);

% Voxel dimensions
dx = xv(2) - xv(1);
dy = yv(2) - yv(1);
dz = zv(2) - zv(1);

% Default arguments
if nargin < 3 theta = [0]; end
if nargin < 4 r0 = 1; end
if nargin < 5 CofM = [nx/2 ny/2 nz/2]; end
if nargin < 6 np = []; end
if nargin < 7 pm = 0; end
if nargin < 8 mode = 'lat'; end

% Create an xyz coordinate mesh in the egg frame for the extracted planes
switch lower(mode)
  
  case 'lat'
    
    % XY planes at Z values derived from theta
    xrng = xv;
    yrng = yv;
    zrng = r0 * sin(theta * pi/180); 
    [xs_e,ys_e,zs_e] = ndgrid(xrng, yrng, zrng);
    
    % Record extracted plane dimensions
    nxp = nx;
    nyp = ny;
    
  case 'lon'

    % XZ planes at Y values derived from theta
    xrng = xv;
    yrng = r0 * sin(theta * pi/180); 
    zrng = zv;
    [xs_e,ys_e,zs_e] = ndgrid(xrng, yrng, zrng);
    
    % Record extracted plane dimensions
    nxp = nx;
    nyp = nz;
    
  case 'fan'

    % Planes through the poles and a given set of longitudes
    
    % Create a basic sagittal plane
    xrng = xv;
    yrng = 0;
    zrng = zv;
    [xs_e0,ys_e0,zs_e0] = ndgrid(xrng, yrng, zrng);

    % Record extracted plane dimensions
    nxp = nx;
    nyp = nz;
    
    % Create a 3 x nsamp sample point matrix for application of rotations
    xyz_e0 = [xs_e0(:)';ys_e0(:)';zs_e0(:)'];
    
    % Make space for plane sample coordinates
    xs_e = zeros(nxp, nyp, nth);
    ys_e = zeros(nxp, nyp, nth);
    zs_e = zeros(nxp, nyp, nth);
  
    % Rotate plane to each longitude
    for lc = 1:nth
      
      th = theta(lc) * pi / 180;;
      ct = cos(th); st = sin(th);
      Rz = [ct -st 0; st ct 0; 0 0 1];

      % Apply rotations transforming egg frame to sampling frame
      % NOTE: Apply prime meridian rotation first prior to NP Euler rotations
      xyz_e = Rz * xyz_e0;
    
      % Extract and reshape x, y and z coordinate vectors
      xs_e(:,:,lc) = reshape(xyz_e(1,:), [nxp nyp 1]);
      ys_e(:,:,lc) = reshape(xyz_e(2,:), [nxp nyp 1]);
      zs_e(:,:,lc) = reshape(xyz_e(3,:), [nxp nyp 1]);
    
    end
    
  otherwise
    fprintf('Unknown eggplane mode: %s\n', mode);
    return
end

% Make space for extracted planes
lplanes = zeros(nxp, nyp, nth);

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
    
  % Rotate sample plane parallel to prime meridian
  pm = pm * pi /180;
  cpm = cos(pm); spm = sin(pm);
  Rzpm = [cpm -spm 0; spm cpm 0; 0 0 1];
  
  % Create a 3 x nsamp sample point matrix for application of rotations
  xyz_e = [xs_e(:)';ys_e(:)';zs_e(:)'];
  
  % Apply rotations transforming egg frame to sampling frame
  % NOTE: Apply prime meridian rotation first prior to NP Euler rotations
  xyz = Rz * Ry * Rzpm * xyz_e;
    
  % Extract and reshape x, y and z coordinate vectors
  switch mode
    case {'lat','fan'}
      xs = reshape(xyz(1,:), [nxp nyp nth]);
      ys = reshape(xyz(2,:), [nxp nyp nth]);
      zs = reshape(xyz(3,:), [nxp nyp nth]);
    case 'lon'
      xs = permute(reshape(xyz(1,:), [nxp nth nyp]),[1 3 2]);
      ys = permute(reshape(xyz(2,:), [nxp nth nyp]),[1 3 2]);
      zs = permute(reshape(xyz(3,:), [nxp nth nyp]),[1 3 2]);
  end  
end

% Convert from real-world to voxel index coordinates
xi = interp1(xv,1:nx,xs,'linear','extrap');
yi = interp1(yv,1:ny,ys,'linear','extrap');
zi = interp1(zv,1:nz,zs,'linear','extrap');

% Extract the planes from the dataset
lplanes = tricubic8(s,xi,yi,zi);