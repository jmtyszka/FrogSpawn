function eggshowbfs(sf,r,segrng)
% eggshowbfs(sf,r,segrng)
%
% Extract a mid-sagittal slice and overlay BFS at one or more normalized
% radii (eg [0.80 0.85 0.90]). Optional degree range allows circle
% segment overlay (eg [-80 80] gives segment from -80 to 80 degrees
% relative to equator, [-180 180] would give default full circle). 
%
% ARGS:
% sf = upsampling scale factor [0.5]
% r  = displayed normalized radius - can be a vector of values [1.0]
% segrng = segment range in degrees [-180 180]
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 10/18/2004 JMT From scratch
%          01/06/2005 JMT Add circle segment range
%
% Copyright 2004-2005 California Institute of Technology.
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

% Default arguments
if nargin < 1 sf = 1; end
if nargin < 2 r = 1.0; end
if nargin < 3 segrng = [-180 180]; end

% Allow user to select the egg data file from datadir
[eggname, eggdir] = uigetfile('*.mat');
if isequal(eggname,0) | isequal(eggdir,0)
  fprintf('No egg data selected - exiting\n');
  return
end

% Check for existance of best fit sphere parameters
eggstub = eggname(1:(length(eggname)-4));
parfile = [eggdir eggstub '_bfs.mat'];
if exist(parfile) ~= 2
  fprintf('*** Could not find %s\n', parfile);
  fprintf('*** Run egg_bfs to generate this file\n');
  return
end

% Load BFS parameters
fprintf('Loading BFS parameters from %s\n', parfile);
load(parfile);

% Load uint8 egg data as s[]
fprintf('Loading data from %s\n', eggname);
s = []; % Default
eggmat = [eggdir filesep eggname];
load(eggmat);

% Check for existance of s[]
if isempty(s)
  fprintf('No egg data (s matrix) in %s - exiting\n', eggname);
  return
end

% Grab dataset dimensions
[nx,ny,nz] = size(s);
nth = 1;

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

% Midsagittal plane
xrng = xv;
yrng = 0;
zrng = zv;
[xs_e,ys_e,zs_e] = ndgrid(xrng, yrng, zrng);
    
% Record extracted plane dimensions
nxp = nx;
nyp = nz;

% Make space for extracted planes
lplanes = zeros(nxp, nyp, nth);

% Don't rotate if no north pole is specified (for test shells, etc)
if isempty(npole)
  
  xs = xs_e;
  ys = ys_e;
  zs = zs_e;
  
else
  
  % Calculate rotations required to place supplied north pole at [90 0] in sampling frame
  dth = npole(2);
  dph = 90 - npole(1);
  
  % Convert to radians
  dth = dth * pi / 180;
  dph = dph * pi / 180;
  
  % Rotate sampling grid in cartesian space using rotation matrices about z and y
  ct = cos(dth); st = sin(dth);
  cp = cos(dph); sp = sin(dph);
  Rz = [ct -st 0; st ct 0; 0 0 1];
  Ry = [cp 0 sp; 0 1 0; -sp 0 cp];
    
  % Rotate sample plane parallel to prime meridian
  ae = ae * pi /180;
  cpm = cos(ae); spm = sin(ae);
  Rzpm = [cpm -spm 0; spm cpm 0; 0 0 1];
  
  % Create a 3 x nsamp sample point matrix for application of rotations
  xyz_e = [xs_e(:)';ys_e(:)';zs_e(:)'];
  
  % Apply rotations transforming egg frame to sampling frame
  % NOTE: Apply prime meridian rotation first prior to NP Euler rotations
  xyz = Rz * Ry * Rzpm * xyz_e;
    
  % Extract and reshape x, y and z coordinate vectors
  xs = permute(reshape(xyz(1,:), [nxp nth nyp]),[1 3 2]);
  ys = permute(reshape(xyz(2,:), [nxp nth nyp]),[1 3 2]);
  zs = permute(reshape(xyz(3,:), [nxp nth nyp]),[1 3 2]);
  
end

% Convert from real-world to voxel index coordinates
xi = interp1(xv,1:nx,xs,'linear','extrap');
yi = interp1(yv,1:ny,ys,'linear','extrap');
zi = interp1(zv,1:nz,zs,'linear','extrap');

% Extract the planes from the dataset
lplanes = tricubic8(s,xi,yi,zi);

% Setup figure
figure(1); clf; colormap(gray);
set(gcf,'Position',[0 0 1024 800]);

% Display sagittal slice with BFS overlayed
pcolor(squeeze(xs_e),squeeze(zs_e),lplanes); shading interp; axis equal off;

% Overlay BFS in this slice
% Coordinate grid should be in BFS frame
hold on;
for rs = r
  circleseg(0,0,rs*r0,segrng,'g','linewidth',2);
end
hold off;

% Add title showing radius
minr = min(r(:));
maxr = max(r(:));
title(sprintf('%s  r = %0.3f to %0.3f',eggstub,minr,maxr), 'fontsize', 16);
  
% Save figure
pname = sprintf('%sBFS_%0.3f_%0.3f.png',eggstub,minr,maxr);
fprintf('Printing figure to %s\n', pname);
print('-dpng','-r300',fullfile(eggdir,pname));