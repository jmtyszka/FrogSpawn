function simview(simdir,r)
% simview(simdir,r)
%
% Extract a prime meridian sagittal and south pole orthographic view
%
% ARGS :
% simdir = data directory containing SIM BFS data [pwd]
% r      = normalized radius (0..1) for orthographic view
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 07/24/2003 JMT Adapt from eggview.m (JMT)
%
% Copyright 2003-2005 California Institute of Technology.
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

% Various parameters
fver = 1.0;                % Software version number

if nargin < 1 simdir = pwd; end
if nargin < 2 r = 0.9; end

% Splash lines
fprintf('\n');
fprintf('SIM VIEW EXTRACTION\n');
fprintf('-------------------\n');
fprintf('Version     : %0.1f\n', fver);
fprintf('Copyright   : 2003 California Institute of Technology\n');
fprintf('Programming : Mike Tyszka\n');
fprintf('Design      : Mike Tyszka and Andy Ewald\n');
fprintf('\n');

% Allow user to select the SIM data file from datadir
[simname, simdir] = uigetfile(fullfile(simdir,'*.mat'));
if isequal(simname,0) | isequal(simdir,0)
  fprintf('No SIM data selected - exiting\n');
  return
end

% Check for existance of best fit sphere parameters
simstub = simname(1:(length(simname)-4));
parfile = [simdir simstub '_bfs.mat'];
if exist(parfile) ~= 2
  fprintf('*** Could not find %s\n', parfile);
  fprintf('*** Run egg_bfs to generate this file\n');
  return
end

% Load BFS parameters
fprintf('Loading BFS parameters from %s\n', parfile);
load(parfile);

% Load uint8 SIM data as s[]
fprintf('Loading data from %s\n', simname);
s = []; % Default
simmat = fullfile(simdir,simname);
load(simmat);

% Check for existance of s[]
if isempty(s)
  fprintf('No sim data (s matrix) in %s - exiting\n', simname);
  return
end

%----------------------------------------------------
% Extract orthographic blastopore view
%----------------------------------------------------

% Scale BFS parameters from um to samples
r0v = r0 / mean(vsize);
CofMv = CofM ./ vsize + size(s)/2;

sf = 1; % Scale factor for voxel sampling
drv = 1/r0v;
lims = [-90 0 -180 180];
ds = 2 / (mean(r) * r0v) * 180 / (2 * pi) / sf;

% Estimate final shell volume size
nxs = round((lims(2) - lims(1)) / ds);
nys = round((lims(4) - lims(3)) / ds);
nzs = length(r);

fprintf('Extracting orthographic blastopore view :\n');
fprintf('  Base radius      : %0.3f um\n', r0);
fprintf('  Fractional Radii : %0.3f to %0.3f\n', min(r), max(r));
fprintf('  Latitudes        : %0.3f to %0.3f\n', lims(1), lims(2));
fprintf('  Longitudes       : %0.3f to %0.3f\n', lims(3), lims(4));
fprintf('  Degrees/sample   : %0.3f\n', ds);
fprintf('  Shell size       : %d x %d x %d\n', nxs, nys, nzs);

% Extract shell
% NOTE: All dimensions in voxels, not um
[map, maplegend] = eggshell_bfs(s, r0v * r, lims, ds, CofMv, npole, ae);

% Project the shell orthographically centered at the south pole
figure(1); clf; colormap(gray)
axesm ortho
setm(gca,'origin',[-90 0 0]);
meshm(map,maplegend);

% Export orthographic view as a TIF image
orthoname = fullfile(simdir,[simstub '_ortho.tif']);
fprintf('  Exporting orthographic view to %s\n', orthoname);
print('-r300','-dtiff',orthoname);

%----------------------------------------------------
% Extract prime meridian sagittal view
%----------------------------------------------------

fprintf('Extracting prime meridian sagittal plane :\n');

% Extract sagittal prime meridian plane
% Note all dimensions in um, not voxels for this function
[sagplane, xs_e,ys_e,zs_e] = eggplane_bfs(xc, yc, zc, s, 0, r0, CofM, npole, ae, 'lon');

% Project the shell orthographically centered at the south pole
ncol = 256;
sagplane = sagplane / max(sagplane(:)) * ncol + 1;
figure(2); clf; colormap(gray(ncol))
pcolor(squeeze(xs_e),squeeze(zs_e),sagplane);
axis equal off;
shading interp;

% Export orthographic view as a TIF image
sagname = fullfile(simdir,[simstub '_sag.tif']);
fprintf('  Exporting sagittal view to %s\n', sagname);
print('-r300','-dtiff',sagname);