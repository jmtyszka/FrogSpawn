function eggseg(sf)
% eggseg(sf)
%
% Extract a stack of shell segments from an existing dataset to which 
% a sphere has been previously fitted. Best fit sphere data should be
% present in the directory selected by the user at the start of this function.
%
% Allow user to specify shell segment limits (radius, latitude, etc)
%
% ARGS:
% sf = upsampling scale factor (magnification of sampling mesh)
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 01/31/2002 JMT Split out from egg_bfs.m
%          05/30/2002 JMT Simplify interface
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

% Default arguments
if nargin < 1 sf = 1; end

% Default data directory
eggdir = 'C:\Shared Data\Frog Eggs\Frog739';
eggver = 1.3;

% Splash text
fprintf('\n');
fprintf('EGG SHELL SEGMENT EXTRACTOR\n');
fprintf('---------------------------\n');
fprintf('Version     : %0.1f\n', eggver);
fprintf('Copyright   : 2002 California Institute of Technology\n');
fprintf('Programming : Mike Tyszka\n');
fprintf('Design      : Mike Tyszka and Andy Ewald\n');
fprintf('\n');

% Allow user to select the egg data file from datadir
[eggname, eggdir] = uigetfile(fullfile(eggdir,'*.mat'));
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
eggmat = fullfile(eggdir,eggname);
load(eggmat);

% Check for existance of s[]
if isempty(s)
  fprintf('No egg data (s matrix) in %s - exiting\n', eggname);
  return
end

%----------------------------------------------------
% Extract standard DMZ shells
% ASSUMPTION: Isotropic voxel dimensions
%----------------------------------------------------

% Ask user for shell segment parameters
fprintf('\nUSER SHELL SEGMENT PARAMETERS\n');
fprintf('-------------------------------\n');
rA  = input('Start radius         (0 .. 1): ');
rB  = input('End radius           (0 .. 1): ');
laA = input('Start latitude    (-90 .. 90): ');
laB = input('End latitude      (-90 .. 90): ');
loA = input('Start longitude (-180 .. 180): ');
loB = input('End longitude   (-180 .. 180): ');

% Setup range vectors and angular sampling size
% Radial step is 1 voxel / sf
dr = (mean(vsize(:)) / r0) / sf;
r = rA:dr:rB;
lims = [laA laB loA loB];
ds = 2 / (mean(r) * r0) * 180 / (2 * pi) / sf;

% Estimate final shell volume size
nxs = round((lims(2) - lims(1)) / ds);
nys = round((lims(4) - lims(3)) / ds);
nzs = length(r);

fprintf('Extracting shells :\n');
fprintf('  Base radius      : %0.3f voxels\n', r0);
fprintf('  Fractional Radii : %0.3f to %0.3f\n', min(r), max(r));
fprintf('  Latitudes        : %0.3f to %0.3f\n', lims(1), lims(2));
fprintf('  Longitudes       : %0.3f to %0.3f\n', lims(3), lims(4));
fprintf('  Degrees/sample   : %0.3f\n', ds);
fprintf('  Shell size       : %d x %d x %d\n', nxs, nys, nzs);

% Extract shell using specified segment parameters
[map, maplegend] = eggshell_bfs(s,vsize,r0 * r,lims,ds,CofM,npole,ae);

% Save 3D map segment
mapname = fullfile(eggdir,[eggstub '_maps.mat']);
fprintf('  Saving map data in %s\n', mapname);
save(mapname,'map','maplegend','r','lims','ds');