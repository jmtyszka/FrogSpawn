function eggkeller(sf)
% eggkeller(sf)
%
% Extract a "virtual Keller explant" from SIM data.
% Use the lat range [-90 90], lon range [-45 45] and radii [0.60 1.00]
%
% ARGS:
% sf = upsampling scale factor [0.5]
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 03/16/2004 JMT Adapt from eggseg.m (JMT)
%          11/15/2004 JMT Extend explant to poles
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
if nargin < 1 sf = 0.5; end

% Default data directory
fver = 3.1;

% Splash lines
fprintf('\n');
fprintf('EGG SHELL SEGMENT EXTRACTOR\n');
fprintf('---------------------------\n');
fprintf('Version     : %0.1f\n', fver);
fprintf('Copyright   : 2004 California Institute of Technology\n');
fprintf('Programming : Mike Tyszka\n');
fprintf('Design      : Mike Tyszka and Andy Ewald\n');
fprintf('\n');

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

%----------------------------------------------------
% Extract standard DMZ shells
% ASSUMPTION: Isotropic voxel dimensions
%----------------------------------------------------

dr = 1/r0 / sf;

% Keller shell limits
r = 0.60:dr:1.00;
lims = [-90 90 -45 45];

ds = 2 / (mean(r) * r0) * 180 / (2 * pi) / sf;

% Estimate final shell volume size
nxs = round((lims(2) - lims(1)) / ds);
nys = round((lims(4) - lims(3)) / ds);
nzs = length(r);

fprintf('Extracting shells :\n');
fprintf('  Base radius      : %0.3f um\n', r0);
fprintf('  Fractional Radii : %0.3f to %0.3f\n', min(r), max(r));
fprintf('  Latitudes        : %0.3f to %0.3f\n', lims(1), lims(2));
fprintf('  Longitudes       : %0.3f to %0.3f\n', lims(3), lims(4));
fprintf('  Degrees/sample   : %0.3f\n', ds);
fprintf('  Shell size       : %d x %d x %d\n', nxs, nys, nzs);

for rc = 1:length(r)
  fprintf('Extracting radius %0.3f\n',r(rc));
  [es, maplegend] = eggshell_bfs(s, vsize, r0 * r(rc), lims, ds, CofM, npole, ae);
  if rc == 1
    [nxs,nys] = size(es);
    map = zeros(nxs,nys,nzs);
  end
  map(:,:,rc) = es;
end

  % Save Keller explant in ANALYZE 7.5 format
  mapname = [eggstub 'Keller'];
  fprintf('  Saving map data in %s\n', mapname);
  hdr.vsize = vsize;
  anwrite(fullfile(eggdir,mapname), map, hdr);