function eggfan(lon)
% eggfan(lon)
%
% Extract a fan of sagittal slices over a range of longitudes from a SIM or
% other volumetric dataset.
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 03/16/2004 JMT Adapt from arch.m (JMT)
%
% Copyright (c) 2004-2005 California Institute of Technology.
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

if nargin < 1 lon = -45:5:45; end

% Number of sections in fan
nlon = length(lon);

% Default data directory
eggdir = pwd;
fver = 3.1;

% Splash lines
fprintf('\n');
fprintf('XENOPUS FAN EXTRACTION\n');
fprintf('----------------------\n');
fprintf('Version     : %0.1f\n', fver);
fprintf('Copyright   : 2003 California Institute of Technology\n');
fprintf('Programming : Mike Tyszka\n');
fprintf('\n');

% Allow user to select the egg data file from datadir
[eggname, eggdir] = uigetfile(fullfile(eggdir,'*.mat'));
if isequal(eggname,0) | isequal(eggdir,0)
  fprintf('No egg data selected - exiting\n');
  return
end

% Check for existance of best fit sphere parameters
eggstub = eggname(1:(length(eggname)-4));
parfile = fullfile(eggdir,[eggstub '_bfs.mat']);
if exist(parfile) ~= 2
  fprintf('*** Could not find %s\n', parfile);
  fprintf('*** Run egg_bfs to generate this file\n');
  return
end

% Load BFS parameters
fprintf('Loading BFS parameters from %s\n', parfile);
load(parfile);

% Summarize BFS parameters
fprintf('  Radius : %0.1f um\n', r0);
fprintf('  CofM   : (%0.1f,%0.1f,%0.1f) um\n', CofM(1), CofM(2), CofM(3));

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
% Extract fan of sagittal slices
%----------------------------------------------------

fprintf('Extracting slice fan\n');

fanstub = sprintf('%sFan_%d_%d_%d',eggstub,min(lon),lon(2)-lon(1),max(lon));

for lc = 1:nlon
  
  this_lon = lon(lc);
  
  fprintf('  Slice %d: %0.1f deg\n', lc, this_lon);
  ss = eggplane_bfs(s,vsize,this_lon,r0,CofM,npole,ae,'fan');
  ss = ss / max(ss(:));
  
  fname = sprintf('%s_%03d.tif',fanstub,lc);
  fprintf('  Writing to %s\n', fname);
  imwrite(rot90(ss),fullfile(eggdir,fname),'tif');
  
end