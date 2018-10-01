function eggmaps(sf)
% eggmaps(sf)
%
% Extract shells and save in map,maplegend .MAT files
% Cover radii from 0.75:0.025:0.90
%
% ARGS:
% sf = upsampling scale factor [0.5]
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 10/18/2004 JMT Adapt from eggsnake.m (JMT)
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

% Default shell parameters
r = 0.70:0.02:0.92;
lims = [-90 90 -180 180];

% Default arguments
if nargin < 1 sf = 0.5; end

% Default data directory
fver = 1.0;

% Splash lines
fprintf('\n');
fprintf('EMBRYO MAP GENERATOR\n');
fprintf('--------------------\n');
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
% Extract shells and render orthographic south
% polar view. Save to file
%----------------------------------------------------

% Sampling size in degrees/sample
ds = 2 / (mean(r) * r0) * 180 / (2 * pi) / sf;

nr = length(r);

for rc = 1:nr
  
  % Calculate current r in um
  this_r = r0 * r(rc);

  % Extract shell
  fprintf('Extracting shell radius : %0.3f\n', r(rc));
  [map, maplegend] = eggshell_bfs(s, vsize, this_r, lims, ds, CofM, npole, ae);
  
  % Save figure
  pname = sprintf('%s_%0.3f.mat',eggstub,r(rc));
  fprintf('Saving map to %s\n', pname);
  save(pname,'map','maplegend');

end