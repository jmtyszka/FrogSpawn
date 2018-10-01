function eggdmz(sf,r,lims)
% eggdmz(sf,r,lims)
%
% Extract dorsal mesendoderm shells and render as globes with viewpoint
% over 60S 0E.
%
% ARGS:
% sf   = upsampling scale factor [0.5]
% r    = normalized radii to extract [0.80:0.05:0.90]
% lims = lat and lon limits of extracted DMZ [-90 90 -180 180]
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 10/20/2004 1.0 JMT Adapt from eggkeller.m (JMT)
%          11/09/2004 1.1 JMT Add r and lims args
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
if nargin < 2 r = 0.80:0.05:0.90; end
if nargin < 3 lims = [-90 90 -180 180]; end

% Default data directory
fver = 1.1;

% Splash lines
fprintf('\n');
fprintf('DORSAL GLOBE EXTRACTOR\n');
fprintf('----------------------\n');
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

  % Render in orthographic projection
  fprintf('Rendering Dorsal Globe\n');
  figure(1); clf; colormap(gray);
  set(gcf,'Position',[0 0 1024 800]);
  
  axesm('mapprojection','ortho','origin',[-60 0 0]);
  meshm(map, maplegend);
  title(sprintf('%s  r = %0.3f',eggstub, r(rc)), 'fontsize', 16);
  drawnow
  
  % Save figure
  pname = sprintf('%sDMZ_%0.3f.png',eggstub,r(rc));
  fprintf('Printing figure to %s\n', pname);
  print('-dpng','-r300',fullfile(eggdir,pname));

end