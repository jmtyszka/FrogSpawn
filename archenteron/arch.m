function arch
% arch
%
% Perform standard quantitation on archenteron angular extent and other
% morphometry of the gastrula stage xenopus embryo.
%
% AUTHOR : Mike Tyszka, Ph.D.
% DESIGN : Mike Tyszka, Ph.D. and Andy Ewald, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 01/31/2002 JMT From scratch
%          03/15/2004 JMT Thicken lines
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

% Default data directory
eggdir = pwd;
fver = 3.1;

% Splash lines
fprintf('\n');
fprintf('XENOPUS ARCHENTERON MORPHOMETRY\n');
fprintf('-------------------------------\n');
fprintf('Version     : %0.1f\n', fver);
fprintf('Copyright   : 2003 California Institute of Technology\n');
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
fprintf('  Radius : %0.1f voxels\n', r0);
fprintf('  CofM   : (%0.1f,%0.1f,%0.1f)\n', CofM(1), CofM(2), CofM(3));

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
% Extract midline sagittal slice through prime meridian
%----------------------------------------------------

theta = 0; % Prime meridian longitude
[ss,xs,ys,zs] = eggplane_bfs(xc,yc,zc,s,theta,r0,CofM,npole,ae,'lon');

% Remove singlet dimensions
xs = squeeze(xs);
zs = squeeze(zs);
ss = squeeze(ss);

figure(1); clf; colormap(gray);
pcolor(xs,zs,ss);
shading interp;
axis equal xy;
xlabel('x'); ylabel('z');

% Draw BFS (NB Image was drawn BFS frame of reference
hold on;
circle(0,0,r0,'g','linewidth',3);
line([0 0],[0 -r0],'color','g','linewidth',3);
hold off;

%----------------------------------------------------
% Ask user for two points within embryo
% Calculate subtended angles assuming archeteron
% to the left of center, blastopore to bottom of frame.
% All angles measured from south pole (ie lat + 90)
%----------------------------------------------------

% Get first point
[x0 z0] = ginput(1)

% Draw angle
line([0 x0],[0 z0],'color','r');

% Report subtended angle
phi0 = atan2(-x0,-z0) * 180 / pi;

% Get second point
[x1 z1] = ginput(1)

% Draw angle
line([0 x1],[0 z1],'color','b');

% Report subtended angle
phi1 = atan2(-x1,-z1) * 180 / pi;
fprintf('Blastopore lip limit : %0.1f degrees\n', phi0);
fprintf('Archenteron limit    : %0.1f degrees\n', phi1);

%----------------------------------------------------
% Write results to file
%----------------------------------------------------
resfile = fullfile(eggdir,[eggstub '_angles.txt']);
fd = fopen(resfile,'w');

if fd < 0
  fprintf('Problem writing to %s\n', resfile);
  return
end

fprintf(fd, 'SIM Study            : %s\n', eggstub);
fprintf(fd, 'Analysis data        : %s\n', datestr(now));
fprintf(fd, 'Blastopore lip angle : %0.1f degrees\n', phi0);
fprintf(fd, 'Archenteron angle    : %0.1f degrees\n', phi1);

fclose(fd);
