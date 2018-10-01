function arch2d
% arch2d
%
% Perform standard quantitation on archenteron angular extent and other
% morphometry of the gastrula stage xenopus embryo.
%
% Takes 3D optical images as data source.
%
% AUTHOR : Mike Tyszka, Ph.D.
% DESIGN : Mike Tyszka, Ph.D. and Andy Ewald, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 06/19/2003 JMT Adapt from arch3d.m (JMT)
%          03/15/2004 JMT Thicken lines and restrict plot to DBP lip, XNOT
%                         limit and AE limit
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

% Default data directory
eggdir = pwd;
fver = 0.1;

% Splash lines
fprintf('\n');
fprintf('XENOPUS ARCHENTERON MORPHOMETRY\n');
fprintf('-------------------------------\n');
fprintf('Dimensions  : 3D\n');
fprintf('Version     : %0.1f\n', fver);
fprintf('Copyright   : 2003 California Institute of Technology\n');
fprintf('Programming : Mike Tyszka\n');
fprintf('Design      : Mike Tyszka and Andy Ewald\n');
fprintf('\n');

% Allow user to select the egg image file
[eggname, eggdir] = uigetfile({'*.tif;*.jpg','Image Files (*.tif,*.jpg)'},'Select an image file');
if isequal(eggname,0) | isequal(eggdir,0)
  fprintf('No image selected - exiting\n');
  return
end

% Report filename to screen
fprintf('------------------------------------------\n');
fprintf('Image filename : %s\n', eggname);
fprintf('------------------------------------------\n');
fprintf('\n');

% Load the image
try
  ss = imread(fullfile(eggdir,eggname));
catch
  fprintf('Problem loading image - exiting\n');
  return
end

%----------------------------------------------------
% Ask user for 3 points which define the BFC
%----------------------------------------------------

keep_going = 1;

figure(1); clf; colormap(gray);
set(gcf,'Renderer','zbuffer','NumberTitle','off');

while keep_going
  
  % (Re)display sagittal midline image

  imagesc(ss);
  shading interp;
  axis equal xy off;
  xlabel('x'); ylabel('y');

  set(gcf,'Name','Select three points on embryo boundary');
  
  % Initialize point vectors
  x = zeros(3,1);
  y = zeros(3,1);
  
  for p = 1:3
  
    % Get (x,y) coordinate from user
    [x(p),y(p),btn] = ginput(1);
    
    hold on;
    plot(x,y,'og','markersize',12,'linewidth',3);
    hold off;
    
  end
  
  % Calculate specified circle from the three points
  [a0,x0,y0] = circle3pnt(x,y);
  
  % Overlay the circle
  hold on;
  circle(x0,y0,a0,'g','linewidth',3);
  hold off;

  % Ask user to confirm BFC
  ans = questdlg('Confirmation','Accept Best Fit Circle?','Yes','No','Yes');
  
  switch ans
    case 'Yes'
      keep_going = 0;
    case 'No'
      keep_going = 1;
  end
  
end
  
%----------------------------------------------------
% Ask user for five points within embryo:
%
% 1. Blastopore center
% 2. Dorsal blastopore lip
% 3. Anterior limit of XNot
% 4. Anterior limit of archenteron
% 5. Leading mesendoderm edge
% 6. Anterior limit of ventral archenteron
%
% Record angles relative to blastopore center with
% postive angles defined by archenteron.
% Ventral BP limit should be a positive acute angle.
%----------------------------------------------------

%----------------------------------------------------
% 1. Blastopore Center
%----------------------------------------------------
set(gcf,'Name','Define blastopore center');
[x_bp y_bp] = ginput(1);

% Draw line from BFC center to blastopore center
line([x0 x_bp],[y0 y_bp],'color',[1 0 0],'linewidth',3);

% By definition
th_bp = 0.0;
fprintf('Blastopore Center                  : %0.1f deg\n', th_bp);

%----------------------------------------------------
% 2. Dorsal Blastopore Lip
%----------------------------------------------------

set(gcf,'Name','Define Dorsal Blastopore Lip');
[x_dbp y_dbp] = ginput(1);

% Draw line from BFC center to dorsal blastopore lip
line([x0 x_dbp],[y0 y_dbp],'color',[0 0 1],'linewidth',3);

% Calculate relative angle
th_dbp = dangle([x0 y0], [x_bp y_bp], [x_dbp y_dbp]);

% Angle from BP to DBP Lip defines positive angular direction.
% If angle is > 180, need to adjust angle to (-angle + 360)
if th_dbp > 180
  ang_sign = -1;
  ang_offset = 360;
else
  ang_sign = 1;
  ang_offset = 0;
end

% Apply sign change
th_dbp = th_dbp * ang_sign + ang_offset;

fprintf('Dorsal Blastopore Lip              : %0.1f deg\n', th_dbp);

%----------------------------------------------------
% 3. Anterior limit of XNot
%----------------------------------------------------

set(gcf,'Name','Define Anterior Limit of XNot Domain');
[x_axn y_axn] = ginput(1);

% Draw line from BFC center to anterior XNot
line([x0 x_axn],[y0 y_axn],'color',[1 1 0],'linewidth',3);

% Calculate relative angle
th_axn = dangle([x0 y0], [x_bp y_bp], [x_axn y_axn]) * ang_sign + ang_offset;
fprintf('Anterior XNot Domain Limit         : %0.1f deg\n', th_axn);

%----------------------------------------------------
% 4. Anterior limit of archenteron
%----------------------------------------------------

set(gcf,'Name','Define Anterior Limit of Archenteron');
[x_ae y_ae] = ginput(1);

% Draw line from BFC center to Anterior Limit of Archenteron
line([x0 x_ae],[y0 y_ae],'color',[1 1 1],'linewidth',3);

% Calculate relative angle
th_ae = dangle([x0 y0], [x_bp y_bp], [x_ae y_ae]) * ang_sign + ang_offset;
fprintf('Anterior Limit of Archenteron      : %0.1f deg\n', th_ae);

%----------------------------------------------------
% 5. Leading Mesendoderm Edge
%----------------------------------------------------

set(gcf,'Name','Define Leading Mesoderm Edge');
[x_me y_me] = ginput(1);

% Draw line from BFC center to leading mesendoderm edge
% line([x0 x_me],[y0 y_me],'color',[0 1 1]);

% Calculate relative angle
th_me = dangle([x0 y0], [x_bp y_bp], [x_me y_me]) * ang_sign + ang_offset;
fprintf('Leading Mesendoderm Edge           : %0.1f deg\n', th_me);

%----------------------------------------------------
% 6. Anterior limit of ventral archenteron
%----------------------------------------------------

set(gcf,'Name','Define Anterior Limit of Ventral Archenteron');
[x_vae y_vae] = ginput(1);

% Draw line from BFC center to anterior ventral archenteron limit
% line([x0 x_vae],[y0 y_vae],'color',[1 1 1]);

% Calculate relative angle
th_vae = dangle([x0 y0], [x_bp y_bp], [x_vae y_vae]) * ang_sign + ang_offset;

% BP to VAE angle will be > 180 and < 360 using this approach
% Remap to positive 0 - 180 domain
if th_vae > 180
  th_vae = abs(th_vae - 360);
end

fprintf('Anterior Ventral Archenteron Limit : %0.1f deg\n', th_vae);

%----------------------------------------------------
% Angle subtractions
%   1. Anterior XNOT - Dorsal BP Lip 
%   2. Anterior AE Limit - Dorsal BP Lip
%----------------------------------------------------

axnot_dbp = th_axn - th_dbp;
ae_dbp    = th_ae - th_dbp;

fprintf('------------------------------------------\n');
fprintf(' Anterior XNOT - Dorsal BP Lip     : %0.1f deg\n', axnot_dbp);
fprintf(' Anterior Arch Lim - Dorsal BP Lip : %0.1f deg\n', ae_dbp);
fprintf('------------------------------------------\n');

%----------------------------------------------------
% Write results to file
%----------------------------------------------------
resfile = fullfile(eggdir,[eggname '.angles.txt']);
fd = fopen(resfile,'w');

if fd < 0
  fprintf('Problem writing to %s\n', resfile);
  return
end

% 1. Blastopore center
% 2. Dorsal blastopore lip
% 3. Anterior limit of XNot
% 4. Anterior limit of archenteron
% 5. Leading mesendoderm edge
% 6. Anterior limit of ventral archenteron

fprintf(fd, 'Image file                      | %s\n', fullfile(eggdir, eggname));
fprintf(fd, 'Analysis data                   | %s\n', datestr(now));
fprintf(fd, 'Best fit circle points (x,y)    | %7.1f | %7.1f | %7.1f | %7.1f | %7.1f | %7.1f\n', ...
  x(1), y(1), x(2), y(2), x(3), y(3));
fprintf(fd, 'BFC (x0,y0,a0)                  | %7.1f | %7.1f | %7.1f\n', x0, y0, a0);
fprintf(fd, 'Blastopore Center (x,y,th)      | %7.1f | %7.1f | %7.1f\n', x_bp,  y_bp,  th_bp);
fprintf(fd, 'Dorsal Blastopore Lip (x,y,th)  | %7.1f | %7.1f | %7.1f\n', x_dbp, y_dbp, th_dbp);
fprintf(fd, 'Anterior XNot (x,y,th)          | %7.1f | %7.1f | %7.1f\n', x_axn, y_axn, th_axn);
fprintf(fd, 'Anterior Archenteron (x,y,th)   | %7.1f | %7.1f | %7.1f\n', x_ae,  y_ae,  th_ae);
fprintf(fd, 'Leading Mesendoderm (x,y,th)    | %7.1f | %7.1f | %7.1f\n', x_me,  y_me,  th_me);
fprintf(fd, 'Anterior Ventral Arch (x,y,th)  | %7.1f | %7.1f | %7.1f\n', x_vae, y_vae, th_vae);
fprintf(fd, 'AXNOT - DBP                     | %7s | %7s | %7.1f\n', '-', '-', axnot_dbp);
fprintf(fd, 'Ant Arch - DBP                  | %7s | %7s | %7.1f\n', '-', '-', ae_dbp);

fclose(fd);
