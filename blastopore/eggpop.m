function eggpop
% eggpop
%
% Render the best-fit sphere reconstruction from a 2D image of
% a spherical embryo.
%
% Allow user to define the BFC for a plane field embryo image
% then render the image with implied surface topography.
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 10/13/2003 JMT From scratch
%          09/15/2004 JMT Add support for RGB base images
%          10/16/2004 JMT Correct axis orientation in eggpop
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

% User chooses image file
[fname,dname] = uigetfile(...
  {'*.tif;*.png;*.jpg','Image Files'},...
  'Pick an image');

if isequal(fname,0) | isequal(dname,0)
  return
end

% Construct full path to image file
imname = fullfile(dname,fname);

% Initialize BFS parameters
x0 = [];
y0 = [];
a0 = [];

% Create a new color figure
figure(1); clf;
set(gcf,'NumberTitle','off');
set(gcf,'Units','Centimeters','Position',[1 1 15 15]);
set(gcf,'Renderer','zbuffer');

% Load image
s = double(imread(fname));

% Normalize RGB scaling
s = s / max(s(:));

%----------------------------------------------------
% Define embryo border
%----------------------------------------------------
  
keep_going = 1;
    
while keep_going
      
  % Redraw initial image to clear any previous BFS
  figure(1);
  imshow(s);
  
  set(gcf,'Name','Select three points on embryo boundary');
      
  % Initialize point vectors
  x = zeros(3,1);
  y = zeros(3,1);
      
  for p = 1:3
    
    % Get (x,y) coordinate from user
    [x(p),y(p),btn] = ginput(1);
    
    hold on;
    plot(x,y,'og');
    hold off;
    
  end
      
  % Calculate specified circle from the three points
  [a0,x0,y0] = circle3pnt(x,y);
      
  % Overlay the circle
  hold on;
  circle(x0,y0,a0,'g');
  hold off;

  pause;
  
  % Ask user to confirm BFC
  ans = questdlg('Confirmation','Accept Best Fit Circle?','Yes','No','Yes');
      
  switch ans
    case 'Yes'
      keep_going = 0;
    case 'No'
      keep_going = 1;
  end
  
end

%--------------------------------------------------------------
% Create a surface mesh for the image
%--------------------------------------------------------------

% Note axis swap - x is the second coordinate (columns)
[ny,nx,nc] = size(s);
[xc,yc] = meshgrid((1:nx)-x0-1,(1:ny)-y0-1);
r2 = xc.*xc + yc.*yc;
a2 = a0 * a0;
r2(r2 > a2) = a2;
zc = sqrt(a2 - r2);

%--------------------------------------------------------------
% Render the embryo surface in 3D
%--------------------------------------------------------------

figure(2); clf;
h = surface(xc,yc,zc);
set(h,'CData',s,'FaceColor','texturemap','EdgeColor','none');

axis equal vis3d ij off;
camproj('perspective');

set(gcf,'color',[0.8 0.8 0.8]);