function [x3,y3,z3] = egg3plane(xc, yc, zc, s, r0, CofM)
% [x3,y3,z3] = egg3plane(xc, yc, zc, s, r0, CofM)
%
% Interactive 3-plane view of egg shell with BFS overlaid
% User can adjust the 3-plane intersection with left mouse button
% clicks. A right mouse button click returns the intersection point.
%
% ARGS :
% xc,yc,zc = voxel center ndgrids 
% s        = 3D egg data (suggest 256x256x256 or smaller)
% r0       = radius of best-fit sphere in s[] (um)
% CofM     = Center of mass of BFS in s[] (um)
%
% RETURNS :
% x3,y3,z3 = final intersection point of three-plane view in BFS frame
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 03/06/2002 JMT From scratch
%          05/09/2003 JMT Replace vsize with voxel center grids
%          03/16/2004 JMT Set initial cursor to mean coordinate
%          10/20/2004 JMT Correct plane labeling (Cor <-> Sag)
%
% Copyright (c) 2002-2005 California Insitute of Technology.
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

% For image display:
% x = 1st index = matrix rows = y in display axes
% y = 2nd index = matrix cols = x in display axes

% Set background to black
set(gcf,'Color','k','Position',[200 50 720 720],'Name','Locate Blastopore Center');

% Grab data dimensions
[nx,ny,nz] = size(s);

% Find extrema
minx = min(xc(:)); maxx = max(xc(:));
miny = min(yc(:)); maxy = max(yc(:));
minz = min(zc(:)); maxz = max(zc(:));

% Create voxel center axis vectors
xcv = linspace(minx,maxx,nx);
ycv = linspace(miny,maxy,ny);
zcv = linspace(minz,maxz,nz);

% Initialize intersection in real space
x3 = (minx + maxx)/2;
y3 = (miny + maxy)/2;
z3 = (minz + maxz)/2;

% Calculate indices into image data
x3i = round(interp1(xcv,1:nx,x3,'nearest','extrap'));
y3i = round(interp1(ycv,1:ny,y3,'nearest','extrap'));
z3i = round(interp1(zcv,1:nz,z3,'nearest','extrap'));

% Set z-buffer rendering
set(gcf,'Renderer','zbuffer');

keepgoing = 1;

while keepgoing
  
  %--------------------------------------
  % Axial Red : matrix(x,y) -> plot(x,y)
  %--------------------------------------
  
  % Extract coordinate and signal slices
  x_ax = squeeze(xc(:,:,z3i));
  y_ax = squeeze(yc(:,:,z3i));
  s_ax = squeeze(s(:,:,z3i));

  subplot(221), pcolor(x_ax,y_ax,s_ax);
  shading flat;
  axis equal xy;
  title('Axial','Color','w');
  xlabel('x (um)'); ylabel('y (um)');
  set(gca,'XLim',[minx maxx],'YLim',[miny maxy],'XColor','r','YColor','r');
  hax = gca;

  % Draw in sagittal and coronal plane intersections
  line([x3 x3],[miny maxy],'color','b'); % Coronal
  line([minx maxx],[y3 y3],'color','g'); % Sagittal
  
  % Overlay BFS in this slice
  rax = sqrt(r0^2 - (z3-CofM(3))^2);
  if isreal(rax)
    hold on;
    circle(CofM(1),CofM(2),rax,'w');
    hold off;
  end
  
  %------------------------------------------
  % Sagittal Green : matrix(x,z) -> plot(x,y)
  %------------------------------------------

  % Extract coordinate and signal slices
  x_cor = squeeze(xc(:,y3i,:));
  z_cor = squeeze(zc(:,y3i,:));
  s_cor = squeeze(s(:,y3i,:));

  subplot(223), pcolor(x_cor,z_cor,s_cor);
  shading flat;
  axis equal xy;
  title('Sagittal','Color','w');
  xlabel('x (um)'); ylabel('z (um)');
  set(gca,'XLim',[minx maxx],'YLim',[minz maxz],'XColor','g','YColor','g');
  hcor = gca;
  
  % Draw in axial and coronal plane intersections
  line([minx maxx],[z3 z3],'color','r'); % Axial
  line([x3 x3],[minz maxz],'color','b'); % Sagittal
  
  % Overlay BFS in this slice
  rsag = sqrt(r0^2 - (y3-CofM(2))^2);
  if isreal(rsag)
    hold on;
    circle(CofM(1),CofM(3),rsag,'w');
    hold off;
  end
  
  %------------------------------------------
  % Coronal Blue : matrix(y,z) -> plot(x,y)
  %------------------------------------------

  % Extract coordinate and signal slices
  y_sag = squeeze(yc(x3i,:,:));
  z_sag = squeeze(zc(x3i,:,:));
  s_sag = squeeze(s(x3i,:,:));

  subplot(224), pcolor(y_sag,z_sag,s_sag);
  shading flat;
  axis equal xy;
  title('Coronal','Color','w');
  xlabel('y (um)'); ylabel('z (um)');
  set(gca,'XLim',[miny maxy],'YLim',[minz maxz],'XColor','b','YColor','b');
  hsag = gca;

  % Draw in axial and sagittal plane intersections
  line([y3 y3],[minz maxz],'color','g'); % Coronal
  line([miny maxy],[z3 z3],'color','r'); % Axial

  % Overlay BFS in this slice
  rcor = sqrt(r0^2 - (x3-CofM(1))^2);
  if isreal(rcor)
    hold on;
    circle(CofM(2),CofM(3),rcor,'w');
    hold off;
  end
  
  drawnow
  
  % Wait for button press and recenter
  [xb,yb,button] = ginput(1);
  xb = round(xb);
  yb = round(yb);
  
  switch button
    
  case 1 % Left mouse button
    
    switch gca
      
    case hax
      x3 = clamp(xb,minx,maxx);
      y3 = clamp(yb,miny,maxy);
      
    case hcor
      x3 = clamp(xb,minx,maxx);
      z3 = clamp(yb,minz,maxz);
      
    case hsag
      y3 = clamp(xb,miny,maxy);
      z3 = clamp(yb,minz,maxz);
      
    otherwise
      
      % Do nothing
      
    end
    
    % Recalculate indices into image data
    x3i = round(interp1(xcv,1:nx,x3,'nearest','extrap'));
    y3i = round(interp1(ycv,1:ny,y3,'nearest','extrap'));
    z3i = round(interp1(zcv,1:nz,z3,'nearest','extrap'));
  
  case 3 % Right mouse button
    
    % Allow exit from while loop
    keepgoing = 0;
    
  end
    
end